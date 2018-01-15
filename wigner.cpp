#include "wigner.h"


Wigner::Wigner(WaveFunction *_wav, int _N): wav(_wav), N(_N), imagMax(0){
    values = wav->on_grid(N);
    interval = wav->interval;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(2*N));
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(2*N));
    p = fftw_plan_dft_1d(2*N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
}

void Wigner::updateGrid(){
    values = wav->on_grid(N);
}

Wigner::~Wigner(){
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
}

std::vector<double> Wigner::val(int n){
    int toEdge = std::min(n, N-n-1);
    for(int i = -N; i < N; i++){
        if(i < -toEdge || i > toEdge){
            in[N+i][0] = 0;
            in[N+i][1] = 0;
        }
        else{
            std::complex<double> value = std::conj(values[n-i])*values[n+i];
            in[N+i][0] = value.real();
            in[N+i][1] = value.imag();
        }
    }

    fftw_execute(p);
    std::vector<double> result(N);

    for(int i = 0; i < (N+1)/2; i++){
        result[i+N/2] = interval/N/M_PI*out[2*i][0];
        if(fabs(out[2*i][1]) > imagMax) imagMax = fabs(out[2*i][1]);
    }
    for(int i = (N+1)/2; i < N; i++){
        result[i-(N+1)/2] = interval/N/M_PI*out[2*i][0];
        if(fabs(out[2*i][1]) > imagMax) imagMax = fabs(out[2*i][1]);
    }

    if(imagMax > 1e-10){
        ABORT("GETTING IMAGINARY VALUES FOR THE WIGNER DISTRIBUTION");
    }

    return result;
}

void Wigner::writeFile(std::string filename){
    std::ofstream file;
    file.open(filename);

    file << interval << std::endl;

    for(int i = 0; i < N; ++i){
        std::vector<double> vals = val(i);
        for(int j = 0; j < N; ++j){
            file << vals[j] << " ";
        }
    }
    file.close();
}

void Wigner::writeFileSum(std::string filename, std::vector<Wigner*> wigs){
    Wigner* wig0 = wigs[0];
    std::ofstream file;
    file.open(filename);

    file << wig0->interval << std::endl;

    for(int i = 0; i < wig0->N; ++i){
        std::vector<double> vals = wigs[0]->val(i);
        for(int m = 1; m < wigs.size(); ++m){
            std::vector<double> vals2 = wigs[m]->val(i);
            for(int j = 0; j < wig0->N; ++j){
                vals[j] += vals2[j];
            }
        }
        for(int j = 0; j < wig0->N; ++j){
            file << vals[j] << " ";
        }
    }
    file.close();
}

// a simple test, with this wavefunction the Wigner function should be 1/pi*exp(-x^2-p^2), so the contours must be circles
void Wigner::test(){
    double L = 5;
    LegendreScaled basis(60,-L,L);
    WaveFunction wav(&basis);
    wav.set([](double x){ return std::complex<double>(exp(-x*x/2)/pow(M_PI,1.0/4), 0);});
    int n = 101;
    Wigner wig(&wav, n);

    wig.writeFile("test1.txt");
}

// the Hamiltonian is x^2/2+p^2/2 = a^{+}a + 1/2
void Wigner::oscEigenStates(){
    double L = 5;
    LegendreScaled basis(20,-L,L);
    WaveFunction wav(&basis);
    Wigner wig(&wav, 50);

    wav.set([](double x){ return std::complex<double>(exp(-x*x/2)/pow(M_PI,1.0/4), 0);});
    wig.updateGrid();
    wig.writeFile("oscEigen0.txt");

    wav.set([](double x){ return std::complex<double>(sqrt(2)*x*exp(-x*x/2)/pow(M_PI,1.0/4), 0);});
    wig.updateGrid();
    wig.writeFile("oscEigen1.txt");

    wav.set([](double x){ return std::complex<double>((2*x*x-1)*exp(-x*x/2)/pow(M_PI,1.0/4), 0);});
    wig.updateGrid();
    wig.writeFile("oscEigen2.txt");
}

// coherent states are for any complex \alpha: 1/pow(M_PI,0.25)*exp[-(x-\sqrt{2}\alpha)^2/2-Im(\alpha)^2]
// or by changing global phase can also be written as 1/pow(M_PI,0.25)*exp[-(x-\sqrt{2}Re(\alpha))^2/2+i\sqrt{2}Im(\alpha)x]
// then <psi|x|psi>=\sqrt{2}Re(\alpha) and <psi|p|psi>=\sqrt{2}Im(\alpha) and time evolution just rotates \alpha->exp(i\phi)\alpha
void Wigner::oscCoherent(){
    double L = 7;
    LegendreScaled basis(50,-L,L);
    WaveFunction wav(&basis);
    Wigner wig(&wav, 80);
    std::complex<double> alpha;

    //plot for various alpha values
    struct dataSet{
        std::complex<double> alpha;
        std::string filename;
    };

    std::vector<dataSet> plots = {
        {std::complex<double>(0,0), "oscCoh00.txt"},
        {std::complex<double>(0,1), "oscCoh01.txt"},
        {std::complex<double>(1,0), "oscCoh10.txt"},
        {std::complex<double>(1,1), "oscCoh11.txt"},
    };

    for(int i = 0; i < plots.size(); i++){
        dataSet d = plots[i];
        std::complex<double> alpha = d.alpha;
        wav.set([&](double x){ return exp(-pow(x-sqrt(2)*alpha,2)/2.0-pow(alpha.imag(),2))/pow(M_PI,1.0/4);});
        wig.updateGrid();
        wig.writeFile(d.filename);
    }

    //plot the pure state (|a1> + |a2>)/sqrt(2)
    std::complex<double> alpha1 = std::complex<double>(-1.5,0);
    std::complex<double> alpha2 = std::complex<double>(1.5,0);
    wav.set([&](double x){ return 1/sqrt(2)*(exp(-pow(x-sqrt(2)*alpha1,2)/2.0-pow(alpha1.imag(),2))/pow(M_PI,1.0/4) +
                                                 exp(-pow(x-sqrt(2)*alpha2,2)/2.0-pow(alpha2.imag(),2))/pow(M_PI,1.0/4));});
    wig.updateGrid();
    wig.writeFile("oscPure.txt");

    //plot the mixed state with the same weights
    wav.set([&](double x){ return 1/sqrt(2)*exp(-pow(x-sqrt(2)*alpha1,2)/2.0-pow(alpha1.imag(),2))/pow(M_PI,1.0/4);});
    wig.updateGrid();
    WaveFunction wav2(&basis);
    Wigner wig2(&wav2, 80);
    wav2.set([&](double x){ return 1/sqrt(2)*exp(-pow(x-sqrt(2)*alpha2,2)/2.0-pow(alpha2.imag(),2))/pow(M_PI,1.0/4);});
    wig2.updateGrid();
    std::vector<Wigner*> wigs = {&wig,&wig2};
    Wigner::writeFileSum("oscMixed.txt", wigs);
}

//calculate that the integral of W over x or p gives the probabilities for p and x respectively
void Wigner::checkProb(){
    double L = 14;
    LegendreScaled basis(50,-L,L);
    WaveFunction wav(&basis);
    wav.set([](double x){ return std::complex<double>(exp(-x*x/2)/pow(M_PI,1.0/4), 0);});
    Wigner wig(&wav, 250);

    int N = wig.N;

    Eigen::MatrixXd grid(N,N);
    for(int i = 0; i < N; i++){
        std::vector<double> vals = wig.val(i);
        for(int j = 0; j < N; j++){
            grid(i,j) = vals[j];
        }
    }

    std::cout << "Ratio of x probabilities from Wigner divided by the wavefunction amplitude squared:" << std::endl;
    for(int i = 0; i < N; i++){
        double prob = 0;
        for(int j = 0; j < N; j++){
            prob += M_PI/wig.interval*grid(i,j);
        }
        std::cout << std::setprecision(20) << "Ratio: " <<  prob/pow(fabs(wig.values[i]),2) << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Ratio of p probabilities from Wigner divided by the wavefunction amplitude squared:" << std::endl;
    for(int i = 0; i < N; i++){
        double prob = 0;
        for(int j = 0; j < N; j++){
            prob += wig.interval/wig.N*grid(j,i);
        }
        std::cout << i << " Wigner: " << prob << ", DFT: " << pow(fabs(wig.wav->p_on_grid(N)[i]),2) << ", theor:" << exp(-pow(2*M_PI*(i-(N+1)/2)/wig.interval,2))/sqrt(M_PI) <<std::endl;
        std::cout << i << std::setprecision(20) << " Ratio: " <<  prob/pow(fabs(wig.wav->p_on_grid(N)[i]),2) << std::endl;
    }
    std::cout << std::endl;

    double prob = 0;
    for(int i = 0; i < N; i++){
        prob += wig.interval/N*pow(fabs(wig.values[i]),2);
    }
    std::cout << "Wavefunction sum: " << prob << std::endl;

    wig.writeFile("checkProb.txt");
}
