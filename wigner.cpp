#include "wigner.h"

using namespace std;

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

void Wigner::newWaveFunc(WaveFunction *_wav){
    wav = _wav;
    updateGrid();
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

template<typename T>
void Wigner::writeFileVector(std::string filename, std::vector<T> vector){
    std::ofstream file;
    file.open(filename);

    file << interval << std::endl;

    for(int j = 0; j < vector.size(); ++j){
        file << vector[j] << " ";
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

// the Hamiltonian is x^2/2+p^2/2 = a^{+}a + 1/2, write data for the first three eigenstates
void Wigner::oscEigenStates(){
    double L = 5;
    BasisBoundary basis(20,-L,L);
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
    double L = 10;
    LegendreScaled basis(100,-L,L);
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
    std::complex<double> alpha1 = std::complex<double>(-2,0);
    std::complex<double> alpha2 = std::complex<double>(2,0);
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
    double L = 10;
    LegendreScaled basis(50,-L,L);
    WaveFunction wav(&basis);
    wav.set([](double x){ return std::complex<double>(exp(-x*x/2)/pow(M_PI,1.0/4), 0);});
    Wigner wig(&wav, 201);

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
        std::cout << i << " Wigner: " << prob << ", Wavefunc: " << pow(fabs(wig.values[i]),2) << ", theor:" << exp(-pow(wig.interval/(N-1)*(i-N/2),2))/sqrt(M_PI) <<std::endl;
        std::cout << i << std::setprecision(20) << " Ratio: " <<  prob/pow(fabs(wig.values[i]),2) << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Ratio of p probabilities from Wigner divided by the wavefunction amplitude squared:" << std::endl;
    for(int i = 0; i < N; i++){
        double prob = 0;
        if(2*i-N/2 >= 0 && 2*i-N/2 < N){
        for(int j = 0; j < N; j++){
            prob += wig.interval/wig.N*grid(j,2*i-N/2);
        }
        }
        std::cout << i << " Wigner: " << prob << ", DFT: " << pow(fabs(wig.wav->p_on_grid(N)[i]),2) << ", theor:" << exp(-pow(M_PI*N*(2.0*i/(N-1)-1)/(wig.interval),2))/sqrt(M_PI) <<std::endl;
        std::cout << i << std::setprecision(20) << " p:" << M_PI*N*(2.0*i/(N-1)-1)/(wig.interval) << " Ratio: " <<  prob/pow(fabs(wig.wav->p_on_grid(N)[i]),2) << std::endl;
    }
    std::cout << wig.interval << std::endl;

    wig.writeFileVector<std::complex<double>>("pongrid.txt", wig.wav->p_on_grid(N));

    double prob = 0;
    for(int i = 0; i < N; i++){
        prob += wig.interval/N*pow(fabs(wig.values[i]),2);
    }
    std::cout << "Wavefunction sum: " << prob << std::endl;

    wig.writeFile("checkProb.txt");
}

//calculate the time evolution of the wave function wavFunc and plot the Wigner distributions along the way, eps is the x^4 coefficient, T is the final time, N is number of plots
//for the time evolution we find the eigenvectors and eigenvalues of the Hamiltonian operator on the given basis set and then just exponentiate them as usual v(t)=exp(
void Wigner::timeEvo(std::complex<double> (*wavFunc)(double), double perturb, string fileNameStart, double T, int N){
    double L = 7;
    BasisBoundary basis(70,-L,L);
    WaveFunction wav(&basis);
    std::complex<double> alpha(1,1);
    wav.set(wavFunc);
    Wigner wig(&wav, 80);
    wig.writeFile(fileNameStart+".txt");

    //the Hamiltonian matrix in the current basis, we get it from the basis class
    Eigen::MatrixXd ham = basis.matrix("-(d/dx)^2", basis, basis)/2 + basis.matrix("x^2", basis, basis)/2 + perturb*basis.matrix("x^4", basis, basis);

    //find eigenvectors and eigenvalues
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solv(ham, basis.overlap());
//    cout << "Eigenvalues:" << endl;
//    cout << solv.eigenvalues() << endl << endl;

//    cout << "Ham:" << endl;
//    cout << ham << endl << endl;

//    cout << "Overlap matrix: " << endl;
//    cout << basis.overlap() << endl << endl;

//    cout << "Eigenvectors: " << endl;
//    cout << solv.eigenvectors() << endl << endl;

//    cout << "Eigenvector products: " << endl;
//    cout << solv.eigenvectors().adjoint()*solv.eigenvectors() << endl << endl;

    //advance the wavefunction to time t in n steps
//    double T = M_PI/0.5;
//    int N = 20;
    //store the original state, we will just be exponentiating the eigenvector components to find the state at time t
    Eigen::VectorXcd vec = wav.vector;
    //find what combination of eigenstates was the original state
    Eigen::VectorXcd eigVComponents = solv.eigenvectors().fullPivLu().solve(vec);
    for(int n = 0; n < N; n++){
        double t = N > 1 ? T*n/(N) : 0;
        wav.vector = Eigen::VectorXcd::Zero(vec.size());
        for(int i = 0; i < wav.vector.size(); i++){
            //each eigenvector v_n gets just an exponential factor depending on the energy: v_n -> exp(-iE_n t)*v_n
            wav.vector += exp(-std::complex<double>(0,1)*t*solv.eigenvalues()(i)) * eigVComponents(i)*solv.eigenvectors().col(i);
        }
        //calculate the Wigner distribution based on the new basis components calculated above and write the data to a file
        wig.updateGrid();
        wig.writeFile(fileNameStart + to_string(n) + ".txt");
    }
}


//double x1, sig1;
void Wigner::timeEvoTest(double x0, double sig){
    double T = M_PI/0.5;
    double perturb;
    int N = 20;
    complex<double> (*wavFunc)(double) = [](double x) -> complex<double> { complex<double> alpha = complex<double>(1,1); return exp(-pow(x-sqrt(2)*alpha,2)/2.0-pow(alpha.imag(),2))/pow(M_PI,1.0/4);};
    ostringstream out;

    goto next;
    //data with various x^4 perturbation coefficients
    perturb = 0.;
    out.str(string()); out << "timeEvo_eps_" << perturb << "_";
    timeEvo(wavFunc, perturb, out.str(), T, N);

    perturb = 0.01;
    out.str(string()); out << "timeEvo_eps_" << perturb << "_";
    timeEvo(wavFunc, perturb, out.str(), T, N);

    perturb = 0.1;
    out.str(string()); out << "timeEvo_eps_" << perturb << "_";
    timeEvo(wavFunc, perturb, out.str(), T, N);

    perturb = 1.;
    out.str(string()); out << "timeEvo_eps_" << perturb << "_";
    timeEvo(wavFunc, perturb, out.str(), T, N);
next:
    //data for gaussian initial state to observe squeezing
//    x1 = x0;
//    sig1 = sig;
    wavFunc = [](double x) -> complex<double> { double x0 = 1; double sig = 1; return exp(-pow((x-x0)/(sig*sig),2))/pow(M_PI,1.0/4);};
    timeEvo(wavFunc, 0, "timeEvo_x0_1_sig_1_", T, N);

    wavFunc = [](double x) -> complex<double> { double x0 = 1.5; double sig = 1.5; return exp(-pow((x-x0)/(sig*sig),2))/pow(M_PI,1.0/4);};
    timeEvo(wavFunc, 0, "timeEvo_x0_1.5_sig_1.5_", T, N);

    wavFunc = [](double x) -> complex<double> { double x0 = 1.5; double sig = 1.8; return exp(-pow((x-x0)/(sig*sig),2))/pow(M_PI,1.0/4);};
    timeEvo(wavFunc, 0, "timeEvo_x0_1.5_sig_1.8_", T, N);

}
