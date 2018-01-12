#include "wigner.h"


Wigner::Wigner(WaveFunction *wav, int _N): N(_N), imagMax(0){
    values = wav->on_grid(N);
    interval = wav->interval;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(2*N));
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(2*N));
    p = fftw_plan_dft_1d(2*N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
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
