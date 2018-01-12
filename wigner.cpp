#include "wigner.h"


Wigner::Wigner(WaveFunction *wav, int _N): N(_N){
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

std::vector<std::complex<double>> Wigner::val(int n){
    int toEdge = std::min(n, N-n-1);
    for(int i = -N; i < N; i++){
        if(i < -toEdge || i > toEdge){
            in[N+i][0] = 0;
            in[N+i][1] = 0;
//            std::cout << std::complex<double>(0,0);
        }
        else{
            std::complex<double> value = std::conj(values[n-i])*values[n+i];
            in[N+i][0] = value.real();
            in[N+i][1] = value.imag();
//            std::cout << value;
        }
    }
    std::cout << std::endl;

    for(int i = 0; i < 2*N; i++){
        std::cout << "(" << in[i][0] << "," << in[i][1] << ")";
    }

//    std::cout << std::endl;
    fftw_execute(p);
    std::vector<std::complex<double>> result(N);
//    std::cout << "Full DFT:" << std::endl;
//    for(int i = 0; i < 2*N; i++){
//        std::cout << std::complex<double>(out[i][0], out[i][1]);
//    }
//    std::cout << std::endl;

    std::cout << "F:" << std::endl;
    for(int i = 0; i < (N+1)/2; i++){
        result[i+N/2] = interval/N/M_PI*std::complex<double>(out[2*i][0], out[2*i][1]);
//        std::cout << std::complex<double>(out[2*i][0], out[2*i][1]);
    }
    for(int i = (N+1)/2; i < N; i++){
        result[i-(N+1)/2] = interval/N/M_PI*std::complex<double>(out[2*i][0], out[2*i][1]);
//        std::cout << std::complex<double>(out[2*i][0], out[2*i][1]);
    }
    std::cout << std::endl;

    std::cout << "Result for p:" << std::endl;
    for(int i = 0; i < N; i++){
        std::cout << result[i];
    }

    std::cout << std::endl;
    return result;
}

void Wigner::writeFile(std::string filename){
    std::ofstream file;
    file.open(filename);

    file << interval << std::endl;

    for(int i = 0; i < N; ++i){
        std::vector<std::complex<double>> vals = val(i);
        for(int j = 0; j < N; ++j){
            file << vals[j] << " ";
            std::cout << vals[j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout<< std::endl;
    file.close();
}
