#include "wigner.h"


Wigner::Wigner(const std::function<std::complex<double>(double)>& Psi, int N) : basis(N,0,1), _Psi(Psi), _N(N), wav(&basis){

    values = wav.on_grid(0,1,1.0/_N);

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    std::complex<double> mycomplex(1,2);
    fftw_complex c = {1,2};

    fftw_execute(p); /* repeat as needed */
}

Wigner::~Wigner(){
    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
}

std::vector<std::complex<double>> Wigner::val(int n){
    for(int i = 0; i < _N; i++){
        std::complex<double> val1 = n+i < _N ? values[n+i] : 0;
        std::complex<double> val2 = n-i > 0 ? values[n-i] : 0;
        std::complex<double> val = std::conj(val1)*val2;
        in[i][0] = val.real();
        in[i][0] = val.imag();
    }
    fftw_execute(p);
    std::vector<std::complex<double>> result(_N);
    for(int i = 0; i < _N; i++){
        result[i] = std::complex<double>(out[i][0], out[i][1]);
    }
    return result;
}
