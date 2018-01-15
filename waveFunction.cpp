#include "waveFunction.h"

WaveFunction::WaveFunction(const BasisSet* Basis): basis(Basis), vector(Basis->size()){}
WaveFunction::WaveFunction(const BasisSet* Basis, const Eigen::VectorXcd& Vector): basis(Basis), vector(Vector){}

//here vector gives the coefficient of the i-th basis function, so that psi= v_i*psi_i
void WaveFunction::set(const std::function<std::complex<double>(double)>& Psi){
    std::vector<double> x, w;
    basis->quadrature(x, w);

    Eigen::MatrixXcd F(basis->size(), x.size());
    for(int k=0; k<x.size(); k++){
        std::vector<double> val = basis->val(x[k]);
        for(int i=0; i<val.size(); i++){
            F(i, k) = std::conj(val[i]);
        }
    }

    Eigen::VectorXcd psi(x.size());
    for(int k=0; k<x.size(); k++){
        psi(k) = w[k]*Psi(x[k]);
    }

    vector = basis->overlap().inverse()*F*psi;
}

std::vector<std::complex<double>> WaveFunction::on_grid(int N, double left, double right){
    std::vector<std::complex<double>> res;
    if(left == 0 && right == 0){
        left = basis->lowerBoundary();
        right = basis->upperBoundary();
    }
    interval = right-left;
    double x = left;
    for(int i=0; i < N; i++){
//        std::cout << x << std::endl;
        std::vector<double> val = basis->val(x);
        res.push_back(0);
        for(int i=0; i<val.size(); i++) res.back()+=val[i]*vector(i);
        x += (right-left)/(N-1);
    }
    return res;
}

std::vector<std::complex<double>> WaveFunction::p_on_grid(int N){
    std::vector<std::complex<double>> values = on_grid(N);
    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *N);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *N);
    fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    for(int i = 0; i < N; ++i){
        in[i][0] = values[i].real();
        in[i][1] = values[i].imag();
    }
    fftw_execute(p);
    for(int i = 0; i < (N+1)/2; i++){
        values[i+N/2] = interval/N/sqrt(2*M_PI)*std::complex<double>(out[i][0], out[i][1]);
    }
    for(int i = (N+1)/2; i < N; i++){
        values[i-(N+1)/2] = interval/N/sqrt(2*M_PI)*std::complex<double>(out[i][0], out[i][1]);
    }
    return values;
}


//gives the integral of the wavefunction squared
double WaveFunction::norm() const{
    return vector.dot(basis->overlap()*vector).real();
}

