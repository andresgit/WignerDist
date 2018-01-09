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

std::vector<std::complex<double>> WaveFunction::on_grid(double left, double right, double spacing){
    std::vector<std::complex<double>> res;
    for(double x=left; x<=right; x+=spacing){
        std::vector<double> val = basis->val(x);
        res.push_back(0);
        for(int i=0; i<val.size(); i++) res.back()+=val[i]*vector(i);
    }
    return res;
}


//gives the integral of the wavefunction squared
double WaveFunction::norm() const{
    return vector.dot(basis->overlap()*vector).real();
}

