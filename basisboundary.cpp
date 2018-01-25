#include "basisboundary.h"

BasisBoundary::BasisBoundary(int N, const double X0, const double X1, bool Zero) :  BasisSet("Boundary", N), _Zero(Zero), x0(X0), x1(X1), _nPhi(N), legScaled(Zero ? N+2:N, X0, X1), quadOrder(Zero ? N+2:N)
{
}

void BasisBoundary::valDer(const double X, std::vector<double>& Val, std::vector<double>& Der) const{
    legScaled.valDer(X, Val, Der);
    if(_Zero){
        double V0 = Val[0];
        double V1 = Val[1];
        double D0 = Der[0];
        double D1 = Der[1];

        for(int i = 0; i < _nPhi; i+=2){
            Val[i] = Val[i+2]-V0;
            Val[i+1] = Val[i+3]-V1;
            Der[i] = Der[i+2]-D0;
            Der[i+1] = Der[i+3]-D1;
        }
        if(_nPhi%2 == 0){
            Val[_nPhi-1] = Val[_nPhi+1]- V1;
            Der[_nPhi-1] = Der[_nPhi+1] -D1;
        }
        Val.pop_back();
        Val.pop_back();
        Der.pop_back();
        Der.pop_back();
    }
}

Eigen::MatrixXd BasisBoundary::matrix(const std::string & Operator, const BasisSet & LeftBas, const BasisSet & RightBas) const{
    std::vector<double> Points;
    std::vector<double> Weights;

    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(_nPhi, _nPhi);
    if(Operator == "-(d/dx)^2"){
        quadrature(Points, Weights, quadOrder-1);
        for(int k = 0; k < Points.size(); ++k){
            std::vector<double> der1 = LeftBas.der(Points[k]);
            std::vector<double> der2 = RightBas.der(Points[k]);
            for(int m = 0; m < _nPhi; ++m){
                for(int n = 0; n < _nPhi; ++n){
                    res(m,n) += Weights[k]*der1[m]*der2[n];
//                    res(n,m) = res(m,n);
                }
            }
        }
    }
    else if(Operator == "x^2"){
        quadrature(Points, Weights, quadOrder+1);
        for(int k = 0; k < Points.size(); ++k){
            std::vector<double> val1 = LeftBas.val(Points[k]);
            std::vector<double> val2 = RightBas.val(Points[k]);
//            for(auto value: val1) std::cout << value << " ";
//            std::cout << std::endl;
            for(int m = 0; m < _nPhi; ++m){
                for(int n = 0; n < _nPhi; ++n){
                    res(m,n) += Weights[k]*val1[m]*val2[n]*pow(Points[k],2);
//                    res(n,m) = res(m,n);
                }
            }
        }
    }
    else if(Operator == "x^4"){
        quadrature(Points, Weights, quadOrder+2);
        for(int k = 0; k < Points.size(); ++k){
            std::vector<double> val1 = LeftBas.val(Points[k]);
            std::vector<double> val2 = RightBas.val(Points[k]);
//            for(auto value: val1) std::cout << value << " ";
//            std::cout << std::endl;
            for(int m = 0; m < _nPhi; ++m){
                for(int n = 0; n < _nPhi; ++n){
                    res(m,n) += Weights[k]*val1[m]*val2[n]*pow(Points[k],4);
//                    res(n,m) = res(m,n);
                }
            }
        }
    }
    return res;
}

void BasisBoundary::quadrature(std::vector<double>& Points, std::vector<double>& Weights, int Kpoints) const{
    legScaled.quadrature(Points, Weights, Kpoints);
}
bool BasisBoundary::test(bool Print, bool Plot) const{ //false, false by default

}
