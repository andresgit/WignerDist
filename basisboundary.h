#ifndef BASISBOUNDARY_H
#define BASISBOUNDARY_H

#include <basisSet.h>

class BasisBoundary: public BasisSet
{
private:
    bool _Zero;
    double x0,x1;
    int _nPhi;
    LegendreScaled legScaled;
    int quadOrder;
public:
    BasisBoundary(int N, const double X0, const double X1, bool Zero=false);
    void valDer(const double X, std::vector<double>& Val, std::vector<double>& Der) const;
    void quadrature(std::vector<double>& Points, std::vector<double>& Weights, int Kpoints=-1) const;
    bool test(bool Print=false, bool Plot=false) const;
    Eigen::MatrixXd matrix(const std::string & Operator, const BasisSet & LeftBas, const BasisSet & RightBas) const;

    inline double lowerBoundary() const {return x0;}
    inline double upperBoundary() const {return x1;}
};

#endif // BASISBOUNDARY_H
