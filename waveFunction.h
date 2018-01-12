#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <vector>
#include <complex>
#include <functional>

#include <Eigen/Dense>
#include "basisSet.h"

class WaveFunction{
    const BasisSet* basis;

public:
    Eigen::VectorXcd vector;
    double interval;

    WaveFunction(const BasisSet* Basis);
    WaveFunction(const BasisSet* Basis, const Eigen::VectorXcd& Vector);

    void set(const std::function<std::complex<double>(double)>& Psi);
    std::vector<std::complex<double>> on_grid(int N, double left=0, double right=0);

    double norm() const;

    WaveFunction& operator/=(const std::complex<double>& scalar){ vector/=scalar; return *this; }
    WaveFunction& operator*=(const std::complex<double>& scalar){ vector*=scalar; return *this; }
    WaveFunction& operator+=(const WaveFunction& other){ vector+=other.vector; return *this; }
    WaveFunction& operator-=(const WaveFunction& other){ vector-=other.vector; return *this; }
    
    WaveFunction operator/(const std::complex<double>& scalar) const { return WaveFunction{basis, vector/scalar}; }
    WaveFunction operator*(const std::complex<double>& scalar) const { return WaveFunction{basis, vector*scalar}; }
    WaveFunction operator+(const WaveFunction& other) const { return WaveFunction{basis, vector+other.vector}; }
    WaveFunction operator-(const WaveFunction& other) const { return WaveFunction{basis, vector-other.vector}; }

    friend WaveFunction operator*(const std::complex<double>& scalar, const WaveFunction& wf){
        return wf*scalar;
    }

};



#endif
