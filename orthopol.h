#ifndef ORTHOPOL_H // "include guard": suppresses multiple includes the header file
#define ORTHOPOL_H

#include <vector> // the standard library vector and string
#include <string> // class std::string for manipulating strings
#include <cfloat> // defines floating point constants valid for the given C++ implementation
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <iostream> // standard I/O routines, prefix with "std::"
#include "abort.h"


class OrthogonalPolynomial{

    static const double tolerance;  // tolerance for tests

protected:
    static const double infty;      //infinty from cfloat

    const std::string _name;
    /// coefficients for recurrence relations
    virtual double a(int N) const =0;
    virtual double b(int N) const =0;
    virtual double c(int N) const =0;
    virtual double normsq0() const =0; // this is required for setting up the quadrature rules
    virtual double val0() const {return 1.;}

    /// directly evaluate the first few polynomials
    virtual std::vector<double> directVal(double X) const
    {
        std::cout<<" --- WARNING: no direct evaluation defined, cannot check "
                   +name()+" ---"<<std::endl;
        return std::vector<double>(0);
    }
    // return true, if first few values agree with literature definition within tolerance
    bool compareLiterature(bool Verbose=false);

public:

    OrthogonalPolynomial(std::string Name):_name(Name){}
    virtual ~OrthogonalPolynomial(){}

    std::string name() const {return _name;}

    /// return values/derivatives up to degree N-1 (order N)
    std::vector<double> val(int N, const double X) const;
    std::vector<double> der(int N, const double X) const;
    virtual void valDer(int N, const double X, std::vector<double>& Val, std::vector<double>& Der) const;

    /// N-point gauss quadrature rule
    void quadrature(int N, std::vector<double>& Points, std::vector<double>& Weights) const;

    /// lower boundary of integration range for polynomial
    virtual double lowerBoundary() const {std::cout<<"no lowerBoundary() specified for "+name()<<std::endl; exit(1);}

    /// upper boundary of integration range for polynomial
    virtual double upperBoundary() const {std::cout<<"no upperBoundary() specified for "+name()<<std::endl; exit(1);}

    /// true if orthogonal up to Order within tolerance, optionally print and generate plot files
    bool test(int Order, bool Print=false, bool Plot=false);

    /// run a range of tests (returns true if all tests pass)
    static bool Test();
    /* another meaning of "static":
     *  - Test() must be called *without* an object
     *  - it is public OrthogonalPolynomial::Test()
     *    will work wherever "orthogonalPolynomial.h" is included
     */

};
// ==== end of abstract base class header, below there are specifications ====

// basic example - Legendre polynomials
class OrthogonalLegendre:public OrthogonalPolynomial{

    inline double a(int N) const {return (2.*double(N)-1.)/double(N);}
    inline double b(int N) const {return 0.;}
    inline double c(int N) const {return   (-double(N)+1.)/double(N);}
    inline double normsq0() const {return 2.;} // normalization

    std::vector<double> directVal(double X) const {
        std::vector<double> val;
        val.push_back(1.);
        val.push_back(X);
        val.push_back( 1.5*X*X-0.5);
        val.push_back((2.5*X*X-1.5)*X);
        return val;
    }

public:
    OrthogonalLegendre():OrthogonalPolynomial("legendre"){}
    inline double lowerBoundary() const {return -1.;}
    inline double upperBoundary() const {return  1.;}
};


///
#endif
