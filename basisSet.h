#ifndef BASISSET_H
#define BASISSET_H

#include <vector>       // the standard library vector
#include <string>       // class std::string for manipulating strings
#include <cfloat>       // defines floating point constants valid for the given C++ implementation
#include <iostream> // standard I/O routines, prefix with "std::"
#include "abort.h"
#include "Eigen/Dense"



class BasisSet{    
protected:
    // here not private, but protected, as we want them available in all derived classes
    static const double tolerance;  // tolerance for tests
    static const double infty;      // a large number indicating "infinity"

    int _nPhi;            // number of functions in basis set
    int _quadOrder;       // number of points needed for exact quadrature for all overlaps
    std::string _name;    // name for debugging reasons

    /// evaluate using explicit formula for the first few functions (for tests)
    virtual std::vector<double> directVal(double X) const
    {
        std::cout<<" --- WARNING: no direct evaluation defined, cannot check "
                   +name()+" ---"<<std::endl;
        return std::vector<double>(0);
    }
    /* [STYLE] comment on how to put comments:
     *
     * - for multi-line comments use the present format
     *   it can be collapsed easily in QTcreator
     *
     * - for single-line comments, the format
     *   // this is a single line comment
     *   is fine.
     *
     * - do NOT use the beautiful format of
     *   //========to access private data===========//
     *   general rule: no redundand info (such as ==== )
     *
     * - for any info that describes the purpose of functions
     *   or the use of arguments (i.e. what a "user" would be interested in)
     *   use the format suitable for the documentation utility "Doxygen":
     */
    /// this is a single line Doxygen comment
    /** and this
     * is a multi-line Doxygen comment
     */
    /* Doxygen comments will show in automatically generated code documentation,
     *   extremely useful for the dummy user, but, by own experience,
     *   also for the smart developper
     */

public:
    // virtual destructor is needed for ensuring
    // correct destruction of derived class if polymorphism is used
    virtual ~BasisSet(){}

    BasisSet(std::string Name, int N):_name(Name), _nPhi(N),_quadOrder(N){}

    /* [STYLE] as a rule, no data of a class should be public
     *         provide access through functions as below
     *  reasons:
     *         - no accidental change of internal data
     *         - controlled interface with constant meaning
     *           even if the form of information storage changes,
     *           e.g. a size is evaluated each time, rather than just
     *           read from some stored value
     */
    int size() const {return _nPhi;}
    int quadOrder() const {return _quadOrder;}
    std::string name() const {return _name;}

    // we found that this member function is not very useful and will be eleminated
    bool compareLiterature(bool verbose);

    /// standard sanity test optionally print and generate plot files
    virtual bool test(bool Print, bool Plot) const =0;

    std::vector<double> val(const double X) const;  ///< values for all basis functions
    std::vector<double> der(const double X) const;  ///< derivatives for all basis functions

    /// values and derivatives for all basis functions
    virtual void valDer(const double X, std::vector<double>& Val, std::vector<double>& Der) const =0;
    // pure virtual as BasisSet has no way of evaluating it as it stands

    /** Kpoints - points quadrature Points[k], Weights[k], k<Kpoints;
     * Kpoints=-1 returns rule such that
     * int dx Phi_m(x) Phi_n(x) is exact for all functions in basis set
     */
    virtual void quadrature(std::vector<double>& Points, std::vector<double>& Weights, int Kpoints=-1) const =0;
    // pure virtual as there is no general recipe for quadratures

    /// lower boundary of function argument
    virtual double lowerBoundary() const {
        ABORT("no lowerBoundary() specified for "+name());
    }
    /// upper boundary of function argument
    virtual double upperBoundary() const {
        ABORT("no upperBoundary() specified for "+name());
    }
    /* [STYLE] this construction does not force you to define a boundary,
     *         but when you try to use it, the program tells you it is missing
     *  this is not terribly good style, but a legitimate way of lazy programming
     */

    /// run a set of tests
    static bool Test();

    /// file with columns x Phi0[x] Phi0'[x] Phi1[x] Phi1'[x] ... Phi'(N-1)[x]
    void graph() const;

    /// returns overlap matrix <Phi_i|Phi_j>
    Eigen::MatrixXd overlap() const;

};

/// shifted and scale Legendre polynomials
class LegendreScaled:public BasisSet{
private:
    // defined on [x0,x1]
    double x0,x1;

    /// directly evaluate the first 4 functions
    std::vector<double> directVal(double X) const {
        std::vector<double> val;
        double y=2*(X-x0)/(x1-x0)-1.;
        val.push_back(1.);
        val.push_back(y);
        val.push_back( 1.5*y*y-0.5);
        val.push_back((2.5*y*y-1.5)*y);
        return val;
    }
public:
    /// create polynomials up to degree N-1 in [X0,X1]
    LegendreScaled(int N, double X0, double X1);

    void valDer(const double X, std::vector<double>& Val, std::vector<double>& Der) const;
    void quadrature(std::vector<double>& Points, std::vector<double>& Weights, int Kpoints=-1) const;

    inline double lowerBoundary() const {return x0;}
    inline double upperBoundary() const {return x1;}
    // [STYLE] "inline" encourages the compiler to
    // replace the function calls with direct insertion of x0 or x1

    bool test(bool Print=false, bool Plot=false) const;

};

#endif // BASISSET_H
