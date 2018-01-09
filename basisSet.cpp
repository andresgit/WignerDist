#include "basisSet.h"

#include "ap.h"                     //alglib to return quadrature
#include "orthopol.h"               //Our orthopogonal polynomial class in TOOLS
#include <iostream>
#include <vector>
#include <fstream>


using namespace std; // [STYLE] this is common practice, although sometimes discouraged

const double BasisSet::tolerance=1.e-13;
const double BasisSet::infty=DBL_MAX;

static const int nPlot=101;


vector<double> BasisSet::val(const double X) const{
    vector<double> _val,dum;
    valDer(X,_val,dum);
    return _val;
}

vector<double> BasisSet::der(const double X) const{
    vector<double> _der,dum;
    valDer(X,dum,_der);
    return _der;
}

void BasisSet::graph() const {

    // [STYLE] here we DO create a filename variable
    //         as it can be unambigously re-used later for info
    string filnam=_name+".dat";

    // open file
    std::ofstream datfil(filnam.c_str());

    // get range (with a reasonable guess if boundaries are too large)
    double plot0=max(lowerBoundary(),-3.);
    double plot1=min(upperBoundary(), 3.);

    vector<double> val,der;
    for(int k=0;k<nPlot;k++){
        double y=plot0+(plot1-plot0)*k/double(nPlot-1);
        // first column
        datfil<<y;
        // values and deriviatives of set
        valDer(y,val,der);
        for(int n=0;n<val.size();n++) {datfil<<" "<<val[n]<<" "<<der[n];}
        datfil<<endl;
    }
    // close file
    datfil.close();
    std::cout<<"\n--- function values and derivatives in "+filnam+" ---"<<std::endl;
}

Eigen::MatrixXd BasisSet::overlap() const{
    /* [STYLE] as a rule we use
     *         captital first letters for function arguments,
     *         lower case first letters for local variables
     * (or any other rule the is adhered to within one code project)
     */
    std::vector<double> val,der,points,weights;

    // by definition, quadrature rule should
    // return just the right number of points to be exact for overlap integrals
    quadrature(points,weights);

    Eigen::MatrixXd phi(size(),points.size());   //Phi(i,j)=Phi_i(x_j)
    for (int i=0; i<points.size(); i++) {
        valDer(points[i],val,der);
        for (int k=0; k<size(); k++)
            phi(k,i)=val[k];
    }

    // the beauty of Eigen:
    // we can create new matrices on the spot and without any extra computational cost
    return phi
            *Eigen::Map<Eigen::VectorXd>(weights.data(),weights.size(),1).asDiagonal()
            *phi.transpose();
    /* [STYLE] the code above may not be beautiful to look at but
     *         it is good style and good performance because
     * - we have avoided creating any auxiliary variables
     * - we do not need to copy data around
     * - we use the numbers in weights in place and interprete them
     *   as the values of a diagonal matrix
     *   (no useless storage of zeros)
     */
}

bool BasisSet::Test(){

    bool good=true;

    // run standard tests
    good=good and LegendreScaled(12,-2.,3.).test(true,true);

    return good;
}

bool LegendreScaled::test(bool Print, bool Plot) const {

    bool good=true;


    // verify for 11 points
    for(int k=0;k<11;k++){
        // boundaries can be infty, restrict range within [-5,5]
        double X=max(lowerBoundary(),-5.)+
                0.1*k*(min(upperBoundary(),5.)-max(lowerBoundary(),-5.));

        // direct and recursive evaluation
        vector<double>valRec,derRec,valDir;
        valDir=directVal(X);
        valDer(X,valRec,derRec);

        // [C++] unfortunately, std::min does not automaticall convert types
        //       we must make sure we have both "int"
        for(int n=0;n<min(size(),(int)valDir.size());n++){
            good=good and (abs(valDir[n]-valRec[n])<tolerance);
        }
    }

    if(Print){
        cout<<"Overlap matrix\n"<<overlap()<<endl;
        cout<<"\n";
        if(directVal(0).size()==0)cout<<"No ";
        else if(good)cout<<"OK ";
        else         cout<<"BAD";
        cout<<" compare with explict formulae for "<<min((int)directVal(0).size(),size())<<" functions \""<<name()<<"\""<<endl;
    }

    if(Plot)graph();
    return good;
}

LegendreScaled::LegendreScaled(int N, double X0, double X1):
    BasisSet("LegendreScaled",N),x0(X0),x1(X1){_quadOrder=N;}

void LegendreScaled::valDer(const double X, std::vector<double> &Val, std::vector<double> &Der) const{

    // [STYLE] a simple an error-save way of writing this (fast may be later)
    double y=2*(X-x0)/(x1-x0)-1.;
    OrthogonalLegendre().valDer(_nPhi,y,Val,Der);
    // remember that the derivative becomes scaled with the interval
    for(int k=0;k<Der.size();k++) {Der[k]/=(x1-x0)/2.;}
}


// [STYLE] clearly, we  only need to scale the quadrature rule from the Legendre polynomials
//         again, we start from "simple-and-save"
void LegendreScaled::quadrature(std::vector<double> &Points, std::vector<double> &Weights, int Kpoints) const{

    if(Kpoints<1)OrthogonalLegendre().quadrature(_quadOrder,Points,Weights);
    else         OrthogonalLegendre().quadrature(   Kpoints,Points,Weights);

    for(unsigned int k=0;k<Points.size();k++)
    {
        // shift and scale from [-1,1] to [x0,x1]
        Points[k]=x0+(Points[k]+1)*(x1-x0)/2.;
        Weights[k]*=(x1-x0)/2.;
    }

}
