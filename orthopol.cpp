#include "orthopol.h"
#include "orthogonalChebyshev.h"
#include <iostream>          // standard I/O routines, prefix with "std::"
#include <vector>
#include "ap.h"
#include "integration.h"
#include "orthogonalLaguerre.h"

using namespace std; // [STYLE] this is common practice, although sometimes discouraged

const double OrthogonalPolynomial::tolerance=1.e-13;
const double OrthogonalPolynomial::infty=DBL_MAX;

void OrthogonalPolynomial::valDer(int N, const double X, std::vector<double> &Val, std::vector<double> &Der) const {
    // reserve storage
    Val.resize(N);
    Der.resize(N);

    // degree = 0
    if(N>0){
        Val[0]=val0();
        Der[0]=0;
    }
    // degree = 1
    if (N>1) {
        Val[1]=a(1)* X*Val[0]        + b(1);
        Der[1]=a(1)*(X*Der[0]+Val[0]);
    }
    // degree > 1: use two-point recurrence
    for (int n=2; n<N; n++){
        Val[n] = (a(n)*X + b(n)) * Val[n-1]                 + c(n)*Val[n-2];
        Der[n] = (a(n)*X + b(n)) * Der[n-1] + a(n)*Val[n-1] + c(n)*Der[n-2];
    }
}

/// values at X up to degree N-1
vector<double> OrthogonalPolynomial::val(int N, const double X) const{
    vector<double> _val,dum;
    valDer(N,X,_val,dum);
    return _val;
}

/// derivatives at X up to degree N-1
vector<double> OrthogonalPolynomial::der(int N, const double X) const{
    vector<double> _der,dum;
    valDer(N,X,dum,_der);
    return _der;
}

bool OrthogonalPolynomial::Test(){ 
    bool good=true;
    good=good and OrthogonalLegendre().test(40,true,true);
    good=good and OrthogonalChebyshev().test(40,true,true);
    good=good and OrthogonalLaguerre().test(40,true,true);

    return good;
}

bool OrthogonalPolynomial::compareLiterature(bool Verbose){

    bool good=true;
    // verify for 11 points
    for(int k=0;k<11;k++){
        // boundaries can be infty, restrict range within [-5,5]
        double X=max(lowerBoundary(),-5.)+
                0.1*k+(min(upperBoundary(),5.)-max(lowerBoundary(),-5.));

        // direct and recursive evaluation
        vector<double>valRec,derRec,valDir;
        valDir=directVal(X);
        valDer(valDir.size(),X,valRec,derRec);

        for(int n=0;n<valDir.size();n++){
            good=good and (abs(valDir[n]-valRec[n])<tolerance);
        }
    }
    if(Verbose){
        if(directVal(0).size()==0)
            cout<<"no comparison literature: directVal() not defined for "<<name()<<endl;
        else if(good)
            cout<<"OK compares literature: "<<name()<<" up to degree "<<directVal(0).size()-1<<endl;
        else
            cout<<"BAD comparison to literature failed: "<<name()<<" up to degree "<<directVal(0).size()-1<<endl;
    }
    return good;
}
bool OrthogonalPolynomial::test(int Order, bool Print, bool Plot){

    // get quadrature points and weights
    vector<double> pts,wgs;
    quadrature(Order,pts,wgs);

    // integrals int dx w(x) Qn(x) Qm(x) for m<=n
    vector<double> ints((Order*(Order+1))/2,0.); // reserve storage and initialize =0

    // sum over quadrature points
    for(int k=0;k<pts.size();k++)
    {
        // store values of all polynomials for argument x=pts[k]
        vector<double> vals(val(Order,pts[k]));

        // add into integrals
        for(int n=0,mn=0;n<Order;n++){
            double valnw=vals[n]*wgs[k]; // [PERFORMANCE] pull multiplication out of loop
            for(int m=0;m<=n;m++,mn++)
                ints[mn]+=valnw*vals[m];
            /* [PERFORMANCE] if you have more than one multiplication in innermost loop
             *               you probably overlooked something (can cost a factor ~2 in time)
             */
        }
    }

    // check orthogonality
    int nCheck=0;
    double epsMin=DBL_MAX,epsMax=0;
    int mMin,nMin,mMax,nMax;
    for(int n=0,mn=0;n<Order;n++){
        for(int m=0;m<=n;m++,mn++){
            if(n==m)continue; // do not check diagonal
            nCheck++;
            double eps=ints[mn]/sqrt(ints[(n*(n+1))/2+n]*ints[(m*(m+1))/2+m]);
            if(eps>tolerance){
                if(eps<epsMin){epsMin=eps; mMin=m; nMin=n;}
                if(eps>epsMax){epsMax=eps; mMax=m; nMax=n;}
            }
        }
    }

    // report if desired
    if(Print){
        cout<<endl;
        if(epsMax>tolerance)
            cout<<"FAILED "<<name()<<" at order="<<Order<<": epsMin="<<epsMin<<"("<<mMin<<","<<nMin<<"), epsMax="<<epsMax<<"("<<mMax<<","<<nMax<<")"<<endl;
        else
            cout<<"OK "<<name()<<" orthogonal at order="<<Order<<", "<<nCheck<<" integrals checked"<<endl;
    }

    // here we should open a file and write base-points and values to it
    if(Plot)
        cout<<"Plot==true not implemented yet"<<endl;

    // return true if error below tolerance
    return epsMax<tolerance and compareLiterature(Print);
    ;
}


void OrthogonalPolynomial::quadrature(int N, vector<double>& Points, vector<double>& Weights) const {
    alglib::real_1d_array alpha;
    alglib::real_1d_array beta;
    alglib::ae_int_t info;
    alglib::real_1d_array xq;
    alglib::real_1d_array wq;

    // convert to general recurrence coefficients as used in alglib
    alpha.setlength(N);
    beta.setlength(N);
    if(N>0){
        alpha[0]=-b(1)/a(1);
        beta[0] =0.;
    }
    for(unsigned int n=1;n<N;n++){
        alpha[n]=-b(n+1)/a(n+1);
        beta[n] =-c(n+1)/(a(n+1)*a(n));
    }

    // generate the rules
    alglib::gqgeneraterec(alpha,beta,normsq0(),N,info,xq,wq);
    switch(info){
    case  1: break; // success
    case -3: cout<<"QuadratureRule: internal eigenproblem solver hasn't converged"<<endl;exit(1);
    case -2: cout<<"QuadratureRule: Beta[i]<=0"<<endl;exit(1);
    case -1: cout<<"QuadratureRule: incorrect N was passed"<<endl;exit(1);
    default: cout<<"QuadratureRule: error value info = "<<info<<endl;exit(1);
    }

    // transfer to std::vector
    Points.resize(N);
    Weights.resize(N);
    for (int i=0; i<Weights.size(); i++) {
        Points[i]=xq[i];
        Weights[i]=wq[i];
    };
}
