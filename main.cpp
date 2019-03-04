#include <iostream>
#include "wigner.h"
#include <math.h>
#include <basisSet.h>

using namespace std;

std::complex<double> f1(double x){
    return std::complex<double>(exp(-x*x/2), 0);
}


int main(int argc, char *argv[])
{
    Wigner::test();
//    Wigner::oscEigenStates();
//    Wigner::oscCoherent();
//    Wigner::checkProb();

//    double x0, sig;
//    x0 = 1;
//    sig = 2;
//    if(argc == 3){
//        x0 = atof(argv[1]);
//        sig = atof(argv[2]);
//    }

//    Wigner::timeEvoTest(0, 0);

    return 0;
}

