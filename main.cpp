#include <iostream>
#include "wigner.h"
#include <math.h>

using namespace std;

std::complex<double> f1(double x){
    return std::complex<double>(exp(-x*x/2), 0);
}


int main()
{
//    Wigner::test();
//    Wigner::oscEigenStates();
//    Wigner::oscCoherent();
    Wigner::checkProb();

    return 0;
}

