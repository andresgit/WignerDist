#include <iostream>
#include "wigner.h"
#include <math.h>

using namespace std;

std::complex<double> f1(double x){
    return std::complex<double>(exp(-x*x/2), 0);
}


int main()
{
    double L = 5;
    LegendreScaled basis(20,-L,L);
    WaveFunction wav(&basis);
    wav.set(f1);
    Wigner wig(&wav, int(50));

    wig.writeFile("test1.txt");

    return 0;
}

