#ifndef WIGNER
#define WIGNER

#include "basisSet.h"
#include "waveFunction.h"
#include "fftw3.h"
#include <fstream>
#include <iostream>

class Wigner{
public:
    LegendreScaled basis;
    WaveFunction wav;
    fftw_complex *in, *out;
    fftw_plan p;
    const std::function<std::complex<double>(double)>& _Psi;
    int _N;
    std::vector<std::complex<double>> values;


    Wigner(const std::function<std::complex<double>(double)>& Psi, int N = 100);
    ~Wigner();

    std::vector<std::complex<double>> val(int n);

    void writeFile(std::string filename);

};

#endif // WIGNER
