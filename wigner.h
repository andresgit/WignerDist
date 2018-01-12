#ifndef WIGNER
#define WIGNER

#include "basisSet.h"
#include "waveFunction.h"
#include <fftw3.h>
#include <fstream>
#include <iostream>
#include <math.h>

class Wigner{
public:
    fftw_complex *in, *out;
    fftw_plan p;
    int N;
    std::vector<std::complex<double>> values;
    double interval;


    Wigner(WaveFunction *wav, int N = 100);
    ~Wigner();

    std::vector<std::complex<double>> val(int n);

    void writeFile(std::string filename);

};

#endif // WIGNER
