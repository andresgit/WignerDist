#ifndef WIGNER
#define WIGNER

#include "basisSet.h"
#include "waveFunction.h"
#include "fftw3.h"

class Wigner{

    LegendreScaled basis;
    WaveFunction wav;
    fftw_complex *in, *out;
    fftw_plan p;
    const std::function<std::complex<double>(double)>& _Psi;
    int _N;
    std::vector<std::complex<double>> values;

public:
    Wigner(const std::function<std::complex<double>(double)>& Psi, int N = 100);
    Wigner::~Wigner();

    std::vector<std::complex<double>> val(int n);

};

#endif // WIGNER
