#ifndef WIGNER
#define WIGNER

#include "basisSet.h"
#include "waveFunction.h"
#include <fftw3.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <abort.h>
#include <iomanip>

class Wigner{
public:
    fftw_complex *in, *out;
    fftw_plan p;
    WaveFunction *wav;
    int N;
    std::vector<std::complex<double>> values;
    double interval;
    double imagMax;


    Wigner(WaveFunction *wav, int N = 100);
    ~Wigner();

    std::vector<double> val(int n);
    void updateGrid();

    void writeFile(std::string filename);

    template<typename T>
    void writeFileVector(std::string filename, std::vector<T> vector);

    static void checkProb();

    static void writeFileSum(std::string filename, std::vector<Wigner*> wigs);

    static void test();

    static void oscEigenStates();

    static void oscCoherent();

};

#endif // WIGNER
