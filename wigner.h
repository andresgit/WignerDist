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
#include "basisboundary.h"
#include <sstream>

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
    void newWaveFunc(WaveFunction *_wav);

    void writeFile(std::string filename);

    template<typename T>
    void writeFileVector(std::string filename, std::vector<T> vector);

    static void checkProb();

    static void writeFileSum(std::string filename, std::vector<Wigner*> wigs);

    static void test();

    static void oscEigenStates();

    static void oscCoherent();

    static void timeEvo(std::complex<double> (*wavFunc)(double), double eps, std::string fileNameStart, double T, int N);

    static void timeEvoTest(double x0, double sig);
};

#endif // WIGNER
