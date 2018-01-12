#include <iostream>
#include "wigner.h"
#include <math.h>

using namespace std;

std::complex<double> f1(double x){
    return std::complex<double>(exp(-x*x/2), 0);
}

//void fftwTest(){
//    std::vector<double> reals = {0,0,0,0,0,0,6.1349e-5,0,0,0,0,0};
////    std::vector<std::complex<double>> vec = {std::complex<double>(0,0), std::complex<double>(1,0), std::complex<double>(0,0)};

//    int N = reals.size();
//    fftw_complex *in, *out;
//    fftw_plan p;

//    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
//    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
//    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

//    for(int i = 0; i < N; i++){
//        in[i][0] = reals[i];
//        in[i][1] = 0;
//    }

//    for(int i = 0; i < N; i++){
//        std::cout << "(" << in[i][0] << "," << in[i][1] << ")";
//    }
//    std::cout << std::endl;

//    fftw_execute(p);

//    for(int i = 0; i < N; i++){
//        std::cout << std::complex<double>(out[i][0], out[i][1]);
//    }
//    std::cout << std::endl;
//}


int main()
{
//    fftwTest();
    double L = 5;
    LegendreScaled basis(20,-L,L);
    WaveFunction wav(&basis);
    wav.set(f1);
    Wigner wig(&wav, 100);

//    std::cout << "Val at 0:" << std::endl;
//    wig.val(0);
//    std::cout << "Val at 3:" << std::endl;
//    wig.val(3);
//    std::cout << std::endl;

//    std::vector<std::complex<double>> vals = wig.val(2);
//    std::cout << vals.size() << std::endl;

    wig.writeFile("test1.txt");

    std::cout << "Test value:" << std::endl;
    wig.val(7);

//    Eigen::VectorXcd vec = wav.vector;
//    for(int i = 0; i < vec.size(); ++i){
//        cout << vec(i) << " ";
//    }

//    std::cout << "Values: " << std::endl;
//    std::vector<std::complex<double>> vec = wig.values;
//    for(int i = 0; i < vec.size(); ++i){
//        cout << vec[i] << " ";
//    }

//    std::cout << std::endl;

    return 0;
}

