#include <iostream>
#include "wigner.h"

using namespace std;

std::complex<double> f1(double x){
    return std::complex<double>(x*x, x*x*x);
}

int main()
{
    Wigner wig(&f1, 11);

    wig.writeFile("test1.txt");

//    Eigen::VectorXcd vec = wig.wav.vector;
//    for(int i = 0; i < vec.size(); ++i){
//        cout << vec(i) << " ";
//    }

//    std::cout << "Values: " << std::endl;
//    std::vector<std::complex<double>> vec = wig.values;
//    for(int i = 0; i < vec.size(); ++i){
//        cout << vec[i] << " ";
//    }

    cout << wig.values.size() << " " << wig._N << endl;
    return 0;
}

