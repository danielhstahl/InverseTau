#define _USE_MATH_DEFINES
#include "FunctionalUtilities.h"
#include "FangOost.h"
#include "ODESolver.h"
#include <iostream>
#include <chrono>
int main(){
    constexpr int discreteX=298; //discrete X
    constexpr int discreteTau=1024; //discrete t
    constexpr int discreteU=512; 
    auto delta=.8;
    auto alpha=.1;
    auto sigma=.3;
    //auto S0=50.0;

    auto M=5.0; //hitting boundary
    auto xMin=0.0;
    auto initLower=0.0;//0 probability at 0
    auto initHigher=1.0; //probabilty of 1 at M
    auto xMax=M;
    auto tMin=0.0;
    auto tMax=176.7939;

    auto alphaFn=[&](const auto& x){
        return alpha*x;
    };
    auto sigmaFn=[&](const auto& x){
        return pow(x, 2.0*delta)*sigma*sigma*.5;
    };
    auto uFn=[&](const auto& u){
        return std::complex<double>(0.0, -u);
    };
    auto du=fangoost::computeDU(tMin, tMax);
    auto started = std::chrono::high_resolution_clock::now();
    auto matrixOfXByU=futilities::for_each_parallel(0, discreteU, [&](const auto& index){
        auto fnConstant=uFn(fangoost::getU(du, index));
        return odesolver::solveODE_diff(sigmaFn, alphaFn, [&](const auto& x){
            return fnConstant;
        }, initLower, initHigher, xMin, xMax, discreteX);
    });
    auto dT=fangoost::computeDX(discreteTau, tMin, tMax);
    auto vectorOfDistributions=futilities::for_each_parallel(0, discreteX, [&](const auto& index){
        return fangoost::computeInv(discreteTau, discreteU, tMin, tMax, [&](const auto& u){
            return matrixOfXByU[round(u.imag()/du)][index];
        });
    });
    auto done = std::chrono::high_resolution_clock::now();
    std::cout << "Total time: "<<std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count()<<std::endl;
    std::cout<<"x, density"<<std::endl;
    for(int i=0;i<discreteTau;++i){
        std::cout<<fangoost::getX(tMin, dT, i)<<", "<<vectorOfDistributions[60][i]<<std::endl;
    }
}