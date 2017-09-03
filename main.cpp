#define _USE_MATH_DEFINES
#include "FunctionalUtilities.h"
#include "FangOost.h"
#include "ODESolver.h"
#include <iostream>
int main(){
    constexpr int discreteX=100; //discrete X
    constexpr int discreteU=64; 
    auto delta=.8;
    auto alpha=.1;
    auto sigma=.3;
    //auto S0=50.0;

    auto M=100.0;
    auto xMin=0.0;
    auto initLower=0.0;//0 probability at 0
    auto initHigher=1.0; //probabilty of 1 at M
    auto xMax=M;
    auto alphaFn=[&](const auto& x){
        return std::complex<double>(x*alpha, 0.0);
    };
    auto sigmaFn=[&](const auto& x){
        return std::complex<double>(x*sigma*sigma*x*.5, 0.0);
    };
    auto uFn=[&](const auto& u){
        return std::complex<double>(0.0, -u);
    };
    auto matrixOfXByU=futilities::for_each_parallel(0, discreteU, [&](const auto& index){
        return odesolver::solveODE<discreteX>(alphaFn, sigmaFn, [&](const auto& x){
            return uFn(fangoost::getU(fangoost::computeDU(xMin, xMax), index));
        }, xMin, xMax, initLower, initHigher);
    });

    auto vectorOfDistributions=futilities::for_each_parallel(0, discreteU, [&](const auto& index){
        return fangoost::computeInvDiscrete(discreteX, 0.0, 10.0, matrixOfXByU[index]);
    });
    for(auto& val: vectorOfDistributions[50]){
        std::cout<<val<<std::endl;
    }
}