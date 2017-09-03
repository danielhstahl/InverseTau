#define _USE_MATH_DEFINES
#include "FunctionalUtilities.h"
#include "FangOost.h"
#include "ODESolver.h"
#include <iostream>

std::vector<std::vector<double> > transpose(const std::vector<std::vector<std::complex<double> > >& data, double cp) {
    std::vector<std::vector<double> > result(data[0].size(),
                                          std::vector<double>(data.size()));
    for (std::vector<double>::size_type i = 0; i < data[0].size(); i++) 
        for (std::vector<double>::size_type j = 0; j < data.size(); j++) {
            result[i][j] = data[j][i].real()*cp;
        }
    return result;
}
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
    auto cp=fangoost::computeCP(du);
    auto matrixOfXByU=futilities::for_each_parallel(0, discreteU, [&](const auto& index){
        auto fnConstant=uFn(fangoost::getU(du, index));
        return odesolver::solveODE_diff(sigmaFn, alphaFn, [&](const auto& x){
            return fnConstant;
        }, initLower, initHigher, xMin, xMax, discreteX);
    });
    auto dT=fangoost::computeDX(discreteTau, tMin, tMax);
    auto transposeVector=transpose(matrixOfXByU, cp);
    auto vectorOfDistributions=futilities::for_each_parallel(0, discreteX, [&](const auto& index){
        return fangoost::computeInvDiscrete(discreteTau, tMin, tMax, transposeVector[index]);
    });
    
    std::cout<<"x, density"<<std::endl;
    for(int i=0;i<discreteTau;++i){
        std::cout<<fangoost::getX(tMin, dT, i)<<", "<<vectorOfDistributions[60][i]<<std::endl;
    }
}