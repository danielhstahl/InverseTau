## Inverse tau

This is the c++ implementation of the paper on generating distributions for first hitting times from GBMs.  See [Tau Distribution](https://github.com/phillyfan1138/TauDistribution).

However, it doesn't work yet.  Something is wrong when it inverts the matrix using Thomas' algorithm (see ODESolver).  I've checked that the R code and C++ match perfectly going into an inversion algorithm, but not coming out.