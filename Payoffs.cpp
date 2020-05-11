#include "Payoffs.h"

//Max
template <class T>
T myMax(T a, T b){return (a<b)?b:a;}

double putPayoff(double spot, double *K){return myMax(*K - spot, 0.);}

double callPayoff(double spot, double *K){return myMax(spot - *K, 0.);}
