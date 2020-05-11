#ifndef PAYOFFS_H
#define PAYOFFS_H

/**
    *File to declare payoffs that can be used with the class MC_engine_GBM to estimate the time-0 price.
    *The examiner can implement its own payoff and test them directly calling the function computeEuropeanPrice of the class MC_engine_GBM.
*/


double callPayoff(double spot, double *K); /**<Payoff of a Call option*/
double putPayoff(double spot, double* K); /**<Payoff of a Put option*/

#endif // PAYOFF_H
