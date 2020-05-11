#ifndef MONTECARLOENGINE_H
#define MONTECARLOENGINE_H

/**
    *A class to estimate the price of European options on a multivariate geometric Brownian motion using a Monte Carlo estimator.
*/

#include "GBM.h"

class MC_engine_GBM{
private:

    GBM model; /**<The geometric Brownian motion on which we want to compute the price of options.*/
    double T; /**<The maturity of the options.*/
    int nbSimuls; /**<The number of simulations for the Monte Carlo estimator.*/

public:
    /**
        *Constructor.
        *@param S The geometric Brownian motion on which we want to compute the price of options.
        *@param T The maturity of the options.
        *@param nbSimuls The number of simulations.
    */
    MC_engine_GBM(GBM S, double T, int nbSimuls);


    /**
        *Estimate the price of an European option using a simple Monte Carlo estimator.
        *The payoff can have as many parameters as one wants, but should be a function of S_T ONLY (and NOT of function of St
        *for t<T).
        *@param Payoff a pointer to a function returning the payoff. The function function should take S_T as first parameter.
        *The other parameters of the function should be of type double, and be passed using a pointer to an array containing all the arguments
        *@see Payoff.h
        *@param params_payoff a pointer to an array of doubles that contain all the parameters needed to compute the Payoff function.
        *@param weights a pointer to an array containing the weights of each stock in the final payoff.
    */
    double computeEuropeanPrice( double (*Payoff)(double spot, double* params), double* params_payoff, double* weights );

    /**
        *Compute the price of an European option using a Monte Carlo estimator with a control variate.
        *The payoff can have as many parameters as one wants, but should be a function of ST ONLY (and NOT of function of St
        *for t<T).
        *@see computeEuropeanPrice( double (*Payoff)(double spot, double* params), double* params_payoff, double* weights ).
        *@param Payoff a pointer to a function returning the payoff. The function function should take ST as first parameter
        *and as many parameters as one would like after ST.
        *@param params_payoff a pointer to an array of doubles that contain all the parameters needed to compute the Payoff function.
        *@param weights a pointer to an array containing the weights of each stock in the final payoff.
    */
    double computeEuropeanPrice_controlVariate(double (*Payoff)(double spot, double* params), double* params_payoff, double* weights);


    /**
        *This function computes and prints the estimate of the price of a call, a put, and call-put options, with confidences intervals.
        *Note: this function has been implemented only for the numerical experiments discussed in the report.
        *This will most likely not be useful to the examiner.
        *The results will be printed on the screen in a row, in this order:
        *nbSimuls, Call price, call lower bound, call upper bound, put price, put lb, put up, c-p, c-p lb, c-p ub.
    */
    void testCallPut(double* K, double* weights, double* res );

    /**
        *Same as testCallPut() but for the control variate estimator.
        *@see testCallPut.
    */
    void testCallPut_control(double* K, double* weights);

};

#endif // MONTECARLOENGINE_H
