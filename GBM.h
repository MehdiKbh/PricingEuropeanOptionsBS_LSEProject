#ifndef GBM_H
#define GBM_H

/**
    *Class to generate sample path from a multivariate geometric Brownian Motion.
    *The class relies on Sampler.h for generating random numbers from a multivariate normal distribution.
    *This class is inspired by the programming session 3.
*/

#include "matrix.h"
#include "MatrixOperations.h"
#include "Sampler.h"

const double X0[]={100.};

class GBM{

public:

    /**
        *Constructor.
        *@param n the number of stocks/Brownian Motions.
        *@param x0 a pointer to an array of the initial stock prices (size n).
        *@param R the risk-free instantaneous interest rate.
        *@param VarCovar the Variance/Covariance Matrix of the log returns.
        *@see Sampler.
    */
    GBM(int n=1, const double* x0=X0, double R=0, Matrix VarCovar=UNIT);

    void operator=(const GBM& model); /**<Overloading of the assignment operator.*/

    double* get_drift() const; /**<Get the drift mu-0.5*sigma^2 of the log returns.*/
    Matrix get_S0() const; /**<Get the initial stocks prices.*/
    int get_N() const; /**<Get the number of stocks.*/
    Matrix get_var() const; /**<Get the variance matrix of the log returns.*/
    double get_r() const; /**<Get the risk-free rate.*/

    //Simulations

    /**
        *Sample one path of the GBM.
        *@param T the maturity (in years).
        *@param nbSteps the number of temporal steps between 0 and T.
        *@return a matrix of shape (N,nbSteps+1), with one sample path (including S0) for each stock.
    */
    Matrix samplepath(double T, int nbSteps);

    /**
        *Compute the sample means for each time step and the sample variance for the terminal time step.
        *Note that this has been implemented only to test empirically the sample paths generator.
        *@see Matrix samplepath(double T, int nbSteps).
        *@param T Maturity.
        *@param nbSteps the number of time steps between 0 and T.
        *@param nbSimuls the number of simulations to compute the sample means and variance.
    */
    Matrix samplepath_mean_var(double T, int nbSteps, int nbSimuls);



private:
    int N; /**<The number of stocks.*/
    Matrix S0; /**<The column vector of initial stocks prices (size nx1).*/
    double r; /**<The risk-free rate*/
    Matrix Var; /**<The variance matrix of the log returns.*/

    /**
        *Generate one step ahead of the GBM.
        *Starting from the stocks prices currentX, this generates one step ahead, of temporal size h, of the GBM.
        *@param currentX the initial stocks prices.
        *@param h the temporal step.
        *@return matrix of shape (N,1) of one-step ahead simulation for each stock.
    */
    Matrix step(Matrix currentX, double h);

};
#endif // GBM_H

