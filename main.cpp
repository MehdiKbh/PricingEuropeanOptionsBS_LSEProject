#include <iostream>
#include <math.h>
#include <cstdlib>
#include <time.h>
#include <fstream>
using namespace std;

#include "matrix.h"
#include "MatrixOperations.h"
#include "Sampler.h"
#include "GBM.h"
#include "MonteCarloEngine.h"
#include "Payoffs.h"

int main(){

    //Parameters of the stocks
    int nbStocks = 3;
    double sigma[] = {.1, .2, .3};
    Matrix vols(sigma, nbStocks);
    double rho[] = {.2, .8, .5};
    Matrix correls(rho, nbStocks);
    double x0[] = {100, 110, 120};
    Matrix Sinit(x0, nbStocks);
    double r = 0.01;

    //Parameters for sample paths
    int nbSteps = 10; //Temporal steps for the sample paths

    //Parameters of the options
    double weights[] = {1./3., 1./3., 1./3.};
    Matrix w(weights, nbStocks);
    double T=1.; //Maturity
    int n_MC = 50000; //number of Monte Carlo simulations
    double K = 100.; //Strike of the option


    //PB1
    cout << " ############ PB1 ############\n";

    cout << nbStocks << " Stocks\n\n";
    cout << "Volatilities:\n";
    vols.T().print_matrix();
    cout << "Correlations:\n";
    correls.T().print_matrix();
    //Get the covariance Matrix from the values of sigmas and rho and print it
    Matrix covariance = BuildCovarianceMatrix(sigma, rho, nbStocks);
    cout << "Covariance Matrix: \n";
    covariance.print_matrix();
    //Compute the Cholesky Factor A
    Matrix A = covariance.CholeskyDecomp();
    cout << "Cholesky Factor of the covariance matrix (matrix A): \n";
    A.print_matrix();
    cout << "A is indeed lower triangular.\n\n";
    //Check if A is indeed the Cholesky Factor
    cout << "A*A_T=\n",
    (A*A.T()).print_matrix();
    cout << "We find exactly the covariance matrix \n";
    cout << endl;


    //PB2
    cout << "\n\n ############ PB2 ############\n";
     // Set the seed of the random number generator using time
    srand(time(0));

    //Generate samples paths
    GBM S(nbStocks, x0, r, covariance); //Define a object of type GBM (geometric brownian motion)
    cout << "Sample paths with " << nbSteps << " time steps\n";
    S.samplepath(T,nbSteps).print_matrix();
    cout << "Another sample path\n";
    S.samplepath(T, nbSteps).print_matrix();

    //Test the sample means and variance of the sample paths generator
    //This part is not of the most importance and can be commented
    int n=10000;
    int n_steps = 4;
    cout << "Test of the sample means for each time step, and the sample terminal variance, for " << n << " simulations and " << n_steps << " time steps:" << endl;
    S.samplepath_mean_var(T, n_steps, n);
    cout << endl;



    //PB3
    cout << "\n\n ############ PB3 ############\n";
    //Test of the Monte Carlo estimators
    MC_engine_GBM simulator(S, T, n_MC); //Declare a mc estimator
    double* pt_K = &K; //Declare a pointer to the strike K

    cout << "Estimation of European options prices with " << n_MC << " simulations\n" << endl;
    cout << "Call option strike " << K << " maturity " << T << " years with the classical MC estimator:\n";
    double c = simulator.computeEuropeanPrice(&callPayoff, pt_K, weights);
    cout << "Price: " << c << endl;
    cout << endl;
    cout << "With the Control Variate estimator:" << endl;
    double cc = simulator.computeEuropeanPrice_controlVariate(&callPayoff, pt_K, weights);
    cout << "Price: " << cc << endl;
    cout << endl;

    cout << "Put option strike " << K << " maturity " << T << " years with the classical MC estimator:\n";
    double p = simulator.computeEuropeanPrice(&putPayoff, pt_K, weights);
    cout << "Price: " << p << endl;
    cout << endl;
    cout << "With the Control variate estimator:" << endl;
    double pc = simulator.computeEuropeanPrice_controlVariate(&putPayoff, pt_K, weights);
    cout << "Price: " << pc << endl;
    cout << endl;

    cout << "C - P with the classical MC estimator: " << c - p << endl;
    cout << "C - P with the Control Variate estimator: " << cc-pc << endl;
    cout << "Theoretical (correct) value (S-PV(K)): " << (w.T()*Sinit).get_element(0,0) - exp(-r*T)*K << endl;

    return 0;

}
