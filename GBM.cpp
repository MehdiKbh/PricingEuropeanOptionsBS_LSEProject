#include "GBM.h"
#include "matrix.h"
#include "MatrixOperations.h"

#include <iostream>
#include <cmath>
using namespace std;


//Constructor
GBM::GBM(int n, const double* x0, double R, Matrix VarCovar){
    N=n;
    Var = VarCovar;
    r = R;

    Matrix tmp_x0(x0, n);
    S0 = tmp_x0;
}


//Getters

//Drift mu-0.5*sigma^2
double* GBM::get_drift()const {
    double* drift = new double[this->N];
    for(int i=0; i<this->N; i++){
        drift[i] = r - 0.5 * Var.get_element(i,i);
    }
    return drift;}

Matrix GBM::get_S0() const {return S0;}
int GBM::get_N() const {return N;}
Matrix GBM::get_var() const {return Var;}
double GBM::get_r() const{return r;}


//Assignment operator
void GBM::operator=(const GBM& model){
    if(this!=&model){
        N = model.get_N();
        S0 = model.get_S0();
        r = model.get_r();
        Var = model.get_var();
    }
}


//Compute one step ahead
Matrix GBM::step(Matrix currentX, double h){
    double* drift = this->get_drift();
    for(int i=0;i<this->N;i++){drift[i]*=h;} //mu*Delta_t
    NormalSampler_multi Z(drift, h*Var);
    Matrix expo = exp(Z.getnumber());
    return HadamardProd(currentX, expo);
}


Matrix GBM::samplepath(double T, int nbSteps){
    Matrix path(N, nbSteps+1); //Return matrix
    for(int k=0; k<N; k++){path(k, 0)=S0.get_element(k,0);} //The first column is S0

    double h = T/(float)nbSteps; //Time step

    for(int i=1; i<=nbSteps; i++){
        Matrix currentS = path.get_exctract(0, N-1, i-1,i-1); //Take the previous stock price as the currentX
        path.modif_matrix(this->step(currentS, h), 0, i); //Compute one step ahead and put it in the newe column
    }
    return path;
}


Matrix GBM::samplepath_mean_var(double T, int nbSteps, int nbSimuls){

    Matrix means(this->N, nbSteps+1); //Sample Means for each time step tj
    Matrix cov_terminal(this->N, this->N); //Sample variance for the terminal time step tM=T

    Matrix simuls_terminal(this->N, nbSimuls); //Keep the last time step simulations in memory in order to compute the variance

    //Simulation and computation of the sample mean for each tj
    for(int j=0; j<nbSimuls; j++){
        Matrix paths = this->samplepath(T, nbSteps); //Generate Sample paths
        means = means + paths; //Update the mean
        simuls_terminal.modif_matrix(paths.get_exctract(0, this->N-1, nbSteps, nbSteps), 0, j); //Keep the last time step in memory
    }
    means = means/(float)nbSimuls;


    //Computation of the sample covariance for the terminal time step T
    for(int j=0; j<nbSimuls; j++){
        //cov = (x_i-m)*(x_i-m)^T
        cov_terminal = cov_terminal + (simuls_terminal.get_exctract(0,this->N-1, j,j) - means.get_exctract(0, this->N-1, nbSteps,nbSteps))*(simuls_terminal.get_exctract(0,this->N-1, j,j)-means.get_exctract(0, this->N-1, nbSteps,nbSteps)).T();
    }
    cov_terminal = cov_terminal/(float)(nbSimuls-1);


    //Compute theoretical mean E[S_t] = S_0exp(rt)
    Matrix theoretical_mean(this->N, nbSteps+1);
    double h=T/nbSteps;
    double t=0;
    for(int j=0; j<=nbSteps; j++){
        theoretical_mean.modif_matrix(S0*exp(r*t), 0, j); //E[S_j*h] = S_0exp(r(j*h))
        t+=h;
    }

    //Compute theoretical variance Cov(S(t)_i, S(t)_j) = S(0)_iS(0)_j*exp(2rt)*exp(sigma_i*sigma_j*rho_ij) - 1
    Matrix theoretical_var(this->N, this->N);
    for(int i=0; i<this->N; i++){
        for(int j=0; j<this->N; j++){
            theoretical_var(i,j) = S0.get_element(i,0)*S0.get_element(j,0)*exp(2*r*T)*(exp( Var.get_element(i,j)*T ) - 1);
        }
    }

    //Print results
    cout << "Sample mean for the whole path :\n";
    means.print_matrix();
    cout << "vs theoretical mean:\n";
    theoretical_mean.print_matrix();
    cout << "Absolute differences (%)\n";
    (100*(abs(means-theoretical_mean)/theoretical_mean)).print_matrix();
    cout << "Sample variance terminal value:\n";
    cov_terminal.print_matrix();
    cout << "vs theoretical variance:\n";
    theoretical_var.print_matrix();
    cout << "1 - ratio theoretical/sample variances (%):\n";
    ( 100*(ones(this->N, this->N) - (theoretical_var/cov_terminal))).print_matrix();

}

