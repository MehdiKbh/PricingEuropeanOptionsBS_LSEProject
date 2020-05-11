#include "Sampler.h"
#include"matrix.h"
#include "MatrixOperations.h"

#include<iostream>
#include<cstdlib>
#include<cmath>
using namespace std;


//Uniform
UniformSampler::UniformSampler(double a, double b)
{
    low_bound=a;
    up_bound=b;
}

//Generate random number using the function rand()
//Inspired by the programming session 2
double UniformSampler::getnumber()
{
       int myrand=rand();
       while((myrand==0)||(myrand==RAND_MAX)){myrand = rand(); } //Number in (0, RAND_MAX).

       double myuni = myrand/static_cast<double>(RAND_MAX); //Create a number in (0,1).
       return low_bound + myuni*(up_bound-low_bound); //Rescale
}


//Normal unidimensional
NormalSampler_unidim::NormalSampler_unidim(double mean, double stdev)
{
   mu=mean;
   sigma=stdev;
}

//Generate random numbers using the Box-Muller method
//Inspired by the programming session 2
double NormalSampler_unidim::getnumber()
{
        UniformSampler U1;
        double R = -2*log(U1.getnumber());
        UniformSampler U2(0., 2*M_PI);

        double X1 = sqrt(R)*cos(U2.getnumber());
        return mu + sigma*X1;
}



//Normal multidimensional
NormalSampler_multi::NormalSampler_multi(const double* m, Matrix VarCovar){

    Matrix temp(m, VarCovar.get_nbRows());
    mu = temp;
    var = VarCovar;
}


//Uses the standard_normal_multi() method (see above)
Matrix NormalSampler_multi::getnumber(){
    Matrix Z = this->standard_normal_multi();
    Matrix A = var.CholeskyDecomp();
    Matrix res = mu + A*Z;
    return res;
};

//Generate n independent rv following a normal distribution
Matrix NormalSampler_multi::standard_normal_multi(){
    Matrix Z(var.get_nbRows(),1); //Column vector

    //Generate independent normal rv
    for(int i=0; i<var.get_nbRows(); i++){
        Z(i,0) = NormalSampler_unidim().getnumber();
    }
    return Z;
}
