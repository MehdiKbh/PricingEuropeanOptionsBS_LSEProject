#ifndef SAMPLER_H
#define SAMPLER_H

#include "matrix.h"
#include "MatrixOperations.h"

/**
    *Template virtual class to handle a sampler of random variable.
    *A template has been used to be able to handle univariate rv (double) and multivariate rv (matrix).
    *This class is inspired by the programming session 2.
*/

template <class T>
class Sampler
{
      public:
             virtual T getnumber()=0; /**<Return one sample of the random variable. Pure virtual function.*/
};

/**
    *A class to generate random numbers following a univariate uniform distribution on the interval (a,b).
*/
class UniformSampler:public Sampler<double>
{
      public:
              /**
                *Constructor.
                *@param a lower-bound of the interval on which the rv follows a uniform distribution.
                *@param b upper-bound of the interval on which the rv follows a uniform distribution.
              */
              UniformSampler(double a=0.0, double b=1.0);
              virtual double getnumber(); /**<Return one sample of the random variable.*/

      private:
              double low_bound; /**<Lower-bound of the interval on which the rv follows a uniform distribution.*/
              double up_bound;  /**<Upper-bound of the interval on which the rv follows a uniform distribution.*/
};


/**
    *A class to generate random numbers following a univariate normal distribution.
*/
class NormalSampler_unidim:public Sampler<double>
{
   public:
           /**
                *Constructor.
                *@param mean the mean of the normal rv.
                *@param stdev the standard deviation of the normal rv.
            */
           NormalSampler_unidim(double mean=0.0, double stdev=1.0);
           virtual double getnumber(); /**<Return one sample of the random variable.*/

   private:
            double mu; /**<Mean*/
            double sigma; /**<Standard deviation*/
};


//Matrix of size 1x1, containing 1.
const Matrix UNIT = unit_mat(1);
const double ZERO_MEAN[] = {0.};

/**
    *A class to generate random numbers following a multivariate normal distribution.
*/
class NormalSampler_multi:public Sampler<Matrix>
{
   public:

           /**
                *Constructor.
                *@param m a pointer to an array containing the mean vector of the normal rv.
                *@param VarCovar the Variance-Covariance matrix of the norml rv.
            */
           NormalSampler_multi(const double* m=ZERO_MEAN, Matrix VarCovar=UNIT);

           /**
                *Return one sample of the random variable.
                *Note that here the sample is of size n.
            */
           virtual Matrix getnumber();

   private:
            Matrix mu; /**<Mean*/
            Matrix var; /**<Variance matrix*/
            Matrix standard_normal_multi(); /**<A function to generate random numbers from a standard normal multivariate distribution.*/


};

#endif

