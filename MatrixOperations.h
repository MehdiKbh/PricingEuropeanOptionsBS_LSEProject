#ifndef MATRIXOPERATIONS_H
#define MATRIXOPERATIONS_H

/**
    *Operators and functions that can be used on a matrix.
    *This should always be included with the matrix.h file.

*/


#include "matrix.h"

/**
    *Takes volatilities and correlations for n stocks, and computes the covariance matrix Mij = sigma_i*sigma_j*rho_ij.
    *Caution: The volatilities should be ordered sigma_1, sigma_2, ...
    *The correlations should be ordered rho_12, rho_13, ..., rho_1n, rho_23, rho_24, ...
    *@param n the number of stocks
    *@param sigma a pointer to an array containing the n volatilities of each stock
    *@param rho a pointer to an array containing the (n choose 2) correlations between the stocks
*/
Matrix BuildCovarianceMatrix(const double* sigma, const double* rho, int n);


//Operations between matrices
Matrix operator+(const Matrix, const Matrix); /**<Add matrices elementwise*/
Matrix operator-(const Matrix, const Matrix); /**<Subtract matrices elementwise*/
Matrix HadamardProd(const Matrix A, const Matrix B); /**<Multiply matrices elementwise. This is NOT the matrices product. */

/**
    *Matrices product.
    *This is NOT the elementwise multiplication (Hadamard product), or whatever other product between matrices.
    *@see HadamardProd(const Matrix A, const Matrix B).
*/
Matrix operator*(const Matrix, const Matrix);
Matrix operator/(const Matrix, const Matrix); /**<Divide matrices elementwise.*/


//Operations with a scalar
Matrix operator*(const double, const Matrix); /**<Left-multiplication -- Multiply all elements of a matrix by a scalar.*/
Matrix operator*(const Matrix, const double); /**<Right-multiplication -- Multiply all elements of a matrix by a scalar.*/
Matrix operator/(const Matrix, const double); /**<Divide all elements of a matrix by a scalar.*/


Matrix unit_mat(int n); /**<Return the identity matrix of size n/*/
Matrix ones(int n, int p); /**<Return a matrix which all elements are 1.*/


//Maths functions
Matrix exp(const Matrix &A); /**<Return the matrix of exponential of all elements*/
Matrix abs(const Matrix &A); /**<Return the matrix of absolute value of all elements*/


#endif // MATRIXOPERATIONS_H
