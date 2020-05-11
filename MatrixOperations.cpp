#include "MatrixOperations.h"
#include "matrix.h"

#include <iostream>
#include <cstdlib>
#include <math.h>
using namespace std;


//Scalar operators

//Right-Multiplication by a scalar
Matrix operator*(const Matrix A, const double c){
    int n = A.get_nbRows(); int p = A.get_nbColumns(); //Size
    Matrix prod_mat(n, p); //Result matrix
    for(int i=0; i<n; i++){
        for(int j=0; j<p; j++){
            prod_mat(i,j) = c*A.get_element(i,j);
        }
    }
    return prod_mat;
}

//Left-Multiplication by a scalar
Matrix operator*(const double c, const Matrix A){return A*c;}

//Division of all elements by a scalar
Matrix operator/(const Matrix A, const double c){return (1./c)*A;}


//Matrix operators

//Addition elementwise
Matrix operator+(const Matrix A, const Matrix B){

    //Check if matrices are of same size
    if(A.get_nbRows()!=B.get_nbRows()|| A.get_nbColumns()!= B.get_nbColumns()){
        cout << "The matrices should be of same size to add them \n";
        exit(EXIT_FAILURE);
    }


    int n = A.get_nbRows(); int p = A.get_nbColumns(); //Size
    Matrix sum_mat(n, p); //Result matrix
    for(int i=0; i<n; i++){
        for(int j=0; j<p; j++){
            sum_mat(i,j) = A.get_element(i,j) + B.get_element(i,j);
        }
    }
    return sum_mat;
}

//Subtraction elementwise
Matrix operator-(const Matrix A, const Matrix B){return (-1.)*B + A;}


//Matrices product
Matrix operator*(const Matrix A, const Matrix B){

    //Check the sizes=
    if(A.get_nbColumns()!=B.get_nbRows()){
        cout << "The number of columns of A should be equal to the number of rows of B to take the matrix product";
        exit(EXIT_FAILURE);
    }

    Matrix prod_mat(A.get_nbRows(), B.get_nbColumns()); //Result matrix
    for(int i=0; i<A.get_nbRows(); i++){
        for(int j=0; j<B.get_nbColumns(); j++){
            for(int k=0; k<A.get_nbColumns(); k++){
                prod_mat(i,j)+= A.get_element(i,k)*B.get_element(k,j);
            }
        }
    }
    return prod_mat;
}

//Divide elementwise
Matrix operator/(const Matrix A, const Matrix B){
    //Check if matrices are of same size
    if(A.get_nbRows()!=B.get_nbRows()|| A.get_nbColumns()!= B.get_nbColumns()){
        cout << "The matrices should be of same size to divide them\n";
        exit(EXIT_FAILURE);
    }

    int n = A.get_nbRows(); int p = A.get_nbColumns(); //Size
    Matrix div_mat(n, p); //Return matrix
    for(int i=0; i<n; i++){
        for(int j=0; j<p; j++){
            div_mat(i,j) = A.get_element(i,j) / B.get_element(i,j);
        }
    }
    return div_mat;
}

//Multiply elementwise
Matrix HadamardProd(const Matrix A, const Matrix B){
    //Check if matrices are of same size
    if(A.get_nbRows()!=B.get_nbRows()|| A.get_nbColumns()!= B.get_nbColumns()){
        cout << "The matrices should be of same size to multiply them \n";
        exit(EXIT_FAILURE);
    }

    int n = A.get_nbRows(); int p = A.get_nbColumns(); //Size
    Matrix prod(n, p); //Return matrix
    for(int i=0; i<n; i++){
        for(int j=0; j<p; j++){
            prod(i,j) = A.get_element(i,j) * B.get_element(i,j);
        }
    }
    return prod;
}

//Build and return the covariance matrix, given the volatilities and correlations
Matrix BuildCovarianceMatrix(const double* sigma, const double* rho, int n){

    //Convert rho into a Matrix to make the assignment easier
    Matrix correl(n,n); //Correlation matrix
    int k=0;
    for(int i=0; i<n; i++){
        correl(i,i) = 1;
        for(int j=i+1; j<n; j++){
            correl(i,j) = rho[k];
            correl(j,i) = rho[k];
            k++;
        }
    }

    //Compute the covariance matrix
    Matrix covariance(n,n);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            covariance(i,j) = sigma[i]*sigma[j]*correl(i,j);
        }
    }
    return covariance;
}

//Exponential elementwise
Matrix exp(const Matrix &A){
    Matrix expo(A.get_nbRows(), A.get_nbColumns());
    for(int i=0; i<A.get_nbRows(); i++){
        for(int j=0; j<A.get_nbColumns(); j++){
            expo(i,j) = exp(A.get_element(i,j));
        }
    }
    return expo;
}

//Absolute value elementwise
Matrix abs(const Matrix &A){
    Matrix abs(A.get_nbRows(), A.get_nbColumns());
     for(int i=0; i<A.get_nbRows(); i++){
        for(int j=0; j<A.get_nbColumns(); j++){
            abs(i,j) = (A.get_element(i,j)>=0)?A.get_element(i,j):-A.get_element(i,j);
        }
    }
    return abs;
}


//Return the identity matrix of size n
Matrix unit_mat(int n){
    Matrix unit(n,n);
    for(int i=0;i<n;i++){
        unit(i,i) = 1;
    }
    return unit;
}

Matrix ones(int n, int p){
    Matrix one(n,p);
    for(int i=0;i<n;i++){
        for(int j=0; j<p; j++){
            one(i,j) = 1;
        }
    }
    return one;
}

