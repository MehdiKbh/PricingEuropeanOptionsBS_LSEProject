#include "matrix.h"
#include "MatrixOperations.h"
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <cmath>
using namespace std;

//Cholesky Factorisation for PB1
Matrix Matrix::CholeskyDecomp(void){

    //Test if the Matrix is symmetric, and exit if it is not
    if(this->isSymmetric()==false){
        cout << "The matrix should be symmetric";
        exit(EXIT_FAILURE);
    }


    int n = this->n; //Size of the square matrix
    Matrix A(n,n); //Return Matrix

    for(int j=0; j<n; j++){
        for(int i=j; i<n; i++){
            //Diagonal Elements
            if(i==j){
                if(j>0){
                    Matrix row_i = A.get_exctract(i,i,0,j-1); //Get the i-th row of the matrix
                    double ps = (row_i*row_i.T())(0,0); //Compute the sum of Mij^2 for the i-th line
                    A(i,j) = sqrt(this->get_element(i,j) - ps);
                }
                else{
                    //First diagonal element
                    A(i,j) = sqrt(this->get_element(i,j));
                }
            }
            //Non-diagonal elements
            else{
                if(j>0){
                    Matrix row_i = A.get_exctract(i,i,0,j-1); //i-th row of the matrix
                    Matrix row_j = A.get_exctract(j,j,0,j-1).T(); //j-th row of the matrix
                    double ps = (row_i*row_j)(0,0); //Sum of the Mik*Mjk for i-th row and j-th row
                    A(i,j) = (this->get_element(i,j)-ps)/(A.get_element(j,j));
                }
                else{
                    A(i,j) = (this->get_element(i,j))/(A.get_element(j,j));
                }

            }
        }
    }
    return A;
}

//Get access to a sub-matrix -- Used in the Cholesky Decomposition member function above
Matrix Matrix::get_exctract(int i_start, int i_stop, int j_start, int j_stop) const{

    //Return a null matrix if indexes are not consistent
    if(i_start>i_stop || j_start>j_stop){
        Matrix exctr(0,0);
        return exctr;
    }

    //Check if we can indeed extract such a sub-matrix
    if(i_start<0 || j_start<0 || i_stop>=this->n || j_stop>=this->p){
        cout << "Index out of range\n";
        exit(EXIT_FAILURE);
    }

    //If indexes are ok and not too large
    else{
        int n_extr = (i_stop - i_start)+1; //The last row number of the matrix that we want to extract
        int p_extr = (j_stop - j_start)+1; //The last column number of the matrix that we want to extract

        Matrix exctr(n_extr, p_extr); //Return Matrix
        for(int i=0; i<n_extr; i++){
            for(int j=0; j<p_extr; j++){
                exctr(i,j) = this->get_element(i+i_start, j+j_start); //Do not modify the initial matrix
            }
        }
        return exctr;
    }
}


//Default Constructor
Matrix::Matrix(int nbRows, int nbColumns){
    n = nbRows; //Number of rows
    p = nbColumns; //Number of columns
    pt = new double*[n]; //Dynamic memory allocation of the elements of the matrix
    //pt is a pointer to pointers
    for(int i=0; i<n; i++){
        pt[i] = new double[p];}

    //Initialize the array to 0
    for(int i=0; i<n; i++){
        for(int j=0; j<p; j++){
            pt[i][j] = 0;
        }
    }
}

//Copy Constructor
Matrix::Matrix(const Matrix& a){
    n = a.n;
    p = a.p;
    pt = new double*[n];
    for(int i=0; i<n; i++){
        pt[i] = new double[p];
        for(int j=0; j<p; j++){
            pt[i][j] = a.pt[i][j];
        }
    }
}

//Create a 1d-matrix using the elements of an array of doubles
Matrix::Matrix(const double* a, int d){
    //This convert an array of doubles into a column vector - NOT INTO A ROW VECTOR!!
    n = d; //Number of lines
    p = 1; //Only one columns -- column vector
    pt = new double*[n];
    for(int i=0; i<n; i++){
        pt[i] = new double[p];
        for(int j=0; j<p; j++){
            pt[i][j] = a[i];
        }
    }
}

//Destructor
Matrix::~Matrix(){delete[] pt;}

//Overloading assignment operation
void Matrix::operator=(const Matrix& a){
    if(this!=&a){
        delete[] pt; //Delete existing memory
        //New number of rows and columns
        n = a.get_nbRows();
        p = a.get_nbColumns();
        //Create memory with the size of a
        pt = new double*[n];
        for(int i=0; i<n; i++){
            pt[i] = new double[p];}

        //Initialize the elements of the array to the elements of a
        for(int i=0; i<n; i++){
            for(int j=0; j<p; j++){
                pt[i][j] = a.get_element(i,j);
            }
        }
    }
}


//Getters
int Matrix::get_nbRows(void) const {return n;}
int Matrix::get_nbColumns(void) const {return p;}


//Access to the elements of the matrix by reference
//Since it returns elements by reference, it can be used to get elements as well as to modify them
double &Matrix::operator()(int i, int j){
    if ((i < 0)||(i>=n)||(j < 0)||(j>=p))
    {
        cout << "Index out of bounds! \n";
        exit(EXIT_FAILURE);
    }
    return pt[i][j];
}


//Get elements without modifying them
//Different of the () operator, since it does NOT modify the elements
double Matrix::get_element(int i, int j) const
{   if ((i < 0)||(i>=n)||(j < 0)||(j>=p))
    {
        cout << "Index out of bounds! \n";
        exit(EXIT_FAILURE);
    }
    return pt[i][j];
    }


//Modify only a part of the matrix
//The modification will delete the considered elements of the matrix, and set the matrix A instead
//Caution: this will modify a part of the matrix, and it CAN NOT be undone!
void Matrix::modif_matrix(const Matrix& A, int i_start, int j_start){

    //Check if A can fit into the matrix, starting at position (i_start, j_start)
    if(A.get_nbRows()+i_start>this->get_nbRows() || A.get_nbColumns()+j_start>this->get_nbColumns()){
        cout << "The matrix is too large to fit in this position\n";
        exit(EXIT_FAILURE);
    }
    //Put the matrix A
    for(int i=i_start; i<i_start+A.get_nbRows(); i++){
        for(int j=j_start; j<j_start+A.get_nbColumns(); j++){
            pt[i][j] = A.get_element(i-i_start,j-j_start);
        }
    }
}


//Print a matrix on the screen
void Matrix::print_matrix(void) const{
    cout << "[";
    for(int i=0; i<n; i++){
            cout << "[";
        for(int j=0; j<p; j++){
            cout << this->get_element(i,j);
            if(j<p-1){cout  << ",";}
        }
        cout << "]";
        if(i<n-1){cout << ",\n";}
        else{}
    }
    cout << "]\n";
    cout << endl;
}


//Compute and return the transpose of a Matrix
Matrix Matrix::T(void){
    Matrix Transpose(this->p, this->n);
    for(int i=0; i<this->n; i++){
        for(int j=0; j<this->p; j++){
            Transpose(j,i) = (*this)(i,j);
        }
    }
    return Transpose;
}


//True if a  square matrix is symmetric
bool Matrix::isSymmetric(void){
    //Check if the matrix is indeed a square matrix
    if(this->n!=this->p){
        cout << "The matrix should be square \n";
        exit(EXIT_FAILURE);
    }

    for(int i=0; i<this->n; i++){
        for(int j=i+1; j<this->n; j++){
            //Note that we avoid to use != essentially to be able to handle matrix that have small numerical errors
            if( abs((*this)(i,j)-(*this)(j,i))>1e-8 ){
                return false;
            }
        }
    }
    return true;
}
