#ifndef MATRIX_H
#define MATRIX_H


/**
    *A class to handle two-dimensional arrays (Matrices) efficiently.
    *This class is fundamental in our work, since we will use it for the Cholesky Decomposition, and for generating
    *multivariate random variables.
    *The Cholesky Decompostion of a Matrix is directly implemented as a member function of the class.

*/

class Matrix{
public:

    /**
        *Compute the Cholesky Decomposition of the matrix.
        *This will return a matrix A, lower-triangular such that AA^T=M.
        *The matrix M considered need to be squared and positive semi-definite for this to work.
        *@return The Cholesky Decomposition of the matrix.
    */
    Matrix CholeskyDecomp(void); //Cholesky factorisation


    /**
        *Constructor.
        *The memory of the 2d-array will be allocated dynamically.
        *@param  nbRows the number of rows in the array (0 by default).
        *@param nbColumns the number of columns in the array (0 by default).
    */
    Matrix(int nbRows=0, int nbColumns=0);

    /**
        *Copy constructor.
        *@param a the matrix to be copied.
    */
    Matrix(const Matrix& a);

    /**
        *Constructor using an existing pointer to an array of doubles.
        *This will create a 1d-matrix (1*n), with the elements equal to the elements in the array of doubles.
        *Note that this will create a column vector, NOT a row vector.
        *@param  a a pointer to an array of doubles.
        *@param d the number of elements in the array of doubles to be copied.
    */
    Matrix(const double* a, int d);

    ~Matrix();/**< Destructor.*/

    void operator=(const Matrix& A);/**< Overloading of the assignment operator.*/

    /**
        *Overloading of the operator ().
        *This will return the (i,j) element of the matrix by reference.
        *This can be used to modify/set the element (i,j).
        *@param i the row position of the needed element.
        *@param i the column position of the needed element.
        *@return the reference to the (i,j)-element in the matrix.
    */
    double& operator()(int i, int j);


    //getters
    int get_nbRows(void) const; /**< Get access of the (i,j)-element without modifying it.*/
    int get_nbColumns(void) const; /**< Get the number of columns.*/
    double get_element(int i, int j) const; /**< Get the number of rows.*/


    /**
        *Get access to a sub-matrix of the matrix.
        *This will return the sub-matrix starting at (i_start, j_start) and stopping at (i_stop-1, j_stop-1).
        *@param i_start the position of the top-left element of the sub-matrix.
        *@param i_stop the position of the bottom-left element of the sub-matrix.
        *@param j_start the position of the top-right element of the sub-matrix.
        *@param j_stop the position of the bottom-right element of the sub-matrix.
        *@return the sub-matrix.

    */
    Matrix get_exctract(int i_start, int i_stop, int j_start, int j_stop) const;


    /**
        *Modify a sub-matrix of the matrix.
        *This will put the matrix A into the initial matrix, starting at position (i_start, j_start).
        *Hence, the number of rows of A should be smaller than the number of rows of A - i_start, and same for
        *the number of columns.
        *@param i_start the row position at which we want to start modifying the matrix.
        *@param j_start the column position at which we want to start modifying the matrix.
        *@param A the matrix that will be put into the matrix we want to modify, with the top-left at position (i_start, j_start).

    */
    void modif_matrix(const Matrix& A, int i_start, int j_start);


    //Useful functions
    void print_matrix(void) const; /**< A function to print the matrix on the screen*/
    Matrix T(void); /**< A function that computes and return the transpose of the matrix.*/
    bool isSymmetric(void); /**< A function that tests if the matrix is symmetric.*/



protected:
    int n; /**< Number of rows in the matrix.*/
    int p; /**< Number of columns in the matrix.*/
    double **pt; /**< Pointer to elements of the matrix.*/

};

#endif // MATRIX_H
