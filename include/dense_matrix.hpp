/*
    dense_matrix.hpp

    General matrix class intended to work with
    (i) Intel® Math Kernel Library
    (ii) openBLAS
    in terms of utilizing  matrix storage formats,
    memory functions and inversion & multiplication routines.
*/

#ifndef cs_matrix_hpp__
#define cs_matrix_hpp__

#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include <algorithm>
#include <utility>
#include <iterator>
#include <cstring>
#include <vector>
#include <typeinfo>
#include <random>
#include <stdlib.h>

#define MKL_INT size_t

#ifdef intelmkl
#include "mkl.h"
#else
#include <cblas.h>
#include <mkl_lapacke.h> // TEMPORARY FIX ( "added mkl_")
#endif

#ifndef _min
#define _min(x, y) (((x) < (y)) ? (x) : (y))
#endif

#ifndef _max
#define _max(x, y) (((x) > (y)) ? (x) : (y))
#endif

#ifndef worksize
#define worksize 10000
#endif

namespace evolm
{

    template <typename T>
    class matrix
    {

    public:
        /* CONSTRUCTORS & DESTRUCTOR */

        matrix(size_t row, size_t col); /* Constructor for rectangular matrix. */
        matrix(size_t lda);             /* Constructor for square symmetric half-store (compact) matrix. */
        matrix(const char *type);       /* Default constructor: rectangular (type = 'r'|| 'R')/symmetric (type = 's'|| 'S') matrix, no memmory allocated. */
        matrix();                       /* Default constructor: rectangular matrix, no memmory allocated. */
        matrix(const matrix &obj);      /* Copy constructor. */
        ~matrix();                      /* Destructor. */

        /* METHODS */

        T &at(size_t atRow, size_t atCol);    /* Get/put element from/in a matrix. */
        void resize(size_t row, size_t col);  /* Resizes/allocates memmory for A. */
        void resize(size_t lda);              /* Resizes/allocates memmory for A; overloaded method for symmetrical matrix. */
        void print(std::string whiichMatrix); /* Prints part of a matrix into a LOG file. */
        void printf(std::string whiichMatrix,
                    bool append);                     /* Prints all of a matrix into a specified file file name. */
        void scale(T val);                            /* Scaling matrix by scalar: A = A*val. */
        size_t size() const;                          /* Gives total number of elements in a matrix. */
        size_t capacity();                            /* Gives total number of allocated elements in a matrix. */
        matrix<size_t> shape() const;                 /* Returns a vector with number of rows and columns. */
        void clear();                                 /* Frees an allocated memmory of A. */
        void fclear();                                /* Removes (if exists) a binary file associated with A. */
        void fclear(const std::string &fname);        /* Removes (if exists) a named binary file. */
        bool empty();                                 /* Checks if memmory is allocated to A. */
        void symtorec();                              /* Transform symmetric matrix in compact form to rectangular form. */
        void rectosym();                              /* Transform square matrix to triangular compact form (only for symmetric matrices). */
        void transpose();                             /* Transpose matrix. */
        void fwrite(const std::string &fname);        /* Save matrix to disk, and not clear memory. */
        void fwrite(const std::string &fname, std::vector<std::int64_t> &vect_data); /* Save matrix to disk with additional data in vector, and not clear memory. */
        void fwrite(const std::string &fname, std::vector<std::string> &vect_data);  /* Save matrix to disk with additional data in vector, and not clear memory. */
        void fwrite();                                /* Move matrix to the disk and clear memory. */
        void invert();                                /* Matrix inversion. */
        void lchol();                                 /* Cholesky factorisation, gives lower triangular outpur. */
        void linsolve(matrix<T> &rhs);                /* Solving vinear system A * x = rhs, result: x overwrites rhs; A overwrited by LU decomposition. This method is not making copies of arrays. */
        void linsolve_fact(matrix<T> &rhs);           /* Solving vinear system A * x = rhs, result: x overwrites rhs; A is already factorized. This method is not making copies of arrays. */
        void fread();                                 /* Restore matrix saved on disk into the memory. */
        void fread(const std::string &fname);         /* Restore matrix from the disk into the memory. */
        void fread(const std::string &fname, std::vector<std::int64_t> &vect_data); /* Restore matrix from the disk into the memory with additional data into vector. */
        void fread(const std::string &fname, std::vector<std::string> &vect_data); /* Restore matrix from the disk into the memory with additional data into vector. */
        bool is_ondisk();                             /* Checking if data is located on disk or in memory. */
        matrix<T> fget(size_t irow[], size_t icol[]); /* Reads and returns just part of data from a file.*/
        matrix<T> fget();                             /* Overloaded. Reads and returns all of data from a file.*/
        bool eq(const matrix &rhs);                   /* Compare dimensins and shapes of two matrix objects. */
        bool is_symmetric();                          /* Checks is matrix is symmetric */

        T *return_array()
        {
            try
            {
                return A;
            }
            catch (const std::exception &e)
            {
                std::cerr << e.what() << " in matrix<T>::return_array()" << '\n';
                throw e;
            }
            catch (...)
            {
                std::cerr << "Exception in matrix<T>::return_array()" << '\n';
                throw;
            }
        };

        void insert_array(T *B)
        {
            try
            {
                if (!allocated)
                {
                    throw std::string("The memory for Matrix is not allocated. Use resize(). matrix<T>::insert_array(T *B)");
                }

                size_t sz;

                if (!compact)
                    sz = numRow * numCol;
                else
                    sz = (numCol * numCol + numCol) / 2;

                for (size_t i = 0; i < sz; i++)
                    A[i] = B[i];
            }
            catch (const std::exception &e)
            {
                std::cerr << e.what() << " in matrix<T>::insert_array(T *)" << '\n';
                throw e;
            }
            catch (...)
            {
                std::cerr << "Exception in matrix<T>::insert_array(T *)" << '\n';
                throw;
            }
        };

        void to_vector(std::vector<T> &vect)
        {
            try
            {
                size_t sz;

                if (!compact)
                    sz = numRow * numCol;
                else
                    sz = (numCol * numCol + numCol) / 2;

                for (size_t i = 0; i < sz; i++)
                    vect.push_back(A[i]);
            }
            catch (const std::exception &e)
            {
                std::cerr << e.what() << " in matrix<T>::to_vector(std::vector<T> &)" << '\n';
                throw e;
            }
            catch (...)
            {
                std::cerr << "Exception in matrix<T>::to_vector(std::vector<T> &)" << '\n';
                throw;
            }
        };

        void from_vector(std::vector<T> &vect)
        {
            try
            {
                size_t sz_vect = vect.size();

                size_t sz;

                if (allocated)
                {
                    if (!compact)
                        sz = numRow * numCol;
                    else
                        sz = (numCol * numCol + numCol) / 2;

                    if (sz != sz_vect)
                        throw std::string("The size of pre-allocated Matrix is not the same as the size of the input vector. matrix<T>::from_vect(std::vector<T> &)");

                    for (size_t i = 0; i < sz; i++)
                        A[i] = vect[i];
                }
                else
                {
                    resize(sz_vect, 1);

                    sz = numRow;

                    if (sz != sz_vect)
                        throw std::string("The size of post-allocated Matrix is not the same as the size of the input vector. matrix<T>::from_vect(std::vector<T> &)");

                    for (size_t i = 0; i < sz; i++)
                        A[i] = vect[i];
                }
            }
            catch (const std::exception &e)
            {
                std::cerr << e.what() << " in matrix<T>::from_vector(std::vector<T> &)" << '\n';
                throw e;
            }
            catch (const std::string &e)
            {
                std::cerr << e << " in matrix<T>::from_vector(std::vector<T> &)" << '\n';
                throw e;
            }
            catch (...)
            {
                std::cerr << "Exception in matrix<T>::from_vector(std::vector<T> &)" << '\n';
                throw;
            }
        };

        void from_vector2d(std::vector<std::vector<T>> &vect)
        {
            try
            {
                size_t sz_vect = vect.size();
                size_t sz_vect2 = vect[0].size();

                size_t sz;

                if (allocated)
                {
                    if (!compact)
                        sz = numRow * numCol;
                    else
                    {
                        symtorec();
                        resize(numRow, numCol);
                        sz = numRow * numCol;
                    }

                    if (sz != sz_vect * sz_vect2)
                        throw std::string("The size of allocated Matrix is not the same as the size of the input vector. matrix<T>::from_vect2d(std::vector<std::vector<T>> &)");

                    size_t l = 0;
                    for (size_t i = 0; i < sz_vect; i++)
                    {
                        for (size_t j = 0; j < sz_vect2; j++)
                        {
                            A[l] = vect[i][j];
                            l = l + 1;
                        }
                    }
                }
                else
                {
                    resize(sz_vect, sz_vect2);

                    sz = numRow * numCol;

                    if (sz != sz_vect * sz_vect2)
                        throw std::string("The size of allocated Matrix is not the same as the size of the input vector. matrix<T>::from_vect2d(std::vector<std::vector<T>> &)");

                    size_t l = 0;
                    for (size_t i = 0; i < sz_vect; i++)
                    {
                        for (size_t j = 0; j < sz_vect2; j++)
                        {
                            A[l] = vect[i][j];
                            l = l + 1;
                        }
                    }
                }
            }
            catch (const std::exception &e)
            {
                std::cerr << e.what() << " in matrix<T>::from_vector2d(std::vector<std::vector<T>> &)" << '\n';
                throw e;
            }
            catch (const std::string &e)
            {
                std::cerr << e << " in matrix<T>::from_vector2d(std::vector<std::vector<T>> &)" << '\n';
                throw e;
            }
            catch (...)
            {
                std::cerr << "Exception in matrix<T>::from_vector2d(std::vector<std::vector<T>> &)" << '\n';
                throw;
            }
        };

        void to_vector2d(std::vector<std::vector<T>> &vect)
        {
            try
            {
                if (allocated)
                {
                    if (compact)
                    {
                        symtorec();
                        resize(numRow, numCol);
                    }

                    size_t l = 0;
                    for (size_t i = 0; i < numRow; i++)
                    {
                        std::vector<T> v0;
                        for (size_t j = 0; j < numCol; j++)
                        {
                            v0.push_back(A[l]);
                            l = l + 1;
                        }
                        vect.push_back(v0);
                    }
                }
                else
                {
                    throw std::string("The array to be copied to 2d std::vector is not allocated!");
                }
            }
            catch (const std::exception &e)
            {
                std::cerr << e.what() << " in matrix<T>::to_vector2d(std::vector<std::vector<T>> &)" << '\n';
                throw e;
            }
            catch (const std::string &e)
            {
                std::cerr << e << " in matrix<T>::to_vector2d(std::vector<std::vector<T>> &)" << '\n';
                throw e;
            }
            catch (...)
            {
                std::cerr << "Exception in matrix<T>::to_vector2d(std::vector<std::vector<T>> &)" << '\n';
                throw;
            }
        };

        /* TYPE CONVERTORS */

        matrix<float> _float();
        matrix<double> _double();
        matrix<int> _int();
        matrix<size_t> _size_t();
        matrix<bool> _bool();
        void cast_fget(size_t irow[], size_t icol[], matrix<float> &out);
        void cast_fget(size_t irow[], size_t icol[], matrix<double> &out);
        void cast_fget(size_t irow[], size_t icol[], matrix<int> &out);
        void cast_fget(size_t irow[], size_t icol[], matrix<size_t> &out);
        void cast_fget(size_t irow[], size_t icol[], matrix<bool> &out);

        void cast_fget(size_t irow[], size_t icol[], std::vector<std::vector<float>> &out);
        void cast_fget(size_t irow[], size_t icol[], float **out);

        /* OPERATORS */

        matrix operator+(const matrix &rhs); /* Overloaded '+' operator to add two matrix objects. */
        matrix operator+(const T &val);
        matrix operator<<(const matrix &rhs); /* Overloaded '<<' operator to combine two matrix objects, without transpose. */
        matrix operator>>(const matrix &rhs); /* Overloaded '<<' operator to combine two matrix objects, with transpose. */
        matrix operator-(const matrix &rhs);  /* Overloaded '-' operator to substract two matrix objects. */
        matrix operator-(const T val);        /* Overloaded '-' operator to substitute a scalar from a matrix object. */
        matrix operator^(const int val);      /* Overloaded '^' operator to multiply matrix by itself and find inversion. */
        matrix operator^(const char *val);    /* Overloaded '^' operator to transpose matrix. */
        matrix operator*(const matrix &rhs);  /* Overloaded '*' operator to multiply two matrix objects. */
        matrix operator*(const T &val);
        bool operator==(const matrix &rhs) const;       /* Compare complete equality of two matrix objects. */
        matrix &operator=(const matrix &rhs);           /* Overloaded assignment '=' operator. */
        T &operator()(size_t atRow, size_t atCol);      /* Access element of a matrix, memory reallocation is allowed. */
        T operator()(size_t atRow, size_t atCol) const; /* Access element of a matrix. */
        T &operator[](size_t i);                        /* Access element of a matrix. */
        T operator[](size_t i) const;                   /* Access element of a matrix. */

        /* VARIABLES & CONSTANTS */

        bool failbit; /* Error state flag. */
        int failinfo;

    private:

        lapack_int *ipiv_A; // pointer to an array that may contain pivoting indices due to linsolve_LU(...) method
        mutable bool is_allocated_ipv_A = false; // flag pointing to if ipiv_A is associated to an allocated memory

        /* METHODS */

        int allocate(size_t row, size_t col); /* Allocate memory for a rectangular matrix. */
        int allocate(size_t lda);             /* Allocate memory for a half-store (compact) symmetric matrix. */
        void resize();                        /* Resizes memmory allocated for A. */

        /* Interfaces to MKL/OpenBLAS routines */

        void dotprod(double *_A, double *B, double *C, MKL_INT rowA/*, MKL_INT rowB*/, MKL_INT colA, MKL_INT colB);
        void dotprod(float *_A, float *B, float *C, MKL_INT rowA/*, MKL_INT rowB*/, MKL_INT colA, MKL_INT colB);
        void inv_rec(double *_A, MKL_INT rowA, MKL_INT colA);
        void inv_rec(float *_A, MKL_INT rowA, MKL_INT colA);
        void linsolve_lu(double *_A, MKL_INT rowA, MKL_INT colA, double *_B, MKL_INT colB);
        void linsolve_lu(float *_A, MKL_INT rowA, MKL_INT colA, float *_B, MKL_INT colB);
        void linsolve_lu_fact(double *_A, MKL_INT rowA, /*MKL_INT colA,*/ double *_B, MKL_INT colB);
        void linsolve_lu_fact(float *_A, MKL_INT rowA, /*MKL_INT colA,*/ float *_B, MKL_INT colB);
        void linsolve_chol(double *_A, MKL_INT rowA, MKL_INT colA, double *_B, MKL_INT colB);
        void linsolve_chol(float *_A, MKL_INT rowA, MKL_INT colA, float *_B, MKL_INT colB);
        void linsolve_chol_fact(double *_A, MKL_INT rowA, /*MKL_INT colA,*/ double *_B, MKL_INT colB);
        void linsolve_chol_fact(float *_A, MKL_INT rowA, /*MKL_INT colA,*/ float *_B, MKL_INT colB);
        void inv_sym(double *_A, MKL_INT colA);
        void inv_sym(float *_A, MKL_INT colA);
        void get_lchol(double *_A, MKL_INT colA);
        void get_lchol(float *_A, MKL_INT colA);
        void get_lchol(double *_A, MKL_INT rowA, MKL_INT colA);
        void get_lchol(float *_A, MKL_INT rowA, MKL_INT colA);
        void gemmt_intrf(double *_A, double *B, MKL_INT rowA, MKL_INT colA, MKL_INT colB);
        void gemmt_intrf(float *_A, float *B, MKL_INT rowA, MKL_INT colA, MKL_INT colB);

        /* VARIABLES & CONSTANTS */

        T *A;                   /* The main matrix container. */
        size_t numRow;          /* Number of rows in A. */
        size_t numCol;          /* Number of columns in A. */
        size_t resizedElements; /* Number of allocated/resized elements in A. */
        bool allocated;         /* Flag indicated that the memory for A have been allocated. */
        bool rectangular;       /* TRUE if the matrix is rectangular. */
        bool symetric;          /* TRUE if the matrix is symmetric. */
        bool compact;           /* TRUE if for the symmetric matrix only a lower triangular part is stored. */
        bool ondisk;
        std::string debug_file = "MATRIX.log";
        std::string binFilename; /* Name of binary file to store A on disck. */
        std::fstream fA;

        const size_t matrix_type = 1001;

    public:
        /* FRIEND OPERATORS */

        /*
         Such operator's definitions are inlined in the class definition because operators don't exist outside the class definition
         until a template instantiation of the class is generated (which happens during the compilation process).
         Here we perform implicit type conversion on non-templated non-member functions.
         We create the non-template friend functions which are automatically created for each template instantiation of the class.
        */
        friend matrix operator*(const T val, const matrix &rhs)
        {
            /*
                Overloaded '*' operator to multiply matrix by scalar.

                Return value: matrix object.

                Example:

                    matrix <double> obj;    // empty matrix
                    matrix <double> M(n,m); // M is (n,m) matrix initialized by 0.

                    for (auto i = 0; i < M.size(); i++)
                        M[i] = 1.0;

                    obj = -2.0 * M;  // obj become (n,m) matrix where all elements are -2.0;
                                     // M remains unchanged (all elements of the matrix are 1.0)
            */

            if (rhs.ondisk)
                throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::operator*");

            return rhs * val;
        }

        friend matrix operator*(const matrix &lhs, const T val)
        {
            /*
                Overloaded '*' operator to multiply matrix by scalar.

                Return value: matrix object.

                Example:

                    matrix <double> obj;    // empty matrix
                    matrix <double> M(n,m); // M is (n,m) matrix initialized by 0.

                    for (auto i = 0; i < M.size(); i++)
                        M[i] = 1.0;

                    obj = M * (-2.0); // obj become (n,m) matrix where all elements are -2.0;
                                      // M remains unchanged (all elements of the matrix are 1.0)
            */

            if (lhs.ondisk)
                throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::operator*");

            matrix<T> C;
            if (!lhs.compact)
            {
                int status = C.allocate(lhs.numRow, lhs.numCol);
                if (status != 0)
                {
                    C.failbit = true;
                    throw std::string("Memory allocation error: matrix<T>::operator*");
                }

                C.allocated = true;
                C.resizedElements = lhs.numRow * lhs.numCol;
                C.numCol = lhs.numCol;
                C.numRow = lhs.numRow;
            }
            else
            {
                int status = C.allocate(lhs.numRow);
                if (status != 0)
                {
                    C.failbit = true;
                    throw std::string("Memory allocation error: matrix<T>::operator*");
                }

                C.allocated = true;
                C.resizedElements = (lhs.numRow * lhs.numRow + lhs.numRow) / 2;
                C.numCol = lhs.numRow;
                C.numRow = lhs.numRow;
            }

            auto n_threads = std::thread::hardware_concurrency();
            auto block_size = static_cast<unsigned int>(C.size() / (n_threads));

            if (block_size < worksize)
            {
                block_size = static_cast<unsigned int>(C.size());
                n_threads = 1;
            }

            // #pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < C.size(); i++)
                C.A[i] = lhs.A[i] * val;

            return matrix(C);
        }

        friend matrix operator+(const T lhs, const matrix &rhs)
        {
            /*
                Overloaded '+' operator to add a scalar to a matrix.

                Return value: matrix object.

                Example:

                    matrix <double> obj;    // empty matrix
                    matrix <double> M(n,m); // M is (n,m) matrix initialized by 0.

                    obj = -2.0 + M; // obj become (n,m) matrix where all elements are -2.0; M remains unchanged
            */

            if (rhs.ondisk)
                throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::operator+");

            return rhs + lhs;
        }

        friend matrix operator+(const matrix &lhs, const T rhs)
        {
            /*
                Overloaded '+' operator to add a scalar to a matrix.

                Return value: matrix object.

                Example:

                    matrix <double> obj;    // empty matrix
                    matrix <double> M(n,m); // M is (n,m) matrix initialized by 0.

                    obj = M + (-2.0); // obj become (n,m) matrix where all elements are -2.0; M remains unchanged
            */

            if (lhs.ondisk)
                throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::operator+");

            matrix<T> C;
            if (!lhs.compact)
            {
                int status = C.allocate(lhs.numRow, lhs.numCol);
                if (status != 0)
                {
                    C.failbit = true;
                    throw std::string("Memory allocation error: matrix<T>::operator+");
                }

                C.allocated = true;
                C.resizedElements = lhs.numRow * lhs.numCol;
                C.numCol = lhs.numCol;
                C.numRow = lhs.numRow;
            }
            else
            {
                int status = C.allocate(lhs.numRow);
                if (status != 0)
                {
                    C.failbit = true;
                    throw std::string("Memory allocation error: matrix<T>::operator+");
                }

                C.allocated = true;
                C.resizedElements = (lhs.numRow * lhs.numRow + lhs.numRow) / 2;
                C.numCol = lhs.numRow;
                C.numRow = lhs.numRow;
            }

            auto n_threads = std::thread::hardware_concurrency();
            auto block_size = static_cast<unsigned int>(C.size() / (n_threads));

            if (block_size < worksize)
            {
                block_size = static_cast<unsigned int>(C.size());
                n_threads = 1;
            }

            // #pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < C.size(); i++)
                C.A[i] = lhs.A[i] + rhs;

            return matrix(C);
        }
    };

    //}

    //===============================================================================================================

    template <typename T>
    void matrix<T>::transpose()
    {
        /*
           Handles only rectangular matrices (including full symmetric).
           No transpose of triangular matrix since we do not use 'U' format.

           Return value: none; modifies the container A of the calling object.

           Example:

                matrix <double> M(n,m); // M is (n,m) matrix initialized by 0.
                matrix <double> res;    // empty matrix

                for (auto i = 0; i < M.size(); i++)
                    M[i] = i;

                M.transpose(); // now M is (m,n) matrix
                res = M;       // now res is (m,n) matrix

        */

        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::transpose()");

        if (!compact)
        {
            matrix<T> C(numCol, numRow);

            /*auto n_threads = std::thread::hardware_concurrency();
            auto block_size = static_cast<unsigned int>(size() / (n_threads));

            if (block_size < worksize)
            {
                block_size = size();
                n_threads = 1;
            }*/

            // #pragma omp parallel for // schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < numRow; i++)
            {
                for (size_t j = 0; j < numCol; j++)
                {
                    C.A[j * numRow + i] = A[i * numCol + j];
                }
            }
            /*for (size_t n = 0; n < numRow*numCol; n++)
            {
                size_t i = n/numRow;
                size_t j = n%numRow;
                C.A[n] = A[j * numCol + i];
            }*/

            resize(numCol, numRow);
            *this = C;
            C.clear();
        }
        else
        {

            /*We also have to accept symmetric matrices in compact form.
            Due to the reasons described in the upper comment, we return the same matrix
            but in restored (not compact) format.*/

            symtorec();
        }
    }

    //===============================================================================================================

    template <typename T>
    T &matrix<T>::operator()(size_t atRow, size_t atCol)
    {
        /*
            Fast access operator (array-like).

            Return value: the element at the specified position in the container.

            Example:
                     matrix <double> M(2,2);
                     M(1,1) = 0.8;
                     double val = M(1,1);   // val = 0.8.
        */

#ifdef UTEST // use this only when debugging and testing
        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::operator()");
#endif
        if (!compact)
            return A[atRow * numCol + atCol];
        else
        {
#ifdef UTEST // use this only when debugging and testing
            if ( atCol > atRow )
                throw std::string("The column value is greater than the row value. This is not allowed for symmetric matrix in compact (L-store) format!");
#endif
            return A[atRow * (atRow + 1) / 2 + atCol];
        }
    }

    //===============================================================================================================

    template <typename T>
    T matrix<T>::operator()(size_t atRow, size_t atCol) const
    {
        /*
            Fast access operator (array-like).

            Return value: the element at the specified position in the container.

            Example:
                     matrix <double> M(2,2);
                     M(1,1) = 0.8;
                     double val = M(1,1);   // val = 0.8.
        */

#ifdef UTEST // use this only when debugging and testing
        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::operator()");
#endif
        if (!compact)
            return A[atRow * numCol + atCol];
        else
            return A[atRow * (atRow + 1) / 2 + atCol];
    }

    //===============================================================================================================

    template <typename T>
    T &matrix<T>::operator[](size_t i)
    {
        /*
            Fast access operator (direct).

            Return value: the element at the specified position in the container.

            Example:

                matrix <double> M(n,m);
                double val = M[ M.size()-1 ]; // val = 0.0

                for (auto i = 0; i < M.size(); i++)
                    M[ i ] = 1.0;

                val = M[ M.size()-1 ]; // val = 1.0

        */

#ifdef UTEST // use this only when debugging and testing
        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::operator[]");
#endif
        return A[i];
    }

    //===============================================================================================================

    template <typename T>
    T matrix<T>::operator[](size_t i) const
    {
        /*
            Fast access operator (direct).

            Return value: the element at the specified position in the container.

            Example:

                matrix <double> M(n,m);
                double val = M[ M.size()-1 ]; // val = 0.0

                for (auto i = 0; i < M.size(); i++)
                    M[ i ] = 1.0;

                val = M[ M.size()-1 ]; // val = 1.0

        */

#ifdef UTEST // use this only when debugging and testing
        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::operator[]");
#endif
        return A[i];
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::resize()
    {
        /*
            Private member. Reallocates memory for the container A.
            Note: row and column of a matrix have to be modified before calling the resize() method.

            Return value: none.
        */

        size_t sz;
        if (!compact)
            sz = numRow * numCol;
        else
        {
            size_t lda = _max(numRow, numCol);
            sz = (lda * lda + lda) / 2;
        }

#ifdef intelmkl
        A = (T *)mkl_realloc(A, sz * sizeof(T));
#else
        A = (T *)realloc(A, sz * sizeof(T));
#endif
        if (A == NULL)
        {
#ifdef intelmkl
            mkl_free(A);
#else
            free(A);
#endif
            allocated = false;
            failbit = true;
            throw std::string("Memory allocation error. matrix<T>::resize()");
        }
        allocated = true;
        resizedElements = sz;
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::resize(size_t row, size_t col)
    {
        /*
            Reallocates memmorey for the member array A.

            Return value: none.
            Example:

                matrix <double> M(3,5);
                size_t val;
                val = M.size(); // val = 15.
                M.resize(2,4);
                val = M.size(); // val = 8.
        */

        size_t sz;
        size_t lda;
        if (!compact)
        {
            sz = row * col;
        }
        else
        {
            lda = _max(row, col);
            sz = (lda * lda + lda) / 2;
        }
        if (allocated)
        {
#ifdef intelmkl
            A = (T *)mkl_realloc(A, sz * sizeof(T));
#else
            A = (T *)realloc(A, sz * sizeof(T));
#endif
            if (A == NULL)
            {
#ifdef intelmkl
                mkl_free(A);
#else
                free(A);
#endif
                allocated = false;
                failbit = true;
                throw std::string("Memory reallocation error. matrix<T>::resize(size_t, size_t)");
            }
        }
        else
        {
            if (!compact)
            {
                int status = allocate(row, col);
                if (status != 0)
                {
                    failbit = true;
                    throw std::string("Memory allocation error. matrix<T>::resize(size_t, size_t)");
                }
            }
            else
            {
                int status = allocate(lda);
                if (status != 0)
                {
                    failbit = true;
                    throw std::string("Memory allocation error. matrix<T>::resize(size_t, size_t)");
                }
            }
            allocated = true;
        }

        resizedElements = sz;
        if (!compact)
        {
            numRow = row;
            numCol = col;
            rectangular = true;
            symetric = false;
        }
        else
        {
            numRow = numCol = _max(row, col);
            rectangular = false;
            symetric = true;
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::resize(size_t lda)
    {
        /*
            Reallocates memmorey for the member array A.
            Overloaded method for symmetrical matrix in compact form.

            Return value: none.
            Example:

                matrix <double> M(3);
                size_t val;
                val = M.size(); // val = 6.
                M.resize(4);
                val = M.size(); // val = 10.
        */

        size_t sz = (lda * lda + lda) / 2;

        if (allocated)
        {
#ifdef intelmkl
            A = (T *)mkl_realloc(A, sz * sizeof(T));
#else
            A = (T *)realloc(A, sz * sizeof(T));
#endif
            if (A == NULL)
            {
#ifdef intelmkl
                mkl_free(A);
#else
                free(A);
#endif
                allocated = false;
                failbit = true;
                throw std::string("Memory reallocation error. matrix<T>::resize(size_t)");
            }
        }
        else
        {
            compact = true;
            int status = allocate(lda);
            if (status != 0)
            {
                failbit = true;
                throw std::string("Memory allocation error. matrix<T>::resize(size_t)");
            }
            allocated = true;
            compact = true;
        }

        resizedElements = sz;
        numRow = numCol = lda;
        rectangular = false;
        symetric = true;
    }

    //===============================================================================================================

    template <typename T>
    int matrix<T>::allocate(size_t lda)
    {
        /*
            Private member.
            Allocates memmory for a symmetric matrix in compact form.

            Return value: integer value; if 0 - allocation is successfull, otherwise - error.
        */

        size_t sz = static_cast<size_t>((lda * lda + lda) / 2);
#ifdef intelmkl
        A = (T *)mkl_malloc(sz * sizeof(T), sizeof(T) * 8);
#else
        A = (T *)malloc(sz * sizeof(T));
#endif
        if (A == NULL)
        {
#ifdef intelmkl
            mkl_free(A);
#else
            free(A);
#endif
            allocated = false;
            failbit = true;
            throw std::string("Memory allocation error. matrix<T>::allocate(size_t)");
        }
        for (size_t i = 0; i < sz; i++)
        {
            A[i] = 0.0;
        }

        return 0;
    }

    //===============================================================================================================

    template <typename T>
    int matrix<T>::allocate(size_t row, size_t col)
    {
        /*
            Private member.
            Allocates memmory for an arbitrary rectangular matrix.

            Return value: integer value; if 0 - allocation is successfull, otherwise - error.
        */

#ifdef intelmkl
        A = (T *)mkl_malloc(row * col * sizeof(T), sizeof(T) * 8);
#else
        // A = (T *)aligned_alloc(sizeof(T) * 8, row * col * sizeof(T));
        A = (T *)malloc(row * col * sizeof(T));
#endif
        if (A == NULL)
        {
#ifdef intelmkl
            mkl_free(A);
#else
            free(A);
#endif
            allocated = false;
            failbit = true;
            throw std::string("Memory allocation error. matrix<T>::allocate(size_t, size_t)");
        }
        for (size_t i = 0; i < (row * col); i++)
        {
            A[i] = 0.0;
        }

        return 0;
    }

    //===============================================================================================================

    template <typename T>
    matrix<T>::matrix(size_t row, size_t col)
    {
        /*
           Constructor for regular rectangular matrix
        */

        failbit = false;
        failinfo = 0;

        std::random_device rd;
        srand(rd());

        int iNum = rand() % 100000;
        binFilename = "matrix_" + std::to_string(iNum);

        rectangular = true;
        compact = false;
        symetric = false;
        allocated = true;
        ondisk = false;

        int status = allocate(row, col);
        if (status != 0)
        {
            failbit = true;
            allocated = false;
            failinfo = -1;
            throw std::string("Memory allocation error. matrix<T>::matrix(size_t, size_t)");
        }

        numRow = row;
        numCol = col;
        resizedElements = row * col;
    }

    //===============================================================================================================

    template <typename T>
    matrix<T>::matrix(size_t lda)
    {
        /*
           Constructor for symmetrical matrix in a compact form
        */

        failbit = false;
        failinfo = 0;

        std::random_device rd;
        srand(rd());

        int iNum = rand() % 100000;
        binFilename = "matrix_" + std::to_string(iNum);

        rectangular = false;
        symetric = true;
        allocated = true;
        compact = true;
        ondisk = false;

        int status = allocate(lda);
        if (status != 0)
        {
            failbit = true;
            allocated = false;
            failinfo = -1;
            throw std::string("Memory allocation error. matrix<T>::matrix(size_t)");
        }

        numRow = numCol = lda;
        resizedElements = (lda * lda + lda) / 2;
    }

    //===============================================================================================================

    template <typename T>
    matrix<T>::matrix()
    {
        /*
           Default constructor for rectangular matrix;
           no memmory allocation.
        */

        failbit = false;
        failinfo = 0;

        std::random_device rd;
        srand(rd());
        int iNum = rand() % 100000;

        binFilename = "matrix_" + std::to_string(iNum);

        allocated = false;
        rectangular = true;
        symetric = false;
        compact = false;
        ondisk = false;
        numRow = numCol = 0;
        resizedElements = 0;
    }

    //===============================================================================================================

    template <typename T>
    matrix<T>::matrix(const char *type)
    {
        /*
           Constructor: explicite determination of a matrix type;
           no memmory allocation.
        */

        failbit = false;
        failinfo = 0;

        std::random_device rd;
        srand(rd());

        int iNum = rand() % 100000;
        binFilename = "matrix_" + std::to_string(iNum);

        allocated = false;
        ondisk = false;

        if (*type == 's' || 'S')
        {
            rectangular = false;
            symetric = true;
            compact = true;
        }
        else if (*type == 'r' || 'R')
        {
            rectangular = true;
            symetric = false;
            compact = false;
        }
        numRow = numCol = 0;
        resizedElements = 0;
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::scale(T val)
    {
        /*
            Matrix scaling by scalar.

            Return value: none.

            Example:

                matrix <double> M(3);
                double val;
                for(auto i = 0; i < M.size(); i++)
                    M[i] = 2.0;

                M.scale(0.5);
                val = M[ 0 ]; // val = 1.0
                val = M[ 1 ]; // val = 1.0
                ...
                val = M[ M.size()-1 ]; // val = 1.0

        */

        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::scale(T)");

        if (allocated)
        {

            auto n_threads = std::thread::hardware_concurrency();
            auto block_size = static_cast<unsigned int>(size() / (n_threads));

            if (block_size < worksize)
            {
                block_size = static_cast<unsigned int>(size());
                n_threads = 1;
            }

            // #pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < size(); i++)
                A[i] = A[i] * val;
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::symtorec()
    {
        /*
            Trnsforms symmetric matrix in compact form to regular rectangular matrix (not compact).

            Return value: none.

            Example:
                    matrix <mytype> M(3); // symmetric matrix in compact form
                    matrix <mytype> res; // empty matrix
                    M.symtorec(); // M now is regular square (symmetric) matrix
                    res = M; // res now is regular square (symmetric) matrix
        */

        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::symtorec()");

        if (allocated && compact)
        {

            matrix<T> C(numRow, numRow);

            /*auto n_threads = std::thread::hardware_concurrency();
            auto block_size = static_cast<unsigned int>(size() / (n_threads));

            if (block_size < worksize)
            {
                block_size = size();
                n_threads = 1;
            }

    //#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < numRow; i++)
            {
                for (size_t j = 0; j <= i; j++)
                {
                    C.A[i * numRow + j] = C.A[j * numRow + i] = A[i * (i + 1) / 2 + j];
                }
            }*/
            // #pragma omp parallel sections
            {
                // #pragma omp section
                {
                    for (size_t i = 0; i < numRow; i++)
                    {
                        for (size_t j = 0; j <= i; j++)
                            C.A[i * numRow + j] = A[i * (i + 1) / 2 + j];
                    }
                }
                // #pragma omp section
                {
                    for (size_t i = 0; i < numRow; i++)
                    {
                        for (size_t j = 0; j < i; j++)
                            C.A[j * numRow + i] = A[i * (i + 1) / 2 + j];
                    }
                }
            }

            compact = false;
            rectangular = true;
            resize(numRow, numRow);
            fclear();
            *this = C;
        }
        else
        {
            failbit = true;
            if (!allocated)
                throw std::string("Matrix is empty. matrix<T>::symtorec()");
            else
                throw std::string("Matrix is not in a compact form. matrix<T>::symtorec()");
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::rectosym()
    {
        /*
            Trnsforms regular rectangular (symmetric) matrix into compact form.

            Return value: none.

            Example:
                    matrix <mytype> M(3); // symmetric matrix in compact form
                    matrix <mytype> res; // empty matrix
                    M.symtorec();
                    res = M; // res now is regular square (symmetric) matrix
                    res.rectosym(); // res now is symmetric matrix in a compact form
        */

        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::rectosym()");

        if (allocated && !compact)
        {
            matrix<T> C(numCol);

            auto n_threads = std::thread::hardware_concurrency();
            auto block_size = static_cast<unsigned int>(size() / (n_threads));

            if (block_size < worksize)
            {
                block_size = static_cast<unsigned int>(size());
                n_threads = 1;
            }

            // #pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < numCol; i++)
            {
                for (size_t j = 0; j <= i; j++)
                {
                    C.A[i * (i + 1) / 2 + j] = A[i * numCol + j];
                }
            }

            compact = true;
            rectangular = false;
            resize(numCol);
            fclear();
            *this = C;
        }
        else
        {
            failbit = true;
            if (!allocated)
                throw std::string("Matrix is empty. matrix<T>::rectosym()");
            else
                std::cout<<"Warning matrix<T>::rectosym(): Matrix is already in a compact form."<<'\n';
        }
    }

    //===============================================================================================================

    template <typename T>
    T &matrix<T>::at(size_t atRow, size_t atCol)
    {
        /*
            Safe access to the elements of a matrix.
            Indexes atRow and atCol caunted starting from 0.

            Return value: the element at the specified position in the container.
            Example:
                     matrix <double> M; // M is empty matrix.
                     M.at(10,10) = 0.8;          // M now is (11,11) matrix.
                     double val = M.at(10,10);   // val = 0.8.
        */

        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::at(size_t, size_t)");

        size_t index = 0;

        if (!compact)
        {

            if ((atRow + 1 > numRow) || (atCol + 1 > numCol))
            {
                numRow = _max(numRow, atRow + 1);
                numCol = _max(numCol, atCol + 1);
                if (allocated)
                {
                    resize();
                }
                else
                {
                    int status = allocate(numRow, numCol);
                    if (status != 0)
                    {
                        failbit = true;
                        throw std::string("Memory allocation error. matrix<T>::at(size_t, size_t)");
                    }
                    allocated = true;
                    resizedElements = numRow * numCol;
                }
            }
            index = atRow * numCol + atCol;
        }
        else
        {
            if (atRow < atCol)
            {
                failbit = true;
                throw std::string("For a matrix in a compact form 'columns > rows' is not allowed. matrix<T>::at(size_t, size_t)");
            }
            if (atRow + 1 > numRow)
            {
                numRow = numCol = atRow + 1;
                if (allocated)
                {
                    resize();
                }
                else
                {
                    int status = allocate(numRow);
                    if (status != 0)
                    {
                        failbit = true;
                        throw std::string("Memory allocation error. matrix<T>::at(size_t, size_t)");
                    }
                    allocated = true;
                    resizedElements = (numRow * numRow + numRow) / 2;
                }
            }
            index = atRow * (atRow + 1) / 2 + atCol;
        }

        return A[index];
    }

    //===============================================================================================================

    template <typename T>
    size_t matrix<T>::size() const
    {
        /*
            Returns number of elements allocated for A.

            Return value: size_t value, number of elements in a matrix.
            Example:
                matrix <double> M(3); // symmetric matrix in a compact form.
                double val = M.size(); // val = 6.

        */

        if (!compact)
            return numRow * numCol;
        else
            return (numCol * numCol + numCol) / 2;
    }
    //===============================================================================================================

    template <typename T>
    matrix<size_t> matrix<T>::shape() const
    {
        /*
            Returns the number of rows and columns of A.
        */

        matrix<size_t> shp(2, 1);
        shp[0] = numRow;
        shp[1] = numCol;

        return shp;
    }

    //===============================================================================================================

    template <typename T>
    matrix<T> matrix<T>::operator+(const matrix<T> &rhs)
    {
        /*
            Matrix addition.
            Accepts a general rectangular matrix as well as symmetric matrix in a compact form.

            Return value: matrix.
            Example:
                     matrix <double> M1(2,3);
                     matrix <double> M2(2,3);
                     matrix <double> res;
                     res = M1 + M2;
        */

        if (ondisk || rhs.ondisk)
            throw std::string("One of the matrices is empty. Use fread() to relocate data to memory. matrix<T>::operator+");

        matrix<T> C;

        if (!eq(rhs))
        {
            C.failbit = true;
            throw std::string("Matrices dimensions are not equal. matrix<T>::operator+");
        }

        if (!rhs.compact)
        {
            int status = C.allocate(rhs.numRow, rhs.numCol);
            if (status != 0)
            {
                C.failbit = true;
                throw std::string("Memory allocation error. matrix<T>::operator+");
            }

            C.allocated = true;
            resizedElements = rhs.numRow * rhs.numCol;
            C.numCol = rhs.numCol;
            C.numRow = rhs.numRow;
        }
        else
        {
            int status = C.allocate(rhs.numRow);
            if (status != 0)
            {
                C.failbit = true;
                throw std::string("Memory allocation error. matrix<T>::operator+");
            }

            C.allocated = true;
            resizedElements = (rhs.numRow * rhs.numRow + rhs.numRow) / 2;
            C.numCol = rhs.numRow;
            C.numRow = rhs.numRow;
        }

        auto n_threads = std::thread::hardware_concurrency();
        auto block_size = static_cast<unsigned int>(C.size() / (n_threads));

        if (block_size < worksize)
        {
            block_size = static_cast<unsigned int>(C.size());
            n_threads = 1;
        }

        // #pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
        for (size_t i = 0; i < C.size(); i++)
            C.A[i] = A[i] + rhs.A[i];

        return matrix(C);
    }

    //===============================================================================================================

    template <typename T>
    matrix<T> matrix<T>::operator<<(const matrix<T> &rhs)
    {
        /*
            Matrix concatenation, 'horizontally'.
            Accepts a general rectangular matrix as well as symmetric matrix in a compact form.

            Return value: matrix.
            Example:
                     matrix <double> M1(2,3);
                     matrix <double> M2(2,3);
                     matrix <double> res;
                     res = M1 << M2;

                     return: res(2,6).
        */

        if (ondisk || rhs.ondisk)
            throw std::string("One of the matrices is empty. Use fread() to relocate data to memory. matrix<T>::operator<<");

        matrix<T> C;

        // checking: the rows should be equal
        if ((rhs.numRow != numRow))
        {
            C.failbit = true;
            throw std::string("Matrices dimensions (rows) are not equal. matrix<T>::operator<<");
        }

        if ((rhs.compact && !compact) || (!rhs.compact && compact))
        {
            C.failbit = true;
            throw std::string("Matrices types are not the same. matrix<T>::operator<<");
        }

        // allocating memory
        int status = C.allocate(rhs.numRow, rhs.numCol + numCol);
        if (status != 0)
        {
            C.failbit = true;
            throw std::string("Memory allocation error. matrix<T>::operator<<");
        }

        C.allocated = true;
        resizedElements = rhs.numRow * (rhs.numCol + numCol);
        C.numCol = rhs.numCol + numCol;
        C.numRow = rhs.numRow;

        if (rhs.compact)
        {
            // #pragma omp parallel sections
            {
                // #pragma omp section
                {
                    for (size_t i = 0; i < C.numRow; i++)
                    {
                        for (size_t j = 0; j <= i; j++)
                            C(i, j) = A[i * (i + 1) / 2 + j];
                    }
                }
                // #pragma omp section
                {
                    for (size_t i = 0; i < C.numRow; i++)
                    {
                        for (size_t j = 0; j < i; j++)
                            C(j, i) = A[i * (i + 1) / 2 + j];
                    }
                }
                // #pragma omp section
                {
                    for (size_t i = 0; i < C.numRow; i++)
                    {
                        for (size_t j = numCol; j <= numCol + i; j++)
                            C(i, j) = rhs.A[i * (i + 1) / 2 + (j - numCol)];
                    }
                }
                // #pragma omp section
                {
                    for (size_t i = 0; i < C.numRow; i++)
                    {
                        for (size_t j = numCol; j < numCol + i; j++)
                            C(j - numCol, i + numCol) = rhs.A[i * (i + 1) / 2 + (j - numCol)];
                    }
                }
            }
        }
        else
        {
            // #pragma omp parallel sections
            {
                // #pragma omp section
                {
                    for (size_t i = 0; i < C.numRow; i++)
                    {
                        for (size_t j = 0; j < numCol; j++)
                            C(i, j) = A[i * numCol + j];
                    }
                }
                // #pragma omp section
                {
                    for (size_t i = 0; i < C.numRow; i++)
                    {
                        for (size_t j = numCol; j < numCol + rhs.numCol; j++)
                            C(i, j) = rhs.A[i * rhs.numCol + (j - numCol)];
                    }
                }
            }

        } // end of if-else section

        return matrix(C);
    }

    //===============================================================================================================

    template <typename T>
    matrix<T> matrix<T>::operator>>(const matrix<T> &rhs)
    {
        /*
            Matrix concatenation, 'vertically'.
            Accepts a general rectangular matrix as well as symmetric matrix in a compact form.

            Return value: matrix.
            Example:
                     matrix <double> M1(2,3);
                     matrix <double> M2(5,3);
                     matrix <double> res;
                     res = M1 >> M2;

                     return: res(7,3).
        */

        if (ondisk || rhs.ondisk)
            throw std::string("One of the matrices is empty. Use fread() to relocate data to memory. matrix<T>::operator>>");

        matrix<T> C;

        // checking: the cols should be equal
        if ((rhs.numCol != numCol))
        {
            C.failbit = true;
            throw std::string("Matrices dimensions (columns) are not equal. matrix<T>::operator>>");
        }

        if ((rhs.compact && !compact) || (!rhs.compact && compact))
        {
            C.failbit = true;
            throw std::string("Matrices types are not the same. matrix<T>::operator>>");
        }

        // allocating memory
        int status = C.allocate(rhs.numRow + numRow, rhs.numCol);
        if (status != 0)
        {
            C.failbit = true;
            throw std::string("Memory allocation error. matrix<T>::operator>>");
        }

        C.allocated = true;
        resizedElements = (rhs.numRow + numRow) * rhs.numCol;
        C.numCol = rhs.numCol;
        C.numRow = rhs.numRow + numRow;

        if (rhs.compact)
        {
            // #pragma omp parallel sections
            {
                // #pragma omp section
                {
                    for (size_t i = 0; i < numRow; i++)
                    {
                        for (size_t j = 0; j <= i; j++)
                            C(i, j) = A[i * (i + 1) / 2 + j];
                    }
                }
                // #pragma omp section
                {
                    for (size_t i = 0; i < numRow; i++)
                    {
                        for (size_t j = 0; j < i; j++)
                            C(j, i) = A[i * (i + 1) / 2 + j];
                    }
                }
                // #pragma omp section
                {
                    for (size_t i = numRow; i < C.numRow; i++)
                    {
                        for (size_t j = 0; j <= (i - numRow); j++)
                            C(i, j) = rhs.A[(i - numRow) * ((i - numRow) + 1) / 2 + j];
                    }
                }
                // #pragma omp section
                {
                    for (size_t i = numRow; i < C.numRow; i++)
                    {
                        for (size_t j = 0; j <= (i - numRow); j++)
                            C(j + numRow, i - numRow) = rhs.A[(i - numRow) * ((i - numRow) + 1) / 2 + j];
                    }
                }
            }
        }
        else
        {
            // #pragma omp parallel sections
            {
                // #pragma omp section
                {
                    for (size_t i = 0; i < numRow; i++)
                    {
                        for (size_t j = 0; j < numCol; j++)
                            C(i, j) = A[i * numCol + j];
                    }
                }
                // #pragma omp section
                {
                    for (size_t i = numRow; i < C.numRow; i++)
                    {
                        for (size_t j = 0; j < numCol; j++)
                            C(i, j) = rhs.A[(i - numRow) * numCol + j];
                    }
                }
            }

        } // end of if-else section

        return matrix(C);
    }

    //===============================================================================================================

    template <typename T>
    matrix<T> matrix<T>::operator-(const matrix<T> &rhs)
    {
        /*
            Matrix substitution.
            Accepts a general rectangular matrix as well as symmetric matrix in a compact form.

            Return value: matrix.
            Example:
                     matrix <double> M1(2,3);
                     matrix <double> M2(2,3);
                     matrix <double> res;
                     res = M1 - M2;
        */

        if (ondisk || rhs.ondisk)
            throw std::string("One of the matrices is empty. Use fread() to relocate data to memory. matrix<T>::operator-");

        matrix<T> C;

        if (!eq(rhs))
        {
            C.failbit = true;
            throw std::string("Matrices dimensions are not equal. matrix<T>::operator-");
        }

        if (!rhs.compact)
        {
            int status = C.allocate(rhs.numRow, rhs.numCol);
            if (status != 0)
            {
                C.failbit = true;
                throw std::string("Memory allocation error. matrix<T>::operator-");
            }

            C.allocated = true;
            resizedElements = rhs.numRow * rhs.numCol;
            C.numCol = rhs.numCol;
            C.numRow = rhs.numRow;
        }
        else
        {
            int status = C.allocate(rhs.numRow);
            if (status != 0)
            {
                C.failbit = true;
                throw std::string("Memory allocation error. matrix<T>::operator-");
            }

            C.allocated = true;
            resizedElements = (rhs.numRow * rhs.numRow + rhs.numRow) / 2;
            C.numCol = rhs.numRow;
            C.numRow = rhs.numRow;
        }

        auto n_threads = std::thread::hardware_concurrency();
        auto block_size = static_cast<unsigned int>(C.size() / (n_threads));

        if (block_size < worksize)
        {
            block_size = static_cast<unsigned int>(C.size());
            n_threads = 1;
        }

        // #pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
        for (size_t i = 0; i < C.size(); i++)
            C.A[i] = A[i] - rhs.A[i];

        return matrix(C);
    }

    //===============================================================================================================

    template <typename T>
    matrix<T> matrix<T>::operator-(const T val)
    {
        /*
            Overloaded '-' operator to substract a scalar from a matrix.

            Return value: matrix object.
            Example:

                matrix <double> res;    // empty matrix
                matrix <double> M(n,m); // M is (n,m) matrix initialized by 0.

                res = M - 2.0; // res become (n,m) matrix where all elements are -2.0;
                               // M remains unchanged
        */

        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::operator-");

        matrix<T> C;

        if (!compact)
        {
            int status = C.allocate(numRow, numCol);
            if (status != 0)
            {
                C.failbit = true;
                throw std::string("Memory allocation error. matrix<T>::operator-");
            }

            C.allocated = true;
            resizedElements = numRow * numCol;
            C.numCol = numCol;
            C.numRow = numRow;
        }
        else
        {
            int status = C.allocate(numRow);
            if (status != 0)
            {
                C.failbit = true;
                throw std::string("Memory allocation error. matrix<T>::operator-");
            }

            C.allocated = true;
            resizedElements = (numRow * numRow + numRow) / 2;
            C.numCol = numRow;
            C.numRow = numRow;
        }

        auto n_threads = std::thread::hardware_concurrency();
        auto block_size = static_cast<unsigned int>(C.size() / (n_threads));

        if (block_size < worksize)
        {
            block_size = static_cast<unsigned int>(C.size());
            n_threads = 1;
        }

        // #pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
        for (size_t i = 0; i < C.size(); i++)
            C.A[i] = A[i] - val;

        return matrix(C);
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::gemmt_intrf(double *_A, double *B, MKL_INT rowA, MKL_INT colA, MKL_INT colB)
    {
        /*
            Interface to cblas_dgemmt routine which computes a matrix-matrix product with general matrices
            but updates only the upper or lower triangular part of the result matrix.
        */
#ifdef intelmkl
        cblas_dgemmt(CblasRowMajor, CblasLower, CblasNoTrans, CblasTrans, rowA, colA, 1.0, _A, colA, _A, colA, 0.0, B, colB);
#else
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, rowA, colB, colA, 1.0, _A, colA, _A, colA, 0.0, B, colB);
#endif
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::gemmt_intrf(float *_A, float *B, MKL_INT rowA, MKL_INT colA, MKL_INT colB)
    {
        /*
            Interface to cblas_sgemmt routine which computes a matrix-matrix product with general matrices
            but updates only the upper or lower triangular part of the result matrix.
        */
#ifdef intelmkl
        cblas_sgemmt(CblasRowMajor, CblasLower, CblasNoTrans, CblasTrans, rowA, colA, 1.0, _A, colA, _A, colA, 0.0, B, colB);
#else
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, rowA, colB, colA, 1.0, _A, colA, _A, colA, 0.0, B, colB);
#endif
    }
    //===============================================================================================================

    template <typename T>
    matrix<T> matrix<T>::operator^(const int val)
    {
        /*
            Overloaded matrix power operator which computes:

            1)  matrix product: A * A' if val == 2.
                Return value: rectangular matrix.

            2)  inverse of matrix A: A^-1 if val == -1.
                Return value: (a) general rectangular matrix if A is not symmetric;
                              (b) matrix in compact form if A is symmetric.

            3)  inversion of matrix product: (A * A')^-1 if val == -2.
                Return value: rectangular matrix.

            Original matrix A remains unchanged.

            Example:
                matrix <double> M(n,n);
                matrix <double> res;

                for (auto i = 0; i < M.size(); i++){
                    M[i] = static_cast <double> (i);

                res = M^(-2);
        */

        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::operator^");

        matrix<T> C(*this);

        if (val == 2)
        {

            matrix<T> tmpM;

            if (C.compact)
                C.symtorec();

            int status = tmpM.allocate(numRow, numRow);
            if (status != 0)
            {
                C.failbit = true;
                throw std::string("Memory allocation error. matrix<T>::operator^");
            }

            tmpM.allocated = true;
            tmpM.compact = false;
            tmpM.numRow = numRow;
            tmpM.numCol = numRow;
            tmpM.resizedElements = numRow * numRow;

            gemmt_intrf(C.A, tmpM.A, C.numRow, C.numCol, tmpM.numCol);
            
            /*
                Because C is symmetrical, it consists of values only for 'L' part, 'U' part are zeros.
                Therefore we restore and further return regular rectangular matrix.
            */
            for (size_t i = 0; i < tmpM.numRow; i++)
            {
                for (size_t j = 0; j <= i; j++)
                {
                    tmpM(j, i) = tmpM(i, j);
                }
            }

            C = tmpM;
        }
        else if (val == -1)
        {

            if (numRow == numCol)
            {

                if (!C.compact)
                {
                    try
                    {
                        inv_rec(C.A, numRow, numCol);
                    }
                    catch (std::string err)
                    {
                        C.failbit = true;
                        throw err;
                    }
                }
                else
                {
                    try
                    {
                        inv_sym(C.A, numCol);
                    }
                    catch (std::string err)
                    {
                        C.failbit = true;
                        throw err;
                    }
                }
            }
            else
            {
                C.failbit = true;
                throw std::string("The matrix is not square. matrix<T>::operator^");
            }
        }
        else if (val == -2)
        {

            /* First do A*A' */

            matrix<T> tmpM;

            if (C.compact)
                C.symtorec();

            int status = tmpM.allocate(numRow, numRow);
            if (status != 0)
            {
                C.failbit = true;
                throw std::string("Memory allocation error. matrix<T>::operator^");
            }

            tmpM.allocated = true;
            tmpM.compact = false;
            tmpM.numRow = numRow;
            tmpM.numCol = numRow;
            tmpM.resizedElements = numRow * numRow;

            gemmt_intrf(C.A, tmpM.A, C.numRow, C.numCol, tmpM.numCol);

            /*
                Because C is symmetrical, it consists of values only for 'L' part, 'U' part are zeros.
                Therefore we restore and further return regular rectangular matrix.
            */
            for (size_t i = 0; i < tmpM.numRow; i++)
            {
                for (size_t j = 0; j <= i; j++)
                {
                    tmpM(j, i) = tmpM(i, j);
                }
            }

            C = tmpM;
            tmpM.clear();

            /* Then do inv(A*A') */

            try
            {
                inv_rec(C.A, C.numRow, C.numCol);
            }
            catch (std::string err)
            {
                C.failbit = true;
                throw err;
            }
        }
        else
        {
            C.failbit = true;
            throw std::string("Not supported operation. matrix<T>::operator^");
        }

        return matrix(C);
    }

    //===============================================================================================================

    template <typename T>
    matrix<T> matrix<T>::operator^(const char *val)
    {
        /*
            Overloaded matrix power operator which computes matrix transpose.
            Original matrix A remains unchanged.

            Return value: rectangular matrix.

            Example:
                matrix <double> M(n,m);
                matrix <double> res;

                res = M^"T"; // res now is (m,n) matrix while M is (n,m) matrix.
        */

        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::operator^");

        matrix<T> C(*this);
        if (*val == 'T' || 't')
        {
            C.transpose();
        }
        else
        {
            C.failbit = true;
            throw "Invalid operation. matrix<T>::operator^";
        }
        return matrix(C);
    }

    //===============================================================================================================

    template <typename T>
    matrix<T> matrix<T>::operator*(const matrix<T> &rhs)
    {
        /*
            Matrix dot product:
            C = A * B;
            where A & B are rectangular matrices, C is always rectangular.

            Return value: rectangular matrix.
        */

        if (ondisk || rhs.ondisk)
            throw std::string("One of the matrices is empty. Use fread() to relocate data to memory. matrix<T>::operator*");

        matrix<T> C;

        /* Check if matrices can be multiplied. */
        if (numCol != rhs.numRow)
        {
            C.failbit = true;
            throw std::string("Matrices are not consistent for multiplication. matrix<T>::operator*");
        }

        int status = C.allocate(numRow, rhs.numCol);
        if (status != 0)
        {
            C.failbit = true;
            throw std::string("Memory allocation error. matrix<T>::operator*");
        }

        C.allocated = true;
        C.numRow = numRow;
        C.numCol = rhs.numCol;
        C.resizedElements = C.numRow * C.numCol;

        dotprod(A, rhs.A, C.A, numRow/*, rhs.numRow*/, numCol, rhs.numCol);

        return C;//matrix(C);
    }

    //===============================================================================================================

    template <typename T>
    matrix<T> matrix<T>::operator*(const T &val)
    {
        /*
            Overloaded '*' operator to multiply matrix by scalar.

            Return value: matrix object.

            Example:

                matrix <double> obj;    // empty matrix
                matrix <double> M(n,m); // M is (n,m) matrix initialized by 0.

                for (auto i = 0; i < M.size(); i++)
                    M[i] = 1.0;

                obj = M * (-2.0); // obj become (n,m) matrix where all elements are -2.0;
                                  // M remains unchanged (all elements of the matrix are 1.0)
        */

        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::operator*");

        matrix<T> C;
        if (!compact)
            C.resize(numRow, numCol);
        else
            C.resize(numRow);

        // #pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
        for (size_t i = 0; i < C.size(); i++)
            C.A[i] = A[i] * val;

        return matrix(C);
    }

    //===============================================================================================================

    template <typename T>
    matrix<T> matrix<T>::operator+(const T &val)
    {
        /*
            Overloaded '+' operator to add scalar to matrix.

            Return value: matrix object.

            Example:

                matrix <double> obj;    // empty matrix
                matrix <double> M(n,m); // M is (n,m) matrix initialized by 0.

                for (auto i = 0; i < M.size(); i++)
                    M[i] = 1.0;

                obj = M + (-2.0); // obj become (n,m) matrix where all elements are -2.0;
                                  // M remains unchanged (all elements of the matrix are 1.0)
        */

        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::operator*");

        matrix<T> C;
        if (!compact)
        {
            int status = C.allocate(numRow, numCol);
            if (status != 0)
            {
                C.failbit = true;
                throw std::string("Memory allocation error: matrix<T>::operator*");
            }

            C.allocated = true;
            C.resizedElements = numRow * numCol;
            C.numCol = numCol;
            C.numRow = numRow;
        }
        else
        {
            int status = C.allocate(numRow);
            if (status != 0)
            {
                C.failbit = true;
                throw std::string("Memory allocation error: matrix<T>::operator*");
            }

            C.allocated = true;
            C.resizedElements = (numRow * numRow + numRow) / 2;
            C.numCol = numRow;
            C.numRow = numRow;
        }

        auto n_threads = std::thread::hardware_concurrency();
        auto block_size = static_cast<unsigned int>(C.size() / (n_threads));

        if (block_size < worksize)
        {
            block_size = static_cast<unsigned int>(C.size());
            n_threads = 1;
        }

        // #pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
        for (size_t i = 0; i < C.size(); i++)
            C.A[i] = A[i] + val;

        return matrix(C);
    }

    //===============================================================================================================

    template <typename T>
    bool matrix<T>::operator==(const matrix<T> &rhs) const
    {
        /*
            Checks matrices complete equality: A == B.

            Return value: logical TRUE if two matrices are equal and FALSE otherwise.
        */

        if (ondisk || rhs.ondisk)
            throw std::string("One of the matrices is empty. Use fread() to relocate data to memory. matrix<T>::operator==");

        if ((rhs.numRow != numRow) || (rhs.numCol != numCol))
            return false;

        if ((rhs.compact && !compact) || (!rhs.compact && compact))
            return false;

        for (size_t i = 0; i < size(); i++)
            if ( A[i] != rhs.A[i] )
                return false;
        
        return true;
        //return std::equal(A, A + sizeof A / sizeof *A, rhs.A);
    }

    //===============================================================================================================

    template <typename T>
    bool matrix<T>::eq(const matrix<T> &rhs)
    {
        /*
            Checks matrices shapes equality: A == B.

            Return value: logical TRUE if shapes of two matrices are equal and FALSE otherwise.
                          Values of matrices elements are not compared.
        */

        if ((rhs.numRow != numRow) || (rhs.numCol != numCol))
            return false;

        if ((rhs.compact && !compact) || (!rhs.compact && compact))
            return false;

        return true;
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::dotprod(double *_A, double *B, double *C, MKL_INT rowA/*, MKL_INT rowB*/, MKL_INT colA, MKL_INT colB)
    {

        /*
         * Matrix product (rectangular matrices):
         * C = A * B;
         *
         * void cblas_dgemm (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa,...
         * 					const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n,...
         * 					const MKL_INT k, const double alpha, const double *a,...
         * 					const MKL_INT lda, const double *b, const MKL_INT ldb,...
         * 					const double beta, double *c, const MKL_INT ldc);
         */

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rowA, colB, colA, 1.0, _A, colA, B, colB, 0.0, C, colB);
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::dotprod(float *_A, float *B, float *C, MKL_INT rowA/*, MKL_INT rowB*/, MKL_INT colA, MKL_INT colB)
    {

        /*
         * Matrix product (rectangular matrices):
         * C = A * B;
         *
         * void cblas_sgemm (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa,...
         * 					const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n,...
         * 					const MKL_INT k, const float alpha, const float *a,...
         * 					const MKL_INT lda, const float *b, const MKL_INT ldb,...
         * 					const float beta, float *c, const MKL_INT ldc);
         */

        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rowA, colB, colA, 1.0, _A, colA, B, colB, 0.0, C, colB);
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::clear()
    {
        /*
            Clears a memory allocated for the container A.

            Return value: none.
        */

        if (allocated)
        {
#ifdef intelmkl
            mkl_free(A);
#else
            free(A);
#endif
            allocated = false;
            failbit = false;
            failinfo = 0;
            resizedElements = 0;
            numRow = numCol = 0;
        }

        if ( is_allocated_ipv_A )
        {
#ifdef intelmkl
            mkl_free(ipiv_A);
#else
            free(ipiv_A);
#endif
            is_allocated_ipv_A = false;
        }

        // rename file name:
        binFilename.clear();
        std::random_device rd;
        srand(rd());
        int iNum = rand() % 100000;
        binFilename = "matrix_" + std::to_string(iNum);

        ondisk = false;
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::fclear()
    {
        /*
            Removes a binary file (if exists) associated with A.
            Note: this shouldbe called before clear() becaause
            the last call changes the name of associated binary file!

            Return value: none.
        */

        std::ifstream f(binFilename.c_str());
        if (f.good())
            remove(binFilename.c_str());
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::fclear(const std::string &fname)
    {
        /*
            Removes a binary file (if exists) associated with A.
            Note: this shouldbe called before clear() becaause
            the last call changes the name of associated binary file!

            Return value: none.
        */

        std::ifstream f(fname.c_str());
        if (f.good())
            remove(fname.c_str());
    }

    //===============================================================================================================

    template <typename T>
    bool matrix<T>::empty()
    {
        /*
            Checks if a matrix is empty.

            Return value: logical TRUE if matrix is empty, otherwise returns FALSE.
        */

        if (allocated)
            return false;

        return true;
    }

    //===============================================================================================================

    template <typename T>
    bool matrix<T>::is_symmetric()
    {
        /*
            Checks if a matrix is symmetric.

            Return value: logical TRUE if matrix is symmetric, otherwise returns FALSE.
        */

        if (symetric)
            return true;
        else
            return false;
    }

    //===============================================================================================================

    template <typename T>
    matrix<T>::matrix(const matrix<T> &obj)
    {
        /*
            Copy constructor.
        */
        compact = obj.compact;
        rectangular = obj.rectangular;
        symetric = obj.symetric;
        failbit = obj.failbit;
        failinfo = obj.failinfo;
        numCol = obj.numCol;
        numRow = obj.numRow;
        resizedElements = obj.resizedElements;
        binFilename = obj.binFilename;
        allocated = false;
        ondisk = obj.ondisk;

        if (compact)
            allocate(numRow);
        else
            allocate(numRow, numCol);
        
        allocated = true;

        if (obj.allocated)
        {
            std::memcpy(A, obj.A, sizeof(T) * this->size());

            if (obj.is_allocated_ipv_A)
            {
                ipiv_A = obj.ipiv_A;
                is_allocated_ipv_A = obj.is_allocated_ipv_A;
                obj.is_allocated_ipv_A = false;
            }
        }
        else
        {
#ifdef intelmkl
            mkl_free(A);
#else
            free(A);
#endif
            allocated = false;
            failbit = false;
            failinfo = 0;
            resizedElements = 0;
        }
    }

    //===============================================================================================================

    template <typename T>
    matrix<T> &matrix<T>::operator=(const matrix<T> &rhs)
    {
        /*
            Overloaded assignment operator.

            Procedure: (1) Copy rhs to tmpObj; (2) Swaps lhs with tmpObj.
        */

        matrix<T> tmpObj(rhs);

        std::swap(compact, tmpObj.compact);
        std::swap(ondisk, tmpObj.ondisk);
        std::swap(rectangular, tmpObj.rectangular);
        std::swap(symetric, tmpObj.symetric);
        std::swap(failbit, tmpObj.failbit);
        std::swap(failinfo, tmpObj.failinfo);
        std::swap(numCol, tmpObj.numCol);
        std::swap(numRow, tmpObj.numRow);
        std::swap(resizedElements, tmpObj.resizedElements);

        if ( ondisk ) // exchange file names only if rhs is not on memory!
            std::swap(binFilename, tmpObj.binFilename);

        std::swap(allocated, tmpObj.allocated);
        std::swap(A, tmpObj.A);
        
        std::swap(ipiv_A, tmpObj.ipiv_A);
        std::swap(is_allocated_ipv_A, tmpObj.is_allocated_ipv_A);

        return *this;
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::fwrite()
    {
        /*
            Moves matrix to a DISK and clears memory allocated for the container A.

            Return value: none.
        */

        fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        fA.open(binFilename, fA.binary | fA.trunc | fA.out);

        if (!fA.is_open())
        {
            failbit = true;
            throw std::string("Error while opening a binary file. matrix<T>::fwrite()");
        }

        size_t sz;
        if (!compact)
            sz = numRow * numCol;
        else
            sz = (numCol * numCol + numCol) / 2;

        fA.write(reinterpret_cast<char *>(A), sz * sizeof(T));

        fA.close();

        ondisk = true;

        //if (allocated)
        //{
#ifdef intelmkl
            mkl_free(A);
#else
            free(A);
#endif
            allocated = false;
        //}
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::fwrite(const std::string &fname)
    {
        /*
            Saves matrix to a DISK and not clears memory allocated for the container A.

            Return value: none.
        */
        fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        fA.open(fname, fA.binary | fA.trunc | fA.out);

        if (!fA.is_open())
        {
            failbit = true;
            throw std::string("Error while opening a binary file. matrix<T>::fwrite(const std::string &)");
        }

        size_t sz;
        if (!compact)
        {
            sz = numRow * numCol;
        }
        else
        {
            sz = (numCol * numCol + numCol) / 2;
        }

        size_t B[6];
        B[0] = matrix_type;
        B[1] = (size_t)compact;
        B[2] = numRow;
        B[3] = numCol;
        B[4] = (size_t)sizeof(T); // points to a type of matrix
        B[5] = 0;
        
        fA.write(reinterpret_cast<char *>(&B), 6 * sizeof(size_t));

        fA.write(reinterpret_cast<char *>(A), sz * sizeof(T));

        fA.close();
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::fwrite(const std::string &fname, std::vector<std::int64_t> &vect_data)
    {
        /*
            Saves matrix to a DISK and not clears memory allocated for the container A.

            Return value: none.
        */
        fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        fA.open(fname, fA.binary | fA.trunc | fA.out);

        if (!fA.is_open())
        {
            failbit = true;
            throw std::string("Error while opening a binary file. matrix<T>::fwrite(const std::string &, std::vector<std::int64_t> &)");
        }

        size_t sz;
        if (!compact)
        {
            sz = numRow * numCol;
        }
        else
        {
            sz = (numCol * numCol + numCol) / 2;
        }

        size_t B[6];
        B[0] = matrix_type;
        B[1] = (size_t)compact;
        B[2] = numRow;
        B[3] = numCol;
        B[4] = (size_t)sizeof(T); // points to a type of matrix
        B[5] = 0; // type of vect_data: 1 - string; 0 - int64_t

        size_t vect_data_size = vect_data.size();
        
        fA.write(reinterpret_cast<char *>(&B), 6 * sizeof(size_t));

        fA.write(reinterpret_cast<char *>(A), sz * sizeof(T));

        fA.write( reinterpret_cast<const char *>( &vect_data_size ), sizeof(size_t) ); // 5. writing size of vect_data
        fA.write( reinterpret_cast<const char *>( vect_data.data() ), vect_data_size * sizeof(std::int64_t) ); // 6. writing all vect_data

        fA.close();
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::fwrite(const std::string &fname, std::vector<std::string> &vect_data)
    {
        /*
            Saves matrix to a DISK and not clears memory allocated for the container A.

            Return value: none.
        */
        fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        fA.open(fname, fA.binary | fA.trunc | fA.out);

        if (!fA.is_open())
        {
            failbit = true;
            throw std::string("Error while opening a binary file. matrix<T>::fwrite(const std::string &, std::vector<std::int64_t> &)");
        }

        size_t sz;
        if (!compact)
        {
            sz = numRow * numCol;
        }
        else
        {
            sz = (numCol * numCol + numCol) / 2;
        }

        size_t B[6];
        B[0] = matrix_type;
        B[1] = (size_t)compact;
        B[2] = numRow;
        B[3] = numCol;
        B[4] = (size_t)sizeof(T); // points to a type of matrix
        B[5] = 1; // type of vect_data: 1 - string; 0 - int64_t
       
        size_t vect_data_size = vect_data.size();
        
        fA.write(reinterpret_cast<char *>(&B), 6 * sizeof(size_t));

        fA.write(reinterpret_cast<char *>(A), sz * sizeof(T));

        fA.write( reinterpret_cast<const char *>( &vect_data_size ), sizeof(size_t) ); // 5. writing size of vect_data
        for (size_t i = 0; i < vect_data_size; i++)
        {
            std::string str = vect_data[i];
            size_t str_size = str.size();
            fA.write( reinterpret_cast<const char *>( &str_size ), sizeof(size_t) ); // 5. writing size of vect_data
            fA.write( reinterpret_cast<const char *>( &str[0] ), str_size ); // 6. writing all vect_data
        }

        fA.close();
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::fread()
    {
        /*
            Moves saved on DISK matrix back to the memory.

            Return value: none.
        */
        if (!ondisk)
            throw std::string("The data was not moved to a binary file. matrix<T>::fread()");

        std::ifstream f(binFilename.c_str());
        if (f.good())
        {

            fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            fA.open(binFilename, fA.binary | fA.in);

            if (!fA.is_open())
            {
                failbit = true;
                throw std::string("Error while opening a binary file. matrix<T>::fread()");
            }

            size_t sz;
            if (!compact)
            {

                sz = numRow * numCol;

                int status = allocate(numRow, numCol);
                if (status != 0)
                    throw std::string("Memory allocation error. matrix<T>::fread()");

                allocated = true;
            }
            else
            {
                sz = (numCol * numCol + numCol) / 2;

                int status = allocate(numCol);
                if (status != 0)
                    throw std::string("Memory allocation error. matrix<T>::fread()");

                allocated = true;
            }

            fA.read(reinterpret_cast<char *>(A), sz * sizeof(T));

            fA.close();

            ondisk = false;
        }
        else
            throw std::string("matrix<T>::fread() => Cannot open the binary file to read from!");
    }
    //===============================================================================================================

    template <typename T>
    void matrix<T>::fread(const std::string &fname)
    {
        /*
            Moves saved on DISK matrix back to the memory.

            Return value: none.
        */

        std::ifstream f( fname.c_str() );
        if ( f.good() )
        {

            fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            fA.open(fname, fA.binary | fA.in);

            if (!fA.is_open())
            {
                failbit = true;
                throw std::string("Error while opening a binary file. matrix<T>::fread(const std::string &)");
            }

            size_t B[6];

            fA.read(reinterpret_cast<char *>(&B), 6 * sizeof(size_t));

            size_t matr_type = B[0];
            compact = (bool)B[1];
            numRow = B[2];
            numCol = B[3];
            size_t var_inbytes = B[4];
            
            if ( matr_type != matrix_type )
                throw std::string("The wrong kind of matrix trying to read the data from file! matrix<T>::fread(const std::string &)");

            if ( var_inbytes != sizeof(T) )
                throw std::string("The data type stored in file ("+ std::to_string(var_inbytes) + "-bytes data) is not consistent with type of the matrix trying to read the data. matrix<T>::fread(const std::string &)");

            size_t sz;
            if (!compact)
            {
                rectangular = true;
                symetric = false;

                sz = numRow * numCol;

                int status = allocate(numRow, numCol);
                if (status != 0)
                    throw std::string("Memory allocation error. matrix<T>::fread(const std::string &)");

                allocated = true;
            }
            else
            {
                rectangular = false;
                symetric = true;

                sz = (numCol * numCol + numCol) / 2;

                int status = allocate(numCol);
                if (status != 0)
                    throw std::string("Memory allocation error. matrix<T>::fread(const std::string &fname)");

                allocated = true;
            }

            fA.read(reinterpret_cast<char *>(A), sz * sizeof(T));

            fA.close();

            ondisk = false;

            resizedElements = sz;
        }
        else
            throw std::string("matrix<T>::fread(const std::string &) => Cannot open the binary file: " + fname + " to read from!");
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::fread(const std::string &fname, std::vector<std::int64_t> &vect_data)
    {
        /*
            Moves saved on DISK matrix back to the memory.

            Return value: none.
        */

        std::ifstream f( fname.c_str() );
        if ( f.good() )
        {

            fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            fA.open(fname, fA.binary | fA.in);

            if (!fA.is_open())
            {
                failbit = true;
                throw std::string("Error while opening a binary file. matrix<T>::fread(const std::string &, std::vector<std::int64_t> &)");
            }

            size_t B[6];

            fA.read(reinterpret_cast<char *>(&B), 6 * sizeof(size_t));

            size_t matr_type = B[0];
            compact = (bool)B[1];
            numRow = B[2];
            numCol = B[3];
            size_t var_inbytes = B[4];

            if ( matr_type != matrix_type )
                throw std::string("The wrong kind of matrix trying to read the data from file! matrix<T>::fread(const std::string &)");

            if ( var_inbytes != sizeof(T) )
                throw std::string("The data type stored in file ("+ std::to_string(var_inbytes) + "-bytes data) is not consistent with type of the matrix trying to read the data. matrix<T>::fread(const std::string &)");

            size_t sz;
            if (!compact)
            {
                rectangular = true;
                symetric = false;

                sz = numRow * numCol;

                int status = allocate(numRow, numCol);
                if (status != 0)
                    throw std::string("Memory allocation error. matrix<T>::fread(const std::string &, std::vector<std::string> &)");

                allocated = true;
            }
            else
            {
                rectangular = false;
                symetric = true;

                sz = (numCol * numCol + numCol) / 2;

                int status = allocate(numCol);
                if (status != 0)
                    throw std::string("Memory allocation error. matrix<T>::fread(const std::string &fname, std::vector<std::string> &)");

                allocated = true;
            }

            size_t vect_data_size;

            fA.read(reinterpret_cast<char *>(A), sz * sizeof(T));

            fA.read( reinterpret_cast<char *>( &vect_data_size ), sizeof(size_t) ); // 3. reading size of keys
            vect_data.resize(vect_data_size);
            fA.read( reinterpret_cast<char *>( vect_data.data() ), vect_data_size * sizeof(std::int64_t) ); // 4. reading all keys

            fA.close();

            ondisk = false;

            resizedElements = sz;
        }
        else
            throw std::string("matrix<T>::fread(const std::string &, std::vector<std::string> &) => Cannot open the binary file: " + fname +" to read from!");
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::fread(const std::string &fname, std::vector<std::string> &vect_data)
    {
        /*
            Moves saved on DISK matrix back to the memory.

            Return value: none.
        */

        std::ifstream f( fname.c_str() );
        if ( f.good() )
        {

            fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            fA.open(fname, fA.binary | fA.in);

            if (!fA.is_open())
            {
                failbit = true;
                throw std::string("Error while opening a binary file. matrix<T>::fread(const std::string &, std::vector<std::int64_t> &)");
            }

            size_t B[6];

            fA.read(reinterpret_cast<char *>(&B), 6 * sizeof(size_t));

            size_t matr_type = B[0];
            compact = (bool)B[1];
            numRow = B[2];
            numCol = B[3];
            size_t var_inbytes = B[4];
            
            if ( matr_type != matrix_type )
                throw std::string("The wrong kind of matrix trying to read the data from file! matrix<T>::fread(const std::string &)");

            if ( var_inbytes != sizeof(T) )
                throw std::string("The data type stored in file ("+ std::to_string(var_inbytes) + "-bytes data) is not consistent with type of the matrix trying to read the data. matrix<T>::fread(const std::string &)");

            size_t sz;
            if (!compact)
            {
                rectangular = true;
                symetric = false;

                sz = numRow * numCol;

                int status = allocate(numRow, numCol);
                if (status != 0)
                    throw std::string("Memory allocation error. matrix<T>::fread(const std::string &)");

                allocated = true;
            }
            else
            {
                rectangular = false;
                symetric = true;

                sz = (numCol * numCol + numCol) / 2;

                int status = allocate(numCol);
                if (status != 0)
                    throw std::string("Memory allocation error. matrix<T>::fread(const std::string &, std::vector<std::string> &)");

                allocated = true;
            }

            size_t vect_data_size;

            fA.read(reinterpret_cast<char *>(A), sz * sizeof(T));

            fA.read( reinterpret_cast<char *>( &vect_data_size ), sizeof(size_t) ); // 3. reading size of keys
            vect_data.resize(vect_data_size);
            for (size_t i = 0; i < vect_data_size; i++)
            {
                size_t str_size;
                std::string str;
                fA.read( reinterpret_cast<char *>( &str_size ), sizeof(size_t) ); // 5. writing size of ids
                str.resize(str_size);
                fA.read( reinterpret_cast<char *>( &str[0] ), str_size ); // 6. reading specific id
                vect_data[i] = str;
            }

            fA.close();

            ondisk = false;

            resizedElements = sz;
        }
        else
            throw std::string("matrix<T>::fread(const std::string &, , std::vector<std::string> &) => Cannot open the binary file: " + fname + " to read from!");
    }

    //===============================================================================================================

    template <typename T>
    bool matrix<T>::is_ondisk()
    {
        /*
            Checks if data saved on DISK or is in memory.

            Return value: bool.
        */

        if (ondisk)
            return true;
        else
            return false;
    }

    //===============================================================================================================

    template <typename T>
    matrix<T> matrix<T>::fget()
    {
        /*
            Moves all of a saved on DISK matrix back to the memory.
            Differs to fread() in data being copied not to the original matrix itself but
            to another matrix. The original matrix remains on DISK.
            Operates on both, rectangular and symmetric matrix in compact format (lower triangular part).

            Return value: matrix, the type of original matrix.

            Example:
                    int dim = 6000;
                    matrix <mytype> M1( dim );
                    matrix <mytype> res1;

                    ...

                    res1 = M1.fget();
        */

        matrix<T> C;

        if (!ondisk)
            throw std::string("The data was not moved to a binary file. matrix<T>::fget()");

        std::ifstream f(binFilename.c_str());

        if (f.good())
        {
            fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            fA.open(binFilename, fA.binary | fA.in);

            if (!fA.is_open())
            {
                failbit = true;
                throw std::string("Error while opening a binary file. matrix<T>::fget()");
            }

            size_t sz;

            if (!compact)
            {
                sz = numRow * numCol;

                C.resize(numRow, numCol);
            }
            else
            {
                sz = (numCol * numCol + numCol) / 2;

                C.resize(numCol);
            }

            fA.read(reinterpret_cast<char *>(C.A), sz * sizeof(T));

            fA.close();
        }

        return C;
    }

    //===============================================================================================================

    template <typename T>
    matrix<float> matrix<T>::_float()
    {
        /*
            Type conversion method.
            Copy the original matrix of type <T> to a new matrix of type <float>.

            In the case of matrix stored on disk and where only part of it needs to be casted,
            it is recomended to use cast_fget( ... ) method since it not increases memorey
            required for copying.

            Return value: matrix <float>. The original matrix remains unchanged.

            Example:
                    matrix <T> M1( dim );
                    matrix <float> res;

                    ...

                    res1 = M1._float();
        */

        /*if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::_float()");*/

        typedef float cast_type;

        matrix<cast_type> C;

        size_t sz;

        if (!compact)
        {
            C.resize(numRow, numCol);
            sz = numRow * numCol;
        }
        else
        {
            C.resize(numCol);
            sz = (numCol * numCol + numCol) / 2;
        }

        bool write_it = false;

        if (ondisk)
        {
            write_it = true;

            std::ifstream f(binFilename.c_str());

            if (f.good())
            {
                fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
                fA.open(binFilename, fA.binary | fA.in);

                if (!fA.is_open())
                {
                    failbit = true;
                    throw std::string("Error while opening a binary file. matrix<T>::_float()");
                }

                if (!compact)
                {

                    sz = numRow * numCol;

                    int status = allocate(numRow, numCol);
                    if (status != 0)
                        throw std::string("Memory allocation error. matrix<T>::_float()");

                    allocated = true;
                }
                else
                {

                    sz = (numCol * numCol + numCol) / 2;

                    int status = allocate(numCol);
                    if (status != 0)
                        throw std::string("Memory allocation error. matrix<T>::_float()");

                    allocated = true;
                }

                fA.read(reinterpret_cast<char *>(A), sz * sizeof(T));

                fA.close();

                ondisk = false;
            }
        }

        // #pragma omp parallel for
        for (size_t i = 0; i < sz; i++)
        {
            C[i] = static_cast<cast_type>(A[i]);
        }

        if (write_it)
        {
            fwrite();
            C.fwrite();
        }

        return C;
    }

    //===============================================================================================================

    template <typename T>
    matrix<double> matrix<T>::_double()
    {
        /*
            Type conversion method.
            Copy the original matrix of type <T> to a new matrix of type <double>.

            In the case of matrix stored on disk and where only part of it needs to be casted,
            it is recomended to use cast_fget( ... ) method since it not increases memorey
            required for copying.

            Return value: matrix <double>. The original matrix remains unchanged.

            Example:
                    matrix <T> M1( dim );
                    matrix <double> res;

                    ...

                    res1 = M1._double();
        */

        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::_double()");

        typedef double cast_type;

        matrix<cast_type> C;

        size_t sz;

        if (!compact)
        {
            C.resize(numRow, numCol);
            sz = numRow * numCol;
        }
        else
        {
            C.resize(numCol);
            sz = (numCol * numCol + numCol) / 2;
        }

        bool write_it = false;

        if (ondisk)
        {
            write_it = true;

            std::ifstream f(binFilename.c_str());

            if (f.good())
            {
                fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
                fA.open(binFilename, fA.binary | fA.in);

                if (!fA.is_open())
                {
                    failbit = true;
                    throw std::string("Error while opening a binary file. matrix<T>::_double()");
                }

                if (!compact)
                {

                    sz = numRow * numCol;

                    int status = allocate(numRow, numCol);
                    if (status != 0)
                        throw std::string("Memory allocation error. matrix<T>::_double()");

                    allocated = true;
                }
                else
                {

                    sz = (numCol * numCol + numCol) / 2;

                    int status = allocate(numCol);
                    if (status != 0)
                        throw std::string("Memory allocation error. matrix<T>::_double()");

                    allocated = true;
                }

                fA.read(reinterpret_cast<char *>(A), sz * sizeof(T));

                fA.close();

                ondisk = false;
            }
        }

        // #pragma omp parallel for
        for (size_t i = 0; i < sz; i++)
        {
            C[i] = static_cast<cast_type>(A[i]);
        }

        if (write_it)
        {
            fwrite();
            C.fwrite();
        }

        return C;
    }

    //===============================================================================================================

    template <typename T>
    matrix<int> matrix<T>::_int()
    {
        /*
            Type conversion method.
            Copy the original matrix of type <T> to a new matrix of type <int>.

            In the case of matrix stored on disk and where only part of it needs to be casted,
            it is recomended to use cast_fget( ... ) method since it not increases memorey
            required for copying.

            Return value: matrix <int>. The original matrix remains unchanged.

            Example:
                    matrix <T> M1( dim );
                    matrix <int> res;

                    ...

                    res1 = M1._int();
        */

        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::_int()");

        typedef int cast_type;

        matrix<cast_type> C;

        size_t sz;

        if (!compact)
        {
            C.resize(numRow, numCol);
            sz = numRow * numCol;
        }
        else
        {
            C.resize(numCol);
            sz = (numCol * numCol + numCol) / 2;
        }

        bool write_it = false;

        if (ondisk)
        {
            write_it = true;

            std::ifstream f(binFilename.c_str());

            if (f.good())
            {
                fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
                fA.open(binFilename, fA.binary | fA.in);

                if (!fA.is_open())
                {
                    failbit = true;
                    throw std::string("Error while opening a binary file. matrix<T>::_int()");
                }

                if (!compact)
                {

                    sz = numRow * numCol;

                    int status = allocate(numRow, numCol);
                    if (status != 0)
                        throw std::string("Memory allocation error. matrix<T>::_int()");

                    allocated = true;
                }
                else
                {

                    sz = (numCol * numCol + numCol) / 2;

                    int status = allocate(numCol);
                    if (status != 0)
                        throw std::string("Memory allocation error. matrix<T>::_int()");

                    allocated = true;
                }

                fA.read(reinterpret_cast<char *>(A), sz * sizeof(T));

                fA.close();

                ondisk = false;
            }
        }

        // #pragma omp parallel for
        for (size_t i = 0; i < sz; i++)
        {
            C[i] = static_cast<cast_type>(A[i]);
        }

        if (write_it)
        {
            fwrite();
            C.fwrite();
        }

        return C;
    }

    //===============================================================================================================

    template <typename T>
    matrix<size_t> matrix<T>::_size_t()
    {
        /*
            Type conversion method.
            Copy the original matrix of type <T> to a new matrix of type <size_t>.

            In the case of matrix stored on disk and where only part of it needs to be casted,
            it is recomended to use cast_fget( ... ) method since it not increases memorey
            required for copying.

            Return value: matrix <size_t>. The original matrix remains unchanged.

            Example:
                    matrix <T> M1( dim );
                    matrix <size_t> res;

                    ...

                    res1 = M1._size_t();
        */

        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::_size_t()");

        typedef size_t cast_type;

        matrix<cast_type> C;

        size_t sz;

        if (!compact)
        {
            C.resize(numRow, numCol);
            sz = numRow * numCol;
        }
        else
        {
            C.resize(numCol);
            sz = (numCol * numCol + numCol) / 2;
        }

        bool write_it = false;

        if (ondisk)
        {
            write_it = true;

            std::ifstream f(binFilename.c_str());

            if (f.good())
            {
                fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
                fA.open(binFilename, fA.binary | fA.in);

                if (!fA.is_open())
                {
                    failbit = true;
                    throw std::string("Error while opening a binary file. matrix<T>::_size_t()");
                }

                if (!compact)
                {

                    sz = numRow * numCol;

                    int status = allocate(numRow, numCol);
                    if (status != 0)
                        throw std::string("Memory allocation error. matrix<T>::_size_t()");

                    allocated = true;
                }
                else
                {

                    sz = (numCol * numCol + numCol) / 2;

                    int status = allocate(numCol);
                    if (status != 0)
                        throw std::string("Memory allocation error. matrix<T>::_size_t()");

                    allocated = true;
                }

                fA.read(reinterpret_cast<char *>(A), sz * sizeof(T));

                fA.close();

                ondisk = false;
            }
        }

        // #pragma omp parallel for
        for (size_t i = 0; i < sz; i++)
        {
            C[i] = static_cast<cast_type>(A[i]);
        }

        if (write_it)
        {
            fwrite();
            C.fwrite();
        }

        return C;
    }

    //===============================================================================================================

    template <typename T>
    matrix<bool> matrix<T>::_bool()
    {
        /*
            Type conversion method.
            Copy the original matrix of type <T> to a new matrix of type <bool>.

            In the case of matrix stored on disk and where only part of it needs to be casted,
            it is recomended to use cast_fget( ... ) method since it not increases memorey
            required for copying.

            Return value: matrix <bool>. The original matrix remains unchanged.

            Example:
                    matrix <T> M1( dim );
                    matrix <bool> res;

                    ...

                    res1 = M1._bool();
        */

        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::_bool()");

        typedef bool cast_type;

        matrix<cast_type> C;

        size_t sz;

        if (!compact)
        {
            C.resize(numRow, numCol);
            sz = numRow * numCol;
        }
        else
        {
            C.resize(numCol);
            sz = (numCol * numCol + numCol) / 2;
        }

        bool write_it = false;

        if (ondisk)
        {
            write_it = true;

            std::ifstream f(binFilename.c_str());

            if (f.good())
            {
                fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
                fA.open(binFilename, fA.binary | fA.in);

                if (!fA.is_open())
                {
                    failbit = true;
                    throw std::string("Error while opening a binary file. matrix<T>::_bool()");
                }

                if (!compact)
                {

                    sz = numRow * numCol;

                    int status = allocate(numRow, numCol);
                    if (status != 0)
                        throw std::string("Memory allocation error. matrix<T>::_bool()");

                    allocated = true;
                }
                else
                {

                    sz = (numCol * numCol + numCol) / 2;

                    int status = allocate(numCol);
                    if (status != 0)
                        throw std::string("Memory allocation error. matrix<T>::_bool()");

                    allocated = true;
                }

                fA.read(reinterpret_cast<char *>(A), sz * sizeof(T));

                fA.close();

                ondisk = false;
            }
        }

        // #pragma omp parallel for
        for (size_t i = 0; i < sz; i++)
        {
            C[i] = static_cast<cast_type>(A[i]);
        }

        if (write_it)
        {
            fwrite();
            C.fwrite();
        }

        return C;
    }

    //===============================================================================================================

    template <typename T>
    matrix<T> matrix<T>::fget(size_t irow[], size_t icol[])
    {
        /*
            Moves part of a saved on DISK matrix back to the memory.
            Operates on both, rectangular and symmetric matrix in compact format (lower triangular part).

            Return value: always rectangular matrix.

            irow := two element c-style array, where
                    irow[0] is lower row bound of extracted matrix,
                    irow[1] is upper row bound of extracted matrix.
            icol := two element c-style array, where
                    icol[0] is lower column bound of extracted matrix,
                    icol[1] is upper column bound of extracted matrix.

            Example:
                    int dim = 6000;
                    matrix <mytype> M1( dim );
                    matrix <mytype> res1;

                    ...

                    size_t irow[] = {0, dim-1};
                    size_t icol[] = {0, floor(dim/2)};

                    res1 = M1.fget( irow, icol );
        */

        matrix<T> C;

        if (!ondisk)
            throw std::string("The data was not moved to a binary file. matrix<T>::fget(size_t, size_t)");

        std::ifstream f(binFilename.c_str());

        if (f.good())
        {

            fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            fA.open(binFilename, fA.binary | fA.in);

            if (!fA.is_open())
            {
                failbit = true;
                throw std::string("Error while opening a binary file. matrix<T>::fget(size_t, size_t)");
            }

            size_t rows = irow[1] - irow[0] + 1;
            size_t cols = icol[1] - icol[0] + 1;

            matrix<T> a(cols, 1);

            C.resize(rows, cols);

            size_t pos = 0;
            size_t c_rows = 0;

            if (!compact)
            {
                for (size_t i = irow[0]; i <= irow[1]; i++)
                {
                    pos = i * numCol + icol[0];

                    fA.seekg(pos * sizeof(T));
                    fA.read(reinterpret_cast<char *>(a.A), cols * sizeof(T));

                    for (size_t j = 0; j < cols; j++)
                        C(c_rows, j) = a.A[j];

                    c_rows++;
                }
            }
            else
            {
                size_t a_cols = 0;
                size_t k = 0;
                size_t j1 = 0;

                for (size_t i = irow[0]; i <= irow[1]; i++)
                {
                    if (i >= icol[0])
                    {
                        pos = i * (i + 1) / 2 + icol[0];

                        if (a_cols < cols)
                            a_cols = a_cols + 1;

                        fA.seekg(pos * sizeof(T));
                        fA.read(reinterpret_cast<char *>(a.A), a_cols * sizeof(T));

                        for (size_t l = 0; l < a_cols; l++)
                            C(c_rows, l) = a.A[l];
                    }

                    c_rows++;

                    if (i >= icol[1])
                        continue;

                    k = 0;

                    if (i >= icol[0])
                        k = i - icol[0] + 1;

                    j1 = icol[0] + k;

                    for (size_t j = j1; j <= icol[1]; j++)
                    {

                        pos = j * (j + 1) / 2 + i;

                        fA.seekg(pos * sizeof(T));
                        fA.read(reinterpret_cast<char *>(a.A), 1 * sizeof(T));

                        C(i - irow[0], j - icol[0]) = a.A[0];
                    }
                }
            }

            fA.close();
            a.clear();
        }

        return C;
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::cast_fget(size_t irow[], size_t icol[], std::vector<std::vector<float>> &out)
    {
        /*
            Moves part of a saved on DISK matrix back to the memory,
            copy it to a new matrix and
            cast to a new type.

            The mehtod is similar to fget( ... ). The difference is in performed type conversion.

            Operates on both, rectangular and symmetric matrix in compact format (lower triangular part).

            Return value: Data copied to the matrix out which is always rectangular matrix.
                          Memory preallocation for out is not required.

            irow := two element c-style array, where
                    irow[0] is lower row bound of extracted matrix,
                    irow[1] is upper row bound of extracted matrix.
            icol := two element c-style array, where
                    icol[0] is lower column bound of extracted matrix,
                    icol[1] is upper column bound of extracted matrix.
            out  := matrix where the method copies data.

            Example:
                    int dim = 6000;
                    matrix <T> M1( dim );
                    matrix <float> res;

                    ...

                    size_t irow[] = {0, dim-1};
                    size_t icol[] = {0, floor(dim/2)};

                    M1.cast_fget( irow, icol, res );


            UPDATE: also include the case for the matrix on memory.
        */

        typedef float cast_type;

        if (ondisk == false)
        {
            // throw std::string("The data was not moved to a binary file. matrix<T>::cast_fget(size_t, size_t, std::vector<std::vector<float>> &)");

            // size_t rows = irow[1] - irow[0] + 1;
            // size_t cols = icol[1] - icol[0] + 1;

            // out.resize(rows, cols);

            if (!compact)
            {
                for (size_t i = irow[0]; i <= irow[1]; i++)
                {
                    for (size_t j = icol[0]; j <= icol[1]; j++)
                    {
                        // out(i - irow[0], j - icol[0]) = static_cast<cast_type>(A[i * numCol + j]);
                        out[i - irow[0]][j - icol[0]] = static_cast<cast_type>(A[i * numCol + j]);
                    }
                }
            }
            else
            {
                throw std::string("This method is not implemented for compact format matrix. matrix<T>::cast_fget(size_t, size_t, std::vector<std::vector<float>> &)");
            }
        }
        else
        {
            std::ifstream f(binFilename.c_str());

            if (f.good())
            {

                fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
                fA.open(binFilename, fA.binary | fA.in);

                if (!fA.is_open())
                {
                    failbit = true;
                    throw std::string("Error while opening a binary file. matrix<T>::cast_fget(size_t, size_t, std::vector<std::vector<float>> &)");
                }

                // size_t rows = irow[1] - irow[0] + 1;
                size_t cols = icol[1] - icol[0] + 1;

                matrix<T> a(cols, 1);

                // out.resize(rows, cols);

                size_t pos = 0;
                size_t c_rows = 0;

                if (!compact)
                {
                    for (size_t i = irow[0]; i <= irow[1]; i++)
                    {
                        pos = i * numCol + icol[0];

                        fA.seekg(pos * sizeof(T));
                        fA.read(reinterpret_cast<char *>(a.A), cols * sizeof(T));

                        for (size_t j = 0; j < cols; j++)
                            out[c_rows][j] = static_cast<cast_type>(a.A[j]);

                        c_rows++;
                    }
                }
                else
                {
                    size_t a_cols = 0;
                    size_t k = 0;
                    size_t j1 = 0;

                    for (size_t i = irow[0]; i <= irow[1]; i++)
                    {
                        if (i >= icol[0])
                        {
                            pos = i * (i + 1) / 2 + icol[0];

                            if (a_cols < cols)
                                a_cols = a_cols + 1;

                            fA.seekg(pos * sizeof(T));
                            fA.read(reinterpret_cast<char *>(a.A), a_cols * sizeof(T));

                            for (size_t l = 0; l < a_cols; l++)
                                out[c_rows][l] = static_cast<cast_type>(a.A[l]);
                        }

                        c_rows++;

                        if (i >= icol[1])
                            continue;

                        k = 0;

                        if (i >= icol[0])
                            k = i - icol[0] + 1;

                        j1 = icol[0] + k;

                        for (size_t j = j1; j <= icol[1]; j++)
                        {

                            pos = j * (j + 1) / 2 + i;

                            fA.seekg(pos * sizeof(T));
                            fA.read(reinterpret_cast<char *>(a.A), 1 * sizeof(T));

                            out[i - irow[0]][j - icol[0]] = static_cast<cast_type>(a.A[0]);
                        }
                    }
                }

                fA.close();
                a.clear();
            }
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::cast_fget(size_t irow[], size_t icol[], float **out)
    {
        /*
            Moves part of a saved on DISK matrix back to the memory,
            copy it to a new matrix and
            cast to a new type.

            The mehtod is similar to fget( ... ). The difference is in performed type conversion.

            Operates on both, rectangular and symmetric matrix in compact format (lower triangular part).

            Return value: Data copied to the matrix out which is always rectangular matrix.
                          Memory preallocation for out is not required.

            irow := two element c-style array, where
                    irow[0] is lower row bound of extracted matrix,
                    irow[1] is upper row bound of extracted matrix.
            icol := two element c-style array, where
                    icol[0] is lower column bound of extracted matrix,
                    icol[1] is upper column bound of extracted matrix.
            out  := matrix where the method copies data.

            Example:
                    int dim = 6000;
                    matrix <T> M1( dim );
                    matrix <float> res;

                    ...

                    size_t irow[] = {0, dim-1};
                    size_t icol[] = {0, floor(dim/2)};

                    M1.cast_fget( irow, icol, res );


            UPDATE: also include the case for the matrix on memory.
        */

        typedef float cast_type;

        if (ondisk == false)
        {
            // throw std::string("The data was not moved to a binary file. matrix<T>::cast_fget(size_t, size_t, float **)");

            // size_t rows = irow[1] - irow[0] + 1;
            // size_t cols = icol[1] - icol[0] + 1;

            // out.resize(rows, cols);

            if (!compact)
            {
                for (size_t i = irow[0]; i <= irow[1]; i++)
                {
                    for (size_t j = icol[0]; j <= icol[1]; j++)
                    {
                        // out(i - irow[0], j - icol[0]) = static_cast<cast_type>(A[i * numCol + j]);
                        out[i - irow[0]][j - icol[0]] = static_cast<cast_type>(A[i * numCol + j]);
                    }
                }
            }
            else
            {
                throw std::string("This method is not implemented for compact format matrix. matrix<T>::cast_fget(size_t, size_t, float **)");
            }
        }
        else
        {
            std::ifstream f(binFilename.c_str());

            if (f.good())
            {

                fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
                fA.open(binFilename, fA.binary | fA.in);

                if (!fA.is_open())
                {
                    failbit = true;
                    throw std::string("Error while opening a binary file. matrix<T>::cast_fget(size_t, size_t, float **)");
                }

                // size_t rows = irow[1] - irow[0] + 1;
                size_t cols = icol[1] - icol[0] + 1;

                matrix<T> a(cols, 1);

                // out.resize(rows, cols);

                size_t pos = 0;
                size_t c_rows = 0;

                if (!compact)
                {
                    for (size_t i = irow[0]; i <= irow[1]; i++)
                    {
                        pos = i * numCol + icol[0];

                        fA.seekg(pos * sizeof(T));
                        fA.read(reinterpret_cast<char *>(a.A), cols * sizeof(T));

                        for (size_t j = 0; j < cols; j++)
                            out[c_rows][j] = static_cast<cast_type>(a.A[j]);

                        c_rows++;
                    }
                }
                else
                {
                    size_t a_cols = 0;
                    size_t k = 0;
                    size_t j1 = 0;

                    for (size_t i = irow[0]; i <= irow[1]; i++)
                    {
                        if (i >= icol[0])
                        {
                            pos = i * (i + 1) / 2 + icol[0];

                            if (a_cols < cols)
                                a_cols = a_cols + 1;

                            fA.seekg(pos * sizeof(T));
                            fA.read(reinterpret_cast<char *>(a.A), a_cols * sizeof(T));

                            for (size_t l = 0; l < a_cols; l++)
                                out[c_rows][l] = static_cast<cast_type>(a.A[l]);
                        }

                        c_rows++;

                        if (i >= icol[1])
                            continue;

                        k = 0;

                        if (i >= icol[0])
                            k = i - icol[0] + 1;

                        j1 = icol[0] + k;

                        for (size_t j = j1; j <= icol[1]; j++)
                        {

                            pos = j * (j + 1) / 2 + i;

                            fA.seekg(pos * sizeof(T));
                            fA.read(reinterpret_cast<char *>(a.A), 1 * sizeof(T));

                            out[i - irow[0]][j - icol[0]] = static_cast<cast_type>(a.A[0]);
                        }
                    }
                }

                fA.close();
                a.clear();
            }
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::cast_fget(size_t irow[], size_t icol[], matrix<float> &out)
    {
        /*
            Moves part of a saved on DISK matrix back to the memory,
            copy it to a new matrix and
            cast to a new type.

            The mehtod is similar to fget( ... ). The difference is in performed type conversion.

            Operates on both, rectangular and symmetric matrix in compact format (lower triangular part).

            Return value: Data copied to the matrix out which is always rectangular matrix.
                          Memory preallocation for out is not required.

            irow := two element c-style array, where
                    irow[0] is lower row bound of extracted matrix,
                    irow[1] is upper row bound of extracted matrix.
            icol := two element c-style array, where
                    icol[0] is lower column bound of extracted matrix,
                    icol[1] is upper column bound of extracted matrix.
            out  := matrix where the method copies data.

            Example:
                    int dim = 6000;
                    matrix <T> M1( dim );
                    matrix <float> res;

                    ...

                    size_t irow[] = {0, dim-1};
                    size_t icol[] = {0, floor(dim/2)};

                    M1.cast_fget( irow, icol, res );


            UPDATE: also include the case for the matrix on memory.
        */

        typedef float cast_type;

        // matrix<cast_type> C;

        if (ondisk == false)
        {
            // throw std::string("The data was not moved to a binary file. matrix<T>::cast_fget(size_t, size_t, matrix<float> &)");

            size_t rows = irow[1] - irow[0] + 1;
            size_t cols = icol[1] - icol[0] + 1;

            out.resize(rows, cols);

            if (!compact)
            {
                for (size_t i = irow[0]; i <= irow[1]; i++)
                {
                    for (size_t j = icol[0]; j <= icol[1]; j++)
                    {
                        out(i - irow[0], j - icol[0]) = static_cast<cast_type>(A[i * numCol + j]);
                    }
                }
            }
            else
            {
                throw std::string("This method is not implemented for compact format matrix. matrix<T>::cast_fget(size_t, size_t, matrix<float> &)");
            }
        }
        else
        {
            std::ifstream f(binFilename.c_str());

            if (f.good())
            {

                fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
                fA.open(binFilename, fA.binary | fA.in);

                if (!fA.is_open())
                {
                    failbit = true;
                    throw std::string("Error while opening a binary file. matrix<T>::cast_fget(size_t, size_t, matrix<float> &)");
                }

                size_t rows = irow[1] - irow[0] + 1;
                size_t cols = icol[1] - icol[0] + 1;

                matrix<T> a(cols, 1);

                out.resize(rows, cols);

                size_t pos = 0;
                size_t c_rows = 0;

                if (!compact)
                {
                    for (size_t i = irow[0]; i <= irow[1]; i++)
                    {
                        pos = i * numCol + icol[0];

                        fA.seekg(pos * sizeof(T));
                        fA.read(reinterpret_cast<char *>(a.A), cols * sizeof(T));

                        for (size_t j = 0; j < cols; j++)
                            out(c_rows, j) = static_cast<cast_type>(a.A[j]);

                        c_rows++;
                    }
                }
                else
                {
                    size_t a_cols = 0;
                    size_t k = 0;
                    size_t j1 = 0;

                    for (size_t i = irow[0]; i <= irow[1]; i++)
                    {
                        if (i >= icol[0])
                        {
                            pos = i * (i + 1) / 2 + icol[0];

                            if (a_cols < cols)
                                a_cols = a_cols + 1;

                            fA.seekg(pos * sizeof(T));
                            fA.read(reinterpret_cast<char *>(a.A), a_cols * sizeof(T));

                            for (size_t l = 0; l < a_cols; l++)
                                out(c_rows, l) = static_cast<cast_type>(a.A[l]);
                        }

                        c_rows++;

                        if (i >= icol[1])
                            continue;

                        k = 0;

                        if (i >= icol[0])
                            k = i - icol[0] + 1;

                        j1 = icol[0] + k;

                        for (size_t j = j1; j <= icol[1]; j++)
                        {

                            pos = j * (j + 1) / 2 + i;

                            fA.seekg(pos * sizeof(T));
                            fA.read(reinterpret_cast<char *>(a.A), 1 * sizeof(T));

                            out(i - irow[0], j - icol[0]) = static_cast<cast_type>(a.A[0]);
                        }
                    }
                }

                fA.close();
                a.clear();
            }
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::cast_fget(size_t irow[], size_t icol[], matrix<double> &out)
    {
        /*
            Moves part of a saved on DISK matrix back to the memory,
            copy it to a new matrix and
            cast to a new type.

            The mehtod is similar to fget( ... ). The difference is in performed type conversion.

            Operates on both, rectangular and symmetric matrix in compact format (lower triangular part).

            Return value: Data copied to the matrix out which is always rectangular matrix.
                          Memory preallocation for out is not required.

            irow := two element c-style array, where
                    irow[0] is lower row bound of extracted matrix,
                    irow[1] is upper row bound of extracted matrix.
            icol := two element c-style array, where
                    icol[0] is lower column bound of extracted matrix,
                    icol[1] is upper column bound of extracted matrix.
            out  := matrix where the method copies data.

            Example:
                    int dim = 6000;
                    matrix <T> M1( dim );
                    matrix <float> res;

                    ...

                    size_t irow[] = {0, dim-1};
                    size_t icol[] = {0, floor(dim/2)};

                    M1.cast_fget( irow, icol, res );
        */

        typedef double cast_type;

        // matrix<cast_type> C;

        if (ondisk == false)
        {
            // throw std::string("The data was not moved to a binary file. matrix<T>::cast_fget(size_t, size_t, matrix<double> &)");

            size_t rows = irow[1] - irow[0] + 1;
            size_t cols = icol[1] - icol[0] + 1;

            out.resize(rows, cols);

            if (!compact)
            {
                for (size_t i = irow[0]; i <= irow[1]; i++)
                {
                    for (size_t j = icol[0]; j <= icol[1]; j++)
                    {
                        out(i - irow[0], j - icol[0]) = static_cast<cast_type>(A[i * numCol + j]);
                    }
                }
            }
            else
            {
                throw std::string("This method is not implemented for compact format matrix. matrix<T>::cast_fget(size_t, size_t, matrix<double> &)");
            }
        }
        else
        {
            std::ifstream f(binFilename.c_str());

            if (f.good())
            {

                fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
                fA.open(binFilename, fA.binary | fA.in);

                if (!fA.is_open())
                {
                    failbit = true;
                    throw std::string("Error while opening a binary file. matrix<T>::cast_fget(size_t, size_t, matrix<double> &)");
                }

                size_t rows = irow[1] - irow[0] + 1;
                size_t cols = icol[1] - icol[0] + 1;

                matrix<T> a(cols, 1);

                out.resize(rows, cols);

                size_t pos = 0;
                size_t c_rows = 0;

                if (!compact)
                {
                    for (size_t i = irow[0]; i <= irow[1]; i++)
                    {
                        pos = i * numCol + icol[0];

                        fA.seekg(pos * sizeof(T));
                        fA.read(reinterpret_cast<char *>(a.A), cols * sizeof(T));

                        for (size_t j = 0; j < cols; j++)
                            out(c_rows, j) = static_cast<cast_type>(a.A[j]);

                        c_rows++;
                    }
                }
                else
                {
                    size_t a_cols = 0;
                    size_t k = 0;
                    size_t j1 = 0;

                    for (size_t i = irow[0]; i <= irow[1]; i++)
                    {
                        if (i >= icol[0])
                        {
                            pos = i * (i + 1) / 2 + icol[0];

                            if (a_cols < cols)
                                a_cols = a_cols + 1;

                            fA.seekg(pos * sizeof(T));
                            fA.read(reinterpret_cast<char *>(a.A), a_cols * sizeof(T));

                            for (size_t l = 0; l < a_cols; l++)
                                out(c_rows, l) = static_cast<cast_type>(a.A[l]);
                        }

                        c_rows++;

                        if (i >= icol[1])
                            continue;

                        k = 0;

                        if (i >= icol[0])
                            k = i - icol[0] + 1;

                        j1 = icol[0] + k;

                        for (size_t j = j1; j <= icol[1]; j++)
                        {

                            pos = j * (j + 1) / 2 + i;

                            fA.seekg(pos * sizeof(T));
                            fA.read(reinterpret_cast<char *>(a.A), 1 * sizeof(T));

                            out(i - irow[0], j - icol[0]) = static_cast<cast_type>(a.A[0]);
                        }
                    }
                }

                fA.close();
                a.clear();
            }
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::cast_fget(size_t irow[], size_t icol[], matrix<int> &out)
    {
        /*
            Moves part of a saved on DISK matrix back to the memory,
            copy it to a new matrix and
            cast to a new type.

            The mehtod is similar to fget( ... ). The difference is in performed type conversion.

            Operates on both, rectangular and symmetric matrix in compact format (lower triangular part).

            Return value: Data copied to the matrix out which is always rectangular matrix.
                          Memory preallocation for out is not required.

            irow := two element c-style array, where
                    irow[0] is lower row bound of extracted matrix,
                    irow[1] is upper row bound of extracted matrix.
            icol := two element c-style array, where
                    icol[0] is lower column bound of extracted matrix,
                    icol[1] is upper column bound of extracted matrix.
            out  := matrix where the method copies data.

            Example:
                    int dim = 6000;
                    matrix <T> M1( dim );
                    matrix <float> res;

                    ...

                    size_t irow[] = {0, dim-1};
                    size_t icol[] = {0, floor(dim/2)};

                    M1.cast_fget( irow, icol, res );
        */

        typedef int cast_type;

        // matrix<cast_type> C;

        if (ondisk == false)
        {
            // throw std::string("The data was not moved to a binary file. matrix<T>::cast_fget(size_t, size_t, matrix<int> &)");

            size_t rows = irow[1] - irow[0] + 1;
            size_t cols = icol[1] - icol[0] + 1;

            out.resize(rows, cols);

            if (!compact)
            {
                for (size_t i = irow[0]; i <= irow[1]; i++)
                {
                    for (size_t j = icol[0]; j <= icol[1]; j++)
                    {
                        out(i - irow[0], j - icol[0]) = static_cast<cast_type>(A[i * numCol + j]);
                    }
                }
            }
            else
            {
                throw std::string("This method is not implemented for compact format matrix. matrix<T>::cast_fget(size_t, size_t, matrix<int> &)");
            }
        }
        else
        {
            std::ifstream f(binFilename.c_str());

            if (f.good())
            {

                fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
                fA.open(binFilename, fA.binary | fA.in);

                if (!fA.is_open())
                {
                    failbit = true;
                    throw std::string("Error while opening a binary file. matrix<T>::cast_fget(size_t, size_t, matrix<int> &)");
                }

                size_t rows = irow[1] - irow[0] + 1;
                size_t cols = icol[1] - icol[0] + 1;

                matrix<T> a(cols, 1);

                out.resize(rows, cols);

                size_t pos = 0;
                size_t c_rows = 0;

                if (!compact)
                {
                    for (size_t i = irow[0]; i <= irow[1]; i++)
                    {
                        pos = i * numCol + icol[0];

                        fA.seekg(pos * sizeof(T));
                        fA.read(reinterpret_cast<char *>(a.A), cols * sizeof(T));

                        for (size_t j = 0; j < cols; j++)
                            out(c_rows, j) = static_cast<cast_type>(a.A[j]);

                        c_rows++;
                    }
                }
                else
                {
                    size_t a_cols = 0;
                    size_t k = 0;
                    size_t j1 = 0;

                    for (size_t i = irow[0]; i <= irow[1]; i++)
                    {
                        if (i >= icol[0])
                        {
                            pos = i * (i + 1) / 2 + icol[0];

                            if (a_cols < cols)
                                a_cols = a_cols + 1;

                            fA.seekg(pos * sizeof(T));
                            fA.read(reinterpret_cast<char *>(a.A), a_cols * sizeof(T));

                            for (size_t l = 0; l < a_cols; l++)
                                out(c_rows, l) = static_cast<cast_type>(a.A[l]);
                        }

                        c_rows++;

                        if (i >= icol[1])
                            continue;

                        k = 0;

                        if (i >= icol[0])
                            k = i - icol[0] + 1;

                        j1 = icol[0] + k;

                        for (size_t j = j1; j <= icol[1]; j++)
                        {

                            pos = j * (j + 1) / 2 + i;

                            fA.seekg(pos * sizeof(T));
                            fA.read(reinterpret_cast<char *>(a.A), 1 * sizeof(T));

                            out(i - irow[0], j - icol[0]) = static_cast<cast_type>(a.A[0]);
                        }
                    }
                }

                fA.close();
                a.clear();
            }
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::cast_fget(size_t irow[], size_t icol[], matrix<size_t> &out)
    {
        /*
            Moves part of a saved on DISK matrix back to the memory,
            copy it to a new matrix and
            cast to a new type.

            The mehtod is similar to fget( ... ). The difference is in performed type conversion.

            Operates on both, rectangular and symmetric matrix in compact format (lower triangular part).

            Return value: Data copied to the matrix out which is always rectangular matrix.
                          Memory preallocation for out is not required.

            irow := two element c-style array, where
                    irow[0] is lower row bound of extracted matrix,
                    irow[1] is upper row bound of extracted matrix.
            icol := two element c-style array, where
                    icol[0] is lower column bound of extracted matrix,
                    icol[1] is upper column bound of extracted matrix.
            out  := matrix where the method copies data.

            Example:
                    int dim = 6000;
                    matrix <T> M1( dim );
                    matrix <float> res;

                    ...

                    size_t irow[] = {0, dim-1};
                    size_t icol[] = {0, floor(dim/2)};

                    M1.cast_fget( irow, icol, res );
        */

        typedef size_t cast_type;

        // matrix<cast_type> C;

        if (ondisk == false)
        {
            // throw std::string("The data was not moved to a binary file. matrix<T>::cast_fget(size_t, size_t, matrix<size_t> &)");

            size_t rows = irow[1] - irow[0] + 1;
            size_t cols = icol[1] - icol[0] + 1;

            out.resize(rows, cols);

            if (!compact)
            {
                for (size_t i = irow[0]; i <= irow[1]; i++)
                {
                    for (size_t j = icol[0]; j <= icol[1]; j++)
                    {
                        out(i - irow[0], j - icol[0]) = static_cast<cast_type>(A[i * numCol + j]);
                    }
                }
            }
            else
            {
                throw std::string("This method is not implemented for compact format matrix. matrix<T>::cast_fget(size_t, size_t, matrix<size_t> &)");
            }
        }
        else
        {

            std::ifstream f(binFilename.c_str());

            if (f.good())
            {

                fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
                fA.open(binFilename, fA.binary | fA.in);

                if (!fA.is_open())
                {
                    failbit = true;
                    throw std::string("Error while opening a binary file. matrix<T>::cast_fget(size_t, size_t, matrix<size_t> &)");
                }

                size_t rows = irow[1] - irow[0] + 1;
                size_t cols = icol[1] - icol[0] + 1;

                matrix<T> a(cols, 1);

                out.resize(rows, cols);

                size_t pos = 0;
                size_t c_rows = 0;

                if (!compact)
                {
                    for (size_t i = irow[0]; i <= irow[1]; i++)
                    {
                        pos = i * numCol + icol[0];

                        fA.seekg(pos * sizeof(T));
                        fA.read(reinterpret_cast<char *>(a.A), cols * sizeof(T));

                        for (size_t j = 0; j < cols; j++)
                            out(c_rows, j) = static_cast<cast_type>(a.A[j]);

                        c_rows++;
                    }
                }
                else
                {
                    size_t a_cols = 0;
                    size_t k = 0;
                    size_t j1 = 0;

                    for (size_t i = irow[0]; i <= irow[1]; i++)
                    {
                        if (i >= icol[0])
                        {
                            pos = i * (i + 1) / 2 + icol[0];

                            if (a_cols < cols)
                                a_cols = a_cols + 1;

                            fA.seekg(pos * sizeof(T));
                            fA.read(reinterpret_cast<char *>(a.A), a_cols * sizeof(T));

                            for (size_t l = 0; l < a_cols; l++)
                                out(c_rows, l) = static_cast<cast_type>(a.A[l]);
                        }

                        c_rows++;

                        if (i >= icol[1])
                            continue;

                        k = 0;

                        if (i >= icol[0])
                            k = i - icol[0] + 1;

                        j1 = icol[0] + k;

                        for (size_t j = j1; j <= icol[1]; j++)
                        {

                            pos = j * (j + 1) / 2 + i;

                            fA.seekg(pos * sizeof(T));
                            fA.read(reinterpret_cast<char *>(a.A), 1 * sizeof(T));

                            out(i - irow[0], j - icol[0]) = static_cast<cast_type>(a.A[0]);
                        }
                    }
                }

                fA.close();
                a.clear();
            }
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::cast_fget(size_t irow[], size_t icol[], matrix<bool> &out)
    {
        /*
            Moves part of a saved on DISK matrix back to the memory,
            copy it to a new matrix and
            cast to a new type.

            The mehtod is similar to fget( ... ). The difference is in performed type conversion.

            Operates on both, rectangular and symmetric matrix in compact format (lower triangular part).

            Return value: Data copied to the matrix out which is always rectangular matrix.
                          Memory preallocation for out is not required.

            irow := two element c-style array, where
                    irow[0] is lower row bound of extracted matrix,
                    irow[1] is upper row bound of extracted matrix.
            icol := two element c-style array, where
                    icol[0] is lower column bound of extracted matrix,
                    icol[1] is upper column bound of extracted matrix.
            out  := matrix where the method copies data.

            Example:
                    int dim = 6000;
                    matrix <T> M1( dim );
                    matrix <float> res;

                    ...

                    size_t irow[] = {0, dim-1};
                    size_t icol[] = {0, floor(dim/2)};

                    M1.cast_fget( irow, icol, res );
        */

        typedef bool cast_type;

        // matrix<cast_type> C;

        if (ondisk == false)
        {
            // throw std::string("The data was not moved to a binary file. matrix<T>::cast_fget(size_t, size_t, matrix<bool> &)");

            size_t rows = irow[1] - irow[0] + 1;
            size_t cols = icol[1] - icol[0] + 1;

            out.resize(rows, cols);

            if (!compact)
            {
                for (size_t i = irow[0]; i <= irow[1]; i++)
                {
                    for (size_t j = icol[0]; j <= icol[1]; j++)
                    {
                        out(i - irow[0], j - icol[0]) = static_cast<cast_type>(A[i * numCol + j]);
                    }
                }
            }
            else
            {
                throw std::string("This method is not implemented for compact format matrix. matrix<T>::cast_fget(size_t, size_t, matrix<bool> &)");
            }
        }
        else
        {
            std::ifstream f(binFilename.c_str());

            if (f.good())
            {

                fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
                fA.open(binFilename, fA.binary | fA.in);

                if (!fA.is_open())
                {
                    failbit = true;
                    throw std::string("Error while opening a binary file. matrix<T>::cast_fget(size_t, size_t, matrix<bool> &)");
                }

                size_t rows = irow[1] - irow[0] + 1;
                size_t cols = icol[1] - icol[0] + 1;

                matrix<T> a(cols, 1);

                out.resize(rows, cols);

                size_t pos = 0;
                size_t c_rows = 0;

                if (!compact)
                {
                    for (size_t i = irow[0]; i <= irow[1]; i++)
                    {
                        pos = i * numCol + icol[0];

                        fA.seekg(pos * sizeof(T));
                        fA.read(reinterpret_cast<char *>(a.A), cols * sizeof(T));

                        for (size_t j = 0; j < cols; j++)
                            out(c_rows, j) = static_cast<cast_type>(a.A[j]);

                        c_rows++;
                    }
                }
                else
                {
                    size_t a_cols = 0;
                    size_t k = 0;
                    size_t j1 = 0;

                    for (size_t i = irow[0]; i <= irow[1]; i++)
                    {
                        if (i >= icol[0])
                        {
                            pos = i * (i + 1) / 2 + icol[0];

                            if (a_cols < cols)
                                a_cols = a_cols + 1;

                            fA.seekg(pos * sizeof(T));
                            fA.read(reinterpret_cast<char *>(a.A), a_cols * sizeof(T));

                            for (size_t l = 0; l < a_cols; l++)
                                out(c_rows, l) = static_cast<cast_type>(a.A[l]);
                        }

                        c_rows++;

                        if (i >= icol[1])
                            continue;

                        k = 0;

                        if (i >= icol[0])
                            k = i - icol[0] + 1;

                        j1 = icol[0] + k;

                        for (size_t j = j1; j <= icol[1]; j++)
                        {

                            pos = j * (j + 1) / 2 + i;

                            fA.seekg(pos * sizeof(T));
                            fA.read(reinterpret_cast<char *>(a.A), 1 * sizeof(T));

                            out(i - irow[0], j - icol[0]) = static_cast<cast_type>(a.A[0]);
                        }
                    }
                }

                fA.close();
                a.clear();
            }
        }
    }

    //===============================================================================================================
    
    template <typename T>
    void matrix<T>::linsolve(matrix<T> &rhs)
    {
        // NOTE: in this setup A will be overwriten by LU factorization!
        //       rhs will be overwriten by solution.

        if ( numRow != numCol )
            throw std::string("A square matrix is required for matrix<T>::linsolve(matrix<T> &rhs) method!");

        if ( is_symmetric() )
            linsolve_chol(A, numRow, numCol, rhs.A, rhs.numCol);
        else
            linsolve_lu(A, numRow, numCol, rhs.A, rhs.numCol);
    }

    //===============================================================================================================
    
    template <typename T>
    void matrix<T>::linsolve_fact(matrix<T> &rhs)
    {
        // NOTE: in this setup A is already factorized!
        //       rhs will be overwriten by solution.

        if ( numRow != numCol )
            throw std::string("A square matrix is required for matrix<T>::linsolve(matrix<T> &rhs) method!");

        if ( is_symmetric() )
            linsolve_chol_fact(A, numRow/*, numCol*/, rhs.A, rhs.numCol);
        else
            linsolve_lu_fact(A, numRow/*, numCol*/, rhs.A, rhs.numCol);
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::linsolve_lu(double *_A, MKL_INT rowA, MKL_INT colA, double *_B, MKL_INT colB)
    {
        // 1. First, make LU factorizatioon of A: A_fact = P*L*U; ipiv = P.

        lapack_int info = 0;
        lapack_int row = rowA;
        lapack_int col = colA;
        lapack_int col_b = colB;
        int matrix_order = LAPACK_ROW_MAJOR;

        lapack_int *ipiv;
        lapack_int sz_ipiv = std::max(row, col);

#ifdef intelmkl
        ipiv = (lapack_int *)mkl_malloc(sz_ipiv * sizeof(lapack_int), sizeof(T) * 8);
#else
        ipiv = (lapack_int *)malloc(sz_ipiv * sizeof(lapack_int));
#endif
        if (ipiv == NULL)
        {
#ifdef intelmkl
            mkl_free(ipiv);
#else
            free(ipiv);
#endif
            failbit = true;
            throw std::string("Memory allocation error. matrix<T>::linsolve_lu(...)");
        }

        for (lapack_int i = 0; i < sz_ipiv; i++)
            ipiv[i] = 1;
        
        info = LAPACKE_dgetrf(matrix_order, row, col, _A, col, ipiv);
    
        if (info != 0)
        {
#ifdef intelmkl
            mkl_free(ipiv);
#else
            free(ipiv);
#endif
            throw std::string("Error during computation of the LU factorization of a general m-by-n matrix. matrix<T>::linsolve_lu(...)");
        }

        // saving ipiv to the class member ipiv_A
        if ( is_allocated_ipv_A )
        {
            throw std::string("The matrix is already factorized! matrix<T>::linsolve_lu(...)");
        }
        else
        {
            ipiv_A = ipiv;
            is_allocated_ipv_A = true;
        }
        
        // 2. Solves a SQUARE system of linear equations with an LU-factored square coefficient matrix

        info = LAPACKE_dgetrs( matrix_order, 'N', row, col_b, _A, row, ipiv, _B, col_b );

        if (info != 0)
        {
#ifdef intelmkl
            mkl_free(ipiv);
#else
            free(ipiv);
#endif
            is_allocated_ipv_A = false;
            throw std::string("Error during solving a linear system with LU-factored general matrix. matrix<T>::linsolve_lu(...)");
        }
// #ifdef intelmkl
//         mkl_free(ipiv);
// #else
//         free(ipiv);
// #endif

    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::linsolve_lu(float *_A, MKL_INT rowA, MKL_INT colA, float *_B, MKL_INT colB)
    {
        // 1. First, make LU factorizatioon of A: A_fact = P*L*U; ipiv = P.

        lapack_int info = 0;
        lapack_int row = rowA;
        lapack_int col = colA;
        lapack_int col_b = colB;
        int matrix_order = LAPACK_ROW_MAJOR;

        lapack_int *ipiv;
        lapack_int sz_ipiv = std::max(row, col);

#ifdef intelmkl
        ipiv = (lapack_int *)mkl_malloc(sz_ipiv * sizeof(lapack_int), sizeof(T) * 8);
#else
        ipiv = (lapack_int *)malloc(sz_ipiv * sizeof(lapack_int));
#endif
        if (ipiv == NULL)
        {
#ifdef intelmkl
            mkl_free(ipiv);
#else
            free(ipiv);
#endif
            failbit = true;
            throw std::string("Memory allocation error. matrix<T>::linsolve_lu(...)");
        }
        for (lapack_int i = 0; i < sz_ipiv; i++)
            ipiv[i] = 1;

        info = LAPACKE_sgetrf(matrix_order, row, col, _A, col, ipiv);

        if (info != 0)
        {
#ifdef intelmkl
            mkl_free(ipiv);
#else
            free(ipiv);
#endif
            throw std::string("Error during computation of the LU factorization of a general m-by-n matrix. matrix<T>::linsolve_lu(...)");
        }

        // saving ipiv to the class member ipiv_A
        if ( is_allocated_ipv_A )
        {
            throw std::string("The matrix is already factorized! matrix<T>::linsolve_lu(...)");
        }
        else
        {
            ipiv_A = ipiv;
            is_allocated_ipv_A = true;
        }

        // 2. Solves a SQUARE system of linear equations with an LU-factored square coefficient matrix

        info = LAPACKE_sgetrs( matrix_order, 'N', row, col_b, _A, row, ipiv, _B, col_b );

        if (info != 0)
        {
#ifdef intelmkl
            mkl_free(ipiv);
#else
            free(ipiv);
#endif
            is_allocated_ipv_A = false;
            throw std::string("Error during solving a linear system with LU-factored general matrix. matrix<T>::linsolve_lu(...)");
        }
// #ifdef intelmkl
//         mkl_free(ipiv);
// #else
//         free(ipiv);
// #endif
    }
    
    //===============================================================================================================

    template <typename T>
    void matrix<T>::linsolve_lu_fact(double *_A, MKL_INT rowA, /*MKL_INT colA,*/ double *_B, MKL_INT colB)
    {
        // Solving linear system using already factorised A.

        lapack_int info = 0;
        lapack_int row = rowA;
        //lapack_int col = colA;
        lapack_int col_b = colB;
        int matrix_order = LAPACK_ROW_MAJOR;

        if ( !is_allocated_ipv_A )
        {
            throw std::string("The matrix is not factorized! matrix<T>::linsolve_lu(...)");
        }
        
        // Solves a SQUARE system of linear equations with an LU-factored square coefficient matrix

        info = LAPACKE_dgetrs( matrix_order, 'N', row, col_b, _A, row, ipiv_A, _B, col_b );

        if (info != 0)
        {
#ifdef intelmkl
            mkl_free(ipiv_A);
#else
            free(ipiv_A);
#endif
            is_allocated_ipv_A = false;
            throw std::string("Error during solving a linear system with LU-factored general matrix. matrix<T>::linsolve_lu_fact(...)");
        }
// #ifdef intelmkl
//         mkl_free(ipiv);
// #else
//         free(ipiv);
// #endif

    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::linsolve_lu_fact(float *_A, MKL_INT rowA, /*MKL_INT colA,*/ float *_B, MKL_INT colB)
    {
        // Solving linear system using already factorised A.

        lapack_int info = 0;
        lapack_int row = rowA;
        //lapack_int col = colA;
        lapack_int col_b = colB;
        int matrix_order = LAPACK_ROW_MAJOR;
        
        if ( !is_allocated_ipv_A )
        {
            throw std::string("The matrix is not factorized! matrix<T>::linsolve_lu(...)");
        }

        // Solves a SQUARE system of linear equations with an LU-factored square coefficient matrix

        info = LAPACKE_sgetrs( matrix_order, 'N', row, col_b, _A, row, ipiv_A, _B, col_b );

        if (info != 0)
        {
#ifdef intelmkl
            mkl_free(ipiv_A);
#else
            free(ipiv_A);
#endif
            is_allocated_ipv_A = false;
            throw std::string("Error during solving a linear system with LU-factored general matrix. matrix<T>::linsolve_lu_fact(...)");
        }
// #ifdef intelmkl
//         mkl_free(ipiv);
// #else
//         free(ipiv);
// #endif
std::cout<<"using ipiv_A, completed and exited. "<<'\n';
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::linsolve_chol(double *_A, MKL_INT rowA, MKL_INT colA, double *_B, MKL_INT colB)
    {
        // 1. First, make Cholesky factorizatioon of symmetric (compact) A: A_fact = L*L^T.

        lapack_int info = 0;
        lapack_int row = rowA;
        lapack_int col = colA;
        lapack_int col_b = colB;
        int matrix_order = LAPACK_ROW_MAJOR;
        
        info = LAPACKE_dpptrf(matrix_order, 'L', col, _A);
        if (info != 0)
        {
            failbit = true;
            failinfo = info;
            std::string err = "Error during computation of  the Cholesky factorization of a symmetric (Hermitian) positive-definite matrix using packed storage. matrix<T>::linsolve_chol(...). Error = " + std::to_string(info);
            throw std::string(err);
        }
        
        // 2. Solves a system of linear equations with a Cholesky-factored symmetric (compact) positive-definite coefficient matrix

        info = LAPACKE_dpptrs( matrix_order, 'L', row, col_b, _A, _B, col_b );

        if (info != 0)
        {
            std::string err = "Error during solving a linear system with Cholesky-factored general matrix. matrix<T>::linsolve_chol(...)" + std::to_string(info);
            throw std::string(err);
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::linsolve_chol(float *_A, MKL_INT rowA, MKL_INT colA, float *_B, MKL_INT colB)
    {
        // 1. First, make Cholesky factorizatioon of symmetric (compact) A: A_fact = L*L^T.

        lapack_int info = 0;
        lapack_int row = rowA;
        lapack_int col = colA;
        lapack_int col_b = colB;
        int matrix_order = LAPACK_ROW_MAJOR;
        
        info = LAPACKE_spptrf(matrix_order, 'L', col, _A);
        if (info != 0)
        {
            failbit = true;
            failinfo = info;
            std::string err = "Error during computation of  the Cholesky factorization of a symmetric (Hermitian) positive-definite matrix using packed storage. matrix<T>::linsolve_chol(...). Error = " + std::to_string(info);
            throw std::string(err);
        }
        
        // 2. Solves a system of linear equations with a Cholesky-factored symmetric (compact) positive-definite coefficient matrix

        info = LAPACKE_spptrs( matrix_order, 'L', row, col_b, _A, _B, col_b );

        if (info != 0)
        {
            std::string err = "Error during solving a linear system with Cholesky-factored general matrix. matrix<T>::linsolve_chol(...)" + std::to_string(info);
            throw std::string(err);
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::linsolve_chol_fact(double *_A, MKL_INT rowA, /*MKL_INT colA,*/ double *_B, MKL_INT colB)
    {
        // Solving linear system using already factorised A.

        lapack_int info = 0;
        lapack_int row = rowA;
        //lapack_int col = colA;
        lapack_int col_b = colB;
        int matrix_order = LAPACK_ROW_MAJOR;

        // Solves a system of linear equations with a Cholesky-factored symmetric (compact) positive-definite coefficient matrix

        info = LAPACKE_dpptrs( matrix_order, 'L', row, col_b, _A, _B, col_b );

        if (info != 0)
        {
            std::string err = "Error during solving a linear system with Cholesky-factored general matrix. matrix<T>::linsolve_chol_fact(...)" + std::to_string(info);
            throw std::string(err);
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::linsolve_chol_fact(float *_A, MKL_INT rowA, /*MKL_INT colA,*/ float *_B, MKL_INT colB)
    {
        // Solving linear system using already factorised A.

        lapack_int info = 0;
        lapack_int row = rowA;
        //lapack_int col = colA;
        lapack_int col_b = colB;
        int matrix_order = LAPACK_ROW_MAJOR;
        
        // Solves a system of linear equations with a Cholesky-factored symmetric (compact) positive-definite coefficient matrix

        info = LAPACKE_spptrs( matrix_order, 'L', row, col_b, _A, _B, col_b );

        if (info != 0)
        {
            std::string err = "Error during solving a linear system with Cholesky-factored general matrix. matrix<T>::linsolve_chol_fact(...)" + std::to_string(info);
            throw std::string(err);
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::inv_rec(double *_A, MKL_INT rowA, MKL_INT colA)
    {
        /*
            Matrix inversion. Interface to mkl routines.
            Here:
                A - matrix to invert;
                row - number of rows in A;
                col - number of columns in A;
                ipiv - empty vector of size lda.

            1. First, make LU factorizatioon of A: A_fact = P*L*U; ipiv = P.
            2. Second, make inversion of factorized A: A = inv(A_fact).
            The result of inversion is in matrix A.

            Return value: none.
        */

        lapack_int info = 0;
        lapack_int row = rowA;
        lapack_int col = colA;
        int matrix_order = LAPACK_ROW_MAJOR;

        lapack_int *ipiv;

#ifdef intelmkl
        ipiv = (lapack_int *)mkl_malloc(row * sizeof(lapack_int), sizeof(T) * 8);
#else
        ipiv = (lapack_int *)malloc(row * sizeof(lapack_int));
#endif
        if (ipiv == NULL)
        {
#ifdef intelmkl
            mkl_free(ipiv);
#else
            free(ipiv);
#endif
            failbit = true;
            throw std::string("Memory allocation error. matrix<T>::inv_rec(...)");
        }

        for (lapack_int i = 0; i < (row); i++)
            ipiv[i] = 1;

        info = LAPACKE_dgetrf(matrix_order, row, col, _A, col, ipiv);

        //info = LAPACKE_dsytrf( matrix_order, 'L', row, _A, col, ipiv ); // Bunch-Kaufman factorization
    
        if (info != 0)
        {
#ifdef intelmkl
            mkl_free(ipiv);
#else
            free(ipiv);
#endif
            throw std::string("Error during computation of the LU factorization of a general m-by-n matrix. matrix<T>::inv_rec(double *, MKL_INT, MKL_INT)");
        }

        info = LAPACKE_dgetri(matrix_order, row, _A, row, ipiv);

        //info = LAPACKE_dsytri ( matrix_order, 'L', row, _A, col, ipiv ); // due to Bunch-Kaufman factorization

        if (info != 0)
        {
#ifdef intelmkl
            mkl_free(ipiv);
#else
            free(ipiv);
#endif
            throw std::string("Error during computation the inverse of an LU-factored general matrix. matrix<T>::inv_rec(double *, MKL_INT, MKL_INT)");
        }
#ifdef intelmkl
        mkl_free(ipiv);
#else
        free(ipiv);
#endif
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::inv_rec(float *_A, MKL_INT rowA, MKL_INT colA)
    {
        /*
            Matrix inversion. Interface to mkl routines.
            Here:
                A - matrix to invert;
                row - number of rows in A;
                col - number of columns in A;
                ipiv - empty vector of size lda.

            1. First, make LU factorizatioon of A: A_fact = P*L*U; ipiv = P.
            2. Second, make inversion of factorized A: A = inv(A_fact).
            The result of inversion is in matrix A.

            Return value: none.
        */
        lapack_int info = 0;
        lapack_int row = rowA;
        lapack_int col = colA;
        int matrix_order = LAPACK_ROW_MAJOR;

        lapack_int *ipiv;

#ifdef intelmkl
        ipiv = (lapack_int *)mkl_malloc(row * sizeof(lapack_int), sizeof(T) * 8);
#else
        // ipiv = (lapack_int *)aligned_alloc(sizeof(T) * 8, row * sizeof(lapack_int));
        ipiv = (lapack_int *)malloc(row * sizeof(lapack_int));
#endif
        if (ipiv == NULL)
        {
#ifdef intelmkl
            mkl_free(ipiv);
#else
            free(ipiv);
#endif
            failbit = true;
            throw std::string("Memory allocation error. matrix<T>::inv_rec(...)");
        }
        for (lapack_int i = 0; i < (row); i++)
            ipiv[i] = 1;

        info = LAPACKE_sgetrf(matrix_order, row, col, _A, col, ipiv);

        if (info != 0)
        {
#ifdef intelmkl
            mkl_free(ipiv);
#else
            free(ipiv);
#endif
            throw std::string("Error during computation of the LU factorization of a general m-by-n matrix. matrix<T>::inv_rec(float *, MKL_INT, MKL_INT)");
        }

        info = LAPACKE_sgetri(matrix_order, row, _A, row, ipiv);

        if (info != 0)
        {
#ifdef intelmkl
            mkl_free(ipiv);
#else
            free(ipiv);
#endif
            throw std::string("Error during computation the inverse of an LU-factored general matrix. matrix<T>::inv_rec(float *, MKL_INT, MKL_INT)");
        }
#ifdef intelmkl
        mkl_free(ipiv);
#else
        free(ipiv);
#endif
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::inv_sym(double *_A, MKL_INT colA)
    {
        /*
            Symetric matrix inversion. Interface to mkl routines.

            Return value: none.
        */

        lapack_int info = 0;

        int matrix_order = LAPACK_ROW_MAJOR;

        info = LAPACKE_dpptrf(matrix_order, 'L', colA, _A);
        if (info != 0)
        {
            failbit = true;
            failinfo = (int)info;
            std::cerr<<"info output from LAPACKE_dpptrf: "<<info<<"\n";
            throw std::string("Error during computation of the Cholesky factorization of a symmetric (Hermitian) positive-definite matrix using packed storage. matrix<T>::inv_sym(double *, MKL_INT)");
        }

        info = LAPACKE_dpptri(matrix_order, 'L', colA, _A);
        if (info != 0)
        {
            failbit = true;
            failinfo = (int)info;
            throw std::string("Error during computation the inverse of a packed symmetric (Hermitian) positive-definite matrix after the Cholesky factorization. matrix<T>::inv_sym(double *, MKL_INT)");
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::inv_sym(float *_A, MKL_INT colA)
    {
        /*
            Symetric matrix inversion. Interface to mkl routines.

            Return value: none.
        */

        lapack_int info = 0;

        int matrix_order = LAPACK_ROW_MAJOR;

        info = LAPACKE_spptrf(matrix_order, 'L', colA, _A);

        if (info != 0)
        {
            std::cout<<"info = "<<info<<"\n";
            throw std::string("Error during computation of the Cholesky factorization of a symmetric (Hermitian) positive-definite matrix using packed storage. matrix<T>::inv_sym(float *, MKL_INT)");
        }

        info = LAPACKE_spptri(matrix_order, 'L', colA, _A);
        if (info != 0)
        {
            throw std::string("Error during computation the inverse of a packed symmetric (Hermitian) positive-definite matrix after the Cholesky factorization. matrix<T>::inv_sym(float *, MKL_INT)");
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::invert()
    {
        /*
            Symetric/General matrix inversion.

            Return value: none.

            Example:

                matrix <double> M(n,n);
                matrix <double> res;    // empty matrix

                for (auto i = 0; i < M.size(); i++)
                    M[i] = i;

                M.invert(); // now M is (n,n) inverted matrix
                res = M;       // now res is (n,n) inverted matrix

        */

        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::invert()");

        if (numRow == numCol)
        {
            if (!compact)
            {
                try
                {
                    inv_rec(A, numRow, numCol);
                }
                catch (std::string err)
                {
                    failbit = true;
                    throw err;
                }
            }
            else
            {
                try
                {
                    inv_sym(A, numCol);
                }
                catch (std::string err)
                {
                    failbit = true;
                    throw err;
                }
            }
        }
        else
        {
            failbit = true;
            throw std::string("Matrix is not square. matrix<T>::invert()");
        }
    }
    //===============================================================================================================

    template <typename T>
    void matrix<T>::lchol()
    {
        /*
            Symetric/General matrix Cholesky.

            Return value: none.

            Example:

                matrix <double> M(n,n);

                for (auto i = 0; i < M.size(); i++)
                    M[i] = i;

                M.lchol(); // now M is (n,n) lower triangular matrix, so L*L' -> original M

        */

        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::lchol()");

        if (numRow == numCol)
        {
            if (!compact)
            {
                try
                {
                    get_lchol(A, numRow, numCol);

                    // put zeros on the upper part of matrix:
                    for (size_t i = 0; i < numRow - 1; i++)
                    {
                        for (size_t j = i + 1; j < numCol; j++)
                            A[i * numRow + j] = 0.0;
                    }
                }
                catch (std::string err)
                {
                    failbit = true;
                    throw err;
                }
            }
            else
            {
                try
                {
                    get_lchol(A, numCol);
                }
                catch (std::string err)
                {
                    failbit = true;
                    throw err;
                }
            }
        }
        else
        {
            failbit = true;
            throw std::string("Matrix is not square. matrix<T>::lchol()");
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::get_lchol(double *_A, MKL_INT colA)
    {
        /*
            Cholesky decomposition of symmetric matrix in compact format. Interface to mkl routines.

            Return value: none.
        */

        lapack_int info = 0;

        int matrix_order = LAPACK_ROW_MAJOR;

        info = LAPACKE_dpptrf(matrix_order, 'L', colA, _A);
        if (info != 0)
        {
            failbit = true;
            failinfo = (int)info;
            throw std::string("Error during computation of  the Cholesky factorization of a symmetric (Hermitian) positive-definite matrix using packed storage. matrix<T>::get_lchol(double *, MKL_INT)");
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::get_lchol(float *_A, MKL_INT colA)
    {
        /*
            Cholesky decomposition of symmetric matrix in compact format. Interface to mkl routines.

            Return value: none.
        */

        lapack_int info = 0;

        int matrix_order = LAPACK_ROW_MAJOR;
        
        info = LAPACKE_spptrf(matrix_order, 'L', colA, _A);

        if (info != 0)
        {
            failbit = true;
            failinfo = (int)info;
            throw std::string("Error during computation of  the Cholesky factorization of a symmetric (Hermitian) positive-definite matrix using packed storage. matrix<T>::get_lchol(double *, MKL_INT)");
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::get_lchol(double *_A, MKL_INT colA, MKL_INT lda)
    {
        /*
            Cholesky decomposition of symmetric matrix in rectangular format. Interface to mkl routines.

            Return value: none.
        */

        lapack_int info = 0;

        int matrix_order = LAPACK_ROW_MAJOR;

        info = LAPACKE_dpotrf(matrix_order, 'L', colA, _A, lda);
        if (info != 0)
        {
            failbit = true;
            failinfo = (int)info;
            throw std::string("Error during computationof  the Cholesky factorization of a symmetric (Hermitian) positive-definite matrix using packed storage. matrix<T>::get_lchol(double *, MKL_INT, MKL_INT)");
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::get_lchol(float *_A, MKL_INT colA, MKL_INT lda)
    {
        /*
            Cholesky decomposition of symmetric matrix in rectangular format. Interface to mkl routines.

            Return value: none.
        */

        lapack_int info = 0;

        int matrix_order = LAPACK_ROW_MAJOR;

        info = LAPACKE_spotrf(matrix_order, 'L', colA, _A, lda);
        if (info != 0)
        {
            failbit = true;
            failinfo = (int)info;
            throw std::string("Error during computationof  the Cholesky factorization of a symmetric (Hermitian) positive-definite matrix using packed storage. matrix<T>::get_lchol(double *, MKL_INT, MKL_INT)");
        }
    }

    //===============================================================================================================

    template <typename T>
    size_t matrix<T>::capacity()
    {
        /*
            Returns number of allocated elements for A.
        */

        return resizedElements;
    }

    //===============================================================================================================

    template <typename T>
    matrix<T>::~matrix()
    {
        /*
            Class destructor.
        */

        if (allocated)
        {
#ifdef intelmkl
            mkl_free(A);
#else
            free(A);
#endif
            allocated = false;
        }

        if ( is_allocated_ipv_A )
        {
#ifdef intelmkl
            mkl_free(ipiv_A);
#else
            free(ipiv_A);
#endif
        }
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::print(std::string whiichMatrix)
    {
        /*
            Prints part of a matrix into a LOG file.

            Return value: none.
        */

        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::print(std::string)");

        int integer;
        size_t linteger;
        const std::type_info &ti1 = typeid(integer);
        const std::type_info &ti2 = typeid(linteger);
        const std::type_info &ti3 = typeid(A[0]);
        bool isInt = false;

        if (ti3 == ti1 || ti3 == ti2)
            isInt = true;

        FILE *dbgFile;
        dbgFile = fopen(debug_file.c_str(), "a");

        size_t maxRows = 50;
        fprintf(dbgFile, "%s", whiichMatrix.c_str());
        // fprintf (dbgFile, "\n");
        if (rectangular)
        {
            fprintf(dbgFile, "%s%s", ", Rectangular matrix, of type ", typeid(A[0]).name());
            fprintf(dbgFile, "\n\n");
            for (size_t i = 0; i < _min(maxRows, numRow); i++)
            {
                for (size_t j = 0; j < _min(maxRows, numCol); j++)
                {
                    if (isInt)
                        fprintf(dbgFile, "%12d", (int)A[i * numCol + j]);
                    else
                        fprintf(dbgFile, "%12.5G", (double)A[i * numCol + j]);
                }
                fprintf(dbgFile, "\n");
            }
        }
        else if (symetric)
        {
            fprintf(dbgFile, "%s%s", ", symetric matrix, of type ", typeid(A[0]).name());
            fprintf(dbgFile, "\n\n");
            for (size_t i = 0; i < _min(maxRows, numRow); i++)
            {
                for (size_t j = 0; j <= i; j++)
                {
                    if (isInt)
                        fprintf(dbgFile, "%12d", (int)A[i * (i + 1) / 2 + j]);
                    else
                        fprintf(dbgFile, "%12.5G", (double)A[i * (i + 1) / 2 + j]);
                }
                fprintf(dbgFile, "\n");
            }
        }
        fprintf(dbgFile, "\n");
        fprintf(dbgFile, "\n");

        fclose(dbgFile);
    }

    //===============================================================================================================

    template <typename T>
    void matrix<T>::printf(std::string whiichMatrix, bool append)
    {
        /*
            Prints all the matrix into the 'whiichMatrix' file.

            Return value: none.
        */
       
        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::printf(std::string)");

        int integer;
        size_t linteger;
        const std::type_info &ti1 = typeid(integer);
        const std::type_info &ti2 = typeid(linteger);
        const std::type_info &ti3 = typeid(A[0]);
        bool isInt = false;

        if (ti3 == ti1 || ti3 == ti2)
            isInt = true;

        FILE *dbgFile;
        if (append)
            dbgFile = fopen(whiichMatrix.c_str(), "a");
        else
            dbgFile = fopen(whiichMatrix.c_str(), "w");

        if (dbgFile == NULL)
            throw std::string("There is problem with opening file in void matrix<T>::printf(std::string)!");

        if (rectangular)
        {
            for (size_t i = 0; i < numRow; i++)
            {
                for (size_t j = 0; j < numCol; j++)
                {
                    if (isInt)
                        fprintf(dbgFile, "%12d", (int)A[i * numCol + j]);
                    else
                        fprintf(dbgFile, "%17.10G", (double)A[i * numCol + j]);
                }
                fprintf(dbgFile, "\n");
            }
        }
        else if (symetric)
        {
            for (size_t i = 0; i < numRow; i++)
            {
                for (size_t j = 0; j <= i; j++)
                {
                    if (isInt)
                        fprintf(dbgFile, "%12d", (int)A[i * (i + 1) / 2 + j]);
                    else
                        fprintf(dbgFile, "%17.10G", (double)A[i * (i + 1) / 2 + j]);
                }
                fprintf(dbgFile, "\n");
            }
        }

        fclose(dbgFile);
    }

} // namespace evolm

//===============================================================================================================

#endif // cs_matrix_hpp__
