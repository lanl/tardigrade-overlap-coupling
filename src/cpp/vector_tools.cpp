/******************************************************************************
*                               vector_tools.cpp                              *
===============================================================================
* A collection of functions and related utilities intended to help perform    *
* vector operations in cpp.                                                   *
*******************************************************************************
*/

#include "vector_tools.h"


//Operator overloading
template<typename T>
std::vector<T>& operator+=(std::vector<T> &lhs, const std::vector<T> &rhs){
    /*!
     * Overload the += operator for vectors
     * 
     * :param std::vector<T> &lhs: The left-hand side vector
     * :param std::vector<T> &rhs: The right-hand side vector
     */

    if (lhs.size() != rhs.size()){
        throw std::length_error("vectors must be the same size to add");
    }

    for (size_type i=0; i<lhs.size(); i++){
        lhs[i] += rhs[i];
    }
    return lhs;
}

template<typename T>
std::vector<T> operator+(std::vector<T> lhs, const std::vector<T> &rhs){
    /*!
     * Overload the + operator for vectors
     * 
     * :param std::vector<T> &lhs: The left-hand side vector
     * :param std::vector<T> &rhs: The right-hand side vector
     */

    if (lhs.size() != rhs.size()){
        throw std::length_error("vectors must be the same size to add");
    }
    return lhs += rhs;
}

template<typename T>
std::vector<T> operator-(std::vector<T> v){
    /*!
     * Overload the negative operator for vectors
     * 
     * :param std::vector<T> &v: The vector in question
     */

    for (size_type i=0; i<v.size(); i++){
        v[i] = -v[i];
    }
    return v;
}

template<typename T>
std::vector<T>& operator-=(std::vector<T> &lhs, const std::vector<T> &rhs){
    /*!
     * Overload the -= operator for vectors
     * 
     * :param std::vector<T> &lhs: The left-hand side vector
     * :param std::vector<T> &rhs: The right-hand side vector
     */
    return lhs += -rhs;
}

template<typename T>
std::vector<T> operator-(std::vector<T> lhs, const std::vector<T> &rhs){
    /*!
     * Overload the subtraction operator for vectors
     * 
     * :param std::vector<T> &lhs: The left-hand side vector
     * :param std::vector<T> &rhs: The right-hand side vector
     */
    
    return lhs -= rhs;
}


template<typename T, typename t>
std::vector<T>& operator*=(std::vector<T> &lhs, const t rhs){
    /*!
     * Overload the *= operator for vectors
     * 
     * :param std::vector<T> lhs: The left-hand side vector
     * :param const t rhs: The right-hand side scalar
     */
    for (size_type i=0; i<lhs.size(); i++){
        lhs[i] *= rhs;
    }
    return lhs;
}

template<typename T, typename t>
std::vector<T> operator*(const t lhs, std::vector<T> rhs){
    /*!
     * Overload the / operator for vectors
     * 
     * :param const t lhs: The left-hand side scalar
     * :param std::vector<T> rhs: The right-hand side vector
     */
    return rhs*=lhs;
}


template<typename T, typename t>
std::vector<T> operator*(std::vector<T> lhs, const t rhs){
    /*!
     * Overload the / operator for vectors
     * 
     * :param std::vector<T> lhs: The left-hand side vector
     * :param const t rhs: The right-hand side scalar
     */
    return lhs*=rhs;
}

template<typename T, typename t>
std::vector<T>& operator/=(std::vector<T> &lhs, const t rhs){
    /*!
     * Overload the /= operator for vectors
     * 
     * :param std::vector<T> lhs: The left-hand side vector
     * :param const t rhs: The right-hand side scalar
     */
    return lhs*=(1./rhs);
}

template<typename T, typename t>
std::vector<T> operator/(std::vector<T> lhs, const t rhs){
    /*!
     * Overload the / operator for vectors
     * 
     * :param std::vector<T> lhs: The left-hand side vector
     * :param const t rhs: The right-hand side scalar
     */
    return lhs/=rhs;
}

namespace vectorTools{

    //Computation Utilities
    template<typename T>
    int computeMean(const std::vector< std::vector< T > > &A, std::vector< T > &v){
        /*!
         * Compute the column-wise mean of A
         * 
         * :param matrixType &A: The matrix of vectors
         * :param vectorType &v: The resulting mean
         */

        if (A.size() == 0){
            std::cerr << "Error: Matrix must have a size greater than zero\n";
            return 1;
        }

        //Size the output vector
        v = std::vector<T>(A[0].size(), 0);

        //Compute the mean
        for (auto it = A.begin(); it!=A.end(); it++){
            v += *it/A.size();
        }
        return 0;
    }

    template<typename T>
    int cross(const std::vector< T > &a, const std::vector< T > &b, std::vector< T > &c){
        /*!
         * Compute the cross product of two vectors i.e. a x b
         * Note that if a and b are 2D vectors a 3D vector for c will be returned.
         * 
         * TODO: Generalize this to n dimensions.
         * 
         * :param std::vector< T > &a: The first vector
         * :param std::vector< T > &b: The second vector
         */

        size_type size = a.size();
        c = std::vector< T >(size, 0);

        if (size == 2){
            c.resize(3);
            c[2] =  a[0]*b[1] - a[1]*b[0];
        }
        else if (size == 3){
            c[0] =  a[1]*b[2] - a[2]*b[1];
            c[1] = -a[0]*b[2] + a[2]*b[0];
            c[2] =  a[0]*b[1] - a[1]*b[0];
        }
        else{
            throw std::length_error("Only 2D and 3D vectors are accepted");
        }

        return 0;
    }

    template<typename T>
    std::vector< T > cross(const std::vector< T > &a, const std::vector< T > &b){
        /*!
         * Compute the cross product of two vectors i.e. a x b
         * Note that if a and b are 2D vectors a 3D vector for c will be returned
         * 
         * TODO: Generalize this to n dimensions
         * 
         * :param std::vector< T > &a: The first vector
         * :param std::vector< T > &b: The second vector
         */

         std::vector< T > c;
         cross(a, b, c);
         return c;
    }

    template<typename T>
    int dot(const std::vector< T > &a, const std::vector< T > &b, T &v){
        /*!
         * Compute the dot product of two vectors i.e. v = a_i b_i
         * 
         * :param std::vector< T > &a: The first vector
         * :param std::vector< T > &b: The second vector
         * :param T &v: The output quantity
         */

        //Get the size and perform error handing
        size_type size = a.size();
        if (size != b.size()){
            throw std::length_error("vectors must be the same size to add");
        }

        //Set v to null
        v = 0;

        for (size_type i=0; i<size; i++){
            v += a[i]*b[i];
        }
        return 0;
    }

    template<typename T>
    T dot(const std::vector< T > &a, const std::vector< T > &b){
        /*!
         * Compute the dot product of two vectors i.e. v = a_i b_i
         * 
         * :param std::vector< T > &a: The first vector
         * :param std::vector< T > &b: The second vector
         */

        T v;
        dot(a, b, v);
        return v;
    }

    template<typename T>
    std::vector< T > dot(const std::vector< std::vector< T > > &A, const std::vector< T > &b){
        /*!
         * Compute the dot product between a matrix and a vector resulting i.e. c_i = A_ij b_j
         * 
         * :param std::vector< std::vector< T > > &A: The matrix
         * :param std::vector< T > &b: The vector
         */

        size_type size = A.size();

        std::vector< T > c(size);

        unsigned int i=0;
        for (auto A_i=A.begin(); A_i!=A.end(); A_i++, i++){
            c[i] = dot(*A_i, b);
        }
        return c;
    }

    //Printing Utilities
    template<typename T>
    int print(std::vector< T > &v){
        /*!
         * Print the contents of the vector to the terminal assuming << has been defined for each component
         * 
         * :param std::vector< T > &v: The vector to be displayed
         */

        for (auto it = v.begin(); it!=v.end(); it++){
            std::cout << *it << " ";
        }
        std::cout << "\n";
        return 1;
    }
}
