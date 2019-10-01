/******************************************************************************
*                              vector_tools.h                                 *
===============================================================================
* A collection of functions and related utilities intended to help perform    *
* vector operations in cpp.                                                   *
*******************************************************************************
*/

#ifndef VECTOR_TOOLS_H
#define VECTOR_TOOLS_H

#include<stdio.h>
#include<iostream>
#include<exception>
#include<fstream>
#include<vector>
#include<map>
#include<math.h>
#include<assert.h>
#include<string.h>


//Type definitions
typedef unsigned int size_type;

//Operator overloading
template<typename T>
std::vector<T>& operator+=(std::vector<T> &lhs, const std::vector<T> &rhs);

template<typename T>
std::vector<T> operator+(std::vector<T> lhs, const std::vector<T> &rhs);

template<typename T>
std::vector<T> operator-(std::vector<T> v);

template<typename T>
std::vector<T>& operator-=(std::vector<T> &lhs, const std::vector<T> &rhs);

template<typename T>
std::vector<T> operator-(std::vector<T> lhs, const std::vector<T> &rhs);

template<typename T, typename t>
std::vector<T>& operator*=(std::vector<T> &lhs, const t rhs);

template<typename T, typename t>
std::vector<T> operator*(const t lhs, std::vector<T> rhs);

template<typename T, typename t>
std::vector<T> operator*(std::vector<T> lhs, const t rhs);

template<typename T, typename t>
std::vector<T>& operator/=(std::vector<T> &lhs, const t rhs);

template<typename T, typename t>
std::vector<T> operator/(std::vector<T> lhs, const t rhs);

namespace vectorTools{

    //Computation Utilities
    template<typename T>
    int computeMean(const std::vector<std::vector< T > > &A, std::vector< T > &v);

    template<typename T>
    int cross(const std::vector< T > &a, const std::vector< T > &b, std::vector< T > &c);

    template<typename T>
    std::vector< T > cross(const std::vector< T > &a, const std::vector< T > &b);

    template<typename T>
    int dot(const std::vector< T > &a, const std::vector< T > &b, T &c);

    template<typename T>
    T dot(const std::vector< T > &a, const std::vector< T > &b);

    template<typename T>
    std::vector< T > dot(const std::vector< std::vector< T > > &A, const std::vector< T > &b);

    template<typename T>
    double l2norm(const std::vector< T > &v);

    template<typename T>
    double l2norm(const std::vector< std::vector< T > > &A);

    //Printing Utilities
    template<typename T>
    int print(std::vector< T > &v);

}

#include "vector_tools.cpp"
#endif
