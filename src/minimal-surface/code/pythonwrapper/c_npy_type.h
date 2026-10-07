#pragma once
#include <numpy/arrayobject.h>

template<typename CType>
struct npy_type{};

template<>
struct npy_type<int>{
    inline static const NPY_TYPES value = NPY_INT;
};

template<>
struct npy_type<double>{
    inline static const NPY_TYPES value = NPY_DOUBLE;
};

template<>
struct npy_type<float>{
    inline static const NPY_TYPES value = NPY_FLOAT;
};

template<>
struct npy_type<unsigned char>{
    inline static const NPY_TYPES value = NPY_UBYTE;
};


template<NPY_TYPES NPYType>
struct c_type {};

template<>
struct c_type<NPY_INT> {
    typedef int Type;
};

template<>
struct c_type<NPY_DOUBLE> {
    typedef double Type;
};

template<>
struct c_type<NPY_FLOAT> {
    typedef float Type;
};

template<>
struct c_type<NPY_UBYTE> {
    typedef unsigned char Type;
};



