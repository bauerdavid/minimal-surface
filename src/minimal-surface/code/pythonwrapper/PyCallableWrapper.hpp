#include <Python.h>
#include <iostream>
#include "python_utils.h"
#include <SimpleITK.h>
#include <vector>
#include <numpy/arrayobject.h>
#include <stdexcept>

namespace sitk = itk::simple;


template<typename RetVal, class ...Args>
class PyCallableWrapperBase {
protected:
    PyObject* callable;
public:
    typedef RetVal result_type;
    PyCallableWrapperBase(PyObject* obj): callable(obj){
		Py_INCREF(obj);
    };

    PyCallableWrapperBase(const PyCallableWrapperBase& other): callable(other.callable) {
     	Py_INCREF(callable);
    }

    PyCallableWrapperBase(PyCallableWrapperBase&& other): callable(other.callable) {
        Py_INCREF(callable);
        other.callable = 0;
    }

    PyCallableWrapperBase(): callable(nullptr) { }

    ~PyCallableWrapperBase() {
        Py_XDECREF(callable);
    }

    PyCallableWrapperBase& operator=(const PyCallableWrapperBase& other) {
        Py_XDECREF(callable);
        PyCallableWrapperBase tmp = other;
        *this = std::move(tmp);
        Py_INCREF(this->callable);
        return *this;
    }

    PyCallableWrapperBase& operator=(PyCallableWrapperBase&& other) {
        Py_XDECREF(callable);
        callable = other.callable;
        Py_INCREF(callable);
        other.callable = 0;
        return *this;
    }
    RetVal operator()(Args... args){
        if(callable)
            return this->call_pyobject(args...);
        throw std::runtime_error("Callable was not set!");
    }


    RetVal call_pyobject(Args... args) {
        // The GIL must be held for the whole call, including the conversion of
        // the return value by cast_from_python, and must be released on every
        // exit path - including the throws below.
        GILGuard gil;
        PyObject* arg_list = BuildArgs<Args...>(args...);
        if(arg_list == NULL){
            PyErr_Print();
            throw std::runtime_error("Failed to build arguments for Python callable.");
        }
        PyRef arg_list_ref(arg_list);
        PyObject* retval = PyObject_CallObject(callable, arg_list);
        if(retval == NULL){
            PyErr_Print();
            throw std::runtime_error("Failed to execute Python callable.");
        }
        PyRef retval_ref(retval);
        return cast_from_python<RetVal>(retval);
    }
};

template<class RetVal, class ...Args>
class PyCallableWrapper: public PyCallableWrapperBase<RetVal, Args...> {
    using PyCallableWrapperBase<RetVal, Args...>::PyCallableWrapperBase;
};


template<typename ...Args>
using CallbackWrapper = PyCallableWrapper<void, Args...>;

typedef CallbackWrapper<> CallbackWrapper_0A;

template<typename T>
using CallbackWrapper_1A = CallbackWrapper<T>;

template<typename T1, typename T2>
using CallbackWrapper_2A = CallbackWrapper<T1, T2>;

typedef PyCallableWrapper<sitk::Image, const sitk::Image&, const sitk::Image&> InitialContourCalculatorWrapper;

