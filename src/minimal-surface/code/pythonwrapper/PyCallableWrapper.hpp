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
    // Takes a new reference to obj, which may be NULL.
    PyCallableWrapperBase(PyObject* obj): callable(obj){
        Py_XINCREF(callable);
    };

    PyCallableWrapperBase(const PyCallableWrapperBase& other): callable(other.callable) {
        Py_XINCREF(callable);
    }

    // Steals the reference rather than taking a new one. Taking a reference
    // here and then clearing other.callable - so that the moved-from object
    // no longer releases anything - leaked one reference per move.
    PyCallableWrapperBase(PyCallableWrapperBase&& other) noexcept: callable(other.callable) {
        other.callable = nullptr;
    }

    PyCallableWrapperBase(): callable(nullptr) { }

    ~PyCallableWrapperBase() {
        Py_XDECREF(callable);
    }

    // Acquires the new reference before releasing the old one, so that
    // self-assignment through an alias cannot free the object being assigned.
    // The previous body released this->callable and then let the delegated
    // move-assignment release it a second time.
    PyCallableWrapperBase& operator=(const PyCallableWrapperBase& other) {
        PyObject* old = callable;
        callable = other.callable;
        Py_XINCREF(callable);
        Py_XDECREF(old);
        return *this;
    }

    PyCallableWrapperBase& operator=(PyCallableWrapperBase&& other) noexcept {
        if (this != &other) {
            Py_XDECREF(callable);
            callable = other.callable;
            other.callable = nullptr;
        }
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

