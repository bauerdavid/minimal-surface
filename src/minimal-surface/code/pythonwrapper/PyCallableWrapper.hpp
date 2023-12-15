#include <Python.h>
#include <iostream>
#include "python_utils.h"
#include <SimpleITK.h>
#include <vector>
#include <numpy/arrayobject.h>
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
        *this = std::move(tmp)
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
        PyGILState_STATE gstate;
		gstate = PyGILState_Ensure();
        PyObject* arg_list = BuildArgs<Args...>(args...);
        if(PyErr_Occurred() != NULL){
            PyErr_Print();
        }
        PyObject* retval = PyObject_CallObject(callable, arg_list);
        if(PyErr_Occurred() != NULL){
            PyErr_Print();
        }
        PyGILState_Release(gstate);
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

