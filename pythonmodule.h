#define PY_SSIZE_T_CLEAN
#pragma once
#include <Python.h>
#include <numpy/arrayobject.h>
#include "MinimalSurfaceEstimator.h"
#include "SimpleITK.h"
#include "structmember.h"
#include <iostream>
#include <fstream>

static PyObject* compute_minimal_surface(PyObject* self, PyObject* args);

PyMODINIT_FUNC PyInit_MinArea(void);