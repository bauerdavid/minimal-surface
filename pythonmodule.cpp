#include "pythonmodule.h"
using namespace std;

#define PYDEFINE_CONST(module_, macro) PyModule_AddIntConstant(module_, #macro, macro)
namespace sitk = itk::simple;

typedef struct {
	PyObject_HEAD
	MinimalSurfaceEstimator* mMinimalSurfaceEstimator;
} Estimator;

static void
Estimator_dealloc(Estimator* self)
{
	delete self->mMinimalSurfaceEstimator;
	Py_TYPE(self)->tp_free((PyObject*)self);
}


static PyObject* Estimator_new(PyTypeObject* type, PyObject* args, PyObject* kwargs) {
	Estimator* self;
	self = (Estimator*)type->tp_alloc(type, 0);

	if (self != NULL) {
		self->mMinimalSurfaceEstimator = new MinimalSurfaceEstimator();
	}
	return (PyObject*)self;
}

static int Estimator_init(Estimator* self, PyObject* args, PyObject* kwargs) {
	return 0;
}

static PyMemberDef Estimator_members[] = {
	//{"first", T_OBJECT_EX, offsetof(Estimator, mMinimalSurfaceEstimator), 0, "first name"},
	{NULL}
};

static PyObject* Estimator_HookStageInitializationEvent(Estimator* self, PyObject* args) {
	int stage;
	PyObject* callback;
	if(!PyArg_ParseTuple(args, "iO", &stage, &callback))
		return NULL;
	if (!PyCallable_Check(callback))
		return NULL;
	Py_INCREF(callback);
	EventSource::basic_callback_type handler = [callback, stage]() {
		PyGILState_STATE gstate;
		gstate = PyGILState_Ensure();
		static PyObject* cb = callback;
		PyObject_CallNoArgs(callback);
		PyEval_CallFunction(callback, "()");
		PyGILState_Release(gstate);
	};
	self->mMinimalSurfaceEstimator->HookStageInitializationEvent(eStageEnum(stage), handler);
	Py_RETURN_NONE;
}

static PyObject* Estimator_HookStageFinishedEvent(Estimator* self, PyObject* args) {
	int stage;
	PyObject* callback;
	if (!PyArg_ParseTuple(args, "iO", &stage, &callback))
		return NULL;
	if (!PyCallable_Check(callback))
		return NULL;
	Py_INCREF(callback);
	EventSource::basic_callback_type handler = [callback]() {
		PyGILState_STATE gstate;
		gstate = PyGILState_Ensure(); 
		PyObject_CallNoArgs(callback);
		PyGILState_Release(gstate);
	};
	self->mMinimalSurfaceEstimator->HookStageFinishedEvent(eStageEnum(stage), handler);
	Py_RETURN_NONE;
}

static PyObject* Estimator_HookStageIterationEvent(Estimator* self, PyObject* args) {
	int stage;
	PyObject* callback;
	if (!PyArg_ParseTuple(args, "iO", &stage, &callback))
		return NULL;
	if (!PyCallable_Check(callback))
		return NULL;
	Py_INCREF(callback);
	EventSource::iter_callback_type handler = [callback](int iteration) {
		PyGILState_STATE gstate;
		gstate = PyGILState_Ensure(); 
		PyObject* args = Py_BuildValue("(i)", iteration);
		PyObject_Call(callback, args, NULL);
		Py_DECREF(args);
		PyGILState_Release(gstate);
	};
	self->mMinimalSurfaceEstimator->HookStageIterationEvent(eStageEnum(stage), handler);
	Py_RETURN_NONE;
}

static PyObject* Estimator_HookStageDataInitializedEvent(Estimator* self, PyObject* args) {
	int stage;
	PyObject* callback;
	int data_type = 0;
	if (!PyArg_ParseTuple(args, "iO|i", &stage, &callback, &data_type)) {
		PyErr_Format(PyExc_TypeError, "Couldn't parse arguments in %s", __FUNCTION__);
		return NULL;
	}
	if (!PyCallable_Check(callback)) {
		PyErr_Format(PyExc_TypeError, "Second argument is not callable");
		return NULL;
	}
	Py_INCREF(callback);
	EventSource::data_callback_type handler = [callback](sitk::Image& map, int idx) {
		vector<unsigned> size = map.GetSize();
		int ndims = size.size();
		if (ndims > 0) {
			cout << size[0] << ", " << size[1];
			if (ndims > 2)
				cout << ", " << size[2];
			cout << endl;
		}
		else {
			cout << "empty image" << endl;
		}int n_pixels = 1;
		npy_intp* image_dims = (npy_intp*)malloc(3 * sizeof(npy_intp));
		if (image_dims == NULL) {
			PyErr_NoMemory();
			return;
		}
		for (int i = 0; i < ndims; i++) {
			image_dims[i] = size[ndims - i - 1];
			n_pixels *= image_dims[i];
		}
		cout << "ndims: " << ndims << endl;
		if (ndims < 3)
			image_dims[2] = 1;
		double* map_buffer = map.GetBufferAsDouble();
		cout << "got buffer" << endl;
		PyGILState_STATE gstate;
		gstate = PyGILState_Ensure(); 
		PyObject* out = PyArray_SimpleNewFromData(3, image_dims, NPY_DOUBLE, map_buffer);
		if (out == NULL) {
			cout << "Couldn't create numpy array" << endl;
		}
		PyObject* maxval = PyArray_Max((PyArrayObject*)out, NPY_MAXDIMS, NULL);
		double val;
		PyArray_ScalarAsCtype(maxval, &val);
		cout << "max: " << val << endl;
		Py_DECREF(maxval);
		//PyObject* out = Py_None; Py_INCREF(Py_None);
		PyObject* args = Py_BuildValue("(Oi)", out, idx);
		if (args == NULL) {
			cout << "args were not built" << endl;
			goto closure;
		}
		if (PyObject_Call(callback, args, NULL) == NULL) {
			cout << "There was an error" << endl;
		}
closure:
		Py_XDECREF(args);
		free(image_dims);
		Py_XDECREF(out);
		PyGILState_Release(gstate);
	};
	self->mMinimalSurfaceEstimator->HookStageDataInitializedEvent(eStageEnum(stage), handler, eDataEnum(data_type));
	Py_RETURN_NONE;
}

static PyObject* Estimator_HookStageUpdatedEvent(Estimator* self, PyObject* args) {
	int stage;
	PyObject* callback;
	int data_type = 0;
	if (!PyArg_ParseTuple(args, "iO|i", &stage, &callback, &data_type)) {
		PyErr_Format(PyExc_TypeError, "Couldn't parse arguments in %s", __FUNCTION__);
		return NULL;
	}
	if (!PyCallable_Check(callback)) {
		PyErr_Format(PyExc_TypeError, "Second argument is not callable");
		return NULL;
	}
	Py_INCREF(callback);
	EventSource::data_callback_type handler = [callback](sitk::Image& map, int idx) {
		vector<unsigned> size = map.GetSize();
		if (size.size() > 0) {
			cout << size[0] << ", " << size[1];
			if (size.size() > 2)
				cout << ", " << size[2];
			cout << endl;
		}
		else {
			cout << "empty image" << endl;
		}
		int ndims = size.size();

		int n_pixels = 1;
		npy_intp image_dims[3];
		for (int i = 0; i < ndims; i++) {
			image_dims[i] = size[ndims - i - 1];
			n_pixels *= image_dims[i];
		}
		if (ndims < 3)
			image_dims[2] = 1; 
		double* map_buffer = map.GetBufferAsDouble();
		PyGILState_STATE gstate;
		gstate = PyGILState_Ensure(); 
		PyObject* out = PyArray_SimpleNewFromData(3, image_dims, NPY_DOUBLE, map_buffer);
		PyObject* args = Py_BuildValue("(Oi)", out, idx);
		if (PyObject_Call(callback, args, NULL) == NULL) {
			PyErr_Print();
		}
		Py_DECREF(args);
		PyGILState_Release(gstate);
	};
	EventSource::data_callback_type handler2 = [callback](sitk::Image& map, int idx) {
		cout << "HANDLER2" << endl;
		PyGILState_STATE gstate = PyGILState_Ensure();
		if (PyObject_CallFunction(callback, "(i)", idx) == NULL) {
		}
		PyGILState_Release(gstate);
	};
	self->mMinimalSurfaceEstimator->HookStageUpdatedEvent(eStageEnum(stage), handler, eDataEnum(data_type));
	Py_RETURN_NONE;
}

static PyObject* Estimator_HookImageTransformCalculatedEvent(Estimator* self, PyObject* args) {
	PyObject* callback;
	if (!PyArg_ParseTuple(args, "O", &callback))
		return NULL;
	if (!PyCallable_Check(callback))
		return NULL;
	Py_INCREF(callback);
	MinimalSurfaceEstimator::image_transform_callback_type handler = [callback](const std::vector<double>& rotation, const std::vector<double>& translation) {
		int ndims = 1;
		npy_intp rot_mat_dims = 9;
		double* rot_mat_buffer = (double*)malloc(9 * sizeof(double));
		std::copy(rotation.begin(), rotation.end(), rot_mat_buffer);

		npy_intp transl_mat_dims = 3;
		double* transl_mat_buffer = (double*)malloc(3 * sizeof(double));
		std::copy(translation.begin(), translation.end(), transl_mat_buffer);
		PyGILState_STATE gstate;
		gstate = PyGILState_Ensure();
		PyObject* rotation_obj = PyArray_SimpleNewFromData(ndims, &rot_mat_dims, NPY_DOUBLE, rot_mat_buffer);
		PyObject* translation_obj = PyArray_SimpleNewFromData(ndims, &transl_mat_dims, NPY_DOUBLE, transl_mat_buffer);
		PyObject* args = Py_BuildValue("(OO)", rotation_obj, translation_obj);
		PyObject_Call(callback, args, NULL);
		Py_DECREF(args);
		Py_DECREF(rotation_obj);
		Py_DECREF(translation_obj);
		PyGILState_Release(gstate);
		free(rot_mat_buffer);
		free(transl_mat_buffer);
	};
	self->mMinimalSurfaceEstimator->HookImageTransformCalculatedEvent(handler);
	Py_RETURN_NONE;
}

static PyObject* Estimator_HookPlaneCenterCalculatedEvent(Estimator* self, PyObject* args) {
	PyObject* callback;
	if (!PyArg_ParseTuple(args, "O", &callback))
		return NULL;
	if (!PyCallable_Check(callback))
		return NULL;
	Py_INCREF(callback);
	MinimalSurfaceEstimator::vector_callback_type handler = [callback](const std::vector<double>& center) {
		int ndims = 1;
		npy_intp center_mat_dims = 3;
		double* center_mat_buffer = (double*)malloc(3 * sizeof(double));
		std::copy(center.begin(), center.end(), center_mat_buffer);

		PyGILState_STATE gstate;
		gstate = PyGILState_Ensure();
		PyObject* center_obj = PyArray_SimpleNewFromData(ndims, &center_mat_dims, NPY_DOUBLE, center_mat_buffer);
		PyObject* args = Py_BuildValue("(O)", center_obj);
		PyObject_Call(callback, args, NULL);
		Py_DECREF(args);
		Py_DECREF(center_obj);
		PyGILState_Release(gstate);
		free(center_mat_buffer);
	};
	self->mMinimalSurfaceEstimator->HookPlaneCenterCalculatedEvent(handler);
	Py_RETURN_NONE;
}

static PyObject* Estimator_ResolvePath(Estimator* self, PyObject* args) {
	if (PyArray_API == NULL)
	{
		import_array();
	}
	PyObject* point_o;
	PyObject* data_o;
	if (!PyArg_ParseTuple(args, "OO", &point_o, &data_o)) {
		return NULL;
	}
	PyObject* point_array;
	if ((point_array = PyArray_FROM_OTF(point_o, NPY_INT, NPY_IN_ARRAY)) == NULL) {
		return NULL;
	}
	int point_ndim = PyArray_NDIM(point_array);
	if (point_ndim != 1)
		return NULL;
	npy_intp* point_dims = PyArray_DIMS(point_array);
	if (*point_dims != 3)
		return NULL;
	int* point_buffer = (int*)PyArray_DATA(point_array);
	Vec3<int> point = { point_buffer[0], point_buffer[1], point_buffer[2] };
	
	PyObject* data_array;
	if ((data_array = PyArray_FROM_OTF(data_o, NPY_DOUBLE, NPY_IN_ARRAY)) == NULL) {
		return NULL;
	}
	int data_ndim = PyArray_NDIM(data_array);
	if (data_ndim != 3)
		return NULL;
	npy_intp* data_dims = PyArray_DIMS(data_array);
	double* data_buffer = (double*)PyArray_DATA(data_array);
	sitk::Image data({ (unsigned)data_dims[2], (unsigned)data_dims[1], (unsigned)data_dims[0] }, sitk::sitkFloat64);
	memcpy(data.GetBufferAsDouble(), data_buffer, data_dims[2] * data_dims[1] * data_dims[0] * sizeof(double));
	vector<Vec3<int>> path;
	Py_BEGIN_ALLOW_THREADS
	path = ResolvePath(point, data);
	Py_END_ALLOW_THREADS
	npy_intp path_dims[] = { path.size(), 3 };
	PyObject* path_o = PyArray_SimpleNew(2, path_dims, NPY_INT);
	int* path_buffer = (int*)PyArray_DATA(path_o);
	for (int i = 0; i < path.size(); i++) {
		Vec3<int> point = path[i];
		path_buffer[i * 3] = point.x();
		path_buffer[i * 3 + 1] = point.y();
		path_buffer[i * 3 + 2] = point.z();
	}
	return path_o;
}

static PyObject* Estimator_calculate(Estimator* self, PyObject* args) {
	if (PyArray_API == NULL)
	{
		import_array();
	}
	PyObject* image_obj;
	PyObject* point1_obj;
	PyObject* point2_obj;
	int use_correction;
	double beta;
	double alpha;
	if (!PyArg_ParseTuple(args, "OOOpdd", &image_obj, &point1_obj, &point2_obj, &use_correction, &beta, &alpha)) {
		return NULL;
	}
	//initialize image info
	PyObject* image_array;
	if ((image_array = PyArray_FROM_OTF(image_obj, NPY_DOUBLE, NPY_IN_ARRAY)) == NULL) {
		return NULL;
	}
	int image_ndim = PyArray_NDIM(image_array);
	npy_intp* image_dims = PyArray_DIMS(image_array);
	double* image_data = (double*)PyArray_DATA(image_array);
	sitk::Image image({ (unsigned)image_dims[2], (unsigned)image_dims[1], (unsigned)image_dims[0] }, sitk::sitkFloat64);
	memcpy(image.GetBufferAsDouble(), image_data, image_dims[2] * image_dims[1] * image_dims[0] * sizeof(double));
	//initialize points info
	PyObject* point1_array;
	if ((point1_array = PyArray_FROM_OTF(point1_obj, NPY_DOUBLE, NPY_IN_ARRAY)) == NULL) {
		return NULL;
	}
	int point1_ndim = PyArray_NDIM(point1_array);
	npy_intp* point1_dims = PyArray_DIMS(point1_array);
	double* point1_data = (double*)PyArray_DATA(point1_array);
	Vec3<double> point1 = { point1_data[0], point1_data[1], point1_data[2] };

	PyObject* point2_array;
	if ((point2_array = PyArray_FROM_OTF(point2_obj, NPY_DOUBLE, NPY_IN_ARRAY)) == NULL) {
		return NULL;
	}
	int point2_ndim = PyArray_NDIM(point2_array);
	npy_intp* point2_dims = PyArray_DIMS(point2_array);
	double* point2_data = (double*)PyArray_DATA(point2_array);
	Vec3<double> point2 = { point2_data[0], point2_data[1], point2_data[2] };

	self->mMinimalSurfaceEstimator->SetUsesCorrection(use_correction);
	Py_BEGIN_ALLOW_THREADS
	self->mMinimalSurfaceEstimator->Calculate(image, point1, point2, beta, alpha);
	Py_END_ALLOW_THREADS
	const double* transport_buffer = self->mMinimalSurfaceEstimator->GetTransportFunctionCalculator().GetTransportFunction().GetBufferAsDouble();
	PyObject* out = PyArray_SimpleNew(3, image_dims, NPY_DOUBLE);
	double* outdata = (double*)PyArray_DATA(out);
	memcpy(outdata, transport_buffer, image_dims[0] * image_dims[1] * image_dims[2] * sizeof(double));
	return out;
}

static PyMethodDef Estimator_methods[] = {
	{"hook_stage_init_event", (PyCFunction)Estimator_HookStageInitializationEvent, METH_VARARGS, ""},
	{"hook_stage_finished_event", (PyCFunction)Estimator_HookStageFinishedEvent, METH_VARARGS, ""},
	{"hook_stage_iteration_event", (PyCFunction)Estimator_HookStageIterationEvent, METH_VARARGS, ""},
	{"hook_stage_updated_event", (PyCFunction)Estimator_HookStageUpdatedEvent, METH_VARARGS, ""},
	{"hook_stage_data_init_event", (PyCFunction)Estimator_HookStageDataInitializedEvent, METH_VARARGS, ""},
	{"hook_transform_calculated_event", (PyCFunction)Estimator_HookImageTransformCalculatedEvent, METH_VARARGS, ""},
	{"hook_plane_center_calculated_event", (PyCFunction)Estimator_HookPlaneCenterCalculatedEvent, METH_VARARGS, ""},
	{"resolve_shortest_paths", (PyCFunction)Estimator_ResolvePath, METH_VARARGS, ""},
	{"calculate", (PyCFunction)Estimator_calculate, METH_VARARGS, ""},
	{nullptr, nullptr, 0, nullptr}/* Sentinel */
};

static PyTypeObject MinArea_EstimatorType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	"MinArea.MinimalSurfaceEstimator",             /* tp_name */
	sizeof(Estimator), /* tp_basicsize */
	0,                         /* tp_itemsize */
	(destructor)Estimator_dealloc,                         /* tp_dealloc */
	0,                         /* tp_print */
	0,                         /* tp_getattr */
	0,                         /* tp_setattr */
	0,                         /* tp_compare */
	0,                         /* tp_repr */
	0,                         /* tp_as_number */
	0,                         /* tp_as_sequence */
	0,                         /* tp_as_mapping */
	0,                         /* tp_hash */
	0,                         /* tp_call */
	0,                         /* tp_str */
	0,                         /* tp_getattro */
	0,                         /* tp_setattro */
	0,                         /* tp_as_buffer */
	Py_TPFLAGS_DEFAULT,        /* tp_flags */
	"This will be the documentation of MinimalSurfaceEstimator",           /* tp_doc */
	0,                         /* tp_traverse */
	0,                         /* tp_clear */
	0,                         /* tp_richcompare */
	0,                         /* tp_weaklistoffset */
	0,                         /* tp_iter */
	0,                         /* tp_iternext */
	Estimator_methods,             /* tp_methods */
	Estimator_members,             /* tp_members */
	0,                         /* tp_getset */
	0,                         /* tp_base */
	0,                         /* tp_dict */
	0,                         /* tp_descr_get */
	0,                         /* tp_descr_set */
	0,                         /* tp_dictoffset */
	(initproc)Estimator_init,      /* tp_init */
	0,                         /* tp_alloc */
	Estimator_new,                 /* tp_new */
};


static PyMethodDef minimalsurface_methods[] = {
	{nullptr, nullptr, 0, nullptr}/* Sentinel */
};

static PyModuleDef minimalsurface_module = {
	PyModuleDef_HEAD_INIT,
	"MinArea",
	"Some random documentation",
	100,
	minimalsurface_methods
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
PyInit_MinArea(void)
{
	PyObject* m;

	if (PyType_Ready(&MinArea_EstimatorType) < 0)
		return NULL;

	m = PyModule_Create(&minimalsurface_module);
	if (m == NULL)
		return NULL;
	PYDEFINE_CONST(m, AREA_EIKONAL_STAGE);
	PYDEFINE_CONST(m, ROTATED_AREA_EIKONAL_STAGE);
	PYDEFINE_CONST(m, PLANE_PHASEFIELD_STAGE);
	PYDEFINE_CONST(m, TRANSPORT_FUNCTION_STAGE);

	PYDEFINE_CONST(m, DEFAULT_DATA);
	PYDEFINE_CONST(m, DISTANCE_DATA);
	PYDEFINE_CONST(m, PHASEFIELD_DATA);
	PYDEFINE_CONST(m, MEETING_POINTS_DATA);
	PYDEFINE_CONST(m, CURVATURE_DATA);
	Py_INCREF(&MinArea_EstimatorType);
	PyModule_AddObject(m, "Estimator", (PyObject*)&MinArea_EstimatorType);
	return m;
}