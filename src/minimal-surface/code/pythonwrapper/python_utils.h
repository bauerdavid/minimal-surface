#include <Python.h>
#include "arg_format_string.h"
#include "c_npy_type.h"
#include <SimpleITK.h>
#include <numpy/arrayobject.h>
#include <vector>
#include <type_traits>
#include <typeinfo>
#ifndef _MSC_VER
#include <cxxabi.h>
#endif
#include <memory>
#include <string>
#include <cstdlib>

namespace sitk = itk::simple;


template <class T>
std::string
type_name()
{
    typedef typename std::remove_reference<T>::type TR;
    std::unique_ptr<char, void(*)(void*)> own
           (
#ifndef _MSC_VER
                abi::__cxa_demangle(typeid(TR).name(), nullptr,
                                           nullptr, nullptr),
#else
                nullptr,
#endif
                std::free
           );
    std::string r = own != nullptr ? own.get() : typeid(TR).name();
    if (std::is_const<TR>::value)
        r += " const";
    if (std::is_volatile<TR>::value)
        r += " volatile";
    if (std::is_lvalue_reference<T>::value)
        r += "&";
    else if (std::is_rvalue_reference<T>::value)
        r += "&&";
    return r;
}


template<class... T>
std::vector<std::string> getTypeNames() {
    return {type_name<T>()...};
}

PyObject* sitk_2_np(sitk::Image& img){
	if (PyArray_API == NULL)
	{
		import_array();
	}
    std::vector<unsigned> size = img.GetSize();
    int ndims = size.size();
    if (ndims == 0) {
        std::cout << "empty image" << std::endl;
        return NULL;
    }
    int n_pixels = 1;
    npy_intp* image_dims = (npy_intp*)malloc(3 * sizeof(npy_intp));
    if (image_dims == NULL) {
        PyErr_NoMemory();
        return NULL;
    }
    for (int i = 0; i < ndims; i++) {
        image_dims[i] = size[ndims - i - 1];
        n_pixels *= image_dims[i];
    }
    if (ndims < 3)
        image_dims[2] = 1;
    double* img_buffer = img.GetBufferAsDouble();
    PyObject* out = PyArray_SimpleNewFromData(3, image_dims, NPY_DOUBLE, img_buffer);
    if (out == NULL) {
        std::cout << "Couldn't create numpy array" << std::endl;
    }
    free(image_dims);
    return out;
}

PyObject* sitk_2_np(const sitk::Image& img){
	if (PyArray_API == NULL)
	{
		import_array();
	}
    std::vector<unsigned> size = img.GetSize();
    int ndims = size.size();
    if (ndims == 0) {
        std::cout << "empty image" << std::endl;
        return NULL;
    }
    int n_pixels = 1;
    npy_intp* image_dims = (npy_intp*)malloc(3 * sizeof(npy_intp));
    if (image_dims == NULL) {
        PyErr_NoMemory();
        return NULL;
    }
    for (int i = 0; i < ndims; i++) {
        image_dims[i] = size[ndims - i - 1];
        n_pixels *= image_dims[i];
    }
    if (ndims < 3)
        image_dims[2] = 1;
    const double* img_buffer = img.GetBufferAsDouble();
    PyObject* out = PyArray_SimpleNew(3, image_dims, NPY_DOUBLE);
    void* arr_data = PyArray_DATA(out);
    memcpy(arr_data, img_buffer, n_pixels*sizeof(double));
    if (out == NULL) {
        std::cout << "Couldn't create numpy array" << std::endl;
    }
    free(image_dims);
    return out;
}

template<sitk::PixelIDValueEnum pixelID>
sitk::Image np_2_sitk(PyObject* arr_obj){
    std::cout << "getting array info" << std::endl;
    int ndim = PyArray_NDIM(arr_obj);
    npy_intp* dims = PyArray_DIMS(arr_obj);
    int type = PyArray_TYPE(arr_obj);
    void* data = PyArray_DATA(arr_obj);
    std::vector<unsigned> im_size;
    int n_pixels = 1;
    for(int i=0; i< ndim; i++){
        im_size.push_back(dims[ndim-1-i]);
        n_pixels *= dims[ndim-1-i];
    }
    std::cout << "constructing image" << std::endl;
    sitk::Image img(im_size, pixelID);
    std::cout << "getting buffer" << std::endl;
    CType<pixelID>::Type* buffer = PixelManagerTrait<pixelID>::GetBuffer(img);
    std::cout << "copying data" << std::endl;
    memcpy(buffer, data, n_pixels*sizeof(CType<pixelID>::Type));
    std::cout << "array converted into image" << std::endl;
    return img;
}

template<typename T>
PyObject* vector_2_np(const std::vector<T>& vec){
	if (PyArray_API == NULL)
	{
		import_array();
	}
    int ndims = 1;
    npy_intp dims = vec.size();
    PyGILState_STATE gstate;
    gstate = PyGILState_Ensure();
    PyObject* out = PyArray_SimpleNew(ndims, &dims, npy_type<T>::value);
    void* data = PyArray_DATA(out);
    memcpy(data, &vec[0], vec.size()*sizeof(T));
//    Py_XDECREF(out);
    if (out == NULL) {
        std::cout << "Couldn't create numpy array" << std::endl;
    }
    PyGILState_Release(gstate);
    return out;
}

template<typename T>
struct processed_type{
    typedef T type;
};

template<>
struct processed_type<sitk::Image&> {
    typedef PyObject* type;
};

template<>
struct processed_type<const sitk::Image&> {
    typedef PyObject* type;
};

template<>
struct processed_type<std::vector<double>&> {
    typedef PyObject* type;
};

template<typename InType>
typename processed_type<InType>::type process_arg(InType object){
    return object;
}

template<>
typename processed_type<sitk::Image&>::type process_arg<sitk::Image&>(sitk::Image& img){
    return sitk_2_np(img);
}

template<>
typename processed_type<const sitk::Image&>::type process_arg<const sitk::Image&>(const sitk::Image& img){
    return sitk_2_np(img);
}

template<>
typename processed_type<std::vector<double>&>::type process_arg<std::vector<double>&>(std::vector<double>& vec){
    return vector_2_np<double>(vec);
}


template<typename ...Args>
PyObject* BuildArgs(Args... args){
    PyObject* out_args = Py_BuildValue(args_format_string<processed_type<Args>::type...>(), process_arg<Args>(args)...);
    return out_args;
}

template<>
PyObject* BuildArgs<>(){
    PyObject* args = Py_BuildValue("()");
    return args;
}

template<typename T>
T cast_from_python(PyObject* object);

template<>
long cast_from_python<long>(PyObject* object){
    long val = PyLong_AsLong(object);
    return val;
}

template<>
int cast_from_python<int>(PyObject* object){
    int val = PyLong_AsLong(object);
    return val;
}

template<>
sitk::Image cast_from_python<sitk::Image>(PyObject* object){
    return np_2_sitk<sitk::sitkFloat64>(object);
}

template<>
void cast_from_python<void>(PyObject* object){}