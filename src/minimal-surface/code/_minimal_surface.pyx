# distutils: language = c++
from libcpp.functional cimport function
from libcpp.vector cimport vector
from libc.string cimport memcpy
from libcpp cimport bool
cimport numpy as np
import numpy as np


DISTANCE_ITERATION = 1
DISTANCE_MEAN_ITERATION = 2
PLANE_PHASEFIELD_ITERATION = 3
TRANSPORT_FUNCTION_ITERATION = 4
DONE_ITERATION = 5

AREA_EIKONAL_STAGE = 0
ROTATED_AREA_EIKONAL_STAGE = 1
PLANE_PHASEFIELD_STAGE = 2
TRANSPORT_FUNCTION_STAGE = 3

DEFAULT_DATA = 0
DISTANCE_DATA = 0
PHASEFIELD_DATA = 1
MEETING_POINTS_DATA = 2
CURVATURE_DATA = 3

cdef extern from *:
    ctypedef int Image_ref "itk::simple::Image&" #hack
    ctypedef int init_image_func_header "itk::simple::Image(itk::simple::Image&, itk::simple::Image&)" #hack No. 2

cdef extern from "Vec.h":
    cdef cppclass Vec3[T]:
        Vec3()
        Vec3(T x, T y, T z)
        T& x()
        T& y()
        T& z()
        T* begin()

cdef extern from "SimpleITK.h" namespace "itk::simple":
    cdef enum PixelIDValueEnum:
        sitkFloat64 = 9
    ctypedef PixelIDValueEnum sitkFloat64_ "sitk::sitkFloat64"
    cdef cppclass Image:
        Image()
        Image(int, int, PixelIDValueEnum)
        Image(vector[unsigned int], PixelIDValueEnum)
        double* GetBufferAsDouble()
    cdef Image Cast(Image, PixelIDValueEnum)

cdef extern from "python_utils.h":
    cdef Image np_2_sitk[PID](object)
    cdef object sitk_2_np(Image)

cdef extern from "pythonwrapper/PyCallableWrapper.hpp":
    cdef cppclass CallbackWrapper_0A:
        CallbackWrapper_0A(object)
        CallbackWrapper_0A()
        void operator()()

    cdef cppclass CallbackWrapper_1A[T]:
        CallbackWrapper_1A(object)
        CallbackWrapper_1A()
        void operator()(T)

    cdef cppclass CallbackWrapper_2A[T1, T2]:
        CallbackWrapper_2A(object)
        CallbackWrapper_2A()
        void operator()(T1, T2)

    cdef cppclass InitialContourCalculatorWrapper:
        InitialContourCalculatorWrapper(object)
        InitialContourCalculatorWrapper()

ctypedef function[void()] basic_callback_type
ctypedef function[void(int)] iter_callback_type
ctypedef function[void(Image_ref, int)] data_callback_type
ctypedef function[void(vector[double]&, vector[double]&)] image_transform_callback_type
ctypedef function[void(vector[double]&)] vector_callback_type
ctypedef function[init_image_func_header] init_contour_callback_type

ctypedef CallbackWrapper_0A basic_callback_wrapper
ctypedef CallbackWrapper_1A[int] iter_callback_wrapper
ctypedef CallbackWrapper_2A[Image_ref, int] data_callback_wrapper
ctypedef CallbackWrapper_2A[vector[double]&, vector[double]&] image_transform_callback_wrapper
ctypedef CallbackWrapper_1A[vector[double]&] vector_callback_wrapper

cdef extern from "MinimalSurfaceEstimator.h":
    cdef cppclass TransportFunctionES:
        const Image& GetTransportFunction() nogil const
    cdef enum eStageEnum:
        pass

    cdef enum eDataEnum:
        pass

    cdef cppclass MinimalSurfaceEstimator:
        void HookStageInitializationEvent(eStageEnum, basic_callback_type&&)
        void HookStageFinishedEvent(eStageEnum, basic_callback_type&&)
        void HookStageIterationEvent(eStageEnum, iter_callback_type&&)
        void HookStageDataInitializedEvent(eStageEnum, data_callback_type&&, eDataEnum)
        void HookStageUpdatedEvent(eStageEnum, data_callback_type&&, eDataEnum)
        void HookImageTransformCalculatedEvent(image_transform_callback_type&&)
        void HookPlaneCenterCalculatedEvent(vector_callback_type&&)
        void HookIterationEvent(iter_callback_type&&)
        void SetUsesCorrection(bool)
        void SetUsingMeetingPoints(bool)
        void Calculate(Image, Image, Vec3[double], Vec3[double], int) nogil
        void Calculate(Image, Vec3[double], Vec3[double], int) nogil
        void SetInitialContourCalculatorFunc(init_contour_callback_type)
        const TransportFunctionES& GetTransportFunctionCalculator() nogil const
        Image GetTransportSliceFromPoints(Image, Vec3[double], Vec3[double]) nogil
        Image CalculateEikonalAndTransportInit(Image, Image, Vec3[double], Vec3[double]) nogil
        void SetTransportInitSlice(Image) nogil
        Image GetCombinedDistanceMap() nogil
        Image GetTempInitContour() nogil
    cdef vector[Vec3[int]] ResolvePath(Vec3[int], const Image&) nogil

cdef class MinimalSurfaceCalculator:
    cdef MinimalSurfaceEstimator calculator
    def __cinit__(self):
        ...

    cpdef void hook_iteration_event(self, object callback):
        if not callable(callback):
            print("callback is not callable!")
            return
        cdef iter_callback_wrapper wrapper = iter_callback_wrapper(callback)
        self.calculator.HookIterationEvent(<iter_callback_type>wrapper)

    cpdef void hook_stage_init_event(self, int stage, object callback):
        if not callable(callback):
            print("callback is not callable!")
            return
        cdef basic_callback_wrapper wrapper = basic_callback_wrapper(callback)
        cdef basic_callback_type cast_wrapper = <basic_callback_type>wrapper
        self.calculator.HookStageInitializationEvent(<eStageEnum>stage, cast_wrapper)

    cpdef void hook_stage_finished_event(self, int stage, object callback):
        if not callable(callback):
            print("callback is not callable!")
            return
        cdef basic_callback_wrapper wrapper = basic_callback_wrapper(callback)
        self.calculator.HookStageFinishedEvent(<eStageEnum>stage, <basic_callback_type>wrapper)

    cpdef void hook_stage_iteration_event(self, int stage, object callback):
        if not callable(callback):
            print("callback is not callable!")
            return
        cdef iter_callback_wrapper wrapper = iter_callback_wrapper(callback)
        self.calculator.HookStageIterationEvent(<eStageEnum> stage, <iter_callback_type> wrapper)

    cpdef void hook_stage_updated_event(self, int stage, object callback):
        if not callable(callback):
            print("callback is not callable!")
            return
        cdef iter_callback_wrapper wrapper = iter_callback_wrapper(callback)
        self.calculator.HookStageIterationEvent(<eStageEnum> stage, <iter_callback_type> wrapper)

    cpdef void hook_stage_data_init_event(self, int stage, object callback, int data_type=0):
        if not callable(callback):
            print("callback is not callable!")
            return
        cdef data_callback_wrapper wrapper = data_callback_wrapper(callback)
        self.calculator.HookStageDataInitializedEvent(<eStageEnum> stage, <data_callback_type> wrapper, <eDataEnum> data_type)

    cpdef void hook_transform_calculated_event(self, object callback):
        if not callable(callback):
            print("callback is not callable!")
            return
        cdef image_transform_callback_wrapper wrapper = image_transform_callback_wrapper(callback)
        self.calculator.HookImageTransformCalculatedEvent(<image_transform_callback_type> wrapper)

    cpdef void hook_plane_center_calculated_event(self, object callback):
        if not callable(callback):
            print("callback is not callable!")
            return
        cdef vector_callback_wrapper wrapper = vector_callback_wrapper(callback)
        self.calculator.HookPlaneCenterCalculatedEvent(<vector_callback_type> wrapper)

    cpdef void set_initial_plane_calculator(self, object func):
        if not callable(func):
            print("func is not callable!")
            return
        cdef InitialContourCalculatorWrapper wrapper = InitialContourCalculatorWrapper(func)
        self.calculator.SetInitialContourCalculatorFunc(<init_contour_callback_type> wrapper)

    cpdef np.ndarray[np.int_t, ndim=2] resolve_shortest_paths(self, np.ndarray[np.int_t, ndim=1] point, np.ndarray[np.float_t, ndim=3] data):
        if len(point) != 3:
            print("point should be a size 3 array")
        cdef Vec3[int] point_vec = Vec3[int](point[0], point[1], point[2])
        cdef vector[unsigned int] im_size = [data.shape[0], data.shape[1], data.shape[2]]
        cdef Image img = Image(im_size, sitkFloat64)
        cdef int[::1] temp_view = data
        memcpy(img.GetBufferAsDouble(), &temp_view[0], im_size[2] * im_size[1] * im_size[0] * sizeof(double))
        cdef vector[Vec3[int]] path
        with nogil:
            path = ResolvePath(point_vec, img)
        cdef np.ndarray[np.int_t, ndim=2] path_arr = np.ndarray((path.size(), 3), dtype=int)
        cdef int i
        for i in range(path.size()):
            path_arr[i, 0] = path[i].x()
            path_arr[i, 1] = path[i].y()
            path_arr[i, 2] = path[i].z()
        return path_arr

    cpdef np.ndarray[np.float_t, ndim=3] calculate(
            self,
            np.ndarray[np.float_t, ndim=3] phi,
            np.ndarray[np.float_t, ndim=3] image,
            np.ndarray[np.float_t, ndim=1] point1,
            np.ndarray[np.float_t, ndim=1] point2,
            bool use_correction,
            int max_iterations=10000,
    ):
        cdef int i
        cdef Image sitk_image = np_2_sitk[sitkFloat64_](image)
        cdef Image sitk_phi = np_2_sitk[sitkFloat64_](phi)
        cdef Vec3[double] point1_vec
        cdef Vec3[double] point2_vec
        cdef double* point1_data = point1_vec.begin()
        cdef double* point2_data = point2_vec.begin()
        for i in range(3):
            point1_data[i] = point1[i]
            point2_data[i] = point2[i]
        self.calculator.SetUsesCorrection(use_correction)
        with nogil:
            self.calculator.Calculate(sitk_phi, sitk_image, point1_vec, point2_vec, max_iterations)
        cdef const double* transport_buffer = self.calculator.GetTransportFunctionCalculator().GetTransportFunction().GetBufferAsDouble()
        cdef np.ndarray[np.float_t, ndim=3] output = np.empty((image.shape[0], image.shape[1], image.shape[2]), dtype=float)
        cdef double[:, :, :] out_view = output
        memcpy(&out_view[0, 0, 0], transport_buffer, image.size * sizeof(double))
        return output

    cpdef void init_transport_slice(self,
            np.ndarray[np.float_t, ndim=3] image,
            np.ndarray[np.float_t, ndim=1] point1,
            np.ndarray[np.float_t, ndim=1] point2):
        cdef int i
        cdef vector[unsigned int] im_size
        for i in range(3):
            im_size.push_back(<unsigned int>image.shape[2-i])
        cdef Image sitk_image = Image(im_size, sitkFloat64)
        cdef double[:, :, :] image_data = image
        memcpy(sitk_image.GetBufferAsDouble(), &image_data[0, 0, 0], image.size*sizeof(double))
        cdef Vec3[double] point1_vec
        cdef Vec3[double] point2_vec
        cdef double* point1_data = point1_vec.begin()
        cdef double* point2_data = point2_vec.begin()
        for i in range(3):
            point1_data[i] = point1[i]
            point2_data[i] = point2[i]
        cdef Image transport_slice
        with nogil:
            transport_slice = self.calculator.GetTransportSliceFromPoints(sitk_image, point1_vec, point2_vec)
            self.calculator.SetTransportInitSlice(transport_slice)

    cpdef void calc_eikonal_and_transport_init(
            self,
            np.ndarray[np.float_t, ndim=3] phi,
            np.ndarray[np.float_t, ndim=3] image,
            np.ndarray[np.float_t, ndim=1] point1,
            np.ndarray[np.float_t, ndim=1] point2,
            bool use_correction,
    ):
        cdef int i
        cdef Image sitk_phi = np_2_sitk[sitkFloat64_](phi)
        cdef Image sitk_image = np_2_sitk[sitkFloat64_](image) if image is not None else sitk_phi
        cdef Vec3[double] point1_vec
        cdef Vec3[double] point2_vec
        cdef double* point1_data = point1_vec.begin()
        cdef double* point2_data = point2_vec.begin()
        for i in range(3):
            point1_data[i] = point1[i]
            point2_data[i] = point2[i]
        cdef Image transport_slice
        self.calculator.SetUsesCorrection(use_correction)
        with nogil:
            transport_slice = self.calculator.CalculateEikonalAndTransportInit(sitk_phi, sitk_image, point1_vec, point2_vec)
            self.calculator.SetTransportInitSlice(transport_slice)

    cpdef void set_using_meeting_points(self, bool use_meeting_points):
        self.calculator.SetUsingMeetingPoints(use_meeting_points)

    cpdef object get_distance_map(self):
        cdef Image distance_map = self.calculator.GetCombinedDistanceMap()
        cdef object arr = sitk_2_np(distance_map).copy()
        return arr

    cpdef object get_init_plane(self):
        cdef Image init_plane = Cast(self.calculator.GetTempInitContour(), sitkFloat64)
        print("init_plane")
        cdef object arr = sitk_2_np(init_plane).copy()
        print("arr")
        return arr
