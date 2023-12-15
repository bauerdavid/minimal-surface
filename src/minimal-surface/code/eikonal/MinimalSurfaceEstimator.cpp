#include "MinimalSurfaceEstimator.h"
#include "Utils.h"
#include <utility>
#include <cmath>
#include <omp.h>
#define LOG_RUN
#ifdef LOG_RUN
#include <fstream>
#endif
using namespace std;

MinimalSurfaceEstimator::MinimalSurfaceEstimator()
{
    using namespace std::placeholders;
	omp_set_nested(1);
	omp_set_dynamic(1);
	mInitialContourCalculatorFunc = [this](sitk::Image& slice_image, sitk::Image& slice_distance){return this->CalculateInitialContourFromPhaseField(slice_image, slice_distance);};
}

MinimalSurfaceEstimator::~MinimalSurfaceEstimator()
{
}

sitk::Image MinimalSurfaceEstimator::GetTransportSliceFromPoints(sitk::Image image, Vec3<double> point1, Vec3<double> point2){
    Vec3<double> center = (point1 + point2)/2;
    mAreaEikonal.SetMeetingPlaneCenter(center);
    Vec3<double> normal_vec = (point2-point1);
    normal_vec /= normal_vec.Norm();
    normal_vec = StandardizeVector(normal_vec);
    mAreaEikonal.SetMeetingPlaneNormal(normal_vec);
    mRotationMatrix = mAreaEikonal.GetRotationMatrix();
    vector<double> translation = calculateOffsetFromRotation(mRotationMatrix, image.GetSize());

    sitk::Image rotated_image, sample_image;
    sample_image = resample_img(image, rotated_image, mRotationMatrix);
    vector<unsigned> rotated_size = rotated_image.GetSize();

    vector<double> plane_center_transformed = sample_image.TransformPhysicalPointToContinuousIndex(center);
	mTransportInitPlaneSlice = (int)round(plane_center_transformed[0]);
	sitk::Image plane_sample_image({1, rotated_size[1], rotated_size[2]}, sitk::sitkFloat64);

	vector<double> origin;
	rotate(mRotationMatrix, { plane_center_transformed[0], 0, 0}, origin);
    std::transform(origin.begin(), origin.end(), translation.begin(), origin.begin(), std::plus<double>());
    plane_sample_image.SetDirection(mRotationMatrix);
    plane_sample_image.SetOrigin(origin);
	plane_sample_image = sitk::Resample(image, plane_sample_image, sitk::Transform(), sitk::sitkLinear, 0., sitk::sitkUnknown, true);
	sitk::Image plane_image_2d({rotated_size[1], rotated_size[2]}, sitk::sitkFloat64);
	memcpy(plane_image_2d.GetBufferAsDouble(), plane_sample_image.GetBufferAsDouble(), plane_image_2d.GetNumberOfPixels() * sizeof(double));
	return sitk::Cast(mInitialContourCalculatorFunc(plane_image_2d, plane_image_2d), sitk::sitkUInt8);
}

void MinimalSurfaceEstimator::CalculateEikonal(sitk::Image phi, Vec3<double> point1, Vec3<double> point2){
    mAreaEikonal.Calculate(phi, point1, point2);
    if(mAreaEikonal.IsUsingMeetingPoints())
	    mAreaEikonal.UpdateMeetingPlane();
	mAreaEikonal.CombineDistance();
	mRotationMatrix = mAreaEikonal.GetRotationMatrix();
	vector<double> translation = calculateOffsetFromRotation(mRotationMatrix, mAreaEikonal.GetSampleImage().GetSize());
	ImageTransformCalculatedEvent(mRotationMatrix, translation);
	mRotatedAreaEikonal = mAreaEikonal.Rotate(mRotationMatrix);
	mAreaEikonal.SmoothDistances();
	vector<double> plane_center_physical = mAreaEikonal.GetMeetingPlaneCenter();
	vector<double> plane_center_transformed = mRotatedAreaEikonal.TransformPhysicalPointToContinuousIndex(plane_center_physical);
	mTransportInitPlaneSlice = (int)round(plane_center_transformed[0]);
	PlaneCenterCalculatedEvent(plane_center_transformed);
#ifndef ROTATED_DISTMAP_W_FM
	mRotatedAreaEikonal.Calculate();
#endif
	vector<double> start_point = mAreaEikonal.TransformContinuousIndexToPhysicalPoint(point1);
	start_point = mRotatedAreaEikonal.TransformPhysicalPointToContinuousIndex(start_point);
	vector<double> end_point = mAreaEikonal.TransformContinuousIndexToPhysicalPoint(point2);
	end_point = mRotatedAreaEikonal.TransformPhysicalPointToContinuousIndex(end_point);
	mRotatedAreaEikonal.SmoothDistances();
	mRotatedAreaEikonal.CombineDistance(mTransportInitPlaneSlice, start_point[0], end_point[0]);
	IterationEvent(DISTANCE_MEAN_ITERATION);
}

sitk::Image MinimalSurfaceEstimator::CalculateEikonalAndTransportInit(sitk::Image phi, sitk::Image image, Vec3<double> point1, Vec3<double> point2){
    CalculateEikonal(phi, point1, point2);
	//calculate slice distance
	sitk::Image distanceSlice = GetImageSlice<sitk::sitkFloat64>(mRotatedAreaEikonal.GetCombinedDistanceMap(), 0, mTransportInitPlaneSlice);
	vector<unsigned> slice_size = distanceSlice.GetSize();
	sitk::Image plane_sample_image({1, slice_size[0], slice_size[1]}, sitk::sitkFloat64);
    mRotationMatrix = mAreaEikonal.GetRotationMatrix();
    vector<double> translation = calculateOffsetFromRotation(mRotationMatrix, mAreaEikonal.GetSampleImage().GetSize());
    vector<double> plane_center_physical = mAreaEikonal.GetMeetingPlaneCenter();
	vector<double> plane_center_transformed = mRotatedAreaEikonal.TransformPhysicalPointToContinuousIndex(plane_center_physical);
	plane_sample_image.SetDirection(mRotationMatrix);
	vector<double> origin;
	rotate(mRotationMatrix, { plane_center_transformed[0], 0, 0}, origin);
	std::transform(origin.begin(), origin.end(), translation.begin(), origin.begin(), std::plus<double>());
	plane_sample_image.SetOrigin(origin);
	plane_sample_image = sitk::Resample(image, plane_sample_image, sitk::Transform(), sitk::sitkLinear, 0., sitk::sitkUnknown, true);
	sitk::Image plane_image_2d(slice_size, sitk::sitkFloat64);
	memcpy(plane_image_2d.GetBufferAsDouble(), plane_sample_image.GetBufferAsDouble(), plane_image_2d.GetNumberOfPixels() * sizeof(double));
	sitk::Image init_contour = sitk::Cast(mInitialContourCalculatorFunc(plane_image_2d, distanceSlice), sitk::sitkUInt8);
    return init_contour;
}

sitk::Image MinimalSurfaceEstimator::CalculateEikonalAndTransportInit(sitk::Image phi, Vec3<double> point1, Vec3<double> point2){
    return CalculateEikonalAndTransportInit(phi, phi, point1, point2);
}

void MinimalSurfaceEstimator::CalculateTransportFunction(sitk::Image& initial_contour, int maxIterations){
    sitk::Image initial_slice = sitk::SignedMaurerDistanceMap(initial_contour, true, false);
    initial_slice = sitk::Clamp(initial_slice, sitk::sitkFloat64, -25, 25);
    realnum maxdist = mRotatedAreaEikonal.GetCurrentDistance();
    IterationEvent(TRANSPORT_FUNCTION_ITERATION);
    auto size = mRotatedAreaEikonal.GetCombinedDistanceMap().GetSize();
    mTransportFunctionCalculator.Calculate(mRotatedAreaEikonal.GetCombinedDistanceMap(), initial_slice, mTransportInitPlaneSlice, maxdist, maxIterations);

	//rotate transport function to original position
	vector<unsigned int> sample_size = mAreaEikonal.GetPhiMap().GetSize();
	mTransportFunctionCalculator.RotateTransportFunction(mRotationMatrix, sample_size);
}

void MinimalSurfaceEstimator::Calculate(sitk::Image phi, sitk::Image image, Vec3<double> point1, Vec3<double> point2, int maxIterations) {
	sitk::Image initial_contour;
	if(mTempInitContour == nullptr){
        if(mAreaEikonal.IsUsingMeetingPoints()) {
            SetTransportInitSlice(CalculateEikonalAndTransportInit(phi, image, point1, point2));
        }
        else {
            SetTransportInitSlice(GetTransportSliceFromPoints(phi, point1, point2));
        }
	}
	if(!mAreaEikonal.IsInitialized()){
        CalculateEikonal(phi, point1, point2);
    }
	initial_contour = *mTempInitContour;
    CalculateTransportFunction(initial_contour, maxIterations);
}

void MinimalSurfaceEstimator::Calculate(sitk::Image phi, Vec3<double> point1, Vec3<double> point2, int maxIterations){
    Calculate(phi, phi, point1, point2, maxIterations);
}

void MinimalSurfaceEstimator::SetTransportInitSlice(const sitk::Image& image){
    mTempInitContour.reset(new sitk::Image(image));
    IterationEvent(PLANE_PHASEFIELD_ITERATION);
}


void MinimalSurfaceEstimator::HookStageInitializationEvent(eStageEnum stage, basic_callback_type& handler)
{
	switch (stage) {
	case AreaEikonalStage:
		mAreaEikonal.HookInitializationEvent(handler);
		break;
	case RotatedAreaEikonalStage:
		mRotatedAreaEikonal.HookInitializationEvent(handler);
		break;
	case PlanePhaseFieldStage:
		mInitialContourCalculator.HookInitializationEvent(handler);
		break;
	case TransportFunctionStage:
		mTransportFunctionCalculator.HookInitializationEvent(handler);
		break;
	}
}

void MinimalSurfaceEstimator::HookStageFinishedEvent(eStageEnum stage, basic_callback_type& handler)
{
	switch (stage) {
	case AreaEikonalStage:
		mAreaEikonal.HookFinishedEvent(handler);
		break;
	case RotatedAreaEikonalStage:
		mRotatedAreaEikonal.HookFinishedEvent(handler);
		break;
	case PlanePhaseFieldStage:
		mInitialContourCalculator.HookFinishedEvent(handler);
		break;
	case TransportFunctionStage:
		mTransportFunctionCalculator.HookFinishedEvent(handler);
		break;
	}
}

void MinimalSurfaceEstimator::HookStageIterationEvent(eStageEnum stage, iter_callback_type& handler)
{
	switch (stage) {
	case AreaEikonalStage:
		mAreaEikonal.HookIterationEvent(handler);
		break;
	case RotatedAreaEikonalStage:
		mRotatedAreaEikonal.HookIterationEvent(handler);
		break;
	case PlanePhaseFieldStage:
		mInitialContourCalculator.HookIterationEvent(handler);
		break;
	case TransportFunctionStage:
		mTransportFunctionCalculator.HookIterationEvent(handler);
		break;
	}
}

void MinimalSurfaceEstimator::HookStageUpdatedEvent(eStageEnum stage, data_callback_type& handler, eDataEnum dataType)
{
	switch (stage) {
	case AreaEikonalStage:
		switch (dataType) {
		case Distance:
			mAreaEikonal.HookUpdatedDistanceMapEvent(handler);
			break;
		case PhaseField:
			mAreaEikonal.HookUpdatedPhaseFieldMapEvent(handler);
			break;
		case Curvature:
			mAreaEikonal.HookUpdatedCurvatureMapEvent(handler);
			break;
		case MeetingPoints:
			mAreaEikonal.HookUpdatedMeetingPointsMapEvent(handler);
			break;
		}
		break;
	case RotatedAreaEikonalStage:
		switch (dataType) {
		case Distance:
			mRotatedAreaEikonal.HookUpdatedDistanceMapEvent(handler);
			break;
		case PhaseField:
			mRotatedAreaEikonal.HookUpdatedPhaseFieldMapEvent(handler);
			break;
		case Curvature:
			mRotatedAreaEikonal.HookUpdatedCurvatureMapEvent(handler);
			break;
		case MeetingPoints:
			mRotatedAreaEikonal.HookUpdatedMeetingPointsMapEvent(handler);
			break;
		}
		break;
	case PlanePhaseFieldStage:
		mInitialContourCalculator.HookUpdatedPhaseFieldMapEvent(handler);
		break;
	case TransportFunctionStage:
		mTransportFunctionCalculator.HookUpdatedTransportFunctionMapEvent(handler);
		break;
	}
}

void MinimalSurfaceEstimator::HookStageDataInitializedEvent(eStageEnum stage, data_callback_type& handler, eDataEnum dataType)
{
	switch (stage) {
	case AreaEikonalStage:
		switch (dataType) {
		case Distance:
			mAreaEikonal.HookInitializedDistanceMapEvent(handler);
			break;
		case PhaseField:
			mAreaEikonal.HookInitializedPhaseFieldMapEvent(handler);
			break;
		case Curvature:
			mAreaEikonal.HookInitializedCurvatureMapEvent(handler);
			break;
		case MeetingPoints:
			mAreaEikonal.HookInitializedMeetingPointsMapEvent(handler);
			break;
		}
		break;
	case RotatedAreaEikonalStage:
		switch (dataType) {
		case Distance:
			mRotatedAreaEikonal.HookInitializedDistanceMapEvent(handler);
			break;
		case PhaseField:
			mRotatedAreaEikonal.HookInitializedPhaseFieldMapEvent(handler);
			break;
		case Curvature:
			mRotatedAreaEikonal.HookInitializedCurvatureMapEvent(handler);
			break;
		case MeetingPoints:
			mRotatedAreaEikonal.HookInitializedMeetingPointsMapEvent(handler);
			break;
		}
		break;
	case PlanePhaseFieldStage:
		mInitialContourCalculator.HookInitializedPhaseFieldMapEvent(handler);
		break;
	case TransportFunctionStage:
		mTransportFunctionCalculator.HookInitializedTransportFunctionMapEvent(handler);
		break;
	}
}

void MinimalSurfaceEstimator::HookCalculatedPhiMapEvent(data_callback_type& handler)
{
	mAreaEikonal.HookCalculatedPhiMapEvent(handler);
}

void MinimalSurfaceEstimator::HookUpdateMeetingPlaneEvent(AreaEikonalES::plane_update_callback_type& handler)
{
	mAreaEikonal.HookUpdateMeetingPlaneEvent(handler);
	mRotatedAreaEikonal.HookUpdateMeetingPlaneEvent(handler);
}

void MinimalSurfaceEstimator::SetUsesCorrection(bool useCorrection)
{
	mAreaEikonal.SetUsesCorrection(useCorrection);
	mRotatedAreaEikonal.SetUsesCorrection(useCorrection);
}

void MinimalSurfaceEstimator::SetUsingMeetingPoints(bool useMeetingPoints){
    mAreaEikonal.SetUsingMeetingPoints(useMeetingPoints);
}


const AreaEikonalES& MinimalSurfaceEstimator::GetAreaEikonal() const
{
	return mAreaEikonal;
}

const AreaEikonalES& MinimalSurfaceEstimator::GetRotatedAreaEikonal() const
{
	return mRotatedAreaEikonal;
}

const PlanePhaseFieldES& MinimalSurfaceEstimator::GetInitialContourCalculator() const
{
	return mInitialContourCalculator;
}

const TransportFunctionES& MinimalSurfaceEstimator::GetTransportFunctionCalculator() const
{
	return mTransportFunctionCalculator;
}

sitk::Image MinimalSurfaceEstimator::CalculateInitialContourFromPhaseField(sitk::Image& plane_image_2d, sitk::Image& distanceSlice) {
    mInitialContourCalculator.Calculate(plane_image_2d, distanceSlice);
    return sitk::Less(mInitialContourCalculator.GetPhaseField(), 0);
}

void MinimalSurfaceEstimator::SetInitialContourCalculatorFunc(std::function<sitk::Image(sitk::Image&, sitk::Image&)> func){
    mInitialContourCalculatorFunc = func;
}

const sitk::Image& MinimalSurfaceEstimator::GetCombinedDistanceMap() const{
    return mRotatedAreaEikonal.GetCombinedDistanceMap();
}

sitk::Image MinimalSurfaceEstimator::GetTempInitContour() const {
    return *mTempInitContour;
}

