#include "MinimalSurfaceEstimator.h"
#include "Utils.h"
#include <utility>
#include <omp.h>
#define LOG_RUN
#ifdef LOG_RUN
#include <fstream>
#endif
using namespace std;

MinimalSurfaceEstimator::MinimalSurfaceEstimator()
{
	omp_set_nested(1);
	omp_set_dynamic(1);
}

MinimalSurfaceEstimator::~MinimalSurfaceEstimator()
{
}


void MinimalSurfaceEstimator::Calculate(sitk::Image image, Vec3<double> point1, Vec3<double> point2, double beta, double alpha) {
	mAreaEikonal.Calculate(image, point1, point2, beta, alpha);

	mAreaEikonal.UpdateMeetingPlane();
	mAreaEikonal.CombineDistance();
	vector<double> meeting_plane_normal(mAreaEikonal.GetMeetingPlaneNormal());
	vector<double> rotation_matrix = rotation_matrix_from_vectors<double>(vector<double>({ 1, 0, 0 }), meeting_plane_normal);
	vector<double> translation = calculateOffsetFromRotation(rotation_matrix, mAreaEikonal.GetSampleImage().GetSize());
	ImageTransformCalculatedEvent(rotation_matrix, translation);
	mRotatedAreaEikonal = mAreaEikonal.Rotate(rotation_matrix);
	mAreaEikonal.SmoothDistances();
	vector<double> plane_center_physical = mAreaEikonal.GetMeetingPlaneCenter();
	vector<double> plane_center_transformed = mRotatedAreaEikonal.TransformPhysicalPointToContinuousIndex(plane_center_physical);
	double plane_slice = plane_center_transformed[0];
	PlaneCenterCalculatedEvent(plane_center_transformed);
	mRotatedAreaEikonal.Calculate();

	vector<double> start_point = mAreaEikonal.TransformContinuousIndexToPhysicalPoint(point1);
	start_point = mRotatedAreaEikonal.TransformPhysicalPointToContinuousIndex(start_point);
	vector<double> end_point = mAreaEikonal.TransformContinuousIndexToPhysicalPoint(point2);
	end_point = mRotatedAreaEikonal.TransformPhysicalPointToContinuousIndex(end_point);
	mRotatedAreaEikonal.CombineDistance(plane_slice, start_point[0], end_point[0]);
	mRotatedAreaEikonal.SmoothDistances();

	IterationEvent(DISTANCE_MEAN_ITERATION);
	//calculate slice distance
	sitk::Image distanceSlice = GetImageSlice<sitk::sitkFloat64>(mRotatedAreaEikonal.GetCombinedDistanceMap(), 0, plane_slice);
	vector<unsigned> slice_size = distanceSlice.GetSize();
	sitk::Image plane_sample_image({1, slice_size[0], slice_size[1]}, sitk::sitkFloat64);

	plane_sample_image.SetDirection(rotation_matrix);
	vector<double> origin;
	rotate(rotation_matrix, { plane_center_transformed[0], 0, 0}, origin);
	std::transform(origin.begin(), origin.end(), translation.begin(), origin.begin(), std::plus<double>());
	plane_sample_image.SetOrigin(origin);
	plane_sample_image = sitk::Resample(image, plane_sample_image, sitk::Transform(), sitk::sitkLinear, 0., sitk::sitkUnknown, true);
	sitk::Image plane_image_2d(slice_size, sitk::sitkFloat64);
	memcpy(plane_image_2d.GetBufferAsDouble(), plane_sample_image.GetBufferAsDouble(), plane_image_2d.GetNumberOfPixels() * sizeof(double));
	save_image("Y:/BIOMAG/shortest path/dist_plane_slice.tif", distanceSlice);
	save_image("Y:/BIOMAG/shortest path/img_plane_slice.tif", plane_image_2d);
	IterationEvent(PLANE_PHASEFIELD_ITERATION);
	mInitialContourCalculator.Calculate(plane_image_2d, distanceSlice);

	//mInitialContourCalculator.Calculate(plane_image_2d, distanceSlice);
	//calculate transport function
	sitk::Image initial_slice = sitk::SignedMaurerDistanceMap(sitk::Greater(mInitialContourCalculator.GetPhaseField(), 0), false, false);
	initial_slice = sitk::Clamp(initial_slice, sitk::sitkFloat64, -25, 25);
	realnum maxdist = mRotatedAreaEikonal.GetCurrentDistance();
	IterationEvent(TRANSPORT_FUNCTION_ITERATION);

	mTransportFunctionCalculator.Calculate(mRotatedAreaEikonal.GetCombinedDistanceMap(), initial_slice, plane_slice, maxdist);

	//rotate transport function to original position
	vector<unsigned int> sample_size = mAreaEikonal.GetPhiMap().GetSize();
	mTransportFunctionCalculator.RotateTransportFunction(rotation_matrix, sample_size);
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
