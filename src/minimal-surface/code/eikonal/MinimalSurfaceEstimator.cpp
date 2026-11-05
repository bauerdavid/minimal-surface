#include "MinimalSurfaceEstimator.h"
#include "Utils.h"
#include "PlaneContour.h"
#include <utility>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#define LOG_RUN
#ifdef LOG_RUN
#include <fstream>
#endif
using namespace std;

namespace {
// The calculation stops getting faster at about 6 to 8 threads (64^3 and 96^3 volumes, on an
// 8-core / 16-thread Windows machine and on 8 Linux CPUs); more threads are not slower any more,
// but not faster either. It used to be much slower: with all 16 threads the evolution took up to
// 2.8 times as long as without OpenMP, because of nested parallel regions. Unless the user chose
// a number with OMP_NUM_THREADS, a calculation uses at most kDefaultMaxThreads, so that a machine
// with many cores does not start dozens of threads for loops over a few thousand points (not
// measured on such a machine). The previous setting is put back afterwards so that other OpenMP
// code in the same process is not affected.
const int kDefaultMaxThreads = 8;

bool UserChoseThreadCount() {
#ifdef _MSC_VER
	size_t length = 0;
	getenv_s(&length, nullptr, 0, "OMP_NUM_THREADS");
	return length > 0;
#else
	return std::getenv("OMP_NUM_THREADS") != nullptr;
#endif
}

class DefaultThreadLimit {
	int mPrevious;
public:
	DefaultThreadLimit() : mPrevious(omp_get_max_threads()) {
		if (mPrevious > kDefaultMaxThreads && !UserChoseThreadCount())
			omp_set_num_threads(kDefaultMaxThreads);
	}
	~DefaultThreadLimit() {
		omp_set_num_threads(mPrevious);
	}
	DefaultThreadLimit(const DefaultThreadLimit&) = delete;
	DefaultThreadLimit& operator=(const DefaultThreadLimit&) = delete;
};
}

MinimalSurfaceEstimator::MinimalSurfaceEstimator()
{
	using namespace std::placeholders;
	omp_set_dynamic(1);
	mInitialContourCalculatorFunc = [this](sitk::Image& slice_image, sitk::Image& slice_distance){return this->CalculateInitialContourFromPhaseField(slice_image, slice_distance);};
}

MinimalSurfaceEstimator::~MinimalSurfaceEstimator()
{
}

sitk::Image MinimalSurfaceEstimator::GetTransportSliceFromPoints(sitk::Image image, Vec3<double> point1, Vec3<double> point2){
	DefaultThreadLimit thread_limit;
	_SECTION("est.GetTransportSliceFromPoints");
	mSurfaceOffsets.clear();
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
	_SECTION("est.CalculateEikonal");
	_SECTION_BEGIN(first_timer, "est.CalculateEikonal.1 first evolution");
	mAreaEikonal.Calculate(phi, point1, point2);
	_SECTION_END(first_timer);
	_SECTION_BEGIN(between_timer, "est.CalculateEikonal.2 plane, combine, rotate");
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
	_SECTION_END(between_timer);
	_SECTION_BEGIN(rotated_timer, "est.CalculateEikonal.3 rotated evolution");
#ifndef ROTATED_DISTMAP_W_FM
	mRotatedAreaEikonal.Calculate();
#endif
	_SECTION_END(rotated_timer);
	vector<double> start_point = mAreaEikonal.TransformContinuousIndexToPhysicalPoint(point1);
	start_point = mRotatedAreaEikonal.TransformPhysicalPointToContinuousIndex(start_point);
	vector<double> end_point = mAreaEikonal.TransformContinuousIndexToPhysicalPoint(point2);
	end_point = mRotatedAreaEikonal.TransformPhysicalPointToContinuousIndex(end_point);
	mRotatedPoints[0] = start_point;
	mRotatedPoints[1] = end_point;
	mRotatedAreaEikonal.SmoothDistances();
	mRotatedAreaEikonal.CombineDistance(mTransportInitPlaneSlice, start_point[0], end_point[0]);
	IterationEvent(DISTANCE_MEAN_ITERATION);
}

sitk::Image MinimalSurfaceEstimator::CalculateEikonalAndTransportInit(sitk::Image phi, sitk::Image image, Vec3<double> point1, Vec3<double> point2){
	DefaultThreadLimit thread_limit;
	CalculateEikonal(phi, point1, point2);
	_SECTION("est.initial contour");
	mSurfaceOffsets.clear();
	//calculate slice distance
	sitk::Image distanceSlice = GetImageSlice<sitk::sitkFloat64>(mRotatedAreaEikonal.GetCombinedDistanceMap(), 0, mTransportInitPlaneSlice);
	sitk::Image plane_image_2d = SampleInitPlane(image, sitk::sitkLinear);
	sitk::Image init_contour;
	mContourFromDistance = false;
	if (!mCustomInitialContourFunc && mInitialContourMethod == InitialContourMethod::Distance) {
		// At the edge of the object. The propagation settles at the bottom of the valley instead, which
		// CalculateTransportFunction takes care of.
		init_contour = CalculateInitialContourFromDistance(distanceSlice, false);
		if (init_contour.GetNumberOfPixels() == 0) {
			cout << "The distances on the meeting plane give no usable contour, using the phase field" << endl;
			init_contour = sitk::Cast(CalculateInitialContourFromPhaseField(plane_image_2d, distanceSlice), sitk::sitkUInt8);
		}
		else {
			mContourFromDistance = true;
		}
	}
	else {
		init_contour = sitk::Cast(mInitialContourCalculatorFunc(plane_image_2d, distanceSlice), sitk::sitkUInt8);
	}
	return init_contour;
}

sitk::Image MinimalSurfaceEstimator::SampleInitPlane(const sitk::Image& image, sitk::InterpolatorEnum interpolator){
	vector<unsigned> rotated_size = mRotatedAreaEikonal.GetCombinedDistanceMap().GetSize();
	vector<unsigned> slice_size = { rotated_size[1], rotated_size[2] };
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
	plane_sample_image = sitk::Resample(image, plane_sample_image, sitk::Transform(), interpolator, 0., sitk::sitkUnknown, true);
	sitk::Image plane_image_2d(slice_size, sitk::sitkFloat64);
	memcpy(plane_image_2d.GetBufferAsDouble(), plane_sample_image.GetBufferAsDouble(), plane_image_2d.GetNumberOfPixels() * sizeof(double));
	return plane_image_2d;
}

const sitk::Image& MinimalSurfaceEstimator::GetRotatedDistanceMap(int i) const {
	return mRotatedAreaEikonal.GetDistanceMap(i);
}

std::vector<double> MinimalSurfaceEstimator::GetRotatedPoint(int i) const {
	return mRotatedPoints[i];
}

int MinimalSurfaceEstimator::GetTransportInitPlaneSlice() const {
	return mTransportInitPlaneSlice;
}

sitk::Image MinimalSurfaceEstimator::CalculateEikonalAndTransportInit(sitk::Image phi, Vec3<double> point1, Vec3<double> point2){
	return CalculateEikonalAndTransportInit(phi, phi, point1, point2);
}

void MinimalSurfaceEstimator::PropagateFrom(const sitk::Image& contour, int maxIterations, bool exact){
	sitk::Image initial_slice = sitk::SignedMaurerDistanceMap(contour, true, false);
	// The distance is 0 on the voxels of the contour at its edge, so its zero level is half a voxel
	// inside of the contour. That is what it always was; only when the surface is to be at the
	// contour it is corrected.
	if (exact)
		initial_slice = sitk::Add(initial_slice, 0.5);
	initial_slice = sitk::Clamp(initial_slice, sitk::sitkFloat64, -25, 25);
	realnum maxdist = mRotatedAreaEikonal.GetCurrentDistance();
	mTransportFunctionCalculator.Calculate(mRotatedAreaEikonal.GetCombinedDistanceMap(), initial_slice, mTransportInitPlaneSlice, maxdist, maxIterations);
}

void MinimalSurfaceEstimator::CalculateTransportFunction(sitk::Image& initial_contour, int maxIterations){
	_SECTION("est.CalculateTransportFunction");
	IterationEvent(TRANSPORT_FUNCTION_ITERATION);
	mSurfaceOffsets.clear();
	PropagateFrom(initial_contour, maxIterations, mAlignSurface);
	// The propagation settles at one surface (the bottom of the valley of the distances) whatever the
	// contour in the plane was, and forgets the contour within two slices. A contour that is not on
	// that surface leaves a notch in the result at the plane. So the surface it settled at is read
	// off, extrapolated to the plane, and the propagation is done again from a contour there:
	// - for the contour that was found from the distances: that is the result;
	// - with SetAlignSurface, for any contour: the surface is then moved by the difference between
	//   the two contours, so that it passes through the given one.
	double seed_x, seed_y, chord;
	vector<double> settled;
	if ((mAlignSurface || mContourFromDistance) && PlaneSeed(seed_x, seed_y, chord) && SettledRadii(seed_x, seed_y, chord, settled)) {
		sitk::Image settled_contour = ContourFromRadii(initial_contour, seed_x, seed_y, settled, chord);
		vector<double> offsets;
		bool usable = settled_contour.GetNumberOfPixels() > 0;
		if (usable && mAlignSurface) {
			vector<double> radii = ContourRadii(initial_contour, seed_x, seed_y);
			usable = !radii.empty();
			if (usable) {
				offsets = SmoothRayValues(settled, radii);
				for (double offset : offsets)
					if (fabs(offset) > 0.25 * chord)
						usable = false;
			}
		}
		if (usable) {
			PropagateFrom(settled_contour, maxIterations, true);
			if (mAlignSurface) {
				mSurfaceOffsets = offsets;
				mPlaneSeedX = seed_x;
				mPlaneSeedY = seed_y;
				ApplySurfaceOffsets();
			}
		}
	}

	//rotate transport function to original position
	vector<unsigned int> sample_size = mAreaEikonal.GetPhiMap().GetSize();
	mTransportFunctionCalculator.RotateTransportFunction(mRotationMatrix, sample_size);
}

void MinimalSurfaceEstimator::Calculate(sitk::Image phi, sitk::Image image, Vec3<double> point1, Vec3<double> point2, int maxIterations) {
	DefaultThreadLimit thread_limit;
	_SECTION("est.Calculate (total)");
	sitk::Image initial_contour;
	if(mTempInitContour == nullptr){
		cout <<"Calculating initial contour..." <<endl;
		if(mAreaEikonal.IsUsingMeetingPoints()) {
			SetTransportInitSlice(CalculateEikonalAndTransportInit(phi, image, point1, point2));
		}
		else {
			SetTransportInitSlice(GetTransportSliceFromPoints(phi, point1, point2));
		}
	}
	if(!mAreaEikonal.IsInitialized()){
		cout <<"Calculating Eikonal" <<endl;
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


void MinimalSurfaceEstimator::HookStageInitializationEvent(eStageEnum stage, const basic_callback_type& handler)
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

void MinimalSurfaceEstimator::HookStageInitializationEvent(eStageEnum stage, basic_callback_type&& handler)
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

void MinimalSurfaceEstimator::HookStageFinishedEvent(eStageEnum stage, const basic_callback_type& handler)
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

void MinimalSurfaceEstimator::HookStageFinishedEvent(eStageEnum stage, basic_callback_type&& handler)
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

void MinimalSurfaceEstimator::HookStageIterationEvent(eStageEnum stage, const iter_callback_type& handler)
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

void MinimalSurfaceEstimator::HookStageIterationEvent(eStageEnum stage, iter_callback_type&& handler)
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

void MinimalSurfaceEstimator::HookStageUpdatedEvent(eStageEnum stage, const data_callback_type& handler, eDataEnum dataType)
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

void MinimalSurfaceEstimator::HookStageUpdatedEvent(eStageEnum stage, data_callback_type&& handler, eDataEnum dataType)
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

void MinimalSurfaceEstimator::HookStageDataInitializedEvent(eStageEnum stage, const data_callback_type& handler, eDataEnum dataType)
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

void MinimalSurfaceEstimator::HookStageDataInitializedEvent(eStageEnum stage, data_callback_type&& handler, eDataEnum dataType)
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

void MinimalSurfaceEstimator::HookCalculatedPhiMapEvent(const data_callback_type& handler)
{
	mAreaEikonal.HookCalculatedPhiMapEvent(handler);
}

void MinimalSurfaceEstimator::HookCalculatedPhiMapEvent(data_callback_type&& handler)
{
	mAreaEikonal.HookCalculatedPhiMapEvent(handler);
}

void MinimalSurfaceEstimator::HookUpdateMeetingPlaneEvent(const AreaEikonalES::plane_update_callback_type& handler)
{
	mAreaEikonal.HookUpdateMeetingPlaneEvent(handler);
	mRotatedAreaEikonal.HookUpdateMeetingPlaneEvent(handler);
}

void MinimalSurfaceEstimator::HookUpdateMeetingPlaneEvent(AreaEikonalES::plane_update_callback_type&& handler)
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
	mCustomInitialContourFunc = true;
}

void MinimalSurfaceEstimator::SetInitialContourMethod(InitialContourMethod method){
	mInitialContourMethod = method;
	mCustomInitialContourFunc = false;
	mInitialContourCalculatorFunc = [this](sitk::Image& slice_image, sitk::Image& slice_distance){return this->CalculateInitialContourFromPhaseField(slice_image, slice_distance);};
}

void MinimalSurfaceEstimator::SetInitialContourFromDistance(bool from_distance){
	SetInitialContourMethod(from_distance ? InitialContourMethod::Distance : InitialContourMethod::PhaseField);
}

bool MinimalSurfaceEstimator::PlaneSeed(double& seed_x, double& seed_y, double& chord) const {
	// The point where the line between the two points crosses the meeting plane, in the plane's
	// pixels. The plane lies between the points, so it is inside of the object, unlike the
	// midpoint projected onto a plane that is not perpendicular to the line.
	const vector<double>& p0 = mRotatedPoints[0];
	const vector<double>& p1 = mRotatedPoints[1];
	if (p0.size() != 3 || p1.size() != 3)
		return false;
	double dx = p1[0] - p0[0];
	double t = fabs(dx) > 1e-6 ? (mTransportInitPlaneSlice - p0[0]) / dx : 0.5;
	t = min(max(t, 0.0), 1.0);
	seed_x = p0[1] + t * (p1[1] - p0[1]);
	seed_y = p0[2] + t * (p1[2] - p0[2]);
	chord = sqrt((p1[0] - p0[0]) * (p1[0] - p0[0]) + (p1[1] - p0[1]) * (p1[1] - p0[1]) + (p1[2] - p0[2]) * (p1[2] - p0[2]));
	return true;
}

sitk::Image MinimalSurfaceEstimator::CalculateInitialContourFromDistance(const sitk::Image& distance_slice, bool at_floor) {
	double seed_x, seed_y, chord;
	if (!PlaneSeed(seed_x, seed_y, chord))
		return sitk::Image();
	return ContourFromDistanceSlice(distance_slice, seed_x, seed_y, chord, at_floor);
}

void MinimalSurfaceEstimator::SetAlignSurface(bool align) {
	mAlignSurface = align;
}

bool MinimalSurfaceEstimator::SettledRadii(double seed_x, double seed_y, double chord, std::vector<double>& radii) const {
	const TransportFunctionES& calculator = mTransportFunctionCalculator;
	const sitk::Image& transport = calculator.GetTransportFunction();
	vector<unsigned> size = transport.GetSize();
	const int xs = size[0], ys = size[1], zs = size[2];
	const double* buffer = transport.GetBufferAsDouble();
	const int rays = PlaneContourRayCount();
	// Planes this far from the contour's: out of its reach (it is forgotten within two slices), and
	// not so far that the object has changed much. The radius of an object changes with the square
	// of the distance from the middle, so a + b * t^2 is fitted and a is the radius at the plane.
	const int distances[] = { 3, 4, 5 };
	vector<double> sum_r(rays, 0.0), sum_rt(rays, 0.0);
	double sum_t = 0, sum_tt = 0;
	int count = 0;
	for (int sign = -1; sign <= 1; sign += 2) {
		for (int distance : distances) {
			const int x = mTransportInitPlaneSlice + sign * distance;
			if (x < 1 || x > xs - 2)
				return false;
			vector<double> values((size_t)ys * zs);
			for (int zz = 0; zz < zs; zz++)
				for (int yy = 0; yy < ys; yy++)
					values[(size_t)zz * ys + yy] = buffer[x + (size_t)xs * (yy + (size_t)ys * zz)];
			vector<double> r = FieldRadii(values, ys, zs, seed_x, seed_y, 0.0);
			if (r.empty())
				return false;
			const double tt = (double)distance * distance;
			for (int i = 0; i < rays; i++) {
				if (r[i] < 0)
					return false;
				sum_r[i] += r[i];
				sum_rt[i] += r[i] * tt;
			}
			sum_t += tt;
			sum_tt += tt * tt;
			count++;
		}
	}
	const double denominator = count * sum_tt - sum_t * sum_t;
	radii.assign(rays, 0.0);
	for (int i = 0; i < rays; i++) {
		double b = (count * sum_rt[i] - sum_t * sum_r[i]) / denominator;
		radii[i] = max((sum_r[i] - b * sum_t) / count, 1.0);
		if (fabs(radii[i] - sum_r[i] / count) > 0.2 * chord)
			return false;
	}
	return true;
}

void MinimalSurfaceEstimator::ApplySurfaceOffsets() {
	// The transport function is at the bottom of the valley away from the plane, and the contour of
	// the plane there. Make it the signed distance to its zero level, and shift that by the offset
	// of the ray it is on (the half is the distance to the edge of a voxel, see CalculateTransportFunction).
	// The shift is the same all the way along the object, so that the surface
	// is as smooth as it was, passes through the given contour in the plane, and has no notch.
	const TransportFunctionES& calculator = mTransportFunctionCalculator;
	sitk::Image transport = calculator.GetTransportFunction();
	sitk::Image inside = sitk::GreaterEqual(transport, 0.0);
	sitk::Image distance = sitk::Cast(sitk::SignedMaurerDistanceMap(inside, true, false, false, 0.0), sitk::sitkFloat64);
	vector<unsigned> size = distance.GetSize();
	const int xs = size[0], ys = size[1], zs = size[2];
	double* buffer = distance.GetBufferAsDouble();
	for (int zz = 0; zz < zs; zz++) {
		for (int yy = 0; yy < ys; yy++) {
			const double shift = PlaneContourRayValue(mSurfaceOffsets, yy - mPlaneSeedX, zz - mPlaneSeedY);
			double* row = buffer + (size_t)xs * (yy + (size_t)ys * zz);
			for (int xx = 0; xx < xs; xx++)
				row[xx] = min(max(row[xx] + 0.5 - shift, -25.0), 25.0);
		}
	}
	mTransportFunctionCalculator.SetTransportFunction(distance);
}

const sitk::Image& MinimalSurfaceEstimator::GetCombinedDistanceMap() const{
	return mRotatedAreaEikonal.GetCombinedDistanceMap();
}

sitk::Image MinimalSurfaceEstimator::GetTempInitContour() const {
	return *mTempInitContour;
}

