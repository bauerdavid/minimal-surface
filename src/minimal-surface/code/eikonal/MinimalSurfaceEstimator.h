#pragma once
#include "PlanePhaseField.h"
#include "AreaEikonal.h"
#include "Transport.h"
#include "Vec.h"
#include "SimpleITK.h"
#include "szevent.h"
#include <functional>
#include <memory>



#define DISTANCE_ITERATION 1
#define DISTANCE_MEAN_ITERATION 2
#define PLANE_PHASEFIELD_ITERATION 3
#define TRANSPORT_FUNCTION_ITERATION 4
#define DONE_ITERATION 5

#define AREA_EIKONAL_STAGE 0
#define ROTATED_AREA_EIKONAL_STAGE 1
#define PLANE_PHASEFIELD_STAGE 2
#define TRANSPORT_FUNCTION_STAGE 3

#define DEFAULT_DATA 0
#define DISTANCE_DATA 0
#define PHASEFIELD_DATA 1
#define MEETING_POINTS_DATA 2
#define CURVATURE_DATA 3

enum eStageEnum { AreaEikonalStage, RotatedAreaEikonalStage, PlanePhaseFieldStage, TransportFunctionStage
};
enum eDataEnum {Default, Distance=0, PhaseField, MeetingPoints, Curvature };
// How the contour of the meeting plane is found when the user has not set a function for it.
// Distance: from the distances of the two points on the plane (needs the meeting plane to come
// from the eikonal; falls back to PhaseField when they give nothing usable). PhaseField: evolves a
// phase field in the plane, which is much slower.
enum class InitialContourMethod { Distance, PhaseField };

namespace sitk = itk::simple;
class MinimalSurfaceEstimator: public EventSource
{
	AreaEikonalES mAreaEikonal;
	AreaEikonalES mRotatedAreaEikonal;
	PlanePhaseFieldES mInitialContourCalculator;
	TransportFunctionES mTransportFunctionCalculator;
	std::function<sitk::Image(sitk::Image&, sitk::Image&)> mInitialContourCalculatorFunc;
	std::unique_ptr<sitk::Image> mTempInitContour;
	InitialContourMethod mInitialContourMethod = InitialContourMethod::Distance;
	bool mCustomInitialContourFunc = false;
	// See SetAlignSurface. The offsets are empty when there is nothing to align.
	bool mAlignSurface = false;
	// Whether the contour of the plane was found from the distances (and so is not the user's).
	bool mContourFromDistance = false;
	std::vector<double> mSurfaceOffsets;
	double mPlaneSeedX = 0, mPlaneSeedY = 0;
	int mTransportInitPlaneSlice = -1;
	std::vector<double> mRotatedPoints[2];
	std::vector<double> mRotationMatrix;
public:
	DEFINE_EVENT_TYPE(image_transform, std::vector<double>&, std::vector<double>&);
	DEFINE_EVENT_TYPE(vector, std::vector<double>&);
	image_transform_event_type ImageTransformCalculatedEvent;
	vector_event_type PlaneCenterCalculatedEvent;
	MinimalSurfaceEstimator();
	~MinimalSurfaceEstimator();
	std::vector<std::vector<Vec3<int>>> mMinimalPaths[2];
	sitk::Image GetTransportSliceFromPoints(sitk::Image, Vec3<double>, Vec3<double>);
	void CalculateEikonal(sitk::Image, Vec3<double>, Vec3<double>);
	sitk::Image CalculateEikonalAndTransportInit(sitk::Image, sitk::Image, Vec3<double>, Vec3<double>);
    sitk::Image CalculateEikonalAndTransportInit(sitk::Image, Vec3<double>, Vec3<double>);
    void CalculateTransportFunction(sitk::Image& initial_contour, int maxIterations);
	void SetTransportInitSlice(const sitk::Image&);
    void Calculate(sitk::Image, sitk::Image, Vec3<double> point1, Vec3<double> point2, int maxIterations=10000);
    void Calculate(sitk::Image, Vec3<double> point1, Vec3<double> point2, int maxIterations=10000);
	void HookStageInitializationEvent(eStageEnum stage, const basic_callback_type& handler);
	void HookStageInitializationEvent(eStageEnum stage, basic_callback_type&& handler);
	void HookStageFinishedEvent(eStageEnum stage, const basic_callback_type& handler);
	void HookStageFinishedEvent(eStageEnum stage, basic_callback_type&& handler);
	void HookStageIterationEvent(eStageEnum stage, const iter_callback_type& handler);
	void HookStageIterationEvent(eStageEnum stage, iter_callback_type&& handler);
	void HookStageUpdatedEvent(eStageEnum stage, const data_callback_type& handler, eDataEnum dataType = Default);
	void HookStageUpdatedEvent(eStageEnum stage, data_callback_type&& handler, eDataEnum dataType = Default);
	void HookStageDataInitializedEvent(eStageEnum stage, const data_callback_type& handler, eDataEnum dataType = Default);
	void HookStageDataInitializedEvent(eStageEnum stage, data_callback_type&& handler, eDataEnum dataType = Default);
	void HookCalculatedPhiMapEvent(const data_callback_type& handler);
	void HookCalculatedPhiMapEvent(data_callback_type&& handler);
	void HookUpdateMeetingPlaneEvent(const AreaEikonalES::plane_update_callback_type& handler);
	void HookUpdateMeetingPlaneEvent(AreaEikonalES::plane_update_callback_type&& handler);
	void SetUsesCorrection(bool useCorrection);
	void SetUsingMeetingPoints(bool);
	const AreaEikonalES& GetAreaEikonal() const;
	const AreaEikonalES& GetRotatedAreaEikonal() const;
	const PlanePhaseFieldES& GetInitialContourCalculator() const;
	const TransportFunctionES& GetTransportFunctionCalculator() const;
	const sitk::Image& GetCombinedDistanceMap() const;
	sitk::Image GetTempInitContour() const;
	// The distance maps of the two points in the rotated volume, where the meeting plane is the
	// slice GetTransportInitPlaneSlice() along the first axis.
	const sitk::Image& GetRotatedDistanceMap(int i) const;
	int GetTransportInitPlaneSlice() const;
	// The point (0 or 1) as a continuous index in the rotated volume.
	std::vector<double> GetRotatedPoint(int i) const;
	// `image` (in the space of the input volume) resampled onto the meeting plane, as a 2D image.
	// Only valid after CalculateEikonal.
	sitk::Image SampleInitPlane(const sitk::Image& image, sitk::InterpolatorEnum interpolator);
	DEFINE_HOOK(ImageTransformCalculatedEvent, image_transform);
	DEFINE_HOOK(PlaneCenterCalculatedEvent, vector);
	sitk::Image CalculateInitialContourFromPhaseField(sitk::Image&, sitk::Image&);
	void SetInitialContourCalculatorFunc(std::function<sitk::Image(sitk::Image&, sitk::Image&)>);
	void SetInitialContourMethod(InitialContourMethod);
	void SetInitialContourFromDistance(bool from_distance);
	// The propagation of the transport function settles at the bottom of the valley of the distances,
	// whatever the contour in the meeting plane was, so a contour that is not there (a segmentation
	// of the edge of the object in the plane, say) leaves a notch in the surface at the plane. With
	// this on, the propagation starts from the contour at the bottom of the valley instead, and the
	// surface is then moved by the distance between that and the given contour, on each ray from the
	// middle of the plane, so that it passes exactly through the given contour in the plane and is
	// smooth through it. The transport function is then a signed distance (clamped to +-25).
	void SetAlignSurface(bool align);
	// The shift of the surface on each ray from the middle of the plane; empty if nothing is aligned.
	const std::vector<double>& GetSurfaceOffsets() const { return mSurfaceOffsets; }
	bool PlaneSeed(double& seed_x, double& seed_y, double& chord) const;
	// The transport function, from the contour in the meeting plane (with its zero level at the contour).
	// exact: the zero level of the plane is at the edge of the contour, not half a voxel inside of it.
	void PropagateFrom(const sitk::Image& contour, int maxIterations, bool exact);
	// The radius, on each ray from the seed, of the surface that the propagation settled at away from
	// the plane, extrapolated to the plane; false if it cannot be told.
	bool SettledRadii(double seed_x, double seed_y, double chord, std::vector<double>& radii) const;
	void ApplySurfaceOffsets();
	// The contour from the distance slice on the meeting plane; empty if there is none to be found.
	sitk::Image CalculateInitialContourFromDistance(const sitk::Image& distance_slice, bool at_floor);
};