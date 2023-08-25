#pragma once
#include "PlanePhaseField.h"
#include "AreaEikonal.h"
#include "Vec.h"
#include "SimpleITK.h"
#include "szevent.h"
#include <functional>



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

namespace sitk = itk::simple;
class MinimalSurfaceEstimator: public EventSource
{
	AreaEikonalES mAreaEikonal;
	AreaEikonalES mRotatedAreaEikonal;
	PlanePhaseFieldES mInitialContourCalculator;
	TransportFunctionES mTransportFunctionCalculator;
	std::function<sitk::Image(sitk::Image&, sitk::Image&)> mInitialContourCalculatorFunc;
public:
	DEFINE_EVENT_TYPE(image_transform, std::vector<double>&, std::vector<double>&);
	DEFINE_EVENT_TYPE(vector, std::vector<double>&);
	image_transform_event_type ImageTransformCalculatedEvent;
	vector_event_type PlaneCenterCalculatedEvent;
	MinimalSurfaceEstimator();
	~MinimalSurfaceEstimator();
	std::vector<std::vector<Vec3<int>>> mMinimalPaths[2];
	void Calculate(sitk::Image, Vec3<double> point1, Vec3<double> point2, double beta, double alpha = 0.01, int maxIterations=10000);
	void HookStageInitializationEvent(eStageEnum stage, basic_callback_type& handler);
	void HookStageFinishedEvent(eStageEnum stage, basic_callback_type& handler);
	void HookStageIterationEvent(eStageEnum stage, iter_callback_type& handler);
	void HookStageUpdatedEvent(eStageEnum stage, data_callback_type& handler, eDataEnum dataType = Default);
	void HookStageDataInitializedEvent(eStageEnum stage, data_callback_type& handler, eDataEnum dataType = Default);
	void HookCalculatedPhiMapEvent(data_callback_type& handler);
	void HookUpdateMeetingPlaneEvent(AreaEikonalES::plane_update_callback_type& handler);
	void SetUsesCorrection(bool useCorrection);
	const AreaEikonalES& GetAreaEikonal() const;
	const AreaEikonalES& GetRotatedAreaEikonal() const;
	const PlanePhaseFieldES& GetInitialContourCalculator() const;
	const TransportFunctionES& GetTransportFunctionCalculator() const;
	DEFINE_HOOK(ImageTransformCalculatedEvent, image_transform);
	DEFINE_HOOK(PlaneCenterCalculatedEvent, vector);
	sitk::Image CalculateInitialContourFromPhaseField(sitk::Image&, sitk::Image&);
	void SetInitialContourCalculatorFunc(std::function<sitk::Image(sitk::Image&, sitk::Image&)>);
};