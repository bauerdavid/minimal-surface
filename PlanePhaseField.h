#pragma once
#include "SimpleITK.h"
#include <unordered_set>
#include "szevent.h"
#include "EventSource.h"
namespace sitk = itk::simple;

class PlanePhaseField
{
protected:
	int mPositivePhaseCount = 0;
	int mIterationCount = 0;
	int mIterationMax = 20000;
	int mLastUpdate = 0;
	int mUpdateTolerance = 200;
	bool mDone = false;
	sitk::Image mGradX;
	sitk::Image mGradY;
	sitk::Image mPhaseField;
	sitk::Image mVelocity;
public:
	PlanePhaseField() {}

	// get the distance gradients along the selected x slice from the two neighboring slices
	void Initialize(const sitk::Image& imageSlice, const sitk::Image& distanceSlice);
	void Iterate();
	void Calculate(const sitk::Image& imageSlice, const sitk::Image& distanceSlice);
	bool IsDone() const;
	const sitk::Image& GetGradX() const;
	const sitk::Image& GetGradY() const;
	const sitk::Image& GetPhaseField() const;
	// collects bounding points into a set
	std::unordered_set<unsigned long> RetrieveBound() const;
};

class PlanePhaseFieldES : public PlanePhaseField, public EventSource {
	EVENT_TYPE(data) InitializedPhaseFieldMapEvent;
	EVENT_TYPE(data) UpdatedPhaseFieldMapEvent;
	EVENT_TYPE(data) InitializedGradXMapEvent;
	EVENT_TYPE(data) InitializedGradYMapEvent;
public:
	void Initialize(const sitk::Image& imageSlice, const sitk::Image& distanceSlice) {
		PlanePhaseField::Initialize(imageSlice, distanceSlice);
		InitializedPhaseFieldMapEvent(mPhaseField, 0);
		InitializedGradXMapEvent(mGradX, 0);
		InitializedGradYMapEvent(mGradY, 0);
		InitializationEvent();
	}
	void Iterate() {
		PlanePhaseField::Iterate();
		UpdatedPhaseFieldMapEvent(mPhaseField, 0);
		IterationEvent(mIterationCount);
	}
	void Calculate(const sitk::Image& imageSlice, const sitk::Image& distanceSlice) {
		Initialize(imageSlice, distanceSlice);
		while (!IsDone())
			Iterate();
		FinishedEvent();
	}
	/*
	void Calculate(const sitk::Image& imageSlice, const sitk::Image& distanceSlice) {
		PlanePhaseField::Calculate(imageSlice, distanceSlice);
		InitializationEvent();
		InitializedPhaseFieldMapEvent(mPhaseField, 0);
		IterationEvent(0);
		UpdatedPhaseFieldMapEvent(mPhaseField, 0);
		FinishedEvent();
	}*/

	DEFINE_HOOK(InitializedPhaseFieldMapEvent, data);
	DEFINE_HOOK(UpdatedPhaseFieldMapEvent, data);
	DEFINE_HOOK(InitializedGradXMapEvent, data);
	DEFINE_HOOK(InitializedGradYMapEvent, data);
};