#pragma once
#include "commontype.h"
#include "ImageOp.h"

#include <vector>
#include <sitkImage.h>
#include <stack>
#include "PlanePhaseField.h"
#include "Vec.h"
#include "szevent.h"
#include "Utils.h"


//#define DEBUG_CURVATURE
#define DEBUG_STATES
#define ITERATE_ROTATED

namespace sitk = itk::simple;

#define MAXMINPATH 500 
#define INIT_DIST init_dist

#ifndef PI
#define PI 3.1415926536
#endif PI


class AreaEikonal
{
protected:
	Vec3<double> mStartPoint;
	Vec3<double> mEndPoint;

	sitk::Image mDistanceMap[2];
	sitk::Image mCombinedDistanceMap;
	sitk::Image mPhiMap;
	sitk::Image mMeetingPointsMap;
	sitk::Image mPhaseFieldMap[2];
	sitk::Image mCurvatureMap[2];

	sitk::Image mSampleImage;
	
#ifdef USE_VECTOR_AS_SET
	std::vector<POINT3D> mActivePoints[2];
	std::vector<POINT3D> mInactivePoints[2];
	std::vector<std::pair<POINT3D, double>> mVelocities[2];
#else
	std::unordered_set<POINT3D_SET> mActivePoints[2];
	std::unordered_set<POINT3D_SET> mInactivePoints[2];
	std::unordered_map<POINT3D_MAP(double)> mVelocities[2];
#endif

	std::unordered_set<POINT3D_SET> mMeetingPoints;
	int mMeetingPointsCount;
	Vec3<double> mMeetingPlaneCenter;
	Vec3<double> mMeetingPlaneNormal;
	int mIterationCount;
	int mMeetingPointsMin;
	bool mPlaneFinalized;
	double mMeetingPlaneOffset;
	bool mUsesCorrection;
	double mCurrentDistance;
#ifdef DEBUG_CURVATURE
	sitk::Image m_new_update[2];
#endif
	int GetSmallestPlaneSize();
	void InitializeContainers(const sitk::Image& image);
	void InitializeDistanceMap(int idx, Vec3<double> center, int initialRadius);
	void InitializeNeighbors();
	virtual void UpdateCurvature(int i);
	double GetMinS(int i) const ;
	double UpdateVelo(int i, double S);
	virtual void UpdateField(int i, double maxv);
	virtual void UpdateDistance(int i, double current_distance);
	virtual void UpdateMeetPoints();
	virtual void FitMeetingPlane();

public:
	virtual void Initialize(const sitk::Image& image, Vec3<double>& start_point, Vec3<double>& end_point, double beta, double alpha);
	AreaEikonal Rotate(std::vector<double>& rotation_matrix, bool inverse = false) const;
	bool IsDone() const;
	virtual void UpdateMeetingPlane();
	void CombineDistance();
	void CombineDistance(int slice, double p1_x, double p2_x);
	void SmoothDistances();
	virtual void Iterate();
	void Calculate(const sitk::Image& image, Vec3<double>& point1, Vec3<double>& point2, double beta, double alpha);
	virtual void Calculate();
	Vec3<double> GetMeetingPlaneCenter() const;
	Vec3<double> GetMeetingPlaneNormal() const;
	double GetMeetingPlaneOffset() const;
	double GetCurrentDistance() const;
	const sitk::Image& GetDistanceMap(int i) const;
	const sitk::Image& GetCombinedDistanceMap() const;
	const sitk::Image& GetPhiMap() const;
	const sitk::Image& GetPhaseFieldMap(int i) const;
	const sitk::Image& GetMeetingPointsMap() const;
	const sitk::Image& GetSampleImage() const;
	std::vector<double> TransformPhysicalPointToContinuousIndex(const std::vector<double>& point) const;
	std::vector<double> TransformContinuousIndexToPhysicalPoint(const std::vector<double>& point) const;
	void SetUsesCorrection(bool useCorrection);
	std::vector<Vec3<int>> ResolvePath(Vec3<int> point, int idx) const;
};



class AreaEikonalES : public AreaEikonal, public EventSource {
public:
	DEFINE_EVENT_TYPE(plane_update, Vec3<double>, Vec3<double>);
private:
	EVENT_TYPE(plane_update) UpdateMeetingPlaneEvent;
	EVENT_TYPE(data) UpdatedDistanceMapEvent;//CALLED
	EVENT_TYPE(data) CombinedDistanceMapEvent;
	EVENT_TYPE(data) CalculatedPhiMapEvent;
	EVENT_TYPE(data) UpdatedMeetingPointsMapEvent;//CALLED
	EVENT_TYPE(data) UpdatedPhaseFieldMapEvent;//CALLED
	EVENT_TYPE(data) UpdatedCurvatureMapEvent;//CALLED
	EVENT_TYPE(data) InitializedMeetingPointsMapEvent;
	EVENT_TYPE(data) InitializedPhaseFieldMapEvent;
	EVENT_TYPE(data) InitializedDistanceMapEvent;
	EVENT_TYPE(data) InitializedCurvatureMapEvent;

	void UpdateCurvature(int i) override {
		AreaEikonal::UpdateCurvature(i);
		UpdatedCurvatureMapEvent(mCurvatureMap[i], i);
	}

	void UpdateField(int i, double maxv) override {
		AreaEikonal::UpdateField(i, maxv);
		UpdatedPhaseFieldMapEvent(mPhaseFieldMap[i], i);
	}

	void UpdateDistance(int i, double current_distance) override {
		AreaEikonal::UpdateDistance(i, current_distance);
	}

	void UpdateMeetPoints() override {
		AreaEikonal::UpdateMeetPoints();
		UpdatedMeetingPointsMapEvent(mMeetingPointsMap, 0);
	}
public:
	using AreaEikonal::Calculate;

	AreaEikonalES(AreaEikonal eikonal) : AreaEikonal(eikonal) {
		InitializedCurvatureMapEvent(mCurvatureMap[0], 0);
		InitializedCurvatureMapEvent(mCurvatureMap[1], 1);
		InitializedDistanceMapEvent(mDistanceMap[0], 0);
		InitializedDistanceMapEvent(mDistanceMap[1], 1);
		InitializedMeetingPointsMapEvent(mMeetingPointsMap, 0);
		InitializedPhaseFieldMapEvent(mPhaseFieldMap[0], 0);
		InitializedPhaseFieldMapEvent(mPhaseFieldMap[1], 1);
	}

	AreaEikonalES& operator=(const AreaEikonal&& areaEikonal) {
		AreaEikonal::operator=(areaEikonal);
		InitializedCurvatureMapEvent(mCurvatureMap[0], 0);
		InitializedCurvatureMapEvent(mCurvatureMap[1], 1);
		InitializedDistanceMapEvent(mDistanceMap[0], 0);
		InitializedDistanceMapEvent(mDistanceMap[1], 1);
		InitializedMeetingPointsMapEvent(mMeetingPointsMap, 0);
		InitializedPhaseFieldMapEvent(mPhaseFieldMap[0], 0);
		InitializedPhaseFieldMapEvent(mPhaseFieldMap[1], 1);
		return *this;
	}
	AreaEikonalES() {}
	void UpdateMeetingPlane() override {
		int prev_point_count = mMeetingPointsCount >= mMeetingPointsMin ? mMeetingPointsCount : mMeetingPointsMin;
		UpdateMeetPoints();
		if (mMeetingPointsCount > prev_point_count) {
			FitMeetingPlane();
			UpdateMeetingPlaneEvent(mMeetingPlaneCenter, mMeetingPlaneNormal);
		}
		else if (mMeetingPointsCount == prev_point_count && mMeetingPointsCount > 0) {
			mPlaneFinalized = true;
		}
	}
	void Initialize(const sitk::Image& image, Vec3<double>& start_point, Vec3<double>& end_point, double beta, double alpha) override {
		AreaEikonal::Initialize(image, start_point, end_point, beta, alpha);
		CalculatedPhiMapEvent(mPhiMap, 0);
		InitializedPhaseFieldMapEvent(mPhaseFieldMap[0], 0);
		InitializedPhaseFieldMapEvent(mPhaseFieldMap[1], 1);
		InitializedDistanceMapEvent(mDistanceMap[0], 0);
		InitializedDistanceMapEvent(mDistanceMap[1], 1);
		InitializedCurvatureMapEvent(mCurvatureMap[0], 0);
		InitializedCurvatureMapEvent(mCurvatureMap[1], 1);
		InitializedMeetingPointsMapEvent(mMeetingPointsMap, 0);
		InitializationEvent();
	}
	void Iterate() override {
		AreaEikonal::Iterate();
		UpdatedDistanceMapEvent(mDistanceMap[0], 0);
		UpdatedDistanceMapEvent(mDistanceMap[1], 1);
		IterationEvent(mIterationCount);
	}

	void Calculate() override {
		AreaEikonal::Calculate();
		FinishedEvent();
	}

	void CombineDistance() {
		AreaEikonal::CombineDistance();
		CombinedDistanceMapEvent(mCombinedDistanceMap, 0);
	}

	void CombineDistance(int slice, double p1_x, double p2_x) {
		AreaEikonal::CombineDistance(slice, p1_x, p2_x);
		CombinedDistanceMapEvent(mCombinedDistanceMap, 0);
	}

	DEFINE_HOOK(InitializedPhaseFieldMapEvent, data);
	DEFINE_HOOK(InitializedDistanceMapEvent, data);
	DEFINE_HOOK(InitializedMeetingPointsMapEvent, data);
	DEFINE_HOOK(InitializedCurvatureMapEvent, data);
	DEFINE_HOOK(UpdateMeetingPlaneEvent, plane_update);
	DEFINE_HOOK(UpdatedDistanceMapEvent, data);
	DEFINE_HOOK(UpdatedPhaseFieldMapEvent, data);
	DEFINE_HOOK(UpdatedCurvatureMapEvent, data);
	DEFINE_HOOK(UpdatedMeetingPointsMapEvent, data);
	DEFINE_HOOK(CalculatedPhiMapEvent, data);
	DEFINE_HOOK(CombinedDistanceMapEvent, data);
};
