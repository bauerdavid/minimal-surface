#pragma once
#include "commontype.h"
#define _NO_BOU_ 1e11
#include "SimpleITK.h"
#include "EventSource.h"
#include "Utils.h"
namespace sitk = itk::simple;

enum {Talox, Taloy, Taloz};

class TransportFunction
{
protected:
	int mIterationsMax = 100000;
	int mIterationCount = 0;
	/*
	Boundary state
	 1: read-only
	 0: calculation is needed
	 -1: not initialized
	*/
	sitk::Image mTransportFunctionMap[2];
	sitk::Image mReadOnlyMap;
	sitk::Image mGradX;
	sitk::Image mGradY;
	sitk::Image mGradZ;
	sitk::Image& GetTransportFunction();
public:
	TransportFunction():m_active(0) {}

	int m_active;
	// Initialize gradients fro, distance map, transform function from inimap, and boundaries (considering inimap, and setting data edges as boundaries)
	void Initialize(sitk::Image distanceMap, const sitk::Image& initialSlice, int initialSliceIdx, realnum distanceMax);
	void Calculate(sitk::Image distanceMap, const sitk::Image& initialSlice, int initialSliceIdx, realnum distanceMax);
	void RotateTransportFunction(const std::vector<double>& rotation_matrix, std::vector<unsigned> sample_size);
	const sitk::Image& GetTransportFunction() const;
	
	const sitk::Image& GetReadOnlyMap() const;
	const sitk::Image& GetGradX() const;
	const sitk::Image& GetGradY() const;
	const sitk::Image& GetGradZ() const;
protected:
	bool Iterate();
};

class TransportFunctionES : public TransportFunction, public EventSource {
public:
	EVENT_TYPE(data) UpdatedTransportFunctionMapEvent;
	EVENT_TYPE(data) InitializedTransportFunctionMapEvent;
	void Initialize(sitk::Image distanceMap, const sitk::Image& initialSlice, int initialSliceIdx, realnum distanceMax) {
		TransportFunction::Initialize(distanceMap, initialSlice, initialSliceIdx, distanceMax);
		InitializedTransportFunctionMapEvent(mTransportFunctionMap[0], 0);
		InitializedTransportFunctionMapEvent(mTransportFunctionMap[1], 1);
		InitializationEvent();
	}
	bool Iterate() {
		bool updated = TransportFunction::Iterate();
		UpdatedTransportFunctionMapEvent(GetTransportFunction(), mIterationCount%2);
		IterationEvent(mIterationCount);
		return updated;
	}
	void Calculate(sitk::Image distanceMap, const sitk::Image& initialSlice, int initialSliceIdx, realnum distanceMax)
	{
		_PROFILING;
		Initialize(distanceMap, initialSlice, initialSliceIdx, distanceMax);
		if (!m_active) return;

		for (int ii = 0; ii < mIterationsMax; ++ii) {
			bool updated = Iterate();
			if (!updated)
				break;
		}
		FinishedEvent();
	}

	DEFINE_HOOK(UpdatedTransportFunctionMapEvent, data);
	DEFINE_HOOK(InitializedTransportFunctionMapEvent, data);
};

