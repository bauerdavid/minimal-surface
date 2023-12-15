#pragma once
#include "commontype.h"
#define _NO_BOU_ 1e11
#include "SimpleITK.h"
#include "EventSource.h"
#include "Utils.h"
#include "IndexedPriorityQueue.h"
//#define ORDERED_TRAVERSAL

namespace sitk = itk::simple;

enum {Talox, Taloy, Taloz};


template <class _Ty = void, class Comparator = std::greater<_Ty>>
struct abs_comparator {
	typedef _Ty _FIRST_ARGUMENT_TYPE_NAME;
	typedef _Ty _SECOND_ARGUMENT_TYPE_NAME;
	typedef bool _RESULT_TYPE_NAME;
	static constexpr Comparator comp{};
	_NODISCARD constexpr bool operator()(const _Ty& _Left, const _Ty& _Right) const {
		return abs_comparator::comp(abs(_Left), abs(_Right));
	}
};

class TransportFunction
{
protected:
	int mIterationsMax = 10000;
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
	sitk::Image& GetTransportFunction(int);
	sitk::Image& GetTransportFunction();
#ifdef ORDERED_TRAVERSAL
	IndexedPriorityQueue<POINT3D_MAP(double), std::greater<double>> mNarrowBand;
	double ComputeValue(std::vector<int> p);
	void Iterate();
	bool Finished();
#else
	bool Iterate();
#endif
public:
	TransportFunction():m_active(0) {}
	realnum GetTransportValueByDistGrad(int, int, int);
	int m_active;
	// Initialize gradients fro, distance map, transform function from inimap, and boundaries (considering inimap, and setting data edges as boundaries)
	void Initialize(sitk::Image distanceMap, const sitk::Image& initialSlice, int initialSliceIdx, realnum distanceMax);
	void Calculate(sitk::Image distanceMap, const sitk::Image& initialSlice, int initialSliceIdx, realnum distanceMax, int maxIterations);
	void RotateTransportFunction(const std::vector<double>& rotation_matrix, std::vector<unsigned> sample_size);
	const sitk::Image& GetTransportFunction(int) const;
	const sitk::Image& GetTransportFunction() const;

	const sitk::Image& GetReadOnlyMap() const;
	const sitk::Image& GetGradX() const;
	const sitk::Image& GetGradY() const;
	const sitk::Image& GetGradZ() const;
protected:
};

class TransportFunctionES : public TransportFunction, public EventSource {
public:
	EVENT_TYPE(data) UpdatedTransportFunctionMapEvent;
	EVENT_TYPE(data) InitializedTransportFunctionMapEvent;
#ifdef ORDERED_TRAVERSAL
	void Initialize(sitk::Image distanceMap, const sitk::Image& initialSlice, int initialSliceIdx, realnum distanceMax) {
		TransportFunction::Initialize(distanceMap, initialSlice, initialSliceIdx, distanceMax);
		InitializedTransportFunctionMapEvent(mTransportFunctionMap[0], 0);
		InitializationEvent();
	}
	void Iterate() {
		TransportFunction::Iterate();
		UpdatedTransportFunctionMapEvent(GetTransportFunction(0), 0);
		IterationEvent(mIterationCount);
	}
	void Calculate(sitk::Image distanceMap, const sitk::Image& initialSlice, int initialSliceIdx, realnum distanceMax) {
		Initialize(distanceMap, initialSlice, initialSliceIdx, distanceMax);
		while (!Finished()) {
			Iterate();
		}
		FinishedEvent();
	}
#else
	void Initialize(sitk::Image distanceMap, const sitk::Image& initialSlice, int initialSliceIdx, realnum distanceMax) {
		TransportFunction::Initialize(distanceMap, initialSlice, initialSliceIdx, distanceMax);
		InitializedTransportFunctionMapEvent(mTransportFunctionMap[0], 0);
		InitializedTransportFunctionMapEvent(mTransportFunctionMap[1], 1);
		InitializationEvent();
	}
	bool Iterate() {
		bool updated = TransportFunction::Iterate();
		UpdatedTransportFunctionMapEvent(GetTransportFunction(0), 0);
		IterationEvent(mIterationCount);
		return updated;
	}

	void Calculate(sitk::Image distanceMap, const sitk::Image& initialSlice, int initialSliceIdx, realnum distanceMax, int maxIterations)
	{
		_PROFILING;
		Initialize(distanceMap, initialSlice, initialSliceIdx, distanceMax);
		if (!m_active) return;
		for (int ii = 0; ii < maxIterations; ++ii) {
			bool updated = Iterate();
			if (!updated) {
				FinishedEvent();
				return;
			}
		}
		FinishedEvent();
		std::cout << "transport function did not converge" << std::endl;
	}
#endif
	


	DEFINE_HOOK(UpdatedTransportFunctionMapEvent, data);
	DEFINE_HOOK(InitializedTransportFunctionMapEvent, data);
};