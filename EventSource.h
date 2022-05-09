#pragma once
#include "szevent.h"
#include "SimpleITK.h"
namespace sitk = itk::simple;
class EventSource {
public:
	//DEFINE_EVENT_TYPE(data, const sitk::Image&);
	DEFINE_EVENT_TYPE(data, sitk::Image&, int);
	DEFINE_EVENT_TYPE(iter, int);
	DEFINE_EVENT_TYPE(basic);
protected:
	iter_event_type IterationEvent;
	basic_event_type InitializationEvent;
	basic_event_type FinishedEvent;
public:
	DEFINE_HOOK(IterationEvent, iter);
	DEFINE_HOOK(InitializationEvent, basic);
	DEFINE_HOOK(FinishedEvent, basic);
};