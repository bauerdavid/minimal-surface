#pragma once
#include "commontype.h"
#include <vector>
#include <unordered_set>
#include "Transport.h"
#include "Vec.h"
#include <SimpleITK.h>

#define XS_ 220 
#define YS_ 155 
#define ZS_ 110 
namespace sitk = itk::simple;
class CImageOp
{
public:
	CImageOp(void);
	~CImageOp(void);
public:

	sitk::Image CreateTestImage(); // 9,19
	sitk::Image CreateTestImage2(); // 9,19
	sitk::Image Blur(const sitk::Image& image);
};
