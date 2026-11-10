#include "PlanePhaseField.h"
#include "sitkImageOperators.h"

using namespace std;

void PlanePhaseField::Initialize(const sitk::Image& imageSlice, const sitk::Image& distanceSlice)
{
	mPositivePhaseCount = 0;
	mIterationCount = 0;
	mIterationMax = 20000;
	mLastUpdate = 0;
	mUpdateTolerance = 200;
	mDone = false;
	sitk::Image normalized_distance = SmartSigmoid(distanceSlice);

	vector<unsigned> size = normalized_distance.GetSize();
	int xs = size[0], ys = size[1];
	sitk::Image img_gradient = sitk::GradientMagnitudeRecursiveGaussian(imageSlice);
	sitk::Image imgGradientInv = SmartSigmoid(sitk::Cast(1 - img_gradient, sitk::sitkFloat64));
	mGradX = 1000 * sitk::Derivative(normalized_distance, 0, 1, false) * imgGradientInv;
	mGradY = 1000 * sitk::Derivative(normalized_distance, 1, 1, false) * imgGradientInv;
	mPhaseField = sitk::Image(xs, ys, sitk::sitkFloat64) + 1;
	DrawEllipse<sitk::sitkFloat64>(mPhaseField, xs / 2 - 5, ys / 2 - 5, xs / 2, ys / 2, -1.);
	double* field_buffer = mPhaseField.GetBufferAsDouble();
	mVelocity = sitk::Image(xs, ys, sitk::sitkFloat64);


	for (int yy = 0; yy < ys; ++yy) {
		for (int xx = 0; xx < xs; ++xx) {
			//if current point is not on the edge of the data
			if(BUF_IDX2D(field_buffer, xs, ys, xx, yy) > 0)
				++mPositivePhaseCount;
		}
	}
}

void PlanePhaseField::Iterate()
{
	_PROFILING;

	// One step of the phase field. The derivatives of the field are taken in place from its
	// neighbours (with the border repeated, like the zero-flux boundary of the SimpleITK filters
	// that this used to call), so nothing is allocated: a step used to allocate three images,
	// which cost more than the arithmetic.
	int xs = mPhaseField.GetWidth(), ys = mPhaseField.GetHeight();
	double* velo_buffer = mVelocity.GetBufferAsDouble();
	double* phasefield_buffer = mPhaseField.GetBufferAsDouble();
	const double* gx_buffer = mGradX.GetBufferAsDouble();
	const double* gy_buffer = mGradY.GetBufferAsDouble();

	realnum maxv(0);
OMP_PARALLEL
	{
		double temp_max = 0;
OMP_FOR
		for (int y = 0; y < ys; y++) {
			const double* row = phasefield_buffer + (size_t)y * xs;
			const double* up = phasefield_buffer + (size_t)(y > 0 ? y - 1 : 0) * xs;
			const double* down = phasefield_buffer + (size_t)(y < ys - 1 ? y + 1 : ys - 1) * xs;
			for (int x = 0; x < xs; x++) {
				const double phase = row[x];
				const double left = row[x > 0 ? x - 1 : 0], right = row[x < xs - 1 ? x + 1 : xs - 1];
				const double lapl = left + right + up[x] + down[x] - 4 * phase;
				const double field_gx = 0.5 * (right - left), field_gy = 0.5 * (down[x] - up[x]);
				const size_t yx = (size_t)y * xs + x;
				const double velo = lapl * 4 + phase - phase * phase * phase + (gx_buffer[yx] * field_gx + gy_buffer[yx] * field_gy);
				velo_buffer[yx] = velo;
				if (fabs(velo) > temp_max)
					temp_max = fabs(velo);
			}
		}
OMP_CRITICAL_NO_TAG
		{
			if (temp_max > maxv)
				maxv = temp_max;
		}
	}
	if (maxv <= 1e-5f)
		maxv = 1;
	int positive_count(0);

OMP_PARALLEL_FOR_REDUCTION(+, positive_count)
	for (int yx = 0; yx < ys * xs; ++yx) {
		phasefield_buffer[yx] += velo_buffer[yx] / (2 * maxv);
		if (phasefield_buffer[yx] < -1) phasefield_buffer[yx] = -1;
		else if (phasefield_buffer[yx] > 1) phasefield_buffer[yx] = 1;
		if (phasefield_buffer[yx] > 0) ++positive_count;
	}
	if (++mIterationCount >= mIterationMax) {
		mDone = true;
	}
	else if (positive_count - mPositivePhaseCount > 5) {
		mLastUpdate = mIterationCount;
		mPositivePhaseCount = positive_count;
	}
	else if (mIterationCount - mLastUpdate >= 200) {
		mDone = true;
	}
}


void PlanePhaseField::Calculate(const sitk::Image& imageSlice, const sitk::Image& distanceSlice) {
	Initialize(imageSlice, distanceSlice);
	while (!IsDone())
		Iterate();
}

bool PlanePhaseField::IsDone() const {
	return mDone;
}

const sitk::Image& PlanePhaseField::GetGradX() const {
	return mGradX;
}
const sitk::Image& PlanePhaseField::GetGradY() const {
	return mGradY;
}
const sitk::Image& PlanePhaseField::GetPhaseField() const {
	return mPhaseField;
}

std::unordered_set<unsigned long> PlanePhaseField::RetrieveBound() const
{
	_PROFILING;
	int xs = mPhaseField.GetWidth(), ys = mPhaseField.GetHeight();
	const double* field_buffer = mPhaseField.GetBufferAsDouble();

	std::unordered_set<unsigned long> bound;
	for (int yy = 1; yy < ys - 1; ++yy) {
		for (int xx = 1; xx < xs - 1; ++xx) {

			if (BUF_IDX2D(field_buffer, xs, ys, xx, yy) > 0) {
				if (BUF_IDX2D(field_buffer, xs, ys, xx, yy + 1) <= 0 || BUF_IDX2D(field_buffer, xs, ys, xx, yy - 1) <= 0
					|| BUF_IDX2D(field_buffer, xs, ys, xx + 1, yy) <= 0 || BUF_IDX2D(field_buffer, xs, ys, xx - 1, yy) <= 0)
					bound.emplace((yy << 16) + xx);
			}

		}
	}
	return bound;
}