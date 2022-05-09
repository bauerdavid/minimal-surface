#include "PlanePhaseField.h"
#include "Utils.h"
#include "sitkImageOperators.h"

using namespace std;

void PlanePhaseField::Initialize(const sitk::Image& distanceSlice)
{
	mPositivePhaseCount = 0;
	mIterationCount = 0;
	mIterationMax = 20000;
	mLastUpdate = 0;
	mUpdateTolerance = 200;
	mDone = false;
	vector<unsigned> size = distanceSlice.GetSize();
	int xs = size[0], ys = size[1];
	const double* distance_buffer = distanceSlice.GetBufferAsDouble();

	mGradX = sitk::Derivative(distanceSlice, 0, 1, false);
	mGradY = sitk::Derivative(distanceSlice, 1, 1, false);
	mPhaseField = sitk::Image(xs, ys, sitk::sitkFloat64) + 1;
	double* field_buffer = mPhaseField.GetBufferAsDouble();
	mVelocity = sitk::Image(xs, ys, sitk::sitkFloat64);


	for (int yy = 0; yy < ys; ++yy) {
		for (int xx = 0; xx < xs; ++xx) {
			//if current point is not on the edge of the data
			if (!(yy < 5 || yy >= ys - 5 || xx < 5 || xx >= xs - 5))
				BUF_IDX2D(field_buffer, xs, ys, xx, yy) = -1.0f;
			else ++mPositivePhaseCount;
		}
	}
}

void PlanePhaseField::Iterate()
{
	_PROFILING;

	int xs = mPhaseField.GetWidth(), ys = mPhaseField.GetHeight();
	sitk::Image field_grad_x = sitk::Derivative(mPhaseField, 0, 1, false);
	sitk::Image field_grad_y = sitk::Derivative(mPhaseField, 1, 1, false);
	sitk::Image laplacian = sitk::Laplacian(mPhaseField, false);
	double* velo_buffer = mVelocity.GetBufferAsDouble();
	double* laplacian_buffer = laplacian.GetBufferAsDouble();
	double* phasefield_buffer = mPhaseField.GetBufferAsDouble();
	double* field_gx_buffer = field_grad_x.GetBufferAsDouble();
	double* field_gy_buffer = field_grad_y.GetBufferAsDouble();
	double* gx_buffer = mGradX.GetBufferAsDouble();
	double* gy_buffer = mGradY.GetBufferAsDouble();
	//sitk::Image velocity = sitk::Laplacian(mPhaseField, false) * 4 + mPhaseField - sitk::Pow(mPhaseField, 3) + 10 * (mGradX * field_grad_x + mGradY * field_grad_y);

	realnum maxv(0);
OMP_PARALLEL
	{
		double temp_max = 0;
OMP_FOR
		for (int yx = 0; yx < xs * ys; yx++) {
			int y = yx / xs;
			int x = yx % xs;
			double& phase = BUF_IDX2D(phasefield_buffer, xs, ys, x, y);
			double& lapl = BUF_IDX2D(laplacian_buffer, xs, ys, x, y);
			double& velo = BUF_IDX2D(velo_buffer, xs, ys, x, y);
			double& field_gx = BUF_IDX2D(field_gx_buffer, xs, ys, x, y);
			double& field_gy = BUF_IDX2D(field_gy_buffer, xs, ys, x, y);
			double& gx = BUF_IDX2D(gx_buffer, xs, ys, x, y);
			double& gy = BUF_IDX2D(gy_buffer, xs, ys, x, y);
			velo = lapl * 4 + phase - pow(phase, 3) + 10 * (gx * field_gx + gy * field_gy);
			if (velo > temp_max)
				temp_max = velo;
		}
		if (temp_max > maxv) {
			maxv = temp_max;
		}
	}

	if (maxv > 1e-5f)
		maxv = 0.025f / maxv;
	int positive_count(0);

OMP_PARALLEL_FOR_REDUCTION(+, positive_count)
	for (int yx = 0; yx < ys * xs; ++yx) {
		phasefield_buffer[yx] += velo_buffer[yx] * maxv;
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


/*void PlanePhaseField::Iterate() {
	_PROFILING;
	vector<unsigned> size = mPhaseField.GetSize();
	int xs(size[0]), ys(size[1]);
	mVelocity = sitk::Laplacian(mPhaseField, false);
	mVelocity *= 4.0f;
	double* phase_buffer = mPhaseField.GetBufferAsDouble();
	double* velo_buffer = mPhaseField.GetBufferAsDouble();
	double* gx_buffer = mGradX.GetBufferAsDouble();
	double* gy_buffer = mGradY.GetBufferAsDouble();
	realnum fac = 1.0f;
	realnum maxv(0.0f);
#pragma omp parallel
	{
		double temp_maxv(0);
#pragma omp for
		for (int yx = 0; yx < xs * ys; yx++) {
			int yy = yx / xs;
			if (yy > 0 && yy < ys - 1) {
				int xx = yx % xs;
				if (xx > 0 && xx < xs - 1) {
					realnum fmm = BUF_IDX2D(phase_buffer, xs, ys, xx, yy);
					BUF_IDX2D(velo_buffer, xs, ys, xx, yy) -= fac * fmm * fmm * fmm;
					BUF_IDX2D(velo_buffer, xs, ys, xx, yy) += fac * fmm;
					realnum gfx = 0.5f * (BUF_IDX2D(phase_buffer, xs, ys, xx + 1, yy) - BUF_IDX2D(phase_buffer, xs, ys, xx - 1, yy));
					realnum gfy = 0.5f * (BUF_IDX2D(phase_buffer, xs, ys, xx, yy+1) - BUF_IDX2D(phase_buffer, xs, ys, xx, yy));

					BUF_IDX2D(velo_buffer, xs, ys, xx, yy) += 10.0f * (gfx * BUF_IDX2D(gx_buffer, xs, ys, xx, yy) + gfy * BUF_IDX2D(gy_buffer, xs, ys, xx, yy));
					realnum cv = BUF_IDX2D(velo_buffer, xs, ys, xx, yy);
					if (cv < 0)
						cv *= -1;
					if (cv > temp_maxv)
						temp_maxv = cv;
				}
			}
		}
#pragma omp critical
		if (temp_maxv > maxv)
			maxv = temp_maxv;
	}

	if (maxv > 1e-5f)
		maxv = 0.025f / maxv;
	int positive_count(0);

#pragma omp parallel for reduction(+: positive_count)
	for (int yy = 1; yy < ys - 1; ++yy) {
		for (int xx = 1; xx < xs - 1; ++xx) {
			BUF_IDX2D(phase_buffer, xs, ys, xx, yy) += BUF_IDX2D(velo_buffer, xs, ys, xx, yy) * maxv;
			if (BUF_IDX2D(phase_buffer, xs, ys, xx, yy) > 0) ++positive_count;
		}
	}
	if (++mIterationCount >= mIterationMax)
		mDone = true;
	else if (positive_count > mPositivePhaseCount)
		mLastUpdate = mIterationCount;
	else if (mIterationCount - mLastUpdate >= 200)
		mDone = true;

	mPositivePhaseCount = positive_count;
}*/

void PlanePhaseField::Calculate(const sitk::Image& distanceSlice) {
	Initialize(distanceSlice);
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