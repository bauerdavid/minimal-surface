#include "PlanePhaseField.h"
#include "Utils.h"
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
			velo = lapl * 4 + phase - pow(phase, 3) + (gx * field_gx + gy * field_gy);
			if (abs(velo) > temp_max)
				temp_max = abs(velo);
		}
		if (temp_max > maxv) {
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

void PlanePhaseField::Calculate(const sitk::Image& imageSlice, const sitk::Image& distanceSlice) {
	Initialize(imageSlice, distanceSlice);
	while (!IsDone())
		Iterate();
}
/*
void PlanePhaseField::Calculate(const sitk::Image& imageSlice, const sitk::Image& distanceSlice) {
	int n_pixels = distanceSlice.GetNumberOfPixels();
	const double initial_distance = 5.0;
	const double sigma = 1.0;
	const double propagation_scaling = 0.5;
	double new_min = ImageQuantile(distanceSlice, 0.01);
	sitk::MinimumMaximumImageFilter min_max_filter;
	min_max_filter.Execute(distanceSlice);
	sitk::Image&& rescaled_distance_slice = sitk::Clamp((distanceSlice - new_min) / (min_max_filter.GetMaximum() - new_min), sitk::sitkUnknown, 0., 1.);
	
	const double Q = 0.95;
	const double eps = 0.1;
	double q = ImageQuantile(rescaled_distance_slice, Q);
	double m = ImageQuantile(rescaled_distance_slice, 0.01);
	double sigm_alpha = (q - m) / log((1 / eps - 1) / (1 / Q - 1));
	double sigm_beta = log(1 / eps - 1) * sigm_alpha + m;
	
	sitk::Image sigmoid = sitk::Cast(sitk::Sigmoid(rescaled_distance_slice, sigm_alpha, sigm_beta, 1., 0.), sitk::sitkFloat64) - 0.1;
	
	sitk::FastMarchingImageFilter fast_marching;
	unsigned seed_x = distanceSlice.GetWidth() / 2, seed_y = distanceSlice.GetHeight() / 2;
	fast_marching.AddTrialPoint({ seed_x, seed_y, 0 });
	sitk::Image&& fast_marching_output = fast_marching.Execute(sitk::Image(distanceSlice.GetSize(), sitk::sitkUInt16) + 1) - 9;
	sitk::Image&& initial_contour = sitk::GeodesicActiveContourLevelSet(fast_marching_output, sigmoid, 0.02, propagation_scaling);
	sitk::Image&& initial_mask = sitk::BinaryThreshold(initial_contour, -1000, 2.99);
	
	//sitk::Image&& smooth_image_slice = sitk::CurvatureAnisotropicDiffusion(imageSlice, 0.125, 9., 1, 10);
	sitk::Image gradient_magnitude = sitk::Cast(sitk::GradientMagnitudeRecursiveGaussian(imageSlice, sigma), sitk::sitkFloat64);
	
	double q2 = ImageQuantile(gradient_magnitude, 0.99);
	double m2 = ImageQuantile(gradient_magnitude, 0.01);
	double sigm_alpha2 = (q2 - m2) / log((1 / eps - 1) / (1 / 0.99 - 1));
	double sigm_beta2 = log(1 / eps - 1) * sigm_alpha2 + m;
	sitk::Image&& sigmoid2 = 1 - sitk::Cast(sitk::Sigmoid(gradient_magnitude, sigm_alpha2, sigm_beta2, 1., 0.), sitk::sitkFloat64);
	sitk::Image&& contour_dist = sitk::Cast(sitk::SignedMaurerDistanceMap(initial_mask, false, false), sitk::sitkFloat64);
	sitk::Image&& contour = sitk::GeodesicActiveContourLevelSet(fast_marching_output, sigmoid2, 0.02, propagation_scaling);
	mPhaseField = sitk::Cast(sitk::BinaryThreshold(contour, -1000, 0, 0, 1), sitk::sitkFloat64);
}*/

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