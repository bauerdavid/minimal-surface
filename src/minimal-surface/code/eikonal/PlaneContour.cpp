#include "PlaneContour.h"
#include "SimpleITK.h"
#include <algorithm>
#include <cmath>
#include <vector>

using namespace std;

namespace {
const double kPi = 3.14159265358979323846;
// Number of rays from the seed. More rays do not change the result visibly.
const int kRays = 96;
// Longest ray, in pixels. The rays stop at the border of the image anyway.
const int kMaxRadius = 200;
// The curve prefers a larger radius, up to this many half chords: the same dark valley in the
// texture inside of the object would otherwise win over the surface.
const double kOuterPreference = 0.2;
const double kOuterPreferenceLimit = 2.0;
// Price of the curve changing its radius by a pixel from one ray to the next.
const double kStepPenalty = 0.02;
// The valley's foot is looked for this far (pixels) around the curve ...
const int kValleyWindow = 4;
// ... and the edge of the object this far (pixels) inside of the valley's bottom, at this fraction
// of the height of the inner flank.
const int kFlankWidth = 6;
const double kFlankLevel = 0.05;
// Smoothing of the radius over the rays, in rays.
const double kRaySigma = 1.5;
// A contour with the area of a circle smaller than this fraction of the chord is the collapsed
// curve of a seed on the surface, not the object.
const double kMinRadiusPerChord = 0.1;
// Cost of the samples outside of the image.
const double kOutsideCost = 5.0;

double SampleBilinear(const vector<double>& image, int w, int h, double x, double y) {
	x = min(max(x, 0.0), w - 1.0);
	y = min(max(y, 0.0), h - 1.0);
	int x0 = (int)floor(x), y0 = (int)floor(y);
	int x1 = min(x0 + 1, w - 1), y1 = min(y0 + 1, h - 1);
	double fx = x - x0, fy = y - y0;
	double top = image[(size_t)y0 * w + x0] * (1 - fx) + image[(size_t)y0 * w + x1] * fx;
	double bottom = image[(size_t)y1 * w + x0] * (1 - fx) + image[(size_t)y1 * w + x1] * fx;
	return top * (1 - fy) + bottom * fy;
}

// The radius of each ray of the closed path with the least cost, where cost[ray][radius] is the
// cost of the path being at that radius on that ray, and the radius changes by at most one from a
// ray to the next, at kStepPenalty a pixel. A closed path is found by going around three times
// and keeping the middle one.
vector<int> ClosedPath(const vector<vector<double>>& cost) {
	const int rays = (int)cost.size();
	const int radii = (int)cost[0].size();
	const int total = 3 * rays;
	vector<vector<double>> acc(total, vector<double>(radii));
	vector<vector<int>> from(total, vector<int>(radii, 0));
	acc[0] = cost[0];
	for (int i = 1; i < total; i++) {
		const vector<double>& prev = acc[i - 1];
		const vector<double>& here = cost[i % rays];
		for (int r = 0; r < radii; r++) {
			double best = 0;
			int best_from = -1;
			// the radius on the previous ray: r + 1, r, r - 1
			for (int d = -1; d <= 1; d++) {
				int p = r - d;
				if (p < 0 || p >= radii) continue;
				double c = prev[p] + kStepPenalty * abs(d);
				if (best_from < 0 || c < best) {
					best = c;
					best_from = p;
				}
			}
			acc[i][r] = best + here[r];
			from[i][r] = best_from;
		}
	}
	vector<int> path(total);
	path[total - 1] = (int)(min_element(acc[total - 1].begin(), acc[total - 1].end()) - acc[total - 1].begin());
	for (int i = total - 1; i > 0; i--)
		path[i - 1] = from[i][path[i]];
	return vector<int>(path.begin() + rays, path.begin() + 2 * rays);
}

// Gaussian smoothing of a closed sequence.
vector<double> SmoothClosed(const vector<double>& values, double sigma) {
	const int n = (int)values.size();
	const int radius = (int)(4.0 * sigma + 0.5);
	vector<double> kernel(2 * radius + 1);
	double sum = 0;
	for (int k = -radius; k <= radius; k++) {
		kernel[k + radius] = exp(-0.5 * k * k / (sigma * sigma));
		sum += kernel[k + radius];
	}
	vector<double> smoothed(n, 0.0);
	for (int i = 0; i < n; i++) {
		for (int k = -radius; k <= radius; k++)
			smoothed[i] += kernel[k + radius] * values[((i + k) % n + n) % n];
		smoothed[i] /= sum;
	}
	return smoothed;
}
}

int PlaneContourRayCount() {
	return kRays;
}

double PlaneContourRayValue(const vector<double>& per_ray, double dx, double dy) {
	double angle = atan2(dy, dx);
	if (angle < 0) angle += 2 * kPi;
	double position = angle / (2 * kPi) * kRays;
	int i0 = (int)floor(position) % kRays;
	int i1 = (i0 + 1) % kRays;
	double f = position - floor(position);
	return per_ray[i0] * (1 - f) + per_ray[i1] * f;
}

sitk::Image ContourFromRadii(const sitk::Image& like, double seed_x, double seed_y, const vector<double>& radii, double chord) {
	vector<unsigned> size = like.GetSize();
	const int w = size[0], h = size[1];
	sitk::Image contour(size, sitk::sitkUInt8);
	unsigned char* out = contour.GetBufferAsUInt8();
	size_t area = 0;
	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			double dx = x - seed_x, dy = y - seed_y;
			bool inside = hypot(dx, dy) <= PlaneContourRayValue(radii, dx, dy);
			out[(size_t)y * w + x] = inside ? 1 : 0;
			if (inside) area++;
		}
	}
	if (area == 0 || sqrt(area / kPi) < kMinRadiusPerChord * chord)
		return sitk::Image();
	return contour;
}

bool FindValley(const sitk::Image& distance, double seed_x, double seed_y, double chord, PlaneValley& valley) {
	vector<unsigned> size = distance.GetSize();
	const int w = size[0], h = size[1];
	const size_t pixel_count = (size_t)w * h;
	if (w < 3 || h < 3 || !(chord > 1))
		return false;

	// Where either point has not been reached (outside of the data the maps are -1) there is no
	// distance, only a placeholder that would be the lowest value: use the highest instead.
	vector<double> values(distance.GetBufferAsDouble(), distance.GetBufferAsDouble() + pixel_count);
	double highest = 0;
	for (double v : values)
		highest = max(highest, v);
	if (!(highest > 0))
		return false;
	for (double& v : values)
		if (!(v > 0)) v = highest;
	sitk::Image smoothed_image(size, sitk::sitkFloat64);
	copy(values.begin(), values.end(), smoothed_image.GetBufferAsDouble());
	smoothed_image = sitk::DiscreteGaussian(smoothed_image, 1.0, 32, 0.01, false);
	const double* smoothed = smoothed_image.GetBufferAsDouble();
	vector<double> D(smoothed, smoothed + pixel_count);
	const double lowest_d = *min_element(D.begin(), D.end());
	const double range = *max_element(D.begin(), D.end()) - lowest_d;
	if (!(range > 0))
		return false;

	const int radii = min((int)hypot((double)w, (double)h), kMaxRadius);
	if (radii < 4)
		return false;
	vector<vector<double>> profile(kRays, vector<double>(radii));
	vector<vector<double>> cost(kRays, vector<double>(radii));
	for (int i = 0; i < kRays; i++) {
		double angle = 2 * kPi * i / kRays;
		for (int r = 0; r < radii; r++) {
			double x = seed_x + cos(angle) * r, y = seed_y + sin(angle) * r;
			profile[i][r] = SampleBilinear(D, w, h, x, y);
			bool outside = y < -0.5 || y > h - 0.5 || x < -0.5 || x > w - 0.5;
			cost[i][r] = outside ? kOutsideCost :
				(profile[i][r] - lowest_d) / range - kOuterPreference * min(r / (0.5 * chord), kOuterPreferenceLimit);
		}
	}
	vector<int> path = ClosedPath(cost);

	// The bottom of the valley near the path, and the edge of the object: the point of the inner
	// flank, between the bottom and the crest inside of it, where the distance has risen by
	// kFlankLevel of the way to the crest.
	vector<double> floor_radius(kRays), edge(kRays);
	for (int i = 0; i < kRays; i++) {
		const vector<double>& p = profile[i];
		int lo = max(1, path[i] - kValleyWindow), hi = min(radii - 1, path[i] + kValleyWindow + 1);
		int bottom = lo;
		for (int r = lo; r < hi; r++)
			if (p[r] < p[bottom]) bottom = r;
		floor_radius[i] = bottom;
		int from = max(1, bottom - kFlankWidth);
		int crest = from;
		for (int r = from; r <= bottom; r++)
			if (p[r] > p[crest]) crest = r;
		double level = p[bottom] + kFlankLevel * (p[crest] - p[bottom]);
		int r = crest;
		while (r <= bottom && p[r] > level) r++;
		edge[i] = r <= bottom ? r : crest;
	}
	valley.floor_radius = SmoothClosed(floor_radius, kRaySigma);
	valley.edge_radius = SmoothClosed(edge, kRaySigma);
	for (int i = 0; i < kRays; i++) {
		valley.floor_radius[i] = max(valley.floor_radius[i], 1.0);
		valley.edge_radius[i] = max(valley.edge_radius[i], 1.0);
	}
	return true;
}

sitk::Image ContourFromDistanceSlice(const sitk::Image& distance, double seed_x, double seed_y, double chord, bool at_floor) {
	PlaneValley valley;
	if (!FindValley(distance, seed_x, seed_y, chord, valley))
		return sitk::Image();
	return ContourFromRadii(distance, seed_x, seed_y, at_floor ? valley.floor_radius : valley.edge_radius, chord);
}

vector<double> FieldRadii(const vector<double>& values, int w, int h, double seed_x, double seed_y, double level) {
	// The first place along the ray where the field (bilinear) drops below the level.
	if (SampleBilinear(values, w, h, seed_x, seed_y) < level)
		return vector<double>();
	const double step = 0.1;
	const int steps = (int)(hypot((double)w, (double)h) / step);
	vector<double> radii(kRays, -1.0);
	for (int i = 0; i < kRays; i++) {
		double angle = 2 * kPi * i / kRays;
		double previous = SampleBilinear(values, w, h, seed_x, seed_y);
		for (int k = 1; k <= steps; k++) {
			double x = seed_x + cos(angle) * k * step, y = seed_y + sin(angle) * k * step;
			if (x < -0.5 || x > w - 0.5 || y < -0.5 || y > h - 0.5)
				break;
			double v = SampleBilinear(values, w, h, x, y);
			if (v < level) {
				double f = previous - v > 1e-12 ? (previous - level) / (previous - v) : 0.0;
				radii[i] = (k - 1 + f) * step;
				break;
			}
			previous = v;
		}
	}
	return radii;
}

vector<double> ContourRadii(const sitk::Image& contour, double seed_x, double seed_y) {
	vector<unsigned> size = contour.GetSize();
	const int w = size[0], h = size[1];
	// Only the object around the seed, without holes: a contour from a threshold or a segmentation
	// has holes and bits that are not connected to it, and a ray would leave at the first of them.
	sitk::Image filled = sitk::BinaryFillhole(sitk::Cast(contour, sitk::sitkUInt8));
	sitk::Image labels = sitk::ConnectedComponent(filled);
	const int seed_pixel_x = min(max((int)lround(seed_x), 0), w - 1), seed_pixel_y = min(max((int)lround(seed_y), 0), h - 1);
	const uint32_t* label_buffer = labels.GetBufferAsUInt32();
	const uint32_t seed_label = label_buffer[(size_t)seed_pixel_y * w + seed_pixel_x];
	if (seed_label == 0)
		return vector<double>();
	vector<double> values((size_t)w * h);
	for (size_t i = 0; i < values.size(); i++)
		values[i] = label_buffer[i] == seed_label ? 1.0 : 0.0;
	// The bilinear mask is a half where a pixel centre inside meets one outside, which is the
	// edge of the pixel. The border of the image is outside.
	vector<double> radii = FieldRadii(values, w, h, seed_x, seed_y, 0.5);
	for (double& r : radii)
		if (r < 0) r = hypot((double)w, (double)h);
	return radii;
}

vector<double> SmoothRayValues(const vector<double>& a, const vector<double>& b) {
	vector<double> difference(a.size());
	for (size_t i = 0; i < a.size(); i++)
		difference[i] = a[i] - b[i];
	return SmoothClosed(difference, kRaySigma);
}
