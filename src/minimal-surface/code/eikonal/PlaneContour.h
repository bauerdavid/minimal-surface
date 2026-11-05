#pragma once
#include <sitkImage.h>
#include <vector>

namespace sitk = itk::simple;

// The initial contour of the transport function, from the distance values on the meeting plane.
//
// `distance` is the sum of the distances of the two points on the plane: the two fronts have met
// there, so it is the length of the cheapest path from one point through the pixel to the other.
// Along the surface of the object that path is cheap, so the surface is a valley of this map, with
// the (more expensive) inside of the object enclosed by it and a long monotone ascent outside.
// The valley is a closed curve around (`seed_x`, `seed_y`), the pixel where the line between the
// two points crosses the plane, which lies inside of the object.
//
// Texture inside of the object makes valleys of its own, so the curve is not the lowest value on
// each ray from the seed. It is found with dynamic programming over all of the rays at once, which
// needs a closed, smooth curve and prefers the outer one, and is then moved to the foot of the
// inner flank of the valley, which is where the edge of the object is.
//
// `chord` is the distance of the two points. Returns an 8 bit image of the size of `distance`, or
// an empty (default constructed) image if there is nothing to find: the distances are constant, or
// the curve is much smaller than the object has to be (the seed lies on the surface and the curve
// collapsed onto it). The caller then needs another way of getting the contour.
//
// The contour is at the edge of the object, at the foot of the inner flank of the valley, or, with
// `at_floor`, at the bottom of the valley. The propagation of the transport function settles at the
// bottom of the valley whatever the contour in the plane was (a contour 4 pixels off is forgotten
// within two slices), so a contour at the edge leaves a notch in the surface at the plane: the
// propagated surface, further away, is outside of it. A contour at the bottom does not.
sitk::Image ContourFromDistanceSlice(const sitk::Image& distance, double seed_x, double seed_y, double chord, bool at_floor = false);

// The radius, on each of PlaneContourRayCount() rays from the seed, of the bottom of the valley and
// of the edge of the object, as ContourFromDistanceSlice finds them. False if there is no valley.
struct PlaneValley {
	std::vector<double> floor_radius;
	std::vector<double> edge_radius;
};
bool FindValley(const sitk::Image& distance, double seed_x, double seed_y, double chord, PlaneValley& valley);

int PlaneContourRayCount();
// The value of a quantity given on each ray, in the direction (dx, dy) from the seed.
double PlaneContourRayValue(const std::vector<double>& per_ray, double dx, double dy);
// The contour with the given radius on each ray; empty if it is degenerate (see above).
sitk::Image ContourFromRadii(const sitk::Image& like, double seed_x, double seed_y, const std::vector<double>& radii, double chord);
// The radius at which each ray from the seed leaves `contour`; empty if the seed is not inside it.
std::vector<double> ContourRadii(const sitk::Image& contour, double seed_x, double seed_y);
// The first radius on each ray from the seed at which `values` (w x h, row by row, bilinear) falls below
// `level`; -1 for a ray that leaves the image before; empty if the seed is below the level.
std::vector<double> FieldRadii(const std::vector<double>& values, int w, int h, double seed_x, double seed_y, double level);
// a - b on each ray, smoothed over the rays like the radii are.
std::vector<double> SmoothRayValues(const std::vector<double>& a, const std::vector<double>& b);
