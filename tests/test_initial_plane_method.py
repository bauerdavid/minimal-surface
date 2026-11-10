"""Tests for how the initial contour in the meeting plane is found.

By default it comes from the distances of the two points on the plane: the surface of the object is a
valley of their sum, and the contour follows it. That needs no phase field and no function from the
user. set_initial_plane_method chooses the phase field instead; a function given to
set_initial_plane_calculator wins over both.

The volume is a ball with a soft edge. The speed function is made like the one in the README: slow
inside of the object and in flat background, fast along the edge. The distance method relies on that:
if the surface were no cheaper to travel along than the inside of the object, the surface would not
be a valley of the distances and there would be nothing to follow.
"""
import numpy as np
import pytest

minimal_surface = pytest.importorskip("minimal_surface")

SHAPE = (30, 40, 44)
CENTER = (15, 20, 22)
RADIUS = 10.0
POINTS = ((15., 20., 12.), (15., 20., 32.))


def make_volumes():
    zz, yy, xx = np.mgrid[0:SHAPE[0], 0:SHAPE[1], 0:SHAPE[2]]
    dist = np.sqrt((zz - CENTER[0]) ** 2 + (yy - CENTER[1]) ** 2 + (xx - CENTER[2]) ** 2)
    image = 1.0 / (1.0 + np.exp((dist - RADIUS) / 1.2))
    delta = 0.1
    gradient = np.sqrt(sum(g ** 2 for g in np.gradient(image)))
    features = gradient * (1 + (delta - 1) * image ** 2) + delta * (1 - image ** 2)
    features = (features - features.min()) / (features.max() - features.min())
    alpha = 3e-3
    phi = alpha + (1 - alpha) * np.exp(-5 * features)
    return image, phi / phi.max(), dist <= RADIUS


def threshold_slice(image_slice, distance_slice):
    s = np.squeeze(np.asarray(image_slice))
    return (s >= 0.5 * s.max()).astype(float)


def iou(a, b):
    return np.logical_and(a, b).sum() / np.logical_or(a, b).sum()


def init_plane(method=None, callback=None):
    image, phi, _ = make_volumes()
    calculator = minimal_surface.MinimalSurfaceCalculator()
    if method is not None:
        calculator.set_initial_plane_method(method)
    if callback is not None:
        calculator.set_initial_plane_calculator(callback)
    calculator.calc_eikonal_and_transport_init(phi, image, np.array(POINTS[0]), np.array(POINTS[1]), True)
    return calculator, calculator.get_init_plane()


def test_default_contour_is_the_object_in_the_plane():
    # The ball is cut by the plane through its centre in a disc of RADIUS.
    _, plane = init_plane()
    assert plane.sum() > 0
    area = plane.sum()
    assert abs(area - np.pi * RADIUS ** 2) < 0.15 * np.pi * RADIUS ** 2
    _, reference = init_plane(callback=threshold_slice)
    assert iou(plane > 0, reference > 0) > 0.85


def test_default_is_deterministic():
    _, first = init_plane()
    _, second = init_plane()
    np.testing.assert_array_equal(first, second)


def test_default_is_the_distance_method():
    _, default = init_plane()
    _, explicit = init_plane("distance")
    np.testing.assert_array_equal(default, explicit)


def test_phase_field_method_gives_a_contour_too():
    _, plane = init_plane("phase_field")
    assert abs(plane.sum() - np.pi * RADIUS ** 2) < 0.3 * np.pi * RADIUS ** 2


def test_user_function_wins_over_the_method():
    _, with_function = init_plane("distance", callback=threshold_slice)
    _, only_function = init_plane(callback=threshold_slice)
    np.testing.assert_array_equal(with_function, only_function)


def test_method_name_is_checked():
    calculator = minimal_surface.MinimalSurfaceCalculator()
    with pytest.raises(ValueError):
        calculator.set_initial_plane_method("circle")


def test_the_result_with_the_default_method_is_the_ball():
    image, phi, ball = make_volumes()
    calculator = minimal_surface.MinimalSurfaceCalculator()
    result = calculator.calculate(phi, image, np.array(POINTS[0]), np.array(POINTS[1]))
    assert iou(result >= 0, ball) > 0.8


def test_rotated_distance_maps_and_plane_sampling():
    calculator, _ = init_plane()
    d1, d2, plane_slice = calculator.get_rotated_distance_maps()
    assert d1.shape == d2.shape
    assert 0 <= plane_slice < d1.shape[2]
    p1, p2 = calculator.get_rotated_points()
    # the plane lies between the two points
    assert min(p1[2], p2[2]) < plane_slice < max(p1[2], p2[2])
    image, _, ball = make_volumes()
    sampled = calculator.sample_init_plane(image)
    assert sampled.size == d1.shape[0] * d1.shape[1]
    # through the centre of the ball: the middle of the plane is inside, the rest of it mostly out
    assert sampled.max() > 0.9


# ---- the notch at the meeting plane -------------------------------------------------------------
# The propagation settles at the bottom of the valley of the distances whatever the contour in the
# plane was and forgets the contour within two slices. A contour that is not on that surface left a
# notch in the result at the plane. Here the plane is x = CENTER[2], perpendicular to the first
# axis of the volume's last axis, and the ball is a sphere, so the radius of a section of the result
# at an offset k from the plane is known.

def section_radius(result, offset):
    return np.sqrt((result[:, :, CENTER[2] + offset] >= 0).sum() / np.pi)


def disc_of_radius(radius):
    """A function for set_initial_plane_calculator: the disc in the middle of the object of the plane."""
    def disc(image_slice, distance_slice):
        s = np.squeeze(np.asarray(image_slice))
        inside = s >= 0.5 * s.max()
        rows, cols = np.nonzero(inside)
        yy, xx = np.mgrid[0:s.shape[0], 0:s.shape[1]]
        return ((yy - rows.mean()) ** 2 + (xx - cols.mean()) ** 2 <= radius ** 2).astype(float)
    return disc


def run(callback=None, align=False):
    image, phi, _ = make_volumes()
    calculator = minimal_surface.MinimalSurfaceCalculator()
    if callback is not None:
        calculator.set_initial_plane_calculator(callback)
    calculator.set_align_surface(align)
    return calculator, calculator.calculate(phi, image, np.array(POINTS[0]), np.array(POINTS[1]))


def test_default_has_no_notch_at_the_plane():
    _, result = run()
    around = np.mean([section_radius(result, k) for k in (-2, 2)])
    assert abs(section_radius(result, 0) - around) < 0.5
    # and it is the ball
    for k, radius in ((0, 10.0), (2, np.sqrt(100 - 4)), (4, np.sqrt(100 - 16))):
        assert abs(section_radius(result, k) - radius) < 0.6


def test_a_contour_that_is_too_small_leaves_a_notch_unless_aligned():
    small = disc_of_radius(8.0)
    _, plain = run(small)
    assert section_radius(plain, 3) - section_radius(plain, 0) > 0.5
    calculator, aligned = run(small, align=True)
    # the surface passes through the given contour in the plane (a contour of pixels is only known to
    # about a pixel), and is smooth through it
    assert abs(section_radius(aligned, 0) - 8.0) < 1.0
    around = np.mean([section_radius(aligned, k) for k in (-2, 2)])
    assert abs(section_radius(aligned, 0) - around) < 0.5
    assert len(calculator.get_surface_offsets()) > 0


def test_align_leaves_a_right_contour_alone():
    right = disc_of_radius(10.0)
    _, plain = run(right)
    calculator, aligned = run(right, align=True)
    for k in (-4, -2, 0, 2, 4):
        assert abs(section_radius(aligned, k) - section_radius(plain, k)) < 0.6
    assert abs(np.mean(calculator.get_surface_offsets())) < 0.6


def test_no_offsets_without_align():
    calculator, _ = run(disc_of_radius(8.0))
    assert len(calculator.get_surface_offsets()) == 0
