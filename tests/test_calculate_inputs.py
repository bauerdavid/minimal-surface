"""Tests for the arguments of MinimalSurfaceCalculator.calculate,
calc_eikonal_and_transport_init and init_transport_slice.

Points and arrays use array index order, (z, y, x); the SimpleITK (x, y, z)
order is internal. A point outside the volume used to be an out-of-bounds read
in C++ that killed the host process, and only float64 arrays were accepted.

The shape is non-cubic, so that mixing up the axes shows. The image is a ball
with a soft edge; the initial plane callback thresholds the slice at half its
maximum. SYNTHETIC_RESULTS were recorded with the module as it was before the
input checks were added; float64 input must keep giving them.
"""
import numpy as np
import pytest

minimal_surface = pytest.importorskip("minimal_surface")

SHAPE = (30, 40, 44)
CENTER = (15, 20, 22)
RADIUS = 10.0
MAX_ITERATIONS = 10000

# (point1, point2): along x through the centre, and along z through the centre.
# The two give different results, which is what pins the (z, y, x) point order.
ALONG_X = ((15., 20., 12.), (15., 20., 32.))
ALONG_Z = ((5., 20., 22.), (25., 20., 22.))

# (points, use_meeting_points) -> number of voxels where the result is >= 0, sum, min, max
SYNTHETIC_RESULTS = {
    (ALONG_X, False): (27693, -139393.15450160307, -8.062257766723633, 8.485280990600586),
    (ALONG_X, True): (1951, -676864.082667103, -15.0, 8.9442720413208),
    (ALONG_Z, False): (17033, -151693.1586778948, -6.0, 12.041594505310059),
    (ALONG_Z, True): (1906, -877638.6230107842, -20.51828384399414, 8.059475451519965),
}


def make_volumes():
    zz, yy, xx = np.mgrid[0:SHAPE[0], 0:SHAPE[1], 0:SHAPE[2]]
    dist = np.sqrt((zz - CENTER[0]) ** 2 + (yy - CENTER[1]) ** 2 + (xx - CENTER[2]) ** 2)
    image = 1.0 / (1.0 + np.exp((dist - RADIUS) / 1.2))
    grad = np.sqrt(sum(g ** 2 for g in np.gradient(image)))
    phi = 0.2 + 0.8 * np.exp(-5.0 * grad / grad.max() * 5)
    return image, phi / phi.max()


def threshold_slice(image_slice, distance_slice):
    s = np.squeeze(np.asarray(image_slice))
    return (s >= 0.5 * s.max()).astype(float)


def make_calculator(use_meeting_points):
    calculator = minimal_surface.MinimalSurfaceCalculator()
    calculator.set_initial_plane_calculator(threshold_slice)
    calculator.set_using_meeting_points(use_meeting_points)
    return calculator


def calculate(phi, image, points, use_meeting_points=False, max_iterations=MAX_ITERATIONS):
    calculator = make_calculator(use_meeting_points)
    point1, point2 = (np.array(p) for p in points)
    return calculator.calculate(phi, image, point1, point2, True, max_iterations)


@pytest.fixture(scope="module")
def volumes():
    return make_volumes()


# SYNTHETIC_RESULTS were recorded on Windows (MSVC). On Linux (GCC) the results with meeting
# points are the same to about 1e-10. The results without meeting points differ slightly: they
# resample the volume, and the zero level set that the voxel count depends on moves by a few
# voxels with the compiler's floating-point behaviour (-2 of 27693 and -56 of 17033 voxels, sum
# 5e-6 away at most). The tolerances accept that and still catch a real change of the result.
COUNT_RTOL = 1e-2
SUM_RTOL = 1e-4
EXTREMA_RTOL = 1e-6


def assert_matches_recorded(out, key):
    inside, total, low, high = SYNTHETIC_RESULTS[key]
    assert out.dtype == np.float64 and out.shape == SHAPE
    assert int((out >= 0).sum()) == pytest.approx(inside, rel=COUNT_RTOL)
    assert out.sum() == pytest.approx(total, rel=SUM_RTOL)
    assert out.min() == pytest.approx(low, rel=EXTREMA_RTOL)
    assert out.max() == pytest.approx(high, rel=EXTREMA_RTOL)


# --------------------------------------------------------------------------
# results do not change
# --------------------------------------------------------------------------

@pytest.mark.parametrize("points, meeting", list(SYNTHETIC_RESULTS))
def test_float64_results_are_unchanged(volumes, points, meeting):
    image, phi = volumes
    out = calculate(phi, image, points, meeting)

    assert_matches_recorded(out, (points, meeting))


# --------------------------------------------------------------------------
# any real numeric dtype and memory layout
# --------------------------------------------------------------------------

def test_float32_image_works_and_matches_float64_of_the_same_data(volumes):
    image, phi = volumes
    image32, phi32 = image.astype(np.float32), phi.astype(np.float32)
    expected = calculate(phi32.astype(np.float64), image32.astype(np.float64), ALONG_X)

    out = calculate(phi32, image32, ALONG_X)

    assert out.dtype == np.float64
    np.testing.assert_array_equal(out, expected)


def test_uint16_image_gives_the_same_result_as_float64_of_the_same_data(volumes):
    image, phi = volumes
    image16 = np.round(image * 1000).astype(np.uint16)
    phi16 = np.round(phi * 1000).astype(np.uint16)
    expected = calculate(phi16.astype(np.float64), image16.astype(np.float64), ALONG_X)

    out = calculate(phi16, image16, ALONG_X)

    np.testing.assert_array_equal(out, expected)
    assert int((out >= 0).sum()) > 0


@pytest.mark.parametrize("dtype", [np.float16, np.float32, np.int16, np.int32, np.int64,
                                   np.uint8, np.uint32, np.bool_])
def test_every_real_dtype_gives_the_float64_result_of_the_same_data(volumes, dtype):
    image, phi = volumes
    image_t = (image > 0.5).astype(dtype) if dtype is np.bool_ else np.round(image * 100).astype(dtype)
    phi_t = phi.astype(dtype) if dtype in (np.float16, np.float32) else np.round(phi * 100).astype(dtype)
    # keep phi positive for the integer types
    if dtype not in (np.float16, np.float32):
        phi_t = np.maximum(phi_t, 1)
    expected = calculate(phi_t.astype(np.float64), image_t.astype(np.float64), ALONG_X)

    out = calculate(phi_t, image_t, ALONG_X)

    np.testing.assert_array_equal(out, expected)


def test_fortran_ordered_and_strided_input_give_the_same_result(volumes):
    image, phi = volumes
    expected = calculate(phi, image, ALONG_X)

    fortran_image, fortran_phi = np.asfortranarray(image), np.asfortranarray(phi)
    assert not fortran_image.flags.c_contiguous
    np.testing.assert_array_equal(calculate(fortran_phi, fortran_image, ALONG_X), expected)

    padded = np.zeros(SHAPE[:2] + (SHAPE[2] * 2,))
    padded[..., ::2] = image
    strided = padded[..., ::2]
    assert not strided.flags.c_contiguous
    np.testing.assert_array_equal(calculate(phi, strided, ALONG_X), expected)

    transposed = np.ascontiguousarray(image.transpose(2, 1, 0)).transpose(2, 1, 0)   # same data, F layout
    np.testing.assert_array_equal(calculate(phi, transposed, ALONG_X), expected)


def test_inputs_are_not_modified(volumes):
    image, phi = volumes
    image_before, phi_before = image.copy(), phi.copy()
    calculate(phi, image, ALONG_X)
    calculate(phi.astype(np.float32), image.astype(np.float32), ALONG_X)
    np.testing.assert_array_equal(image, image_before)
    np.testing.assert_array_equal(phi, phi_before)


@pytest.mark.parametrize("points", [ALONG_X, ALONG_Z])
def test_init_transport_slice_accepts_float32_and_gives_the_float64_slice(volumes, points):
    _, phi = volumes
    phi32 = phi.astype(np.float32)
    reference = make_calculator(False)
    reference.init_transport_slice(phi32.astype(np.float64), np.array(points[0]), np.array(points[1]))
    candidate = make_calculator(False)
    candidate.init_transport_slice(phi32, np.array(points[0]), np.array(points[1]))

    np.testing.assert_array_equal(candidate.get_init_plane(), reference.get_init_plane())


def test_init_transport_slice_handles_a_fortran_ordered_array(volumes):
    _, phi = volumes
    reference = make_calculator(False)
    reference.init_transport_slice(phi, np.array(ALONG_X[0]), np.array(ALONG_X[1]))
    candidate = make_calculator(False)
    candidate.init_transport_slice(np.asfortranarray(phi), np.array(ALONG_X[0]), np.array(ALONG_X[1]))

    np.testing.assert_array_equal(candidate.get_init_plane(), reference.get_init_plane())


def test_calc_eikonal_and_transport_init_accepts_float32_and_image_none(volumes):
    image, phi = volumes
    image32, phi32 = image.astype(np.float32), phi.astype(np.float32)
    p1, p2 = np.array(ALONG_X[0]), np.array(ALONG_X[1])

    reference = make_calculator(True)
    reference.calc_eikonal_and_transport_init(
        phi32.astype(np.float64), image32.astype(np.float64), p1, p2, True)
    candidate = make_calculator(True)
    candidate.calc_eikonal_and_transport_init(phi32, image32, p1, p2, True)
    np.testing.assert_array_equal(candidate.get_init_plane(), reference.get_init_plane())
    np.testing.assert_array_equal(candidate.get_distance_map(), reference.get_distance_map())

    no_image = make_calculator(True)
    no_image.calc_eikonal_and_transport_init(phi32, None, p1, p2, True)
    same_image = make_calculator(True)
    same_image.calc_eikonal_and_transport_init(phi32.astype(np.float64), phi32.astype(np.float64), p1, p2, True)
    np.testing.assert_array_equal(no_image.get_init_plane(), same_image.get_init_plane())


def test_points_may_be_any_real_sequence_or_dtype(volumes):
    image, phi = volumes
    expected = calculate(phi, image, ALONG_X)
    calculator = make_calculator(False)

    out = calculator.calculate(phi, image, [15, 20, 12], (15, 20, 32), True, MAX_ITERATIONS)
    np.testing.assert_array_equal(out, expected)
    out = make_calculator(False).calculate(
        phi, image, np.array([15, 20, 12], dtype=np.int16), np.array([15, 20, 32], dtype=np.float32),
        True, MAX_ITERATIONS)
    np.testing.assert_array_equal(out, expected)


# --------------------------------------------------------------------------
# points outside the volume, or not a point
# --------------------------------------------------------------------------

def all_entry_points(volumes):
    """Callables (point1, point2) -> None, one per method that takes points."""
    image, phi = volumes

    def calculate_(p1, p2):
        make_calculator(False).calculate(phi, image, p1, p2, True, 50)

    def calculate_meeting(p1, p2):
        make_calculator(True).calculate(phi, image, p1, p2, True, 50)

    def eikonal_init(p1, p2):
        make_calculator(True).calc_eikonal_and_transport_init(phi, image, p1, p2, True)

    def transport_slice(p1, p2):
        make_calculator(False).init_transport_slice(phi, p1, p2)

    return [calculate_, calculate_meeting, eikonal_init, transport_slice]


OUTSIDE_POINTS = [
    (SHAPE[0], 20., 22.), (15., SHAPE[1], 22.), (15., 20., SHAPE[2]),            # one past the end
    (-1., 20., 22.), (15., -1., 22.), (15., 20., -1.),                           # before the start
    (15., 20., 320.),                                                            # far outside
    (40., 20., 22.),             # fits the x axis (44) but not z (30): (x, y, z) order would accept it
    (SHAPE[0] - 0.5, 20., 22.), (15., 20., SHAPE[2] - 0.5),                      # the voxel edge itself
    (np.nan, 20., 22.), (15., np.inf, 22.), (15., 20., -np.inf),
]


@pytest.mark.parametrize("bad_point", OUTSIDE_POINTS, ids=repr)
@pytest.mark.parametrize("which", [0, 1])
def test_a_point_outside_the_volume_is_rejected_by_every_entry_point(volumes, bad_point, which):
    good_point = (15., 20., 22.)
    for entry in all_entry_points(volumes):
        points = [good_point, good_point]
        points[which] = bad_point
        with pytest.raises(ValueError) as excinfo:
            entry(np.array(points[0]), np.array(points[1]))
        message = str(excinfo.value)
        assert "point%d" % (which + 1) in message
        assert str(SHAPE) in message                       # names the shape ...
        assert "outside" in message
        if np.all(np.isfinite(bad_point)):                 # ... and the point
            assert str(list(map(float, bad_point))) in message


def test_the_message_names_the_axis_in_array_order(volumes):
    image, phi = volumes
    with pytest.raises(ValueError, match=r"coordinate 0 \(axis 0"):
        make_calculator(False).calculate(phi, image, np.array([40., 20., 22.]), np.array(ALONG_X[1]), True, 50)
    with pytest.raises(ValueError, match=r"coordinate 2 \(axis 2"):
        make_calculator(False).calculate(phi, image, np.array(ALONG_X[0]), np.array([15., 20., 44.]), True, 50)


@pytest.mark.parametrize("bad_point", [(1., 2.), (1., 2., 3., 4.), (), [[15., 20., 22.]], 5.,
                                       np.zeros((3, 1))], ids=repr)
def test_a_point_that_is_not_three_values_is_rejected(volumes, bad_point):
    for entry in all_entry_points(volumes):
        with pytest.raises(ValueError) as excinfo:
            entry(bad_point, np.array(ALONG_X[1]))
        assert str(SHAPE) in str(excinfo.value)


@pytest.mark.parametrize("bad_point", [np.array([1j, 2, 3]), ["a", "b", "c"], None], ids=repr)
def test_a_point_that_is_not_real_numbers_is_rejected(volumes, bad_point):
    for entry in all_entry_points(volumes):
        with pytest.raises(TypeError):
            entry(bad_point, np.array(ALONG_X[1]))


@pytest.mark.parametrize("points", [
    ((0., 0., 0.), (29., 39., 43.)),               # opposite corners
    ((15., 20., 0.), (15., 20., 43.)),             # the two ends of the x axis
    ((0., 20., 22.), (29., 20., 22.)),             # the two ends of the z axis
    ((15., 20., -0.4), (15., 20., 43.4)),          # inside the end voxels, off their centres
    ((15., 20., 12.6), (15.25, 20.5, 31.75)),
], ids=repr)
@pytest.mark.parametrize("meeting", [False, True])
def test_points_on_the_border_are_accepted_and_do_not_crash(volumes, points, meeting):
    image, phi = volumes
    out = calculate(phi, image, points, meeting, max_iterations=200)
    assert out.shape == SHAPE and np.isfinite(out).all()


def test_a_rejected_call_leaves_the_calculator_usable(volumes):
    image, phi = volumes
    calculator = make_calculator(False)
    with pytest.raises(ValueError):
        calculator.calculate(phi, image, np.array([15., 20., 12.]), np.array([15., 20., 320.]), True, 50)
    out = calculator.calculate(phi, image, np.array(ALONG_X[0]), np.array(ALONG_X[1]), True, MAX_ITERATIONS)

    assert_matches_recorded(out, (ALONG_X, False))


# --------------------------------------------------------------------------
# input that is not a real numeric 3D array
# --------------------------------------------------------------------------

def test_non_numeric_arrays_are_rejected_with_a_type_error(volumes):
    image, phi = volumes
    p1, p2 = np.array(ALONG_X[0]), np.array(ALONG_X[1])
    calculator = make_calculator(False)
    for bad in (image.astype(np.complex128), image.astype(str), image.astype(object), None):
        with pytest.raises(TypeError, match="image"):
            calculator.calculate(phi, bad, p1, p2, True, 50)
        with pytest.raises(TypeError, match="phi"):
            calculator.calculate(bad, image, p1, p2, True, 50)
        with pytest.raises(TypeError):
            calculator.init_transport_slice(bad, p1, p2)
        with pytest.raises(TypeError):
            calculator.calc_eikonal_and_transport_init(bad, image, p1, p2, True)


def test_arrays_of_the_wrong_dimensionality_or_a_mismatched_shape_are_rejected(volumes):
    image, phi = volumes
    p1, p2 = np.array(ALONG_X[0]), np.array(ALONG_X[1])
    calculator = make_calculator(False)
    for bad in (image[0], image[..., None, None], np.float64(1.0)):
        with pytest.raises(ValueError, match="3D"):
            calculator.calculate(phi, bad, p1, p2, True, 50)
    with pytest.raises(ValueError, match="empty"):
        calculator.calculate(phi[:0], image[:0], p1, p2, True, 50)
    with pytest.raises(ValueError, match="same shape"):
        calculator.calculate(phi, image[:-1], p1, p2, True, 50)
    with pytest.raises(ValueError, match="same shape"):
        calculator.calc_eikonal_and_transport_init(phi, image[:, :-1], p1, p2, True)
