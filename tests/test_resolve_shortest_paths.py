"""Tests for MinimalSurfaceCalculator.resolve_shortest_paths.

The function follows the steepest descent of a distance map from a start point.
Points and the returned path use array index order, (z, y, x). These tests use
the distance from a known source voxel on a non-cubic array, so that mixing up
the axis order, or reading the array with the wrong layout, changes the result.
"""
import numpy as np
import pytest

minimal_surface = pytest.importorskip("minimal_surface")


def distance_from(shape, source):
    zz, yy, xx = np.mgrid[0:shape[0], 0:shape[1], 0:shape[2]]
    return np.sqrt((zz - source[0]) ** 2 + (yy - source[1]) ** 2 + (xx - source[2]) ** 2)


@pytest.fixture
def calculator():
    return minimal_surface.MinimalSurfaceCalculator()


SHAPE = (8, 10, 14)          # all different, so any axis mix-up shows
SOURCE = (2, 3, 4)
START = (6, 8, 12)


def resolve(calculator, point, data):
    return calculator.resolve_shortest_paths(np.asarray(point, dtype=np.int64), data)


def test_path_runs_from_the_point_to_the_source(calculator):
    data = distance_from(SHAPE, SOURCE).astype(np.float64)
    path = resolve(calculator, START, data)

    assert path.dtype == np.int64
    assert path.ndim == 2 and path.shape[1] == 3
    assert tuple(path[0]) == START
    assert tuple(path[-1]) == SOURCE


def test_path_stays_inside_the_array_and_moves_one_voxel_at_a_time(calculator):
    data = distance_from(SHAPE, SOURCE).astype(np.float64)
    path = resolve(calculator, START, data)

    assert (path >= 0).all()
    assert (path < np.asarray(SHAPE)).all()
    steps = np.abs(np.diff(path, axis=0))
    assert steps.max() <= 1
    assert steps.sum(axis=1).min() >= 1          # every step actually moves


def test_path_descends_the_distance_map(calculator):
    data = distance_from(SHAPE, SOURCE).astype(np.float64)
    path = resolve(calculator, START, data)

    values = data[path[:, 0], path[:, 1], path[:, 2]]
    assert (np.diff(values) < 0).all()


@pytest.mark.parametrize("start, source", [
    ((6, 8, 12), (2, 3, 4)),
    ((1, 9, 2), (5, 2, 11)),
    ((7, 1, 13), (3, 7, 5)),
])
def test_axis_order_is_zyx_for_input_and_output(calculator, start, source):
    data = distance_from(SHAPE, source).astype(np.float64)
    path = resolve(calculator, start, data)

    assert tuple(path[0]) == start
    assert tuple(path[-1]) == source


def test_non_contiguous_input_gives_the_same_path(calculator):
    data = distance_from(SHAPE, SOURCE).astype(np.float64)
    expected = resolve(calculator, START, data)

    fortran = np.asfortranarray(data)
    assert not fortran.flags.c_contiguous
    np.testing.assert_array_equal(resolve(calculator, START, fortran), expected)

    padded = np.zeros((SHAPE[0], SHAPE[1], SHAPE[2] * 2))
    padded[..., ::2] = data
    strided = padded[..., ::2]
    assert not strided.flags.c_contiguous
    np.testing.assert_array_equal(resolve(calculator, START, strided), expected)


def test_a_point_that_is_not_three_values_is_rejected(calculator):
    data = distance_from(SHAPE, SOURCE).astype(np.float64)
    with pytest.raises(ValueError):
        resolve(calculator, (1, 2), data)
    with pytest.raises(ValueError):
        resolve(calculator, (1, 2, 3, 4), data)


@pytest.mark.parametrize("point", [
    (SHAPE[0], 0, 0), (0, SHAPE[1], 0), (0, 0, SHAPE[2]),
    (-1, 0, 0), (0, -1, 0), (0, 0, -1),
])
def test_a_point_outside_the_array_is_rejected(calculator, point):
    data = distance_from(SHAPE, SOURCE).astype(np.float64)
    with pytest.raises(ValueError):
        resolve(calculator, point, data)


@pytest.mark.parametrize("shape, start, source", [
    # Each of these used to walk out of the volume: the neighbour scan was
    # centred on a position clamped away from the border, while the step it
    # chose was applied to the real position.
    ((14, 14, 5), (13, 12, 4), (13, 13, 4)),
    ((5, 9, 12), (0, 0, 0), (4, 8, 11)),
    ((12, 5, 9), (11, 4, 8), (0, 0, 0)),
    ((9, 12, 5), (0, 11, 0), (8, 0, 4)),
])
def test_paths_that_start_or_end_at_the_border_stay_inside(calculator, shape, start, source):
    data = distance_from(shape, source).astype(np.float64)
    path = resolve(calculator, start, data)

    assert (path >= 0).all()
    assert (path < np.asarray(shape)).all()
    assert tuple(path[0]) == start
    assert tuple(path[-1]) == source


def test_random_paths_stay_inside_the_array_and_reach_the_source(calculator):
    rng = np.random.default_rng(7)
    checked = 0
    for trial in range(400):
        shape = tuple(int(s) for s in rng.integers(5, 18, 3))
        edge = trial % 2 == 0           # half of the cases hug the border
        def pick():
            return tuple(int(rng.choice([0, s - 1, rng.integers(0, s)])) if edge
                         else int(rng.integers(0, s)) for s in shape)
        start, source = pick(), pick()
        if start == source:
            continue
        data = distance_from(shape, source).astype(np.float64)
        path = resolve(calculator, start, data)
        assert (path >= 0).all() and (path < np.asarray(shape)).all(), (shape, start, source)
        assert tuple(path[-1]) == source, (shape, start, source)
        checked += 1
    assert checked > 300
