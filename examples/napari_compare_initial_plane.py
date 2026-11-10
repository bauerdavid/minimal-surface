"""Compare the ways of finding the initial contour in the meeting plane, on the cells3d nuclei, in napari.

Run it in an environment that has napari and minimal_surface (and scikit-image, scipy, pooch):

    python examples/napari_compare_initial_plane.py

It opens the nuclei channel (the crop of the README example) with two points on the surface of a nucleus.
Move the two points (select the layer, drag them) or add your own, then press "Calculate": the result of
each method is added as a labels layer, and the time it took is printed.

    distance     the contour comes from the distances of the two points on the meeting plane (the default)
    phase_field  the contour is found by evolving a phase field in the plane
    threshold    the contour is the image in the plane above a threshold, like `segment_slice` in the README
                 (the one that needs a person to look at the plane; here a fixed threshold)

With "align_surface" the result passes through the contour of the method in the plane (otherwise only the one
found from the distances is corrected, so that the surface has no notch at the plane), see the README.

The initial plane of each method is added as a hidden labels layer too ("... plane"): it is a 2D image in the
meeting plane, which is tilted, so it does not line up with the volume; show it on its own to see the contour.
"""
import time

import numpy as np
from scipy import ndimage

import minimal_surface

METHODS = ("distance", "phase_field", "threshold")
# README example: the part of the nuclei channel around one nucleus, and two points on its surface (z, y, x)
CROP = (slice(None), slice(60, 130), slice(130, 190))
POINTS = np.array([[35., 60., 20.], [35., 12., 39.]])
THRESHOLD = 0.24


def calc_features(image):
    delta = 0.1
    max_ = np.quantile(image, .95)
    min_ = image.min()
    image = np.clip((image - min_) / (max_ - min_), 0, 1)
    # considers both the image gradient and the intensity: high where the gradient is high and the
    # intensity is low
    gradient = ndimage.gaussian_gradient_magnitude(image, (1., 1., 1.))
    weights = 1 + (delta - 1) * image ** 2
    phi = gradient * weights + delta * (1 - image ** 2)
    return (phi - phi.min()) / (phi.max() - phi.min())


def load_data():
    from skimage.data import cells3d
    data = cells3d()[:, 1][CROP].astype(float)
    data = (data - data.min()) / (data.max() - data.min())
    data = ndimage.gaussian_filter(data, 1.)
    alpha = 3e-3
    phi = alpha + (1 - alpha) * np.exp(-5 * calc_features(data))
    return data, phi / phi.max()


def calculate(method, phi, data, p1, p2, align=False):
    """(transport function >= 0, the initial plane, seconds) for one method."""
    calculator = minimal_surface.MinimalSurfaceCalculator()
    if method == "threshold":
        calculator.set_initial_plane_calculator(lambda plane, distance: (np.squeeze(plane) >= THRESHOLD).astype(float))
    else:
        calculator.set_initial_plane_method(method)
    calculator.set_align_surface(align)
    start = time.time()
    transport = calculator.calculate(phi, data, p1, p2)
    seconds = time.time() - start
    return transport >= 0, np.squeeze(calculator.get_init_plane()) > 0, seconds


def main():
    import napari
    from magicgui import magicgui

    data, phi = load_data()
    viewer = napari.Viewer(ndisplay=2)
    viewer.add_image(data, name="nuclei", colormap="gray")
    points = viewer.add_points(POINTS, name="points", face_color="cyan", size=3)
    viewer.dims.set_point(0, int(POINTS[0, 0]))

    @magicgui(call_button="Calculate", methods={"widget_type": "Select", "choices": METHODS, "value": list(METHODS)})
    def calculate_all(methods=list(METHODS), align_surface=False):
        if len(points.data) != 2:
            print("Put exactly two points on the surface of a nucleus (now %d)." % len(points.data))
            return
        p1, p2 = (np.asarray(p, float) for p in points.data)
        if np.any(p1 < -0.5) or np.any(p2 < -0.5) or np.any(p1 > np.array(data.shape) - 0.5) or np.any(p2 > np.array(data.shape) - 0.5):
            print("The points have to be inside of the volume of shape %s." % (data.shape,))
            return
        for method in methods:
            try:
                mask, plane, seconds = calculate(method, phi, data, p1, p2, align_surface)
            except Exception as error:  # a point pair the algorithm cannot handle should not close the viewer
                print("%-12s failed: %s" % (method, error))
                continue
            print("%-12s %5.2f s   %6d voxels   plane %5d pixels" % (method, seconds, mask.sum(), plane.sum()))
            for name, layer_data in (("%s" % method, mask.astype(np.uint8)), ("%s plane" % method, plane.astype(np.uint8))):
                if name in viewer.layers:
                    viewer.layers[name].data = layer_data
                else:
                    viewer.add_labels(layer_data, name=name, opacity=0.5, visible=not name.endswith("plane"))
        points.visible = True

    viewer.window.add_dock_widget(calculate_all, area="right", name="Initial plane methods")
    napari.run()


if __name__ == "__main__":
    main()
