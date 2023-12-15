# import SimpleITK
import scipy.ndimage
import skimage.measure
from PyQt5.QtCore import QObject, pyqtSignal, Qt, QMutex
from PyQt5 import QtGui
from PyQt5.QtWidgets import QComboBox, QVBoxLayout, QSlider, QSpinBox, QWidget, QPushButton, QCheckBox, QLabel, \
    QDoubleSpinBox, QSizePolicy, QFrame
import numpy as np
# import tifffile
# import sys
import napari
from napari.layers import Image, Points, Shapes, Labels
# from napari.qt.threading import thread_worker
from napari.utils.notifications import show_info
from napari.utils.colormaps.standardize_color import transform_color
# from magicgui.widgets import FunctionGui
from threading import Thread
# from multiprocessing import Lock
from skimage.data import cells3d
from napari_nd_annotator.boundingbox import BoundingBoxLayer
import itertools
import MinArea
import cv2
from copy import deepcopy
from scipy.ndimage.filters import generic_filter
from scipy.signal.windows import gaussian
import SimpleITK as sitk
from typing import Union
c_rot = 0

blur_func = lambda img, val: scipy.ndimage.gaussian_filter(img, val)
# blur_func = lambda img, _: sitk.GetArrayFromImage(sitk.CurvatureAnisotropicDiffusion(sitk.GetImageFromArray(img.astype(float)), 0.0625, 9., numberOfIterations=10))
def bilateral_3d(image: np.ndarray, filter_size, sigmaSpace, sigmaIntensity):
    ndim = image.ndim
    if type(filter_size) is not tuple:
        filter_size = (filter_size,) * ndim
    kernel = None
    def gauss_fun(x):
        return 1/(sigmaIntensity*np.sqrt(2*np.pi)) * np.exp(-(1/2)*(x**2)/sigmaIntensity**2)
    for d, s in enumerate(filter_size):
        new_kernel = gaussian(s, std=sigmaSpace).reshape(tuple(s if i == d else 1 for i in range(ndim)))
        kernel = new_kernel if kernel is None else kernel[..., np.newaxis] @ new_kernel[np.newaxis]
    print(kernel.shape)

    def bilateral_filter(a):
        a = a.reshape(filter_size)
        gauss_term = kernel*gauss_fun(a)
        return a*gauss_term/gauss_term.sum()

    return generic_filter(image, bilateral_filter, size=filter_size)

def layer_to_sitk_image(layer: Union[Image, Labels]):
    if type(layer) not in [Image, Labels]:
        raise TypeError("layer should be either 'Image' or 'Labels'")
    img = sitk.GetImageFromArray(layer.data)
    img.SetOrigin(np.flip(layer.translate))
    img.SetDirection(np.flip(layer.rotate.reshape(-1)))
    return img

def sitk_image_to_layer(img: sitk.Image, cls=Image, **kwargs):
    layer = cls(sitk.GetArrayFromImage(img), **kwargs)
    layer.translate = np.flip(img.GetOrigin())
    layer.rotate = np.flip(img.GetDirection()).reshape(3, 3)
    return layer
def generate_label_colors(n):
    import colorsys
    from random import shuffle, seed
    rgb_tuples = [colorsys.hsv_to_rgb(x * 1.0 / n, 1., 1.) for x in range(n)]
    rgb_tuples = [(int(rgb_tuples[i][0] * 255), int(rgb_tuples[i][1] * 255), int(rgb_tuples[i][2] * 255)) for i in
                  range(n)]
    seed(0)
    shuffle(rgb_tuples)
    return rgb_tuples


def color_to_hex_string(rgb_tuple):
    return '#%.2X%.2X%.2X' % rgb_tuple


def get_bb_corners(bb):
    return np.concatenate([bb.min(0, keepdims=True), bb.max(0, keepdims=True)])


def bb_2_slice(bb):
    bb_corners = get_bb_corners(bb).round().astype(int)
    return tuple(slice(bb_corners[0, i], bb_corners[1, i]) for i in range(3))


class ColorPairsCallback:
    def __init__(self, n_colors=50):
        self.prev_selection = set()
        self.prev_len = 0
        self.prev_data = None
        color_cycle = list(map(lambda color_tuple: color_to_hex_string(color_tuple), generate_label_colors(n_colors)))
        self.color_cycle = itertools.cycle(color_cycle)

    def __call__(self, event):
        points_layer = event.source
        if len(points_layer.selected_data) > 1:
            points_layer.selected_data = self.prev_selection
            points_layer.refresh()
            return
        self.prev_selection = points_layer.selected_data.copy()
        if event.type == "data":
            if len(points_layer.data) > self.prev_len:
                if len(points_layer.data) % 2 == 1:
                    points_layer.current_face_color = next(self.color_cycle)
                    points_layer.face_color[-1] = transform_color(points_layer.current_face_color)
                self.prev_data = points_layer.data.copy()
            else:
                removed_idx = np.squeeze(np.argwhere(
                    np.all(np.any(~np.equal(self.prev_data[np.newaxis], points_layer.data[:, np.newaxis]), -1), axis=0)))
                if removed_idx.ndim == 0:
                    removed_idx = int(removed_idx)
                    points_layer.selected_data.clear()
                    if removed_idx % 2 == 0 and removed_idx < len(points_layer.data):
                        points_layer.selected_data.add(removed_idx)
                    elif removed_idx % 2 == 1:
                        points_layer.selected_data.add(removed_idx - 1)
                    points_layer.remove_selected()
            points_layer.refresh()
            self.prev_data = points_layer.data
            self.prev_len = len(points_layer.data)


center_mat = None
rot_mat = None
transl_mat = None


class BlurSlider(QSlider):
    def __init__(self, parent, layer_name, blur_func=None):
        super().__init__()
        self.layer = None
        self.layer_name = layer_name
        self.parent = parent
        self.image = None
        self.valueChanged.connect(self.update_image)
        self.setOrientation(Qt.Horizontal)
        self.blur_func = blur_func if blur_func is not None else lambda img, val: scipy.ndimage.gaussian_filter(img, val)

    def setMaximum(self, a0: float) -> None:
        super().setMaximum(a0*10)

    def setMinimum(self, a0: float) -> None:
        super().setMinimum(a0*10)

    def setValue(self, a0: float) -> None:
        super().setValue(a0*10)

    def value(self):
        return super().value()/10

    def mousePressEvent(self, ev: QtGui.QMouseEvent) -> None:
        super().mousePressEvent(ev)
        image_layer = self.parent.image_layer
        points_layer = self.parent.points_layer
        if image_layer is None or points_layer is None or len(points_layer.data) < 2:
            return
        else:
            bb = pts_2_bb(points_layer.data[0], points_layer.data[1])
            bb = np.clip(bb, 0, np.asarray(image_layer.data.shape) - 1).round().astype(int)
            bb_slice = bb_2_slice(bb)
            self.image = image_layer.data[bb_slice]
            offset = bb.min(0)

            blurred = self.blur_func(self.image, self.value())
            self.layer = Image(blurred, name=self.layer_name, translate=image_layer.translate, colormap=image_layer.colormap, contrast_limits=image_layer.contrast_limits)
            self.layer.translate = offset
            self.parent.viewer.add_layer(self.layer)

    def mouseReleaseEvent(self, ev: QtGui.QMouseEvent) -> None:
        super().mouseReleaseEvent(ev)
        if self.layer is None:
            return
        self.parent.viewer.layers.remove(self.layer)
        self.layer = None
        self.image = None

    def update_image(self, _):
        if self.image is not None:
            self.layer.data = self.blur_func(self.image, self.value())
            self.layer.refresh()


class EstimatorWidget(QWidget):
    image_data_received = pyqtSignal(str, "PyQt_PyObject", "PyQt_PyObject")
    mask_data_received = pyqtSignal("PyQt_PyObject", "PyQt_PyObject")
    remove_layer = pyqtSignal(str)
    layer_invalidated = pyqtSignal(str)

    def __init__(self, viewer: napari.Viewer):
        super().__init__()
        self.image_data_received.connect(self._add_image)

        self.layer_invalidated.connect(self.refresh_layer)
        self.mask_data_received.connect(self.mask_received)
        self.remove_layer.connect(self.data_remover)
        self.viewer = viewer
        self.points_layer = self.viewer.add_points(ndim=3, size=2, name="[surface points]")
        self.points_layer.events.connect(ColorPairsCallback())
        self.points_layer.mouse_drag_callbacks.append(self.drag_points_callback)
        self.bounding_box_layer = BoundingBoxLayer(ndim=3, face_color="transparent", edge_color="green", name="[bounding boxes]")
        self.mask_layer = None
        self.viewer.add_layer(self.bounding_box_layer)
        self.new_layer = None
        self.previous_index = None
        layout = QVBoxLayout()
        layout.addWidget(QLabel("Image"))
        self.image_layer_dropdown = QComboBox()
        layout.addWidget(self.image_layer_dropdown)

        layout.addWidget(QLabel("Alpha"))
        self.alpha_slider = QSlider()
        self.alpha_slider.setMinimum(0)
        self.alpha_slider.setMaximum(100)
        self.alpha_slider.setTickInterval(1)
        self.alpha_slider.setValue(1)
        self.alpha_slider.setOrientation(Qt.Horizontal)
        layout.addWidget(self.alpha_slider)

        layout.addWidget(QLabel("Beta"))
        self.beta_spinner = QDoubleSpinBox()
        self.beta_spinner.setMinimum(1.)
        self.beta_spinner.setMaximum(1000.)
        self.beta_spinner.setValue(30.)
        self.beta_spinner.setSingleStep(1.)
        layout.addWidget(self.beta_spinner)

        self.blur_label = QLabel("blur sigma: 0.")
        layout.addWidget(self.blur_label)
        self.blur_sigma_spinner = BlurSlider(self, "crop", blur_func)
        self.blur_sigma_spinner.setMinimum(0.)
        self.blur_sigma_spinner.setMaximum(10.)
        self.blur_sigma_spinner.valueChanged.connect(self.blur_changed)
        self.blur_sigma_spinner.setValue(1.5)
        layout.addWidget(self.blur_sigma_spinner)

        self.use_correction_checkbox = QCheckBox("Use correction")
        self.use_correction_checkbox.setChecked(True)
        layout.addWidget(self.use_correction_checkbox)

        self.call_button = QPushButton("Run")
        self.call_button.clicked.connect(lambda: self.call_button.setDisabled(True))
        self.call_button.clicked.connect(self.start_estimation)
        layout.addWidget(self.call_button)
        self.setLayout(layout)
        viewer.layers.events.connect(self.update_layers)
        self.update_layers()
        self.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Maximum)

    @property
    def image_layer(self):
        layer_name = self.image_layer_dropdown.currentText()
        return None if layer_name == '' else self.viewer.layers[layer_name]

    def blur_changed(self, value):
        self.blur_sigma = value/10
        self.blur_label.setText("Blur sigma: %f" % self.blur_sigma)

    def _add_image(self, layer_name, data, layer_args):
        if 'opacity' not in layer_args:
            layer_args['opacity'] = 0.3
        if 'blending' not in layer_args:
            layer_args['blending'] = "additive"
        if 'rgb' not in layer_args:
            layer_args['rgb'] = False
        self.viewer.add_layer(Image(data, name=layer_name, **layer_args))
        print(layer_name, 'added')

    def mask_received(self, data, offset):
        offset = offset.round().astype(int)
        bb_slice = tuple(slice(offs, offs+extent) for offs, extent in zip(offset, data.shape))
        self.mask_layer.data[bb_slice][data > 0] = data[data > 0]
        self.mask_layer.refresh()

    def refresh_layer(self, layer_name):
        if layer_name not in self.viewer.layers:
            return
        layer = self.viewer.layers[layer_name]
        layer.refresh()

    def start_estimation(self):
        if self.image_layer_dropdown.currentText() not in self.viewer.layers:
            self.call_button.setDisabled(False)
            raise ValueError("Missing image layer")
        elif self.points_layer is None:
            self.call_button.setDisabled(False)
            raise ValueError("Missing points layer")
        elif self.bounding_box_layer is None:
            self.call_button.setDisabled(False)
            raise ValueError("Missing bounding box layer")
        elif type(self.viewer.layers.selection.active) != Labels:
            self.call_button.setDisabled(False)
            raise TypeError("Selected layer should be of type 'Labels'")
        self.mask_layer = self.viewer.layers.selection.active
        image = self.viewer.layers[self.image_layer_dropdown.currentText()]
        points = self.points_layer
        alpha = self.alpha_slider.value() / 100
        beta = self.beta_spinner.value()
        use_correction = self.use_correction_checkbox.isChecked()

        # viewer.add_image(data)
        viewer.window.worker = Thread(target=self.estimate, args=[image.data, points.data, use_correction, beta, alpha],
                                      daemon=True)
        viewer.window.worker.start()

    def hook_callbacks(self, estimator, stage, layer_name, layer_params, data_idx=None):

        def initializer(arr, idx):
            self.data_initializer(layer_name, data_idx, layer_params)(arr, idx)
        estimator.hook_stage_data_init_event(stage, initializer)
        estimator.hook_stage_iteration_event(stage, self.data_updater(layer_name))
        estimator.hook_stage_finished_event(stage, self.data_finalizer(layer_name))

    def estimate(self, image, points, use_correction, beta, alpha):
        import time
        # time.sleep(5)
        print("started")
        postscript = ""
        estimator = MinArea.Estimator()
        estimator.hook_stage_data_init_event(
            MinArea.AREA_EIKONAL_STAGE,
            lambda arr, idx: self.data_initializer(
                "Area Eikonal %d%s" % (0, postscript),
                0,
                {
                    'colormap': "plasma",
                    'translate': offset,
                    "visible": False
                }
            )(arr, idx)
        )
        estimator.hook_stage_iteration_event(MinArea.AREA_EIKONAL_STAGE,
                                            lambda idx: self.data_updater("Area Eikonal %d%s" % (0, postscript))(idx))
        estimator.hook_stage_finished_event(MinArea.AREA_EIKONAL_STAGE,
                                           lambda: self.data_finalizer("Area Eikonal %d%s" % (0, postscript))())
        estimator.hook_stage_data_init_event(
            MinArea.AREA_EIKONAL_STAGE,
            lambda arr, idx: self.data_initializer(
                "Area Eikonal %d%s" % (1, postscript),
                1,
                {
                    'colormap': "plasma",
                    'translate': offset,
                    "visible": False
                }
            )(arr, idx)
        )
        estimator.hook_stage_iteration_event(MinArea.AREA_EIKONAL_STAGE,
                                            lambda idx: self.data_updater("Area Eikonal %d%s" % (1, postscript))(idx))
        estimator.hook_stage_finished_event(MinArea.AREA_EIKONAL_STAGE,
                                            lambda: self.data_finalizer("Area Eikonal %d%s" % (1, postscript))())
        #self.hook_callbacks(estimator, MinArea.AREA_EIKONAL_STAGE, "Area Eikonal 0", {'colormap': "plasma", 'translate': offset, "visible": False}, 0)
        #self.hook_callbacks(estimator, MinArea.AREA_EIKONAL_STAGE, "Area Eikonal 1", {'colormap': "plasma", 'translate': offset, "visible": False}, 1)
        # self.hook_callbacks(estimator, MinArea.ROTATED_AREA_EIKONAL_STAGE, "Rotated Area Eikonal", {'colormap': "plasma"}, 0)

        '''estimator.hook_stage_data_init_event(
            MinArea.ROTATED_AREA_EIKONAL_STAGE,
            lambda arr, idx: self.data_initializer("Rotated Area Eikonal %d%s" % (0, postscript),
                                                   0,
                                                   {
                                                       'colormap': "plasma",
                                                       'translate': offset + transl_mat,
                                                       'rotate': rot_mat,
                                                       'visible': False
                                                   })(arr, idx)
        )
        estimator.hook_stage_iteration_event(MinArea.ROTATED_AREA_EIKONAL_STAGE,
                                             lambda idx: self.data_updater("Rotated Area Eikonal %d%s" % (0, postscript))(idx))
        estimator.hook_stage_finished_event(MinArea.ROTATED_AREA_EIKONAL_STAGE,
                                            lambda: self.data_finalizer("Rotated Area Eikonal %d%s" % (0, postscript))())'''
        '''estimator.hook_stage_data_init_event(
            MinArea.ROTATED_AREA_EIKONAL_STAGE,
            lambda arr, idx: self.data_initializer("Rotated Area Eikonal %d%s" % (1, postscript),
                                                   1,
                                                   {
                                                       'colormap': "plasma",
                                                       'translate': offset + transl_mat,
                                                       'rotate': rot_mat,
                                                       'visible': False
                                                   })(arr, idx)
        )
        estimator.hook_stage_iteration_event(MinArea.ROTATED_AREA_EIKONAL_STAGE,
                                             lambda idx: self.data_updater("Rotated Area Eikonal %d%s" % (1, postscript))(idx))
        estimator.hook_stage_finished_event(MinArea.ROTATED_AREA_EIKONAL_STAGE,
                                            lambda: self.data_finalizer("Rotated Area Eikonal %d%s" % (1, postscript))())'''

        estimator.hook_stage_data_init_event(
            MinArea.PLANE_PHASEFIELD_STAGE,
            lambda arr, idx: self.data_initializer("Plane PhaseField%s" % postscript,
                                              layer_args={
                                                  'colormap': "plasma",
                                                  'translate': offset + transl_mat +
                                                               np.asarray([0, 0, center_mat[-1]]) @ rot_mat.T,
                                                  'rotate': rot_mat,
                                                  'opacity': 0.6,
                                                  "visible": False
                                              }
                                              )(arr, idx)
        )
        estimator.hook_stage_iteration_event(MinArea.PLANE_PHASEFIELD_STAGE, lambda idx: self.data_updater("Plane PhaseField%s" % postscript)(idx))
        estimator.hook_stage_finished_event(MinArea.PLANE_PHASEFIELD_STAGE, lambda: self.data_finalizer("Plane PhaseField%s" % postscript)())

        estimator.hook_stage_data_init_event(
            MinArea.TRANSPORT_FUNCTION_STAGE,
            lambda arr, idx: self.data_initializer("Transport Function%s" % postscript,
                                                   0,
                                                   {
                                                       'colormap': "turbo",
                                                       'translate': offset + transl_mat,
                                                       'rotate': rot_mat,
                                                       'opacity': 0.5,
                                                       'rendering': "iso",
                                                       "iso_threshold": 0.
                                                   })(arr, idx)
        )
        estimator.hook_stage_iteration_event(MinArea.TRANSPORT_FUNCTION_STAGE, lambda idx: self.data_updater("Transport Function%s" % postscript)(idx))
        estimator.hook_stage_finished_event(MinArea.TRANSPORT_FUNCTION_STAGE, lambda: self.data_finalizer("Transport Function%s" % postscript)())

        def tform_calculated(rotation, translation):
            global rot_mat, transl_mat
            rot_mat = np.flip(rotation.copy()).reshape(3, 3)
            transl_mat = np.flip(translation.copy())
            print(transl_mat)

        estimator.hook_transform_calculated_event(tform_calculated)

        def center_calculated(center):
            global center_mat
            center_mat = np.flip(center.copy())

        estimator.hook_plane_center_calculated_event(center_calculated)
        for i in range(len(points)//2):
            bounding_box = pts_2_bb(points[2*i], points[2*i+1])
            bounding_box = np.clip(bounding_box, 0, np.asarray(image.shape) - 1).round().astype(int)
            offset = bounding_box.min(0, keepdims=True)
            bb_slice = bb_2_slice(bounding_box)
            point1 = np.flip(points[2*i] - offset)
            point2 = np.flip(points[2*i+1] - offset)

            data = image[bb_slice]
            data = (data - data.min()) / (data.max() - data.min())
            # data = scipy.ndimage.gaussian_filter(data, self.blur_sigma)
            data = blur_func(data, self.blur_sigma)
            offset = np.clip(bounding_box.min(0), 0, np.asarray(self.image_layer.data.shape) - 1)
            # start = time.time()

            output = estimator.calculate(data, point1, point2, use_correction, beta, alpha)
            segmented = (output >=0)
            labelled = skimage.measure.label(segmented)
            obj_pixel = np.argwhere(output == output.max())[0]
            obj_label = labelled[tuple(obj_pixel)]
            mask = (labelled == obj_label)*(i+1)
            self.image_data_received.emit("Result", output, {"colormap": "plasma", "translate": offset})
            self.mask_data_received.emit(mask, offset)
            postscript = " - %d" % (i+1)
            # self.remove_layer.emit("Transport Function")
        show_info("Estimation done")
        self.call_button.setDisabled(False)

    def objectName(self):
        return "estimate widget"

    def update_layers(self, event=None):
        type_ = event.type if event else "reordered"
        index = getattr(event, "index", None)
        current_image = self.image_layer_dropdown.currentText()
        if type_ in ["reordered", "removed"]:
            img_idx = 0
            for layer in self.viewer.layers:
                if type(layer) == Image:
                    if img_idx >= self.image_layer_dropdown.count():
                        self.image_layer_dropdown.addItem(layer.name)
                    else:
                        self.image_layer_dropdown.setItemText(img_idx, layer.name)
                    if self.image_layer is None \
                            or self.image_layer_dropdown.itemText(img_idx) == current_image:
                        self.image_layer_dropdown.setCurrentIndex(img_idx)
                    img_idx += 1
            if type_ == "removed":
                if self.new_layer is not None:
                    if type(self.new_layer) == Points:
                        self.points_layer = self.new_layer
                    elif type(self.new_layer) == BoundingBoxLayer:
                        self.bounding_box_layer = self.new_layer
                    self.viewer.add_layer(self.new_layer)
                    if self.previous_index == len(self.viewer.layers) - 1:
                        self.viewer.layers.move(1, 0)
                        self.viewer.layers.move(1, 0)
                    self.viewer.layers.move(len(self.viewer.layers)-1, self.previous_index)
                    self.new_layer = None
                while img_idx < self.image_layer_dropdown.count():
                    self.image_layer_dropdown.removeItem(img_idx)
        elif type_ == "inserted":
            layer = self.viewer.layers[-1]
            if type(layer) == Image:
                self.image_layer_dropdown.addItem(layer.name)
        elif type_ == "removing":
            layer = self.viewer.layers[index]
            if layer in [self.points_layer, self.bounding_box_layer]:
                self.previous_index = index
                self.new_layer = deepcopy(layer)
                # self.new_layer.name = "[surface points]"

    def data_initializer(self, name, selected_idx=None, layer_args=None):
        if layer_args is None:
            layer_args = {}

        def initialize(arr, idx):
            try:
                print("initializing", name)
                if selected_idx in [None, idx]:
                    print(name, "initialized")
                    self.image_data_received.emit(name, arr, layer_args)
                else:
                    print("%s was not initialized as id (%d) was not %d" % (name, idx, selected_idx))
            except Exception as e:
                print(e)

        return initialize
        # return lambda arr, idx: print("initializing", name)

    def data_updater(self, name):
        def update_viewer(iteration):
            if iteration % 10 == 0:
                self.layer_invalidated.emit(name)

        return update_viewer

    def data_finalizer(self, name):
        def finalize():
            if name not in self.viewer.layers:
                print("%s not in layer list" % name)
                return
            self.viewer.layers[name].data = np.copy(self.viewer.layers[name].data)
            print(name, "finalized")

        return finalize

    def data_remover(self, name):
        self.viewer.layers.remove(name)

    def drag_points_callback(self, layer, event):
        if layer.mode != "select":
            return
        yield
        if event.type == "mouse_move":
            if len(layer.selected_data) > 2 or len(layer.selected_data) == 0:
                return
            if len(layer.selected_data) == 2:
                selection_iter = iter(layer.selected_data)
                idx1 = next(selection_iter)
                idx2 = next(selection_iter)
                if idx1//2 != idx2//2:
                    return
            else:
                index = next(iter(layer.selected_data))
                if index % 2 == 0:
                    idx1 = index
                    idx2 = index+1
                else:
                    idx1 = index - 1
                    idx2 = index
            p1 = layer.data[idx1]
            p2 = layer.data[idx2]
            bb = pts_2_bb(p1, p2)
            bb_slice = bb_2_slice(bb)
            image = self.image_layer.data[bb_slice]
            self.image_layer.visible = False
            cropped_image = self.viewer.add_image(image, colormap=self.image_layer.colormap, translate=bb.min(0),
                                                  rendering="additive", gamma=2.)
            self.viewer.layers.selection.active = layer
        yield
        while event.type == "mouse_move":
            yield
        self.viewer.layers.remove(cropped_image)
        self.image_layer.visible = True


class ShortestPathsWidget(QWidget):
    shapes_data_received = pyqtSignal(str, "PyQt_PyObject", "PyQt_PyObject")

    def __init__(self, viewer: napari.Viewer):
        super().__init__()
        self.shapes_data_received.connect(self.add_shapes)
        self.viewer = viewer
        self.estimator = MinArea.Estimator()
        layout = QVBoxLayout()

        layout.addWidget(QLabel("Extract contour from meeting plane"))
        self.meeting_plane_dropdown = QComboBox()
        layout.addWidget(self.meeting_plane_dropdown)
        plane_2_contour_btn = QPushButton("Extract")
        plane_2_contour_btn.clicked.connect(self.plane_2_contour)
        layout.addWidget(plane_2_contour_btn)
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        line.setFixedHeight(1)
        line.setStyleSheet("background-color: #414851;")
        layout.addWidget(line)
        self.distmap_layer_dropdown = QComboBox()
        layout.addWidget(QLabel("Distance map"))
        layout.addWidget(self.distmap_layer_dropdown)

        self.points_layer_dropdown = QComboBox()
        layout.addWidget(QLabel("Start points"))
        layout.addWidget(self.points_layer_dropdown)

        run_button = QPushButton("Show shortest paths")
        run_button.clicked.connect(self.show_shortest_paths)
        layout.addWidget(run_button)
        viewer.layers.events.connect(self.update_layer_list)
        self.setLayout(layout)
        self.update_layer_list(None)

    @property
    def distance_map(self):
        current_text = self.distmap_layer_dropdown.currentText()
        return self.viewer.layers[current_text] if current_text else None

    @property
    def points(self):
        current_text = self.points_layer_dropdown.currentText()
        return self.viewer.layers[current_text] if current_text else None

    @property
    def meeting_plane(self):
        current_text = self.meeting_plane_dropdown.currentText()
        return self.viewer.layers[current_text] if current_text else None

    def add_shapes(self, layer_name, data, layer_args):
        if 'opacity' not in layer_args:
            layer_args['opacity'] = 0.3
        self.viewer.add_layer(Shapes(data, name=layer_name, shape_type="path", **layer_args))

    def plane_2_contour(self):
        if self.distance_map is None or self.meeting_plane is None:
            return
        meeting_plane_layer = self.meeting_plane
        plane_mask = (meeting_plane_layer.data < 0).astype(np.uint8)
        contour = np.squeeze(cv2.findContours(plane_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)[0][0])
        contour = np.concatenate([np.fliplr(contour), np.ones([len(contour), 1])], axis=1)
        contour = contour @ meeting_plane_layer.rotate.T + meeting_plane_layer.translate
        contour_layer = Points(contour, ndim=self.distance_map.ndim, name="contour points", size=2)
        self.viewer.add_layer(contour_layer)
        for i in range(self.points_layer_dropdown.count()):
            item = self.points_layer_dropdown.itemText(i)
            if item == contour_layer.name:
                self.points_layer_dropdown.setCurrentIndex(i)
                break

    def update_layer_list(self, event):
        current_distance_map = self.distmap_layer_dropdown.currentText()
        self.distmap_layer_dropdown.clear()

        current_points = self.points_layer_dropdown.currentText()
        self.points_layer_dropdown.clear()

        current_meeting_plane = self.meeting_plane_dropdown.currentText()
        self.meeting_plane_dropdown.clear()

        for layer in self.viewer.layers:
            if isinstance(layer, Image):
                self.distmap_layer_dropdown.addItem(layer.name)
                self.meeting_plane_dropdown.addItem(layer.name)
                if layer.name == current_distance_map:
                    self.distmap_layer_dropdown.setCurrentIndex(self.distmap_layer_dropdown.count() - 1)
                if layer.name == current_meeting_plane:
                    self.meeting_plane_dropdown.setCurrentIndex(self.meeting_plane_dropdown.count() - 1)
            elif type(layer) is Points:
                self.points_layer_dropdown.addItem(layer.name)
                if layer.name == current_points:
                    self.points_layer_dropdown.setCurrentIndex(self.points_layer_dropdown.count() - 1)

    def calculate(self, points, distance_map, translate, rotate):
        paths = []
        for point in points:
            shortest_path = self.estimator.resolve_shortest_paths(point, distance_map).copy()
            paths.append(np.fliplr(shortest_path))
        colors = generate_label_colors(len(paths))
        self.shapes_data_received.emit("Shortest path", paths, {
            'translate': translate,
            'rotate': rotate,
            'edge_color': np.concatenate([np.asarray(colors)/255, np.ones([len(colors), 1])], axis=1),
            'edge_width': 0.1
        })

    def show_shortest_paths(self):
        distance_map = self.distance_map
        points = self.points
        if distance_map is None or points is None:
            return
        translation = distance_map.translate
        rotation = distance_map.rotate
        distance_map = distance_map.data
        points = np.fliplr(points.data - translation).astype(int)
        viewer.window.worker2 = Thread(target=self.calculate, args=[points, distance_map, translation, rotation])
        viewer.window.worker2.start()


def pts_2_bb(p1, p2):
    center = (p1 + p2) / 2
    size = np.sqrt(((p1 - p2) ** 2).sum()) * 1.2 + 10
    bb = np.asarray(np.where(list(itertools.product((False, True), repeat=3)), center + size / 2, center - size / 2))
    bb = np.clip(bb, 0, np.asarray(viewer.layers["image"].data.shape) - 1)
    return bb

class InitialContourCalculator(QWidget):
    def __init__(self, viewer):
        self.viewer = viewer


if __name__ == "__main__":
    viewer = napari.Viewer()
    # image = tifffile.imread(r"Y:\BIOMAG\shortest path\interm_imgs\backup\ph0_testimage.tif")
    image = cells3d()[:, 0]

    layer = viewer.add_image(image, colormap="magma")
    mask = viewer.add_labels(np.zeros_like(image, dtype=np.uint8), name="Segmentation")
    viewer.window.add_dock_widget(EstimatorWidget(viewer))
    viewer.window.add_dock_widget(ShortestPathsWidget(viewer))
    napari.run()
