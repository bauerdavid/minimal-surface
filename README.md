# Minimal surface

`minimal-surface` is a Python package which can be utilized for the segmentation of objects in 3D images using as few as 2 points and a 2D slice! This algorithm uses energy minimization to estimate the surface of the object. To achieve fast execution the method was implemented in C++, but can be called from Python through a thin wrapper.

## Citation
**When the pen is mightier than the sword: semi-automatic 2 and 3D image labelling**<br>
Réka Hollandi, David Bauer, Akos Diosdi, Bálint Schrettner, Timea Toth, Dominik Hirling, Gábor Hollandi, Maria Harmati, József Molnár, Peter Horvath<br>
bioRxiv 2024.01.15.575658<br>
doi: https://doi.org/10.1101/2024.01.15.575658

## Installation
It can be installed from PyPI using `pip`:
```
python -m pip install minimal-surface
```

> **Note**: This package is currently available only for Windows. In the future we plan to release it for Linux and Mac systems.

## Usage
The easiest way to utilize our tool is by using [Annotation Toolbox](https://github.com/bauerdavid/napari-nD-annotator), a [napari](https://www.napari.org) plugin created for fast 2 and 3D image annotation.

### Demo


https://github.com/bauerdavid/minimal-surface/assets/36735863/36cc3e18-87c3-4f6d-ae64-fc0e317e1cf4



## Example
In order to run this example, you'll need these packages:
* numpy
* scipy
* scikit-image
* pooch
* matplotlib
```python
import minimal_surface
from skimage.data import cells3d
from scipy import ndimage
import matplotlib.pyplot as plt
import numpy as np

def calc_features(image):
    delta = 0.1
    max_ = np.quantile(image, .95)
    min_ = image.min()
    image = np.clip((image-min_)/(max_-min_), 0, 1)

    # This feature will consider both image gradient and intensity
    # (giving higher values when the gradient is high and the intensity is low)
    gradient = ndimage.gaussian_gradient_magnitude(image, (1., 1., 1.))
    weights = 1+(delta-1)*image**2
    phi = gradient*weights +delta*(1-image**2)
    phi = (phi-(min_:=phi.min()))/(phi.max()-min_)
    return phi

# Crop the part to the smallest size possible to prevent long computation
data = cells3d()[:, 1, 60:130, 130:190] 

# Two points on the surface of the nucleus
p1 = np.asarray([ 20., 60., 35.])
p2 = np.asarray([ 39., 12., 35.])

data = (data - (min_ := data.min())) / (data.max() - min_)
data = ndimage.gaussian_filter(data, 1.)

features = calc_features(data)
alpha = 3e-3
phi = alpha + (1-alpha)*np.exp(-5*features)
phi = phi/phi.max()

calculator = minimal_surface.MinimalSurfaceCalculator()
def segment_slice(image_slice, distance_map):
    # This function segments the input image slice.
    # Asking for user input is also possible for accurate segmentation.
    return image_slice >= 0.24

calculator.set_initial_plane_calculator(segment_slice)
# If you have multiple objects, you can run the first part separately. which requires user input,
# then finish the computationally expensive part together for the whole dataset. 
# calculator.calc_eikonal_and_transport_init(phi, data, p1, p2, True)

transport_function = calculator.calculate(phi, data, p1, p2)
print("calculated")
mask = transport_function >= 0.
# The segmentation is contained in mask
plt.imshow(mask[35], cmap="plasma")
plt.show()
```

## Issues
If you bump into any problems regarding the package, please [file an issue].

[file an issue]: https://github.com/bauerdavid/minimal-surface/issues
