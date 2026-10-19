# Minimal surface

`minimal-surface` is a Python package which can be utilized for the segmentation of objects in 3D images using as few as 2 points and a 2D slice! This algorithm uses energy minimization to estimate the surface of the object. To achieve fast execution the method was implemented in C++, but can be called from Python through a thin wrapper.

The C++ algorithm was implemented by Jozsef Molnar. Python wrapper was implemented by David Bauer.

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

### Supported platforms
Wheels are published to PyPI for:

| | Python 3.10 | 3.11 | 3.12 | 3.13 | 3.14 |
|---|---|---|---|---|---|
| Windows (x86-64) | yes | yes | yes | yes | yes |
| macOS 11 and later, Apple silicon and Intel | yes | yes | yes | yes | yes |
| Linux (x86-64, glibc 2.28 or later: Ubuntu 20.04+, Debian 10+, RHEL 8+) | yes | yes | yes | yes | yes |

Python 3.9 and older are not supported (0.0.2 was the last release for Python 3.9). On other platforms, such as
Linux on ARM, `pip` builds the package from source, which needs the SimpleITK libraries described below.
The supported Python versions are listed in `setup.cfg` (`python_requires` and the classifiers), and
[`.github/workflows/build.yml`](.github/workflows/build.yml) builds and tests every one of them on each platform.

### Building from source
The extension is written in C++ and links against the SimpleITK C++ libraries, which `setup.py` finds through
the `SIMPLEITK_PATH` environment variable. It can point at either

* a prefix that has conda-forge's `libsimpleitk` and `libitk-devel` installed (a conda or micromamba environment),
  which is the quickest way to get a working SDK on Linux, macOS and Windows, or
* a SimpleITK [SuperBuild](https://simpleitk.readthedocs.io/en/master/building.html) directory (static ITK and
  SimpleITK libraries), for example `C:\SimpleITK-build` on Windows.

With micromamba (on Windows leave out `compilers`: Visual Studio is used instead; on macOS add `llvm-openmp`):
```
micromamba create -n minimal-surface -c conda-forge python=3.13 libsimpleitk=2.5 libitk-devel numpy pip pytest compilers
micromamba activate minimal-surface
export SIMPLEITK_PATH=$CONDA_PREFIX      # Windows (cmd): set SIMPLEITK_PATH=%CONDA_PREFIX%
python -m pip install .
python -m pytest tests
```
The build is checked on Linux, macOS (Apple silicon and Intel) and Windows by
[`.github/workflows/build.yml`](.github/workflows/build.yml), which is the reference for the exact steps.
It also repairs each wheel so that it carries the ITK libraries it needs, and tests the result on a plain Python
installation. Publishing a GitHub release uploads the Windows and macOS wheels and the source distribution to PyPI
([`.github/workflows/release.yml`](.github/workflows/release.yml)).

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

# Two points on the surface of the nucleus. Points use the same axis order as
# the array indices, (z, y, x), just like `data[z, y, x]`.
# Both lie in slice z=35, which is the one displayed at the end.
p1 = np.asarray([ 35., 60., 20.])
p2 = np.asarray([ 35., 12., 39.])

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

## Inputs
`calculate`, `calc_eikonal_and_transport_init` and `init_transport_slice` take 3D arrays of any real numeric dtype (bool, integer or floating point) and any memory layout; they are converted to C-contiguous `float64` internally, and your arrays are not modified. Anything else (complex, strings, objects, `None` for `phi`) raises a `TypeError`, and `phi` and `image` must have the same shape.

Points are `(z, y, x)`, like array indices, and may be fractional. A point must lie inside the volume (each coordinate in `[-0.5, n - 0.5)` for an axis of length `n`); a point of the wrong length, or outside the array on any axis, raises a `ValueError` naming the point and the shape.

## Issues
If you bump into any problems regarding the package, please [file an issue].

[file an issue]: https://github.com/bauerdavid/minimal-surface/issues

## License
minimal-surface is released under the BSD 3-Clause license (see [LICENSE](LICENSE)).

The wheels also contain third-party libraries (ITK and what it depends on, the OpenMP runtime) and are built with
SimpleITK and Eigen. Their licenses are collected in `THIRD_PARTY_LICENSES.txt`, which is stored in the wheel's
`.dist-info/licenses` directory.
