from distutils.core import setup, Extension, DEBUG
import glob
import os
import sys
import numpy as np


try:
    from Cython.Build import cythonize
except ImportError:
    cythonize = None


# https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html#distributing-cython-modules
def no_cythonize(extensions, **_ignore):
    for extension in extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = os.path.splitext(sfile)
            if ext in (".pyx", ".py"):
                if extension.language == "c++":
                    ext = ".cpp"
                else:
                    ext = ".c"
                sfile = path + ext
            sources.append(sfile)
        extension.sources[:] = sources
    return extensions


if "SIMPLEITK_PATH" in os.environ:
    sitk_path = os.environ["SIMPLEITK_PATH"]
else:
    if sys.platform == "win32":
        sitk_path = "C:\\SimpleITK-build"
        extra_build_args = ["/Zc:externC"]
    else:
        sitk_path = os.path.expanduser("~/SimpleITK-build")
        os.environ["CC"] = "g++"

sitk_lib_path = os.path.join(sitk_path, "lib")
if sys.platform == "win32":
    itk_lib_path = os.path.join(sitk_path, "ITK-build", "lib", "MinSizeRel")
    itk_libs = glob.glob(os.path.join(itk_lib_path, "*.lib"))
    sitk_libs = glob.glob(os.path.join(sitk_lib_path, "*.lib"))
    os_libs = [
        "rpcrt4.lib",
        "comctl32.lib",
        "wsock32.lib",
        "ws2_32.lib",
        "dbghelp.lib",
        "psapi.lib",
        "kernel32.lib",
        "user32.lib",
        "gdi32.lib",
        "winspool.lib",
        "shell32.lib",
        "ole32.lib",
        "oleaut32.lib",
        "uuid.lib",
        "comdlg32.lib",
        "advapi32.lib",
    ]
    extra_compile_args = ["/std:c++17"]
else:
    itk_lib_path = os.path.join(sitk_path, "ITK-build", "lib")
    itk_libs = glob.glob(os.path.join(itk_lib_path, "*.a"))
    sitk_libs = glob.glob(os.path.join(sitk_lib_path, "*.*"))
    extra_compile_args = ["-std=c++17", "-fopenmp"]
    os_libs = []
libs = itk_libs+sitk_libs
libs = map(os.path.basename, libs)
libs = map(lambda s: s.rsplit(".", 1)[0], libs)
libs = list(map(lambda s: s[3:] if s.startswith("lib") else s, libs))
print("libs:")
for l in libs:
    print(" - %s" % l)
print("extra_objects:")
for l in itk_libs+sitk_libs:
    print(" - %s" % l)
print("library dirs: %s, %s" % (sitk_lib_path, itk_lib_path))
extension = Extension(
    'minimal_surface',
    sources=["src/minimal-surface/code/_minimal_surface.pyx"] + glob.glob("src/minimal-surface/code/eikonal/*.cpp"),
    include_dirs=[
        np.get_include(),
        "src/minimal-surface/code",
        "src/minimal-surface/code/eikonal",
        "src/minimal-surface/code/pythonwrapper",
        "src/minimal-surface/code/eigen-3.4.0",
        glob.glob(os.path.join(sitk_path, "include", "SimpleITK-*"))[0]
    ],
    depends=["MinimalSurfaceEstimator.h", "SimpleITK.h", "sitkImage.h"],
    library_dirs=[sitk_lib_path, itk_lib_path],
    libraries=libs,
    extra_objects=itk_libs+sitk_libs+os_libs,
    extra_compile_args=extra_compile_args,
    language="c++"
)

CYTHONIZE = cythonize is not None

if CYTHONIZE:
    compiler_directives = {"language_level": 3, "embedsignature": True}
    extensions = cythonize([extension], compiler_directives=compiler_directives)
else:
    extensions = no_cythonize([extension])

setup(name='minimal_surface',
    description='Python Package with minimalsurface C++ extension',
    ext_modules=extensions
)

