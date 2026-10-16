from setuptools import setup, Extension
import glob
import os
import subprocess
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


WINDOWS = sys.platform == "win32"
MACOS = sys.platform == "darwin"

# Windows system libraries that the static ITK/SimpleITK libraries depend on.
WINDOWS_SYSTEM_LIBS = [
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


def _glob(directory, *patterns):
    found = set()
    for pattern in patterns:
        found.update(glob.glob(os.path.join(directory, pattern)))
    return sorted(found)


def _link_group(libraries):
    """Wrap static libraries so that GNU ld can resolve references between them.

    ld only looks at an archive once, in command line order, and the ITK and
    SimpleITK archives reference each other in both directions; a group makes it
    rescan them until nothing more resolves. Other linkers do not need (or accept) it.
    """
    if WINDOWS or MACOS:
        return list(libraries)
    return ["-Wl,--start-group", *libraries, "-Wl,--end-group"]


def describe_sdk(root):
    """Work out how to compile and link against the SimpleITK installation at root.

    Two layouts are supported:

    * a SimpleITK SuperBuild tree (ITK-build/ next to lib/ and include/), where SimpleITK
      and ITK are both static libraries. This is what the Windows wheels were built against.
    * an installed prefix, for example a conda environment holding conda-forge's
      libsimpleitk (static SimpleITK libraries) and libitk-devel (shared ITK libraries).
      On Windows such a prefix keeps its headers and libraries under Library/.

    Returns the keyword arguments for the pieces of Extension that depend on the SDK.
    """
    root = os.path.abspath(root)
    static_ext = ".lib" if WINDOWS else ".a"
    sitk = ("" if WINDOWS else "lib") + "SimpleITK*" + static_ext
    sdk = {
        "include_dirs": [],
        "library_dirs": [],
        "extra_objects": [],
        "extra_link_args": [],
    }

    if os.path.isdir(os.path.join(root, "ITK-build")):
        itk_dir = os.path.join(root, "ITK-build", "lib", "MinSizeRel" if WINDOWS else "")
        itk_libs = _glob(itk_dir, "*" + static_ext)
        sitk_libs = _glob(os.path.join(root, "lib"), "*" + static_ext)
        include_roots = [os.path.join(root, "include")]
        lib_dirs = [os.path.join(root, "lib"), itk_dir]
        static_libs = itk_libs + sitk_libs
        shared_libs = []
    else:
        subdirs = ["Library"] if WINDOWS else []
        include_roots = [os.path.join(root, *subdirs, "include")]
        lib_dir = next(
            (d for d in (os.path.join(root, *subdirs, "lib"), os.path.join(root, "lib"))
             if _glob(d, sitk)),
            None,
        )
        if lib_dir is None:
            raise RuntimeError(
                "No SimpleITK libraries (%s) found under %s. Set SIMPLEITK_PATH to a SimpleITK "
                "SuperBuild directory or to a prefix that has libsimpleitk installed." % (sitk, root))
        lib_dirs = [lib_dir]
        static_libs = _glob(lib_dir, sitk)
        if WINDOWS:
            shared_libs = _glob(lib_dir, "ITK*.lib", "itk*.lib")
        elif MACOS:
            shared_libs = _glob(lib_dir, "libITK*.dylib", "libitk*.dylib")
        else:
            shared_libs = _glob(lib_dir, "libITK*.so", "libitk*.so")
        if not shared_libs:
            raise RuntimeError("No ITK libraries found in %s; install libitk-devel next to libsimpleitk." % lib_dir)
        if not WINDOWS:
            # let the extension find the shared ITK libraries again at run time
            sdk["extra_link_args"].append("-Wl,-rpath," + lib_dir)

    include_dirs = [d for r in include_roots for d in sorted(glob.glob(os.path.join(r, "SimpleITK-*")))]
    if not include_dirs:
        raise RuntimeError("SimpleITK headers (SimpleITK-*) not found under %s" % include_roots)
    sdk["include_dirs"] = include_dirs[-1:]
    sdk["library_dirs"] = lib_dirs

    if WINDOWS:
        sdk["extra_objects"] = static_libs + shared_libs + WINDOWS_SYSTEM_LIBS
    else:
        sdk["extra_link_args"] = _link_group(static_libs) + shared_libs + sdk["extra_link_args"]
    return sdk


def openmp_flags(root):
    """(compile arguments, link arguments) that turn on OpenMP, or empty lists if it is unavailable."""
    if WINDOWS:
        # MSVC implements OpenMP 2.0. Without /openmp it ignores every `#pragma omp`, which is
        # how the Windows wheels on PyPI (0.0.1 and 0.0.2) were built: their OpenMP loops run on
        # one thread.
        return ["/openmp"], []
    if not MACOS:
        return ["-fopenmp"], ["-fopenmp"]
    # Apple's clang has no OpenMP. conda's clang, together with llvm-openmp, does.
    if os.path.exists(os.path.join(root, "lib", "libomp.dylib")):
        return ["-fopenmp"], ["-fopenmp"]
    # Otherwise use Homebrew's libomp, as the napari-nD-annotator build does.
    try:
        prefix = subprocess.run(["brew", "--prefix", "libomp"], capture_output=True, text=True,
                                check=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("setup.py: OpenMP not found; building without it")
        return [], []
    return (["-Xpreprocessor", "-fopenmp", "-I%s/include" % prefix],
            ["-L%s/lib" % prefix, "-lomp"])


if "SIMPLEITK_PATH" in os.environ:
    sitk_path = os.environ["SIMPLEITK_PATH"]
else:
    if WINDOWS:
        sitk_path = "C:\\SimpleITK-build"
    else:
        sitk_path = os.path.expanduser("~/SimpleITK-build")
        os.environ["CC"] = "g++"

sdk = describe_sdk(sitk_path)
omp_compile_args, omp_link_args = openmp_flags(sitk_path)
extra_compile_args = (["/std:c++17"] if WINDOWS else ["-std=c++17"]) + omp_compile_args

# Compile against the modern NumPy C API only, so that use of the NumPy 1
# legacy API is a build error rather than a silent numpy 1 pin.
define_macros = [("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]

# Function-level profiling is opt-in. It is only reachable through
# profile_manager::dump(), which nothing calls, so it stays out of normal
# builds. Enable with MINIMAL_SURFACE_PROFILE=1.
if os.environ.get("MINIMAL_SURFACE_PROFILE", "").strip().lower() not in (
    "", "0", "false", "no", "off",
):
    define_macros.append(("PROFILE_FUNCTIONS", "1"))
    print("setup.py: function profiling enabled (PROFILE_FUNCTIONS)")

# delocate (the tool that bundles the ITK libraries into the macOS wheel) rewrites every library
# path to a longer one, and the extension links about a hundred of them. Without this padding
# the Mach-O header has no room for the longer paths and install_name_tool fails.
macos_link_args = ["-Wl,-headerpad_max_install_names"] if MACOS else []

extension = Extension(
    'minimal_surface',
    sources=["src/minimal-surface/code/_minimal_surface.pyx"] + glob.glob("src/minimal-surface/code/eikonal/*.cpp"),
    include_dirs=[
        np.get_include(),
        "src/minimal-surface/code",
        "src/minimal-surface/code/eikonal",
        "src/minimal-surface/code/pythonwrapper",
        "src/minimal-surface/code/eigen-3.4.0",
        *sdk["include_dirs"],
    ],
    depends=["MinimalSurfaceEstimator.h", "SimpleITK.h", "sitkImage.h"],
    define_macros=define_macros,
    library_dirs=sdk["library_dirs"],
    extra_objects=sdk["extra_objects"],
    extra_compile_args=extra_compile_args,
    extra_link_args=sdk["extra_link_args"] + omp_link_args + macos_link_args,
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
