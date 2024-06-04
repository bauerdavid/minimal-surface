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
        sitk_path = "C:/SimpleITK-build"
# extra_build_args = ["/Zc:externC"]
    else:
        sitk_path = os.path.expanduser("~/SimpleITK-build")
        os.environ["CC"] = "g++"


itk_path = os.path.join(sitk_path, "ITK-build")
itk_lib_path = os.path.join(itk_path, "lib")
itk_libs = [
    "libhdf5_hl_cpp-static",
    "libitkdouble-conversion-5.4",
    "libitksys-5.4",
    "libitkvcl-5.4",
    "libitkv3p_netlib-5.4",
    "libitktestlib-5.4",
    "libitkvnl-5.4",
    "libitkvnl_algo-5.4",
    "libITKVNLInstantiation-5.4",
    "libitkNetlibSlatec-5.4",
    "libITKTransform-5.4",
    "libITKFFT-5.4",
    "libITKMesh-5.4",
    "libitkzlib-5.4",
    "libITKMetaIO-5.4",
    "libITKSpatialObjects-5.4",
    "libITKPath-5.4",
    "libITKImageIntensity-5.4",
    "libITKConvolution-5.4",
    "libITKSmoothing-5.4",
    "libITKLabelMap-5.4",
    "libITKMathematicalMorphology-5.4",
    "libITKQuadEdgeMesh-5.4",
    "libITKFastMarching-5.4",
    "libITKImageFeature-5.4",
    "libITKOptimizers-5.4",
    "libITKPolynomials-5.4",
    "libITKBiasCorrection-5.4",
    "libITKColormap-5.4",
    "libITKDICOMParser-5.4",
    "libITKDeformableMesh-5.4",
    "libITKDenoising-5.4",
    "libITKDiffusionTensorImage-5.4",
    "libITKEXPAT-5.4",
    "libitkgdcmCommon-5.4",
    "libitkgdcmDICT-5.4",
    "libitkgdcmDSED-5.4",
    "libitkgdcmIOD-5.4",
    "libitkgdcmMEXD-5.4",
    "libitkgdcmMSFF-5.4",
    "libitkgdcmcharls-5.4",
    "libitkgdcmjpeg12-5.4",
    "libitkgdcmjpeg16-5.4",
    "libitkgdcmjpeg8-5.4",
    "libitkgdcmopenjp2-5.4",
    "libitkgdcmsocketxx-5.4",
    "libitkgdcmuuid-5.4",
    "libITKznz-5.4",
    "libITKniftiio-5.4",
    "libITKgiftiio-5.4",
    "libITKPDEDeformableRegistration-5.4",
    "libitkgtest-5.4",
    "libitkgtest_main-5.4",
    "libitkhdf5-static-5.4",
    "libitkhdf5_cpp-static-5.4",
    "libitkhdf5_hl-static-5.4",
    "libITKIOBMP-5.4",
    "libITKIOBioRad-5.4",
    "libITKIOBruker-5.4",
    "libITKIOCSV-5.4",
    "libITKIOGDCM-5.4",
    "libITKIOGE-5.4",
    "libITKIOGIPL-5.4",
    "libITKIOHDF5-5.4",
    "libitkjpeg-5.4",
    "libITKIOJPEG-5.4",
    "libitkopenjpeg-5.4",
    "libITKIOJPEG2000-5.4",
    "libitktiff-5.4",
    "libITKIOTIFF-5.4",
    "libITKIOLSM-5.4",
    "libitkminc2-5.4",
    "libITKIOMINC-5.4",
    "libITKIOMRC-5.4",
    "libITKIOMeshBase-5.4",
    "libITKIOMeshBYU-5.4",
    "libITKIOMeshFreeSurfer-5.4",
    "libITKIOMeshGifti-5.4",
    "libITKIOMeshOBJ-5.4",
    "libITKIOMeshOFF-5.4",
    "libITKIOMeshVTK-5.4",
    "libITKIOMeta-5.4",
    "libITKIONIFTI-5.4",
    "libITKNrrdIO-5.4",
    "libITKIONRRD-5.4",
    "libitkpng-5.4",
    "libITKIOPNG-5.4",
    "libITKIOSiemens-5.4",
    "libITKIOXML-5.4",
    "libITKIOSpatialObjects-5.4",
    "libITKIOStimulate-5.4",
    "libITKTransformFactory-5.4",
    "libITKIOTransformBase-5.4",
    "libITKIOTransformHDF5-5.4",
    "libITKIOTransformInsightLegacy-5.4",
    "libITKIOTransformMINC-5.4",
    "libITKIOTransformMatlab-5.4",
    "libITKIOVTK-5.4",
    "libITKKLMRegionGrowing-5.4",
    "libitklbfgs-5.4",
    "libITKMarkovRandomFieldsClassifiers-5.4",
    "libITKOptimizersv4-5.4",
    "libITKQuadEdgeMeshFiltering-5.4",
    "libITKRegionGrowing-5.4",
    "libITKRegistrationMethodsv4-5.4",
    "libITKVTK-5.4",
    "libITKWatersheds-5.4",
    "libITKReview-5.4",
    "libITKTestKernel-5.4",
    "libITKVideoCore-5.4",
    "libITKVideoIO-5.4",
    "libitkSimpleITKFilters-5.4",
    "libITKIOIPL-5.4",
    "libITKIOImageBase-5.4",
    "libITKCommon-5.4",
    "libITKStatistics-5.4",
]

itk_cxx_files = glob.glob(os.path.join(sitk_path, "ITK", "Modules", "**", "*.cxx"), recursive=True) 
# [
#     os.path.join("Core", "Common", "src", "itkRandomVariateGeneratorBase.cxx"),
#     os.path.join("Core", "Common", "src", "itkMersenneTwisterRandomVariateGenerator.cxx"),
#     os.path.join("IO", "HDF5", "src","itkHDF5ImageIO.cxx"),
#     os.path.join("ThirdParty","GDCM","src","gdcm","Source","DataStructureAndEncodingDefinition","gdcmByteValue.cxx"),
#     os.path.join("IO","NRRD","src","itkNrrdImageIO.cxx"),
#     os.path.join("IO","GDCM","src","itkGDCMImageIO.cxx"),
#     os.path.join("ThirdParty","VNL","src","vxl","core","vnl","Templates","vnl_c_vector+float-.cxx"),
#     os.path.join("Core","Common","src","itkDataObject.cxx"),
#     os.path.join("Core","Common","src","itkStreamingProcessObject.cxx"),
# ]
# itk_libs = list(map(lambda lib: f"{os.path.join(itk_lib_path, lib)}.so", itk_libs))
itk_libs = glob.glob(os.path.join(itk_lib_path, "*.so*"))+glob.glob(os.path.join(itk_lib_path, "*.lib"))
# for l1, l2 in zip(sorted(itk_libs_old), sorted(itk_libs)):
#     print(f"{l1}\t{l2}")
sitk_lib_path = os.path.join(sitk_path, "lib")
sitk_libs = glob.glob(os.path.join(sitk_lib_path, "*.*"))+glob.glob(os.path.join(sitk_lib_path, "*.lib"))
extra_compile_args = ["/std:c++17"] if sys.platform == "win32" else ["-std=c++17", "-fopenmp"]
extra_link_args = ["-lgomp"]
extra_objects = itk_libs+sitk_libs
libs = map(os.path.basename, extra_objects)
libs = map(lambda s: s.split(".so", 1)[0], libs)
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
    sources=[
        # os.path.join(sitk_path, "ITK", "Modules", itk_cxx_file) for itk_cxx_file in itk_cxx_files
    ]# + glob.glob(os.path.join("src", "minimal-surface", "code", "sitk_helper", "*.c*"))
     + ["src/minimal-surface/code/_minimal_surface.pyx"]
     + glob.glob("src/minimal-surface/code/eikonal/*.cpp"),
    include_dirs=[
        np.get_include(),
        os.path.join("src","minimal-surface","code"),
        glob.glob(os.path.join(sitk_path, "include", "SimpleITK-*"))[0],
        os.path.join(sitk_path, "ITK-prefix", "include", "ITK-5.4"),
        os.path.join(sitk_path, "ITK-prefix", "include", "ITK-5.4", "vnl")
    ]+ glob.glob(os.path.join("src", "minimal-surface", "code", "*", ""))
     + glob.glob(os.path.join(sitk_path, "ITK", "Modules", "**", "include/"), recursive=True),
    #   + glob.glob(os.path.join(sitk_path, "ITK-prefix", "include", "ITK-5.4", "*", "")),
    depends=["MinimalSurfaceEstimator.h"], #, "sitkInclude.h", "sitkImage.h"],
    library_dirs=[sitk_lib_path, itk_lib_path],
    libraries=libs,
    extra_objects=extra_objects,
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
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

