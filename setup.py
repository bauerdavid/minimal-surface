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
    else:
        sitk_path = "/SimpleITK-build"


extension = Extension(
    'MinArea',
    sources=[r"src\minimal-surface\code\_minimal_surface.pyx"] + glob.glob("src\minimal-surface\code\eikonal\*.cpp"),
    include_dirs=[
        np.get_include(),
        r"src\minimal-surface\code",
        r"src\minimal-surface\code\eikonal",
        r"src\minimal-surface\code\pythonwrapper",
        r"src\minimal-surface\code\eigen-3.4.0",
        os.path.join(sitk_path, r"include\SimpleITK-2.2"),
    ],
    depends=["MinimalSurfaceEstimator.h", "SimpleITK.h", "sitkImage.h"],
    library_dirs=[os.path.join(sitk_path, r"lib"), os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel")],
    extra_objects=[
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkdouble-conversion-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itksys-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkvnl_algo-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkvnl-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkv3p_netlib-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkvcl-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKCommon-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkNetlibSlatec-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKStatistics-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKTransform-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKMesh-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkzlib-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKMetaIO-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKSpatialObjects-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKPath-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKLabelMap-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKMathematicalMorphology-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKQuadEdgeMesh-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKFastMarching-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOImageBase-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKSmoothing-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKImageFeature-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKOptimizers-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKPolynomials-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKBiasCorrection-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKColormap-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKFFT-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKConvolution-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKDICOMParser-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKDeformableMesh-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKDenoising-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKDiffusionTensorImage-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKEXPAT-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkgdcmDICT-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkgdcmMSFF-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKznz-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKniftiio-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKgiftiio-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKPDEDeformableRegistration-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkhdf5_cpp-static-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkhdf5-static-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOBMP-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOBioRad-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOBruker-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOCSV-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOGDCM-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOIPL-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOGE-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOGIPL-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOHDF5-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkjpeg-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOJPEG-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkopenjpeg-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOJPEG2000-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itktiff-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOTIFF-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOLSM-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkminc2-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOMINC-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOMRC-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOMeshBase-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOMeshBYU-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOMeshFreeSurfer-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOMeshGifti-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOMeshOBJ-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOMeshOFF-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOMeshVTK-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOMeta-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIONIFTI-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKNrrdIO-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIONRRD-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkpng-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOPNG-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOSiemens-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOXML-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOSpatialObjects-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOStimulate-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKTransformFactory-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOTransformBase-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOTransformHDF5-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOTransformInsightLegacy-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOTransformMINC-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOTransformMatlab-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKIOVTK-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKKLMRegionGrowing-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itklbfgs-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKMarkovRandomFieldsClassifiers-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKOptimizersv4-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKQuadEdgeMeshFiltering-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKRegionGrowing-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKRegistrationMethodsv4-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKVTK-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKWatersheds-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKReview-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKTestKernel-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKVideoCore-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKVideoIO-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkSimpleITKFilters-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkgdcmIOD-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkgdcmDSED-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkgdcmCommon-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkgdcmjpeg8-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkgdcmjpeg12-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkgdcmjpeg16-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkgdcmopenjp2-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\itkgdcmcharls-5.3.lib"),
        os.path.join(sitk_path, r"ITK-build\lib\MinSizeRel\ITKVNLInstantiation-5.3.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKAnisotropicSmoothing-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKAntiAlias-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKBiasCorrection-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKBinaryMathematicalMorphology-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKClassifiers-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKColormap-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKCommon-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKConnectedComponents-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKConvolution-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKCurvatureFlow-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKDeconvolution-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKDenoising-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKDisplacementField-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKDistanceMap-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKFastMarching-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKFFT-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKImageCompare-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKImageCompose-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKImageFeature-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKImageFilterBase-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKImageFunction-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKImageFusion-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKImageGradient-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKImageGrid-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKImageIntensity-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKImageLabel-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKImageNoise-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKImageSources-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKImageStatistics-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKLabelMap-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKLabelVoting-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKLevelSets-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKMathematicalMorphology-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKPDEDeformableRegistration-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKRegionGrowing-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKRegistrationCommon-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKReview-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKSmoothing-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKSuperPixel-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKThresholding-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKTransform-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_ITKWatersheds-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITK_SimpleITKFilters-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITKBasicFilters0-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITKBasicFilters1-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITKCommon-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITKIO-2.2.lib"),
        os.path.join(sitk_path, r"lib\SimpleITKRegistration-2.2.lib"),
        r"crypt32.lib",
        r"rpcrt4.lib",
        r"comctl32.lib",
        r"wsock32.lib",
        r"ws2_32.lib",
        r"dbghelp.lib",
        r"psapi.lib",
        r"kernel32.lib",
        r"user32.lib",
        r"gdi32.lib",
        r"winspool.lib",
        r"shell32.lib",
        r"ole32.lib",
        r"oleaut32.lib",
        r"uuid.lib",
        r"comdlg32.lib",
        r"advapi32.lib",
    ],
    extra_compile_args=["/std:c++17"],
    language="c++"
)

CYTHONIZE = cythonize is not None

if CYTHONIZE:
    compiler_directives = {"language_level": 3, "embedsignature": True}
    extensions = cythonize([extension], compiler_directives=compiler_directives)
else:
    extensions = no_cythonize([extension])

setup(name='MinArea',
    description='Python Package with minimalsurface C++ extension',
    ext_modules=extensions
)

