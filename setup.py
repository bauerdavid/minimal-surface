from distutils.core import setup, Extension, DEBUG
import glob

sfc_module = Extension('MinArea', sources = glob.glob("*.cpp"),
                      include_dirs=[r"Y:\BIOMAG\single-cell-raman\src\venv\Lib\site-packages\numpy\core\include",
                                    r"Y:\BIOMAG\shortest path\MinArea\MinArea",
                                    r"C:\SimpleITK-build\include\SimpleITK-2.2",
                                    r"Y:\BIOMAG\shortest path\MinArea\packages\Eigen3.3.3.9\lib\native\include"],
                      depends=["MinimalSurfaceEstimator.h", "SimpleITK.h", "sitkImage.h"],
                      library_dirs=[r"C:\SimpleITK-build\lib", r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel"],
                      extra_objects=[r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkdouble-conversion-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itksys-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkvnl_algo-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkvnl-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkv3p_netlib-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkvcl-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKCommon-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkNetlibSlatec-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKStatistics-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKTransform-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKMesh-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkzlib-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKMetaIO-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKSpatialObjects-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKPath-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKLabelMap-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKMathematicalMorphology-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKQuadEdgeMesh-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKFastMarching-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOImageBase-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKSmoothing-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKImageFeature-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKOptimizers-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKPolynomials-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKBiasCorrection-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKColormap-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKFFT-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKConvolution-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKDICOMParser-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKDeformableMesh-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKDenoising-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKDiffusionTensorImage-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKEXPAT-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkgdcmDICT-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkgdcmMSFF-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKznz-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKniftiio-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKgiftiio-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKPDEDeformableRegistration-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkhdf5_cpp-static-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkhdf5-static-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOBMP-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOBioRad-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOBruker-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOCSV-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOGDCM-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOIPL-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOGE-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOGIPL-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOHDF5-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkjpeg-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOJPEG-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkopenjpeg-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOJPEG2000-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itktiff-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOTIFF-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOLSM-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkminc2-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOMINC-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOMRC-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOMeshBase-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOMeshBYU-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOMeshFreeSurfer-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOMeshGifti-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOMeshOBJ-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOMeshOFF-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOMeshVTK-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOMeta-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIONIFTI-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKNrrdIO-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIONRRD-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkpng-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOPNG-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOSiemens-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOXML-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOSpatialObjects-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOStimulate-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKTransformFactory-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOTransformBase-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOTransformHDF5-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOTransformInsightLegacy-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOTransformMINC-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOTransformMatlab-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKIOVTK-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKKLMRegionGrowing-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itklbfgs-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKMarkovRandomFieldsClassifiers-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKOptimizersv4-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKQuadEdgeMeshFiltering-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKRegionGrowing-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKRegistrationMethodsv4-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKVTK-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKWatersheds-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKReview-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKTestKernel-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKVideoCore-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKVideoIO-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkSimpleITKFilters-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkgdcmIOD-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkgdcmDSED-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkgdcmCommon-5.3.lib",
                                    r"crypt32.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkgdcmjpeg8-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkgdcmjpeg12-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkgdcmjpeg16-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkgdcmopenjp2-5.3.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\itkgdcmcharls-5.3.lib",
                                    r"rpcrt4.lib",
                                    r"comctl32.lib",
                                    r"wsock32.lib",
                                    r"ws2_32.lib",
                                    r"dbghelp.lib",
                                    r"psapi.lib",
                                    r"C:\SimpleITK-build\ITK-build\lib\MinSizeRel\ITKVNLInstantiation-5.3.lib",
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
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKAnisotropicSmoothing-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKAntiAlias-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKBiasCorrection-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKBinaryMathematicalMorphology-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKClassifiers-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKColormap-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKCommon-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKConnectedComponents-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKConvolution-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKCurvatureFlow-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKDeconvolution-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKDenoising-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKDisplacementField-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKDistanceMap-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKFastMarching-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKFFT-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKImageCompare-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKImageCompose-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKImageFeature-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKImageFilterBase-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKImageFunction-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKImageFusion-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKImageGradient-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKImageGrid-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKImageIntensity-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKImageLabel-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKImageNoise-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKImageSources-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKImageStatistics-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKLabelMap-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKLabelVoting-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKLevelSets-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKMathematicalMorphology-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKPDEDeformableRegistration-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKRegionGrowing-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKRegistrationCommon-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKReview-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKSmoothing-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKSuperPixel-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKThresholding-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKTransform-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_ITKWatersheds-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITK_SimpleITKFilters-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITKBasicFilters0-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITKBasicFilters1-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITKCommon-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITKIO-2.2.lib",
                                    r"C:\SimpleITK-build\lib\SimpleITKRegistration-2.2.lib",],
                      extra_compile_args = ["/std:c++17"])

setup(name = 'MinArea', version = '0.1.9',
    description = 'Python Package with minimalsurface C++ extension',
    ext_modules = [sfc_module]
    )