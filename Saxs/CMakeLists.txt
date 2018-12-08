# Policy
cmake_minimum_required( VERSION 3.0) 


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCE SaxsPadding.cpp BackupFile.cpp LagrangeInterpolation.cpp 
	RhoSaxsBSP.cpp RhoSaxs.cpp RhoSaxsLI.cpp SaxsBSP.cpp Saxs.cpp
	SaxsHistogram.cpp SaxsHistogramSpline.cpp ScatteringFactors.cpp 
	ScatteringFactorsN.cpp SaxsDebye.cpp SaxsDir.cpp Spline1DInterpolant.cpp
	SaxsBSPStatic.cpp SaxsBSPfixed.cpp SaxsFilter.cpp fftw++.cpp)
add_library(saxs STATIC ${SOURCE})
target_include_directories(saxs PUBLIC "${PROJECT_SOURCE_DIR}/Properties" "${PROJECT_SOURCE_DIR}/Molecules" "${PROJECT_SOURCE_DIR}/libtraj")