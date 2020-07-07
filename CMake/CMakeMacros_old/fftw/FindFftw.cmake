
find_path(FFTW_INCLUDE_DIRS fftw3.h
    HINTS ${FFTW_ROOT} ENV FFTW_ROOT
    PATH_SUFFIXES include
    DOC "Directory where the FFTW header files are located"
)



find_library(FFTW_LIBRARIES
    NAMES libfftw3-3.lib libfftw3f-3.lib libfftw3l-3.lib
    HINTS ${FFTW_LIBRARIES_PATH}
    DOC "Directory where the FFTW library is located"
)


# Standard package handling
include(FindPackageHandleStandardArgs)
if(CMAKE_VERSION VERSION_GREATER 2.8.2)
  find_package_handle_standard_args(FFTW
    REQUIRED_VARS FFTW_INCLUDE_DIRS FFTW_LIBRARIES)
else()
  find_package_handle_standard_args(FFTW
    REQUIRED_VARS FFTW_INCLUDE_DIRS FFTW_LIBRARIES)
endif()


mark_as_advanced(FFTW_INCLUDE_DIRS FFTW_LIBRARIES)

