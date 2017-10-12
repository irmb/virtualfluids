
find_path(JSONCPP_INCLUDE_DIRS json/json.h
    HINTS ${JSONCPP_ROOT} ENV JSONCPP_ROOT
    PATH_SUFFIXES include
    DOC "Directory where the JSONCPP header files are located"
)



find_library(JSONCPP_LIBRARIES
    NAMES jsoncpp.lib 
    HINTS ${JSONCPP_LIBRARIES_PATH}
    DOC "Directory where the JSONCPP library is located"
)


# Standard package handling
include(FindPackageHandleStandardArgs)
if(CMAKE_VERSION VERSION_GREATER 2.8.2)
  find_package_handle_standard_args(JSONCPP
    REQUIRED_VARS JSONCPP_INCLUDE_DIRS JSONCPP_LIBRARIES)
else()
  find_package_handle_standard_args(JSONCPP
    REQUIRED_VARS JSONCPP_INCLUDE_DIRS JSONCPP_LIBRARIES)
endif()


mark_as_advanced(JSONCPP_INCLUDE_DIRS JSONCPP_LIBRARIES)