

IF(MSVC)
    set(libsToLink ws2_32 GridGenerator basics) # ws_32 throws an error on Phoenix
ELSE(MSVC)
    set(libsToLink GridGenerator basics)
ENDIF(MSVC)


vf_add_library(BUILDTYPE shared DEPENDS ${libsToLink})

include (${VF_CMAKE_DIR}/3rd/cuda.cmake)

vf_get_library_name(library_name)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/basics)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/gpu)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/gpu/GridGenerator)
target_include_directories(${library_name} PRIVATE ${VF_SRC_DIR}/basics/Core)


#SET(TPN_WIN32 "/EHsc")
#https://stackoverflow.com/questions/6832666/lnk2019-when-including-asio-headers-solution-generated-with-cmake
#https://stackoverflow.com/questions/27442885/syntax-error-with-stdnumeric-limitsmax

IF(MSVC)
    set_target_properties(${library_name} PROPERTIES LINK_FLAGS "/ignore:4251")
ENDIF(MSVC)
