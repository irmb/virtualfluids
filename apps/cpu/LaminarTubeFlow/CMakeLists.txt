
PROJECT(ltf)

vf_add_library(BUILDTYPE binary DEPENDS VirtualFluidsCore muparser basics ${MPI_CXX_LIBRARIES} FILES ltf.cpp )

vf_get_library_name (library_name)
target_include_directories(${library_name} PRIVATE ${APPS_ROOT_CPU})