PROJECT(GyroidsRow)

vf_add_library(BUILDTYPE binary PRIVATE_LINK cpu_core muparser basics ${MPI_CXX_LIBRARIES} FILES GyroidsRow.cpp )

vf_get_library_name (library_name)
target_include_directories(${library_name} PRIVATE ${APPS_ROOT_CPU})

