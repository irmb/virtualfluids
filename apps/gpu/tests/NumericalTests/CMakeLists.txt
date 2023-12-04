project(NumericalTestsGPU)
set(collectTestFiles ON)
vf_add_library(BUILDTYPE binary PRIVATE_LINK gpu_core GridGenerator basics logger fftw3 gtest gmock gmock_main cudart_static)

vf_get_library_name(library_name)
group_target(${library_name} ${testFolder})

target_include_directories(NumericalTests PRIVATE ${CMAKE_CURRENT_LIST_DIR})
