project(NumericalTestsGPU)
set(collectTestFiles ON)

vf_add_library(BUILDTYPE binary PRIVATE_LINK gpu_core GridGenerator basics logger gtest gmock gmock_main cudart_static)

vf_get_library_name(library_name)
group_target(${library_name} ${testFolder})

target_include_directories(NumericalTests PRIVATE ${CMAKE_CURRENT_LIST_DIR})

if (UNIX)
    find_package(PkgConfig REQUIRED)
    pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)

    target_link_libraries(${library_name} PRIVATE PkgConfig::FFTW)
else()
    message(FATAL_ERROR "Currently the numerical tests can only be compiled on unix systems.")
endif()
