project(lbm LANGUAGES CXX)

vf_add_library(NAME lbm PUBLIC_LINK basics)

if(BUILD_VF_GPU)
    add_subdirectory(cuda)
endif()


vf_add_tests()
