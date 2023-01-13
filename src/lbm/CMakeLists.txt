if(BUILD_VF_CPU)
    project(lbm LANGUAGES CXX)

    vf_add_library(NAME lbm PUBLIC_LINK basics)
    target_link_libraries(lbm PRIVATE project_warnings)

    vf_add_tests()
endif()

if(BUILD_VF_GPU OR BUILD_VF_GKS)
    add_subdirectory(cuda)
endif()
