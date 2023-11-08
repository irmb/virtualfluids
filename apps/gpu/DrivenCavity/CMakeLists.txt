PROJECT(DrivenCavity LANGUAGES CXX)

#LIST(APPEND CS_COMPILER_FLAGS_CXX "-DOMPI_SKIP_MPICXX" )

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics gpu_core GridGenerator FILES DrivenCavity.cpp)
