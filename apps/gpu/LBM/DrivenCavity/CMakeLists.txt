PROJECT(DrivenCavity)

#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-DOMPI_SKIP_MPICXX" )

vf_add_library(BUILDTYPE binary PRIVATE_LINK VirtualFluids_GPU basics GridGenerator FILES DrivenCavity.cpp )

linkCUDA()
