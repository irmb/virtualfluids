PROJECT(DrivenCavityMultiGPU LANGUAGES CUDA CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator MPI::MPI_CXX FILES DrivenCavityMultiGPU.cpp)

set_source_files_properties(DrivenCavityMultiGPU.cpp PROPERTIES LANGUAGE CUDA)

set_target_properties(DrivenCavityMultiGPU PROPERTIES 
	CUDA_SEPARABLE_COMPILATION ON)