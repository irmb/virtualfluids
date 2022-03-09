PROJECT(SphereScaling LANGUAGES CUDA CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator MPI::MPI_CXX FILES SphereScaling.cpp)

set_source_files_properties(SphereScaling.cpp PROPERTIES LANGUAGE CUDA)

set_target_properties(SphereScaling PROPERTIES 
	CUDA_SEPARABLE_COMPILATION ON)
	# VS_DEBUGGER_COMMAND "C:/Program Files/Microsoft MPI/Bin/mpiexec.exe"
    # VS_DEBUGGER_COMMAND_ARGUMENTS "-n 2 \"$<TARGET_FILE:SphereScaling>\"")