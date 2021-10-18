PROJECT(MusselOyster LANGUAGES CUDA CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator MPI::MPI_CXX FILES MusselOyster.cpp)

set_source_files_properties(MusselOyster.cpp PROPERTIES LANGUAGE CUDA)

set_target_properties(MusselOyster PROPERTIES 
	CUDA_SEPARABLE_COMPILATION ON
	VS_DEBUGGER_COMMAND "C:/Program Files/Microsoft MPI/Bin/mpiexec.exe"
    VS_DEBUGGER_COMMAND_ARGUMENTS "-n 2 C:/Users/Master/Documents/MasterAnna/VirtualFluids_dev/build/bin/Debug/MusselOyster.exe")