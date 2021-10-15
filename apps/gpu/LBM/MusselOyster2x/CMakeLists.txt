PROJECT(MusselOyster LANGUAGES CUDA CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics VirtualFluids_GPU GridGenerator MPI::MPI_CXX FILES MusselOyster2x.cpp)

set_source_files_properties(MusselOyster2x.cpp PROPERTIES LANGUAGE CUDA)

set_target_properties(MusselOyster2x PROPERTIES 
	CUDA_SEPARABLE_COMPILATION ON
	VS_DEBUGGER_COMMAND "C:/Program Files/Microsoft MPI/Bin/mpiexec.exe"
    VS_DEBUGGER_COMMAND_ARGUMENTS "-n 2 C:/Users/Master/Documents/MasterAnna/VirtualFluids_dev/build/bin/Release/MusselOyster2x.exe")