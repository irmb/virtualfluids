project(VirtualFluids_GPU LANGUAGES CUDA CXX)

set(additional_libraries "")
if(MSVC)
    set(additional_libraries ws2_32 Traffic) # ws_32 throws an error on Phoenix
endif()

vf_add_library(PUBLIC_LINK basics lbmCuda PRIVATE_LINK ${additional_libraries} GridGenerator MPI::MPI_CXX vf_cuda)

#SET(TPN_WIN32 "/EHsc")
#https://stackoverflow.com/questions/6832666/lnk2019-when-including-asio-headers-solution-generated-with-cmake
#https://stackoverflow.com/questions/27442885/syntax-error-with-stdnumeric-limitsmax

set_target_properties(VirtualFluids_GPU PROPERTIES CUDA_SEPARABLE_COMPILATION ON POSITION_INDEPENDENT_CODE ON)


vf_add_tests()

if(BUILD_VF_UNIT_TESTS)
    set_target_properties(VirtualFluids_GPUTests PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    set_source_files_properties(Kernel/Utilities/DistributionHelperTests.cpp PROPERTIES LANGUAGE CUDA)
	set_source_files_properties(DataStructureInitializer/GridReaderGenerator/IndexRearrangementForStreamsTest.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(Communication/ExchangeData27Test.cpp PROPERTIES LANGUAGE CUDA)
    target_include_directories(VirtualFluids_GPUTests PRIVATE "${VF_THIRD_DIR}/cuda_samples/")
endif()
