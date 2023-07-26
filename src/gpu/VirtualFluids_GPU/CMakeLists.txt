project(VirtualFluids_GPU LANGUAGES CUDA CXX)


vf_add_library(PUBLIC_LINK basics lbm parallel PRIVATE_LINK GridGenerator cuda_helper)

#SET(TPN_WIN32 "/EHsc")
#https://stackoverflow.com/questions/6832666/lnk2019-when-including-asio-headers-solution-generated-with-cmake
#https://stackoverflow.com/questions/27442885/syntax-error-with-stdnumeric-limitsmax

set_target_properties(VirtualFluids_GPU PROPERTIES CUDA_SEPARABLE_COMPILATION ON POSITION_INDEPENDENT_CODE ON)


if(MSVC)
    set_target_properties(VirtualFluids_GPU PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()

vf_add_tests()

if(BUILD_VF_UNIT_TESTS)
    # set_target_properties(VirtualFluids_GPUTests PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_include_directories(VirtualFluids_GPUTests PRIVATE "${VF_THIRD_DIR}/cuda_samples/")
    target_include_directories(VirtualFluids_GPUTests PRIVATE "${VF_ROOT_DIR}/src/gpu/GridGenerator/")
    set_source_files_properties(Communication/ExchangeData27Test.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(Factories/BoundaryConditionFactoryTest.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(DataStructureInitializer/GridReaderGenerator/GridGeneratorTest.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(DataStructureInitializer/GridReaderGenerator/IndexRearrangementForStreamsTest.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK17/CumulantK17Test.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(Kernel/Utilities/DistributionHelperTests.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(Parameter/ParameterTest.cpp PROPERTIES LANGUAGE CUDA)
endif()
