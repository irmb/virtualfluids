project(gpu_core LANGUAGES CUDA CXX)


vf_add_library(NAME gpu_core PUBLIC_LINK basics lbm parallel logger PRIVATE_LINK GridGenerator cuda_helper)

#SET(TPN_WIN32 "/EHsc")
#https://stackoverflow.com/questions/6832666/lnk2019-when-including-asio-headers-solution-generated-with-cmake
#https://stackoverflow.com/questions/27442885/syntax-error-with-stdnumeric-limitsmax

set_target_properties(gpu_core PROPERTIES CUDA_SEPARABLE_COMPILATION ON POSITION_INDEPENDENT_CODE ON)


if(MSVC)
    set_target_properties(gpu_core PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()

vf_add_tests(NAME gpu_core)

if(BUILD_VF_UNIT_TESTS)
    # set_target_properties(gpu_coreTests PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_include_directories(gpu_coreTests PRIVATE "${VF_THIRD_DIR}/cuda_samples/")
    target_include_directories(gpu_coreTests PRIVATE "${VF_ROOT_DIR}/src/gpu/GridGenerator/")
    set_source_files_properties(Communication/ExchangeData27Test.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(Factories/BoundaryConditionFactoryTest.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(DataStructureInitializer/GridReaderGenerator/GridGeneratorTest.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(DataStructureInitializer/GridReaderGenerator/IndexRearrangementForStreamsTest.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK17/CumulantK17Test.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(LBM/GPUHelperFunctions/KernelUtilitiesTests.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(Parameter/ParameterTest.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(PreCollisionInteractor/ActuatorFarmInlinesTest.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(GPU/GridScaling/interpolateStaticToRotatingInlinesTest.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(GPU/GridScaling/interpolateRotatingToStaticInlinesTest.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(LBM/GPUHelperFunctions/CoordinateTransformationTest.h PROPERTIES LANGUAGE CUDA)
endif()
