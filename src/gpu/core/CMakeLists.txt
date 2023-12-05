project(gpu_core LANGUAGES CUDA CXX)


vf_add_library(NAME gpu_core PUBLIC_LINK basics lbm parallel logger project_options project_warnings PRIVATE_LINK GridGenerator cuda_helper)

#SET(TPN_WIN32 "/EHsc")
#https://stackoverflow.com/questions/6832666/lnk2019-when-including-asio-headers-solution-generated-with-cmake
#https://stackoverflow.com/questions/27442885/syntax-error-with-stdnumeric-limitsmax

set_target_properties(gpu_core PROPERTIES CUDA_SEPARABLE_COMPILATION ON POSITION_INDEPENDENT_CODE ON)

target_include_directories(gpu_core SYSTEM PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} "${VF_THIRD_DIR}/cuda_samples/")

if(MSVC)
    set_target_properties(gpu_core PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()

vf_add_tests(NAME gpu_core)

if(VF_ENABLE_UNIT_TESTS)
    # set_target_properties(gpu_coreTests PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_include_directories(gpu_coreTests PRIVATE "${VF_THIRD_DIR}/cuda_samples/")
    target_include_directories(gpu_coreTests PRIVATE "${VF_ROOT_DIR}/src/gpu/GridGenerator/")
    set_source_files_properties(Communication/ExchangeData27Test.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(BoundaryConditions/BoundaryConditionFactoryTest.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(DataStructureInitializer/GridReaderGenerator/GridGeneratorTest.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(DataStructureInitializer/GridReaderGenerator/IndexRearrangementForStreamsTest.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK17/CumulantK17Test.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(Utilities/KernelUtilitiesTests.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(Output/DistributionDebugWriterTest.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(Parameter/ParameterTest.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(PreCollisionInteractor/Actuator/ActuatorFarmInlinesTest.cpp PROPERTIES LANGUAGE CUDA)
endif()


## WARNINGS
# TODO: Remove this and fix those warnings: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/163
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    # Using CMakes De-Duplication of Flags "Hack" (https://cmake.org/cmake/help/latest/prop_tgt/COMPILE_OPTIONS.html#option-de-duplication)
    list(APPEND WARNINGS "-Wold-style-cast;-Wno-old-style-cast")
    list(APPEND WARNINGS "-Wsign-conversion;-Wno-sign-conversion")
    list(APPEND WARNINGS "-Wdouble-promotion;-Wno-double-promotion")
    list(APPEND WARNINGS "-Wshadow;-Wno-shadow")

    list(APPEND WARNINGS "-Wextra") # Wno-unused-parameter comes with -Wextra
    list(APPEND WARNINGS "-Wall") # Wno-unused-function comes with -Wall
    list(APPEND WARNINGS "-Wunused") # -Wno-unused-function comes with -Wunused
    list(APPEND WARNINGS "-Wno-unused-parameter")
    list(APPEND WARNINGS "-Wno-unused-function")
    list(APPEND WARNINGS "-Wno-implicit-float-conversion")
    list(APPEND WARNINGS "-Wno-comment")
    list(APPEND WARNINGS "-Wno-reorder-ctor")
    list(APPEND WARNINGS "-Wno-uninitialized")
    list(APPEND WARNINGS "-Wno-reorder")
    list(APPEND WARNINGS "-Wno-strict-aliasing")

    if(VF_ENABLE_UNIT_TESTS)
        target_compile_options(gpu_coreTests PUBLIC "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-Wnull-dereference;-Xcompiler=-Wno-null-dereference>")
        target_compile_options(gpu_coreTests PUBLIC "$<$<COMPILE_LANGUAGE:CXX>:-Wnull-dereference;-Wno-null-dereference>")
    endif()
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    list(APPEND WARNINGS "-Wuseless-cast;-Wno-useless-cast")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    list(APPEND WARNINGS "/wd4244") # conversion
    list(APPEND WARNINGS "/wd4100") # unreferenced formal parameter
    list(APPEND WARNINGS "/wd4701") # potentially uninitialized local variable
    list(APPEND WARNINGS "/wd4703") # potentially uninitialized local pointer
    list(APPEND WARNINGS "/wd4458") # declaration of 'XXX' hides class member
    list(APPEND WARNINGS "/wd4305") # 'initializing': truncation from double to float  
    list(APPEND WARNINGS "/wd4505") # 'XXX': unreferenced local function has been removed
    list(APPEND WARNINGS "/wd4324") # structure was padded due to alignment specifier
    list(APPEND WARNINGS "/w14265;/wd4265")
    list(APPEND WARNINGS "/wd4389")
    list(APPEND WARNINGS "/w14296;/wd4296")

    target_link_options(gpu_core PUBLIC $<DEVICE_LINK:-Xcompiler=/wd4100>)
endif()

foreach(flag IN LISTS WARNINGS)
    target_compile_options(gpu_core PUBLIC "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${flag}>")
    target_compile_options(gpu_core PUBLIC "$<$<COMPILE_LANGUAGE:CXX>:${flag}>")
endforeach()
