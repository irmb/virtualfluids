
set(VF_LIBRARIES)

IF(METIS_RELEASE_LIBRARY AND METIS_DEBUG_LIBRARY)
   list(APPEND VF_LIBRARIES optimized ${METIS_RELEASE_LIBRARY} debug ${METIS_DEBUG_LIBRARY})
ELSE()
   list(APPEND VF_LIBRARIES metis)
ENDIF()

IF(${VF_CPU_ENABLE_VTK})
   list(APPEND VF_LIBRARIES optimized ${VTK_LIBRARIES} debug ${VTK_LIBRARIES})
ENDIF()

IF(${VF_CP_ENABLE_CATALYST})
   list(APPEND VF_LIBRARIES optimized vtkPVPythonCatalyst debug vtkPVPythonCatalyst)
   list(APPEND VF_LIBRARIES optimized vtkParallelMPI debug vtkParallelMPI )
ENDIF()

if(VF_ENABLE_OPENMP)
    list(APPEND VF_LIBRARIES OpenMP::OpenMP_CXX)
endif()

set(library_name cpu_core)

vf_add_library(NAME ${library_name} BUILDTYPE static PUBLIC_LINK basics muparser ${VF_LIBRARIES} parallel project_options PRIVATE_LINK lbm logger project_warnings)

vf_add_tests(NAME ${library_name})

target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryConditions)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Connectors)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Data)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Interactors)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/LBM)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/LBM/Interpolation)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Parallel)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Simulation)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Visitors)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/SimulationObservers)
target_include_directories(${library_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Utilities)


IF(${VF_CPU_ENABLE_VTK})
   target_include_directories(${library_name} PRIVATE ${VTK_INCLUDE_DIRS})
ENDIF()

# TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/139
# if(VF_ENABLE_MPI)
    target_link_libraries(${library_name} PRIVATE MPI::MPI_CXX)
# endif()

## WARNINGS
# TODO: Remove this and fix those warnings
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    # Using CMakes De-Duplication of Flags "Hack" (https://cmake.org/cmake/help/latest/prop_tgt/COMPILE_OPTIONS.html#option-de-duplication)
    # target_compile_options(${library_name} PUBLIC "-Wold-style-cast;-Wno-old-style-cast")
    # target_compile_options(${library_name} PUBLIC "-Wsign-conversion;-Wno-sign-conversion")
    # target_compile_options(${library_name} PUBLIC "-Wdouble-promotion;-Wno-double-promotion")
    # target_compile_options(${library_name} PUBLIC "-Wshadow;-Wno-shadow")

    target_compile_options(${library_name} PUBLIC "-Wextra") # Wno-unused-parameter comes with -Wextra
    target_compile_options(${library_name} PUBLIC "-Wall") # Wno-unused-function comes with -Wall
    target_compile_options(${library_name} PUBLIC "-Wunused") # -Wno-unused-function comes with -Wunused
    target_compile_options(${library_name} PUBLIC "-Wno-unused-parameter")
    target_compile_options(${library_name} PUBLIC "-Wno-unused-function")
    # target_compile_options(${library_name} PRIVATE "-Wno-implicit-float-conversion")
    # target_compile_options(${library_name} PRIVATE "-Wno-comment")
    # target_compile_options(${library_name} PRIVATE "-Wno-reorder-ctor")
    # target_compile_options(${library_name} PRIVATE "-Wno-uninitialized")
    target_compile_options(${library_name} PUBLIC "-Wno-reorder")
    # target_compile_options(${library_name} PRIVATE "-Wno-strict-aliasing")
endif()

# if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
#     target_compile_options(${library_name} PUBLIC "-Wuseless-cast;-Wno-useless-cast")
# endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(${library_name} PUBLIC "/wd4244") # conversion
    target_compile_options(${library_name} PUBLIC "/wd4100") # unreferenced formal parameter
    target_compile_options(${library_name} PUBLIC "/wd4701") # potentially uninitialized local variable
    target_compile_options(${library_name} PUBLIC "/wd4703") # potentially uninitialized local pointer
    target_compile_options(${library_name} PUBLIC "/wd4458") # declaration of 'XXX' hides class member
    target_compile_options(${library_name} PUBLIC "/wd4305") #'initializing': truncation from double to float
    target_compile_options(${library_name} PUBLIC "/wd4505") # C4505: 'XXX': unreferenced local function has been removed
    target_compile_options(${library_name} PUBLIC "/wd4310") # C4310: cast truncates constant value, triggerd by muParserbase.h
    target_compile_options(${library_name} PUBLIC "/wd4005") # disable macro redefinition (triggered by metis.h)
    target_compile_options(${library_name} PUBLIC "/wd4456") # C4456: declaration of hides previous local declaration
endif()