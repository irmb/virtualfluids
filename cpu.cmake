#workaround for machine with mpi compiler wrapper
#it most define before project

#MPI
#set(CMAKE_C_COMPILER mpicc)
#set(CMAKE_CXX_COMPILER mpicxx)

#Intel MPI
#set(CMAKE_C_COMPILER mpiicc)
#set(CMAKE_CXX_COMPILER mpiicpc)

#Cray
#set(CMAKE_C_COMPILER cc)
#set(CMAKE_CXX_COMPILER CC)

#SuperMUC
#set(CMAKE_C_COMPILER mpicc)
#set(CMAKE_CXX_COMPILER mpiCC)

#debug build for unix
#IF(UNIX)
#SET(CMAKE_BUILD_TYPE DEBUG)
#ENDIF()

SET(VFCPU_USE_METIS ON CACHE BOOL "include METIS library support")
SET(VFCPU_USE_VTK OFF CACHE BOOL "include VTK library support")
SET(VFCPU_USE_CATALYST OFF CACHE BOOL "include Paraview Catalyst support")

SET(VFCPU_USE_HLRN_LUSTRE OFF CACHE BOOL "include HLRN Lustre support")
SET(VFCPU_USE_DEM_COUPLING OFF CACHE BOOL "PE plugin")

SET(VFCPU_ENABLE_LiggghtsCoupling ON CACHE BOOL "enable coupling with LIGGGHTS library")
SET(VFCPU_ENABLE_NonNewtonianFluids ON CACHE BOOL "enable non-Newtonian fluids module")
SET(VFCPU_ENABLE_MultiphaseFlow ON CACHE BOOL "enable multiphase flow module")

if(BUILD_VF_ALL_SAMPLES)
    set(VFCPU_ENABLE_NonNewtonianFluids ON)
    set(VFCPU_ENABLE_MultiphaseFlow ON)
endif() 

#MPI
IF((NOT ${CMAKE_CXX_COMPILER} MATCHES mpicxx) AND (NOT ${CMAKE_CXX_COMPILER} MATCHES mpiicpc))# OR NOT ${CMAKE_CXX_COMPILER} MATCHES cc OR NOT ${CMAKE_CXX_COMPILER} MATCHES mpiCC)
    FIND_PACKAGE(MPI REQUIRED)
ENDIF()
#SET(MPI_CXX_LINK_FLAGS -mpe=mpilog)

if(VFCPU_ENABLE_LiggghtsCoupling)
    add_subdirectory(src/cpu/LiggghtsCoupling)
    SET(VFCPU_USE_VTK ON CACHE BOOL "include VTK library support" FORCE)
endif()


#VTK
IF(${VFCPU_USE_VTK})
    FIND_PACKAGE(VTK REQUIRED)
    INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
ENDIF()

IF(${VFCPU_USE_CATALYST})
    find_package(ParaView 4.3 REQUIRED COMPONENTS vtkPVPythonCatalyst)
    include("${PARAVIEW_USE_FILE}")
ENDIF()

IF(${VFCPU_USE_METIS})
    list(APPEND VF_COMPILER_DEFINITION VF_METIS)
ENDIF()
IF(${VFCPU_USE_VTK})
    list(APPEND VF_COMPILER_DEFINITION VF_VTK)
ENDIF()
IF(${VFCPU_USE_CATALYST})
    list(APPEND VF_COMPILER_DEFINITION VF_CATALYST)
ENDIF()

IF(${VFCPU_USE_BOOST})
    list(APPEND VF_COMPILER_DEFINITION VF_BOOST)
ENDIF()

IF(${VFCPU_USE_HLRN_LUSTRE})
    list(APPEND VF_COMPILER_DEFINITION HLRN_LUSTRE)
ENDIF()

# workaround itanium processoren
IF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "ia64")
    LIST(APPEND VF_COMPILER_DEFINITION _M_IA64)
ENDIF()

if(${VFCPU_USE_METIS} AND NOT METIS_INCLUDEDIR)
    add_subdirectory(${VF_THIRD_DIR}/metis/metis-5.1.0)
endif()

add_subdirectory(${VF_THIRD_DIR}/MuParser)

add_subdirectory(src/cpu/VirtualFluidsCore)

if(BUILD_VF_PYTHON_BINDINGS)
    add_subdirectory(src/cpu/simulationconfig)
endif()

if(VFCPU_ENABLE_NonNewtonianFluids)
    add_subdirectory(src/cpu/NonNewtonianFluids)
endif()

if(VFCPU_ENABLE_MultiphaseFlow)
    add_subdirectory(src/cpu/MultiphaseFlow)
endif()

set (APPS_ROOT_CPU "${VF_ROOT_DIR}/apps/cpu/")
include(${APPS_ROOT_CPU}/Applications.cmake)