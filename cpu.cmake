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

PROJECT(VirtualFluids)

#debug build for unix
#IF(UNIX)
#SET(CMAKE_BUILD_TYPE DEBUG)
#ENDIF()

SET(USE_METIS ON CACHE BOOL "include METIS library support")
SET(USE_MPI ON CACHE BOOL "include MPI library support")
SET(USE_VTK OFF CACHE BOOL "include VTK library support")
SET(USE_CATALYST OFF CACHE BOOL "include Paraview Catalyst support")
SET(USE_BOOST OFF CACHE BOOL "include Boost support")
#SET(USE_PYTHON OFF CACHE BOOL "include Python scripting support")

SET(USE_HLRN_LUSTRE OFF CACHE BOOL "include HLRN Lustre support")
SET(USE_DEM_COUPLING OFF CACHE BOOL "PE plugin")

#MPI
IF((NOT ${CMAKE_CXX_COMPILER} MATCHES mpicxx) AND (NOT ${CMAKE_CXX_COMPILER} MATCHES mpiicpc))# OR NOT ${CMAKE_CXX_COMPILER} MATCHES cc OR NOT ${CMAKE_CXX_COMPILER} MATCHES mpiCC)
    FIND_PACKAGE(MPI REQUIRED)
ENDIF()
#SET(MPI_CXX_LINK_FLAGS -mpe=mpilog)

#SET(BOOST_USE_MULTITHREAD ON)
#SET(Boost_USE_STATIC_LIBS ON)
#SET(Boost_DEBUG TRUE)

#SET(bv ${BOOST_VERSION}) #hack for find boost, after next command ${BOOST_VERSION} would be set to 0
#FIND_PACKAGE(Boost ${bv} COMPONENTS system date_time thread serialization chrono regex)
#FIND_PACKAGE(Boost ${BOOST_VERSION} COMPONENTS system date_time thread serialization chrono regex)
#FIND_PACKAGE(Boost ${bv} COMPONENTS system thread serialization date_time)
#SET(BOOST_VERSION ${bv})
#IF(${USE_PYTHON})
#  FIND_PACKAGE(Boost ${BOOST_VERSION} COMPONENTS system date_time thread serialization chrono regex python)
#ELSE(${USE_PYTHON})
#    FIND_PACKAGE(Boost ${BOOST_VERSION} COMPONENTS system date_time thread serialization chrono regex)
#ENDIF()

IF(${USE_BOOST})
    FIND_PACKAGE(Boost ${BOOST_VERSION})
ENDIF()

#VTK
IF(${USE_VTK})
    #find_package(VTK 6.1 NO_MODULE)
    FIND_PACKAGE(VTK REQUIRED)
    INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
ENDIF()

IF(${USE_CATALYST})
    find_package(ParaView 4.3 REQUIRED COMPONENTS vtkPVPythonCatalyst)
    include("${PARAVIEW_USE_FILE}")
ENDIF()

IF(${USE_METIS})
    list(APPEND VF_COMPILER_DEFINITION VF_METIS)
ENDIF()
IF(${USE_MPI})
    list(APPEND VF_COMPILER_DEFINITION VF_MPI)
ENDIF()
IF(${USE_VTK})
    list(APPEND VF_COMPILER_DEFINITION VF_VTK)
ENDIF()
IF(${USE_CATALYST})
    list(APPEND VF_COMPILER_DEFINITION VF_CATALYST)
ENDIF()

IF(${USE_BOOST})
    list(APPEND VF_COMPILER_DEFINITION VF_BOOST)
ENDIF()

IF(${USE_HLRN_LUSTRE})
    list(APPEND VF_COMPILER_DEFINITION HLRN_LUSTRE)
ENDIF()

IF(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    list(APPEND VF_LINK_OPTIONS -parallel)
    list(APPEND VF_LINK_OPTIONS -irc)
ENDIF()

IF(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    list(APPEND VF_LINK_OPTIONS -lgomp)
    list(APPEND VF_LINK_OPTIONS -lrt)
ENDIF()



# IF(${USE_PYTHON})
# FIND_PACKAGE(PythonLibs)
# INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIR})
# LIST(APPEND VF_COMPILER_DEFINITION VF_PYTHON)
# LIST(APPEND VF_COMPILER_DEFINITION BOOST_PYTHON_STATIC_LIB)
# add_subdirectory(python)
# ENDIF()

# IF(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
# LIST(APPEND VF_COMPILER_DEFINITION MPICH_IGNORE_CXX_SEEK)
# LIST(APPEND VF_COMPILER_DEFINITION MPICH_SKIP_MPICXX)
# ENDIF()


add_subdirectory(${VF_THIRD_DIR}/MuParser)

add_subdirectory(src/cpu/VirtualFluidsCore)

set (APPS_ROOT_CPU "${VF_ROOT_DIR}/apps/cpu/")
include(${APPS_ROOT_CPU}/Applications.cmake)