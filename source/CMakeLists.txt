CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

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

PROJECT(VirtualFluids)

set(SOURCE_ROOT "../source")

#debug build for unix
 #IF(UNIX)
  #SET(CMAKE_BUILD_TYPE DEBUG)
 #ENDIF()

SET(USE_ZOLTAN OFF CACHE BOOL "include Zoltan library support")
SET(USE_METIS ON CACHE BOOL "include METIS library support")
SET(USE_MPI ON CACHE BOOL "include MPI library support")
SET(USE_VTK OFF CACHE BOOL "include VTK library support")
SET(USE_CATALYST OFF CACHE BOOL "include Paraview Catalyst support")
SET(USE_PYTHON OFF CACHE BOOL "include Python scripting support")
SET(USE_FETOL OFF CACHE BOOL "include FETOL library support")
SET(USE_INTEL OFF CACHE BOOL "include Intel compiler support")
SET(USE_GCC OFF CACHE BOOL "include gcc compiler support")

#CAB
INCLUDE("CMake/CMakeCABMacros.txt")

#MPI
FIND_PACKAGE(MPI REQUIRED)
#SET(MPI_CXX_LINK_FLAGS -mpe=mpilog)

SET(BOOST_USE_MULTITHREAD ON)
SET(Boost_USE_STATIC_LIBS ON)
#SET(Boost_DEBUG TRUE)

SET(bv ${BOOST_VERSION}) #hack for find boost, after next command ${BOOST_VERSION} would be set to 0
#FIND_PACKAGE(Boost ${bv} COMPONENTS system date_time thread serialization chrono regex)
#FIND_PACKAGE(Boost ${BOOST_VERSION} COMPONENTS system date_time thread serialization chrono regex)
#FIND_PACKAGE(Boost ${bv} COMPONENTS system thread serialization date_time)
#SET(BOOST_VERSION ${bv})
IF(${USE_PYTHON})
  FIND_PACKAGE(Boost ${BOOST_VERSION} COMPONENTS system date_time thread serialization chrono regex python)
ELSE(${USE_PYTHON})
    FIND_PACKAGE(Boost ${BOOST_VERSION} COMPONENTS system date_time thread serialization chrono regex)
ENDIF()

##################################################################################
#  Java
##############################################################################
### FindJNI.cmake
IF(${USE_FETOL})
    find_package(JNI REQUIRED) 
ENDIF()

#VTK
IF(${USE_VTK})
    #find_package(VTK 6.1 NO_MODULE)
    FIND_PACKAGE(VTK REQUIRED)
    INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
ENDIF()
MESSAGE("VTK_INCLUDE_DIRS = " ${VTK_INCLUDE_DIRS})

IF(${USE_CATALYST})
    find_package(ParaView 4.3 REQUIRED COMPONENTS vtkPVPythonCatalyst)
    include("${PARAVIEW_USE_FILE}")
ENDIF()

#LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DCAB_BOOST)
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DNOMINMAX)
#LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DBOOST_SIGNALS_NO_DEPRECATION_WARNING)
#LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DCAB_RUBY)
#LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -mpe=mpilog)
#LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -noshlib)
#LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DSINGLEPRECISION)

IF(${USE_ZOLTAN})
 LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DVF_ZOLTAN)
ENDIF() 
IF(${USE_METIS})
  LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DVF_METIS)
ENDIF()
IF(${USE_MPI})
  LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DVF_MPI)
ENDIF()
IF(${USE_FETOL})
  LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DVF_FETOL)
ENDIF()
IF(${USE_VTK})
  LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DVF_VTK)
ENDIF()
IF(${USE_CATALYST})
  LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DVF_CATALYST)
ENDIF()
IF(${USE_PYTHON})
  FIND_PACKAGE(PythonLibs)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIR})
  LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DVF_PYTHON)
  LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DBOOST_PYTHON_STATIC_LIB)
  add_subdirectory(python)
ENDIF()

IF(${USE_INTEL})
   LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DMPICH_IGNORE_CXX_SEEK)
   LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -DMPICH_SKIP_MPICXX)
ENDIF()
#message("MPI_CXX_LIBRARY: " ${MPI_CXX_LIBRARY})
#IF(MPI_CXX_LIBRARY)
   #SET(MPI_LIBRARY ${MPI_LIBRARY} ${MPI_CXX_LIBRARY})
   #message("MPI_LIBRARY: " ${MPI_LIBRARY})
#ENDIF() 

add_subdirectory(VirtualFluidsCore)

INCLUDE("Applications.cmake")

