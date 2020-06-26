#LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__unix__)
#LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__UNIX__)
#LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__APPLE__)

IF(NOT CMAKE_CXX_COMPILER)
   MESSAGE(FATAL_ERROR "before cmake-config-file can be included -> project must be extecuted")
ENDIF()
  
#################################################################################
# MPI
#################################################################################
#SET(MPI_DIR  "C:/Program Files (x86)/mpich2")
#SET(MPI_DIR  "C:/Program Files/mpich2")
#SET(USE_MPI_CXX_SYNTAX OFF)
#SET(MPI_COMPILER "C:/Program Files/mpich2/bin/mpicxx")
#SET(MPI_INCLUDE_PATH "C:/Program Files (x86)/mpich2/include")
#SET(MPI_LIBRARY "C:/Program Files/mpich2/lib/libmpi.a")
#SET(MPI_CXX_LIBRARY "C:/Program Files/MPICH2/lib/cxx.lib")
#################################################################################
#  BOOST  
#################################################################################
SET(BOOST_VERSION "1.60.0")
SET(BOOST_ROOT "c:/Tools/boost/boost_1_60_0")
SET(BOOST_DIR ${BOOST_ROOT})
SET(BOOST_LIBRARYDIR ${BOOST_ROOT}"/stageMSVC64/lib")  
#################################################################################
#  VTK  
#################################################################################
set(VTK_DIR "E:/Tools/VTK-8.1.2-build")
#################################################################################
#  ZOLTAN  
#################################################################################
IF(${USE_ZOLTAN})
  SET(ZOLTAN_INCLUDEDIR "d:/Tools/zoltan/include")
  SET(ZOLTAN_DEBUG_LIBRARY "d:/Tools/zoltan/lib/Debug/zoltan.lib") 
  SET(ZOLTAN_RELEASE_LIBRARY "d:/Tools/zoltan/lib/Release/zoltan.lib") 
ENDIF()

#################################################################################
#  METIS  
#################################################################################
IF(${USE_METIS})
  SET(METIS_INCLUDEDIR "c:/Tools/metis-5.1.0/include")
  SET(METIS_DEBUG_LIBRARY "c:/Tools/metis-5.1.0/build/libmetis/Debug/metis.lib") 
  SET(METIS_RELEASE_LIBRARY "c:/Tools/metis-5.1.0/build/libmetis/Release/metis.lib") 
ENDIF()

##################################################################################
#  FETOL
##################################################################################
IF(${USE_FETOL})
  SET(FETOL_INCLUDEDIR "d:/Projects/FETOL/dev/CppFETOLlib")
  SET(FETOL_DEBUG_LIBRARY "d:/Projects/FETOL/dev/CppFETOLlib/build/Debug/fetol.lib") 
  SET(FETOL_RELEASE_LIBRARY "d:/Projects/FETOL/dev/CppFETOLlib/build/Release/fetol.lib") 
  
  SET(YAML_INCLUDEDIR "d:/Tools/yaml-cpp/include")
  SET(YAML_DEBUG_LIBRARY "d:/Tools/yaml-cpp/buildVS11/Debug/libyaml-cppmdd.lib") 
  SET(YAML_RELEASE_LIBRARY "d:/Tools/yaml-cpp/buildVS11/Release/libyaml-cppmd.lib") 
  
  SET(BOND_INCLUDEDIR "d:/Projects/FETOL/dev/bond_src/cpp/bond/fetol")
  SET(BOND_DEBUG_LIBRARY "d:/Projects/FETOL/dev/bond_lib/Debug/bond.lib") 
  SET(BOND_RELEASE_LIBRARY "d:/Projects/FETOL/dev/bond_lib/Release/bond.lib")   
ENDIF()

##################################################################################
#  Java
##############################################################################
### FindJNI.cmake
#find_package(JNI REQUIRED) 
#SET(JNI_INCLUDE_DIRS ${JAVA_INCLUDE_PATH} ${JAVA_INCLUDE_PATH2} ${JAVA_AWT_INCLUDE_PATH})
#SET(JNI_LIBRARIES ${JAVA_AWT_LIBRARY} ${JAVA_JVM_LIBRARY})
#SET(JNI_FOUND 1) 