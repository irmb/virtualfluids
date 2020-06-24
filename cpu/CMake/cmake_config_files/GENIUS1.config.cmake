LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__unix__)
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__UNIX__)

IF(NOT CMAKE_CXX_COMPILER)
   MESSAGE(FATAL_ERROR "before cmake-config-file can be included -> project must be extecuted")
ENDIF()

#################################################################################
# MPI
#################################################################################
SET(USE_MPI_CXX_SYNTAX OFF)

#################################################################################
#  BOOST
#################################################################################
SET(BOOST_VERSION "1.44")
SET(BOOST_USE_MULTITHREAD ON)
SET(BOOST_USE_STATIC_LIBS ON)

IF(BOOST_VERSION MATCHES "1.44" )
  SET(BOOST_INCLUDEDIR "/u/lrz102bk/soft/boost_1_44_0/build/include")
  SET(BOOST_LIBRARYDIR "/u/lrz102bk/soft/boost_1_44_0/build/lib")
ENDIF()

IF(${USE_ZOLTAN})
  SET(ZOLTAN_INCLUDEDIR "/u/lrz102bk/BlueGene/soft/Zoltan_v3.1/bin/include")
  SET(ZOLTAN_DEBUG_LIBRARY "/u/lrz102bk/BlueGene/soft/Zoltan_v3.1/bin/lib/libzoltan.a") 
  SET(ZOLTAN_RELEASE_LIBRARY "/u/lrz102bk/BlueGene/soft/Zoltan_v3.1/bin/lib/libzoltan.a") 
ENDIF()

set(SOURCE_ROOT "/u/lrz102bk/BlueGene/projects/pFluid/source")


