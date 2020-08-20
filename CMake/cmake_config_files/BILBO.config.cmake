LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__unix__)
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__UNIX__)
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__APPLE__)

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
#SET(BOOST_VERSION "1.47")
#SET(BOOST_USE_MULTITHREAD ON)
#SET(BOOST_USE_STATIC_LIBS ON)
#
#IF(BOOST_VERSION MATCHES "1.47" )
#  SET(BOOST_ROOT "/host/tools/boost/boost_1_47_0")
#  SET(BOOST_LIBRARYDIR "/host/tools/boost/boost_1_47_0/stageLinux/lib")
#ENDIF()
#################################################################################
#  ZOLTAN
#################################################################################
#IF(${USE_ZOLTAN})
#  SET(ZOLTAN_INCLUDEDIR "/home/kostja/tools/Zoltan_v3.5/bin/include")
#  SET(ZOLTAN_DEBUG_LIBRARY "/home/kostja/tools/Zoltan_v3.5/bin/lib/libzoltan.a")
#  SET(ZOLTAN_RELEASE_LIBRARY "/home/kostja/tools/Zoltan_v3.5/bin/lib/libzoltan.a")
#ENDIF()
#################################################################################
#  METIS
#################################################################################
IF(${USE_METIS})
    SET(METIS_INCLUDEDIR "/usr/local/include")
    SET(METIS_DEBUG_LIBRARY "/usr/local/lib/libmetis.a")
    SET(METIS_RELEASE_LIBRARY "/usr/local/lib/libmetis.a")
ENDIF()
#set(SOURCE_ROOT "/host/Projects/pFluid/source")

SET(GMOCK_ROOT "~/Libraries/googletest-release-1.10.0" CACHE PATH "GMOCK ROOT")