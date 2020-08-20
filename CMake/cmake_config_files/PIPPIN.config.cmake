#################################################################################
# MPI
#################################################################################
IF(CMAKE_CL_64)
  IF(NEED_MPI) 
     MESSAGE(FATAL_ERROR "kein 64 bit mpi installiert")
  ENDIF()
  #SET(MPI_DIR  "c:/Program Files/MPICH2")
ELSE()
  SET(MPI_DIR  "c:/Program Files (x86)/MPICH2")
ENDIF()

#################################################################################
#  BOOST  
#################################################################################
SET(BOOST_VERSION 1.38.0)
SET(BOOST_USE_MULTITHREAD ON)
SET(BOOST_USE_STATIC_LIBS ON)

IF(CMAKE_CL_64)
   MESSAGE(FATAL_ERROR "kein 64 bit boost ${BOOST_VERSION} installiert")
ELSE()
  SET(BOOST_INCLUDEDIR "C:/Program Files (x86)/boost/boost_1_38")
  SET(BOOST_LIBRARYDIR "C:/Program Files (x86)/boost/boost_1_38/lib")
ENDIF()
