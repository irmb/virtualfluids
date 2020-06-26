IF(NOT CMAKE_CXX_COMPILER)
   MESSAGE(FATAL_ERROR "before cmake-config-file can be included -> project must be extecuted")
ENDIF()
  
#################################################################################
# MPI
#################################################################################
SET(MPI_DIR  "C:/Program Files/mpich2")

#################################################################################
#  BOOST  
#################################################################################
SET(BOOST_VERSION "1.38")
SET(BOOST_USE_MULTITHREAD ON)
SET(BOOST_USE_STATIC_LIBS ON)

IF(BOOST_VERSION MATCHES "1.38" )
  SET(BOOST_INCLUDEDIR "C:/Program Files/boost/boost_1_38")
  SET(BOOST_LIBRARYDIR "C:/Program Files/boost/boost_1_38/lib")
ENDIF()
  

