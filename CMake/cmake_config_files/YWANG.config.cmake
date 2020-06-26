IF(NOT CMAKE_CXX_COMPILER)
   MESSAGE(FATAL_ERROR "before cmake-config-file can be included -> project must be extecuted")
ENDIF()
  
#################################################################################
#  BOOST  
#################################################################################
SET(BOOST_VERSION "1.49")
SET(BOOST_USE_MULTITHREAD ON)
SET(BOOST_USE_STATIC_LIBS ON)

IF(BOOST_VERSION MATCHES "1.49" )
  SET(BOOST_ROOT "C:/boost/boost_1_49_0")
  SET(BOOST_LIBRARYDIR "c:/boost/boost_1_49_0/stageMSVC64/lib")  
ENDIF()

#################################################################################
#  METIS  
#################################################################################
IF(${USE_METIS})
  SET(METIS_INCLUDEDIR "c:/metis-5.1.0/include")
  SET(METIS_DEBUG_LIBRARY "c:/metis-5.1.0/build/libmetis/Debug/metis.lib") 
  SET(METIS_RELEASE_LIBRARY "c:/metis-5.1.0/build/libmetis/Release/metis.lib") 
ENDIF()
