#################################################################################
#  BOOST  
#################################################################################
SET(BOOST_VERSION "1.60.0")
SET(BOOST_ROOT "d:/boost/boost_1_60_0")
SET(BOOST_DIR ${BOOST_ROOT})
SET(BOOST_LIBRARYDIR ${BOOST_ROOT}"/stageMSVC64/lib")  
#################################################################################

#################################################################################
#  METIS  
#################################################################################
IF(${USE_METIS})
  SET(METIS_INCLUDEDIR "d:/metis-5.1.0/include")
  SET(METIS_DEBUG_LIBRARY "d:/metis-5.1.0/build/libmetis/Debug/metis.lib") 
  SET(METIS_RELEASE_LIBRARY "d:/metis-5.1.0/build/libmetis/Release/metis.lib") 
ENDIF()