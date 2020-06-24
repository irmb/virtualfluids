#################################################################################
#  BOOST  
#################################################################################
SET(BOOST_VERSION "1.65.1")
SET(BOOST_ROOT "c:/Libraries/boost/boost_1_65_1")
SET(BOOST_DIR ${BOOST_ROOT})
SET(BOOST_LIBRARYDIR ${BOOST_ROOT}"/stage/lib")  
#################################################################################

#################################################################################
#################################################################################
#  METIS  
#################################################################################
IF(${USE_METIS})
  SET(METIS_INCLUDEDIR "c:/Libraries/metis/metis-5.1.0/include")
  SET(METIS_DEBUG_LIBRARY "c:/Libraries/metis/metis-5.1.0/MSVC2015/libmetis/Debug/metis.lib") 
  SET(METIS_RELEASE_LIBRARY "c:/Libraries/metis/metis-5.1.0/MSVC2015/libmetis/Release/metis.lib") 
ENDIF()