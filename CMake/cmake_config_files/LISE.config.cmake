#################################################################################
#  METIS  
#################################################################################
IF(${USE_METIS})
  SET(METIS_INCLUDEDIR "/home/niikonst/metis-5.1.0/include")
  SET(METIS_DEBUG_LIBRARY "/home/niikonst/metis-5.1.0/build/Linux-x86_64/libmetis/libmetis.a") 
  SET(METIS_RELEASE_LIBRARY "/home/niikonst/metis-5.1.0/build/Linux-x86_64/libmetis/libmetis.a")
ENDIF()
#################################################################################
#  BOOST  
#################################################################################
SET(BOOST_VERSION "1.72.0")
SET(BOOST_ROOT "/sw/tools/boost/1.72.0/skl/openmpi.3.1.5-gcc.9.2.0")
SET(BOOST_DIR ${BOOST_ROOT})
SET(BOOST_LIBRARYDIR ${BOOST_ROOT}"/lib")  