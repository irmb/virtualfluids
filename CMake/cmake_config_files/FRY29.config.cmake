LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__unix__)
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__UNIX__)

########################################################################################
##                            BOOST ALLGMEINGUELTIG                                   ##
########################################################################################
SET(BOOST_VERSION "1.60.0" CACHE STRING "std: 1.60.0")


SET(USE_MPI_CXX_SYNTAX OFF)
SET(BOOST_INCLUDEDIR "/home/irmb/kutscher/boost_1_60_0/boost")
SET(BOOST_LIBRARYDIR "/home/irmb/kutscher/boost_1_60_0/stage/lib")
SET(Boost_INCLUDE_DIR "/home/irmb/kutscher/boost_1_60_0")

IF(BOOST_VERSION AND NOT BOOST_INCLUDEDIR)
  MESSAGE("${BOOST_VERSION} not found on ${CAB_MACHINE} for specified compiler")
ENDIF()

#################################################################################
#  METIS  
#################################################################################
IF(${USE_METIS})
  SET(METIS_INCLUDEDIR "/home/irmb/kutscher/metis-5.1.0/include")
  SET(METIS_DEBUG_LIBRARY "/home/irmb/kutscher/metis-5.1.0/build/Linux-x86_64/libmetis/libmetis.a") 
  SET(METIS_RELEASE_LIBRARY "/home/irmb/kutscher/metis-5.1.0/build/Linux-x86_64/libmetis/libmetis.a") 
ENDIF()



