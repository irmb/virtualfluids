LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__unix__)
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__UNIX__)

########################################################################################
##                            BOOST ALLGMEINGUELTIG                                   ##
########################################################################################

SET(BOOST_VERSION "1.54.0" CACHE STRING "std: 1.54.0")

# IF(BOOST_VERSION MATCHES "1.51.0" )
#   SET(BOOST_INCLUDEDIR "/zhome/academic/HLRS/xrm/xrmkuchr/tools/boost_1_51_0")
#   SET(BOOST_LIBRARYDIR "/zhome/academic/HLRS/xrm/xrmkuchr/tools/boost_1_51_0/stageGCC/lib")
# 
# ENDIF()
# 
# IF(BOOST_VERSION AND NOT BOOST_INCLUDEDIR)
#   MESSAGE("${BOOST_VERSION} not found on ${CAB_MACHINE} for specified compiler")
# ENDIF()


#################################################################################
#  METIS  
#################################################################################
IF(${USE_METIS})
  SET(METIS_INCLUDEDIR "/zhome/academic/HLRS/xrm/xrmkuchr/tools/metis-5.0.2/include")
  SET(METIS_DEBUG_LIBRARY "/zhome/academic/HLRS/xrm/xrmkuchr/tools/metis-5.0.2/build/Linux-x86_64/libmetis/libmetis.a") 
  SET(METIS_RELEASE_LIBRARY "/zhome/academic/HLRS/xrm/xrmkuchr/tools/metis-5.0.2/build/Linux-x86_64/libmetis/libmetis.a") 
ENDIF()
