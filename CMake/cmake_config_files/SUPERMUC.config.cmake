LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__unix__)
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__UNIX__)

########################################################################################
##                            BOOST ALLGMEINGUELTIG                                   ##
########################################################################################
#standard boost
SET(BOOST_VERSION "1.61.0" CACHE STRING "std: 1.61.0")
SET(BOOST_INCLUDEDIR "/lrz/sys/libraries/boost/1.61_icc/include")
SET(BOOST_LIBRARYDIR "/lrz/sys/libraries/boost/1.61_icc/lib")

#IF(BOOST_VERSION AND NOT BOOST_INCLUDEDIR)
#  MESSAGE("${BOOST_VERSION} not found on ${CAB_MACHINE} for specified compiler")
#ENDIF()

#################################################################################
#  METIS  
#################################################################################
SET(METIS_INCLUDEDIR "/lrz/sys/libraries/metis/5.1.0/i4r4/include")
SET(METIS_DEBUG_LIBRARY "/lrz/sys/libraries/metis/5.1.0/i4r4/lib/libmetis.a")
SET(METIS_RELEASE_LIBRARY "/lrz/sys/libraries/metis/5.1.0/i4r4/lib/libmetis.a")
