#################################################################################
# COMPILER FLAG VAR
#################################################################################
IF(NOT CAB_COMPILER})
  IF( ${MSVC_VERSION} MATCHES "1400" AND CMAKE_CL_64)
      SET( CAB_COMPILER "msvc8_64" )
  ELSEIF( ${MSVC_VERSION} MATCHES "1400" AND NOT CMAKE_CL_64)
      SET( CAB_COMPILER "msvc8_32" )
  ELSEIF( ${MSVC_VERSION} MATCHES "1500" AND CMAKE_CL_64)
      SET( CAB_COMPILER "msvc9_64" )
  ELSEIF( ${MSVC_VERSION} MATCHES "1500" AND NOT CMAKE_CL_64)
     SET( CAB_COMPILER "msvc9_32" )
  ELSE()
     MESSAGE(ERROR, "unknown ms studio version MSVC_VERSION = "${MSVC_VERSION})  
  ENDIF()
ENDIF()

#################################################################################
# MPI
#################################################################################
SET(MPI_DIR  "c:/programme/mpich2")

#################################################################################
#  BOOST  
#################################################################################
IF(NEED_BOOST)

  SET(BOOST_USE_MULTITHREAD ON)
  SET(BOOST_USE_STATIC_LIBS ON)
  
  IF(BOOST_VERSION MATCHES "1.35.0" )
  
     SET(BOOST_INCLUDEDIR "c:/cpp/boost/boost_1_35_0")
  
     IF(CMAKE_CL_64) 
     #  SET(BOOST_LIBRARYDIR "c:/cpp/boost/boost_1_35_0/lib/x86_64bit")
     ELSE()
       SET(BOOST_LIBRARYDIR "c:/cpp/boost/boost_1_35_0/lib")
     ENDIF()

  ELSEIF(BOOST_VERSION MATCHES "1.38.0" )
     SET(BOOST_INCLUDEDIR "e:/cpp/boost/boost_1_38")
     SET(BOOST_LIBRARYDIR "e:/cpp/boost/boost_1_38/lib")

  ENDIF()
  
  IF(BOOST_VERSION AND NOT BOOST_INCLUDEDIR)
    MESSAGE("${BOOST_VERSION} not found on ${CAB_MACHINE} for specified compiler")
  ENDIF()
ENDIF()
  


