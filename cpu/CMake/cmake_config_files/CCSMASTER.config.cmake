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
IF(CMAKE_CL_64)
  SET(MPI_DIR  "c:/Program Files/MPICH2")
ELSE()
  SET(MPI_DIR  "c:/Program Files (x86)/MPICH2")
ENDIF()

#################################################################################
#  BOOST  
#################################################################################
IF(NEED_BOOST)
  SET(BOOST_USE_MULTITHREAD ON)
  SET(BOOST_USE_STATIC_LIBS ON)
  
  IF(BOOST_VERSION MATCHES "1.35.0" )
  
     SET(BOOST_INCLUDEDIR "Y:/boost_1_35_0")
  
     IF(CMAKE_CL_64) 
       SET(BOOST_LIBRARYDIR "Y:/boost_1_35_0/lib/x86_64bit")
     ELSE()
       SET(BOOST_LIBRARYDIR "Y:/boost_1_35_0/lib/x86_32bit")
     ENDIF()
  
  ELSEIF(BOOST_VERSION MATCHES "1.34.1" )
     SET(BOOST_INCLUDEDIR "H:/boost1.34.1")
          IF(CMAKE_CL_64) 
       SET(BOOST_LIBRARYDIR "H:/boost1.34.1/stage64/lib")
     ELSE()
       SET(BOOST_LIBRARYDIR "H:/boost1.34.1/lib")
     ENDIF()
  ENDIF()
  
  IF(BOOST_VERSION AND NOT BOOST_INCLUDEDIR)
    MESSAGE("${BOOST_VERSION} not found on ${CAB_MACHINE} for specified compiler")
  ENDIF()
ENDIF()

