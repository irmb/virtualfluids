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
  
    IF(BOOST_VERSION MATCHES "1.38.0" )
  
     SET(BOOST_INCLUDEDIR "E:/3rdPartyLibs/boost/boost_1_38_0")
  
     IF(CMAKE_CL_64) 
       SET(BOOST_LIBRARYDIR "E:/3rdPartyLibs/boost/boost_1_38_0/lib/x86_64bit")
     ELSE()
       SET(BOOST_LIBRARYDIR "E:/3rdPartyLibs/boost/boost_1_38_0/lib/x86_32bit")
     ENDIF()


    ELSEIF(BOOST_VERSION MATCHES "1.37.0" )
  
     SET(BOOST_INCLUDEDIR "E:/3rdPartyLibs/boost/boost_1_37_0")
  
     IF(CMAKE_CL_64) 
       SET(BOOST_LIBRARYDIR "E:/3rdPartyLibs/boost/boost_1_37_0/lib/x86_64bit")
     ELSE()
       SET(BOOST_LIBRARYDIR "E:/3rdPartyLibs/boost/boost_1_37_0/lib/x86_32bit")
     ENDIF()

   ELSEIF(BOOST_VERSION MATCHES "1.36.0" )
  
     SET(BOOST_INCLUDEDIR "E:/3rdPartyLibs/boost/boost_1_36_0")
  
     IF(CMAKE_CL_64) 
       SET(BOOST_LIBRARYDIR "E:/3rdPartyLibs/boost/boost_1_36_0/lib/x86_64bit")
     ELSE()
       SET(BOOST_LIBRARYDIR "E:/3rdPartyLibs/boost/boost_1_36_0/lib/x86_32bit")
     ENDIF()

  ELSEIF(BOOST_VERSION MATCHES "1.35.0.svn" )
  
     SET(BOOST_INCLUDEDIR "E:/3rdPartyLibs/boost/boost_1_35_svn")
  
     IF(CMAKE_CL_64) 
       SET(BOOST_LIBRARYDIR "E:/3rdPartyLibs/boost/boost_1_35_svn/lib/x86_64bit")
     ELSE()
       SET(BOOST_LIBRARYDIR "E:/3rdPartyLibs/boost/boost_1_35_svn/lib/x86_32bit")
     ENDIF()

  ELSEIF(BOOST_VERSION MATCHES "1.35.0.std" )
  
     SET(BOOST_INCLUDEDIR "E:/3rdPartyLibs/boost/boost_1_35_0_std")
  
     IF(CMAKE_CL_64) 
       SET(BOOST_LIBRARYDIR "E:/3rdPartyLibs/boost/boost_1_35_0_std/lib/x86_64bit")
     ELSE()
       SET(BOOST_LIBRARYDIR "E:/3rdPartyLibs/boost/boost_1_35_0_std/lib/x86_32bit")
     ENDIF()
  
  ELSEIF(BOOST_VERSION MATCHES "1.35.0" )
  
     SET(BOOST_INCLUDEDIR "E:/3rdPartyLibs/boost/boost_1_35_0")
  
     IF(CMAKE_CL_64) 
       SET(BOOST_LIBRARYDIR "E:/3rdPartyLibs/boost/boost_1_35_0/lib/x86_64bit")
     ELSE()
       SET(BOOST_LIBRARYDIR "E:/3rdPartyLibs/boost/boost_1_35_0/lib/x86_32bit")
     ENDIF()

  ELSEIF(BOOST_VERSION MATCHES "1.34.1" )
     SET(BOOST_INCLUDEDIR "E:/3rdPartyLibs/boost/boost_1_34_1")
     SET(BOOST_LIBRARYDIR "E:/3rdPartyLibs/boost/boost_1_34_1/lib")
  ENDIF()
  
  IF(BOOST_VERSION AND NOT BOOST_INCLUDEDIR)
    MESSAGE("${BOOST_VERSION} not found on ${CAB_MACHINE} for specified compiler")
  ENDIF()
ENDIF()
