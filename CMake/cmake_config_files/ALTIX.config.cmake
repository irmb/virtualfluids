OPTION(CAB_COMPILER_64BIT "64BIT COMPILATION"  ON)

#################################################################################
# COMPILER 
#################################################################################
SET( ALTIX_CXX_COMPILER  "mpiCC"   )
SET( ALTIX_C_COMPILER    "mpicc"   )

#################################################################################
# COMPILER FLAG VAR
#################################################################################
IF(NOT CAB_COMPILER})
  IF(NOT DEFINED CAB_COMPILER)
    SET( CAB_COMPILER "intel9" )
  ENDIF()
ENDIF()

SET(BOOST_DIR  "/home/hlrb2/h005x/h005xac/boost/include/boost-1_34_1")
SET(USE_BOOST_STATIC_LIBS TRUE)

#################################################################################
# MPI
#################################################################################


#################################################################################
#  BOOST  
#################################################################################
IF(NEED_BOOST)
  SET(BOOST_USE_MULTITHREAD ON)
  SET(BOOST_USE_STATIC_LIBS ON)
  
  #IF(BOOST_VERSION MATCHES "1.35.0" )
  #
  #   SET(BOOST_INCLUDEDIR "E:/3rdPartyLibs/boost/boost_1_35")
  #
  #   IF(CMAKE_CL_64) 
  #     SET(BOOST_LIBRARYDIR "E:/3rdPartyLibs/boost/boost_1_35/lib/x86_64bit")
  #   ELSE()
  #     SET(BOOST_LIBRARYDIR "E:/3rdPartyLibs/boost/boost_1_35/lib/x86_32bit")
  #   ENDIF()
  
  #ELSE
  IF(BOOST_VERSION MATCHES "1.34.1" )
     SET(BOOST_INCLUDEDIR "/home/hlrb2/h005x/h005xac/boost/include/boost-1_34_1")
     SET(BOOST_LIBRARYDIR "/home/hlrb2/h005x/h005xac/boost/include/boost-1_34_1/lib")
  ENDIF()
  
  IF(BOOST_VERSION AND NOT BOOST_INCLUDEDIR)
    MESSAGE("${BOOST_VERSION} not found on ${CAB_MACHINE} for specified compiler")
  ENDIF()
ENDIF()
  
