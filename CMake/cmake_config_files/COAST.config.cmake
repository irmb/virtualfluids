
#SET(CAB_COMPILER gcc41 CACHE STRING "gcc41" FORCE) 
CHECK_FOR_VARIABLE(CAB_COMPILER "gcc41, geller_gcc41")

#################################################################################
# COMPILER FLAG VAR
#################################################################################
IF(NOT CAB_COMPILER})
  IF(NOT DEFINED CAB_COMPILER)
    SET(CAB_COMPILER "CAB_COMPILER-NOTFOUND")
    MESSAGE("CHECK_FOR_VARIABLE - error - set CAB_COMPILER")
  ENDIF()
ENDIF()

SET(CPU_TYPE "Opteron")

#standard boost: 1.37.0
SET(BOOST_VERSION "1.37.0" CACHE STRING "std: 1.37.0")

########################################################################################
##                                      gcc41                                         ##
########################################################################################
IF(CAB_COMPILER MATCHES "gcc41")
   SET(BOOST_USE_MULTITHREAD ON)

   SET(BOOST_COMPILER_SUFFIX -gcc41)
   SET(BOOST_USE_STATIC_LIBS ON)
		
   IF(BOOST_VERSION MATCHES "1.37.0" )
      SET(BOOST_INCLUDEDIR "/home/hegewald/bin/boost/boost_1_37_0")
      SET(BOOST_LIBRARYDIR "/home/hegewald/bin/boost/boost_1_37_0/lib")
   ENDIF()

   #SET(MPI_DIR  "/opt/mpich-gm/gcc-4.1")

   SET_CXX_COMPILER("/usr/bin/g++-4.1")
   SET_C_COMPILER("/usr/bin/gcc-4.1")

ELSEIF(CAB_COMPILER MATCHES "geller_gcc41")

   SET(BOOST_USE_MULTITHREAD ON)

   SET(BOOST_COMPILER_SUFFIX -gcc41)
   SET(BOOST_USE_STATIC_LIBS ON)
		
   IF(BOOST_VERSION MATCHES "1.37.0" )
      SET(BOOST_INCLUDEDIR "/home/hegewald/bin/boost/boost_1_37_0")
      SET(BOOST_LIBRARYDIR "/home/hegewald/bin/boost/boost_1_37_0/lib")
   ENDIF()

   SET(MPI_DIR  "/home/geller/bin/mpich2-1.0.8")

   SET_CXX_COMPILER("/usr/bin/g++-4.2")
   SET_C_COMPILER("/usr/bin/gcc-4.2")
ENDIF()


IF(BOOST_VERSION AND NOT BOOST_INCLUDEDIR)
  MESSAGE("${BOOST_VERSION} not found on ${CAB_MACHINE} for specified compiler")
ENDIF()



