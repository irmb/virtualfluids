OPTION(CAB_COMPILER_64BIT "64BIT COMPILATION"  ON)

IF(NOT CAB_COMPILER})
   SET(BOOST_VERSION "1.37.0")
ENDIF()

IF(NOT CAB_COMPILER})
  SET(CAB_COMPILER "intel10")
ENDIF()

SET(CAB_COMPILER ${CAB_COMPILER} CACHE STRING "intel10" FORCE) 
SET(CPU_TYPE "Opteron")

#################################################################################
#  BOOST AND MPI
#################################################################################
IF(CAB_COMPILER MATCHES "intel10" )
########################################################################################
##                                      intel10                                      ##
########################################################################################
   SET(BOOST_USE_MULTITHREAD ON)
   
   IF(BOOST_VERSION MATCHES "1.37.0" )
   
     SET(BOOST_USE_STATIC_LIBS ON)

     SET(BOOST_INCLUDEDIR "/home/freudigx/sw/boost/intel/boost_1_37_0")
     SET(BOOST_LIBRARYDIR "/home/freudigx/sw/boost/intel/boost_1_37_0/lib")

     SET( BOOST_COMPILER_SUFFIX "-il" ) 
   ENDIF()

   SET(MPI_DIR  "/sw/mpi/mvapich-1.0-intel101")

   SET_CXX_COMPILER("icpc")
   SET_C_COMPILER("icc")
ENDIF()

IF(BOOST_VERSION AND NOT BOOST_INCLUDEDIR)
  MESSAGE("${BOOST_VERSION} not found on "${CAB_MACHINE}" for specified compiler")
ENDIF()
