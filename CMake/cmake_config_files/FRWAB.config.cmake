OPTION(CAB_COMPILER_64BIT "64BIT COMPILATION"  ON)

########################################################################################
##                            BOOST ALLGMEINGUELTIG                                   ##
########################################################################################
#standard boost: 1.38.0
SET(BOOST_VERSION "1.38.0" CACHE STRING "std: 1.38.0")

########################################################################################
##                                      gcc41                                         ##
########################################################################################
IF(CAB_COMPILER MATCHES "gcc41")
   SET(BOOST_USE_MULTITHREAD ON)

   SET(BOOST_COMPILER_SUFFIX -gcc41)
   SET(BOOST_USE_STATIC_LIBS ON)
		
   IF(BOOST_VERSION MATCHES "1.38.0" )
      SET(BOOST_INCLUDEDIR "/nfs/HOME/HLRS/iws/iwsirmb/boost/gcc41/boost_1_38_0/include/boost-1_38")
      SET(BOOST_LIBRARYDIR "/nfs/HOME/HLRS/iws/iwsirmb/boost/gcc41/boost_1_38_0/lib")
   ENDIF()

   SET(MPI_DIR  "/opt/bwgrid/mpi/openmpi/1.3-gnu-4.1")
   SET(MPI_CXX_LIBRARY libmpi_cxx.so)
   SET(MPI_C_LIBRARY libmpi.so)
  
########################################################################################
##                                      intel10.1                                     ##
########################################################################################
ELSEIF(CAB_COMPILER MATCHES "intel10" )
  SET(BOOST_USE_MULTITHREAD ON)

  SET(BOOST_COMPILER_SUFFIX -il)
  SET(BOOST_USE_STATIC_LIBS ON)
  
  IF(BOOST_VERSION MATCHES "1.37.0" )
      SET(BOOST_INCLUDEDIR "/nfs/HOME/HLRS/iws/iwsirmb/boost/intel101/boost_1_38_0/include/boost-1_38")
      SET(BOOST_LIBRARYDIR "/nfs/HOME/HLRS/iws/iwsirmb/boost/intel101/boost_1_38_0/lib")
   ENDIF()

   #SET(MPI_DIR  "/opt/mpich/intel-10.1")
   SET(MPI_DIR  "/opt/bwgrid/mpi/openmpi/1.3-intel-10.1")

   SET_CXX_COMPILER("/opt/intel/vve/10.1.018/bin/icpc")
   SET_C_COMPILER("/opt/intel/vve/10.1.018/bin/icc")
ENDIF()


IF(BOOST_VERSION AND NOT BOOST_INCLUDEDIR)
  MESSAGE("${BOOST_VERSION} not found on ${CAB_MACHINE} for specified compiler")
ENDIF()



