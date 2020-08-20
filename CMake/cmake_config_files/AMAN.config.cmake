########################################################################################
##                                      CPU                                           ##
########################################################################################
SET(CPU_TYPE "Opteron")

########################################################################################
##                            BOOST ALLGMEINGUELTIG                                   ##
########################################################################################
#standard boost: 1.37.0
SET(BOOST_VERSION "1.37.0" CACHE STRING "std: 1.37.0")

########################################################################################
##                                      gcc42                                         ##
########################################################################################
IF(CAB_COMPILER MATCHES "gcc42")
   SET(BOOST_USE_MULTITHREAD ON)

   SET(BOOST_COMPILER_SUFFIX -gcc42)
   SET(BOOST_USE_STATIC_LIBS ON)
		
   IF(BOOST_VERSION MATCHES "1.37.0" )
      SET(BOOST_INCLUDEDIR "/opt/boost/gcc-4.2/include/boost-1_37")
      SET(BOOST_LIBRARYDIR "/opt/boost/gcc-4.2/lib")
   ENDIF()

   SET(MPI_DIR  "/opt/mpich-gm/gcc-4.2")

########################################################################################
##                                      gcc41                                         ##
########################################################################################
ELSEIF(CAB_COMPILER MATCHES "gcc41")
   SET(BOOST_USE_MULTITHREAD ON)

   SET(BOOST_COMPILER_SUFFIX -gcc41)
   SET(BOOST_USE_STATIC_LIBS ON)
		
   IF(BOOST_VERSION MATCHES "1.37.0" )
      SET(BOOST_INCLUDEDIR "/opt/boost/gcc-4.1/include/boost-1_37")
      SET(BOOST_LIBRARYDIR "/opt/boost/gcc-4.1/lib")
   ENDIF()

   SET(MPI_DIR  "/opt/mpich-gm/gcc-4.1")

########################################################################################
##                                      gcc34                                         ##
########################################################################################
ELSEIF(CAB_COMPILER MATCHES "gcc34")
   SET(BOOST_USE_MULTITHREAD ON)

   SET(BOOST_COMPILER_SUFFIX -gcc34)
   SET(BOOST_USE_STATIC_LIBS ON)
		
   IF(BOOST_VERSION MATCHES "1.37.0" )
      SET(BOOST_INCLUDEDIR "/opt/boost/gcc-3.4/include/boost-1_37")
      SET(BOOST_LIBRARYDIR "/opt/boost/gcc-3.4/lib")
   ENDIF()

   SET(MPI_DIR  "/opt/mpich-gm/gcc-3.4")

########################################################################################
##                                      intel10.1                                     ##
########################################################################################
ELSEIF(CAB_COMPILER MATCHES "intel10" )
  SET(BOOST_USE_MULTITHREAD ON)

  SET(BOOST_COMPILER_SUFFIX -il)
  SET(BOOST_USE_STATIC_LIBS ON)
  
  IF(BOOST_VERSION MATCHES "1.37.0" )
    #SET(BOOST_INCLUDEDIR "/opt/boost/intel-10.1/include/boost-1_37")
    #SET(BOOST_LIBRARYDIR "/opt/boost/intel-10.1/lib")
   ENDIF()

   #SET(MPI_DIR  "/opt/mpich/intel-10.1")
ENDIF()


IF(NEED_BOOST AND BOOST_VERSION AND NOT BOOST_INCLUDEDIR)
  MESSAGE("${BOOST_VERSION} not found on ${CAB_MACHINE} for specified compiler")
ENDIF()



