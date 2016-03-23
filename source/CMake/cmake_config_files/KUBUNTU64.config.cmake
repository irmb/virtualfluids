#################################################################################
# COMPILER FLAG VAR
#################################################################################
SET(CAB_COMPILER "gcc43" FORCE) 
SET(BOOST_VERSION "1.37.0")

#################################################################################
#  BOOST AND MPI
#################################################################################
IF(CAB_COMPILER MATCHES "gcc43")
	IF(NEED_BOOST)
	  SET(BOOST_USE_MULTITHREAD ON)
	  SET(BOOST_USE_STATIC_LIBS ON)
	  
	  IF(BOOST_VERSION MATCHES "1.37.0" )

	     SET(USER_DEFINED_BOOST_COMPILER -gcc43)
        SET(BOOST_INCLUDEDIR "/opt/boost_1_37_0/gcc4.3.1/include/boost-1_37")
	     SET(BOOST_LIBRARYDIR "/opt/boost_1_37_0/gcc4.3.1/lib")
        SET(BOOST_COMPILER_SUFFIX -gcc43)
   
     ENDIF()
   ENDIF()
ENDIF()

IF(BOOST_VERSION AND NOT BOOST_INCLUDEDIR)
  MESSAGE("${BOOST_VERSION} not found on ${CAB_MACHINE} for specified compiler")
ENDIF()
