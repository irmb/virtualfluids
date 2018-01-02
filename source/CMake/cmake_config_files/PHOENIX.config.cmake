LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__unix__)
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__UNIX__)

SET(BOOST_VERSION "1.63.0" CACHE STRING "std: 1.63.0")

#################################################################################
#  METIS  
#################################################################################
IF(${USE_METIS})
  SET(METIS_INCLUDEDIR "/cluster/lib/metis/5.1.0/intel/include")
  SET(METIS_DEBUG_LIBRARY "/cluster/lib/metis/5.1.0/intel/lib/libmetis.a") 
  SET(METIS_RELEASE_LIBRARY "/cluster/lib/metis/5.1.0/intel/lib/libmetis.a") 
ENDIF()
