LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__unix__)
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__UNIX__)

#################################################################################
#  METIS  
#################################################################################
IF(${USE_METIS})
  SET(METIS_INCLUDEDIR "/home/niikonst/libs/metis-5.1.0/include")
  SET(METIS_DEBUG_LIBRARY "/home/niikonst/libs/metis-5.1.0/build/Linux-x86_64/libmetis/libmetis.a") 
  SET(METIS_RELEASE_LIBRARY "/home/niikonst/libs/metis-5.1.0/build/Linux-x86_64/libmetis/libmetis.a") 
ENDIF()


