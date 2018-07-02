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

#################################################################################
#  PE  
#################################################################################
IF(${USE_DEM_COUPLING})
  SET(PE_BINARY_DIR "/home/irmb/walberla-git/build" CACHE PATH "pe binary dir")
  SET(PE_ROOT "/home/irmb/walberla-git" CACHE PATH "pe root")
 
  SET(PE_DEBUG_LIBRARY ${PE_BINARY_DIR}/src/pe/libpe.a) 
  SET(PE_RELEASE_LIBRARY ${PE_BINARY_DIR}/src/pe/libpe.a)
  SET(BLOCKFOREST_DEBUG_LIBRARY ${PE_BINARY_DIR}/src/blockforest/libblockforest.a) 
  SET(BLOCKFOREST_RELEASE_LIBRARY ${PE_BINARY_DIR}/src/blockforest/libblockforest.a)
  SET(DOMAIN_DECOMPOSITION_DEBUG_LIBRARY ${PE_BINARY_DIR}/src/domain_decomposition/libdomain_decomposition.a) 
  SET(DOMAIN_DECOMPOSITION_RELEASE_LIBRARY ${PE_BINARY_DIR}/src/domain_decomposition/libdomain_decomposition.a)
  SET(CORE_DEBUG_LIBRARY ${PE_BINARY_DIR}/src/core/libcore.a) 
  SET(CORE_RELEASE_LIBRARY ${PE_BINARY_DIR}/src/core/libcore.a)
  
ENDIF()