#################################################################################
# VirtualFluids MACHINE FILE
# Responsible: Konstantin Kutscher
# OS:          CentOS 7.3
#################################################################################

#################################################################################
#  PE (legacy)
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
  SET(GEOMETRY_DEBUG_LIBRARY ${PE_BINARY_DIR}/src/geometry/libgeometry.a)
  SET(GEOMETRY_RELEASE_LIBRARY ${PE_BINARY_DIR}/src/geometry/libgeometry.a)
  SET(CORE_DEBUG_LIBRARY ${PE_BINARY_DIR}/src/core/libcore.a)
  SET(CORE_RELEASE_LIBRARY ${PE_BINARY_DIR}/src/core/libcore.a)
ENDIF()

## nvidia
set(CMAKE_CUDA_ARCHITECTURES 60) # NVIDIA Tesla P100

set(GPU_APP "apps/gpu/LBM/")
list(APPEND USER_APPS 
    "${GPU_APP}DrivenCavityMultiGPU"
    # "${GPU_APP}SphereScaling"
    # "${GPU_APP}MusselOyster"
    )
