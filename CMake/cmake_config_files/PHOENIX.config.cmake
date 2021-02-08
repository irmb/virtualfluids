#################################################################################
# VirtualFluids MACHINE FILE
# Responsible: Konstantin Kutscher
# OS:          CentOS 7.3
#################################################################################

SET(BOOST_VERSION "1.63.0" CACHE STRING "std: 1.63.0")

#################################################################################
#  METIS  
#################################################################################
SET(METIS_INCLUDEDIR "/cluster/lib/metis/5.1.0/gcc/include")
SET(METIS_DEBUG_LIBRARY "/cluster/lib/metis/5.1.0/gcc/lib/libmetis.a")
SET(METIS_RELEASE_LIBRARY "/cluster/lib/metis/5.1.0/gcc/lib/libmetis.a")


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
  SET(GEOMETRY_DEBUG_LIBRARY ${PE_BINARY_DIR}/src/geometry/libgeometry.a) 
  SET(GEOMETRY_RELEASE_LIBRARY ${PE_BINARY_DIR}/src/geometry/libgeometry.a)
  SET(CORE_DEBUG_LIBRARY ${PE_BINARY_DIR}/src/core/libcore.a) 
  SET(CORE_RELEASE_LIBRARY ${PE_BINARY_DIR}/src/core/libcore.a)
  
ENDIF()



SET(BOOST_ROOT  "/cluster/lib/boost/1.63.0/gcc"  CACHE PATH "BOOST_ROOT")
SET(BOOST_LIBRARYDIR  "/cluster/lib/boost/1.63.0/gcc/lib"  CACHE PATH "BOOST_LIBRARYDIR")

#SET(VTK_DIR "/cluster/lib/vtk/8.1.0/lib/cmake/vtk-8.1" CACHE PATH "VTK directory override" FORCE)
#SET(VTK_DIR "/home/irmb/tools/VTK/build/VTK-8.2.0" CACHE PATH "VTK directory override" FORCE)
#SET(VTK_DIR "/home/stelenz/software/vtk/VTK-8.1.0/build" CACHE PATH "VTK directory override" FORCE)

## nvidia
set(NVCUDASAMPLES_ROOT "/cluster/cuda/11.0/samples")
set(CMAKE_CUDA_ARCHITECTURES 60) # NVIDIA Tesla P100