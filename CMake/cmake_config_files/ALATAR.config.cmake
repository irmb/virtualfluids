LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__unix__)
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__UNIX__)

########################################################################################
##                            BOOST ALLGMEINGUELTIG                                   ##
########################################################################################
#standard boost: 1.50.0
SET(BOOST_VERSION "1.51.0" CACHE STRING "std: 1.51.0")

#SET(MPI_DIR  "/hpc3lustre/software/lib/openmpi-1.4.1-gcc41")

SET(USE_MPI_CXX_SYNTAX OFF)
#SET(MPI_COMPILER "/hpc3lustre/software/lib/openmpi-1.4.1-gcc41/bin/mpicxx")

IF(CAB_COMPILER MATCHES "gcc41")
   SET(BOOST_USE_MULTITHREAD ON)

   SET(BOOST_COMPILER_SUFFIX -gcc41)
   SET(BOOST_USE_STATIC_LIBS ON)

   IF(BOOST_VERSION MATCHES "1.51.0" )
      SET(BOOST_INCLUDEDIR "/hpc3lustre/software/irmb/boost/boost_1_51_0")
      SET(BOOST_LIBRARYDIR "/hpc3lustre/software/irmb/boost/boost_1_51_0/stageGCC/lib")

   ENDIF()
ENDIF()

IF(BOOST_VERSION AND NOT BOOST_INCLUDEDIR)
  MESSAGE("${BOOST_VERSION} not found on ${CAB_MACHINE} for specified compiler")
ENDIF()

#################################################################################
#  ZOLTAN  
#################################################################################
IF(${USE_ZOLTAN})
  SET(ZOLTAN_INCLUDEDIR "/hpc3lustre/software/irmb/Trilinos/trilinos-10.6.4-Source/packages/zoltan/src/include")
  SET(ZOLTAN_DEBUG_LIBRARY "/hpc3lustre/software/irmb/Trilinos/trilinos-10.6.4-Build/packages/zoltan/src/libzoltan.a") 
  SET(ZOLTAN_RELEASE_LIBRARY "/hpc3lustre/software/irmb/Trilinos/trilinos-10.6.4-Build/packages/zoltan/src/libzoltan.a")
ENDIF()

#################################################################################
#  METIS  
#################################################################################
IF(${USE_METIS})
  SET(METIS_INCLUDEDIR "/hpc3lustre/software/irmb/metis/include")
  SET(METIS_DEBUG_LIBRARY "/hpc3lustre/software/irmb/metis/lib/libmetis.a") 
  SET(METIS_RELEASE_LIBRARY "/hpc3lustre/software/irmb/metis/lib/libmetis.a") 
ENDIF()

#################################################################################
#  YAML  
#################################################################################
IF(${USE_YAML})
  SET(YAML_INCLUDEDIR "/hpc3lustre/software/irmb/yaml/yaml-cpp/include")
  SET(YAML_DEBUG_LIBRARY "/hpc3lustre/software/irmb/yaml/yaml-cpp/build/libyaml-cpp.a") 
  SET(YAML_RELEASE_LIBRARY "/hpc3lustre/software/irmb/yaml/yaml-cpp/build/libyaml-cpp.a") 
ENDIF()

##################################################################################
#  Bond
##################################################################################
IF(${USE_BOND})
  SET(BOND_INCLUDEDIR "/hpc3lustre/home/koskuche/projects/bond/bond_src/cpp/bond/fetol")
  SET(BOND_DEBUG_LIBRARY "/hpc3lustre/home/koskuche/projects/bond/bin/lib/libbond.a") 
  SET(BOND_RELEASE_LIBRARY "/hpc3lustre/home/koskuche/projects/bond/bin/lib/libbond.a") 
ENDIF()

