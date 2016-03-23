LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__unix__)
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__UNIX__)

########################################################################################
##                            BOOST ALLGMEINGUELTIG                                   ##
########################################################################################
#standard boost
SET(BOOST_VERSION "1.54.0" CACHE STRING "std: 1.54.0")

#SET(USE_MPI_CXX_SYNTAX OFF)

#IF(CAB_COMPILER MATCHES "gcc41")
#   SET(BOOST_USE_MULTITHREAD ON)

   #SET(BOOST_COMPILER_SUFFIX -gcc41)
#   SET(BOOST_USE_STATIC_LIBS ON)

#    IF(BOOST_VERSION MATCHES "1.54.0" )
#       SET(BOOST_INCLUDEDIR "/gfs1/work/niivfcpu/tools/boost_1_54_0")
#       SET(BOOST_LIBRARYDIR "/gfs1/work/niivfcpu/tools/boost_1_54_0/stageGCC/lib")
#    ENDIF()
#ENDIF()

#IF(BOOST_VERSION AND NOT BOOST_INCLUDEDIR)
#  MESSAGE("${BOOST_VERSION} not found on ${CAB_MACHINE} for specified compiler")
#ENDIF()


#################################################################################
#  METIS  
#################################################################################
IF(${USE_METIS})
  SET(METIS_INCLUDEDIR "/gfs1/work/niivfcpu/libs/metis/include")
  SET(METIS_DEBUG_LIBRARY "/gfs1/work/niivfcpu/libs/metis/lib/libmetis.a") 
  SET(METIS_RELEASE_LIBRARY "/gfs1/work/niivfcpu/libs/metis/lib/libmetis.a") 
ENDIF()

# #################################################################################
# #  YAML  
# #################################################################################
# IF(${USE_YAML})
#   SET(YAML_INCLUDEDIR "/hpc3lustre/software/irmb/yaml/yaml-cpp/include")
#   SET(YAML_DEBUG_LIBRARY "/hpc3lustre/software/irmb/yaml/yaml-cpp/build/libyaml-cpp.a") 
#   SET(YAML_RELEASE_LIBRARY "/hpc3lustre/software/irmb/yaml/yaml-cpp/build/libyaml-cpp.a") 
# ENDIF()
# 
# ##################################################################################
# #  Bond
# ##################################################################################
# IF(${USE_BOND})
#   SET(BOND_INCLUDEDIR "/hpc3lustre/home/koskuche/projects/bond/bond_src/cpp/bond/fetol")
#   SET(BOND_DEBUG_LIBRARY "/hpc3lustre/home/koskuche/projects/bond/bin/lib/libbond.a") 
#   SET(BOND_RELEASE_LIBRARY "/hpc3lustre/home/koskuche/projects/bond/bin/lib/libbond.a") 
# ENDIF()

