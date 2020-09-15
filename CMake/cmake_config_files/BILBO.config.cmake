#################################################################################
# VirtualFluids MACHINE FILE
# Responsible: Soeren Peters
# OS:          MacOS X
#################################################################################

LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__unix__)
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__UNIX__)
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__APPLE__)

#################################################################################
#  METIS
#################################################################################
SET(METIS_INCLUDEDIR "/usr/local/include")
SET(METIS_DEBUG_LIBRARY "/usr/local/lib/libmetis.a")
SET(METIS_RELEASE_LIBRARY "/usr/local/lib/libmetis.a")

