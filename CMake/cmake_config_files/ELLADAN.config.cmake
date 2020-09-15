#################################################################################
# VirtualFluids MACHINE FILE
# Responsible: Soeren Peters
# OS:          Ubuntu 20.04
#################################################################################

LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__unix__)
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__UNIX__)

set(NVCUDASAMPLES_ROOT "~/cuda-samples/Common")
#################################################################################
#  METIS
#################################################################################
set(METIS_INCLUDEDIR "/usr/include")
set(METIS_DEBUG_LIBRARY "/usr/lib/x86_64-linux-gnu/libmetis.so")
set(METIS_RELEASE_LIBRARY "/usr/lib/x86_64-linux-gnu/libmetis.so")

