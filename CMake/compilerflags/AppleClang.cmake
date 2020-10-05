#############################################################################################################
# compiler flags
#############################################################################################################

# debug
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_DEBUG "-g")  # generates debug information. Works best with -O0.
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_DEBUG "-O0")

# release
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_RELEASE "-O3") # optimization level (-O3: most optimization which also could result in larger binaries)

# all
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-fPIC") # position independent code for shared libraries

#############################################################################################################
# warnings
#############################################################################################################
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wall")
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wunreachable-code")

list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-unused-function")
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-reorder-ctor")

# temp:
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-sometimes-uninitialized")
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-inconsistent-missing-override")