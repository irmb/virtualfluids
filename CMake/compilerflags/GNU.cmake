#############################################################################################################
# compiler flags
#############################################################################################################

# debug
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_DEBUG "-g")  # generates debug information. Works best with -O0.
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_DEBUG "-O0") # no optimization

# release
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_RELEASE "-O3") # optimization level (-O3: most optimization which also could result in larger binaries)

# all
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-fPIC") # position independent code for shared libraries
LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-funroll-all-loops")


#############################################################################################################
# warnings
#############################################################################################################
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wall")
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-unused-function")
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-reorder")
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-sign-compare")
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-unknown-pragmas")

#############################################################################################################
# linker options
#############################################################################################################
list(APPEND VF_LINK_OPTIONS -lgomp)
list(APPEND VF_LINK_OPTIONS -lrt)
list(APPEND VF_LINK_OPTIONS -ldl)
