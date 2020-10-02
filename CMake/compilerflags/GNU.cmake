
#############################################################################################################
# Flags
#############################################################################################################
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-g")  # generates debug information. Works best with -O0.
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-O0") # optimization level (-O3: most optimization which also could result in larger binaries)
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-fPIC") # position independent code for shared libraries
LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-funroll-all-loops")

list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wall")
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-unused-function")
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-reorder")
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-sign-compare")


#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-fext-numeric-literals")
#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=0")

list(APPEND VF_LINK_OPTIONS -lgomp)
list(APPEND VF_LINK_OPTIONS -lrt)

list(APPEND VF_LINK_OPTIONS -ldl)
