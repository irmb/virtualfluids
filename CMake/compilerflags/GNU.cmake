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

if(NOT BUILD_VF_INCLUDE_WHAT_YOU_USE) # optimization flag '-funroll-all-loops' is not supported for IWYU
    LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-funroll-all-loops")
endif()

# gcov
if (BUILD_VF_COVERAGE)
    list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "--coverage")
    set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} " --coverage")
endif()

#############################################################################################################
# warnings
#############################################################################################################
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wall")
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wextra")
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-pedantic")

if(BUILD_WARNINGS_AS_ERRORS)
    list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS -Werror)
endif()

list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-unused-function")
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-unused-parameter")
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-reorder")
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-unknown-pragmas")
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-cast-function-type")

#############################################################################################################
# linker options
#############################################################################################################
list(APPEND VF_LINK_OPTIONS -lgomp)
list(APPEND VF_LINK_OPTIONS -lrt)
list(APPEND VF_LINK_OPTIONS -ldl)
