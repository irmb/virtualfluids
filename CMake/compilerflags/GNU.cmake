#############################################################################################################
# compiler flags
#############################################################################################################

# debug
list(APPEND CS_COMPILER_FLAGS_CXX_DEBUG "-g")  # generates debug information. Works best with -O0.
list(APPEND CS_COMPILER_FLAGS_CXX_DEBUG "-O0") # no optimization

# release
list(APPEND CS_COMPILER_FLAGS_CXX_RELEASE "-O3") # optimization level (-O3: most optimization which also could result in larger binaries)

# all
list(APPEND CS_COMPILER_FLAGS_CXX "-fPIC") # position independent code for shared libraries

if(NOT BUILD_VF_INCLUDE_WHAT_YOU_USE) # optimization flag '-funroll-all-loops' is not supported for IWYU
    LIST(APPEND CS_COMPILER_FLAGS_CXX "-funroll-all-loops")
endif()

# gcov
if (BUILD_VF_COVERAGE)
    list(APPEND CS_COMPILER_FLAGS_CXX "--coverage")
    set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} " --coverage")
endif()

#############################################################################################################
# warnings
#############################################################################################################
list(APPEND CS_COMPILER_FLAGS_CXX "-Wall")
list(APPEND CS_COMPILER_FLAGS_CXX "-Wextra")
list(APPEND CS_COMPILER_FLAGS_CXX "-pedantic")

if(BUILD_WARNINGS_AS_ERRORS)
    list(APPEND CS_COMPILER_FLAGS_CXX -Werror)
endif()

list(APPEND CS_COMPILER_FLAGS_CXX "-Wno-unused-function")
list(APPEND CS_COMPILER_FLAGS_CXX "-Wno-unused-parameter")
list(APPEND CS_COMPILER_FLAGS_CXX "-Wno-reorder")
list(APPEND CS_COMPILER_FLAGS_CXX "-Wno-unknown-pragmas")
list(APPEND CS_COMPILER_FLAGS_CXX "-Wno-cast-function-type")

#############################################################################################################
# linker options
#############################################################################################################
list(APPEND VF_LINK_OPTIONS -lgomp)
list(APPEND VF_LINK_OPTIONS -lrt)
list(APPEND VF_LINK_OPTIONS -ldl)
