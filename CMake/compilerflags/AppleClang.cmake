#############################################################################################################
# compiler flags
#############################################################################################################

# debug
list(APPEND CS_COMPILER_FLAGS_CXX_DEBUG "-g")  # generates debug information. Works best with -O0.
list(APPEND CS_COMPILER_FLAGS_CXX_DEBUG "-O0")

# release
list(APPEND CS_COMPILER_FLAGS_CXX_RELEASE "-O3") # optimization level (-O3: most optimization which also could result in larger binaries)

# all
list(APPEND CS_COMPILER_FLAGS_CXX "-fPIC") # position independent code for shared libraries

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
