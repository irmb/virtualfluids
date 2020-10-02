#############################################################################################################
# compiler flags
#############################################################################################################
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/bigobj") # increases that address capacity to 4,294,967,296 (2^32).
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-MP")     # enable multi-threaded compiling


#############################################################################################################
# warnings
#############################################################################################################
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/wd4251") # disable needs to have dll interface

#############################################################################################################
# preprocessor definitions
#############################################################################################################
list(APPEND VF_COMPILER_DEFINITION _CRT_SECURE_NO_DEPRECATE) # disable warnings promoting Microsoft's security enhanced CRT

option(USE_UNSECURE_STL_VECTORS_RELEASE "_SECURE_SCL=0" OFF)
if(USE_UNSECURE_STL_VECTORS_RELEASE)
    list(APPEND VF_COMPILER_DEFINITION _SECURE_SCL=0)
    list(APPEND VF_COMPILER_DEFINITION _SCL_SECURE_NO_WARNINGS)
endif()
