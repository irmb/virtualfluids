###############################################################################################################
##  MSVC
###############################################################################################################

###############################################################################################################
## USE_UNSECURE_STL_VECTORS_RELEASE ?
###############################################################################################################
OPTION(USE_UNSECURE_STL_VECTORS_RELEASE "_SECURE_SCL=0" OFF)
IF(USE_UNSECURE_STL_VECTORS_RELEASE)
    # More MSVC specific compilation flags
    LIST(APPEND VF_COMPILER_DEFINITION _SECURE_SCL=0)
    LIST(APPEND VF_COMPILER_DEFINITION _SCL_SECURE_NO_WARNINGS)
ENDIF()

#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_RELEASE "/W1")
#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_DEBUG "/Wall")

###############################################################################################################
## Flags
###############################################################################################################
LIST(APPEND VF_COMPILER_DEFINITION _CRT_SECURE_NO_DEPRECATE) # disable warnings promoting Microsoft's security enhanced CRT
#LIST(APPEND VF_COMPILER_DEFINITION _SCL_SECURE_NO_WARNINGS)  # disable warnings triggered by Microsoft's checked iterators
#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/wd4996") # deprecated strcpy...
#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/wd4800") # forcing value to bool 'true' or 'false' (performance warning)

LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/bigobj") # increases that address capacity to 4,294,967,296 (2^32).
LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-MP") # enable multi-threaded compiling

LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/wd4251") # disable needs to have dll interface