###############################################################################################################
##  MSVC
###############################################################################################################

###############################################################################################################
## USE_UNSECURE_STL_VECTORS_RELEASE ?
###############################################################################################################
OPTION(USE_UNSECURE_STL_VECTORS_RELEASE "_SECURE_SCL=0" OFF)
IF(USE_UNSECURE_STL_VECTORS_RELEASE)
    # More MSVC specific compilation flags
    LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-D_SECURE_SCL=0")
    LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-D_SCL_SECURE_NO_WARNINGS")
ENDIF()

###############################################################################################################
## Flags
###############################################################################################################
LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-D_CRT_SECURE_NO_DEPRECATE")
LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/wd4996") #deprecated strcpy...
LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/wd4800") #forcing value to bool 'true' or 'false' (performance warning)
LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/bigobj") #ansonsten funzt mit boost das compilieren unter windows nimmer

###############################################################################################################
## OpenMP support
###############################################################################################################
IF(USE_OPENMP)
    LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/openmp")
ENDIF()
