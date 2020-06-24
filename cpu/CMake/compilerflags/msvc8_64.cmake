###############################################################################################################
## 
##  MSVC 2005 (Version 8.0) 32bit (msvc8_64)
##
###############################################################################################################
    
MACRO(SET_COMPILER_SPECIFIC_FLAGS_INTERN build_type use64BitOptions)

   IF(NOT ${use64BitOptions})
      MESSAGE(FATAL_ERROR "SET_COMPILER_SPECIFIC_FLAGS: use64BitOptions must be ON for msvc8_64")
   ENDIF()

   ###############################################################################################################
   ## USE_UNSECURE_STL_VECTORS_RELEASE ?
   ###############################################################################################################
   OPTION(USE_UNSECURE_STL_VECTORS_RELEASE "_SECURE_SCL=0" OFF)
   IF(USE_UNSECURE_STL_VECTORS_RELEASE)
      # More MSVC specific compilation flags
      LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-D_SECURE_SCL=0")
      LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-D_SCL_SECURE_NO_WARNINGS")
      LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-D_CRT_SECURE_NO_DEPRECATE")
   ENDIF()

   ###############################################################################################################
   ## Flags
   ###############################################################################################################
   LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D_CRT_SECURE_NO_WARNINGS)
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/wd4996") #deprecated strcpy...
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/wd4800") #forcing value to bool 'true' or 'false' (performance warning)
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/bigobj") #ansonsten funzt mit boost das compilieren unter windows nimmer

#   IF(${build_type} MATCHES BINARY)
#      LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_DEBUG   "/MTd" )
#      LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_RELEASE "/MT"  )
#   ELSEIF(${build_type} MATCHES STATIC)
#      LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_DEBUG   "/MTd" )
#      LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_RELEASE "/MT"  )
#   ELSEIF(${build_type} MATCHES SHARED)
#      LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_DEBUG   "/MDd" )
#      LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_RELEASE "/MD"  )
#   ELSE()
#      MESSAGE(FATAL_ERROR "build_type=${build_type} doesn't match BINARY, SHARED or STATIC")
#   ENDIF()

   ###############################################################################################################
   ## OpenMP support
   ###############################################################################################################
   IF(USE_OPENMP)
   	LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/openmp")
   ENDIF()
ENDMACRO(SET_COMPILER_SPECIFIC_FLAGS_INTERN build_type use64BitOptions)
