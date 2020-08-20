###############################################################################################################
## 
##  MSVC 2010 (Version 10.0) 64bit (msvc10_64)
##
###############################################################################################################

MACRO(SET_COMPILER_SPECIFIC_FLAGS_INTERN build_type use64BitOptions)

   IF(NOT ${use64BitOptions})
      MESSAGE(FATAL_ERROR "SET_COMPILER_SPECIFIC_FLAGS: use64BitOptions must be ON for msvc10_64")
   ENDIF()

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

#folgendes kann man mittlerweile weglassen..
#unser project kompilert einwandfrei durch!

#hack (solange CMAke OMPILER_FLAGS_<CONFIG> nicht supported
#foreach(flag_var
#        CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
#        CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
#   if(${flag_var} MATCHES "/MD")
#      string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
#   endif(${flag_var} MATCHES "/MD")
#endforeach(flag_var)

#  IF(${build_type} MATCHES BINARY)
#     LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_DEBUG   "/MTd" )
#     LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_RELEASE "/MT"  )
#  ELSEIF(${build_type} MATCHES STATIC)
#     LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_DEBUG   "/MTd" )
#     LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_RELEASE "/MT"  )
#  ELSEIF(${build_type} MATCHES SHARED)
#     LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_DEBUG   "/MDd" )
#     LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_RELEASE "/MD"  )
#  ELSE()
#     MESSAGE(FATAL_ERROR "build_type=${build_type} doesn't match BINARY, SHARED or STATIC")
#  ENDIF()

   ###############################################################################################################
   ## OpenMP support
   ###############################################################################################################
   IF(USE_OPENMP)
   	LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/openmp")
   ENDIF()

   ###############################################################################################################
   ## boost + rcf extensions
   ###############################################################################################################
   IF(NEED_BOOST)
      IF(BOOST_VERSION)
         STRING(REGEX REPLACE "(.*)\\.(.*)\\.(.*)" "\\2" BoostMinorVersion "${BOOST_VERSION}")
         IF(NOT BoostMinorVersion LESS 35 )     
         
            #IF(NOT ${build_type} MATCHES SHARED)
               #LIST(APPEND CAB_COMPILER_ADDITIONAL_LINK_PROPS "/NODEFAULTLIB:\"MSVCRT\"" )
            #ENDIF()
         ENDIF() 
      ENDIF()
   ENDIF()

ENDMACRO(SET_COMPILER_SPECIFIC_FLAGS_INTERN build_type use64BitOptions)
