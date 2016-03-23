###############################################################################################################
## 
##  gcc40
##
###############################################################################################################

MACRO(SET_COMPILER_SPECIFIC_FLAGS_INTERN build_type use64BitOptions)

   ###############################################################################################################
   ## Flags
   ###############################################################################################################
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-O3 -fomit-frame-pointer -finline-functions -funroll-all-loops -fPIC")
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-deprecated") #deprecated header warning (jarl benutzt sstream weil schneller und so) 

   LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "-O3 -fomit-frame-pointer -finline-functions -funroll-all-loops -fPIC")

   IF(CPU_TYPE MATCHES "Opteron")
     LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-mtune=opteron")
     LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-mcpu=opteron" )

     LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "-mtune=opteron")
     LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "-mcpu=opteron" )
   ENDIF()


   ###############################################################################################################
   ## 64Bit support
   ###############################################################################################################
   IF( ${use64BitOptions} ) 
     LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-m64" )
     LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS   "-m64" )
   ENDIF()


   ###############################################################################################################
   ## OpenMP support
   ###############################################################################################################
   IF(USE_OPENMP)
   	MESSAGE(STATUS "gcc40 has no OpenMP support -> OpenMP deactivated")
   	SET(USE_OPENMP "OFF")
   ENDIF()


   ###############################################################################################################
   ## mt support
   ###############################################################################################################
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-pthread")
   LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "-pthread")

   IF(NOT APPLE)
      LIST(APPEND CAB_ADDITIONAL_LINK_PROPS "-lrt")
   ENDIF()

ENDMACRO(SET_COMPILER_SPECIFIC_FLAGS_INTERN build_type use64BitOptions)
