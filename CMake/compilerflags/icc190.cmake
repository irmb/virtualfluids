###############################################################################################################
## 
##  intel190 
##
###############################################################################################################

MACRO(SET_COMPILER_SPECIFIC_FLAGS_INTERN build_type use64BitOptions)
   # auto_ptr warning from mu::Parser
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-wd1478") 

   # main intel compiler flags as CMake-list ("xx;yy;")
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-xHOST;-O3;-ip;-fno-alias;-mcmodel=medium;-qopt-streaming-stores=always;-qopt-zmm-usage=high") 
   
   #Debug
   #LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-g -traceback")
   
   ###############################################################################################################
   ## OpenMP support
   ###############################################################################################################
   IF(USE_OPENMP)
   	LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-qopenmp")
      LIST(APPEND CAB_ADDITIONAL_LINK_FLAGS "-parallel")
   ENDIF()

   ###############################################################################################################
   ## mt support
   ###############################################################################################################
   #LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-pthread")
   #LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "-pthread")
   
   #############################################################################################################
   # c++ 11 support
   #############################################################################################################
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-std=c++11")
   LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "-std=c++11")

ENDMACRO(SET_COMPILER_SPECIFIC_FLAGS_INTERN build_type use64BitOptions)
