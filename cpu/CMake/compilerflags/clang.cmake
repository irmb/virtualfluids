###############################################################################################################
## 
##  clang
##
###############################################################################################################

MACRO(SET_COMPILER_SPECIFIC_FLAGS_INTERN build_type use64BitOptions)
   #############################################################################################################
   # Flags
   #############################################################################################################
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-O3 -fomit-frame-pointer -finline-functions -fPIC -Wbackslash-newline-escape")
 
   LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "-O3 -fomit-frame-pointer -finline-functions -fPIC")

   #############################################################################################################
   # 64Bit support
   #############################################################################################################
   IF( ${use64BitOptions} ) 
     LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-m64" )
     LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS   "-m64" )
   ENDIF()

   #############################################################################################################
   # OpenMP support
   #############################################################################################################
   IF(USE_OPENMP)
     #LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-fopenmp")
     #LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "-fopenmp")
   ENDIF()

   #############################################################################################################
   # mt support
   #############################################################################################################
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-pthread")
   LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "-pthread")

   #############################################################################################################
   # c++ 11 support
   #############################################################################################################
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-std=c++11")
   LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "-std=c++11")

   #############################################################################################################
   # disable warning
   #############################################################################################################
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-deprecated") #deprecated header warning
   #LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wbackslash-newline-escape") #backslash and newline separated by space
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wcomment") #'/*' within block comment

   LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "-Wbackslash-newline-escape") #backslash and newline separated by space
   
   #############################################################################################################
   # c++ 17 support
   #############################################################################################################
   #LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-std=c++17")
   #LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "-std=c++17")
   
   #LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-fext-numeric-literals")
   #LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=0")
   #LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wregister")

   

   IF(NOT APPLE)
      LIST(APPEND CAB_ADDITIONAL_LINK_PROPS "-lrt")
   ENDIF()

ENDMACRO(SET_COMPILER_SPECIFIC_FLAGS_INTERN build_type use64BitOptions)
