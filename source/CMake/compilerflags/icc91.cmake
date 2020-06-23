###############################################################################################################
## 
##  intel9 
##
###############################################################################################################

MACRO(SET_COMPILER_SPECIFIC_FLAGS_INTERN build_type use64BitOptions)

   IF( ${use64BitOptions} )
     LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-D__amd64" ) 
   ENDIF()

   #IF( ${compiler} MATCHES "INTEL_91_IPO")
   #	LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "-ipo" ) 
   #	LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-ipo" ) 
   #	LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-opt-mem-bandwidth2" ) 
   #ENDIF( ${compiler} MATCHES "INTEL_91_IPO")

   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-O3 -fast -fomit-frame-pointer -finline-functions")

   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-wd68" ) #integer conversion resulted in a change of sign
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-wd186") # pointless comparison of unsigned integer with zero
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-wd654")
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-wd654")
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-wd1125") #virtual function override intended
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-wd1224") #warning directive: This file includes at least one deprecated or antiquated header
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-wd377")  #class "std::auto_ptr<RCF::I_ClientTransport>" has no suitable copy constructor
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-wd327")  #class "std::auto_ptr<RCF::I_ClientTransport>" has no suitable copy constructor
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-wd327")  #class "std::auto_ptr<RCF::I_ClientTransport>" has no suitable copy constructor

   LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "-wd266")  #function "__GKfree" declared implicitly
   LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "-O3 -fast -fomit-frame-pointer -finline-functions")

   ###############################################################################################################
   ## OpenMP support
   ###############################################################################################################
   IF(USE_OPENMP)
   	LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-openmp")
   ENDIF()

   ###############################################################################################################
   ## mt support
   ###############################################################################################################
   LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-pthread")
   LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "-pthread")

ENDMACRO(SET_COMPILER_SPECIFIC_FLAGS_INTERN build_type use64BitOptions)
