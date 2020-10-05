#############################################################################################################
# compiler flags
#############################################################################################################

#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-O")
#~ LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-wd654")
#~ LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-wd1125") #virtual function override intended
#~ LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-wd1224") #warning directive: This file includes at least one deprecated or antiquated header
#~ LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-wd377")  #class "std::auto_ptr<RCF::I_ClientTransport>" has no suitable copy constructor
#~ LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-wd327")  #class "std::auto_ptr<RCF::I_ClientTransport>" has no suitable copy constructor
#~ LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-wd327")  #class "std::auto_ptr<RCF::I_ClientTransport>" has no suitable copy constructor
#~
#~ LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "-wd266")  #function "__GKfree" declared implicitly
#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-xHOST -O3 -ip -ipo -fno-alias -mcmodel=medium -qopt-streaming-stores=always")

# all
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-xHOST;-O3;-ip;-fno-alias;-mcmodel=medium;-qopt-streaming-stores=always;-xCORE-AVX512;-qopt-zmm-usage=high")

# debug
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_DEBUG "-g -traceback")


#############################################################################################################
# linker options
#############################################################################################################
list(APPEND VF_LINK_OPTIONS -parallel)

#############################################################################################################
# preprocessor definitions
#############################################################################################################
# LIST(APPEND VF_COMPILER_DEFINITION MPICH_IGNORE_CXX_SEEK)
# LIST(APPEND VF_COMPILER_DEFINITION MPICH_SKIP_MPICXX)
