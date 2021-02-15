#############################################################################################################
# compiler flags
#############################################################################################################

#LIST(APPEND CS_COMPILER_FLAGS_CXX "-O")
#LIST(APPEND CS_COMPILER_FLAGS_CXX "-wd654")
#LIST(APPEND CS_COMPILER_FLAGS_CXX "-wd1125") #virtual function override intended
#LIST(APPEND CS_COMPILER_FLAGS_CXX "-wd1224") #warning directive: This file includes at least one deprecated or antiquated header
#LIST(APPEND CS_COMPILER_FLAGS_CXX "-wd377")  #class "std::auto_ptr<RCF::I_ClientTransport>" has no suitable copy constructor
#LIST(APPEND CS_COMPILER_FLAGS_CXX "-wd327")  #class "std::auto_ptr<RCF::I_ClientTransport>" has no suitable copy constructor
#LIST(APPEND CS_COMPILER_FLAGS_CXX "-wd327")  #class "std::auto_ptr<RCF::I_ClientTransport>" has no suitable copy constructor
#
#LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "-wd266")  #function "__GKfree" declared implicitly
#LIST(APPEND CS_COMPILER_FLAGS_CXX "-xHOST -O3 -ip -ipo -fno-alias -mcmodel=medium -qopt-streaming-stores=always")

# all
list(APPEND CS_COMPILER_FLAGS_CXX "-xHOST;-O3;-ip;-fno-alias;-mcmodel=medium;-qopt-streaming-stores=always;-xCORE-AVX512;-qopt-zmm-usage=high")

# debug
list(APPEND CS_COMPILER_FLAGS_CXX_DEBUG "-g -traceback")


#############################################################################################################
# linker options
#############################################################################################################
list(APPEND VF_LINK_OPTIONS -parallel)

#############################################################################################################
# preprocessor definitions
#############################################################################################################
# LIST(APPEND VF_COMPILER_DEFINITION MPICH_IGNORE_CXX_SEEK)
# LIST(APPEND VF_COMPILER_DEFINITION MPICH_SKIP_MPICXX)
