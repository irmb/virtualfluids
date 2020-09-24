
#############################################################################################################
# Flags
#############################################################################################################
LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-O3;-fomit-frame-pointer;-finline-functions;-funroll-all-loops;-fPIC")
LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-deprecated") #deprecated header warning (jarl benutzt sstream weil schneller und so)

#############################################################################################################
# mt support
#############################################################################################################
LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-pthread")



#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-fext-numeric-literals")
#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=0")
#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wregister")


list(APPEND VF_LINK_OPTIONS -lgomp)
list(APPEND VF_LINK_OPTIONS -lrt)

list(APPEND VF_LINK_OPTIONS -ldl)
