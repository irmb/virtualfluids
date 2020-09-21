###############################################################################################################
##
##  apple clang
##
###############################################################################################################

#############################################################################################################
# Flags
#############################################################################################################
LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-O3;-fomit-frame-pointer;-finline-functions;-fPIC;-Wbackslash-newline-escape")

#############################################################################################################
# mt support
#############################################################################################################
LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-pthread")

#############################################################################################################
# c++ 11 support
#############################################################################################################
LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-std=c++11")


# test
#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_RELEASE "-Wall")
#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS_DEBUG "-Werror")

#############################################################################################################
# disable warning
#############################################################################################################
LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wno-deprecated") #deprecated header warning
#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wbackslash-newline-escape") #backslash and newline separated by space
LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wcomment") #'/*' within block comment

#############################################################################################################
# c++ 17 support
#############################################################################################################
#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-std=c++17")
#LIST(APPEND CAB_COMPILER_ADDTIONAL_C_COMPILER_FLAGS "-std=c++17")

#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-fext-numeric-literals")
#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=0")
#LIST(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-Wregister")
