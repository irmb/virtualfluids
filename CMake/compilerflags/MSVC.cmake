#############################################################################################################
# compiler flags
#############################################################################################################
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/bigobj") # increases that address capacity to 4,294,967,296 (2^32).
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "-MP")     # enable multi-threaded compiling


#############################################################################################################
# warnings
#############################################################################################################
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/W4") # highest warning level

# With W4 the following warnings appear many times. As long they are not eliminated they are suppressed:
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/wd4458") # C4458: declaration of 'XXX' hides class member
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/wd4100") # C4100: 'XXX': unreferenced formal parameter
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/wd4505") # C4505: 'XXX': unreferenced local function has been removed
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/wd4244") # C4244: '=': conversion from 'int' to 'char', possible loss of data, triggered by algorithm(2216,24)
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/wd4310") # C4310: cast truncates constant value, triggerd by muParserbase.h
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/wd4127") # C4127: conditional expression is constant: e.g. sizeof(int)

# Urgent FIXME: This warning should be activated and fixed:
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/wd4701") # C4701: potentially uninitialized local variable 'lMaxX3' used


list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/wd4251") # disable needs to have dll interface
list(APPEND CAB_COMPILER_ADDTIONAL_CXX_COMPILER_FLAGS "/wd4005") # disable macro redefinition (triggered by metis.h)

#############################################################################################################
# preprocessor definitions
#############################################################################################################
list(APPEND VF_COMPILER_DEFINITION _CRT_SECURE_NO_DEPRECATE) # disable warnings promoting Microsoft's security enhanced CRT

option(USE_UNSECURE_STL_VECTORS_RELEASE "_SECURE_SCL=0" OFF)
if(USE_UNSECURE_STL_VECTORS_RELEASE)
    list(APPEND VF_COMPILER_DEFINITION _SECURE_SCL=0)
    list(APPEND VF_COMPILER_DEFINITION _SCL_SECURE_NO_WARNINGS)
endif()
