#FILE ENDINGS
resetFileEndingsToCollect()
addCAndCPPFileTypes()

#GLOB SOURCE FILES IN MY_SRCS
unset(MY_SRCS)
includeRecursiveAllFilesFrom(${targetName} ${CMAKE_CURRENT_LIST_DIR})
includeRecursiveTestFilesFrom(VirtualFluidsBasics ${CMAKE_SOURCE_DIR}/src/core)
