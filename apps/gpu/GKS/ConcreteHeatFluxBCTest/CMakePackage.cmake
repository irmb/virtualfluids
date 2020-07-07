#FILE ENDINGS
resetFileEndingsToCollect()
addCAndCPPFileTypes()
addFileEndingToCollect("*.cu")
addFileEndingToCollect("*.cuh")

#GLOB SOURCE FILES IN MY_SRCS
unset(MY_SRCS)
includeRecursiveAllFilesFrom(${targetName} ${CMAKE_CURRENT_LIST_DIR})