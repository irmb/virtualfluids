#FILE ENDINGS
resetFileEndingsToCollect()
addCAndCPPFileTypes()

#GLOB SOURCE FILES IN MY_SRCS
unset(MY_SRCS)
includeRecursiveAllFilesFrom(${targetName} ${CMAKE_CURRENT_LIST_DIR}/src)
includeRecursiveAllFilesFrom(${targetName} ${CMAKE_CURRENT_LIST_DIR}/include)