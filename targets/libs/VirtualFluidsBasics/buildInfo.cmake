set(buildInfoPath ${CMAKE_BINARY_DIR}/buildInfo)
set(buildInfoFile buildInfo.cpp)
set(buildInfoFileHeader buildInfo.h)
set(buildInfoInput ${CMAKE_CURRENT_LIST_DIR}/buildInfoInput/buildInfo.in.cpp)

include(${CMAKE_SOURCE_DIR}/${cmakeMacroPath}/git/GetGitRevisionDescription.cmake)
get_git_head_revision(git_branch git_commit_hash)

set(COMPILER_FLAGS "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}} ${CMAKE_CXX_FLAGS}")

configure_file(${buildInfoInput} ${buildInfoPath}/${buildInfoFile})

include_directories(${buildInfoPath})
set(MY_SRCS ${MY_SRCS} ${buildInfoPath}/${buildInfoFile} ${CMAKE_CURRENT_LIST_DIR}/${buildInfoFileHeader})
source_group("" FILES  ${buildInfoPath}/${buildInfoFile})
source_group("" FILES  ${CMAKE_CURRENT_LIST_DIR}/${buildInfoFileHeader})