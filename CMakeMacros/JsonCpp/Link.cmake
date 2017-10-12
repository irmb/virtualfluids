include(${CMAKE_CURRENT_LIST_DIR}/FindJsonCpp.cmake)

macro(linkJsonCpp targetName)

	include_directories(${JSONCPP_INCLUDE_DIRS})
	target_link_libraries(${targetName} ${JSONCPP_LIBRARIES})

endmacro(linkJsonCpp)
