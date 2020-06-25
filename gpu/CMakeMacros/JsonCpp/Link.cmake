macro(linkJsonCpp targetName)

	target_include_directories(${targetName} PUBLIC ${JSONCPP_ROOT}/include)
	target_link_libraries(${targetName} jsoncpp)

endmacro(linkJsonCpp)
