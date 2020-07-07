macro(linkMetis targetName)

	include_directories(${METIS_ROOT}/include)
	target_link_libraries(${targetName} metis)

endmacro(linkMetis)
