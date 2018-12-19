
macro(linkMathLink targetName)

	include_directories(${MATHLINK_ROOT})
	target_link_libraries(${targetName} wstp64i4m)

endmacro(linkMathLink)