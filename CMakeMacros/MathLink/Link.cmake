
macro(linkMathLink targetName)

	include_directories(${MATHLINK_ROOT})
	target_link_libraries(${targetName} ${MATHLINK_ROOT}\\wstp64i4.lib)

endmacro(linkMathLink)