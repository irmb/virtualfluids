macro(linkCocoa targetName)
	find_library(COCOA_LIB Cocoa)
	target_link_libraries(${targetName} ${COCOA_LIB})
endmacro(linkCocoa)
