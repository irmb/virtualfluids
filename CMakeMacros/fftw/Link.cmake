
macro(linkFftw targetName)

	include_directories(${FFTW_ROOT}/api)
	target_link_libraries(${targetName} fftw3)

endmacro(linkFftw)
