include(${CMAKE_CURRENT_LIST_DIR}/FindFftw.cmake)

macro(linkFftw targetName)

	include_directories(${FFTW_INCLUDE_DIRS})
	target_link_libraries(${targetName} ${FFTW_LIBRARIES})

endmacro(linkFftw)
