macro(sharedLibs)
	if(${BUILD_SHARED_LIBS})
		set(LIB_TYPE SHARED)
	else()
		set(LIB_TYPE STATIC)

		set(CompilerFlags
				CMAKE_CXX_FLAGS
				CMAKE_CXX_FLAGS_DEBUG
				CMAKE_CXX_FLAGS_RELEASE
				CMAKE_C_FLAGS
				CMAKE_C_FLAGS_DEBUG
				CMAKE_C_FLAGS_RELEASE
				)
		foreach(CompilerFlag ${CompilerFlags})
			string(REPLACE "/MD" "/MT" ${CompilerFlag} "${${CompilerFlag}}")
		endforeach()
        
		set(CompilerFlagsCuda
                CMAKE_CUDA_FLAGS
                CMAKE_CUDA_FLAGS_DEBUG
                CMAKE_CUDA_FLAGS_MINSIZEREL
                CMAKE_CUDA_FLAGS_RELEASE
                CMAKE_CUDA_FLAGS_RELWITHDEBINFO
				)
		foreach(CompilerFlag ${CompilerFlagsCuda})
			string(REPLACE "-MD" "-MT" ${CompilerFlag} "${${CompilerFlag}}")
		endforeach()
	endif()
endmacro(sharedLibs)



########################################################################
#                            AllTest Option                            #
########################################################################

macro(activateAllTestOption)
	set(isAllTestSuite ON)
endmacro(activateAllTestOption)


macro(deactivateAllTestOption)
	set(isAllTestSuite OFF)
endmacro(deactivateAllTestOption)


########################################################################
#                           target name                                #
########################################################################

macro(setTargetNameToFolderName path)
     get_filename_component(targetName "${path}" NAME)
endmacro(setTargetNameToFolderName)