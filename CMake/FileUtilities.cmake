#################################################################################
## Helper functions for building source groups
## and extracting test/production files.
##
## After function call the files are stored in: MY_SRCS
#################################################################################

macro(includeAllFiles targetName file_path)
	set(collectTestFiles ON)
	set(collectProductionFiles ON)

	includeFiles(${targetName} "${file_path}")
endmacro(includeAllFiles)


macro(includeProductionFiles targetName file_path)
	set(collectTestFiles OFF)
	set(collectProductionFiles ON)

	includeFiles(${targetName} "${file_path}")
endmacro(includeProductionFiles)



macro(includeTestFiles targetName file_paths)
	set(collectTestFiles ON)
	set(collectProductionFiles OFF)

	includeFiles(${targetName} "${file_paths}")
endmacro(includeTestFiles)




macro(includeFiles targetName file_paths)

	foreach(file ${file_paths})

		get_filename_component(package_dir ${file} DIRECTORY)
		#message("File: " ${file})
		#message("package_dir: " ${package_dir})

		collectFilesFrom(${file})
		if (package_dir)
		   setSourceGroupForFilesIn(${package_dir} ${targetName})
		endif()

	endforeach()

endmacro(includeFiles)



macro(collectFilesFrom path)
	#input: path from files to collect

	get_filename_component(fileName ${path} NAME)
	if(collectTestFiles)
		if(${fileName} MATCHES "Test" OR ${fileName} MATCHES "Mock")
			set(MY_SRCS ${MY_SRCS} ${path})
		endif()
	endif()
	if(collectProductionFiles)
		if(NOT ${fileName} MATCHES "Test" AND NOT ${fileName} MATCHES "Mock")
			set(MY_SRCS ${MY_SRCS} ${path})
		endif()
	endif()

	#output: MY_SRCS
endmacro()




macro(setSourceGroupForFilesIn package_dir targetName)
#input: target_name PACKAGE_SRCS
	buildSourceGroup(${targetName} ${package_dir})

	if(isAllTestSuite)
		source_group(${targetName}\\${SOURCE_GROUP} FILES ${MY_SRCS})
	else()
		source_group(${SOURCE_GROUP} FILES ${MY_SRCS})
	endif()
#output: -
endmacro(setSourceGroupForFilesIn)




macro(buildSourceGroup targetName path)
#input: targetName (e.g. lib name, exe name)

	unset(SOURCE_GROUP)
	string(REPLACE "/" ";" folderListFromPath ${path})
	set(findTargetName 0)

	foreach(folder ${folderListFromPath})
		if(findTargetName)
			set(SOURCE_GROUP ${SOURCE_GROUP}\\${folder})
		endif()

		if(${folder} STREQUAL ${targetName})
			SET(findTargetName 1)
		endif()
	endforeach()

	#message("SOURCE_GROUP: " ${SOURCE_GROUP})

	if(NOT SOURCE_GROUP)
		set(SOURCE_GROUP "general")
	endif()

#output: SOURCE_GROUP
endmacro(buildSourceGroup)


include (GenerateExportHeader)
macro(generateExportHeader libName)
	#if(${BUILD_SHARED_LIBS})
		GENERATE_EXPORT_HEADER	(${libName}
				#BASE_NAME ${libName}
				#EXPORT_MACRO_NAME ${libName}_EXPORT
				EXPORT_FILE_NAME ${CMAKE_BINARY_DIR}/${libName}_export.h
				#STATIC_DEFINE ${libName}_BUILT_AS_STATIC
				)
	#endif()
endmacro(generateExportHeader)