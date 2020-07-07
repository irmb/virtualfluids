macro(includeRecursiveAllFilesFrom targetName path)
	set(collectTestFiles ON)
	set(collectProductionFiles ON)

	includeRecursiveFilesFrom(${targetName} ${path})
endmacro(includeRecursiveAllFilesFrom)



macro(includeRecursiveProductionFilesFrom targetName path)
	set(collectTestFiles OFF)
	set(collectProductionFiles ON)

	includeRecursiveFilesFrom(${targetName} ${path})
endmacro(includeRecursiveProductionFilesFrom)



macro(includeRecursiveTestFilesFrom targetName path)
	set(collectTestFiles ON)
	set(collectProductionFiles OFF)

	includeRecursiveFilesFrom(${targetName} ${path})
endmacro(includeRecursiveTestFilesFrom)




macro(includeRecursiveFilesFrom targetName path)
	file(GLOB_RECURSE includeSourcePaths ${path}/package.include)

	foreach(package ${includeSourcePaths})
		get_filename_component(package_dir ${package} DIRECTORY)
		collectFilesFrom(${package_dir} "${FILES_TO_COLLECT}")
		setSourceGroupForFilesIn(${package_dir} ${targetName})
	endforeach()
endmacro(includeRecursiveFilesFrom)



macro(collectFilesFrom path file_endings)
	#input: path from files to collect
	unset(COLLECTED_FILES_IN_PATH)
	unset(ALL_FILES_IN_PATH)

	foreach(_ending ${file_endings})
		FILE(GLOB filesWithEnding ${path}/${_ending})
		set(ALL_FILES_IN_PATH ${ALL_FILES_IN_PATH} ${filesWithEnding})
	endforeach()

	foreach(_file ${ALL_FILES_IN_PATH})
		get_filename_component(fileName ${_file} NAME)
		if(collectTestFiles)
			if(${fileName} MATCHES "Test" OR ${fileName} MATCHES "Mock")
				set(COLLECTED_FILES_IN_PATH ${COLLECTED_FILES_IN_PATH} ${_file})
			endif()
		endif()
		if(collectProductionFiles)
			if(NOT ${fileName} MATCHES "Test" AND NOT ${fileName} MATCHES "Mock")
				set(COLLECTED_FILES_IN_PATH ${COLLECTED_FILES_IN_PATH} ${_file})
			endif()
		endif()
	endforeach()
	set(MY_SRCS ${MY_SRCS} ${COLLECTED_FILES_IN_PATH})

	#output: MY_SRCS
endmacro(collectFilesFrom)




macro(setSourceGroupForFilesIn package_dir targetName)
#input: target_name PACKAGE_SRCS
	buildSourceGroup(${targetName} ${package_dir})

	if(isAllTestSuite)
		source_group(${targetName}\\${SOURCE_GROUP} FILES ${COLLECTED_FILES_IN_PATH})
	else()
		source_group(${SOURCE_GROUP} FILES ${COLLECTED_FILES_IN_PATH})
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

	if(NOT SOURCE_GROUP)
		set(SOURCE_GROUP "general")
	endif()

#output: SOURCE_GROUP
endmacro(buildSourceGroup)