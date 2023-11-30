#################################################################################
## Helper functions for building source groups
## and extracting test/production files.
##
## After function call the files are stored in: MY_SRCS
#################################################################################

macro(includeAllFiles folderName file_path)
	if(NOT DEFINED collectTestFiles)
	    set(collectTestFiles ON)
	endif()
	
	if(NOT DEFINED collectProductionFiles)
        set(collectProductionFiles ON)
    endif()

	includeFiles(${folderName} "${file_path}")
endmacro(includeAllFiles)


macro(includeProductionFiles folderName file_path)
	if(NOT DEFINED collectTestFiles)
	    set(collectTestFiles OFF)
	endif()
	
	if(NOT DEFINED collectProductionFiles)
        set(collectProductionFiles ON)
    endif()

	includeFiles(${folderName} "${file_path}")
endmacro(includeProductionFiles)



macro(includeTestFiles folderName file_paths)
	if(NOT DEFINED collectTestFiles)
		set(collectTestFiles ON)
	endif()

	if(NOT DEFINED collectProductionFiles)
		set(collectProductionFiles OFF)
	endif()

	includeFiles(${folderName} "${file_paths}")
endmacro(includeTestFiles)




macro(includeFiles folderName file_paths)

	foreach(file ${file_paths})

		get_filename_component(package_dir ${file} DIRECTORY)
		#message("File: " ${file})
		#message("package_dir: " ${package_dir})

		collectFilesFrom(${file})
		if (package_dir)
		   setSourceGroupForFilesIn(${file} ${package_dir} ${folderName})
		endif()

	endforeach()

	unset(collectTestFiles)
	unset(collectProductionFiles)

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




macro(setSourceGroupForFilesIn file package_dir folderName)
#input: target_name PACKAGE_SRCS
	buildSourceGroup(${folderName} ${package_dir})
	source_group(${SOURCE_GROUP} FILES ${file})
#output: -
endmacro(setSourceGroupForFilesIn)




macro(buildSourceGroup folderName path)
#input: folderName (e.g. name of folder after src/)

	unset(SOURCE_GROUP)
	string(REPLACE "/" ";" folderListFromPath ${path})
	set(findFolderName 0)

	foreach(folder ${folderListFromPath})
		if(findFolderName)
			set(SOURCE_GROUP ${SOURCE_GROUP}\\${folder})
		endif()

		if(${folder} STREQUAL ${folderName})
			SET(findFolderName 1)
		endif()
	endforeach()

	#message("SOURCE_GROUP: " ${SOURCE_GROUP})

	if(NOT SOURCE_GROUP)
		set(SOURCE_GROUP "general")
	endif()

#output: SOURCE_GROUP
endmacro(buildSourceGroup)


function(collectFiles source_files ARG_FILES ARG_FOLDER ARG_EXCLUDE)
	set(local_source_files)

	#cmake_print_variables(ARG_FOLDER)
	#cmake_print_variables(ARG_FILES)
	#cmake_print_variables(ARG_EXCLUDE)

	if (ARG_FILES)
		set(local_source_files ${local_source_files} ${ARG_FILES})
	endif()

	if (ARG_FOLDER)
		foreach(folder ${ARG_FOLDER})
			foreach(file ${VIRTUAL_FLUIDS_GLOB_FILES})
				set (filePath ${folder}/${file})
				#message("${filePath}")
				file (GLOB part_files ${filePath} )
				set(local_source_files ${local_source_files} ${part_files})
				#message("${local_source_files}")
			endforeach()
		endforeach()
	endif()

	if (NOT ARG_FILES AND NOT ARG_FOLDER)
		file ( GLOB_RECURSE all_files ${VIRTUAL_FLUIDS_GLOB_FILES} )
		set(local_source_files ${local_source_files} ${all_files})
	endif()
	
	if (ARG_EXCLUDE)
		foreach(file_path ${local_source_files})
		    set(exclude_file OFF)
			foreach(file_exclude ${ARG_EXCLUDE})
				get_filename_component(file_name ${file_path} NAME)
				if (${file_name} STREQUAL ${file_exclude})
					set(exclude_file TRUE)
				endif()
			endforeach()
			if (NOT ${exclude_file})	
				set(new_files ${new_files} ${file_path})
			endif()
		endforeach()
		set(local_source_files ${new_files})
	endif()

	set("${source_files}" "${local_source_files}" PARENT_SCOPE)
endfunction()