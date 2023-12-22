#=======================================================================================
# ____          ____    __    ______     __________   __      __       __        __
# \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
#  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
#   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
#    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
#     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
#      \    \  |    |   ________________________________________________________________
#       \    \ |    |  |  ______________________________________________________________|
#        \    \|    |  |  |         __          __     __     __     ______      _______
#         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
#          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
#           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
#            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
#
#  This file is part of VirtualFluids. VirtualFluids is free software: you can
#  redistribute it and/or modify it under the terms of the GNU General Public
#  License as published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
#
#  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#  for more details.
#
#  SPDX-License-Identifier: GPL-3.0-or-later
#  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
#
#! \author Soeren Peters
#=======================================================================================

# public function collectFiles
# Return: <MY_SRC> writes all found source files into 
# Input: <ARG_FILES> are added to source files
#        <ARG_FOLDER> all files in folder are added to source files
#        <ARG_EXCLUDE> removes files from source files
#        When ARG_FILES and ARG_FOLDER are not defined all files are taken from directory
macro(collectFiles ARG_FILES ARG_FOLDER ARG_EXCLUDE)
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

    # set output variable <MY_SRC>
    set(MY_SRCS ${local_source_files})

    # create source groups for Visual Studio (https://cmake.org/cmake/help/latest/command/source_group.html)
    createSourceGroups("${MY_SRCS}")
endmacro()


function(createSourceGroups file_paths)

    foreach(file_path ${file_paths})

        get_filename_component(path ${file_path} DIRECTORY)

        if (path)
            source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} PREFIX "general" FILES ${file_path})
        endif()

    endforeach()

endfunction(createSourceGroups)
