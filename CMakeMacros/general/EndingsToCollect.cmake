macro(addFileEndingToCollect file_ending)
	#input: file_ending --> appends it to the list of files to collect

	list (FIND FILES_TO_COLLECT ${file_ending} index)
	if (${index} EQUAL -1)
		set(FILES_TO_COLLECT ${FILES_TO_COLLECT} ${file_ending})
	endif()

	#output: files_to_collect
endmacro(addFileEndingToCollect)



macro(resetFileEndingsToCollect)
	unset(FILES_TO_COLLECT)
endmacro(resetFileEndingsToCollect)




macro(addCAndCPPFileTypes)
	addFileEndingToCollect("*.h")
	addFileEndingToCollect("*.c")
	addFileEndingToCollect("*.cpp")
	addFileEndingToCollect("*.cxx")
	addFileEndingToCollect("*.hpp")
	addFileEndingToCollect("*.cu")
	addFileEndingToCollect("*.cuh")
endmacro(addCAndCPPFileTypes)


macro(addObjCAndObjCPPFileTypesToCollect)
	addFileEndingToCollect("*.m")
	addFileEndingToCollect("*.mm")
	addFileEndingToCollect("*.h")
endmacro(addObjCAndObjCPPFileTypesToCollect)
