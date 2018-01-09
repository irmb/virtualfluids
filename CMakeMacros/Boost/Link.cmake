macro(linkBoost targetName components)
  set(Boost_USE_STATIC_LIBS ON)
  set(Boost_USE_MULTITHREADED ON)
  set(Boost_USE_STATIC_RUNTIME OFF)

  if (WIN32)
	add_definitions( -DBOOST_ALL_NO_LIB )
#	add_definitions( -DBOOST_ALL_DYN_LINK )
  endif()
  
  find_package( Boost REQUIRED COMPONENTS ${components})
  target_include_directories(${targetName} PRIVATE ${Boost_INCLUDE_DIR})
  target_link_libraries(${targetName} ${Boost_LIBRARIES})
endmacro(linkBoost)
