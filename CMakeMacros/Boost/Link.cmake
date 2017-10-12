macro(linkBoost targetName components)
  set(Boost_USE_STATIC_LIBS OFF)
  set(Boost_USE_MULTITHREADED ON)
  set(Boost_USE_STATIC_RUNTIME OFF)

  if (WIN32)
	add_definitions( -DBOOST_ALL_NO_LIB )
	add_definitions( -DBOOST_ALL_DYN_LINK )
  endif()
  
  FIND_PACKAGE( Boost REQUIRED COMPONENTS ${components})
  INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
  TARGET_LINK_LIBRARIES( ${targetName} ${Boost_LIBRARIES} )
endmacro(linkBoost)
