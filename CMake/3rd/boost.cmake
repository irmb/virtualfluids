function(linkBoost)

    set( options )
    set( oneValueArgs )
    set( multiValueArgs COMPONENTS)
    cmake_parse_arguments( ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  if(BUILD_SHARED_LIBS)
     if (WIN32)
         set(Boost_USE_STATIC_LIBS ON)
     else()
         set(Boost_USE_STATIC_LIBS OFF)
     endif()
	 set(Boost_USE_STATIC_RUNTIME OFF)
  else()
	 set(Boost_USE_STATIC_LIBS ON)
   if(WIN32)
	  set(Boost_USE_STATIC_RUNTIME ON)
   endif()
  endif()
	  
  set(Boost_USE_MULTITHREADED ON)

  if (WIN32)
	add_definitions( -DBOOST_ALL_NO_LIB )
#	add_definitions( -DBOOST_ALL_DYN_LINK )
  endif()

    vf_get_library_name(library_name)
    if(DEFINED ARG_COMPONENTS)
        find_package( Boost REQUIRED COMPONENTS ${ARG_COMPONENTS})
        target_link_libraries(${library_name} PRIVATE ${Boost_LIBRARIES})
    else()
        find_package( Boost REQUIRED)
    endif()


  target_include_directories(${library_name} PRIVATE ${Boost_INCLUDE_DIR})
endfunction()
