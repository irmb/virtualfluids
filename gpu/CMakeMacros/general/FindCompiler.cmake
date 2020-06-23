#############################################################
###                COMPILER DETECTION                     ###
#############################################################
# Check for intel compiler
if( CMAKE_CXX_COMPILER MATCHES "icpc" OR CMAKE_CXX_COMPILER_ARG1 MATCHES "icpc" )
    option ( VF_CXX_COMPILER_IS_INTEL "Use Intel compiler" ON  )
    # Intel(R) Compiler has its own library archiver,
    # if you build libraries and do not use xiar,
    # the Intel compiler will complain about invalid
    # archives at the link phase.
    # The Intel(R) archiver is "xiar" usually
    # located in the same folder as the compiler,
    FIND_PROGRAM(XIAR xiar)
    IF(XIAR)
        SET(CMAKE_AR "${XIAR}")
    ENDIF(XIAR)
    MARK_AS_ADVANCED(XIAR)

    # Intel(R) Compiler also comes with its own linker
    # which provides a number of additional benefits when
    # linking code compiled with the Intel(R) compiler.
    # Again, usually in the same place as icc itself,
    FIND_PROGRAM(XILD xild)
    IF(XILD)
       SET(CMAKE_LINKER "${XILD}")
    ENDIF(XILD)
    MARK_AS_ADVANCED(XILD)
else()
    option ( VF_CXX_COMPILER_IS_INTEL "Use Intel compiler" OFF  )
endif()
mark_as_advanced ( VF_CXX_COMPILER_IS_INTEL )

# Check for Gnu compiler
if ( CMAKE_COMPILER_IS_GNUCXX  AND NOT VF_CXX_COMPILER_IS_INTEL )
     option ( VF_CXX_COMPILER_IS_GNU "Use gnu compiler" ON  )
else()
     option ( VF_CXX_COMPILER_IS_GNU "Use gnu compiler" OFF  )
endif()
mark_as_advanced ( VF_CXX_COMPILER_IS_GNU )

# Check for Visual Studio
if(MSVC)
     option (VF_CXX_COMPILER_IS_MSVC "Use Visual Studio compiler" ON)
else()
     option (VF_CXX_COMPILER_IS_MSVC "Use Visual Studio compiler" OFF)
endif()
mark_as_advanced(VF_CXX_COMPILER_IS_MSVC)

# Check for IBM compiler
if( CMAKE_CXX_COMPILER MATCHES "xlc" OR CMAKE_CXX_COMPILER_ARG1 MATCHES "xlc" )
    option ( VF_CXX_COMPILER_IS_IBM "Use IBM compiler" ON  )
else()
    option ( VF_CXX_COMPILER_IS_IBM "Use IBM compiler" OFF  )
endif()
mark_as_advanced ( VF_CXX_COMPILER_IS_IBM )

# Check for NEC SX compiler
if( CMAKE_CXX_COMPILER MATCHES "sxc" OR CMAKE_CXX_COMPILER_ARG1 MATCHES "sxc" OR CMAKE_CXX_COMPILER MATCHES "sxmpic" OR CMAKE_CXX_COMPILER_ARG1 MATCHES "sxmpic" )
    option ( VF_CXX_COMPILER_IS_NEC "Use NEC compiler" ON  )
else()
    option ( VF_CXX_COMPILER_IS_NEC "Use NEC compiler" OFF  )
endif()
mark_as_advanced ( VF_CXX_COMPILER_IS_NEC )

# Check for Clang compiler
if( CMAKE_CXX_COMPILER MATCHES "clang" OR CMAKE_CXX_COMPILER_ARG1 MATCHES "clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang" )
    option ( VF_CXX_COMPILER_IS_CLANG "Use clang compiler" ON  )
else()
    option ( VF_CXX_COMPILER_IS_CLANG "Use clang compiler" OFF  )
endif()
mark_as_advanced ( VF_CXX_COMPILER_IS_CLANG )

if( CMAKE_CXX_COMPILER_ID MATCHES Cray )
    option ( VF_CXX_COMPILER_IS_CRAY "Use Cray compiler" ON   )
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.4)
        message( FATAL_ERROR "Insufficient Cray Compiler Environment version" )
    endif()
else()
    option ( VF_CXX_COMPILER_IS_CRAY "Use Cray compiler" OFF  )
endif()
mark_as_advanced ( VF_CXX_COMPILER_IS_CRAY )

# Check for MPI wrapper
get_filename_component( CXX_COMPILER_WITHOUT_PATH ${CMAKE_CXX_COMPILER} NAME )
if( CXX_COMPILER_WITHOUT_PATH MATCHES "mpi" OR CMAKE_CXX_COMPILER_ARG1 MATCHES "mpi" )
    option ( VF_CXX_COMPILER_IS_MPI_WRAPPER "Compiler is MPI wrapper" ON  )
else()
    option ( VF_CXX_COMPILER_IS_MPI_WRAPPER "Compiler is MPI wrapper" OFF  )
endif()
mark_as_advanced ( VF_CXX_COMPILER_IS_MPI_WRAPPER )

############################################################################################################################
