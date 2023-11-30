function(set_project_options project_name)

    set(PROJECT_OPTIONS_GCC "")
    set(PROJECT_OPTIONS_GCC_DEBUG "-g;-O0")
    set(PROJECT_OPTIONS_GCC_RELEASE "-O3")

    if(NOT VF_ENABLE_INCLUDE_WHAT_YOU_USE) # optimization flag '-funroll-all-loops' is not supported for IWYU
        LIST(APPEND PROJECT_OPTIONS_GCC "-funroll-all-loops")
    endif()

    set(PROJECT_OPTIONS_CLANG "")
    set(PROJECT_OPTIONS_CLANG_DEBUG "-g;-O0")
    set(PROJECT_OPTIONS_CLANG_RELEASE "-O3")

    set(PROJECT_OPTIONS_MSVC "-MP;/bigobj")
    set(PROJECT_OPTIONS_MSVC_DEBUG "")
    set(PROJECT_OPTIONS_MSVC_RELEASE "")

    set(PROJECT_OPTIONS "")
    set(PROJECT_OPTIONS_DEBUG "")
    set(PROJECT_OPTIONS_RELEASE "")
    if(MSVC)
        target_compile_definitions(${project_name} INTERFACE _CRT_SECURE_NO_DEPRECATE)  # disable warnings promoting Microsoft's security enhanced CRT

        set(PROJECT_OPTIONS ${PROJECT_OPTIONS_MSVC})
        set(PROJECT_OPTIONS_DEBUG ${PROJECT_OPTIONS_MSVC_DEBUG})
        set(PROJECT_OPTIONS_RELEASE ${PROJECT_OPTIONS_MSVC_RELEASE})
    elseif(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
        set(PROJECT_OPTIONS ${PROJECT_OPTIONS_CLANG})
        set(PROJECT_OPTIONS_DEBUG ${PROJECT_OPTIONS_CLANG_DEBUG})
        set(PROJECT_OPTIONS_RELEASE ${PROJECT_OPTIONS_CLANG_RELEASE})
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        # gcov: According to https://gcovr.com/en/stable/cookbook.html#out-of-source-builds-with-cmake
        # These flags are used if cmake is called with -DCMAKE_BUILD_TYPE=PROFILE
        set(CMAKE_C_FLAGS_PROFILE --coverage)
        set(CMAKE_CXX_FLAGS_PROFILE --coverage)

        set(PROJECT_OPTIONS ${PROJECT_OPTIONS_GCC})
        set(PROJECT_OPTIONS_DEBUG ${PROJECT_OPTIONS_GCC_DEBUG})
        set(PROJECT_OPTIONS_RELEASE ${PROJECT_OPTIONS_GCC_RELEASE})
    else()
    message(AUTHOR_WARNING "No compiler options set for CXX compiler: '${CMAKE_CXX_COMPILER_ID}'")

    # TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/71
    # support Intel compiler
    endif()

    message(DEBUG "${project_name} debug: ${PROJECT_OPTIONS_DEBUG}")
    message(DEBUG "${project_name} release: ${PROJECT_OPTIONS_RELEASE}")
    message(DEBUG "${project_name} ${PROJECT_OPTIONS}")

    target_compile_options(
        ${project_name}
        INTERFACE
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:DEBUG>>:${PROJECT_OPTIONS_DEBUG}>
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELEASE>>:${PROJECT_OPTIONS_RELEASE}>
        $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_OPTIONS}>
    )

    target_include_directories(${project_name} INTERFACE ${CMAKE_BINARY_DIR})
    target_include_directories(${project_name} INTERFACE ${VF_SRC_DIR})
    target_include_directories(${project_name} INTERFACE ${VF_SRC_DIR}/gpu)
    target_include_directories(${project_name} INTERFACE ${VF_SRC_DIR}/cpu)
    target_include_directories(${project_name} INTERFACE ${VF_SRC_DIR}/basics)
    target_include_directories(${project_name} INTERFACE ${VF_SRC_DIR}/gpu/core)
    target_include_directories(${project_name} INTERFACE ${VF_SRC_DIR}/gpu/GridGenerator)
    target_include_directories(${project_name} INTERFACE ${VF_SRC_DIR}/cpu/core)

    target_compile_definitions(${project_name} INTERFACE ${CMAKE_CXX_COMPILER_ID})

    #################################################################
    ###   OS DEFINES                                              ###
    #################################################################
    IF(WIN32)
        target_compile_definitions(${project_name} INTERFACE __WIN__)
    ELSEIF(UNIX)
        target_compile_definitions(${project_name} INTERFACE __unix__)
        IF(APPLE)
            target_compile_definitions(${project_name} INTERFACE __APPLE__)
        endif()
    ENDIF()

    if(VF_ENABLE_DOUBLE_ACCURACY)
        target_compile_definitions(${project_name} INTERFACE VF_DOUBLE_ACCURACY)
        message(STATUS "Configure VirtualFluids with double precision")
    else()
        message(STATUS "Configure VirtualFluids with single precision")
    endif()


endfunction()
