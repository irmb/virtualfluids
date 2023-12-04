
vf_add_library(PUBLIC_LINK logger basics PRIVATE_LINK project_warnings project_options)

if(MSVC)
    target_link_libraries(parallel PRIVATE ws2_32)
endif()

# TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/139
# if(VF_ENABLE_MPI)
    target_link_libraries(parallel PUBLIC MPI::MPI_CXX)
# endif()

## WARNINGS
# TODO: Remove this and fix those warnings: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/163
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    # Using CMakes De-Duplication of Flags "Hack" (https://cmake.org/cmake/help/latest/prop_tgt/COMPILE_OPTIONS.html#option-de-duplication)
    target_compile_options(parallel PRIVATE "-Wold-style-cast;-Wno-old-style-cast")
    target_compile_options(parallel PRIVATE "-Wsign-conversion;-Wno-sign-conversion")
    target_compile_options(parallel PRIVATE "-Wconversion;-Wno-conversion")

    target_compile_options(parallel PRIVATE "-Wextra") # Wno-unused-parameter comes with -Wextra
    target_compile_options(parallel PRIVATE "-Wno-unused-parameter")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(parallel PRIVATE "/wd4100") # warning C4100: unreferenced formal paramete
endif()