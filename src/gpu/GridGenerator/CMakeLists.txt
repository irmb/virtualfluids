vf_add_library(NAME GridGenerator PUBLIC_LINK basics logger OpenMP::OpenMP_CXX PRIVATE_LINK project_options project_warnings)
vf_add_tests()

target_include_directories(GridGenerator PRIVATE "${VF_ROOT_DIR}/src/gpu/GridGenerator/")


if(VF_ENABLE_UNIT_TESTS)
   target_include_directories(GridGeneratorTests PRIVATE "${VF_ROOT_DIR}/src/basics/")
endif()


## WARNINGS
# TODO: Remove this and fix those warnings: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/163
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    # Using CMakes De-Duplication of Flags "Hack" (https://cmake.org/cmake/help/latest/prop_tgt/COMPILE_OPTIONS.html#option-de-duplication)
    target_compile_options(GridGenerator PUBLIC "-Wold-style-cast;-Wno-old-style-cast")
    target_compile_options(GridGenerator PUBLIC "-Wsign-conversion;-Wno-sign-conversion")
    target_compile_options(GridGenerator PUBLIC "-Wdouble-promotion;-Wno-double-promotion")
    target_compile_options(GridGenerator PUBLIC "-Wshadow;-Wno-shadow")
    target_compile_options(GridGenerator PRIVATE "-Wcast-align;-Wno-cast-align")

    target_compile_options(GridGenerator PRIVATE "-Wextra") # Wno-unused-parameter comes with -Wextra
    target_compile_options(GridGenerator PRIVATE "-Wall") # Wno-unused-function comes with -Wall
    target_compile_options(GridGenerator PRIVATE "-Wunused") # -Wno-unused-function comes with -Wunused
    target_compile_options(GridGenerator PRIVATE "-Wno-unused-parameter")
    target_compile_options(GridGenerator PRIVATE "-Wno-unused-function")
    target_compile_options(GridGenerator PRIVATE "-Wno-implicit-float-conversion")
    target_compile_options(GridGenerator PRIVATE "-Wno-comment")
    target_compile_options(GridGenerator PRIVATE "-Wno-reorder-ctor")
    target_compile_options(GridGenerator PRIVATE "-Wno-uninitialized")
    target_compile_options(GridGenerator PRIVATE "-Wno-reorder")
    target_compile_options(GridGenerator PRIVATE "-Wno-strict-aliasing")
    target_compile_options(GridGenerator PRIVATE "-Wformat=2;-Wno-format-nonliteral")
    target_compile_options(GridGenerator PRIVATE "-Wno-null-dereference")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(GridGenerator PUBLIC "-Wuseless-cast;-Wno-useless-cast")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(GridGenerator PRIVATE "/wd4244") # conversion
    target_compile_options(GridGenerator PRIVATE "/wd4100") # unreferenced formal parameter
    target_compile_options(GridGenerator PRIVATE "/wd4701") # potentially uninitialized local variable
    target_compile_options(GridGenerator PRIVATE "/wd4703") # potentially uninitialized local pointer
    target_compile_options(GridGenerator PRIVATE "/wd4458") # declaration of 'XXX' hides class member
    target_compile_options(GridGenerator PRIVATE "/wd4305") #'initializing': truncation from double to float
    target_compile_options(GridGenerator PRIVATE "/wd4389") #''==': signed/unsigned mismatch
endif()