cmake_minimum_required(VERSION 3.16...3.28)
project(fastwinding LANGUAGES CXX)

# -------------------------------
# libigl configuration
# -------------------------------
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"           OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK              "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"         OFF)
option(LIBIGL_WITH_LIM               "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"          OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           OFF)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          OFF)
option(LIBIGL_WITH_PNG               "Use PNG"            OFF)
option(LIBIGL_WITH_PYTHON            "Use Python"         OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       OFF)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  OFF)
option(LIBIGL_WITH_XML               "Use XML"            OFF)

# Add libigl as subdirectory (you already have libigl/ as submodule)
add_subdirectory(libigl)

# If Eigen path was provided, add it
if (EIGEN3_INCLUDE_DIR)
  include_directories(${EIGEN3_INCLUDE_DIR})
endif()

# -------------------------------
# Project sources
# -------------------------------
set(FW_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
)
file(GLOB WNFILES ${CMAKE_CURRENT_SOURCE_DIR}/WindingNumber/*.cpp)

add_executable(${PROJECT_NAME} ${FW_SOURCES} ${WNFILES})

# -------------------------------
# Build TBB
# -------------------------------
set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL " " FORCE)
set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)

add_subdirectory(tbb)
target_include_directories(tbb_static PUBLIC tbb/include)

# -------------------------------
# Linking
# -------------------------------
target_link_libraries(${PROJECT_NAME}
    PRIVATE
      igl::core
      tbb_static
)

add_executable(fastsdf
  ${CMAKE_CURRENT_SOURCE_DIR}/fastsdf.cpp
  ${WNFILES}
)

target_link_libraries(fastsdf
    PRIVATE
      igl::core
      tbb_static
)
