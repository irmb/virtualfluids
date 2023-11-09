
set(LIGGGHTS_LIBRARIES optimized ${LIGGGHTS_RELEASE_LIBRARY} debug ${LIGGGHTS_DEBUG_LIBRARY})

vf_add_library(BUILDTYPE static PUBLIC_LINK basics muparser MPI::MPI_CXX cpu_core ${LIGGGHTS_LIBRARIES})

vf_get_library_name(library_name)

#target_link_directories(${library_name} PUBLIC ${LIGGGHTS_BINARY_DIR})
target_include_directories(${library_name} PUBLIC ${LIGGGHTS_SOURCE_DIR})

