# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 1995-2013, Regents of the University of Minnesota

# Add this directory for internal users.
#include_directories(.)
# Find sources.
file(GLOB metis_sources *.c)
# Build libmetis.
add_library(metis ${METIS_LIBRARY_TYPE} ${GKlib_sources} ${metis_sources})
if(UNIX)
  target_link_libraries(metis m)

  target_compile_options(metis PRIVATE "-Wno-format")
endif()

if(MSVC)
  target_compile_options(metis PRIVATE "/W0")
endif()

if(METIS_INSTALL)
  install(TARGETS metis
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION lib
    ARCHIVE DESTINATION lib)
endif()
