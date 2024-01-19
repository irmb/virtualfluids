# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 1995-2013, Regents of the University of Minnesota

# Where the header files reside
#include_directories(../)

# Build program.
add_executable(strings strings.c)
add_executable(gksort gksort.c)
add_executable(fis fis.c)
add_executable(rw rw.c)
add_executable(gkgraph gkgraph.c)
foreach(prog strings gksort fis rw gkgraph)
  target_link_libraries(${prog} GKlib)
endforeach(prog)

