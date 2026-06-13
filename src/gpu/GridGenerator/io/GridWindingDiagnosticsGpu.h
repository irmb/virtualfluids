//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_io io
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Hussein Alihussein
//! \brief Declares GPU grid-winding diagnostics helpers for missing links and Q lines.
//! \note Generated with the assistance of OpenAI Codex (GPT-5). Reviewed and adapted by author.
//=======================================================================================
#ifndef GRID_WINDING_DIAGNOSTICS_H_
#define GRID_WINDING_DIAGNOSTICS_H_
#include <functional>
#include <string>
#include <vector>

#include "global.h"

namespace vf::grid_winding
{
struct MissingLinkCollection;
struct QLineCollection;
struct SubgridDistanceStats;
} // namespace vf::grid_winding

namespace vf::parallel
{
class Communicator;
}

class GbTriFaceMesh3D;

namespace vf::gpu {
class Grid;
//! \brief GPU-side wrapper that forwards grid-winding diagnostics
//! collections to the common, platform-agnostic writers and provides
//! helpers to collect diagnostics from GPU grids.
class GridWindingWriter
{
public:
    //! \brief Write missing boundary links using the common diagnostics IO.
    static void writeMissingLinks(const vf::grid_winding::MissingLinkCollection& collection,
                                  const std::string&                              basePath);

    //! \brief Write Q-lines using the common diagnostics IO.
    static void writeQLines(const vf::grid_winding::QLineCollection& collection,
                            const std::string&                        basePath);

};

namespace grid_winding
{

//! \brief Collect Q-lines from GPU grids into a QLineCollection.
void collectQLines(const std::vector<SPtr<Grid>>&          grids,
                   const SPtr<vf::parallel::Communicator>& comm,
                   vf::grid_winding::QLineCollection&      collection);

//! \brief Convenience wrapper for grid-winding diagnostics with STL surfaces.
void runBoundaryDiagnosticsForSurfaces(const std::vector<SPtr<Grid>>&          grids,
                                       const SPtr<GbTriFaceMesh3D>&            buildingSurface,
                                       const SPtr<GbTriFaceMesh3D>&            bodenSurface,
                                       bool                                    writeMissingLinks,
                                       bool                                    writeQLinesFlag,
                                       const std::string&                      path,
                                       const SPtr<vf::parallel::Communicator>& comm,
                                       bool                                    isRoot);

#if defined(VF_HAS_FAST_WINDING)
void runBoundaryDiagnostics(const std::vector<SPtr<Grid>>&                          grids,
                            const std::vector<GbTriFaceMesh3D*>&                    surfaces,
                            std::vector<GbTriFaceMesh3D*>&                          kdMeshes,
                            bool                                                    writeMissingLinks,
                            bool                                                    writeQLinesFlag,
                            const std::string&                                      path,
                            const SPtr<vf::parallel::Communicator>&                 comm,
                            bool                                                    isRoot);

namespace visualization
{
struct NodeFlags
{
    double isFluid{ 0.0 };
    double isBoundary{ 0.0 };
};

using NodeFlagProvider = std::function<NodeFlags(char type)>;

void exportGridVisualization(const std::vector<SPtr<Grid>>& grids,
                             const SPtr<GbTriFaceMesh3D>&   surface,
                             const std::string&            outputPath,
                             NodeFlagProvider              flagProvider = NodeFlagProvider{});
} // namespace visualization
#endif
} // namespace grid_winding
} // namespace vf::gpu
#endif
//! \}
