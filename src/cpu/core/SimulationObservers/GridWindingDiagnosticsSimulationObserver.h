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
//! \addtogroup cpu_SimulationObservers SimulationObservers
//! \ingroup cpu_core core
//! \{
//! \author Hussein Alihussein
//! \note CPU grid-winding diagnostics observer, introduced with assistance from GPT-5.1.
//=======================================================================================

#ifndef GRIDWINDINGDIAGNOSTICSSIMULATIONOBSERVER_H
#define GRIDWINDINGDIAGNOSTICSSIMULATIONOBSERVER_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "SimulationObserver.h"

#include <cpu/core/Interactors/D3Q27GridWindingInteractor.h>

class GbTriFaceMesh3D;
class Grid3D;
class BC;
class BCSet;
class UbScheduler;

namespace vf::parallel
{
class Communicator;
} // namespace vf::parallel

class D3Q27GridWindingInteractor;

namespace vf::grid_winding::cpu
{
void finalizeMissingLinkCollection(vf::grid_winding::MissingLinkCollection &collection,
                                   const SPtr<GbTriFaceMesh3D>            &surface,
                                   const SPtr<Grid3D>                     &grid);
} // namespace vf::grid_winding::cpu

//! \brief Simulation observer that runs CPU grid-winding boundary diagnostics.
class GridWindingDiagnosticsSimulationObserver : public SimulationObserver
{
public:
    GridWindingDiagnosticsSimulationObserver(
        SPtr<Grid3D>                                   grid,
        SPtr<UbScheduler>                              scheduler,
        const std::vector<SPtr<GbTriFaceMesh3D>>      &surfaces,
        const SPtr<BCSet>                             &bcSet,
        const std::vector<SPtr<BC>>                   &boundaryAdapters,
        const SPtr<D3Q27GridWindingInteractor>        &interactor,
        const std::string                             &path,
        const SPtr<vf::parallel::Communicator>        &comm,
        const vf::grid_winding::cpu::BoundaryProcessingConfig &config = {});

    GridWindingDiagnosticsSimulationObserver(
        SPtr<Grid3D>                                   grid,
        SPtr<UbScheduler>                              scheduler,
        const SPtr<GbTriFaceMesh3D>                   &surface,
        const SPtr<BCSet>                             &bcSet,
        const std::vector<SPtr<BC>>                   &boundaryAdapters,
        const SPtr<D3Q27GridWindingInteractor>        &interactor,
        const std::string                             &path,
        const SPtr<vf::parallel::Communicator>        &comm,
        const vf::grid_winding::cpu::BoundaryProcessingConfig &config = {});

    ~GridWindingDiagnosticsSimulationObserver() override = default;

    void update(real step) override;

    [[nodiscard]] const vf::grid_winding::SubgridDistanceStats &getSubgridStats() const
    {
        return stats_;
    }

private:

    void ensureOutputDirectory() const;
    void writeMissingLinkReport() const;
    void writeQReports() const;
    void writeWindingVolumeReport(real step) const;

    std::string path_;
    SPtr<vf::parallel::Communicator> comm_;
    SPtr<GbTriFaceMesh3D> surface_;
    SPtr<BCSet> bcSet_;
    std::vector<SPtr<BC>> boundaryAdapters_;
    SPtr<D3Q27GridWindingInteractor> interactor_;
    vf::grid_winding::cpu::BoundaryProcessingConfig config_;

    bool diagnosticsComputed_{ false };

    vf::grid_winding::SubgridDistanceStats  stats_;
    vf::grid_winding::MissingLinkCollection missingLinks_;
    vf::grid_winding::QLineCollection       qLines_;
};

#endif // GRIDWINDINGDIAGNOSTICSSIMULATIONOBSERVER_H

//! \}
