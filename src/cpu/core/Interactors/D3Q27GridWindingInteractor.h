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
//! \addtogroup cpu_Interactors Interactors
//! \ingroup cpu_core core
//! \{
//! \author Hussein Alihussein
//! \note Helper interactor for CPU grid-winding diagnostics, introduced with assistance from GPT-5.1.
//=======================================================================================

#ifndef D3Q27GRIDWINDINGINTERACTOR_H
#define D3Q27GRIDWINDINGINTERACTOR_H

#include <array>
#include <string>

#include <PointerDefinitions.h>

#include "D3Q27Interactor.h"

#include <basics/geometry3d/winding/GridWindingFastDefaults.h>
#include <basics/geometry3d/winding/GridWindingSubgridDistances.h>

class GbTriFaceMesh3D;
class Grid3D;
class BC;
class BCSet;
class Block3D;

// namespace vf::geometry
// {
// class GeneralizedWindingNumber;
// } // namespace vf::geometry

namespace vf::parallel
{
class Communicator;
} // namespace vf::parallel

//! \brief Interactor that can run CPU grid-winding based boundary diagnostics.
class D3Q27GridWindingInteractor : public D3Q27Interactor
{
public:
    using D3Q27Interactor::D3Q27Interactor;

    D3Q27GridWindingInteractor()  = default;
    ~D3Q27GridWindingInteractor() override = default;

    void initInteractor(const real &timeStep = 0) override;
    void updateInteractor(const real &timestep = 0) override;
    bool setDifferencesToGbObject3D(const SPtr<Block3D> block) override;
    void setQs(const real &timeStep);

    //! \brief Compute subgrid distances and Q-values for the current geometry/grid setup.
    //!
    //! The function is a thin wrapper around vf::grid_winding::cpu::processBoundaryData.
    void computeBoundaryData(const SPtr<GbTriFaceMesh3D>           &surface,
                             const SPtr<Grid3D>                    &grid,
                             const SPtr<BCSet>                     &bcSet,
                             const std::vector<SPtr<BC>>           &boundaryAdapters,
                             double                                 timeStep,
                             const SPtr<vf::parallel::Communicator>& comm,
                             const std::string                     &basePath,
                             const vf::grid_winding::cpu::BoundaryProcessingConfig &config);

    //! \brief Access the last boundary-processing result.
    [[nodiscard]] const vf::grid_winding::cpu::BoundaryProcessingResult &getBoundaryResult() const
    {
        return boundaryResult_;
    }

    //! \brief Convenience access to the accumulated subgrid statistics.
    [[nodiscard]] const vf::grid_winding::SubgridDistanceStats &getSubgridStats() const
    {
        return boundaryResult_.stats;
    }

private:
    void rebuildNodeIndexMaps();
    SPtr<BCSet> findBCSet(const SPtr<Grid3D> &grid) const;

    vf::grid_winding::cpu::BoundaryProcessingResult boundaryResult_{};
};

namespace vf::grid_winding::cpu
{
// void markSolidsWithWinding(const SPtr<Grid3D> &grid,
//                            const vf::geometry::GeneralizedWindingNumber &winding,
//                            const std::array<double, 3> &meshMin,
//                            const std::array<double, 3> &meshMax,
//                            double threshold = 0.5,
//                            double tolerance = 1.0e-4);

BoundaryProcessingResult processBoundaryData(
    const SPtr<GbTriFaceMesh3D>               &surface,
    const SPtr<Grid3D>                        &grid,
    const SPtr<BCSet>                         &bcSet,
    const std::vector<SPtr<BC>>               &boundaryAdapters,
    const ::D3Q27Interactor                   *adapterInteractor,
    double                                     timeStep,
    const SPtr<vf::parallel::Communicator>    &comm,
    const std::string                         &basePath,
    const BoundaryProcessingConfig            &config = {});

void computeSubgridDistancesStandalone(const SPtr<GbTriFaceMesh3D> &surface,
                                       const SPtr<Grid3D>          &grid,
                                       const SPtr<BCSet>           &bcSet,
                                       const std::vector<SPtr<BC>> &boundaryAdapters = {},
                                       const ::D3Q27Interactor     *adapterInteractor = nullptr,
                                       double                       timeStep = 0.0,
                                       vf::grid_winding::SubgridDistanceStats *stats = nullptr);

#if defined(VF_HAS_FAST_WINDING)
//! \brief Mark CPU grid nodes as solids using fast winding.
//!
//! \note Defaults are repeated here because this wrapper is a separate public function.
//! C++ does not automatically reuse defaults from another function it calls.
void markSolidsWithFastWinding(const SPtr<GbTriFaceMesh3D> &surface,
                               const SPtr<Grid3D>          &grid,
                               float accuracyScale = vf::grid_winding::FastWindingDefaultAccuracyScale,
                               float threshold = vf::grid_winding::FastWindingDefaultThreshold,
                               float tolerance = vf::grid_winding::FastWindingDefaultTolerance);
#endif
} // namespace vf::grid_winding::cpu

inline void computeSubgridDistancesStandalone(const SPtr<GbTriFaceMesh3D> &surface,
                                              const SPtr<Grid3D>          &grid,
                                              const SPtr<BCSet>           &bcSet,
                                              const std::vector<SPtr<BC>> &boundaryAdapters = {},
                                              const ::D3Q27Interactor     *adapterInteractor = nullptr,
                                              double                       timeStep = 0.0,
                                              vf::grid_winding::SubgridDistanceStats *stats = nullptr)
{
    vf::grid_winding::cpu::computeSubgridDistancesStandalone(surface, grid, bcSet, boundaryAdapters,
                                                             adapterInteractor, timeStep, stats);
}

#endif // D3Q27GRIDWINDINGINTERACTOR_H

//! \}
