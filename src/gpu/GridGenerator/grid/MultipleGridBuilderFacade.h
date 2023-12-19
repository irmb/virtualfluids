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
//! \addtogroup gpu_grid grid
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Anna Wellmann
//! \brief A class that makes the setup of simulations on multiple gpus easier
//! \details Using this class is optional.

//=======================================================================================

#ifndef MULTIGPUGRIDHELPER_H
#define MULTIGPUGRIDHELPER_H

#include <array>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <vector>

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>
#include <basics/geometry3d/Axis.h>

#include "grid/BoundaryConditions/Side.h"
#include "utilities/communication.h"

struct GridDimensions;
class MultipleGridBuilder;
class Object;
class FileCollection;

//! \class MultipleGridBuilderFacade
//! \brief Simplifies the creation of the grids for a multi-gpu simulation
//!
//! \details Steps to set up the grids:
//!
//! - 1. initialize class with a MultipleGridBuilder and the dimensions of the entire domain
//!
//! - 2. Optional steps (their order does not matter):
//!
//!     - a. for multi gpu:
//!
//!         - addDomainSplit() [can be called multiple times]
//!
//!     - b. for grids with refinement:
//!
//!         - addFineGrid() [an be called multiple times]
//!
//!         - setNumberOfLayersForRefinement() [call once (or not at all for default)]
//!
//!     - c. for fluid domains with inserted geometry: call addGeometry()
//!
//! - 3. call createGrids()
//!
//! - 4. set boundary conditions
//!

using namespace vf::basics::constant;

class MultipleGridBuilderFacade
{
public:
    MultipleGridBuilderFacade(SPtr<MultipleGridBuilder> gridBuilder, SPtr<GridDimensions> gridDimensions,
                              std::optional<real> overlapOfSubdomains = std::nullopt);

    MultipleGridBuilderFacade(SPtr<GridDimensions> gridDimensions, std::optional<real> overlapOfSubdomains = std::nullopt);

    //! \brief split the domain in the passed direction at the passed coordinate
    //! \details multiple splits can be added to a domain
    void addDomainSplit(real coordinate, Axis direction);

    //! \brief set the overlap of the subdomains
    void setOverlapOfSubdomains(real overlap);

    //! \brief adds a fine grid
    //! \details multiple fine grids can be added
    void addFineGrid(std::shared_ptr<Object> gridShape, uint levelFine);

    //! \brief adds a geometry to the fluid domain
    //! \details add a primitive or an stl to the fluid domain
    void addGeometry(SPtr<Object> gridShape);

    //! \brief add number of layers for refinements
    //! \details calls gridBuilder->setNumberOfLayers()
    void setNumberOfLayersForRefinement(uint numberOfLayersFine, uint numberOfLayersBetweenLevels) const;

    //! \brief generate the subdomain which is specified by generatePart
    //! \details has to be called before the boundary conditions
    void createGrids(uint generatePart);

    // Boundary conditions, call after createGrids()
    void setSlipBoundaryCondition(SideType sideType, real normalX, real normalY, real normalZ) const;
    void setStressBoundaryCondition(SideType sideType, real normalX, real normalY, real normalZ, uint samplingOffset,
                                    real z0, real dx) const;
    void setVelocityBoundaryCondition(SideType sideType, real vx, real vy, real vz) const;
    void setPressureBoundaryCondition(SideType sideType, real rho) const;
    void setNoSlipBoundaryCondition(SideType sideType) const;
    void setPeriodicBoundaryCondition(bool periodic_X, bool periodic_Y, bool periodic_Z) const;
    void setPeriodicBoundaryCondition(const std::array<bool, 3>& periodicity) const;
    void setPrecursorBoundaryCondition(SideType sideType, SPtr<FileCollection> fileCollection, int timeStepsBetweenReads,
                                       real velocityX = c0o1, real velocityY = c0o1, real velocityZ = c0o1,
                                       std::vector<uint> fileLevelToGridLevelMap = {}) const;

    SPtr<MultipleGridBuilder> getGridBuilder() const;

protected:
    // index calculations
    std::array<uint, 3> convertToIndices3D(uint index1D) const;
    uint getX3D(uint index1D) const;
    uint getY3D(uint index1D) const;
    uint getZ3D(uint index1D) const;
    uint getIndex1D(uint xIndex, uint yIndex, uint zIndex) const;
    uint getIndex1D(const std::array<uint, 3>& index3D) const;

private:
    //! \brief calculate the number of subdomains in all coordinate directions
    void calculateNumberOfSubdomains();

    //! \brief sort the split locations for multiGPU cases
    void sortSplitLocations();

    //! \brief calculate and save the 3d index from the 1d index generatePart
    void calculatedIndexOfPart(uint generatePart);

    //! \brief for each direction, calculate if the current subdomain has a neighbor in this direction
    void checkForNeighbors();

    //! \brief set up coarse grids and subdomain boxes for all subdomains
    void configureSubDomainGrids();

    //! \brief set up the communication to neighboring subdomains
    void setUpCommunicationNeighbors();

    //! \brief check if all locations for domain splits are inside the grid bounds and there are no duplicates.
    void checkSplitLocations(const std::vector<real>& splits, real lowerBound, real upperBound) const;

    //! \brief add fine grids to the grid builder
    void addFineGridsToGridBuilder();

    //! \brief add geometries to the grid builder
    void addGeometriesToGridBuilder();

    //! \brief check whether a subdomain is the last one in a direction
    bool isFinalSubdomainInDirection(CommunicationDirections::CommunicationDirection direction) const;

    //! \brief get 1D index of the final subdomain in a direction, in the other directions it has the same position as the current subdomain
    uint getIndexOfFinalSubdomainInDirection(CommunicationDirections::CommunicationDirection direction) const;

    //! \brief call the grid builder's setter for a boundary condition
    template <typename function>
    void setBoundaryCondition(SideType sideType, function boundaryConditionFunction) const
    {
        if (!createGridsHasBeenCalled) {
            throw std::runtime_error(
                "MultipleGridBuilderFacade::createGrids() has to be called before setting boundary conditions");
        }

        if (sideType == SideType::GEOMETRY ||
            !hasNeighbors.at(static_cast<CommunicationDirections::CommunicationDirection>(static_cast<int>(sideType)))) {
            boundaryConditionFunction();
        }
    }

protected:
    //! \brief total number of subdomains in each coordinate direction (computed)
    std::array<uint, 3> numberOfSubdomains;

private:
    const SPtr<MultipleGridBuilder> gridBuilder;

    // basic grid dimension (set in constructor)
    const SPtr<GridDimensions> gridDimensionsDomain;

    uint numberOfSubdomainsTotal;

    //! \brief coordinates, signifying where the domain is split into subdomains (must be set in a setter)
    std::vector<real> xSplits;
    std::vector<real> ySplits;
    std::vector<real> zSplits;

    //! \brief overlap between subdomains (must be set in a setter)
    std::optional<real> overlapOfSubdomains = std::nullopt;

    //! \brief index of the current subdomains in relation to all subdomains (computed)
    std::array<uint, 3> index;

    //! \brief hasNeighbors, indicates if the current subdomains has a neighbor in a specific direction (computed)
    //! \details use the enum CommunciationDirection to access the data
    std::map<CommunicationDirections::CommunicationDirection, bool> hasNeighbors = {
        { CommunicationDirections::MX, false }, { CommunicationDirections::PX, false },
        { CommunicationDirections::MY, false }, { CommunicationDirections::PY, false },
        { CommunicationDirections::MZ, false }, { CommunicationDirections::PZ, false }
    };

    //! \brief collects fineGrids: uint is the level, Object* the gridShape
    std::vector<std::pair<std::shared_ptr<Object>, uint>> fineGrids;

    //! \brief collects geometries
    std::vector<std::shared_ptr<Object>> geometries;

    bool createGridsHasBeenCalled = false;
};

#endif
//! \}
