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
//! \addtogroup gpu_GridScaling GridScaling
//! \ingroup gpu_core core
//! \{
//! \author Anna Wellmann, Martin Schoenherr
//=======================================================================================
#ifndef GS_FACTORY
#define GS_FACTORY

#include <functional>

#include "Calculation/Calculation.h"
#include "Parameter/Parameter.h"

struct CUstream_st;

namespace vf::gpu {

struct LBMSimulationParameter;
class Parameter;

using gridScaling = std::function<void(LBMSimulationParameter *, LBMSimulationParameter *, ICells *, ICellNeigh&, CUstream_st *stream)>;
using gridScalingAdvectionDiffusion = std::function<void(LBMSimulationParameter*, LBMSimulationParameter*, ICells*, ICellNeigh&, CUstream_st* stream)>;

class GridScalingFactory
{
public:
    //! \brief An enumeration for selecting a scaling function
    enum class GridScaling {
        //! - ScaleCompressible = basic scaling for compressible fluid flow
        ScaleCompressible,
        //! - not specified scaling
        NotSpecified
    };

    //! \brief An enumeration for selecting a scaling function
    enum class GridScalingAdvectionDiffusion {
        //! - ScaleAdvectionDiffusionCompressible = basic scaling for compressible advection diffusion
        ScaleAdvectionDiffusionCompressible,
        //! - not specified scaling
        NotSpecified
    };

    void setScalingFactory(const GridScalingFactory::GridScaling gridScalingType, const GridScalingFactory::GridScalingAdvectionDiffusion gridScalingTypeAdvectionDiffusion = GridScalingAdvectionDiffusion::NotSpecified);

    [[nodiscard]] gridScaling getGridScalingFC(bool hasTurbulentViscosity) const;
    [[nodiscard]] gridScaling getGridScalingCF(bool hasTurbulentViscosity) const;

    [[nodiscard]] gridScalingAdvectionDiffusion getGridScalingAdvectionDiffusionFC(bool hasTurbulentDiffusivity) const;
    [[nodiscard]] gridScalingAdvectionDiffusion getGridScalingAdvectionDiffusionCF(bool hasTurbulentDiffusivity) const;

private:
    GridScaling gridScalingType = GridScaling::NotSpecified;
    GridScalingAdvectionDiffusion gridScalingTypeAdvectionDiffusion = GridScalingAdvectionDiffusion::NotSpecified;
};

}

#endif

//! \}
