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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file ActuatorFarm.cu
//! \ingroup PreCollisionInteractor
//! \author Henrik Asmuth, Henry Korb, Anna Wellmann
//======================================================================================

#ifndef ACTUATOR_FARM_INLINES
#define ACTUATOR_FARM_INLINES

#include "basics/DataTypes.h"

struct TurbineNodeIndex {
    uint turbine;
    uint blade;
    uint bladeNode;
};

__host__ __device__ __inline__ uint calcNodeIndexInBladeArrays(uint bladeNode, uint numberOfNodesPerBlade, uint blade, uint numberOfBlades, uint turbine)
{
    // see https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/merge_requests/248 for visualization
    return bladeNode + numberOfNodesPerBlade * (blade + numberOfBlades * turbine);
}

__host__ __device__ __inline__ uint calcNodeIndexInBladeArrays(const TurbineNodeIndex &turbineNodeIndex, uint numberOfNodesPerBlade, uint numberOfBlades)
{
    return calcNodeIndexInBladeArrays(turbineNodeIndex.bladeNode, numberOfNodesPerBlade, turbineNodeIndex.blade, numberOfBlades, turbineNodeIndex.turbine);
}

__host__ __device__ __inline__ void calcTurbineBladeAndBladeNode(uint node, uint &bladeNode, uint numberOfNodesPerBlade, uint &blade, uint numberOfBlades, uint &turbine)
{
    // see https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/merge_requests/248 for visualization
    turbine = node / (numberOfNodesPerBlade * numberOfBlades);
    uint x_off = turbine * numberOfNodesPerBlade * numberOfBlades;
    blade = (node - x_off) / numberOfNodesPerBlade;
    uint y_off = numberOfNodesPerBlade * blade + x_off;
    bladeNode = node - y_off;
}
__host__ __device__ __inline__ TurbineNodeIndex calcTurbineBladeAndBladeNode(uint node, uint numberOfNodesPerBlade, uint numberOfBlades)
{
    uint turbine;
    uint blade;
    uint bladeNode;
    calcTurbineBladeAndBladeNode(node, bladeNode, numberOfNodesPerBlade, blade, numberOfBlades, turbine);
    return { /*.turbine = */ turbine, /*.blade = */ blade, /*.bladeNode = */ bladeNode }; // Designated initializers are a C++ 20 feature
}

#endif
