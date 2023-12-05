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
//! \author Henrik Asmuth, Henry Korb
//======================================================================================
#include "ActuatorFarmInlines.h"
#include "ActuatorFarmStandalone.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <logger/Logger.h>
#include <cuda_helper/CudaGrid.h>
#include <basics/constants/NumericConstants.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>

#include "Utilities/GeometryUtils.h"
#include "LBM/GPUHelperFunctions/KernelUtilities.h"
#include "Parameter/Parameter.h"
#include "Parameter/CudaStreamManager.h"
#include "DataStructureInitializer/GridProvider.h"
#include "GPU/CudaMemoryManager.h"

using namespace vf::basics::constant;

std::vector<real> ActuatorFarmStandalone::computeBladeRadii(const real diameter, const uint numberOfNodesPerBlade)
{
    const real dr = c1o2 * diameter / numberOfNodesPerBlade;
    std::vector<real> blade_radii(numberOfNodesPerBlade);
    for (uint node = 0; node < numberOfNodesPerBlade; node++)
        blade_radii[node] = dr * (node + 0.5);
    return blade_radii;
}

void ActuatorFarmStandalone::updateForcesAndCoordinates()
{
    const real lift_coefficient = c1o1;
    const real drag_coefficient = c0o1;
    const real c0 = 20 * c1o10;
    const real delta_azimuth = c2Pi / this->numberOfBlades;

    for (uint turbine = 0; turbine < this->numberOfTurbines; turbine++) {
        const real rotor_speed = this->rotorSpeeds[turbine];
        const real azimuth_old = this->azimuths[turbine];
        const real azimuth_new = azimuth_old + deltaT * rotor_speed;
        this->azimuths[turbine] = azimuth_new > c2Pi ? azimuth_new - c2Pi : azimuth_new;

        for (uint blade = 0; blade < this->numberOfBlades; blade++) {
            const real local_azimuth_new = azimuth_new + blade * delta_azimuth;

            real last_node_radius = c0o1;
            real current_node_radius = c0o1;
            real next_node_radius = this->bladeRadii[0];

            for (uint bladeNode = 0; bladeNode < this->numberOfNodesPerBlade; bladeNode++) {
                const uint node = calcNodeIndexInBladeArrays({ turbine, blade, bladeNode }, this->numberOfNodesPerBlade,
                                                             this->numberOfBlades);

                real u_rel, v_rel, w_rel;
                rotateFromGlobalToBlade(u_rel, v_rel, w_rel,
                                        this->bladeVelocitiesXH[node],
                                        this->bladeVelocitiesYH[node],
                                        this->bladeVelocitiesZH[node],
                                        azimuth_old + delta_azimuth);

                last_node_radius = current_node_radius;
                current_node_radius = next_node_radius;
                next_node_radius =
                    bladeNode < this->numberOfNodesPerBlade - 1 ? this->bladeRadii[bladeNode + 1] : this->diameter * c1o2;

                const real dr = c1o2 * (next_node_radius - last_node_radius);

                v_rel += current_node_radius * rotor_speed;
                const real u_rel_sq = u_rel * u_rel + v_rel * v_rel;
                const real phi = atan2(u_rel, v_rel);

                const real tmp = c4o1 * current_node_radius / this->diameter - c1o1;
                const real chord = c0 * sqrt(c1o1 - tmp * tmp);
                const real normal_coefficient = lift_coefficient * cos(phi) + drag_coefficient * sin(phi);
                const real tangential_coefficient = lift_coefficient * sin(phi) - drag_coefficient * cos(phi);
                const real fx = -c1o2 * u_rel_sq * chord * this->density * normal_coefficient * dr;
                const real fy = -c1o2 * u_rel_sq * chord * this->density * tangential_coefficient * dr;

                rotateFromBladeToGlobal(fx, fy, c0o1,
                                        this->bladeForcesXH[node], this->bladeForcesYH[node], this->bladeForcesZH[node],
                                        local_azimuth_new);
                rotateFromBladeToGlobal(c0o1, c0o1, current_node_radius,
                                        this->bladeCoordsXH[node], this->bladeCoordsYH[node], this->bladeCoordsZH[node],
                                        local_azimuth_new);
                bladeCoordsXH[node] += this->turbinePosXH[turbine];
                bladeCoordsYH[node] += this->turbinePosYH[turbine];
                bladeCoordsZH[node] += this->turbinePosZH[turbine];
            }
        }
    }
}
