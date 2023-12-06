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
//=======================================================================================
#ifndef CalcTurbulenceIntensity_H
#define CalcTurbulenceIntensity_H

#include <string>
#include <vector>

#include <basics/DataTypes.h>

class Parameter;
class CudaMemoryManager;

void allocTurbulenceIntensity(Parameter* para, CudaMemoryManager* cudaMemoryManager);
void calcVelocityAndFluctuations(Parameter* para, CudaMemoryManager* cudaMemoryManager, uint tdiff);
void calcTurbulenceIntensity(Parameter* para, CudaMemoryManager* cudaMemoryManager, uint tdiff);
void resetVelocityFluctuationsAndMeans(Parameter* para, CudaMemoryManager* cudaMemoryManager);
void cudaFreeTurbulenceIntensityArrays(Parameter* para, CudaMemoryManager* cudaMemoryManager);

void writeTurbulenceIntensityToFile(Parameter* para, uint timestep);
void writeVeloFluctuationToFile(Parameter* para, uint timeste);
void writeVeloMeansToFile(Parameter* para, uint timestep);
void writeAllTiDatafToFile(Parameter* para, uint timestep);

void writeTiStuffToFile(Parameter* para, uint timestep, unsigned long long sizeOfTiArray, std::vector<real*>& data,
                        std::vector<std::string>& datanames);

#endif
