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
//! \author Martin Schoenherr
//=======================================================================================
#ifndef Cp_H
#define Cp_H

#include "Calculation/Calculation.h"
#include "Cuda/CudaMemoryManager.h"
#include "Parameter/Parameter.h"

void calcCp(Parameter* para, CudaMemoryManager* cudaMemoryManager, int lev);
void printCpTopIntermediateStep(Parameter* para, unsigned int t, int lev);
void printCpTop(Parameter* para, CudaMemoryManager* cudaMemoryManager, int lev);
void printCpBottom(Parameter* para, CudaMemoryManager* cudaMemoryManager);
void printCpBottom2(Parameter* para, CudaMemoryManager* cudaMemoryManager);



void excludeGridInterfaceNodesForMirror(Parameter* para, int lev);
void calcPressForMirror(Parameter* para, CudaMemoryManager* cudaMemoryManager, int lev);
//Ensight Gold
void printCaseFile(Parameter* para);
void printGeoFile(Parameter* para, bool fileFormat);
void printScalars(Parameter* para, bool fileFormat);
//functions to write binary files
void writeIntToFile(const int &i, std::ofstream &ofile);
void writeFloatToFile(const float &f, std::ofstream &ofile);
void writeStringToFile(const std::string &s, std::ofstream &ofile);

#endif
