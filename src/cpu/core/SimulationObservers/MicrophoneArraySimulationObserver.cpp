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
//! \author Konstantin Kutscher
//=======================================================================================
#include "MicrophoneArraySimulationObserver.h"
#include "BCArray3D.h"
#include "BCSet.h"
#include "Block3D.h"
#include <parallel/Communicator.h>
#include "D3Q27System.h"
#include "DataSet3D.h"
#include "DistributionArray3D.h"
#include "Grid3D.h"
#include "ILBMKernel.h"
#include "UbScheduler.h"
#include "Vector3D.h"
#include <sstream>

MicrophoneArraySimulationObserver::MicrophoneArraySimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                                                       std::shared_ptr<vf::parallel::Communicator> comm)
    : SimulationObserver(grid, s), path(path), comm(comm)
{
    count = 0;
    micID = 0;
}

MicrophoneArraySimulationObserver::~MicrophoneArraySimulationObserver() = default;

void MicrophoneArraySimulationObserver::update(real step)
{
    if (microphones.size() > 0) {
        collectData(step);

        if (scheduler->isDue(step))
            writeFile(step);
    }

    UBLOG(logDEBUG3, "MicrophoneArraySimulationObserver::update:" << step);
}

bool MicrophoneArraySimulationObserver::addMicrophone(Vector3D coords)
{
    micID++;
    //   UbTupleInt3 blockIndexes = grid->getBlockIndexes(coords[0], coords[1], coords[2]);

    int minInitLevel = this->grid->getCoarsestInitializedLevel();
    int maxInitLevel = this->grid->getFinestInitializedLevel();

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        UbTupleInt3 blockIndexes = grid->getBlockIndexes(coords[0], coords[1], coords[2], level);
        SPtr<Block3D> block = grid->getBlock(val<1>(blockIndexes), val<2>(blockIndexes), val<3>(blockIndexes), level);
        if (block) {
            SPtr<ILBMKernel> kernel = block->getKernel();
            if (kernel) {
                SPtr<BCArray3D> bcarray = kernel->getBCSet()->getBCArray();
                UbTupleInt3 nodes       = grid->getNodeIndexes(block, coords[0], coords[1], coords[2]);
                if (!bcarray->isUndefined(val<1>(nodes), val<2>(nodes), val<3>(nodes))) {

                    if (kernel->getCompressible()) {
                        calcMacros = &D3Q27System::calcCompMacroscopicValues;
                    } else {
                        calcMacros = &D3Q27System::calcIncompMacroscopicValues;
                    }
                    SPtr<Mic> mic(new Mic);
                    mic->id           = micID;
                    mic->distridution = kernel->getDataSet()->getFdistributions();
                    mic->nodeIndexes  = grid->getNodeIndexes(block, coords[0], coords[1], coords[2]);
                    microphones.push_back(mic);

                    strVector.push_back(SPtr<std::stringstream>(new std::stringstream));

                    std::string fname = path + "/mic/mic_" + UbSystem::toString(micID) + ".csv";
                    std::ofstream ostr;
                    ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
                    if (!ostr) {
                        ostr.clear();
                        std::string path = UbSystem::getPathFromString(fname);
                        if (path.size() > 0) {
                            UbSystem::makeDirectory(path);
                            ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
                        }
                        if (!ostr)
                            throw UbException(UB_EXARGS, "couldn't open file " + fname);
                    }
                    ostr << "#microphone position: " << coords[0] << "; " << coords[1] << "; " << coords[2] << "; "
                         << "\n";
                    ostr.close();
                    return true;
                }
            }
        }
    }
    return false;
}

void MicrophoneArraySimulationObserver::collectData(real step)
{
    for (std::size_t i = 0; i < microphones.size(); i++) {
        real f[D3Q27System::ENDF + 1];
        microphones[i]->distridution->getPreCollisionDistribution(f, val<1>(microphones[i]->nodeIndexes),
                                                      val<2>(microphones[i]->nodeIndexes),
                                                      val<3>(microphones[i]->nodeIndexes));
        real vx1, vx2, vx3, rho;
        calcMacros(f, rho, vx1, vx2, vx3);
        *strVector[i] << step << ';' << rho << '\n';
    }
}

void MicrophoneArraySimulationObserver::writeFile(real /*step*/)
{
    for (std::size_t i = 0; i < microphones.size(); i++) {
        std::string fname = path + "/mic/mic_" + UbSystem::toString(microphones[i]->id) + ".csv";
        std::ofstream ostr;
        ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
        if (!ostr) {
            ostr.clear();
            std::string path = UbSystem::getPathFromString(fname);
            if (path.size() > 0) {
                UbSystem::makeDirectory(path);
                ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
            }
            if (!ostr)
                throw UbException(UB_EXARGS, "couldn't open file " + fname);
        }
        ostr << strVector[i]->str();
        ostr.close();
        strVector[i] = SPtr<std::stringstream>(new std::stringstream);
    }
}
//! \}
