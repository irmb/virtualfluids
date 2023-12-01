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
//! \file DecreaseViscositySimulationObserver.cpp
//! \ingroup SimulationObservers
//! \author Sonja Uphoff
//=======================================================================================


#include "DecreaseViscositySimulationObserver.h"

#include <vector>

#include "Block3D.h"
#include <parallel/Communicator.h>
#include "Grid3D.h"
#include "LBMKernel.h"
#include "UbScheduler.h"

DecreaseViscositySimulationObserver::DecreaseViscositySimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, mu::Parser *nueFunc,
                                                           std::shared_ptr<vf::parallel::Communicator> comm)

    : SimulationObserver(grid, s), nueFunc(nueFunc), comm(comm)
{
    if (comm->getProcessID() == comm->getRoot()) {
    }
}
//////////////////////////////////////////////////////////////////////////
DecreaseViscositySimulationObserver::~DecreaseViscositySimulationObserver() = default;
//////////////////////////////////////////////////////////////////////////
void DecreaseViscositySimulationObserver::update(real step)
{
    if (scheduler->isDue(step))
        setViscosity(step);
}
//////////////////////////////////////////////////////////////////////////
void DecreaseViscositySimulationObserver::setViscosity(real step)
{

    UBLOG(logDEBUG3, "DecreaseViscositySimulationObserver::update:" << step);
    int gridRank     = grid->getRank();
    int minInitLevel = this->grid->getCoarsestInitializedLevel();
    int maxInitLevel = this->grid->getFinestInitializedLevel();

    if (comm->getProcessID() == comm->getRoot()) {

        for (int level = minInitLevel; level <= maxInitLevel; level++) {
            std::vector<SPtr<Block3D>> blockVector;
            grid->getBlocks(level, gridRank, blockVector);
            for (SPtr<Block3D> block : blockVector) {
                SPtr<ILBMKernel> kernel = block->getKernel();
            }
        }

        int istep      = static_cast<int>(step);
        this->timeStep = istep;
        nueFunc->DefineVar("t", &this->timeStep);
        real nue = nueFunc->Eval();

        for (int level = minInitLevel; level <= maxInitLevel; level++) {
            std::vector<SPtr<Block3D>> blockVector;
            grid->getBlocks(level, gridRank, blockVector);
            for (SPtr<Block3D> block : blockVector) {
                SPtr<ILBMKernel> kernel = block->getKernel();
                if (kernel) {
                    real collFactor = LBMSystem::calcCollisionFactor(nue, block->getLevel());
                    kernel->setCollisionFactor(collFactor);
                }
            }
        }
    }
}
