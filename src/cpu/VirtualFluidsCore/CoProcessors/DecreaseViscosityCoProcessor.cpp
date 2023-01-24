/*
 *  DecreaseViscosityCoProcessor
 *
 *  Created on: 10.05.2013
 *  Author: uphoff
 */

#include "DecreaseViscosityCoProcessor.h"

#include <vector>

#include "Block3D.h"
#include <mpi/Communicator.h>
#include "Grid3D.h"
#include "LBMKernel.h"
#include "UbScheduler.h"

DecreaseViscosityCoProcessor::DecreaseViscosityCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, mu::Parser *nueFunc,
                                                           std::shared_ptr<vf::mpi::Communicator> comm)

    : CoProcessor(grid, s), nueFunc(nueFunc), comm(comm)
{
    if (comm->getProcessID() == comm->getRoot()) {
    }
}
//////////////////////////////////////////////////////////////////////////
DecreaseViscosityCoProcessor::~DecreaseViscosityCoProcessor() = default;
//////////////////////////////////////////////////////////////////////////
void DecreaseViscosityCoProcessor::process(double step)
{
    if (scheduler->isDue(step))
        setViscosity(step);
}
//////////////////////////////////////////////////////////////////////////
void DecreaseViscosityCoProcessor::setViscosity(double step)
{

    UBLOG(logDEBUG3, "DecreaseViscosityCoProcessor::update:" << step);
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
        double nue = nueFunc->Eval();

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
