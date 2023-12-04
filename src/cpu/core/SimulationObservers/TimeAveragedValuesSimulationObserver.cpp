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
//! \file TimeAveragedValuesSimulationObserver.cpp
//! \ingroup SimulationObservers
//! \author Konstantin Kutscher
//=======================================================================================
#include "TimeAveragedValuesSimulationObserver.h"

#include "BCSet.h"
#include "LBMKernel.h"

#include "Block3D.h"
#include <parallel/Communicator.h>
#include "DataSet3D.h"
#include "Grid3D.h"
#include "UbScheduler.h"
#include "basics/writer/WbWriterVtkXmlASCII.h"

#include "BCArray3D.h"

TimeAveragedValuesSimulationObserver::TimeAveragedValuesSimulationObserver() = default;
//////////////////////////////////////////////////////////////////////////
TimeAveragedValuesSimulationObserver::TimeAveragedValuesSimulationObserver(SPtr<Grid3D> grid, const std::string &path,
                                                             WbWriter *const writer, SPtr<UbScheduler> s,
                                                             std::shared_ptr<vf::parallel::Communicator> comm, int options)
    : SimulationObserver(grid, s), path(path), writer(writer), comm(comm), options(options)
{
    init();
    planarAveraging = false;
    timeAveraging   = true;
}
//////////////////////////////////////////////////////////////////////////
TimeAveragedValuesSimulationObserver::TimeAveragedValuesSimulationObserver(SPtr<Grid3D> grid, const std::string &path,
                                                             WbWriter *const writer, SPtr<UbScheduler> s,
                                                             std::shared_ptr<vf::parallel::Communicator> comm, int options,
                                                             std::vector<int> levels, std::vector<real> &levelCoords,
                                                             std::vector<real> &bounds, bool timeAveraging)
    : SimulationObserver(grid, s), path(path), writer(writer), comm(comm), options(options), levels(levels),
      levelCoords(levelCoords), bounds(bounds), timeAveraging(timeAveraging)
{
    init();
    planarAveraging = true;
}
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesSimulationObserver::init()
{
    root         = comm->isRoot();
    gridRank     = grid->getRank();
    minInitLevel = this->grid->getCoarsestInitializedLevel();
    maxInitLevel = this->grid->getFinestInitializedLevel();

    withGhostLayer = false;
    iMinC          = 1;

    minStep       = scheduler->getMinBegin();
    maxStep       = scheduler->getMaxEnd();
    numberOfSteps = (maxStep - minStep);

    // function pointer
    using namespace D3Q27System;
    calcMacros = NULL;
    if (compressible) {
        calcMacros = &calcCompMacroscopicValues;
    } else {
        calcMacros = &calcIncompMacroscopicValues;
    }

    real begin        = scheduler->getMinBegin();
    real gridTimeStep = grid->getTimeStep();

    if (gridTimeStep == begin || gridTimeStep == 0) {
        initData();
    } else {
        blockVector.clear();
        blockVector.resize(maxInitLevel + 1);
        for (int level = minInitLevel; level <= maxInitLevel; level++) {
            grid->getBlocks(level, gridRank, true, blockVector[level]);
            if (blockVector[level].size() > 0)
                compressible = blockVector[level][0]->getKernel()->getCompressible();
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesSimulationObserver::initData()
{
    blockVector.clear();
    blockVector.resize(maxInitLevel + 1);

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        grid->getBlocks(level, gridRank, true, blockVector[level]);

        if (blockVector[level].size() > 0)
            compressible = blockVector[level][0]->getKernel()->getCompressible();

        for (SPtr<Block3D> block : blockVector[level]) {
            UbTupleInt3 nx = grid->getBlockNX();

            if ((options & Density) == Density) {
                SPtr<AverageValuesArray3D> ar = SPtr<AverageValuesArray3D>(
                    new AverageValuesArray3D(2, val<1>(nx) + 1, val<2>(nx) + 1, val<3>(nx) + 1, 0.0));
                block->getKernel()->getDataSet()->setAverageDensity(ar);
            }

            if ((options & Velocity) == Velocity) {
                SPtr<AverageValuesArray3D> av = SPtr<AverageValuesArray3D>(
                    new AverageValuesArray3D(3, val<1>(nx) + 1, val<2>(nx) + 1, val<3>(nx) + 1, 0.0));
                block->getKernel()->getDataSet()->setAverageVelocity(av);
            }

            if ((options & Fluctuations) == Fluctuations) {
                SPtr<AverageValuesArray3D> af = SPtr<AverageValuesArray3D>(
                    new AverageValuesArray3D(6, val<1>(nx) + 1, val<2>(nx) + 1, val<3>(nx) + 1, 0.0));
                block->getKernel()->getDataSet()->setAverageFluctuations(af);
            }

            if ((options & Triplecorrelations) == Triplecorrelations) {
                SPtr<AverageValuesArray3D> at = SPtr<AverageValuesArray3D>(
                    new AverageValuesArray3D(10, val<1>(nx) + 1, val<2>(nx) + 1, val<3>(nx) + 1, 0.0));
                block->getKernel()->getDataSet()->setAverageTriplecorrelations(at);
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesSimulationObserver::update(real step)
{
    if (step == minStep) {
        initData();
        numberOfSteps = (maxStep - minStep) + 1;
        // DEBUG/////////////////////
        // UBLOG(logINFO, "process::step = " << step << ", minStep = " << minStep << ", maxStep = " << maxStep << ",
        // numberOfSteps = " << numberOfSteps << " init()");
        ////////////////////////////
    }
    calculateSubtotal(step);

    if (step == maxStep) {
        // DEBUG/////////////////////
        // UBLOG(logINFO, "process::step = " << step << ", minStep = " << minStep << ", maxStep = " << maxStep << ",
        // numberOfSteps = " << numberOfSteps);
        ////////////////////////////

        // calculateAverageValues((double)numberOfFineSteps);
        calculateAverageValues(numberOfSteps);

        if (timeAveraging) {
            collectData(step);
        }

        if (planarAveraging) {
            planarAverage(step);
        }
    }

    UBLOG(logDEBUG3, "AverageValuesSimulationObserver::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesSimulationObserver::collectData(real step)
{
    int istep = int(step);

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (SPtr<Block3D> block : blockVector[level]) {
            if (block) {
                addData(block);
            }
        }
    }

    std::string pfilePath, partPath, subfolder, cfilePath;
    subfolder = "tav" + UbSystem::toString(istep);
    pfilePath = path + "/tav/" + subfolder;
    partPath  = pfilePath + "/tav" + UbSystem::toString(gridRank) + "_" + UbSystem::toString(istep);

    std::string partName = writer->writeOctsWithNodeData(partPath, nodes, cells, datanames, data);
    size_t found         = partName.find_last_of("/");
    std::string piece    = partName.substr(found + 1);
    piece                = subfolder + "/" + piece;

    std::vector<std::string> cellDataNames;
    std::vector<std::string> pieces = comm->gather(piece);
    if (root) {
        std::string pname =
            WbWriterVtkXmlASCII::getInstance()->writeParallelFile(pfilePath, pieces, datanames, cellDataNames);
        UBLOG(logINFO, "TimeAveragedValuesSimulationObserver::collectData() step: " << istep);
    }

    clearData();
}
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesSimulationObserver::clearData()
{
    nodes.clear();
    cells.clear();
    datanames.clear();
    data.clear();
}
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesSimulationObserver::addData(const SPtr<Block3D> block)
{
    UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);
    //   UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
    UbTupleDouble3 nodeOffset = grid->getNodeOffset(block);
    real dx                 = grid->getDeltaX(block);
    int level                 = block->getLevel();

    // Diese Daten werden geschrieben:
    datanames.resize(0);

    datanames.emplace_back("level");
    datanames.emplace_back("Rho");

    if ((options & Density) == Density) {
        datanames.emplace_back("taRho");
        datanames.emplace_back("taRhoF");
    }

    if ((options & Velocity) == Velocity) {
        datanames.emplace_back("taVx");
        datanames.emplace_back("taVy");
        datanames.emplace_back("taVz");
    }

    if ((options & Fluctuations) == Fluctuations) {
        datanames.emplace_back("taVxx");
        datanames.emplace_back("taVyy");
        datanames.emplace_back("taVzz");
        datanames.emplace_back("taVxy");
        datanames.emplace_back("taVxz");
        datanames.emplace_back("taVyz");
        datanames.emplace_back("taVyz");
    }

    if ((options & Triplecorrelations) == Triplecorrelations) {
        datanames.emplace_back("taVxxx");
        datanames.emplace_back("taVxxy");
        datanames.emplace_back("taVxxz");
        datanames.emplace_back("taVyyy");
        datanames.emplace_back("taVyyx");
        datanames.emplace_back("taVyyz");
        datanames.emplace_back("taVzzz");
        datanames.emplace_back("taVzzx");
        datanames.emplace_back("taVzzy");
        datanames.emplace_back("taVxyz");
    }

    // datanames.push_back("AvP");
    // datanames.push_back("AvPrms");

    data.resize(datanames.size());

    SPtr<ILBMKernel> kernel                 = block->getKernel();
    SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
    SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
    SPtr<AverageValuesArray3D> ar           = kernel->getDataSet()->getAverageDensity();
    SPtr<AverageValuesArray3D> av           = kernel->getDataSet()->getAverageVelocity();
    SPtr<AverageValuesArray3D> af           = kernel->getDataSet()->getAverageFluctuations();
    SPtr<AverageValuesArray3D> at           = kernel->getDataSet()->getAverageTriplecorrelations();
    // int ghostLayerWidth = kernel->getGhostLayerWidth();

    int minX1 = iMinC;
    int minX2 = iMinC;
    int minX3 = iMinC;

    int maxX1 = int(distributions->getNX1());
    int maxX2 = int(distributions->getNX2());
    int maxX3 = int(distributions->getNX3());

    // nummern vergeben und node vector erstellen + daten sammeln
    CbArray3D<int> nodeNumbers((int)maxX1, (int)maxX2, (int)maxX3, -1);

    maxX1 -= 2;
    maxX2 -= 2;
    maxX3 -= 2;

    real f[D3Q27System::ENDF + 1];
    real vx1, vx2, vx3, rho;

    // D3Q27BoundaryConditionPtr bcPtr;

    int nr = (int)nodes.size();

    for (int ix3 = minX3; ix3 <= maxX3; ix3++) {
        for (int ix2 = minX2; ix2 <= maxX2; ix2++) {
            for (int ix1 = minX1; ix1 <= maxX1; ix1++) {
                if (!bcArray->isUndefined(ix1, ix2, ix3) && !bcArray->isSolid(ix1, ix2, ix3)) {
                    int index                  = 0;
                    nodeNumbers(ix1, ix2, ix3) = nr++;
                    nodes.push_back(makeUbTuple(float(val<1>(org) - val<1>(nodeOffset) + ix1 * dx),
                                                float(val<2>(org) - val<2>(nodeOffset) + ix2 * dx),
                                                float(val<3>(org) - val<3>(nodeOffset) + ix3 * dx)));

                    data[index++].push_back(level);

                    distributions->getPreCollisionDistribution(f, ix1, ix2, ix3);
                    calcMacros(f, rho, vx1, vx2, vx3);

                    data[index++].push_back(rho);

                    if ((options & Density) == Density) {
                        data[index++].push_back((*ar)(Rho, ix1, ix2, ix3));
                        data[index++].push_back((*ar)(RhoF, ix1, ix2, ix3));
                    }

                    if ((options & Velocity) == Velocity) {
                        data[index++].push_back((*av)(Vx, ix1, ix2, ix3));
                        data[index++].push_back((*av)(Vy, ix1, ix2, ix3));
                        data[index++].push_back((*av)(Vz, ix1, ix2, ix3));
                    }

                    if ((options & Fluctuations) == Fluctuations) {
                        data[index++].push_back((*af)(Vxx, ix1, ix2, ix3));
                        data[index++].push_back((*af)(Vyy, ix1, ix2, ix3));
                        data[index++].push_back((*af)(Vzz, ix1, ix2, ix3));
                        data[index++].push_back((*af)(Vxy, ix1, ix2, ix3));
                        data[index++].push_back((*af)(Vxz, ix1, ix2, ix3));
                        data[index++].push_back((*af)(Vyz, ix1, ix2, ix3));
                    }

                    if ((options & Triplecorrelations) == Triplecorrelations) {
                        data[index++].push_back((*at)(Vxxx, ix1, ix2, ix3));
                        data[index++].push_back((*at)(Vxxy, ix1, ix2, ix3));
                        data[index++].push_back((*at)(Vxxz, ix1, ix2, ix3));
                        data[index++].push_back((*at)(Vyyy, ix1, ix2, ix3));
                        data[index++].push_back((*at)(Vyyx, ix1, ix2, ix3));
                        data[index++].push_back((*at)(Vyyz, ix1, ix2, ix3));
                        data[index++].push_back((*at)(Vzzz, ix1, ix2, ix3));
                        data[index++].push_back((*at)(Vzzx, ix1, ix2, ix3));
                        data[index++].push_back((*at)(Vzzy, ix1, ix2, ix3));
                        data[index++].push_back((*at)(Vxyz, ix1, ix2, ix3));
                    }
                }
            }
        }
    }

    maxX1 -= 1;
    maxX2 -= 1;
    maxX3 -= 1;

    int SWB, SEB, NEB, NWB, SWT, SET, NET, NWT;

    // cell vector erstellen
    for (int ix3 = minX3; ix3 <= maxX3; ix3++) {
        for (int ix2 = minX2; ix2 <= maxX2; ix2++) {
            for (int ix1 = minX1; ix1 <= maxX1; ix1++) {
                if ((SWB = nodeNumbers(ix1, ix2, ix3)) >= 0 && (SEB = nodeNumbers(ix1 + 1, ix2, ix3)) >= 0 &&
                    (NEB = nodeNumbers(ix1 + 1, ix2 + 1, ix3)) >= 0 && (NWB = nodeNumbers(ix1, ix2 + 1, ix3)) >= 0 &&
                    (SWT = nodeNumbers(ix1, ix2, ix3 + 1)) >= 0 && (SET = nodeNumbers(ix1 + 1, ix2, ix3 + 1)) >= 0 &&
                    (NET = nodeNumbers(ix1 + 1, ix2 + 1, ix3 + 1)) >= 0 &&
                    (NWT = nodeNumbers(ix1, ix2 + 1, ix3 + 1)) >= 0) {
                    cells.push_back(makeUbTuple((unsigned int)SWB, (unsigned int)SEB, (unsigned int)NEB,
                                                (unsigned int)NWB, (unsigned int)SWT, (unsigned int)SET,
                                                (unsigned int)NET, (unsigned int)NWT));
                }
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesSimulationObserver::calculateAverageValues(real timeSteps)
{
    using namespace vf::basics::constant;

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        int i;
        const int block_size = (int) blockVector[level].size();
        //#ifdef _OPENMP
        //   #pragma omp parallel for
        //#endif
        // for(SPtr<Block3D> block : blockVector[level])
        for (i = 0; i < block_size; i++) {
            SPtr<Block3D> block = blockVector[level][i];
            if (block) {
                SPtr<ILBMKernel> kernel                 = block->getKernel();
                SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
                SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
                SPtr<AverageValuesArray3D> ar           = kernel->getDataSet()->getAverageDensity();
                SPtr<AverageValuesArray3D> av           = kernel->getDataSet()->getAverageVelocity();
                SPtr<AverageValuesArray3D> af           = kernel->getDataSet()->getAverageFluctuations();
                SPtr<AverageValuesArray3D> at           = kernel->getDataSet()->getAverageTriplecorrelations();

                int minX1 = iMinC;
                int minX2 = iMinC;
                int minX3 = iMinC;

                int maxX1 = int(distributions->getNX1());
                int maxX2 = int(distributions->getNX2());
                int maxX3 = int(distributions->getNX3());

                maxX1 -= 2;
                maxX2 -= 2;
                maxX3 -= 2;

                real rho {0.}, ux {0.}, uy {0.}, uz {0.}, uxx {0.}, uzz {0.}, uyy {0.}, uxy {0.}, uxz {0.}, uyz {0.}, rhof {0.};

                for (int ix3 = minX3; ix3 <= maxX3; ix3++) {
                    for (int ix2 = minX2; ix2 <= maxX2; ix2++) {
                        for (int ix1 = minX1; ix1 <= maxX1; ix1++) {
                            if (!bcArray->isUndefined(ix1, ix2, ix3) && !bcArray->isSolid(ix1, ix2, ix3)) {
                                //////////////////////////////////////////////////////////////////////////
                                // compute average values
                                //////////////////////////////////////////////////////////////////////////

                                // mean density
                                if ((options & Density) == Density) {
                                    rho                        = (*ar)(Rho, ix1, ix2, ix3) / timeSteps;
                                    rhof                       = (*ar)(RhoF, ix1, ix2, ix3) / timeSteps;
                                    (*ar)(Rho, ix1, ix2, ix3)  = rho;
                                    (*ar)(RhoF, ix1, ix2, ix3) = rhof - rho * rho;
                                }

                                // mean velocity
                                if ((options & Velocity) == Velocity) {
                                    ux = (*av)(Vx, ix1, ix2, ix3) / timeSteps;
                                    uy = (*av)(Vy, ix1, ix2, ix3) / timeSteps;
                                    uz = (*av)(Vz, ix1, ix2, ix3) / timeSteps;

                                    (*av)(Vx, ix1, ix2, ix3) = ux;
                                    (*av)(Vy, ix1, ix2, ix3) = uy;
                                    (*av)(Vz, ix1, ix2, ix3) = uz;
                                }

                                // mean fluctuations
                                if ((options & Fluctuations) == Fluctuations) {
                                    uxx = (*af)(Vxx, ix1, ix2, ix3) / timeSteps;
                                    uyy = (*af)(Vyy, ix1, ix2, ix3) / timeSteps;
                                    uzz = (*af)(Vzz, ix1, ix2, ix3) / timeSteps;
                                    uxy = (*af)(Vxy, ix1, ix2, ix3) / timeSteps;
                                    uxz = (*af)(Vxz, ix1, ix2, ix3) / timeSteps;
                                    uyz = (*af)(Vyz, ix1, ix2, ix3) / timeSteps;

                                    (*af)(Vxx, ix1, ix2, ix3) = uxx - ux * ux;
                                    (*af)(Vyy, ix1, ix2, ix3) = uyy - uy * uy;
                                    (*af)(Vzz, ix1, ix2, ix3) = uzz - uz * uz;
                                    (*af)(Vxy, ix1, ix2, ix3) = uxy - ux * uy;
                                    (*af)(Vxz, ix1, ix2, ix3) = uxz - ux * uz;
                                    (*af)(Vyz, ix1, ix2, ix3) = uyz - uy * uz;
                                }

                                if ((options & Triplecorrelations) == Triplecorrelations) {
                                    // mean triple-correlations
                                    (*at)(Vxxx, ix1, ix2, ix3) =
                                        (*at)(Vxxx, ix1, ix2, ix3) / timeSteps - c3o1 * uxx * ux + c2o1 * ux * ux * ux;
                                    (*at)(Vxxy, ix1, ix2, ix3) = (*at)(Vxxy, ix1, ix2, ix3) / timeSteps -
                                                                 c2o1 * uxy * ux - uxx * uy + c2o1 * ux * ux * uy;
                                    (*at)(Vxxz, ix1, ix2, ix3) = (*at)(Vxxz, ix1, ix2, ix3) / timeSteps -
                                                                 c2o1 * uxz * ux - uxx * uz + c2o1 * ux * ux * uz;
                                    (*at)(Vyyy, ix1, ix2, ix3) =
                                        (*at)(Vyyy, ix1, ix2, ix3) / timeSteps - c3o1 * uyy * uy + c2o1 * uy * uy * uy;
                                    (*at)(Vyyx, ix1, ix2, ix3) = (*at)(Vyyx, ix1, ix2, ix3) / timeSteps -
                                                                 c2o1 * uxy * uy - uyy * ux + c2o1 * uy * uy * ux;
                                    (*at)(Vyyz, ix1, ix2, ix3) = (*at)(Vyyz, ix1, ix2, ix3) / timeSteps -
                                                                 c2o1 * uyz * uy - uyy * uz + c2o1 * uy * uy * uz;
                                    (*at)(Vzzz, ix1, ix2, ix3) =
                                        (*at)(Vzzz, ix1, ix2, ix3) / timeSteps - c3o1 * uzz * uz + c2o1 * uz * uz * uz;
                                    (*at)(Vzzx, ix1, ix2, ix3) = (*at)(Vzzx, ix1, ix2, ix3) / timeSteps -
                                                                 c2o1 * uxz * uz - uzz * ux + c2o1 * uz * uz * ux;
                                    (*at)(Vzzy, ix1, ix2, ix3) = (*at)(Vzzy, ix1, ix2, ix3) / timeSteps -
                                                                 c2o1 * uyz * uz - uzz * uy + c2o1 * uz * uz * uy;
                                    (*at)(Vxyz, ix1, ix2, ix3) = (*at)(Vxyz, ix1, ix2, ix3) / timeSteps - uxy * uz -
                                                                 uxz * uy - uyz * ux + c2o1 * ux * uy * uz;
                                }
                                //////////////////////////////////////////////////////////////////////////
                            }
                        }
                    }
                }
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesSimulationObserver::calculateSubtotal(real step)
{
    if (scheduler->isDue(step)) {

        // DEBUG/////////////////////
        // UBLOG(logINFO, "calculateSubtotal::step = " << step);
        ////////////////////////////
        real f[27];

        //#ifdef _OPENMP
        //#pragma omp parallel private (f)
        //#endif
        {
            for (int level = minInitLevel; level <= maxInitLevel; level++) {
                int i;
                const int block_size = (int) blockVector[level].size();
                //#ifdef _OPENMP
                //#pragma omp for schedule(dynamic)
                //#endif
                // for(SPtr<Block3D> block : blockVector[level])
                for (i = 0; i < block_size; i++) {
                    SPtr<Block3D> block = blockVector[level][i];
                    if (block) {
                        SPtr<ILBMKernel> kernel                 = block->getKernel();
                        SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
                        SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
                        SPtr<AverageValuesArray3D> ar           = kernel->getDataSet()->getAverageDensity();
                        SPtr<AverageValuesArray3D> av           = kernel->getDataSet()->getAverageVelocity();
                        SPtr<AverageValuesArray3D> af           = kernel->getDataSet()->getAverageFluctuations();
                        SPtr<AverageValuesArray3D> at           = kernel->getDataSet()->getAverageTriplecorrelations();

                        int minX1 = iMinC;
                        int minX2 = iMinC;
                        int minX3 = iMinC;

                        int maxX1 = int(distributions->getNX1());
                        int maxX2 = int(distributions->getNX2());
                        int maxX3 = int(distributions->getNX3());

                        maxX1 -= 2;
                        maxX2 -= 2;
                        maxX3 -= 2;

                        for (int ix3 = minX3; ix3 <= maxX3; ix3++) {
                            for (int ix2 = minX2; ix2 <= maxX2; ix2++) {
                                for (int ix1 = minX1; ix1 <= maxX1; ix1++) {
                                    if (!bcArray->isUndefined(ix1, ix2, ix3) && !bcArray->isSolid(ix1, ix2, ix3)) {
                                        //////////////////////////////////////////////////////////////////////////
                                        // read distribution
                                        ////////////////////////////////////////////////////////////////////////////

                                        distributions->getPreCollisionDistribution(f, ix1, ix2, ix3);
                                        //////////////////////////////////////////////////////////////////////////
                                        // compute velocity
                                        //////////////////////////////////////////////////////////////////////////
                                        real vx, vy, vz, rho;
                                        calcMacros(f, rho, vx, vy, vz);
                                        // double press = D3Q27System::calcPress(f, rho, vx, vy, vz);

                                        //////////////////////////////////////////////////////////////////////////
                                        // compute subtotals
                                        //////////////////////////////////////////////////////////////////////////

                                        // mean density
                                        if ((options & Density) == Density) {
                                            (*ar)(0, ix1, ix2, ix3)    = (*ar)(Rho, ix1, ix2, ix3) + rho;
                                            (*ar)(RhoF, ix1, ix2, ix3) = (*ar)(RhoF, ix1, ix2, ix3) + rho * rho;
                                        }

                                        // mean velocity
                                        if ((options & Velocity) == Velocity) {
                                            (*av)(Vx, ix1, ix2, ix3) = (*av)(Vx, ix1, ix2, ix3) + vx;
                                            (*av)(Vy, ix1, ix2, ix3) = (*av)(Vy, ix1, ix2, ix3) + vy;
                                            (*av)(Vz, ix1, ix2, ix3) = (*av)(Vz, ix1, ix2, ix3) + vz;
                                        }

                                        // fluctuations
                                        if ((options & Fluctuations) == Fluctuations) {
                                            (*af)(Vxx, ix1, ix2, ix3) = (*af)(Vxx, ix1, ix2, ix3) + vx * vx;
                                            (*af)(Vyy, ix1, ix2, ix3) = (*af)(Vyy, ix1, ix2, ix3) + vy * vy;
                                            (*af)(Vzz, ix1, ix2, ix3) = (*af)(Vzz, ix1, ix2, ix3) + vz * vz;
                                            (*af)(Vxy, ix1, ix2, ix3) = (*af)(Vxy, ix1, ix2, ix3) + vx * vy;
                                            (*af)(Vxz, ix1, ix2, ix3) = (*af)(Vxz, ix1, ix2, ix3) + vx * vz;
                                            (*af)(Vyz, ix1, ix2, ix3) = (*af)(Vyz, ix1, ix2, ix3) + vy * vz;
                                        }

                                        // triple-correlations
                                        if ((options & Triplecorrelations) == Triplecorrelations) {
                                            (*at)(Vxxx, ix1, ix2, ix3) = (*at)(Vxxx, ix1, ix2, ix3) + vx * vx * vx;
                                            (*at)(Vxxy, ix1, ix2, ix3) = (*at)(Vxxy, ix1, ix2, ix3) + vx * vx * vy;
                                            (*at)(Vxxz, ix1, ix2, ix3) = (*at)(Vxxz, ix1, ix2, ix3) + vx * vx * vz;
                                            (*at)(Vyyy, ix1, ix2, ix3) = (*at)(Vyyy, ix1, ix2, ix3) + vy * vy * vy;
                                            (*at)(Vyyx, ix1, ix2, ix3) = (*at)(Vyyx, ix1, ix2, ix3) + vy * vy * vx;
                                            (*at)(Vyyz, ix1, ix2, ix3) = (*at)(Vyyz, ix1, ix2, ix3) + vy * vy * vz;
                                            (*at)(Vzzz, ix1, ix2, ix3) = (*at)(Vzzz, ix1, ix2, ix3) + vz * vz * vz;
                                            (*at)(Vzzx, ix1, ix2, ix3) = (*at)(Vzzx, ix1, ix2, ix3) + vz * vz * vx;
                                            (*at)(Vzzy, ix1, ix2, ix3) = (*at)(Vzzy, ix1, ix2, ix3) + vz * vz * vy;
                                            (*at)(Vxyz, ix1, ix2, ix3) = (*at)(Vxyz, ix1, ix2, ix3) + vx * vy * vz;
                                        }
                                        //////////////////////////////////////////////////////////////////////////
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesSimulationObserver::planarAverage(real step)
{
    std::ofstream ostr;
    using namespace vf::basics::constant;

    if (root) {
        int istep         = int(step);
        std::string fname = path + "/tav/" + "tav" + UbSystem::toString(istep) + ".csv";

        ostr.open(fname.c_str(), std::ios_base::out);
        if (!ostr) {
            ostr.clear();
            std::string path = UbSystem::getPathFromString(fname);
            if (path.size() > 0) {
                UbSystem::makeDirectory(path);
                ostr.open(fname.c_str(), std::ios_base::out);
            }
            if (!ostr)
                throw UbException(UB_EXARGS, "couldn't open file " + fname);
        }

        ostr << "z";

        if ((options & Density) == Density) {
            ostr << ";Rho;RhoF";
        }
        // mean velocity
        if ((options & Velocity) == Velocity) {
            ostr << ";Vx;Vy;Vz";
        }
        // fluctuations
        if ((options & Fluctuations) == Fluctuations) {
            ostr << ";Vxx;Vyy;Vzz;Vxy;Vxz;Vyz";
        }
        // triple-correlations
        if ((options & Triplecorrelations) == Triplecorrelations) {
            ostr << ";Vxxx;Vxxy;Vxxz;Vyyy;Vyyx;Vyyz;Vzzz;Vzzx;Vzzy;Vxyz";
        }
        ostr << "\n";
    }

    int size              = (int)levels.size();
    int sizeOfLevelCoords = (int)levelCoords.size();

    if (2 * size != sizeOfLevelCoords) {
        UB_THROW(UbException(UB_EXARGS, "Number of levels coordinates don't match number of levels!"));
    }

    int k = 0;

    for (int i = 0; i < size; i++) {
        int level    = levels[i];
        real dx    = grid->getDeltaX(level);
        real start = levelCoords[k];
        real stop  = levelCoords[k + 1];

        for (real j = start; j < stop; j += dx) {
            IntegrateValuesHelper intValHelp(grid, comm, bounds[0], bounds[1], j, bounds[3], bounds[4], j + dx, level);

            std::vector<IntegrateValuesHelper::CalcNodes> cnodes = intValHelp.getCNodes();
            // if (cnodes.size() == 0)
            //{
            //   continue;
            //}
            calculateAverageValuesForPlane(cnodes);

            if (root) {
                real numberOfFluidsNodes = intValHelp.getNumberOfFluidsNodes();
                if (numberOfFluidsNodes > 0) {
                    ostr << j + c1o2 * dx << std::setprecision(15);

                    // mean density
                    if ((options & Density) == Density) {
                        real rho  = saRho / numberOfFluidsNodes;
                        real rhoF = saRhoF / numberOfFluidsNodes;
                        ostr << ";" << rho << ";" << rhoF;
                    }

                    // mean velocity
                    if ((options & Velocity) == Velocity) {
                        real Vx = saVx / numberOfFluidsNodes;
                        real Vy = saVy / numberOfFluidsNodes;
                        real Vz = saVz / numberOfFluidsNodes;
                        ostr << ";" << Vx << ";" << Vy << ";" << Vz;
                    }
                    // fluctuations
                    if ((options & Fluctuations) == Fluctuations) {
                        real Vxx = saVxx / numberOfFluidsNodes;
                        real Vyy = saVyy / numberOfFluidsNodes;
                        real Vzz = saVzz / numberOfFluidsNodes;
                        real Vxy = saVxy / numberOfFluidsNodes;
                        real Vxz = saVxz / numberOfFluidsNodes;
                        real Vyz = saVyz / numberOfFluidsNodes;
                        ostr << ";" << Vxx << ";" << Vyy << ";" << Vzz << ";" << Vxy << ";" << Vxz << ";" << Vyz;
                    }
                    // triple-correlations
                    if ((options & Triplecorrelations) == Triplecorrelations) {
                        real Vxxx = saVxxx / numberOfFluidsNodes;
                        real Vxxy = saVxxy / numberOfFluidsNodes;
                        real Vxxz = saVxxz / numberOfFluidsNodes;
                        real Vyyy = saVyyy / numberOfFluidsNodes;
                        real Vyyx = saVyyx / numberOfFluidsNodes;
                        real Vyyz = saVyyz / numberOfFluidsNodes;
                        real Vzzz = saVzzz / numberOfFluidsNodes;
                        real Vzzx = saVzzx / numberOfFluidsNodes;
                        real Vzzy = saVzzy / numberOfFluidsNodes;
                        real Vxyz = saVxyz / numberOfFluidsNodes;
                        ostr << ";" << Vxxx << ";" << Vxxy << ";" << Vxxz << ";" << Vyyy << ";" << Vyyx << ";" << Vyyz
                             << ";" << Vzzz << ";" << Vzzx << ";" << Vzzy << ";" << Vxyz;
                    }
                    ostr << "\n";
                }
            }
        }
        k += 2;
    }

    if (root) {
        ostr.close();
        UBLOG(logINFO, "TimeAveragedValuesSimulationObserver::planarAverage() step: " << (int)step);
    }
}
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesSimulationObserver::reset()
{
    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (SPtr<Block3D> block : blockVector[level]) {
            if (block) {
                SPtr<AverageValuesArray3D> arho = block->getKernel()->getDataSet()->getAverageDensity();
                if (arho) {
                    arho->reset(0.0);
                }

                SPtr<AverageValuesArray3D> avel = block->getKernel()->getDataSet()->getAverageVelocity();
                if (avel) {
                    avel->reset(0.0);
                }

                SPtr<AverageValuesArray3D> afl = block->getKernel()->getDataSet()->getAverageFluctuations();
                if (afl) {
                    afl->reset(0.0);
                }

                SPtr<AverageValuesArray3D> atrp = block->getKernel()->getDataSet()->getAverageTriplecorrelations();
                if (atrp) {
                    atrp->reset(0.0);
                }
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesSimulationObserver::setWithGhostLayer(bool val)
{
    withGhostLayer = val;

    if (withGhostLayer) {
        iMinC = 0;
    } else {
        iMinC = 1;
    }
}
//////////////////////////////////////////////////////////////////////////
bool TimeAveragedValuesSimulationObserver::getWithGhostLayer() { return withGhostLayer; }
//////////////////////////////////////////////////////////////////////////
void TimeAveragedValuesSimulationObserver::calculateAverageValuesForPlane(
    std::vector<IntegrateValuesHelper::CalcNodes> &cnodes)
{
    saVx = 0;
    saVy = 0;
    saVz = 0;

    saVxx = 0;
    saVyy = 0;
    saVzz = 0;
    saVxy = 0;
    saVxz = 0;
    saVyz = 0;

    saVxxx = 0;
    saVxxy = 0;
    saVxxz = 0;
    saVyyy = 0;
    saVyyx = 0;
    saVyyz = 0;
    saVzzz = 0;
    saVzzx = 0;
    saVzzy = 0;
    saVxyz = 0;

    saRho  = 0;
    saRhoF = 0;

    real lsaVx = 0;
    real lsaVy = 0;
    real lsaVz = 0;

    real lsaVxx = 0;
    real lsaVyy = 0;
    real lsaVzz = 0;
    real lsaVxy = 0;
    real lsaVxz = 0;
    real lsaVyz = 0;

    real lsaVxxx = 0;
    real lsaVxxy = 0;
    real lsaVxxz = 0;
    real lsaVyyy = 0;
    real lsaVyyx = 0;
    real lsaVyyz = 0;
    real lsaVzzz = 0;
    real lsaVzzx = 0;
    real lsaVzzy = 0;
    real lsaVxyz = 0;

    real lsaRho  = 0;
    real lsaRhoF = 0;

    for (IntegrateValuesHelper::CalcNodes cn : cnodes) {
        SPtr<ILBMKernel> kernel                               = cn.block->getKernel();
        SPtr<AverageValuesArray3D> averagedDensity            = kernel->getDataSet()->getAverageDensity();
        SPtr<AverageValuesArray3D> averagedVelocity           = kernel->getDataSet()->getAverageVelocity();
        SPtr<AverageValuesArray3D> averagedFluctuations       = kernel->getDataSet()->getAverageFluctuations();
        SPtr<AverageValuesArray3D> averagedTriplecorrelations = kernel->getDataSet()->getAverageTriplecorrelations();

        for (UbTupleInt3 node : cn.nodes) {
            real aRho  = (*averagedDensity)(Rho, val<1>(node), val<2>(node), val<3>(node));
            real aRhoF = (*averagedDensity)(RhoF, val<1>(node), val<2>(node), val<3>(node));

            real aVx = (*averagedVelocity)(Vx, val<1>(node), val<2>(node), val<3>(node));
            real aVy = (*averagedVelocity)(Vy, val<1>(node), val<2>(node), val<3>(node));
            real aVz = (*averagedVelocity)(Vz, val<1>(node), val<2>(node), val<3>(node));

            real aVxx = (*averagedFluctuations)(Vxx, val<1>(node), val<2>(node), val<3>(node));
            real aVyy = (*averagedFluctuations)(Vyy, val<1>(node), val<2>(node), val<3>(node));
            real aVzz = (*averagedFluctuations)(Vzz, val<1>(node), val<2>(node), val<3>(node));
            real aVxy = (*averagedFluctuations)(Vxy, val<1>(node), val<2>(node), val<3>(node));
            real aVxz = (*averagedFluctuations)(Vxz, val<1>(node), val<2>(node), val<3>(node));
            real aVyz = (*averagedFluctuations)(Vyz, val<1>(node), val<2>(node), val<3>(node));

            real aVxxx = (*averagedTriplecorrelations)(Vxxx, val<1>(node), val<2>(node), val<3>(node));
            real aVxxy = (*averagedTriplecorrelations)(Vxxy, val<1>(node), val<2>(node), val<3>(node));
            real aVxxz = (*averagedTriplecorrelations)(Vxxz, val<1>(node), val<2>(node), val<3>(node));
            real aVyyy = (*averagedTriplecorrelations)(Vyyy, val<1>(node), val<2>(node), val<3>(node));
            real aVyyx = (*averagedTriplecorrelations)(Vyyx, val<1>(node), val<2>(node), val<3>(node));
            real aVyyz = (*averagedTriplecorrelations)(Vyyz, val<1>(node), val<2>(node), val<3>(node));
            real aVzzz = (*averagedTriplecorrelations)(Vzzz, val<1>(node), val<2>(node), val<3>(node));
            real aVzzx = (*averagedTriplecorrelations)(Vzzx, val<1>(node), val<2>(node), val<3>(node));
            real aVzzy = (*averagedTriplecorrelations)(Vzzy, val<1>(node), val<2>(node), val<3>(node));
            real aVxyz = (*averagedTriplecorrelations)(Vxyz, val<1>(node), val<2>(node), val<3>(node));

            lsaRho += aRho;
            lsaRhoF += aRhoF;

            lsaVx += aVx;
            lsaVy += aVy;
            lsaVz += aVz;

            lsaVxx += aVxx;
            lsaVyy += aVyy;
            lsaVzz += aVzz;
            lsaVxy += aVxy;
            lsaVxz += aVxz;
            lsaVyz += aVyz;

            lsaVxxx += aVxxx;
            lsaVxxy += aVxxy;
            lsaVxxz += aVxxz;
            lsaVyyy += aVyyy;
            lsaVyyx += aVyyx;
            lsaVyyz += aVyyz;
            lsaVzzz += aVzzz;
            lsaVzzx += aVzzx;
            lsaVzzy += aVzzy;
            lsaVxyz += aVxyz;
        }
    }
    std::vector<real> values;
    std::vector<real> rvalues;

    values.push_back(lsaRho);
    values.push_back(lsaRhoF);

    values.push_back(lsaVx);
    values.push_back(lsaVy);
    values.push_back(lsaVz);

    values.push_back(lsaVxx);
    values.push_back(lsaVyy);
    values.push_back(lsaVzz);
    values.push_back(lsaVxy);
    values.push_back(lsaVxz);
    values.push_back(lsaVyz);

    values.push_back(lsaVxxx);
    values.push_back(lsaVxxy);
    values.push_back(lsaVxxz);
    values.push_back(lsaVyyy);
    values.push_back(lsaVyyx);
    values.push_back(lsaVyyz);
    values.push_back(lsaVzzz);
    values.push_back(lsaVzzx);
    values.push_back(lsaVzzy);
    values.push_back(lsaVxyz);

    rvalues = comm->gather(values);
    if (root) {
        for (int i = 0; i < (int)rvalues.size(); i += 21) {
            saRho += rvalues[i];
            saRhoF += rvalues[i + 1];

            saVx += rvalues[i + 2];
            saVy += rvalues[i + 3];
            saVz += rvalues[i + 4];

            saVxx += rvalues[i + 5];
            saVyy += rvalues[i + 6];
            saVzz += rvalues[i + 7];
            saVxy += rvalues[i + 8];
            saVxz += rvalues[i + 9];
            saVyz += rvalues[i + 10];

            saVxxx += rvalues[i + 11];
            saVxxy += rvalues[i + 12];
            saVxxz += rvalues[i + 13];
            saVyyy += rvalues[i + 14];
            saVyyx += rvalues[i + 15];
            saVyyz += rvalues[i + 16];
            saVzzz += rvalues[i + 17];
            saVzzx += rvalues[i + 18];
            saVzzy += rvalues[i + 19];
            saVxyz += rvalues[i + 20];
        }
    }
}
