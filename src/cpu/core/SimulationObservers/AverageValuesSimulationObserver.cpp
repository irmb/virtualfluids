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
#include "AverageValuesSimulationObserver.h"

#include "BCSet.h"
#include "LBMKernel.h"

#include "basics/writer/WbWriterVtkXmlASCII.h"

#include "BCArray3D.h"
#include "Block3D.h"
#include <parallel/Communicator.h>
#include "DataSet3D.h"
#include "Grid3D.h"
#include "UbScheduler.h"
#include "WbWriter.h"

using namespace std;

AverageValuesSimulationObserver::AverageValuesSimulationObserver() = default;
//////////////////////////////////////////////////////////////////////////
AverageValuesSimulationObserver::AverageValuesSimulationObserver(SPtr<Grid3D> grid, const std::string &path, WbWriter *const writer,
                                                   SPtr<UbScheduler> s, SPtr<UbScheduler> Avs,
                                                   SPtr<UbScheduler> rsMeans, SPtr<UbScheduler> rsRMS, bool restart)
    : SimulationObserver(grid, s), averageScheduler(Avs), resetSchedulerMeans(rsMeans), resetSchedulerRMS(rsRMS), path(path),
      writer(writer)
{
    resetStepMeans  = (int)rsMeans->getMinBegin();
    resetStepRMS    = (int)rsRMS->getMinBegin();
    averageInterval = (real)Avs->getMinStep();

    gridRank     = grid->getRank();
    minInitLevel = this->grid->getCoarsestInitializedLevel();
    maxInitLevel = this->grid->getFinestInitializedLevel();

    blockVector.resize(maxInitLevel + 1);

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        grid->getBlocks(level, gridRank, true, blockVector[level]);

        if (blockVector[level].size() > 0)
            compressible = blockVector[level][0]->getKernel()->getCompressible();

        if (!restart) {
            for (SPtr<Block3D> block : blockVector[level]) {
                UbTupleInt3 nx                           = grid->getBlockNX();
                SPtr<AverageValuesArray3D> averageValues = SPtr<AverageValuesArray3D>(
                    new AverageValuesArray3D(11, val<1>(nx) + 1, val<2>(nx) + 1, val<3>(nx) + 1, 0.0));
                block->getKernel()->getDataSet()->setAverageValues(averageValues);
            }
        }
    }

    // for musis special use
    // initPlotDataZ(0.0);
    // restartStep = 0.0;
}
//////////////////////////////////////////////////////////////////////////
void AverageValuesSimulationObserver::update(real step)
{
    // resetRMS(step);
    if (resetSchedulerRMS->isDue(step))
        resetDataRMS(step);

    // reset(step);
    if (resetSchedulerMeans->isDue(step))
        resetDataMeans(step);

    if (averageScheduler->isDue(step)) {
        calculateAverageValues(step);
        // for musis special use
        // collectPlotDataZ(step);
    }
    if (scheduler->isDue(step)) {
        collectData(step);
    }

    UBLOG(logDEBUG3, "AverageValuesSimulationObserver::update:" << step);
}

void AverageValuesSimulationObserver::resetDataRMS(real step)
{
    using namespace vf::basics::constant;
    
    resetStepRMS = (int)step;

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (SPtr<Block3D> block : blockVector[level]) {
            if (block) {
                SPtr<ILBMKernel> kernel                 = block->getKernel();
                SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
                SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
                SPtr<AverageValuesArray3D> av           = kernel->getDataSet()->getAverageValues();

                int minX1 = 0;
                int minX2 = 0;
                int minX3 = 0;

                int maxX1 = int(distributions->getNX1());
                int maxX2 = int(distributions->getNX2());
                int maxX3 = int(distributions->getNX3());

                for (int ix3 = minX3; ix3 < maxX3 - 1; ix3++) {
                    for (int ix2 = minX2; ix2 < maxX2 - 1; ix2++) {
                        for (int ix1 = minX1; ix1 < maxX1 - 1; ix1++) {
                            if (!bcArray->isUndefined(ix1, ix2, ix3) && !bcArray->isSolid(ix1, ix2, ix3)) {
                                //////////////////////////////////////////////////////////////////////////
                                // compute average values
                                //////////////////////////////////////////////////////////////////////////
                                (*av)(AvVxx, ix1, ix2, ix3)  = c0o1;
                                (*av)(AvVyy, ix1, ix2, ix3)  = c0o1;
                                (*av)(AvVzz, ix1, ix2, ix3)  = c0o1;
                                (*av)(AvVxy, ix1, ix2, ix3)  = c0o1;
                                (*av)(AvVxz, ix1, ix2, ix3)  = c0o1;
                                (*av)(AvVyz, ix1, ix2, ix3)  = c0o1;
                                (*av)(AvPrms, ix1, ix2, ix3) = c0o1;
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
void AverageValuesSimulationObserver::resetDataMeans(real step)
{
    using namespace vf::basics::constant;
    
    resetStepMeans = (int)step;

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (SPtr<Block3D> block : blockVector[level]) {
            if (block) {
                SPtr<ILBMKernel> kernel                 = block->getKernel();
                SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
                SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
                SPtr<AverageValuesArray3D> av           = kernel->getDataSet()->getAverageValues();

                int minX1 = 0;
                int minX2 = 0;
                int minX3 = 0;

                int maxX1 = int(distributions->getNX1());
                int maxX2 = int(distributions->getNX2());
                int maxX3 = int(distributions->getNX3());

                for (int ix3 = minX3; ix3 < maxX3 - 1; ix3++) {
                    for (int ix2 = minX2; ix2 < maxX2 - 1; ix2++) {
                        for (int ix1 = minX1; ix1 < maxX1 - 1; ix1++) {
                            if (!bcArray->isUndefined(ix1, ix2, ix3) && !bcArray->isSolid(ix1, ix2, ix3)) {
                                //////////////////////////////////////////////////////////////////////////
                                // compute average values
                                //////////////////////////////////////////////////////////////////////////
                                (*av)(AvVx, ix1, ix2, ix3) = c0o1;
                                (*av)(AvVy, ix1, ix2, ix3) = c0o1;
                                (*av)(AvVz, ix1, ix2, ix3) = c0o1;
                                (*av)(AvP, ix1, ix2, ix3)  = c0o1;
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
void AverageValuesSimulationObserver::collectData(real step)
{
    int istep = int(step);

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (SPtr<Block3D> block : blockVector[level]) {
            if (block) {
                addData(block);
            }
        }
    }

    string pfilePath, partPath, subfolder, cfilePath;
    subfolder = "av" + UbSystem::toString(istep);
    pfilePath = path + "/av/" + subfolder;
    cfilePath = path + "/av/av_collection";
    partPath  = pfilePath + "/av" + UbSystem::toString(gridRank) + "_" + UbSystem::toString(istep);

    string partName = writer->writeOctsWithNodeData(partPath, nodes, cells, datanames, data);
    size_t found    = partName.find_last_of("/");
    string piece    = partName.substr(found + 1);
    piece           = subfolder + "/" + piece;

    vector<string> cellDataNames;
    std::shared_ptr<vf::parallel::Communicator> comm = vf::parallel::Communicator::getInstance();
    vector<string> pieces   = comm->gather(piece);
    if (comm->getProcessID() == comm->getRoot()) {
        string pname =
            WbWriterVtkXmlASCII::getInstance()->writeParallelFile(pfilePath, pieces, datanames, cellDataNames);
        found = pname.find_last_of("/");
        piece = pname.substr(found + 1);

        vector<string> filenames;
        filenames.push_back(piece);
        if (step == SimulationObserver::scheduler->getMinBegin()) {
            WbWriterVtkXmlASCII::getInstance()->writeCollection(cfilePath, filenames, istep, false);
        } else {
            WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(cfilePath, filenames, istep, false);
        }
        UBLOG(logINFO, "AverageValuesSimulationObserver step: " << istep);
    }

    clearData();
}
//////////////////////////////////////////////////////////////////////////
void AverageValuesSimulationObserver::clearData()
{
    nodes.clear();
    cells.clear();
    datanames.clear();
    data.clear();
}
//////////////////////////////////////////////////////////////////////////
void AverageValuesSimulationObserver::addData(const SPtr<Block3D> block)
{
    UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);
    //    UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
    UbTupleDouble3 nodeOffset = grid->getNodeOffset(block);
    real dx                 = grid->getDeltaX(block);

    // Diese Daten werden geschrieben:
    datanames.resize(0);
    datanames.emplace_back("AvVx");
    datanames.emplace_back("AvVy");
    datanames.emplace_back("AvVz");
    datanames.emplace_back("AvVxx");
    datanames.emplace_back("AvVyy");
    datanames.emplace_back("AvVzz");
    datanames.emplace_back("AvVxy");
    datanames.emplace_back("AvVxz");
    datanames.emplace_back("AvVyz");
    datanames.emplace_back("AvP");
    datanames.emplace_back("AvPrms");

    data.resize(datanames.size());

    SPtr<ILBMKernel> kernel                 = block->getKernel();
    SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
    SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
    SPtr<AverageValuesArray3D> av           = kernel->getDataSet()->getAverageValues();
    // int ghostLayerWidth = kernel->getGhostLayerWidth();

    int minX1 = 0;
    int minX2 = 0;
    int minX3 = 0;

    int maxX1 = int(distributions->getNX1());
    int maxX2 = int(distributions->getNX2());
    int maxX3 = int(distributions->getNX3());

    // nummern vergeben und node vector erstellen + daten sammeln
    CbArray3D<int> nodeNumbers((int)maxX1, (int)maxX2, (int)maxX3, -1);

    maxX1 -= 2;
    maxX2 -= 2;
    maxX3 -= 2;

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

                    real vx = (*av)(AvVx, ix1, ix2, ix3);
                    real vy = (*av)(AvVy, ix1, ix2, ix3);
                    real vz = (*av)(AvVz, ix1, ix2, ix3);

                    real vxx = (*av)(AvVxx, ix1, ix2, ix3);
                    real vyy = (*av)(AvVyy, ix1, ix2, ix3);
                    real vzz = (*av)(AvVzz, ix1, ix2, ix3);

                    real vxy = (*av)(AvVxy, ix1, ix2, ix3);
                    real vxz = (*av)(AvVxz, ix1, ix2, ix3);
                    real vyz = (*av)(AvVyz, ix1, ix2, ix3);

                    real vp    = (*av)(AvP, ix1, ix2, ix3);
                    real vprms = (*av)(AvPrms, ix1, ix2, ix3);

                    data[index++].push_back(vx);
                    data[index++].push_back(vy);
                    data[index++].push_back(vz);

                    data[index++].push_back(vxx);
                    data[index++].push_back(vyy);
                    data[index++].push_back(vzz);

                    data[index++].push_back(vxy);
                    data[index++].push_back(vxz);
                    data[index++].push_back(vyz);

                    data[index++].push_back(vp);
                    data[index++].push_back(vprms);
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
void AverageValuesSimulationObserver::calculateAverageValues(real timeStep)
{
    using namespace D3Q27System;
    using namespace vf::basics::constant;

    // Funktionszeiger
    calcMacros = NULL;
    if (compressible) {
        calcMacros = &calcCompMacroscopicValues;
    } else {
        calcMacros = &calcIncompMacroscopicValues;
    }

    real f[27];

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (SPtr<Block3D> block : blockVector[level]) {
            if (block) {
                SPtr<ILBMKernel> kernel                 = block->getKernel();
                SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
                SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
                SPtr<AverageValuesArray3D> av           = kernel->getDataSet()->getAverageValues();

                int minX1 = 0;
                int minX2 = 0;
                int minX3 = 0;

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
                                real press = D3Q27System::calcPress(f, rho, vx, vy, vz);

                                //////////////////////////////////////////////////////////////////////////
                                // compute average values
                                //////////////////////////////////////////////////////////////////////////

                                real timeStepAfterResetRMS =
                                    (real)(timeStep - resetStepRMS) / ((real)averageInterval);
                                real timeStepAfterResetMeans =
                                    (real)(timeStep - resetStepMeans) / ((real)averageInterval);

                                // mean velocity
                                (*av)(AvVx, ix1, ix2, ix3) =
                                    ((*av)(AvVx, ix1, ix2, ix3) * timeStepAfterResetMeans + vx) /
                                    (timeStepAfterResetMeans + c1o1);
                                (*av)(AvVy, ix1, ix2, ix3) =
                                    ((*av)(AvVy, ix1, ix2, ix3) * timeStepAfterResetMeans + vy) /
                                    (timeStepAfterResetMeans + c1o1);
                                (*av)(AvVz, ix1, ix2, ix3) =
                                    ((*av)(AvVz, ix1, ix2, ix3) * timeStepAfterResetMeans + vz) /
                                    (timeStepAfterResetMeans + c1o1);

                                // rms
                                (*av)(AvVxx, ix1, ix2, ix3) =
                                    ((vx - (*av)(AvVx, ix1, ix2, ix3)) * (vx - (*av)(AvVx, ix1, ix2, ix3)) +
                                     (*av)(AvVxx, ix1, ix2, ix3) * timeStepAfterResetRMS) /
                                    (timeStepAfterResetRMS + c1o1);
                                (*av)(AvVyy, ix1, ix2, ix3) =
                                    ((vy - (*av)(AvVy, ix1, ix2, ix3)) * (vy - (*av)(AvVy, ix1, ix2, ix3)) +
                                     (*av)(AvVyy, ix1, ix2, ix3) * timeStepAfterResetRMS) /
                                    (timeStepAfterResetRMS + c1o1);
                                (*av)(AvVzz, ix1, ix2, ix3) =
                                    ((vz - (*av)(AvVz, ix1, ix2, ix3)) * (vz - (*av)(AvVz, ix1, ix2, ix3)) +
                                     (*av)(AvVzz, ix1, ix2, ix3) * timeStepAfterResetRMS) /
                                    (timeStepAfterResetRMS + c1o1);

                                // cross-correlations
                                (*av)(AvVxy, ix1, ix2, ix3) =
                                    ((vx - (*av)(AvVx, ix1, ix2, ix3)) * (vy - (*av)(AvVy, ix1, ix2, ix3)) +
                                     (*av)(AvVxy, ix1, ix2, ix3) * timeStepAfterResetRMS) /
                                    (timeStepAfterResetRMS + c1o1);
                                (*av)(AvVxz, ix1, ix2, ix3) =
                                    ((vx - (*av)(AvVx, ix1, ix2, ix3)) * (vz - (*av)(AvVz, ix1, ix2, ix3)) +
                                     (*av)(AvVxz, ix1, ix2, ix3) * timeStepAfterResetRMS) /
                                    (timeStepAfterResetRMS + c1o1);
                                (*av)(AvVyz, ix1, ix2, ix3) =
                                    ((vy - (*av)(AvVy, ix1, ix2, ix3)) * (vz - (*av)(AvVz, ix1, ix2, ix3)) +
                                     (*av)(AvVyz, ix1, ix2, ix3) * timeStepAfterResetRMS) /
                                    (timeStepAfterResetRMS + c1o1);

                                // mean and rms press
                                (*av)(AvP, ix1, ix2, ix3) =
                                    ((*av)(AvP, ix1, ix2, ix3) * timeStepAfterResetMeans + press) /
                                    (timeStepAfterResetMeans + c1o1);
                                (*av)(AvPrms, ix1, ix2, ix3) =
                                    ((press - (*av)(AvP, ix1, ix2, ix3)) * (press - (*av)(AvP, ix1, ix2, ix3)) +
                                     (*av)(AvPrms, ix1, ix2, ix3) * timeStepAfterResetRMS) /
                                    (timeStepAfterResetRMS + c1o1);

                                //////////////////////////////////////////////////////////////////////////
                            }
                        }
                    }
                }
            }
        }
    }
}
////////////////////////////////////////////////////////////////////////////
// void AverageValuesSimulationObserver::initPlotData(double step)
//{
//   std::shared_ptr<vf::parallel::Communicator> comm = vf::parallel::Communicator::getInstance();
//    if (comm->getProcessID() == comm->getRoot())
//    {
//        std::ofstream ostr;
//        string fname = path + "_PlotData_" + UbSystem::toString(step) + ".txt";
//        ostr.open(fname.c_str(), std::ios_base::out);
//        if(!ostr)
//        {
//            ostr.clear();
//            string path = UbSystem::getPathFromString(fname);
//            if(path.size()>0){ UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out);}
//            if(!ostr) throw UbException(UB_EXARGS,"couldn't open file "+fname);
//        }
//        ostr << "Time"<< "\t" <<"Ref.Time"<<"\t"<< "Z_Coor"<< "\t" << "Pore fraction" << "\t";
//        ostr << "Vx"  << "\t" << "Vy" << "\t" << "Vz" << "\t";
//        ostr << "TSx" << "\t" << "TSy"<< "\t" << "TSz"<< "TSxz";
//        ostr << endl;
//        ostr.close();
//    }
//}
//////////////////////////////////////////////////////////////////////////////
// void AverageValuesSimulationObserver::collectPlotData(double step)
//{
//
//    double hminX1 = 0.9;
//    double hminX2 = 0.0;
//    double hmaxX1 = 0.95;
//    double hmaxX2 = 0.01; //systemabmessungen world units
//
//    // 3 level platte standard:
//    double hX3_level[] = {0.305, 0.309,0.3365,0.35};
//    //0.004, 0,0365,0.045
//    //musis: 3 level refinement
//    //double hX3_level[] = {0.42, 0.28, 0.105, 0.0}; //refinement coords
//                        //bsislevel von 0.42-0.28,... (level 0 bis 2 , 3 insgesamt)
//    //musis: 4 level refinement
//    //double hX3_level[] = {0.42, 0.3, 0.195, 0.078, 0.0};
//    //musis: 5 level refinement
//    //double hX3_level[] = {0.396, 0.28, 0.18, 0.08, 0.006, 0.0};
//
//    ostringstream Str;
//    Str << step;
//    string step2string(Str.str());
//    string fname = path + "_PlotZ_" + step2string + ".txt";
//
//
//    for(int level = minInitLevel; level<=maxInitLevel;level++)
//    {
//        double dx = grid->getDeltaX(level);
//
//        for (double hi =hX3_level[level]; hi >= hX3_level[level+1]; hi=hi-dx ){
//            D3Q27IntegrateValuesHelper h1(grid, comm,
//                hminX1, hminX2, hi,
//                hmaxX1, hmaxX2, hi-dx);
//
//            h1.calculateAV();
//            double nn1 = h1.getNumberOfNodes();
//            double ns1 = h1.getNumberOfSolids();
//            if (nn1 > 0.0){
//                // get data and write into txt files
//                if (comm->getProcessID() == comm->getRoot())
//                {
//                    int istep = static_cast<int>(step);
//                    std::ofstream ostr;
//
//                    double AvVx1 = h1.getAvVx1()/nn1;
//                    double AvVx2 = h1.getAvVx2()/nn1;
//                    double AvVx3 = h1.getAvVx3()/nn1;
//
//                    double AvTSx1 = h1.getTSx1()/nn1;
//                    double AvTSx2 = h1.getTSx2()/nn1;
//                    double AvTSx3 = h1.getTSx3()/nn1;
//
//                    double AvTSx1x3 = h1.getTSx1x3()/nn1;
//
//                    ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
//                    if(!ostr)
//                    {
//                        ostr.clear();
//                        string path = UbSystem::getPathFromString(fname);
//                        if(path.size()>0){ UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out |
//std::ios_base::app);}                         if(!ostr) throw UbException(UB_EXARGS,"couldn't open file "+fname);
//                    }
//                    ostr << istep << "\t" << resetStep << "\t" << hi+0.5*dx << "\t" << nn1/(nn1+ns1)*100.0 << "%\t";
//                    ostr << AvVx1 << "\t" << AvVx2 << "\t" << AvVx3 << "\t";
//                    ostr << AvTSx1<< "\t" << AvTSx2<< "\t" << AvTSx3<< "\t" << AvTSx1x3;
//                    ostr << endl;
//                    ostr.close();
//
//                }
//            }
//        }
//
//    }
//}

//! \}
