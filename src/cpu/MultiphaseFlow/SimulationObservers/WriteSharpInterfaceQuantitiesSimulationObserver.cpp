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
//! \file WriteSharpInterfaceQuantitiesSimulationObserver.cpp
//! \ingroup SimulationObservers
//! \author Konstantin Kutscher
//=======================================================================================

#include "WriteSharpInterfaceQuantitiesSimulationObserver.h"
#include "BCSet.h"
#include "LBMKernel.h"
#include <string>
#include <vector>

#include "BCArray3D.h"
#include "Block3D.h"
#include <mpi/Communicator.h>
#include "DataSet3D.h"
#include "Grid3D.h"
#include "LBMUnitConverter.h"
#include "UbScheduler.h"
#include "basics/writer/WbWriterVtkXmlASCII.h"

WriteSharpInterfaceQuantitiesSimulationObserver::WriteSharpInterfaceQuantitiesSimulationObserver() = default;
//////////////////////////////////////////////////////////////////////////
WriteSharpInterfaceQuantitiesSimulationObserver::WriteSharpInterfaceQuantitiesSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s,
                                                                             const std::string &path,
                                                                             WbWriter *const writer,
                                                                             SPtr<LBMUnitConverter> conv,
                                                                             std::shared_ptr<vf::mpi::Communicator> comm)
        : SimulationObserver(grid, s), path(path), writer(writer), conv(conv), comm(comm)
{
    gridRank = comm->getProcessID();
    minInitLevel = this->grid->getCoarsestInitializedLevel();
    maxInitLevel = this->grid->getFinestInitializedLevel();

    blockVector.resize(maxInitLevel + 1);

    for (int level = minInitLevel; level <= maxInitLevel; level++)
    {
        grid->getBlocks(level, gridRank, true, blockVector[level]);
    }
}

//////////////////////////////////////////////////////////////////////////
void WriteSharpInterfaceQuantitiesSimulationObserver::init()
{}

//////////////////////////////////////////////////////////////////////////
void WriteSharpInterfaceQuantitiesSimulationObserver::update(double step)
{
    if (scheduler->isDue(step))
        collectData(step);

    UBLOG(logDEBUG3, "WriteSharpInterfaceQuantitiesSimulationObserver::update:" << step);
}

//////////////////////////////////////////////////////////////////////////
void WriteSharpInterfaceQuantitiesSimulationObserver::collectData(double step)
{
    int istep = static_cast<int>(step);

    for (int level = minInitLevel; level <= maxInitLevel; level++)
    {
        for (SPtr<Block3D> block : blockVector[level])
        {
            if (block)
            {
                addDataMQ(block);
            }
        }
    }

    std::string pfilePath, partPath, subfolder, cfilePath;

    subfolder = "mq" + UbSystem::toString(istep);
    pfilePath = path + "/mq/" + subfolder;
    cfilePath = path + "/mq/mq_collection";
    partPath = pfilePath + "/mq" + UbSystem::toString(gridRank) + "_" + UbSystem::toString(istep);

    std::string partName = writer->writeOctsWithNodeData(partPath, nodes, cells, datanames, data);
    size_t found = partName.find_last_of("/");
    std::string piece = partName.substr(found + 1);
    piece = subfolder + "/" + piece;

    std::vector<std::string> cellDataNames;
    std::vector<std::string> pieces = comm->gather(piece);
    if (comm->getProcessID() == comm->getRoot()) {
        std::string pname =
                WbWriterVtkXmlASCII::getInstance()->writeParallelFile(pfilePath, pieces, datanames, cellDataNames);
        found = pname.find_last_of("/");
        piece = pname.substr(found + 1);

        std::vector<std::string> filenames;
        filenames.push_back(piece);
        if (step == SimulationObserver::scheduler->getMinBegin())
        {
            WbWriterVtkXmlASCII::getInstance()->writeCollection(cfilePath, filenames, istep, false);
        } else
        {
            WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(cfilePath, filenames, istep, false);
        }
        UBLOG(logINFO, "WriteSharpInterfaceQuantitiesSimulationObserver step: " << istep);
    }

    clearData();
}

//////////////////////////////////////////////////////////////////////////
void WriteSharpInterfaceQuantitiesSimulationObserver::clearData()
{
    nodes.clear();
    cells.clear();
    datanames.clear();
    data.clear();
}

//////////////////////////////////////////////////////////////////////////
void WriteSharpInterfaceQuantitiesSimulationObserver::addDataMQ(SPtr<Block3D> block)
{
    using namespace D3Q27System;
    //using namespace UbMath;
    using namespace vf::lbm::dir;
    using namespace vf::basics::constant;
    SPtr<LBMKernel> kernel = dynamicPointerCast<LBMKernel>(block->getKernel());
    //double level   = (double)block->getLevel();
    kernel->swapDistributions();
    // Diese Daten werden geschrieben:
    datanames.resize(0);
    datanames.push_back("Phi");
    datanames.push_back("Vx");
    datanames.push_back("Vy");
    datanames.push_back("Vz");
    datanames.push_back("Pressure");

    data.resize(datanames.size());


    SPtr<BCArray3D> bcArray                  = kernel->getBCSet()->getBCArray();
    SPtr<DistributionArray3D> distributionsF = kernel->getDataSet()->getFdistributions();
    SPtr<DistributionArray3D> distributionsH = kernel->getDataSet()->getHdistributions();
    SPtr<DistributionArray3D> distributionsH2 = kernel->getDataSet()->getH2distributions();
    SPtr<PhaseFieldArray3D> divU             = kernel->getDataSet()->getPhaseField();

    real pressure;

    real f[D3Q27System::ENDF + 1];
    real phi;
    real vx1, vx2, vx3, rho;
    real densityRatio = kernel->getDensityRatio();
    real phiL = kernel->getPhiL();
    real phiH = kernel->getPhiH();

    // knotennummerierung faengt immer bei 0 an!
    int SWB, SEB, NEB, NWB, SWT, SET, NET, NWT;

    if (block->getKernel()->getCompressible()) {
        calcMacros = &D3Q27System::calcCompMacroscopicValues;
    } else {
        calcMacros = &D3Q27System::calcIncompMacroscopicValues;
    }

    // int minX1 = 0;
    // int minX2 = 0;
    // int minX3 = 0;

    int maxX1 = (int)(distributionsF->getNX1());
    int maxX2 = (int)(distributionsF->getNX2());
    int maxX3 = (int)(distributionsF->getNX3());

    int minX1 = 0;
    int minX2 = 0;
    int minX3 = 0;
    
    if (kernel->getGhostLayerWidth() == 2)
    {
        minX1 = 1;
        minX2 = 1;
        minX3 = 1;
    }

    // int maxX1 = (int)(distributions->getNX1());
    // int maxX2 = (int)(distributions->getNX2());
    // int maxX3 = (int)(distributions->getNX3());

    // nummern vergeben und node vector erstellen + daten sammeln
    CbArray3D<int> nodeNumbers((int)maxX1, (int)maxX2, (int)maxX3, -1);
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr phaseField(
        new CbArray3D<real, IndexerX3X2X1>(maxX1, maxX2, maxX3, -999.0));


    for (int ix3 = minX3; ix3 < maxX3; ix3++) {
        for (int ix2 = minX2; ix2 < maxX2; ix2++) {
            for (int ix1 = minX1; ix1 < maxX1; ix1++) {
                if (!bcArray->isUndefined(ix1, ix2, ix3) && !bcArray->isSolid(ix1, ix2, ix3)) {
                    distributionsH->getDistributionInv(f, ix1, ix2, ix3);
                    (*phaseField)(ix1, ix2, ix3) =
                        ((f[DIR_PPP] + f[DIR_MMM]) + (f[DIR_PMP] + f[DIR_MPM])) + ((f[DIR_PMM] + f[DIR_MPP]) + (f[DIR_MMP] + f[DIR_PPM])) +
                        (((f[DIR_PP0] + f[DIR_MM0]) + (f[DIR_PM0] + f[DIR_MP0])) + ((f[DIR_P0P] + f[DIR_M0M]) + (f[DIR_P0M] + f[DIR_M0P])) +
                        ((f[DIR_0PM] + f[DIR_0MP]) + (f[DIR_0PP] + f[DIR_0MM]))) +
                            ((f[DIR_P00] + f[DIR_M00]) + (f[DIR_0P0] + f[DIR_0M0]) + (f[DIR_00P] + f[DIR_00M])) + f[DIR_000];
                }
            }
        }
    }

    if (kernel->getGhostLayerWidth() == 1)
    {
        maxX1 -= 2;
        maxX2 -= 2;
        maxX3 -= 2;
    }
    else if (kernel->getGhostLayerWidth() == 2)
    {
        maxX1 -= 3;
        maxX2 -= 3;
        maxX3 -= 3;
    }

    int nr = (int)nodes.size();
    //real dX1_phi;
    //real dX2_phi;
    //real dX3_phi;
    //real mu;

    for (int ix3 = minX3; ix3 <= maxX3; ix3++) {
        for (int ix2 = minX2; ix2 <= maxX2; ix2++) {
            for (int ix1 = minX1; ix1 <= maxX1; ix1++) {
                if (!bcArray->isUndefined(ix1, ix2, ix3) && !bcArray->isSolid(ix1, ix2, ix3)) {
                    int index                  = 0;
                    nodeNumbers(ix1, ix2, ix3) = nr++;
                    Vector3D worldCoordinates  = grid->getNodeCoordinates(block, ix1, ix2, ix3);
                    nodes.push_back(UbTupleFloat3(float(worldCoordinates[0]), float(worldCoordinates[1]),
                                                  float(worldCoordinates[2])));

                    phi = (*phaseField)(ix1, ix2, ix3);


                    distributionsF->getDistributionInv(f, ix1, ix2, ix3);
                    //real dU = (*divU)(ix1, ix2, ix3);

                    real rhoH = 1.0;
                    real rhoL = 1.0 / densityRatio;
                    // real rhoToPhi = (1.0 - 1.0/densityRatio);
                    //real rhoToPhi = (rhoH - rhoL) / (phiH - phiL);

                    // rho = phi[ZERO] + (1.0 - phi[ZERO])*1.0/densityRatio;

                    
                    rho = (phi>c1o2) ? rhoH : rhoL; // rhoH + rhoToPhi * (phi - phiH);

                        vx1 =
                            ((((f[DIR_PPP] - f[DIR_MMM]) + (f[DIR_PMP] - f[DIR_MPM])) + ((f[DIR_PMM] - f[DIR_MPP]) + (f[DIR_PPM] - f[DIR_MMP]))) +
                            (((f[DIR_P0M] - f[DIR_M0P]) + (f[DIR_P0P] - f[DIR_M0M])) + ((f[DIR_PM0] - f[DIR_MP0]) + (f[DIR_PP0] - f[DIR_MM0]))) + (f[DIR_P00] - f[DIR_M00])) ;

                        vx2 =
                            ((((f[DIR_PPP] - f[DIR_MMM]) + (f[DIR_MPM] - f[DIR_PMP])) + ((f[DIR_MPP] - f[DIR_PMM]) + (f[DIR_PPM] - f[DIR_MMP]))) +
                            (((f[DIR_0PM] - f[DIR_0MP]) + (f[DIR_0PP] - f[DIR_0MM])) + ((f[DIR_MP0] - f[DIR_PM0]) + (f[DIR_PP0] - f[DIR_MM0]))) + (f[DIR_0P0] - f[DIR_0M0])) ;

                        vx3 =
                            ((((f[DIR_PPP] - f[DIR_MMM]) + (f[DIR_PMP] - f[DIR_MPM])) + ((f[DIR_MPP] - f[DIR_PMM]) + (f[DIR_MMP] - f[DIR_PPM]))) +
                            (((f[DIR_0MP] - f[DIR_0PM]) + (f[DIR_0PP] - f[DIR_0MM])) + ((f[DIR_M0P] - f[DIR_P0M]) + (f[DIR_P0P] - f[DIR_M0M]))) + (f[DIR_00P] - f[DIR_00M]));



                    pressure = (((f[DIR_PPP] + f[DIR_MMM]) + (f[DIR_PMP] + f[DIR_MPM])) + ((f[DIR_PMM] + f[DIR_MPP]) + (f[DIR_MMP] + f[DIR_PPM])) +
                        (((f[DIR_PP0] + f[DIR_MM0]) + (f[DIR_PM0] + f[DIR_MP0])) + ((f[DIR_P0P] + f[DIR_M0M]) + (f[DIR_P0M] + f[DIR_M0P])) +
                        ((f[DIR_0PM] + f[DIR_0MP]) + (f[DIR_0PP] + f[DIR_0MM]))) +
                            ((f[DIR_P00] + f[DIR_M00]) + (f[DIR_0P0] + f[DIR_0M0]) + (f[DIR_00P] + f[DIR_00M])) + f[DIR_000])*c1o3*rho;

					if (UbMath::isNaN(vx1) || UbMath::isInfinity(vx1))
                        UB_THROW(UbException(
                            UB_EXARGS, "vx1 is not a number (nan or -1.#IND) or infinity number -1.#INF in block=" +
                                           block->toString() + ", node=" + UbSystem::toString(ix1) + "," +
                                           UbSystem::toString(ix2) + "," + UbSystem::toString(ix3)));
                    // vx1=999.0;
                    if (UbMath::isNaN(vx2) || UbMath::isInfinity(vx2))
                        UB_THROW(UbException(
                            UB_EXARGS, "vx2 is not a number (nan or -1.#IND) or infinity number -1.#INF in block=" +
                                           block->toString() + ", node=" + UbSystem::toString(ix1) + "," +
                                           UbSystem::toString(ix2) + "," + UbSystem::toString(ix3)));
                    // vx2=999.0;
                    if (UbMath::isNaN(vx3) || UbMath::isInfinity(vx3))
                        UB_THROW(UbException(
                            UB_EXARGS, "vx3 is not a number (nan or -1.#IND) or infinity number -1.#INF in block=" +
                                           block->toString() + ", node=" + UbSystem::toString(ix1) + "," +
                                           UbSystem::toString(ix2) + "," + UbSystem::toString(ix3)));

                    if (UbMath::isNaN(phi) || UbMath::isInfinity(phi))
                        UB_THROW(UbException(
                            UB_EXARGS, "phi is not a number (nan or -1.#IND) or infinity number -1.#INF in block=" +
                                           block->toString() + ", node=" + UbSystem::toString(ix1) + "," +
                                           UbSystem::toString(ix2) + "," + UbSystem::toString(ix3)));

                     if (UbMath::isNaN(pressure) || UbMath::isInfinity(pressure))
                        UB_THROW( UbException(UB_EXARGS,"pressure is not a number (nan or -1.#IND) or infinity number -1.#INF in block="+block->toString()+
                        ", node="+UbSystem::toString(ix1)+","+UbSystem::toString(ix2)+","+UbSystem::toString(ix3)));

                    //if (UbMath::isNaN(p1) || UbMath::isInfinity(p1))
                    //    UB_THROW( UbException(UB_EXARGS,"p1 is not a number (nan or -1.#IND) or infinity number -1.#INF in block="+block->toString()+
                    //    ", node="+UbSystem::toString(ix1)+","+UbSystem::toString(ix2)+","+UbSystem::toString(ix3)));

                    //if (UbMath::isNaN(mp) || UbMath::isInfinity(mp))
                    //    UB_THROW(UbException(UB_EXARGS, "mp is not a number (nan or -1.#IND) or infinity number -1.#INF in block=" + block->toString() +
                    //        ", node=" + UbSystem::toString(ix1) + "," + UbSystem::toString(ix2) + "," + UbSystem::toString(ix3)));

                    //if (UbMath::isNaN(delmp) || UbMath::isInfinity(delmp))
                    //    UB_THROW(UbException(UB_EXARGS, "delmp is not a number (nan or -1.#IND) or infinity number -1.#INF in block=" + block->toString() +
                    //        ", node=" + UbSystem::toString(ix1) + "," + UbSystem::toString(ix2) + "," + UbSystem::toString(ix3)));

                    data[index++].push_back(phi);
                    data[index++].push_back(vx1);
                    data[index++].push_back(vx2);
                    data[index++].push_back(vx3);
                    //data[index++].push_back(p1);
                    //data[index++].push_back(phi2[DIR_000]);
                    //data[index++].push_back(mp);
                    //data[index++].push_back(delmp);
                    data[index++].push_back(pressure);
                }
            }
        }
    }
    maxX1 -= 1;
    maxX2 -= 1;
    maxX3 -= 1;
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
    kernel->swapDistributions();
}

real WriteSharpInterfaceQuantitiesSimulationObserver::gradX1_phi(const real *const &h)
{
    using namespace D3Q27System;
	using namespace vf::lbm::dir;
    using namespace vf::basics::constant;
    real sum = c0o1;
    for (int k = FSTARTDIR; k <= FENDDIR; k++) {
        sum += WEIGTH[k] * DX1[k] * h[k];
    }
    return 3.0 * sum;
}
real WriteSharpInterfaceQuantitiesSimulationObserver::gradX2_phi(const real *const &h)
{
    using namespace D3Q27System;
	using namespace vf::lbm::dir;
    using namespace vf::basics::constant;
    real sum = c0o1;
    for (int k = FSTARTDIR; k <= FENDDIR; k++) {
        sum += WEIGTH[k] * DX2[k] * h[k];
    }
    return 3.0 * sum;
}

real WriteSharpInterfaceQuantitiesSimulationObserver::gradX3_phi(const real *const &h)
{
    using namespace D3Q27System;
    using namespace vf::lbm::dir;
    using namespace vf::basics::constant;
    real sum = c0o1;
    for (int k = FSTARTDIR; k <= FENDDIR; k++) {
        sum += WEIGTH[k] * DX3[k] * h[k];
    }
    return 3.0 * sum;
}

real WriteSharpInterfaceQuantitiesSimulationObserver::nabla2_phi(const real *const &h)
{
    using namespace D3Q27System;
    using namespace vf::lbm::dir;
    using namespace vf::basics::constant;
    real sum = c0o1;
    for (int k = FSTARTDIR; k <= FENDDIR; k++) {
        sum += WEIGTH[k] * (h[k] - h[DIR_000]);
    }
    return 6.0 * sum;
}