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
//! \file WriteMultiphaseQuantitiesCoProcessor.cpp
//! \ingroup CoProcessors
//! \author Konstantin Kutscher
//=======================================================================================

#include "WriteMultiphaseQuantitiesCoProcessor.h"
#include "BCProcessor.h"
#include "LBMKernel.h"
#include <string>
#include <vector>
#include "MultiphaseTwoPhaseFieldsVelocityCumulantLBMKernel.h"

#include "BCArray3D.h"
#include "Block3D.h"
#include "Communicator.h"
#include "DataSet3D.h"
#include "Grid3D.h"
#include "LBMUnitConverter.h"
#include "UbScheduler.h"
#include "basics/writer/WbWriterVtkXmlASCII.h"

WriteMultiphaseQuantitiesCoProcessor::WriteMultiphaseQuantitiesCoProcessor() = default;
//////////////////////////////////////////////////////////////////////////
WriteMultiphaseQuantitiesCoProcessor::WriteMultiphaseQuantitiesCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s,
                                                                             const std::string &path,
                                                                             WbWriter *const writer,
                                                                             SPtr<LBMUnitConverter> conv,
                                                                             SPtr<Communicator> comm)
        : CoProcessor(grid, s), path(path), writer(writer), conv(conv), comm(comm)
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
void WriteMultiphaseQuantitiesCoProcessor::init()
{}

//////////////////////////////////////////////////////////////////////////
void WriteMultiphaseQuantitiesCoProcessor::process(double step)
{
    if (scheduler->isDue(step))
        collectData(step);

    UBLOG(logDEBUG3, "WriteMultiphaseQuantitiesCoProcessor::update:" << step);
}

//////////////////////////////////////////////////////////////////////////
void WriteMultiphaseQuantitiesCoProcessor::collectData(double step)
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
        if (step == CoProcessor::scheduler->getMinBegin())
        {
            WbWriterVtkXmlASCII::getInstance()->writeCollection(cfilePath, filenames, istep, false);
        } else
        {
            WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(cfilePath, filenames, istep, false);
        }
        UBLOG(logINFO, "WriteMultiphaseQuantitiesCoProcessor step: " << istep);
    }

    clearData();
}

//////////////////////////////////////////////////////////////////////////
void WriteMultiphaseQuantitiesCoProcessor::clearData()
{
    nodes.clear();
    cells.clear();
    datanames.clear();
    data.clear();
}

//////////////////////////////////////////////////////////////////////////
void WriteMultiphaseQuantitiesCoProcessor::addDataMQ(SPtr<Block3D> block)
{
    using namespace D3Q27System;
    using namespace UbMath;
    SPtr<LBMKernel> kernel = dynamicPointerCast<LBMKernel>(block->getKernel());
    //double level   = (double)block->getLevel();

    // Diese Daten werden geschrieben:
    datanames.resize(0);
    datanames.push_back("Phi");
    datanames.push_back("Vx");
    datanames.push_back("Vy");
    datanames.push_back("Vz");
    datanames.push_back("P1");
    datanames.push_back("Phi2");
    if (dynamicPointerCast<MultiphaseTwoPhaseFieldsVelocityCumulantLBMKernel>(kernel)->pressure) datanames.push_back("Pressure");

    data.resize(datanames.size());


    SPtr<BCArray3D> bcArray                  = kernel->getBCProcessor()->getBCArray();
    SPtr<DistributionArray3D> distributionsF = kernel->getDataSet()->getFdistributions();
    SPtr<DistributionArray3D> distributionsH = kernel->getDataSet()->getHdistributions();
    SPtr<DistributionArray3D> distributionsH2 = kernel->getDataSet()->getH2distributions();
    SPtr<PhaseFieldArray3D> divU             = kernel->getDataSet()->getPhaseField();

    LBMReal f[D3Q27System::ENDF + 1];
    LBMReal phi[D3Q27System::ENDF + 1];
    LBMReal phi2[D3Q27System::ENDF + 1];
    LBMReal vx1, vx2, vx3, rho, p1, beta, kappa;
    LBMReal densityRatio = kernel->getDensityRatio();

    kernel->getMultiphaseModelParameters(beta, kappa);
    LBMReal phiL = kernel->getPhiL();
    LBMReal phiH = kernel->getPhiH();

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

    // int maxX1 = (int)(distributions->getNX1());
    // int maxX2 = (int)(distributions->getNX2());
    // int maxX3 = (int)(distributions->getNX3());

    // nummern vergeben und node vector erstellen + daten sammeln
    CbArray3D<int> nodeNumbers((int)maxX1, (int)maxX2, (int)maxX3, -1);
    CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr phaseField(
        new CbArray3D<LBMReal, IndexerX3X2X1>(maxX1, maxX2, maxX3, -999.0));
    CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr phaseField2(
        new CbArray3D<LBMReal, IndexerX3X2X1>(maxX1, maxX2, maxX3, -999.0));

    for (int ix3 = minX3; ix3 < maxX3; ix3++) {
        for (int ix2 = minX2; ix2 < maxX2; ix2++) {
            for (int ix1 = minX1; ix1 < maxX1; ix1++) {
                if (!bcArray->isUndefined(ix1, ix2, ix3) && !bcArray->isSolid(ix1, ix2, ix3)) {
                    distributionsH->getDistribution(f, ix1, ix2, ix3);
                    (*phaseField)(ix1, ix2, ix3) =
                        ((f[TNE] + f[BSW]) + (f[TSE] + f[BNW])) + ((f[BSE] + f[TNW]) + (f[TSW] + f[BNE])) +
                        (((f[NE] + f[SW]) + (f[SE] + f[NW])) + ((f[TE] + f[BW]) + (f[BE] + f[TW])) +
                        ((f[BN] + f[TS]) + (f[TN] + f[BS]))) +
                            ((f[E] + f[W]) + (f[N] + f[S]) + (f[T] + f[B])) + f[REST];
                    if (distributionsH2) {
                    distributionsH2->getDistribution(f, ix1, ix2, ix3);
                    (*phaseField2)(ix1, ix2, ix3) =
                        ((f[TNE] + f[BSW]) + (f[TSE] + f[BNW])) + ((f[BSE] + f[TNW]) + (f[TSW] + f[BNE])) +
                        (((f[NE] + f[SW]) + (f[SE] + f[NW])) + ((f[TE] + f[BW]) + (f[BE] + f[TW])) +
                        ((f[BN] + f[TS]) + (f[TN] + f[BS]))) +
                            ((f[E] + f[W]) + (f[N] + f[S]) + (f[T] + f[B])) + f[REST];
                }
                    else { (*phaseField2)(ix1, ix2, ix3) = 999.0; }
                    
                }
            }
        }
    }

    maxX1 -= 2;
    maxX2 -= 2;
    maxX3 -= 2;

    // maxX1 -= 1;
    // maxX2 -= 1;
    // maxX3 -= 1;

    // D3Q27BoundaryConditionPtr bcPtr;
    int nr = (int)nodes.size();
    LBMReal dX1_phi;
    LBMReal dX2_phi;
    LBMReal dX3_phi;
    LBMReal mu;

    for (int ix3 = minX3; ix3 <= maxX3; ix3++) {
        for (int ix2 = minX2; ix2 <= maxX2; ix2++) {
            for (int ix1 = minX1; ix1 <= maxX1; ix1++) {
                if (!bcArray->isUndefined(ix1, ix2, ix3) && !bcArray->isSolid(ix1, ix2, ix3)) {
                    int index                  = 0;
                    nodeNumbers(ix1, ix2, ix3) = nr++;
                    Vector3D worldCoordinates  = grid->getNodeCoordinates(block, ix1, ix2, ix3);
                    nodes.push_back(UbTupleFloat3(float(worldCoordinates[0]), float(worldCoordinates[1]),
                                                  float(worldCoordinates[2])));

                    phi[REST] = (*phaseField)(ix1, ix2, ix3);
                    phi2[REST] = (*phaseField2)(ix1, ix2, ix3);

                    if ((ix1 == 0) || (ix2 == 0) || (ix3 == 0)) {
                        dX1_phi = 0.0;
                        dX2_phi = 0.0;
                        dX3_phi = 0.0;
                        mu      = 0.0;
                        // vx1 = 0.0;
                        // vx2 = 0.0;
                        // vx3 = 0.0;
                    } else {
                        phi[E]   = (*phaseField)(ix1 + DX1[E], ix2 + DX2[E], ix3 + DX3[E]);
                        phi[N]   = (*phaseField)(ix1 + DX1[N], ix2 + DX2[N], ix3 + DX3[N]);
                        phi[T]   = (*phaseField)(ix1 + DX1[T], ix2 + DX2[T], ix3 + DX3[T]);
                        phi[W]   = (*phaseField)(ix1 + DX1[W], ix2 + DX2[W], ix3 + DX3[W]);
                        phi[S]   = (*phaseField)(ix1 + DX1[S], ix2 + DX2[S], ix3 + DX3[S]);
                        phi[B]   = (*phaseField)(ix1 + DX1[B], ix2 + DX2[B], ix3 + DX3[B]);
                        phi[NE]  = (*phaseField)(ix1 + DX1[NE], ix2 + DX2[NE], ix3 + DX3[NE]);
                        phi[NW]  = (*phaseField)(ix1 + DX1[NW], ix2 + DX2[NW], ix3 + DX3[NW]);
                        phi[TE]  = (*phaseField)(ix1 + DX1[TE], ix2 + DX2[TE], ix3 + DX3[TE]);
                        phi[TW]  = (*phaseField)(ix1 + DX1[TW], ix2 + DX2[TW], ix3 + DX3[TW]);
                        phi[TN]  = (*phaseField)(ix1 + DX1[TN], ix2 + DX2[TN], ix3 + DX3[TN]);
                        phi[TS]  = (*phaseField)(ix1 + DX1[TS], ix2 + DX2[TS], ix3 + DX3[TS]);
                        phi[SW]  = (*phaseField)(ix1 + DX1[SW], ix2 + DX2[SW], ix3 + DX3[SW]);
                        phi[SE]  = (*phaseField)(ix1 + DX1[SE], ix2 + DX2[SE], ix3 + DX3[SE]);
                        phi[BW]  = (*phaseField)(ix1 + DX1[BW], ix2 + DX2[BW], ix3 + DX3[BW]);
                        phi[BE]  = (*phaseField)(ix1 + DX1[BE], ix2 + DX2[BE], ix3 + DX3[BE]);
                        phi[BS]  = (*phaseField)(ix1 + DX1[BS], ix2 + DX2[BS], ix3 + DX3[BS]);
                        phi[BN]  = (*phaseField)(ix1 + DX1[BN], ix2 + DX2[BN], ix3 + DX3[BN]);
                        phi[BSW] = (*phaseField)(ix1 + DX1[BSW], ix2 + DX2[BSW], ix3 + DX3[BSW]);
                        phi[BSE] = (*phaseField)(ix1 + DX1[BSE], ix2 + DX2[BSE], ix3 + DX3[BSE]);
                        phi[BNW] = (*phaseField)(ix1 + DX1[BNW], ix2 + DX2[BNW], ix3 + DX3[BNW]);
                        phi[BNE] = (*phaseField)(ix1 + DX1[BNE], ix2 + DX2[BNE], ix3 + DX3[BNE]);
                        phi[TNE] = (*phaseField)(ix1 + DX1[TNE], ix2 + DX2[TNE], ix3 + DX3[TNE]);
                        phi[TNW] = (*phaseField)(ix1 + DX1[TNW], ix2 + DX2[TNW], ix3 + DX3[TNW]);
                        phi[TSE] = (*phaseField)(ix1 + DX1[TSE], ix2 + DX2[TSE], ix3 + DX3[TSE]);
                        phi[TSW] = (*phaseField)(ix1 + DX1[TSW], ix2 + DX2[TSW], ix3 + DX3[TSW]);
                        dX1_phi  = 0.0 * gradX1_phi(phi);
                        dX2_phi  = 0.0 * gradX2_phi(phi);
                        dX3_phi  = 0.0 * gradX3_phi(phi);
                        mu = 2 * beta * phi[REST] * (phi[REST] - 1) * (2 * phi[REST] - 1) - kappa * nabla2_phi(phi);

                        //phi2[E] = (*phaseField2)(ix1 + DX1[E], ix2 + DX2[E], ix3 + DX3[E]);
                        //phi2[N] = (*phaseField2)(ix1 + DX1[N], ix2 + DX2[N], ix3 + DX3[N]);
                        //phi2[T] = (*phaseField2)(ix1 + DX1[T], ix2 + DX2[T], ix3 + DX3[T]);
                        //phi2[W] = (*phaseField2)(ix1 + DX1[W], ix2 + DX2[W], ix3 + DX3[W]);
                        //phi2[S] = (*phaseField2)(ix1 + DX1[S], ix2 + DX2[S], ix3 + DX3[S]);
                        //phi2[B] = (*phaseField2)(ix1 + DX1[B], ix2 + DX2[B], ix3 + DX3[B]);
                        //phi2[NE] = (*phaseField2)(ix1 + DX1[NE], ix2 + DX2[NE], ix3 + DX3[NE]);
                        //phi2[NW] = (*phaseField2)(ix1 + DX1[NW], ix2 + DX2[NW], ix3 + DX3[NW]);
                        //phi2[TE] = (*phaseField2)(ix1 + DX1[TE], ix2 + DX2[TE], ix3 + DX3[TE]);
                        //phi2[TW] = (*phaseField2)(ix1 + DX1[TW], ix2 + DX2[TW], ix3 + DX3[TW]);
                        //phi2[TN] = (*phaseField2)(ix1 + DX1[TN], ix2 + DX2[TN], ix3 + DX3[TN]);
                        //phi2[TS] = (*phaseField2)(ix1 + DX1[TS], ix2 + DX2[TS], ix3 + DX3[TS]);
                        //phi2[SW] = (*phaseField2)(ix1 + DX1[SW], ix2 + DX2[SW], ix3 + DX3[SW]);
                        //phi2[SE] = (*phaseField2)(ix1 + DX1[SE], ix2 + DX2[SE], ix3 + DX3[SE]);
                        //phi2[BW] = (*phaseField2)(ix1 + DX1[BW], ix2 + DX2[BW], ix3 + DX3[BW]);
                        //phi2[BE] = (*phaseField2)(ix1 + DX1[BE], ix2 + DX2[BE], ix3 + DX3[BE]);
                        //phi2[BS] = (*phaseField2)(ix1 + DX1[BS], ix2 + DX2[BS], ix3 + DX3[BS]);
                        //phi2[BN] = (*phaseField2)(ix1 + DX1[BN], ix2 + DX2[BN], ix3 + DX3[BN]);
                        //phi2[BSW] = (*phaseField2)(ix1 + DX1[BSW], ix2 + DX2[BSW], ix3 + DX3[BSW]);
                        //phi2[BSE] = (*phaseField2)(ix1 + DX1[BSE], ix2 + DX2[BSE], ix3 + DX3[BSE]);
                        //phi2[BNW] = (*phaseField2)(ix1 + DX1[BNW], ix2 + DX2[BNW], ix3 + DX3[BNW]);
                        //phi2[BNE] = (*phaseField2)(ix1 + DX1[BNE], ix2 + DX2[BNE], ix3 + DX3[BNE]);
                        //phi2[TNE] = (*phaseField2)(ix1 + DX1[TNE], ix2 + DX2[TNE], ix3 + DX3[TNE]);
                        //phi2[TNW] = (*phaseField2)(ix1 + DX1[TNW], ix2 + DX2[TNW], ix3 + DX3[TNW]);
                        //phi2[TSE] = (*phaseField2)(ix1 + DX1[TSE], ix2 + DX2[TSE], ix3 + DX3[TSE]);
                        //phi2[TSW] = (*phaseField2)(ix1 + DX1[TSW], ix2 + DX2[TSW], ix3 + DX3[TSW]);

                       // mu = 2 * beta * phi[REST] * (phi[REST] - 1) * (2 * phi[REST] - 1) - kappa * nabla2_phi(phi);



                    }

                    distributionsF->getDistribution(f, ix1, ix2, ix3);
                    //LBMReal dU = (*divU)(ix1, ix2, ix3);

                    LBMReal rhoH = 1.0;
                    LBMReal rhoL = 1.0 / densityRatio;
                    // LBMReal rhoToPhi = (1.0 - 1.0/densityRatio);
                    LBMReal rhoToPhi = (rhoH - rhoL) / (phiH - phiL);

                    // rho = phi[ZERO] + (1.0 - phi[ZERO])*1.0/densityRatio;
                    rho = rhoH + rhoToPhi * (phi[REST] - phiH);

                    if (dynamicPointerCast<MultiphaseTwoPhaseFieldsVelocityCumulantLBMKernel>(kernel)->pressure) {
                        vx1 =
                            ((((f[TNE] - f[BSW]) + (f[TSE] - f[BNW])) + ((f[BSE] - f[TNW]) + (f[BNE] - f[TSW]))) +
                            (((f[BE] - f[TW]) + (f[TE] - f[BW])) + ((f[SE] - f[NW]) + (f[NE] - f[SW]))) + (f[E] - f[W])) ;

                        vx2 =
                            ((((f[TNE] - f[BSW]) + (f[BNW] - f[TSE])) + ((f[TNW] - f[BSE]) + (f[BNE] - f[TSW]))) +
                            (((f[BN] - f[TS]) + (f[TN] - f[BS])) + ((f[NW] - f[SE]) + (f[NE] - f[SW]))) + (f[N] - f[S])) ;

                        vx3 =
                            ((((f[TNE] - f[BSW]) + (f[TSE] - f[BNW])) + ((f[TNW] - f[BSE]) + (f[TSW] - f[BNE]))) +
                            (((f[TS] - f[BN]) + (f[TN] - f[BS])) + ((f[TW] - f[BE]) + (f[TE] - f[BW]))) + (f[T] - f[B]));

                    }
                    else {
                        vx1 =
                            ((((f[TNE] - f[BSW]) + (f[TSE] - f[BNW])) + ((f[BSE] - f[TNW]) + (f[BNE] - f[TSW]))) +
                            (((f[BE] - f[TW]) + (f[TE] - f[BW])) + ((f[SE] - f[NW]) + (f[NE] - f[SW]))) + (f[E] - f[W])) /
                                (rho * c1o3) +
                            mu * dX1_phi / (2 * rho);

                        vx2 =
                            ((((f[TNE] - f[BSW]) + (f[BNW] - f[TSE])) + ((f[TNW] - f[BSE]) + (f[BNE] - f[TSW]))) +
                            (((f[BN] - f[TS]) + (f[TN] - f[BS])) + ((f[NW] - f[SE]) + (f[NE] - f[SW]))) + (f[N] - f[S])) /
                                (rho * c1o3) +
                            mu * dX2_phi / (2 * rho);

                        vx3 =
                            ((((f[TNE] - f[BSW]) + (f[TSE] - f[BNW])) + ((f[TNW] - f[BSE]) + (f[TSW] - f[BNE]))) +
                            (((f[TS] - f[BN]) + (f[TN] - f[BS])) + ((f[TW] - f[BE]) + (f[TE] - f[BW]))) + (f[T] - f[B])) /
                                (rho * c1o3) +
                            mu * dX3_phi / (2 * rho);

                    }

                    p1 = (((f[TNE] + f[BSW]) + (f[TSE] + f[BNW])) + ((f[BSE] + f[TNW]) + (f[TSW] + f[BNE])) +
                          (((f[NE] + f[SW]) + (f[SE] + f[NW])) + ((f[TE] + f[BW]) + (f[BE] + f[TW])) +
                           ((f[BN] + f[TS]) + (f[TN] + f[BS]))) +
                          ((f[E] + f[W]) + (f[N] + f[S]) + (f[T] + f[B])) + f[REST]) +
                         (vx1 * rhoToPhi * dX1_phi * c1o3 + vx2 * rhoToPhi * dX2_phi * c1o3 +
                          vx3 * rhoToPhi * dX3_phi * c1o3) /
                             2.0;
                    p1 = rho * p1 * c1o3;

                    // calcMacros(f,rho,vx1,vx2,vx3);
                    // tempfield(ix1,ix2,ix3)= ; Added by HS
                    // double press = D3Q27System::calcPress(f,rho,vx1,vx2,vx3);

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

                    if (UbMath::isNaN(phi[REST]) || UbMath::isInfinity(phi[REST]))
                        UB_THROW(UbException(
                            UB_EXARGS, "phi is not a number (nan or -1.#IND) or infinity number -1.#INF in block=" +
                                           block->toString() + ", node=" + UbSystem::toString(ix1) + "," +
                                           UbSystem::toString(ix2) + "," + UbSystem::toString(ix3)));

                    if (UbMath::isNaN(p1) || UbMath::isInfinity(p1))
                        UB_THROW( UbException(UB_EXARGS,"p1 is not a number (nan or -1.#IND) or infinity number -1.#INF in block="+block->toString()+
                        ", node="+UbSystem::toString(ix1)+","+UbSystem::toString(ix2)+","+UbSystem::toString(ix3)));

                    data[index++].push_back(phi[REST]);
                    data[index++].push_back(vx1);
                    data[index++].push_back(vx2);
                    data[index++].push_back(vx3);
                    data[index++].push_back(p1);
                    data[index++].push_back(phi2[REST]);
                    if (dynamicPointerCast<MultiphaseTwoPhaseFieldsVelocityCumulantLBMKernel>(kernel)->pressure) {
                        data[index++].push_back((*dynamicPointerCast<MultiphaseTwoPhaseFieldsVelocityCumulantLBMKernel>(kernel)->pressure)(ix1, ix2, ix3));
                    }
                   // else { data[index++].push_back(999); }
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
}

LBMReal WriteMultiphaseQuantitiesCoProcessor::gradX1_phi(const LBMReal *const &h)
{
    using namespace D3Q27System;
    LBMReal sum = 0.0;
    for (int k = FSTARTDIR; k <= FENDDIR; k++) {
        sum += WEIGTH[k] * DX1[k] * h[k];
    }
    return 3.0 * sum;
}
LBMReal WriteMultiphaseQuantitiesCoProcessor::gradX2_phi(const LBMReal *const &h)
{
    using namespace D3Q27System;
    LBMReal sum = 0.0;
    for (int k = FSTARTDIR; k <= FENDDIR; k++) {
        sum += WEIGTH[k] * DX2[k] * h[k];
    }
    return 3.0 * sum;
}

LBMReal WriteMultiphaseQuantitiesCoProcessor::gradX3_phi(const LBMReal *const &h)
{
    using namespace D3Q27System;
    LBMReal sum = 0.0;
    for (int k = FSTARTDIR; k <= FENDDIR; k++) {
        sum += WEIGTH[k] * DX3[k] * h[k];
    }
    return 3.0 * sum;
}

LBMReal WriteMultiphaseQuantitiesCoProcessor::nabla2_phi(const LBMReal *const &h)
{
    using namespace D3Q27System;
    LBMReal sum = 0.0;
    for (int k = FSTARTDIR; k <= FENDDIR; k++) {
        sum += WEIGTH[k] * (h[k] - h[REST]);
    }
    return 6.0 * sum;
}