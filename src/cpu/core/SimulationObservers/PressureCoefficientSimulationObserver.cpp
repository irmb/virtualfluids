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
//! \file PressureCoefficientSimulationObserver.cpp
//! \ingroup SimulationObservers
//! \author Konstantin Kutscher
//=======================================================================================
#include "PressureCoefficientSimulationObserver.h"
#include <WbWriterVtkXmlASCII.h>

#include "BCArray3D.h"
#include "BCSet.h"
#include "Block3D.h"
#include <parallel/Communicator.h>
#include "D3Q27Interactor.h"
#include "DataSet3D.h"
#include "GbCuboid3D.h"
#include "Grid3D.h"
#include "LBMKernel.h"
#include "UbScheduler.h"

PressureCoefficientSimulationObserver::PressureCoefficientSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s,
                                                               GbCuboid3DPtr plane, const std::string &path,
                                                               std::shared_ptr<vf::parallel::Communicator> comm)
    : SimulationObserver(grid, s), plane(plane), path(path), comm(comm)
{
    maxStep       = scheduler->getMaxEnd();
    numberOfSteps = int(maxStep - scheduler->getMinBegin());
}
//////////////////////////////////////////////////////////////////////////
PressureCoefficientSimulationObserver::~PressureCoefficientSimulationObserver() = default;
//////////////////////////////////////////////////////////////////////////
void PressureCoefficientSimulationObserver::update(real step)
{
    if (scheduler->isDue(step))
        collectData(step);

    UBLOG(logDEBUG3, "D3Q27ForcesSimulationObserver::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void PressureCoefficientSimulationObserver::collectData(real step)
{
    calculateRho();

    if (step == maxStep) {
        writeValues((int)step);
    }
}
//////////////////////////////////////////////////////////////////////////
void PressureCoefficientSimulationObserver::calculateRho()
{
    real f[D3Q27System::ENDF + 1];
    real vx1, vx2, vx3, rho;
    std::vector<real> values;
    std::vector<real> rvalues;

    for (SPtr<D3Q27Interactor> interactor : interactors) {
        typedef std::map<SPtr<Block3D>, std::set<std::vector<int>>> TransNodeIndicesMap;
        for (TransNodeIndicesMap::value_type t : interactor->getBcNodeIndicesMap()) {
            SPtr<Block3D> block                          = t.first;
            std::set<std::vector<int>> &bcNodeIndicesSet = t.second;

            SPtr<ILBMKernel> kernel                 = block->getKernel();
            SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
            SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();

            UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);
            //         UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
            UbTupleDouble3 nodeOffset = grid->getNodeOffset(block);
            real dx                 = grid->getDeltaX(block);

            if (kernel->getCompressible()) {
                calcMacros = &D3Q27System::calcCompMacroscopicValues;
            } else {
                calcMacros = &D3Q27System::calcIncompMacroscopicValues;
            }

            int ghostLayerWidth = kernel->getGhostLayerWidth();
            int minX1           = ghostLayerWidth;
            int maxX1           = (int)bcArray->getNX1() - 1 - ghostLayerWidth;
            int minX2           = ghostLayerWidth;
            int maxX2           = (int)bcArray->getNX2() - 1 - ghostLayerWidth;
            int minX3           = ghostLayerWidth;
            int maxX3           = (int)bcArray->getNX3() - 1 - ghostLayerWidth;

            for (std::vector<int> node : bcNodeIndicesSet) {
                int x1 = node[0];
                int x2 = node[1];
                int x3 = node[2];

                // without ghost nodes
                if (x1 < minX1 || x1 > maxX1 || x2 < minX2 || x2 > maxX2 || x3 < minX3 || x3 > maxX3)
                    continue;

                if (bcArray->isFluid(
                        x1, x2,
                        x3)) // es kann sein, dass der node von einem anderen interactor z.B. als solid gemarkt wurde!!!
                {
                    real cx1 = val<1>(org) - val<1>(nodeOffset) + x1 * dx;
                    real cx2 = val<2>(org) - val<2>(nodeOffset) + x2 * dx;
                    real cx3 = val<3>(org) - val<3>(nodeOffset) + x3 * dx;
                    if (plane->isPointInGbObject3D(cx1, cx2, cx3)) {
                        distributions->getPreCollisionDistribution(f, x1, x2, x3);
                        calcMacros(f, rho, vx1, vx2, vx3);
                        values.push_back(cx1);
                        values.push_back(cx2);
                        values.push_back(cx3);
                        values.push_back(rho);
                    }
                }
            }
        }
    }

    comm->allGather(values, rvalues);
    if (comm->getProcessID() == comm->getRoot()) {
        if (outValues.size() == 0) {
            outValues.resize(rvalues.size());
        }
        int size = (int)rvalues.size();
        for (int i = 0; i < size; i += 4) {
            outValues[i]     = rvalues[i];
            outValues[i + 1] = rvalues[i + 1];
            outValues[i + 2] = rvalues[i + 2];
            outValues[i + 3] += rvalues[i + 3];
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void PressureCoefficientSimulationObserver::writeValues(int step)
{
    if (comm->getProcessID() == comm->getRoot()) {
        datanames.resize(0);
        datanames.emplace_back("rho");
        data.resize(datanames.size());

        std::ofstream ostr;
        std::string fname = path + UbSystem::toString(step) + ".csv";
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

        ostr << "x1;x2;x3;rho\n";

        int size = (int)outValues.size();
        for (int i = 0; i < size; i += 4) {
            ostr << outValues[i] << ";" << outValues[i + 1] << ";" << outValues[i + 2] << ";"
                 << outValues[i + 3] / numberOfSteps << std::endl;

            int index = 0;
            nodes.push_back(makeUbTuple(float(outValues[i]), float(outValues[i + 1]), float(outValues[i + 2])));

            data[index++].push_back(outValues[i + 3] / numberOfSteps);
        }

        ostr.close();

        WbWriterVtkXmlASCII::getInstance()->writeNodesWithNodeData(path + UbSystem::toString(step), nodes, datanames,
                                                                   data);

        fname = path + UbSystem::toString(step) + ".bin";
        std::ofstream out(fname.c_str(), std::ios::out | std::ios::binary);
        if (!out) {
            out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
            std::string path = UbSystem::getPathFromString(fname);
            if (path.size() > 0) {
                UbSystem::makeDirectory(path);
                out.open(fname.c_str(), std::ios::out | std::ios::binary);
            }
            if (!out)
                throw UbException(UB_EXARGS, "couldn't open file " + fname);
        }

        out.write((char *)&outValues[0], outValues.size() * sizeof(real));

        out.close();

        UBLOG(logINFO, "PressureCoefficientSimulationObserver::writeValues() step: " << (int)step);
    }
}
void PressureCoefficientSimulationObserver::readValues(int step)
{
    if (comm->isRoot()) {
        std::string fname = path + UbSystem::toString(step) + ".bin";
        std::ifstream in(fname.c_str(), std::ios::in | std::ios::binary);
        if (!in) {
            throw UbException(UB_EXARGS, "couldn't open file " + fname);
        }

        // get length of file:
        in.seekg(0, in.end);
        int length = (int)in.tellg();
        in.seekg(0, in.beg);

        outValues.resize(length / sizeof(real));

        in.read((char *)&outValues[0], length);

        in.close();

        UBLOG(logINFO, "PressureCoefficientSimulationObserver::readValues() step: " << (int)step);
    }
}
//////////////////////////////////////////////////////////////////////////
void PressureCoefficientSimulationObserver::addInteractor(SPtr<D3Q27Interactor> interactor)
{
    interactors.push_back(interactor);
}
