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
//! \file ShearStressSimulationObserver.cpp
//! \ingroup SimulationObservers
//! \author Konstantin Kutscher, S. Uphoff, M. Geier, E. Goraki Fard
//=======================================================================================
#include "ShearStressSimulationObserver.h"
#include "BCSet.h"
#include "WbWriterVtkXmlASCII.h"

#include "BCArray3D.h"
#include "Block3D.h"
#include <parallel/Communicator.h>
#include "D3Q27Interactor.h"
#include "DataSet3D.h"
#include "Grid3D.h"
#include "Interpolator.h"
#include "LBMKernel.h"
#include "UbScheduler.h"

ShearStressSimulationObserver::ShearStressSimulationObserver(SPtr<Grid3D> grid, const std::string &path, WbWriter *const writer,
                                               SPtr<UbScheduler> s, SPtr<UbScheduler> rs)
    : SimulationObserver(grid, s), Resetscheduler(rs), path(path), writer(writer)
{
    std::shared_ptr<vf::parallel::Communicator> comm = vf::parallel::Communicator::getInstance();
    normals.push_back(0);
    normals.push_back(0);
    normals.push_back(1);
    gridRank     = grid->getRank();
    minInitLevel = this->grid->getCoarsestInitializedLevel();
    maxInitLevel = this->grid->getFinestInitializedLevel();

    blockVector.resize(maxInitLevel + 1);
    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        grid->getBlocks(level, gridRank, true, blockVector[level]);
        for (SPtr<Block3D> block : blockVector[level]) {
            UbTupleInt3 nx                                   = grid->getBlockNX();
            SPtr<ShearStressValuesArray3D> shearStressValues = SPtr<ShearStressValuesArray3D>(
                new ShearStressValuesArray3D(14, val<1>(nx) + 1, val<2>(nx) + 1, val<3>(nx) + 1, vf::basics::constant::c0o1));
            block->getKernel()->getDataSet()->setShearStressValues(shearStressValues);
        }
    }
}
//////////////////////////////////////////////////////////////////////////
ShearStressSimulationObserver::~ShearStressSimulationObserver() = default;
//////////////////////////////////////////////////////////////////////////
void ShearStressSimulationObserver::update(real step)
{
    if (step == 0) {
        initDistance();
    }
    calculateShearStress(step);
    if (scheduler->isDue(step))
        collectData(step);
    UBLOG(logDEBUG3, "D3Q27ShearStressSimulationObserver::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void ShearStressSimulationObserver::collectData(real step)
{
    using namespace std;

    int istep = int(step);
    addData();

    // string partName = writer->writeNodesWithNodeData(path+ UbSystem::toString(gridRank)+ "_" +
    // UbSystem::toString(istep),nodes,datanames,data); size_t found=partName.find_last_of("//"); string piece =
    // partName.substr(found+1);

    // vector<string> cellDataNames;

    // std::shared_ptr<vf::parallel::Communicator> comm = vf::parallel::Communicator::getInstance();
    // vector<string> pieces = comm->gatherStrings(piece);
    // if (comm->getProcessID() == comm->getRoot())
    //{
    //   string pname =
    //   WbWriterVtkXmlASCII::getInstance()->writeParallelFile(path+"_"+UbSystem::toString(istep),pieces,datanames,cellDataNames);

    //   vector<string> filenames;
    //   filenames.push_back(pname);
    //   if (step == SimulationObserver::scheduler->getMinBegin())
    //   {
    //      WbWriterVtkXmlASCII::getInstance()->writeCollection(path+"__Shear_collection",filenames,istep,false);
    //   }
    //   else
    //   {
    //      WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(path+"__Shear_collection",filenames,istep,false);
    //   }
    //   UBLOG(logINFO,"D3Q27ShearStressSimulationObserver step: " << istep);
    //}

    string pfilePath, partPath, subfolder, cfilePath;
    subfolder = "shs" + UbSystem::toString(istep);
    pfilePath = path + "/shs/" + subfolder;
    cfilePath = path + "/shs/shs_collection";
    partPath  = pfilePath + "/shs" + UbSystem::toString(gridRank) + "_" + UbSystem::toString(istep);

    string partName = writer->writeNodesWithNodeData(partPath, nodes, datanames, data);
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
        UBLOG(logINFO, "D3Q27ShearStressSimulationObserver step: " << istep);
    }

    clearData();
}
//////////////////////////////////////////////////////////////////////////
void ShearStressSimulationObserver::clearData()
{
    nodes.clear();
    datanames.clear();
    data.clear();
}
//////////////////////////////////////////////////////////////////////////
void ShearStressSimulationObserver::calculateShearStress(real timeStep)
{
    using namespace vf::lbm::dir;
    using namespace D3Q27System;
    using namespace vf::basics::constant;

    real f[27];
    real vx, vy, vz, sxx, syy, szz, sxy, syz, sxz;

    for (SPtr<D3Q27Interactor> interactor : interactors) {
        typedef std::map<SPtr<Block3D>, std::set<std::vector<int>>> TransNodeIndicesMap;
        for (TransNodeIndicesMap::value_type t : interactor->getBcNodeIndicesMap()) {
            SPtr<Block3D> block                             = t.first;
            std::set<std::vector<int>> &transNodeIndicesSet = t.second;

            SPtr<ILBMKernel> kernel                 = block->getKernel();
            SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
            SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
            SPtr<ShearStressValuesArray3D> ssv      = kernel->getDataSet()->getShearStressValues();

            int ghostLayer     = kernel->getGhostLayerWidth();
            real collFactor = kernel->getCollisionFactor();

            int minX1 = ghostLayer;
            int maxX1 = (int)bcArray->getNX1() - 1 - ghostLayer;
            int minX2 = ghostLayer;
            int maxX2 = (int)bcArray->getNX2() - 1 - ghostLayer;
            int minX3 = ghostLayer;
            int maxX3 = (int)bcArray->getNX3() - 1 - ghostLayer;

            for (std::vector<int> node : transNodeIndicesSet) {
                int ix1 = node[0];
                int ix2 = node[1];
                int ix3 = node[2];

                // without ghost nodes
                if (ix1 < minX1 || ix1 > maxX1 || ix2 < minX2 || ix2 > maxX2 || ix3 < minX3 || ix3 > maxX3)
                    continue;

                if (bcArray->isFluid(ix1, ix2, ix3)) {
                    real q        = (*ssv)(normalq, ix1, ix2, ix3);
                    real numPoint = (*ssv)(numberOfPoint, ix1, ix2, ix3);
                    if (q == 0 || numPoint != 3)
                        continue;
                    // if (q==0)continue;
                    //////////////////////////////////////////////////////////////////////////
                    // read distribution
                    ////////////////////////////////////////////////////////////////////////////
                    distributions->getPreCollisionDistribution(f, ix1, ix2, ix3);
                    //////////////////////////////////////////////////////////////////////////
                    // compute velocity
                    //////////////////////////////////////////////////////////////////////////
                    vx = ((((f[dPPP] - f[dMMM]) + (f[dPMP] - f[dMPM])) + ((f[dPMM] - f[dMPP]) + (f[dPPM] - f[dMMP]))) +
                          (((f[dP0M] - f[dM0P]) + (f[dP0P] - f[dM0M])) + ((f[dPM0] - f[dMP0]) + (f[dPP0] - f[dMM0]))) + (f[dP00] - f[dM00]));

                    vy = ((((f[dPPP] - f[dMMM]) + (f[dMPM] - f[dPMP])) + ((f[dMPP] - f[dPMM]) + (f[dPPM] - f[dMMP]))) +
                          (((f[d0PM] - f[d0MP]) + (f[d0PP] - f[d0MM])) + ((f[dMP0] - f[dPM0]) + (f[dPP0] - f[dMM0]))) + (f[d0P0] - f[d0M0]));

                    vz = ((((f[dPPP] - f[dMMM]) + (f[dPMP] - f[dMPM])) + ((f[dMPP] - f[dPMM]) + (f[dMMP] - f[dPPM]))) +
                          (((f[d0MP] - f[d0PM]) + (f[d0PP] - f[d0MM])) + ((f[dM0P] - f[dP0M]) + (f[dP0P] - f[dM0M]))) + (f[d00P] - f[d00M]));

                    sxy = c3o1 * collFactor / (collFactor - c1o1) *
                          (((f[dPPP] + f[dMMM]) - (f[dPMP] + f[dMPM])) + (-(f[dPMM] + f[dMPP]) + (f[dMMP] + f[dPPM])) +
                           (((f[dPP0] + f[dMM0]) - (f[dPM0] + f[dMP0]))) - vx * vy);

                    sxz = c3o1 * collFactor / (collFactor - c1o1) *
                          (((f[dPPP] + f[dMMM]) + (f[dPMP] + f[dMPM])) + (-(f[dPMM] + f[dMPP]) - (f[dMMP] + f[dPPM])) +
                           ((f[dP0P] + f[dM0M]) - (f[dP0M] + f[dM0P])) - vx * vz);

                    syz = c3o1 * collFactor / (collFactor - c1o1) *
                          (((f[dPPP] + f[dMMM]) - (f[dPMP] + f[dMPM])) + ((f[dPMM] + f[dMPP]) - (f[dMMP] + f[dPPM])) +
                           (-(f[d0PM] + f[d0MP]) + (f[d0PP] + f[d0MM])) - vy * vz);

                    real dxxMyy = c3o1 / c2o1 * collFactor / (collFactor - c1o1) *
                                     (((f[dP0P] + f[dM0M]) + (f[dP0M] + f[dM0P])) - ((f[d0PM] + f[d0MP]) + (f[d0PP] + f[d0MM])) +
                                      ((f[dP00] + f[dM00]) - (f[d0P0] + f[d0M0])) - vx * vx + vy * vy);

                    real dxxMzz = c3o1 / c2o1 * collFactor / (collFactor - c1o1) *
                                     ((((f[dPP0] + f[dMM0]) + (f[dPM0] + f[dMP0])) - ((f[d0PM] + f[d0MP]) + (f[d0PP] + f[d0MM]))) +
                                      ((f[dP00] + f[dM00]) - (f[d00P] + f[d00M])) - vx * vx + vz * vz);

                    // real dyyMzz =3.0/2.0 *collFactor/(collFactor-1.0)*((((f[NE] + f[SW]) + (f[SE] +
                    // f[NW]))-((f[TE] + f[BW])+(f[BE]+ f[TW])))
                    //    +((f[N] + f[S])-(f[T] + f[B])) -vy*vy +vz*vz);

                    sxx = (dxxMyy + dxxMzz) / c3o1; // weil dxxPyyPzz=0

                    syy = (dxxMzz - c2o1 * dxxMyy) / c3o1;

                    szz = (dxxMyy - c2o1 * dxxMzz) / c3o1;

                    //////////////////////////////////////////////////////////////////////////
                    // compute average values
                    //////////////////////////////////////////////////////////////////////////
                    (*ssv)(AvVx, ix1, ix2, ix3) = ((*ssv)(AvVx, ix1, ix2, ix3) * timeStep + vx) / (timeStep + c1o1);
                    (*ssv)(AvVy, ix1, ix2, ix3) = ((*ssv)(AvVy, ix1, ix2, ix3) * timeStep + vy) / (timeStep + c1o1);
                    (*ssv)(AvVz, ix1, ix2, ix3) = ((*ssv)(AvVz, ix1, ix2, ix3) * timeStep + vz) / (timeStep + c1o1);

                    (*ssv)(AvSxx, ix1, ix2, ix3) = ((*ssv)(AvSxx, ix1, ix2, ix3) * timeStep + sxx) / (timeStep + c1o1);
                    (*ssv)(AvSyy, ix1, ix2, ix3) = ((*ssv)(AvSyy, ix1, ix2, ix3) * timeStep + syy) / (timeStep + c1o1);
                    (*ssv)(AvSzz, ix1, ix2, ix3) = ((*ssv)(AvSzz, ix1, ix2, ix3) * timeStep + szz) / (timeStep + c1o1);
                    (*ssv)(AvSxy, ix1, ix2, ix3) = ((*ssv)(AvSxy, ix1, ix2, ix3) * timeStep + sxy) / (timeStep + c1o1);
                    (*ssv)(AvSyz, ix1, ix2, ix3) = ((*ssv)(AvSyz, ix1, ix2, ix3) * timeStep + syz) / (timeStep + c1o1);
                    (*ssv)(AvSxz, ix1, ix2, ix3) = ((*ssv)(AvSxz, ix1, ix2, ix3) * timeStep + sxz) / (timeStep + c1o1);
                }
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void ShearStressSimulationObserver::addData()
{
    using namespace vf::basics::constant;

    // Diese Daten werden geschrieben:
    datanames.resize(0);
    datanames.emplace_back("y^plus");
    datanames.emplace_back("u_tau");
    // datanames.push_back("yPlusFD");

    data.resize(datanames.size());

    for (const auto &interactor : interactors) {
        using TransNodeIndicesMap = std::map<SPtr<Block3D>, std::set<std::vector<int>>>;
        for (TransNodeIndicesMap::value_type t : interactor->getBcNodeIndicesMap()) {
            SPtr<Block3D> block                             = t.first;
            std::set<std::vector<int>> &transNodeIndicesSet = t.second;

            UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);
            //         UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
            UbTupleDouble3 nodeOffset = grid->getNodeOffset(block);
            real dx                 = grid->getDeltaX(block);

            SPtr<ILBMKernel> kernel                 = block->getKernel();
            SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
            SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
            SPtr<ShearStressValuesArray3D> ssv      = kernel->getDataSet()->getShearStressValues();

            int ghostLayer     = kernel->getGhostLayerWidth();
            real collFactor = kernel->getCollisionFactor();

            int minX1 = ghostLayer;
            int maxX1 = (int)bcArray->getNX1() - 1 - ghostLayer;
            int minX2 = ghostLayer;
            int maxX2 = (int)bcArray->getNX2() - 1 - ghostLayer;
            int minX3 = ghostLayer;
            int maxX3 = (int)bcArray->getNX3() - 1 - ghostLayer;

            //         int level=block->getLevel();
            //         if(level==1)
            //         {
            //            int le=0;
            //         }
            for (std::vector<int> node : transNodeIndicesSet) {
                int ix1 = node[0];
                int ix2 = node[1];
                int ix3 = node[2];

                // without ghost nodes
                if (ix1 < minX1 || ix1 > maxX1 || ix2 < minX2 || ix2 > maxX2 || ix3 < minX3 || ix3 > maxX3)
                    continue;

                if (bcArray->isFluid(ix1, ix2, ix3)) {
                    real q        = (*ssv)(normalq, ix1, ix2, ix3);
                    real numPoint = (*ssv)(numberOfPoint, ix1, ix2, ix3);
                    if (q == 0 || numPoint != 3)
                        continue;
                    // if (q==0)continue;

                    int index = 0;
                    nodes.push_back(makeUbTuple(float(val<1>(org) - val<1>(nodeOffset) + ix1 * dx),
                                                float(val<2>(org) - val<2>(nodeOffset) + ix2 * dx),
                                                float(val<3>(org) - val<3>(nodeOffset) + ix3 * dx)));

                    //////get normal and distance//////
                    real A, B, C;
                    A = (*ssv)(normalX1, ix1, ix2, ix3);
                    B = (*ssv)(normalX2, ix1, ix2, ix3);
                    C = (*ssv)(normalX3, ix1, ix2, ix3);

                    ///////////
                    // compute y plus
                    // double vtxSonja, vtySonja, vtzSonja; //tangent velocity
                    // double temp = (*av)(ix1,ix2,ix3,AvVx)*A+(*av)(ix1,ix2,ix3,AvVy)*B+(*av)(ix1,ix2,ix3,AvVz)*C;
                    // vtxSonja = (*av)(ix1,ix2,ix3,AvVx)-normals[0]*temp;
                    // vtySonja = (*av)(ix1,ix2,ix3,AvVy)-normals[1]*temp;
                    // vtzSonja = (*av)(ix1,ix2,ix3,AvVz)-normals[2]*temp;

                    real vtx = (B * B * (*ssv)(AvVx, ix1, ix2, ix3) + C * C * (*ssv)(AvVx, ix1, ix2, ix3) -
                                  A * B * (*ssv)(AvVy, ix1, ix2, ix3) - A * C * (*ssv)(AvVy, ix1, ix2, ix3)) /
                                 (A * A + B * B + C * C);
                    real vty = (-(A * B * (*ssv)(AvVx, ix1, ix2, ix3)) + A * A * (*ssv)(AvVy, ix1, ix2, ix3) +
                                  C * C * (*ssv)(AvVy, ix1, ix2, ix3) - B * C * (*ssv)(AvVz, ix1, ix2, ix3)) /
                                 (A * A + B * B + C * C);
                    real vtz = (-(A * C * (*ssv)(AvVx, ix1, ix2, ix3)) - B * C * (*ssv)(AvVy, ix1, ix2, ix3) +
                                  A * A * (*ssv)(AvVz, ix1, ix2, ix3) + B * B * (*ssv)(AvVz, ix1, ix2, ix3)) /
                                 (A * A + B * B + C * C);

                    real normVt = sqrt(vtx * vtx + vty * vty + vtz * vtz) + 1e-100;
                    real nvtx   = vtx / normVt;
                    real nvty   = vty / normVt;
                    real nvtz   = vtz / normVt;

                    real sx   = c1o2 * ((*ssv)(AvSxx, ix1, ix2, ix3) * nvtx + (*ssv)(AvSxy, ix1, ix2, ix3) * nvty +
                                       (*ssv)(AvSxz, ix1, ix2, ix3) * nvtz);
                    real sy   = c1o2 * ((*ssv)(AvSxy, ix1, ix2, ix3) * nvtx + (*ssv)(AvSyy, ix1, ix2, ix3) * nvty +
                                       (*ssv)(AvSyz, ix1, ix2, ix3) * nvtz);
                    real sz   = c1o2 * ((*ssv)(AvSxz, ix1, ix2, ix3) * nvtx + (*ssv)(AvSyz, ix1, ix2, ix3) * nvty +
                                       (*ssv)(AvSzz, ix1, ix2, ix3) * nvtz);
                    real sabs = sqrt(sx * sx + sy * sy + sz * sz);

                    real viscosity = (c1o1 / c3o1) * (c1o1 / collFactor - c1o2);
                    real rho       = c1o1;
                    real utau      = sqrt(viscosity / rho * sabs);

                    // double q=(*av)(ix1,ix2,ix3,normalq) ;
                    real yPlus = (utau * q) / viscosity;

                    data[index++].push_back(yPlus);
                    data[index++].push_back(utau);
                }
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void ShearStressSimulationObserver::reset(real step)
{
    if (Resetscheduler->isDue(step))
        resetData(step);

    UBLOG(logDEBUG3, "resetSimulationObserver::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void ShearStressSimulationObserver::resetData(real /*step*/)
{
    using namespace vf::basics::constant;

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (const auto &block : blockVector[level]) {
            if (block) {
                //            UbTupleDouble3 org          = grid->getBlockWorldCoordinates(block);
                //            UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
                //            UbTupleDouble3 nodeOffset   = grid->getNodeOffset(block);
                //            double         dx           = grid->getDeltaX(block);

                SPtr<ILBMKernel> kernel                 = block->getKernel();
                SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
                SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
                SPtr<ShearStressValuesArray3D> ssv      = kernel->getDataSet()->getShearStressValues();

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
                                (*ssv)(AvVx, ix1, ix2, ix3) = c0o1;
                                (*ssv)(AvVy, ix1, ix2, ix3) = c0o1;
                                (*ssv)(AvVz, ix1, ix2, ix3) = c0o1;

                                (*ssv)(AvSxx, ix1, ix2, ix3) = c0o1;
                                (*ssv)(AvSyy, ix1, ix2, ix3) = c0o1;
                                (*ssv)(AvSzz, ix1, ix2, ix3) = c0o1;
                                (*ssv)(AvSxy, ix1, ix2, ix3) = c0o1;
                                (*ssv)(AvSyz, ix1, ix2, ix3) = c0o1;
                                (*ssv)(AvSxz, ix1, ix2, ix3) = c0o1;
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
void ShearStressSimulationObserver::addInteractor(SPtr<D3Q27Interactor> interactor) { interactors.push_back(interactor); }
//////////////////////////////////////////////////////////////////////////
void ShearStressSimulationObserver::findPlane(int ix1, int ix2, int ix3, SPtr<Grid3D> grid, SPtr<Block3D> block, real &A,
                                       real &B, real &C, real &D, real &ii)
{
    using namespace vf::lbm::dir;
    using namespace vf::basics::constant;

    real x1plane = c0o1, y1plane = c0o1, z1plane = c0o1;
    real x2plane = c0o1, y2plane = c0o1, z2plane = c0o1;
    real x3plane = c0o1, y3plane = c0o1, z3plane = c0o1;
    SPtr<BoundaryConditions> bcPtr;
    real dx                               = grid->getDeltaX(block);
    SPtr<ILBMKernel> kernel                 = block->getKernel();
    SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
    SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
    bcPtr                                   = bcArray->getBC(ix1, ix2, ix3);
    int x, y, z;

    if (Interpolator::iCellHasSolid(bcArray, ix1, ix2, ix3)) {
        x = ix1;
        y = ix2;
        z = ix3;
    } else if (Interpolator::iCellHasSolid(bcArray, ix1, ix2 - 1, ix3)) {
        x = ix1 + 0;
        y = ix2 - 1;
        z = ix3 + 0;
    } // S
    else if (Interpolator::iCellHasSolid(bcArray, ix1, ix2, ix3 - 1)) {
        x = ix1 + 0;
        y = ix2 + 0;
        z = ix3 - 1;
    } // B
    else if (Interpolator::iCellHasSolid(bcArray, ix1 - 1, ix2, ix3)) {
        x = ix1 - 1;
        y = ix2 + 0;
        z = ix3 + 0;
    } // w
    else if (Interpolator::iCellHasSolid(bcArray, ix1, ix2 - 1, ix3 - 1)) {
        x = ix1 + 0;
        y = ix2 - 1;
        z = ix3 - 1;
    } // BS
    else if (Interpolator::iCellHasSolid(bcArray, ix1 - 1, ix2, ix3 - 1)) {
        x = ix1 - 1;
        y = ix2 + 0;
        z = ix3 - 1;
    } // BW
    else if (Interpolator::iCellHasSolid(bcArray, ix1 - 1, ix2 - 1, ix3)) {
        x = ix1 - 1;
        y = ix2 - 1;
        z = ix3 + 0;
    } // SW

    else if (Interpolator::iCellHasSolid(bcArray, ix1 - 1, ix2 - 1, ix3 - 1)) {
        x = ix1 - 1;
        y = ix2 - 1;
        z = ix3 - 1;
    } // BSW
    else if (Interpolator::iCellHasSolid(bcArray, ix1 + 1, ix2, ix3)) {
        x = ix1 + 1;
        y = ix2 + 0;
        z = ix3 + 0;
    } // E
    else if (Interpolator::iCellHasSolid(bcArray, ix1, ix2 + 1, ix3)) {
        x = ix1 + 0;
        y = ix2 + 1;
        z = ix3 + 0;
    } // N
    else if (Interpolator::iCellHasSolid(bcArray, ix1, ix2, ix3 + 1)) {
        x = ix1 + 0;
        y = ix2 + 0;
        z = ix3 + 1;
    } // T
    else if (Interpolator::iCellHasSolid(bcArray, ix1 + 1, ix2 + 1, ix3)) {
        x = ix1 + 1;
        y = ix2 + 1;
        z = ix3 + 0;
    } // NE
    else if (Interpolator::iCellHasSolid(bcArray, ix1 + 1, ix2, ix3 + 1)) {
        x = ix1 + 1;
        y = ix2 + 0;
        z = ix3 + 1;
    } // TE
    else if (Interpolator::iCellHasSolid(bcArray, ix1, ix2 + 1, ix3 + 1)) {
        x = ix1 + 0;
        y = ix2 + 1;
        z = ix3 + 1;
    } // TN
    else if (Interpolator::iCellHasSolid(bcArray, ix1 + 1, ix2 + 1, ix3 + 1)) {
        x = ix1 + 1;
        y = ix2 + 1;
        z = ix3 + 1;
    } // TNE

    else if (Interpolator::iCellHasSolid(bcArray, ix1 + 1, ix2 - 1, ix3)) {
        x = ix1 + 1;
        y = ix2 - 1;
        z = ix3 + 0;
    } // SE
    else if (Interpolator::iCellHasSolid(bcArray, ix1 - 1, ix2 + 1, ix3)) {
        x = ix1 - 1;
        y = ix2 + 1;
        z = ix3 + 0;
    } // NW
    else if (Interpolator::iCellHasSolid(bcArray, ix1 + 1, ix2, ix3 - 1)) {
        x = ix1 + 1;
        y = ix2 + 0;
        z = ix3 - 1;
    } // BE
    else if (Interpolator::iCellHasSolid(bcArray, ix1 - 1, ix2, ix3 + 1)) {
        x = ix1 - 1;
        y = ix2 + 0;
        z = ix3 + 1;
    } // TW
    else if (Interpolator::iCellHasSolid(bcArray, ix1 + 0, ix2 + 1, ix3 - 1)) {
        x = ix1 + 0;
        y = ix2 + 1;
        z = ix3 - 1;
    } // BN
    else if (Interpolator::iCellHasSolid(bcArray, ix1 + 0, ix2 - 1, ix3 + 1)) {
        x = ix1 + 0;
        y = ix2 - 1;
        z = ix3 + 1;
    } // TS

    else if (Interpolator::iCellHasSolid(bcArray, ix1 - 1, ix2 + 1, ix3 + 1)) {
        x = ix1 - 1;
        y = ix2 + 1;
        z = ix3 + 1;
    } // TNW
    else if (Interpolator::iCellHasSolid(bcArray, ix1 + 1, ix2 - 1, ix3 + 1)) {
        x = ix1 + 1;
        y = ix2 - 1;
        z = ix3 + 1;
    } // TSE
    else if (Interpolator::iCellHasSolid(bcArray, ix1 - 1, ix2 - 1, ix3 + 1)) {
        x = ix1 - 1;
        y = ix2 - 1;
        z = ix3 + 1;
    } // TSW
    else if (Interpolator::iCellHasSolid(bcArray, ix1 + 1, ix2 + 1, ix3 - 1)) {
        x = ix1 + 1;
        y = ix2 + 1;
        z = ix3 - 1;
    } // BNE
    else if (Interpolator::iCellHasSolid(bcArray, ix1 - 1, ix2 + 1, ix3 - 1)) {
        x = ix1 - 1;
        y = ix2 + 1;
        z = ix3 - 1;
    } // BNW
    else if (Interpolator::iCellHasSolid(bcArray, ix1 + 1, ix2 - 1, ix3 - 1)) {
        x = ix1 + 1;
        y = ix2 - 1;
        z = ix3 - 1;
    } // BSE

    else {
        {
            UB_THROW(UbException(UB_EXARGS, "there is no cell  ix1=" + UbSystem::toString(ix1) +
                                                "ix2=" + UbSystem::toString(ix2) + "ix3=" + UbSystem::toString(ix3) +
                                                "GlobalID=" + UbSystem::toString(block->getGlobalID()) +
                                                "dx=" + UbSystem::toString(dx) +
                                                "T=" + UbSystem::toString(bcPtr->getQ(d00P)) +
                                                "B=" + UbSystem::toString(bcPtr->getQ(d00M)) +
                                                "E=" + UbSystem::toString(bcPtr->getQ(dP00)) +
                                                "W=" + UbSystem::toString(bcPtr->getQ(dM00)) +
                                                "N=" + UbSystem::toString(bcPtr->getQ(d0P0)) +
                                                "S=" + UbSystem::toString(bcPtr->getQ(d0M0)) +
                                                "NE=" + UbSystem::toString(bcPtr->getQ(dPP0)) +
                                                "SW=" + UbSystem::toString(bcPtr->getQ(dMM0)) +
                                                "SE=" + UbSystem::toString(bcPtr->getQ(dPM0)) +
                                                "NW=" + UbSystem::toString(bcPtr->getQ(dMP0)) +
                                                "TE=" + UbSystem::toString(bcPtr->getQ(dP0P)) +
                                                "BW=" + UbSystem::toString(bcPtr->getQ(dM0M)) +
                                                "BE=" + UbSystem::toString(bcPtr->getQ(dP0M)) +
                                                "TW=" + UbSystem::toString(bcPtr->getQ(dM0P)) +
                                                "TN=" + UbSystem::toString(bcPtr->getQ(d0PP)) +
                                                "BS=" + UbSystem::toString(bcPtr->getQ(d0MM)) +
                                                "BN=" + UbSystem::toString(bcPtr->getQ(d0PM)) +
                                                "TS=" + UbSystem::toString(bcPtr->getQ(d0MP)) +
                                                "TNE=" + UbSystem::toString(bcPtr->getQ(dPPP)) +
                                                "TNW=" + UbSystem::toString(bcPtr->getQ(dMPP)) +
                                                "TSE=" + UbSystem::toString(bcPtr->getQ(dPMP)) +
                                                "TSW=" + UbSystem::toString(bcPtr->getQ(dMMP)) +
                                                "BNE=" + UbSystem::toString(bcPtr->getQ(dPPM)) +
                                                "BNW=" + UbSystem::toString(bcPtr->getQ(dMPM)) +
                                                "BSE=" + UbSystem::toString(bcPtr->getQ(dPMM)) +
                                                "BSW=" + UbSystem::toString(bcPtr->getQ(dMMM) * dx)));
        }
    }

    if (Interpolator::iCellHasSolid(bcArray, x, y, z)) {
        for (int i = x; i <= x + 1; i++) {
            for (int j = y; j <= y + 1; j++) {
                for (int k = z; k <= z + 1; k++) {
                    Vector3D pointplane1 = grid->getNodeCoordinates(block, i, j, k);

                    real iph = pointplane1[0];
                    real jph = pointplane1[1];
                    real kph = pointplane1[2];

                    if (!bcArray->isSolid(i, j, k)) {
                        SPtr<BoundaryConditions> bcPtrIn = bcArray->getBC(i, j, k);
                        if (bcPtrIn) {
                            for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
                                if (ii <= 2) {
                                    real q = bcPtrIn->getQ(fdir);
                                    if (q != 999.00000) {
                                        if (fdir == dP00) {
                                            // if(!bcArray->isSolid(i, j, k))continue;
                                            if (i + q <= x + 1) {
                                                if (ii == 0) {
                                                    x1plane = iph + q * dx;
                                                    y1plane = jph;
                                                    z1plane = kph;
                                                    ii++;
                                                } else if (ii == 1) {
                                                    x2plane = iph + q * dx;
                                                    y2plane = jph;
                                                    z2plane = kph;
                                                    if (x1plane != x2plane || y1plane != y2plane || z1plane != z2plane)
                                                        ii++;
                                                } else if (ii == 2) {
                                                    x3plane = iph + q * dx;
                                                    y3plane = jph;
                                                    z3plane = kph;
                                                    if ((x3plane != x1plane || y3plane != y1plane ||
                                                         z3plane != z1plane) &&
                                                        (x2plane != x3plane || y2plane != y3plane ||
                                                         z2plane != z3plane))
                                                        ii++;
                                                }
                                            }
                                        }
                                        if (fdir == dM00) {
                                            // if(!bcArray->isSolid(i, j, k))continue;
                                            if (i - q >= x) {
                                                if (ii == 0) {
                                                    x1plane = iph - q * dx;
                                                    y1plane = jph;
                                                    z1plane = kph;
                                                    ii++;
                                                } else if (ii == 1) {
                                                    x2plane = iph - q * dx;
                                                    y2plane = jph;
                                                    z2plane = kph;
                                                    if (x1plane != x2plane || y1plane != y2plane || z1plane != z2plane)
                                                        ii++;
                                                } else if (ii == 2) {
                                                    x3plane = iph - q * dx;
                                                    y3plane = jph;
                                                    z3plane = kph;
                                                    if ((x3plane != x1plane || y3plane != y1plane ||
                                                         z3plane != z1plane) &&
                                                        (x2plane != x3plane || y2plane != y3plane ||
                                                         z2plane != z3plane))
                                                        ii++;
                                                }
                                            }
                                        }
                                        if (fdir == d0P0) {
                                            // if(!bcArray->isSolid(i, j, k))continue;
                                            if (j + q <= y + 1) {
                                                if (ii == 0) {
                                                    x1plane = iph;
                                                    y1plane = jph + q * dx;
                                                    z1plane = kph;
                                                    ii++;
                                                } else if (ii == 1) {
                                                    x2plane = iph;
                                                    y2plane = jph + q * dx;
                                                    z2plane = kph;
                                                    if (x1plane != x2plane || y1plane != y2plane || z1plane != z2plane)
                                                        ii++;
                                                } else if (ii == 2) {
                                                    x3plane = iph;
                                                    y3plane = jph + q * dx;
                                                    z3plane = kph;
                                                    if ((x3plane != x1plane || y3plane != y1plane ||
                                                         z3plane != z1plane) &&
                                                        (x2plane != x3plane || y2plane != y3plane ||
                                                         z2plane != z3plane))
                                                        ii++;
                                                }
                                            }
                                        }
                                        if (fdir == d0M0) {
                                            // if(!bcArray->isSolid(i, j, k))continue;
                                            if (j - q >= y) {
                                                if (ii == 0) {
                                                    x1plane = iph;
                                                    y1plane = jph - q * dx;
                                                    z1plane = kph;
                                                    ii++;
                                                } else if (ii == 1) {
                                                    x2plane = iph;
                                                    y2plane = jph - q * dx;
                                                    z2plane = kph;
                                                    if (x1plane != x2plane || y1plane != y2plane || z1plane != z2plane)
                                                        ii++;
                                                } else if (ii == 2) {
                                                    x3plane = iph;
                                                    y3plane = jph - q * dx;
                                                    z3plane = kph;
                                                    if ((x3plane != x1plane || y3plane != y1plane ||
                                                         z3plane != z1plane) &&
                                                        (x2plane != x3plane || y2plane != y3plane ||
                                                         z2plane != z3plane))
                                                        ii++;
                                                }
                                            }
                                        }

                                        if (fdir == d00P) {
                                            // if(!bcArray->isSolid(i, j, k))continue;
                                            if (k + q <= z + 1) {
                                                if (ii == 0) {
                                                    x1plane = iph;
                                                    y1plane = jph;
                                                    z1plane = kph + q * dx;
                                                    ii++;
                                                } else if (ii == 1) {
                                                    x2plane = iph;
                                                    y2plane = jph;
                                                    z2plane = kph + q * dx;
                                                    if (x1plane != x2plane || y1plane != y2plane || z1plane != z2plane)
                                                        ii++;
                                                } else if (ii == 2) {
                                                    x3plane = iph;
                                                    y3plane = jph;
                                                    z3plane = kph + q * dx;
                                                    if ((x3plane != x1plane || y3plane != y1plane ||
                                                         z3plane != z1plane) &&
                                                        (x2plane != x3plane || y2plane != y3plane ||
                                                         z2plane != z3plane))
                                                        ii++;
                                                }
                                            }
                                        }
                                        if (fdir == d00M) {
                                            // if(!bcArray->isSolid(i, j, k))continue;
                                            if (k - q >= z) {
                                                if (ii == 0) {
                                                    x1plane = iph;
                                                    y1plane = jph;
                                                    z1plane = kph - q * dx;
                                                    ii++;
                                                } else if (ii == 1) {
                                                    x2plane = iph;
                                                    y2plane = jph;
                                                    z2plane = kph - q * dx;
                                                    if (x1plane != x2plane || y1plane != y2plane || z1plane != z2plane)
                                                        ii++;
                                                } else if (ii == 2) {
                                                    x3plane = iph;
                                                    y3plane = jph;
                                                    z3plane = kph - q * dx;
                                                    if ((x3plane != x1plane || y3plane != y1plane ||
                                                         z3plane != z1plane) &&
                                                        (x2plane != x3plane || y2plane != y3plane ||
                                                         z2plane != z3plane))
                                                        ii++;
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
        }

        A = y1plane * (z2plane - z3plane) + y2plane * (z3plane - z1plane) + y3plane * (z1plane - z2plane);
        B = z1plane * (x2plane - x3plane) + z2plane * (x3plane - x1plane) + z3plane * (x1plane - x2plane);
        C = x1plane * (y2plane - y3plane) + x2plane * (y3plane - y1plane) + x3plane * (y1plane - y2plane);
        D = -(x1plane * (y2plane * z3plane - y3plane * z2plane) + x2plane * (y3plane * z1plane - y1plane * z3plane) +
              x3plane * (y1plane * z2plane - y2plane * z1plane));
    }
    if (ii != 3) {

        {
            {
                UB_THROW(UbException(
                    UB_EXARGS, "ii is=" + UbSystem::toString(ii) + "  ix1=" + UbSystem::toString(ix1) +
                                   " ix2=" + UbSystem::toString(ix2) + " ix3=" + UbSystem::toString(ix3) +
                                   " Block3D::GlobalID=" + UbSystem::toString(block->getGlobalID()) + " dx=" +
                                   UbSystem::toString(dx) + " T=" + UbSystem::toString(bcPtr->getQ(d00P)) +
                                   " B=" + UbSystem::toString(bcPtr->getQ(d00M)) +
                                   " E=" + UbSystem::toString(bcPtr->getQ(dP00)) +
                                   " W=" + UbSystem::toString(bcPtr->getQ(dM00)) +
                                   " N=" + UbSystem::toString(bcPtr->getQ(d0P0)) +
                                   " S=" + UbSystem::toString(bcPtr->getQ(d0M0)) +
                                   " NE=" + UbSystem::toString(bcPtr->getQ(dPP0)) +
                                   " SW=" + UbSystem::toString(bcPtr->getQ(dMM0)) +
                                   " SE=" + UbSystem::toString(bcPtr->getQ(dPM0)) +
                                   " NW=" + UbSystem::toString(bcPtr->getQ(dMP0)) +
                                   " TE=" + UbSystem::toString(bcPtr->getQ(dP0P)) +
                                   " BW=" + UbSystem::toString(bcPtr->getQ(dM0M)) +
                                   " BE=" + UbSystem::toString(bcPtr->getQ(dP0M)) +
                                   " TW=" + UbSystem::toString(bcPtr->getQ(dM0P)) +
                                   " TN=" + UbSystem::toString(bcPtr->getQ(d0PP)) +
                                   " BS=" + UbSystem::toString(bcPtr->getQ(d0MM)) +
                                   " BN=" + UbSystem::toString(bcPtr->getQ(d0PM)) +
                                   " TS=" + UbSystem::toString(bcPtr->getQ(d0MP)) +
                                   " TNE=" + UbSystem::toString(bcPtr->getQ(dPPP)) +
                                   " TNW=" + UbSystem::toString(bcPtr->getQ(dMPP)) +
                                   " TSE=" + UbSystem::toString(bcPtr->getQ(dPMP)) +
                                   " TSW=" + UbSystem::toString(bcPtr->getQ(dMMP)) +
                                   " BNE=" + UbSystem::toString(bcPtr->getQ(dPPM)) +
                                   " BNW=" + UbSystem::toString(bcPtr->getQ(dMPM)) +
                                   " BSE=" + UbSystem::toString(bcPtr->getQ(dPMM)) +
                                   " BSW=" + UbSystem::toString(bcPtr->getQ(dMMM))));
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
bool ShearStressSimulationObserver::checkUndefindedNodes(SPtr<BCArray3D> bcArray, int ix1, int ix2, int ix3)
{
    for (int i = ix1; i <= ix1 + 1; i++) {
        for (int j = ix2; j <= ix2 + 1; j++) {
            for (int k = ix3; k <= ix3 + 1; k++) {
                if (bcArray->isUndefined(i, j, k))
                    return true;
            }
        }
    }
    return false;
}
//////////////////////////////////////////////////////////////////////////////////////
void ShearStressSimulationObserver::initDistance()
{
    using namespace vf::lbm::dir;
    using namespace vf::basics::constant;

    for (const auto &interactor : interactors) {
        //      typedef std::map<SPtr<Block3D>, std::set< std::vector<int> > > TransNodeIndicesMap;
        for (const auto &t : interactor->getBcNodeIndicesMap()) {
            SPtr<Block3D> block                                   = t.first;
            const std::set<std::vector<int>> &transNodeIndicesSet = t.second;

            //         UbTupleDouble3 org          = grid->getBlockWorldCoordinates(block);
            //         UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
            //         UbTupleDouble3 nodeOffset   = grid->getNodeOffset(block);
            //         double         dx           = grid->getDeltaX(block);

            SPtr<ILBMKernel> kernel                 = block->getKernel();
            SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
            SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
            SPtr<ShearStressValuesArray3D> ssv      = kernel->getDataSet()->getShearStressValues();

            int ghostLayer = kernel->getGhostLayerWidth();
            //         real collFactor = kernel->getCollisionFactor();

            int minX1 = ghostLayer;
            int maxX1 = (int)bcArray->getNX1() - 1 - ghostLayer;
            int minX2 = ghostLayer;
            int maxX2 = (int)bcArray->getNX2() - 1 - ghostLayer;
            int minX3 = ghostLayer;
            int maxX3 = (int)bcArray->getNX3() - 1 - ghostLayer;

            for (std::vector<int> node : transNodeIndicesSet) {
                int ix1 = node[0];
                int ix2 = node[1];
                int ix3 = node[2];

                // without ghost nodes
                if (ix1 < minX1 || ix1 > maxX1 || ix2 < minX2 || ix2 > maxX2 || ix3 < minX3 || ix3 > maxX3)
                    continue;

                if (bcArray->isFluid(ix1, ix2, ix3)) {
                    auto bc = bcArray->getBC(ix1, ix2, ix3);
                    if (!bc)
                        continue;
                    if ((bc->hasDensityBoundary() || bc->hasVelocityBoundary()))
                        continue;
                    int numberOfCorner = 0;

                    if (bc->getQ(d00P) != 999.000) {
                        numberOfCorner++;
                    }
                    if (bc->getQ(d00M) != 999.000) {
                        numberOfCorner++;
                    }
                    if (bc->getQ(dP00) != 999.000) {
                        numberOfCorner++;
                    }
                    if (bc->getQ(dM00) != 999.000) {
                        numberOfCorner++;
                    }
                    if (bc->getQ(d0P0) != 999.000) {
                        numberOfCorner++;
                    }
                    if (bc->getQ(d0M0) != 999.000) {
                        numberOfCorner++;
                    }
                    // if(bc->hasVelocityBoundary()||bc->hasDensityBoundary())continue;
                    if (numberOfCorner > 1)
                        continue;
                    if (checkUndefindedNodes(bcArray, ix1, ix2, ix3))
                        continue;

                    //////get normal and distance//////
                    real A, B, C, D, ii = c0o1;
                    findPlane(ix1, ix2, ix3, grid, block, A, B, C, D, ii);
                    Vector3D pointplane1 = grid->getNodeCoordinates(block, ix1, ix2, ix3);
                    real ix1ph         = pointplane1[0];
                    real ix2ph         = pointplane1[1];
                    real ix3ph         = pointplane1[2];
                    real normalDis;
                    if (ii != 3) {
                        UB_THROW(UbException(UB_EXARGS, "not enough points to create plane" + UbSystem::toString(ii)));
                    } else {
                        real s = A * ix1ph + B * ix2ph + C * ix3ph +
                                   D; // The sign of s = Ax + By + Cz + D determines which side the point (x,y,z) lies
                                      // with respect to the plane. If s > 0 then the point lies on the same side as the
                                      // normal (A,B,C). If s < 0 then it lies on the opposite side, if s = 0 then the
                                      // point (x,y,z) lies on the plane.
                        if (s > 0) {
                            s = 1;
                        } else if (s < 0) {
                            s = -1;
                        } else {
                            s = 0;
                        }

                        normalDis = ((A * ix1ph + B * ix2ph + C * ix3ph + D) /
                                     sqrt(A * A + B * B + C * C)); /// distance point to plane xp-Xw=distance
                        normalDis *= s;

                        (*ssv)(normalX1, ix1, ix2, ix3)      = A;
                        (*ssv)(normalX2, ix1, ix2, ix3)      = B;
                        (*ssv)(normalX3, ix1, ix2, ix3)      = C;
                        (*ssv)(normalq, ix1, ix2, ix3)       = normalDis;
                        (*ssv)(numberOfPoint, ix1, ix2, ix3) = ii;
                    }
                }
            }
        }
    }
}
