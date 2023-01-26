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
//! \file D3Q27Interactor.cpp
//! \ingroup Interactor
//! \author Sören Freudiger
//! \author Sebastian Geller
//! \author Konstantin Kutscher
//=======================================================================================

#include "D3Q27Interactor.h"
#include <basics/utilities/UbLogger.h>
#include <basics/utilities/UbMath.h>

#include <basics/writer/WbWriterVtkXmlBinary.h>

#include "BCAdapter.h"
#include "BCArray3D.h"
#include "BCProcessor.h"
#include "Block3D.h"
#include "BoundaryConditions.h"
#include "Grid3D.h"
#include "LBMKernel.h"
#include <GbCuboid3D.h>
#include <GbLine3D.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
D3Q27Interactor::D3Q27Interactor() : Interactor3D()
{
    this->reinitWithStoredQsFlag = false;
    this->initRayVectors();
}
//////////////////////////////////////////////////////////////////////////
D3Q27Interactor::D3Q27Interactor(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, int type)
    : Interactor3D(geoObject3D, grid, type), relevantForForces(false)
{
    this->reinitWithStoredQsFlag = false;
    this->initRayVectors();
}
//////////////////////////////////////////////////////////////////////////
D3Q27Interactor::D3Q27Interactor(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, SPtr<BCAdapter> bcAdapter, int type)
    : Interactor3D(geoObject3D, grid, type), relevantForForces(false)
{
    this->reinitWithStoredQsFlag = false;
    this->addBCAdapter(bcAdapter);
    this->initRayVectors();
}
//////////////////////////////////////////////////////////////////////////
D3Q27Interactor::D3Q27Interactor(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, SPtr<BCAdapter> bcAdapter, int type,
                                 Interactor3D::Accuracy a)
    : Interactor3D(geoObject3D, grid, type, a), relevantForForces(false)
{
    this->reinitWithStoredQsFlag = false;
    this->addBCAdapter(bcAdapter);
    this->initRayVectors();
}
//////////////////////////////////////////////////////////////////////////
D3Q27Interactor::~D3Q27Interactor() = default;
//////////////////////////////////////////////////////////////////////////
void D3Q27Interactor::initRayVectors()
{
    using namespace vf::lbm::dir;

    int fdir;
    real c1oS2 = UbMath::one_over_sqrt2;
    real c1oS3 = UbMath::one_over_sqrt3;
    fdir         = DIR_P00;
    rayX1[fdir]  = 1.0;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = 0.0;
    fdir         = DIR_M00;
    rayX1[fdir]  = -1.0;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = 0.0;
    fdir         = DIR_0P0;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = 1.0;
    rayX3[fdir]  = 0.0;
    fdir         = DIR_0M0;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = -1.0;
    rayX3[fdir]  = 0.0;
    fdir         = DIR_00P;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = 1.0;
    fdir         = DIR_00M;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = -1.0;
    fdir         = DIR_PP0;
    rayX1[fdir]  = c1oS2;
    rayX2[fdir]  = c1oS2;
    rayX3[fdir]  = 0.0;
    fdir         = DIR_MM0;
    rayX1[fdir]  = -c1oS2;
    rayX2[fdir]  = -c1oS2;
    rayX3[fdir]  = 0.0;
    fdir         = DIR_PM0;
    rayX1[fdir]  = c1oS2;
    rayX2[fdir]  = -c1oS2;
    rayX3[fdir]  = 0.0;
    fdir         = DIR_MP0;
    rayX1[fdir]  = -c1oS2;
    rayX2[fdir]  = c1oS2;
    rayX3[fdir]  = 0.0;
    fdir         = DIR_P0P;
    rayX1[fdir]  = c1oS2;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = c1oS2;
    fdir         = DIR_M0M;
    rayX1[fdir]  = -c1oS2;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = -c1oS2;
    fdir         = DIR_P0M;
    rayX1[fdir]  = c1oS2;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = -c1oS2;
    fdir         = DIR_M0P;
    rayX1[fdir]  = -c1oS2;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = c1oS2;
    fdir         = DIR_0PP;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = c1oS2;
    rayX3[fdir]  = c1oS2;
    fdir         = DIR_0MM;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = -c1oS2;
    rayX3[fdir]  = -c1oS2;
    fdir         = DIR_0PM;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = c1oS2;
    rayX3[fdir]  = -c1oS2;
    fdir         = DIR_0MP;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = -c1oS2;
    rayX3[fdir]  = c1oS2;

    fdir        = DIR_MPP;
    rayX1[fdir] = -c1oS3;
    rayX2[fdir] = c1oS3;
    rayX3[fdir] = c1oS3;
    fdir        = DIR_PPP;
    rayX1[fdir] = c1oS3;
    rayX2[fdir] = c1oS3;
    rayX3[fdir] = c1oS3;
    fdir        = DIR_MMP;
    rayX1[fdir] = -c1oS3;
    rayX2[fdir] = -c1oS3;
    rayX3[fdir] = c1oS3;
    fdir        = DIR_PMP;
    rayX1[fdir] = c1oS3;
    rayX2[fdir] = -c1oS3;
    rayX3[fdir] = c1oS3;
    fdir        = DIR_MPM;
    rayX1[fdir] = -c1oS3;
    rayX2[fdir] = c1oS3;
    rayX3[fdir] = -c1oS3;
    fdir        = DIR_PPM;
    rayX1[fdir] = c1oS3;
    rayX2[fdir] = c1oS3;
    rayX3[fdir] = -c1oS3;
    fdir        = DIR_MMM;
    rayX1[fdir] = -c1oS3;
    rayX2[fdir] = -c1oS3;
    rayX3[fdir] = -c1oS3;
    fdir        = DIR_PMM;
    rayX1[fdir] = c1oS3;
    rayX2[fdir] = -c1oS3;
    rayX3[fdir] = -c1oS3;
}
//////////////////////////////////////////////////////////////////////////
void D3Q27Interactor::initInteractor(const real &timeStep)
{
    UBLOG(logDEBUG5, "D3Q27Interactor::initInteractor - "
                         << " for timestep = " << timeStep);

    //////////////////////////////////////////////////////////////////////////
    // init bcs
    int nofAdapter = (int)bcAdapters.size();
    if (nofAdapter == 0) {
        UBLOG(logWARNING, "WARNING - D3Q27Interactor::initInteractor Warning - no nodeAdapter available");
    }
    bool needTimeDependence = false;
    for (int pos = 0; pos < nofAdapter; ++pos) {
        bcAdapters[pos]->init(this, timeStep);
        if (bcAdapters[pos]->isTimeDependent())
            needTimeDependence = true;
    }
    if (needTimeDependence)
        this->setTimeDependent();
    else
        this->unsetTimeDependent();

    updateBlocks();
}
//////////////////////////////////////////////////////////////////////////
void D3Q27Interactor::updateInteractor(const real &timestep)
{
    UBLOG(logDEBUG5, "D3Q27Interactor::updateInteractor - for timestep = " << timestep);

    //////////////////////////////////////////////////////////////////////////
    // update bcs
    int nofAdapter = (int)bcAdapters.size();
    if (nofAdapter == 0) {
        UBLOG(logERROR, "WARNING - D3Q27Interactor::updateInteractor Warning - no nodeAdapter available for ");
    }

    bool needTimeDependence = false;

    for (int pos = 0; pos < nofAdapter; ++pos) {
        bcAdapters[pos]->update(this, timestep);
        if (bcAdapters[pos]->isTimeDependent())
            needTimeDependence = true;
    }
    if (needTimeDependence)
        this->setTimeDependent();
    else
        this->unsetTimeDependent();

    for (BcNodeIndicesMap::value_type t : bcNodeIndicesMap) {
        SPtr<Block3D> block                             = t.first;
        std::set<std::vector<int>> &transNodeIndicesSet = t.second;

        if (block->isNotActive() || !block)
            continue;

        SPtr<ILBMKernel> kernel = block->getKernel();
        SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();

        set<std::vector<int>>::iterator setPos;

        for (setPos = transNodeIndicesSet.begin(); setPos != transNodeIndicesSet.end(); ++setPos) {
            int x1          = (*setPos)[0];
            int x2          = (*setPos)[1];
            int x3          = (*setPos)[2];
            Vector3D coords = grid.lock()->getNodeCoordinates(block, x1, x2, x3);
            real worldX1  = coords[0];
            real worldX2  = coords[1];
            real worldX3  = coords[2];

            SPtr<BoundaryConditions> bc = bcArray->getBC(x1, x2, x3);
            if (bc) // may be that the BC has been deleted by the solid setting of another interactor
            {
                for (size_t i = 0; i < bcAdapters.size(); i++)
                    bcAdapters[i]->adaptBC(*this, bc, worldX1, worldX2, worldX3, timestep);
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
// calculation takes place in the real coordinate system !!!
// not normalized!
// x1, x2, x3 are the coordinates at the bottom left of the "system"
// extendedBoundingGeoOfGeoObject MUST already have been magnified by delta_x_level in each direction for SOLID
bool D3Q27Interactor::setDifferencesToGbObject3D(const SPtr<Block3D> block)
{
    using namespace vf::lbm::dir;

    if (!block)
        return false;

    if (block->isNotActive())
        return false; // continue;

    bcNodeIndicesMap[block]                 = set<std::vector<int>>();
    set<std::vector<int>> &transNodeIndices = bcNodeIndicesMap[block];
    solidNodeIndicesMap[block]              = set<UbTupleInt3>();
    set<UbTupleInt3> &solidNodeIndices      = solidNodeIndicesMap[block];

    real timestep    = 0;
    bool oneEntryGotBC = false;
    bool gotQs         = false;
    SPtr<BoundaryConditions> bc;

    SPtr<ILBMKernel> kernel = block->getKernel();
    SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();

    real internX1, internX2, internX3;

    int startIX1 = 0;
    int startIX2 = 0;
    int startIX3 = 0;
    int stopIX1  = (int)bcArray->getNX1();
    int stopIX2  = (int)bcArray->getNX2();
    int stopIX3  = (int)bcArray->getNX3();

    real dx = grid.lock()->getDeltaX(block);

    // other boundingRect than in init, because here the boundrect has to be increased by one dx
    GbCuboid3D extendedBoundingGeoOfGeoObject(
        geoObject3D->getX1Minimum() - 1.02 * dx, geoObject3D->getX2Minimum() - 1.02 * dx,
        geoObject3D->getX3Minimum() - 1.02 * dx, geoObject3D->getX1Maximum() + 1.02 * dx,
        geoObject3D->getX2Maximum() + 1.02 * dx, geoObject3D->getX3Maximum() + 1.02 * dx);

    real deltaX1 = dx, deltaX2 = dx, deltaX3 = dx;

    if (geoObject3D->hasRaytracing() || (this->isInverseSolid() && geoObject3D->raytracingSupportsPointsInside())) {
        // if deltaX1==deltaX2==deltaX3 (must for LB!!)
        if (!UbMath::zero(deltaX1 - deltaX2 + deltaX1 - deltaX3 + deltaX2 - deltaX3))
            throw UbException(
                UB_EXARGS, "fuer den bei LB nicht vorkommenden Fall deltaX1!=deltaX2!=deltaX3  nicht implementiert ");

        vector<real> distNeigh(D3Q27System::FENDDIR + 1, UbMath::sqrt2 * deltaX1);
        distNeigh[DIR_P00] = distNeigh[DIR_M00] = distNeigh[DIR_0P0] = deltaX1;
        distNeigh[DIR_0M0] = distNeigh[DIR_00P] = distNeigh[DIR_00M] = deltaX1;
        distNeigh[DIR_PP0] = distNeigh[DIR_MP0] = distNeigh[DIR_MM0] =
            distNeigh[DIR_PM0]             = UbMath::sqrt2 * deltaX1;
        distNeigh[DIR_P0P] = distNeigh[DIR_0PP] = distNeigh[DIR_M0P] =
            distNeigh[DIR_0MP]             = UbMath::sqrt2 * deltaX1;
        distNeigh[DIR_P0M] = distNeigh[DIR_0PM] = distNeigh[DIR_M0M] =
            distNeigh[DIR_0MM]             = UbMath::sqrt2 * deltaX1;
        distNeigh[DIR_PPP] = distNeigh[DIR_MPP] = distNeigh[DIR_PMP] =
            distNeigh[DIR_MMP]              = UbMath::sqrt3 * deltaX1;
        distNeigh[DIR_PPM] = distNeigh[DIR_MPM] = distNeigh[DIR_PMM] =
            distNeigh[DIR_MMM]              = UbMath::sqrt3 * deltaX1;
        real q;
        bool pointOnBoundary = false;

        //#ifdef _OPENMP
        //      #pragma omp parallel for private(internX1,internX2,internX3,gotQs,bc,q )
        //#endif
        for (int ix3 = startIX3; ix3 < stopIX3; ix3++) {
            for (int ix2 = startIX2; ix2 < stopIX2; ix2++) {
                for (int ix1 = startIX1; ix1 < stopIX1; ix1++) {
                    // TODO: further, investigate if this is not a mistake
                    if (bcArray->isUndefined(ix1, ix2, ix3))
                        continue;

                    Vector3D coords = grid.lock()->getNodeCoordinates(block, ix1, ix2, ix3);
                    internX1        = coords[0];
                    internX2        = coords[1];
                    internX3        = coords[2];

                    // Point in the object test is superfluous, since the start and stop indices already exist
                    // are determined -> only point-in-cube indexes are considered
                    if (extendedBoundingGeoOfGeoObject.isPointInGbObject3D(internX1, internX2, internX3)) {
                        if (this->isSolid()) {
                            if (this->geoObject3D->isPointInGbObject3D(internX1, internX2, internX3)) {
                                {
                                    solidNodeIndices.insert(UbTupleInt3(ix1, ix2, ix3));
                                    bcArray->setSolid(ix1, ix2, ix3);
                                }
                                continue;
                            }
                        } else if (this->isInverseSolid()) {
                            // in inverse solid all nodes are OUTSIDE and on the boundary SOLID
                            if (!this->geoObject3D->isPointInGbObject3D(internX1, internX2, internX3,
                                                                        pointOnBoundary) ||
                                pointOnBoundary == true) {

                                {
                                    solidNodeIndices.insert(UbTupleInt3(ix1, ix2, ix3));
                                    bcArray->setSolid(ix1, ix2, ix3);
                                }
                                continue;
                            }
                        }

                        if (bcArray->isSolid(ix1, ix2, ix3))
                            continue;

                        gotQs = false;

                        for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
                            q = geoObject3D->getIntersectionRaytraceFactor(internX1, internX2, internX3, rayX1[fdir],
                                                                           rayX2[fdir], rayX3[fdir]);
                            q /= distNeigh[fdir];

                            // assert(UbMath::lessEqual(q, 1.0));

                            if (UbMath::inClosedInterval(q, 1.0, 1.0))
                                q = 1.0;
                            if (UbMath::greater(q, 0.0) && UbMath::lessEqual(q, 1.0)) {
                                //#pragma omp critical (BC_CHANGE)
                                {
                                    bc = bcArray->getBC(ix1, ix2, ix3);
                                    if (!bc) {
                                        bc = std::make_shared<BoundaryConditions>();
                                        bcArray->setBC(ix1, ix2, ix3, bc);
                                    }

                                    if (bc->hasNoSlipBoundary()) {
                                        bc->setBoundaryVelocityX1(0.0);
                                        bc->setBoundaryVelocityX2(0.0);
                                        bc->setBoundaryVelocityX3(0.0);
                                    }

                                    for (int index = (int)bcAdapters.size() - 1; index >= 0; --index)
                                        bcAdapters[index]->adaptBCForDirection(*this, bc, internX1, internX2, internX3,
                                                                               q, fdir, timestep);
                                }

                                gotQs = true;
                            }
                        }

                        if (gotQs) {
                            {
                                oneEntryGotBC = true;

                                std::vector<int> p(3);
                                p[0] = ix1;
                                p[1] = ix2;
                                p[2] = ix3;
                                transNodeIndices.insert(p);

                                for (int index = (int)bcAdapters.size() - 1; index >= 0; --index)
                                    bcAdapters[index]->adaptBC(*this, bc, internX1, internX2, internX3, timestep);
                            }
                        }
                    } else if (this->isInverseSolid()) {
                        // bei inverse solid sind alle Knoten AUSSERHALB und auf der boundary SOLID
                        if (!this->geoObject3D->isPointInGbObject3D(internX1, internX2, internX3, pointOnBoundary) ||
                            pointOnBoundary == true) {
                            {
                                solidNodeIndices.insert(UbTupleInt3(ix1, ix2, ix3));
                                bcArray->setSolid(ix1, ix2, ix3);
                            }
                            continue;
                        }
                    }
                }
            }
        }
    } else // clipping -> slower (currently also used for all inverse Solid objects whose raytracing does not work for
           // nodes INSIDE the geo)
    {
        bool pointOnBoundary = false;
        for (int ix1 = startIX1; ix1 < stopIX1; ix1++) {
            for (int ix2 = startIX2; ix2 < stopIX2; ix2++) {
                for (int ix3 = startIX3; ix3 < stopIX3; ix3++) {
                    if (bcArray->isSolid(ix1, ix2, ix3) || bcArray->isUndefined(ix1, ix2, ix3))
                        continue;

                    Vector3D coords = grid.lock()->getNodeCoordinates(block, ix1, ix2, ix3);
                    internX1        = coords[0];
                    internX2        = coords[1];
                    internX3        = coords[2];

                    if (extendedBoundingGeoOfGeoObject.isPointInGbObject3D(internX1, internX2, internX3)) {
                        if (this->isSolid() && this->geoObject3D->isPointInGbObject3D(internX1, internX2, internX3)) {
                            {
                                solidNodeIndices.insert(UbTupleInt3(ix1, ix2, ix3));
                                bcArray->setSolid(ix1, ix2, ix3);
                            }
                            continue;
                        } else if (this->isInverseSolid()) {
                            // bei inverse solid sind alle Knoten AUSSERHALB und auf der boundary SOLID
                            if (!this->geoObject3D->isPointInGbObject3D(internX1, internX2, internX3,
                                                                        pointOnBoundary) ||
                                pointOnBoundary == true) {
                                {
                                    solidNodeIndices.insert(UbTupleInt3(ix1, ix2, ix3));
                                    bcArray->setSolid(ix1, ix2, ix3);
                                }
                                continue;
                            }
                        }

                        gotQs = false;

                        GbPoint3D pointA(internX1, internX2, internX3);
                        for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
                            real x1B = internX1 + D3Q27System::DX1[fdir] * deltaX1;
                            real x2B = internX2 + D3Q27System::DX2[fdir] * deltaX2;
                            real x3B = internX3 + D3Q27System::DX3[fdir] * deltaX3;

                            GbPoint3D pointB(x1B, x2B, x3B);
                            GbLine3D *clippedLine = this->geoObject3D->createClippedLine3D(pointA, pointB);

                            if (clippedLine) {
                                real q = 0.0;
                                if (!this->isInverseSolid()) // A is outside
                                {
                                    real distanceAB = pointA.getDistance(&pointB); // pointA to B
                                    real distanceAP = UbMath::min(pointA.getDistance(clippedLine->getPoint1()),
                                                                    pointA.getDistance(clippedLine->getPoint2()));
                                    q                 = distanceAP / distanceAB;
                                } else {
                                    bool pointIsOnBoundary = false;
                                    if (!clippedLine->getPoint1()->equals(&pointB) &&
                                        !clippedLine->getPoint2()->equals(&pointB)) {
                                        // A is inside, a clipped line must not contain B
                                        real distanceAB = pointA.getDistance(&pointB); // pointA to B
                                        real distanceAP = clippedLine->getLength();
                                        q                 = distanceAP / distanceAB;
                                    } else if (this->geoObject3D->isPointInGbObject3D(
                                                   pointB.getX1Coordinate(), pointB.getX2Coordinate(),
                                                   pointB.getX3Coordinate(), pointIsOnBoundary) &&
                                               pointIsOnBoundary) {
                                        // A is definitely inside, B is exactly on ObjectBoundary => q = 1.0
                                        q = 1.0;
                                    } else {
                                        q = 0.0;
                                    }
                                }

                                if (UbMath::inClosedInterval(q, 1.0, 1.0))
                                    q = 1.0;
                                if (UbMath::lessEqual(q, 1.0) && UbMath::greater(q, 0.0)) {
                                    {
                                        bc = bcArray->getBC(ix1, ix2, ix3);
                                        if (!bc) {
                                            bc = std::make_shared<BoundaryConditions>();
                                            bcArray->setBC(ix1, ix2, ix3, bc);
                                        }
                                        for (int index = (int)bcAdapters.size() - 1; index >= 0; --index)
                                            bcAdapters[index]->adaptBCForDirection(*this, bc, internX1, internX2,
                                                                                   internX3, q, fdir, timestep);
                                    }

                                    gotQs = true;
                                }

                                clippedLine->deletePoint1();
                                clippedLine->deletePoint2();
                                delete clippedLine;
                            }
                        }

                        if (gotQs) {
                            {
                                oneEntryGotBC = true;

                                std::vector<int> p(3);
                                p[0] = ix1;
                                p[1] = ix2;
                                p[2] = ix3;
                                transNodeIndices.insert(p);

                                for (int index = (int)bcAdapters.size() - 1; index >= 0; --index)
                                    bcAdapters[index]->adaptBC(*this, bc, internX1, internX2, internX3, timestep);
                            }
                        }
                    }
                }
            }
        }
    }

    return oneEntryGotBC;
}
//////////////////////////////////////////////////////////////////////////
void D3Q27Interactor::addQsLineSet(std::vector<UbTupleFloat3> &nodes, std::vector<UbTupleInt2> &lines)
{
    using namespace vf::lbm::dir;

    for (SPtr<Block3D> block : bcBlocks) {
        if (!block)
            continue;

        real dx               = grid.lock()->getDeltaX(block);
        UbTupleDouble3 orgDelta = grid.lock()->getNodeOffset(block);

        SPtr<ILBMKernel> kernel = block->getKernel();
        SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();

        map<SPtr<Block3D>, set<std::vector<int>>>::iterator pos = bcNodeIndicesMap.find(block);
        if (pos == bcNodeIndicesMap.end()) {
            UB_THROW(UbException(UB_EXARGS, "block nicht in indizes map!!!"));
        }
        set<std::vector<int>> &transNodeIndicesSet = pos->second;
        set<std::vector<int>>::iterator setPos;

        std::size_t node1Index, node2Index;

        UbTupleDouble3 blockOrg = grid.lock()->getBlockWorldCoordinates(block);

        for (setPos = transNodeIndicesSet.begin(); setPos != transNodeIndicesSet.end(); ++setPos) {
            int ix1 = (*setPos)[0];
            int ix2 = (*setPos)[1];
            int ix3 = (*setPos)[2];

            if (bcArray->isFluid(
                    ix1, ix2,
                    ix3)) // it may be that the node is replaced by another interactor e.g. was marked as solid !!!
            {
                if (!bcArray->hasBC(ix1, ix2, ix3))
                    continue;
                SPtr<BoundaryConditions> bc = bcArray->getBC(ix1, ix2, ix3);

                real x1a = val<1>(blockOrg) - val<1>(orgDelta) + ix1 * dx;
                real x2a = val<2>(blockOrg) - val<2>(orgDelta) + ix2 * dx;
                real x3a = val<3>(blockOrg) - val<3>(orgDelta) + ix3 * dx;
                nodes.push_back(makeUbTuple((float)x1a, (float)x2a, (float)x3a));
                node1Index = nodes.size() - 1;

                for (int dir = D3Q27System::FSTARTDIR; dir <= D3Q27System::FENDDIR; dir++) {
                    if (bc->hasBoundaryConditionFlag(D3Q27System::INVDIR[dir])) {
                        real x1b, x2b, x3b, q = bc->getQ(dir);
                        switch (dir) {
                            case DIR_P00:
                                x1b = x1a + q * dx;
                                x2b = x2a;
                                x3b = x3a;
                                break;
                            case DIR_0P0:
                                x1b = x1a;
                                x2b = x2a + q * dx;
                                x3b = x3a;
                                break;
                            case DIR_M00:
                                x1b = x1a - q * dx;
                                x2b = x2a;
                                x3b = x3a;
                                break;
                            case DIR_0M0:
                                x1b = x1a;
                                x2b = x2a - q * dx;
                                x3b = x3a;
                                break;
                            case DIR_PP0:
                                x1b = x1a + q * dx;
                                x2b = x2a + q * dx;
                                x3b = x3a;
                                break;
                            case DIR_MP0:
                                x1b = x1a - q * dx;
                                x2b = x2a + q * dx;
                                x3b = x3a;
                                break;
                            case DIR_MM0:
                                x1b = x1a - q * dx;
                                x2b = x2a - q * dx;
                                x3b = x3a;
                                break;
                            case DIR_PM0:
                                x1b = x1a + q * dx;
                                x2b = x2a - q * dx;
                                x3b = x3a;
                                break;
                            case DIR_00P:
                                x1b = x1a;
                                x2b = x2a;
                                x3b = x3a + q * dx;
                                break;
                            case DIR_P0P:
                                x1b = x1a + q * dx;
                                x2b = x2a;
                                x3b = x3a + q * dx;
                                break;
                            case DIR_0PP:
                                x1b = x1a;
                                x2b = x2a + q * dx;
                                x3b = x3a + q * dx;
                                break;
                            case DIR_M0P:
                                x1b = x1a - q * dx;
                                x2b = x2a;
                                x3b = x3a + q * dx;
                                break;
                            case DIR_0MP:
                                x1b = x1a;
                                x2b = x2a - q * dx;
                                x3b = x3a + q * dx;
                                break;
                            case DIR_00M:
                                x1b = x1a;
                                x2b = x2a;
                                x3b = x3a - q * dx;
                                break;
                            case DIR_P0M:
                                x1b = x1a + q * dx;
                                x2b = x2a;
                                x3b = x3a - q * dx;
                                break;
                            case DIR_0PM:
                                x1b = x1a;
                                x2b = x2a + q * dx;
                                x3b = x3a - q * dx;
                                break;
                            case DIR_M0M:
                                x1b = x1a - q * dx;
                                x2b = x2a;
                                x3b = x3a - q * dx;
                                break;
                            case DIR_0MM:
                                x1b = x1a;
                                x2b = x2a - q * dx;
                                x3b = x3a - q * dx;
                                break;
                            case DIR_PPP:
                                x1b = x1a + q * dx;
                                x2b = x2a + q * dx;
                                x3b = x3a + q * dx;
                                break;
                            case DIR_MMM:
                                x1b = x1a - q * dx;
                                x2b = x2a - q * dx;
                                x3b = x3a - q * dx;
                                break;
                            case DIR_PPM:
                                x1b = x1a + q * dx;
                                x2b = x2a + q * dx;
                                x3b = x3a - q * dx;
                                break;
                            case DIR_MMP:
                                x1b = x1a - q * dx;
                                x2b = x2a - q * dx;
                                x3b = x3a + q * dx;
                                break;
                            case DIR_PMP:
                                x1b = x1a + q * dx;
                                x2b = x2a - q * dx;
                                x3b = x3a + q * dx;
                                break;
                            case DIR_MPM:
                                x1b = x1a - q * dx;
                                x2b = x2a + q * dx;
                                x3b = x3a - q * dx;
                                break;
                            case DIR_PMM:
                                x1b = x1a + q * dx;
                                x2b = x2a - q * dx;
                                x3b = x3a - q * dx;
                                break;
                            case DIR_MPP:
                                x1b = x1a - q * dx;
                                x2b = x2a + q * dx;
                                x3b = x3a + q * dx;
                                break;
                            default:
                                throw UbException(UB_EXARGS, "unknown direction");
                        }

                        nodes.push_back(makeUbTuple((float)x1b, (float)x2b, (float)x3b));
                        node2Index = nodes.size() - 1;

                        lines.push_back(makeUbTuple((int)node1Index, (int)node2Index));
                    }
                }
            }
        }
    }
}
////////////////////////////////////////////////////////////////////////////
vector<pair<GbPoint3D, GbPoint3D>> D3Q27Interactor::getQsLineSet()
{
    using namespace vf::lbm::dir;

    vector<pair<GbPoint3D, GbPoint3D>> QsLineSet;
    pair<GbPoint3D, GbPoint3D> pointpair;

    UbTupleInt3 blocknx = grid.lock()->getBlockNX();

    int blocknx1 = val<1>(blocknx);
    int blocknx2 = val<2>(blocknx);
    int blocknx3 = val<3>(blocknx);

    for (SPtr<Block3D> block : bcBlocks) {
        SPtr<ILBMKernel> kernel   = block->getKernel();
        SPtr<BCArray3D> bcMatrix  = kernel->getBCProcessor()->getBCArray();
        UbTupleDouble3 nodeOffset = grid.lock()->getNodeOffset(block);

        // Check whether top row is real in the system or not
        bool include_N_Face  = false; // x1=[0..blocknx1[ && x3=[0..blocknx3[
        bool include_E_Face  = false; // x2=[0..blocknx2[ && x3=[0..blocknx3[
        bool include_T_Face  = false; // x1=[0..blocknx1[ && x2=[0..blocknx2[
        bool include_NE_Edge = false; //(x1/x2/x3)=(blocknx1/blocknx2/[0..blocknx3[)
        bool include_TN_Edge = false; //(x1/x2/x3)=([0..blocknx1[/blocknx2/blocknx1)
        bool include_TE_Edge = false; //(x1/x2/x3)=(blocknx1/[0..blocknx2[/blocknx2)
        if (block) {
            if (!block->getConnector(DIR_0P0))
                include_N_Face = true;
            if (!block->getConnector(DIR_P00))
                include_E_Face = true;
            if (!block->getConnector(DIR_00P))
                include_T_Face = true;
            if (!block->getConnector(DIR_PP0) && include_N_Face && include_E_Face)
                include_NE_Edge = true;
            if (!block->getConnector(DIR_0PP) && include_T_Face && include_N_Face)
                include_TN_Edge = true;
            if (!block->getConnector(DIR_P0P) && include_T_Face && include_E_Face)
                include_TE_Edge = true;
        }

        map<SPtr<Block3D>, set<std::vector<int>>>::iterator pos = bcNodeIndicesMap.find(block);
        if (pos == bcNodeIndicesMap.end())
            throw UbException(UB_EXARGS, "block nicht in indizes map!!!" + block->toString());
        set<std::vector<int>> &transNodeIndicesSet = pos->second;
        set<std::vector<int>>::iterator setPos;

        real x1, x2, x3, dx;
        grid.lock()->calcStartCoordinatesAndDelta(block, x1, x2, x3, dx);

        for (setPos = transNodeIndicesSet.begin(); setPos != transNodeIndicesSet.end(); ++setPos) {
            int ix1 = (*setPos)[0];
            int ix2 = (*setPos)[1];
            int ix3 = (*setPos)[2];

            if ((ix1 < blocknx1 && ix2 < blocknx2 && ix3 < blocknx3) ||
                (include_E_Face && ix1 == blocknx1 && ix2 < blocknx2 && ix3 < blocknx3) ||
                (include_N_Face && ix2 == blocknx2 && ix1 < blocknx1 && ix3 < blocknx3) ||
                (include_T_Face && ix3 == blocknx3 && ix1 < blocknx1 && ix2 < blocknx2) ||
                (include_NE_Edge && ix1 == blocknx1 && ix2 == blocknx2) ||
                (include_TN_Edge && ix2 == blocknx2 && ix3 == blocknx3) ||
                (include_TE_Edge && ix1 == blocknx1 && ix3 == blocknx3)) {
                if (bcMatrix->isFluid(
                        ix1, ix2,
                        ix3)) // it may be that the node is replaced by another interactor e.g. was marked as solid !!!
                {
                    if (!bcMatrix->hasBC(ix1, ix2, ix3))
                        continue;
                    SPtr<BoundaryConditions> bc = bcMatrix->getBC(ix1, ix2, ix3);
                    real x1a                  = x1 - val<1>(nodeOffset) + dx * ix1;
                    real x2a                  = x2 - val<2>(nodeOffset) + dx * ix2;
                    real x3a                  = x3 - val<3>(nodeOffset) + dx * ix3;
                    pointpair.first.setX1(x1a);
                    pointpair.first.setX2(x2a);
                    pointpair.first.setX3(x3a);
                    for (int dir = D3Q27System::FSTARTDIR; dir <= D3Q27System::FENDDIR; dir++) {
                        if (bc->hasBoundaryConditionFlag(D3Q27System::INVDIR[dir])) {
                            real x1b, x2b, x3b, q = bc->getQ(dir);
                            switch (dir) {
                                case DIR_P00:
                                    x1b = x1a + q * dx;
                                    x2b = x2a;
                                    x3b = x3a;
                                    break;
                                case DIR_0P0:
                                    x1b = x1a;
                                    x2b = x2a + q * dx;
                                    x3b = x3a;
                                    break;
                                case DIR_M00:
                                    x1b = x1a - q * dx;
                                    x2b = x2a;
                                    x3b = x3a;
                                    break;
                                case DIR_0M0:
                                    x1b = x1a;
                                    x2b = x2a - q * dx;
                                    x3b = x3a;
                                    break;
                                case DIR_PP0:
                                    x1b = x1a + q * dx;
                                    x2b = x2a + q * dx;
                                    x3b = x3a;
                                    break;
                                case DIR_MP0:
                                    x1b = x1a - q * dx;
                                    x2b = x2a + q * dx;
                                    x3b = x3a;
                                    break;
                                case DIR_MM0:
                                    x1b = x1a - q * dx;
                                    x2b = x2a - q * dx;
                                    x3b = x3a;
                                    break;
                                case DIR_PM0:
                                    x1b = x1a + q * dx;
                                    x2b = x2a - q * dx;
                                    x3b = x3a;
                                    break;
                                case DIR_00P:
                                    x1b = x1a;
                                    x2b = x2a;
                                    x3b = x3a + q * dx;
                                    break;
                                case DIR_P0P:
                                    x1b = x1a + q * dx;
                                    x2b = x2a;
                                    x3b = x3a + q * dx;
                                    break;
                                case DIR_0PP:
                                    x1b = x1a;
                                    x2b = x2a + q * dx;
                                    x3b = x3a + q * dx;
                                    break;
                                case DIR_M0P:
                                    x1b = x1a - q * dx;
                                    x2b = x2a;
                                    x3b = x3a + q * dx;
                                    break;
                                case DIR_0MP:
                                    x1b = x1a;
                                    x2b = x2a - q * dx;
                                    x3b = x3a + q * dx;
                                    break;
                                case DIR_00M:
                                    x1b = x1a;
                                    x2b = x2a;
                                    x3b = x3a - q * dx;
                                    break;
                                case DIR_P0M:
                                    x1b = x1a + q * dx;
                                    x2b = x2a;
                                    x3b = x3a - q * dx;
                                    break;
                                case DIR_0PM:
                                    x1b = x1a;
                                    x2b = x2a + q * dx;
                                    x3b = x3a - q * dx;
                                    break;
                                case DIR_M0M:
                                    x1b = x1a - q * dx;
                                    x2b = x2a;
                                    x3b = x3a - q * dx;
                                    break;
                                case DIR_0MM:
                                    x1b = x1a;
                                    x2b = x2a - q * dx;
                                    x3b = x3a - q * dx;
                                    break;
                                case DIR_PPP:
                                    x1b = x1a + q * dx;
                                    x2b = x2a + q * dx;
                                    x3b = x3a + q * dx;
                                    break;
                                case DIR_MMM:
                                    x1b = x1a - q * dx;
                                    x2b = x2a - q * dx;
                                    x3b = x3a - q * dx;
                                    break;
                                case DIR_PPM:
                                    x1b = x1a + q * dx;
                                    x2b = x2a + q * dx;
                                    x3b = x3a - q * dx;
                                    break;
                                case DIR_MMP:
                                    x1b = x1a - q * dx;
                                    x2b = x2a - q * dx;
                                    x3b = x3a + q * dx;
                                    break;
                                case DIR_PMP:
                                    x1b = x1a + q * dx;
                                    x2b = x2a - q * dx;
                                    x3b = x3a + q * dx;
                                    break;
                                case DIR_MPM:
                                    x1b = x1a - q * dx;
                                    x2b = x2a + q * dx;
                                    x3b = x3a - q * dx;
                                    break;
                                case DIR_PMM:
                                    x1b = x1a + q * dx;
                                    x2b = x2a - q * dx;
                                    x3b = x3a - q * dx;
                                    break;
                                case DIR_MPP:
                                    x1b = x1a - q * dx;
                                    x2b = x2a + q * dx;
                                    x3b = x3a + q * dx;
                                    break;
                                default:
                                    throw UbException(UB_EXARGS, "unknown direction");
                            }
                            pointpair.second.setX1(x1b);
                            pointpair.second.setX2(x2b);
                            pointpair.second.setX3(x3b);
                            QsLineSet.push_back(pointpair);
                        }
                    }
                }
            }
        }
    }
    return QsLineSet;
}

void D3Q27Interactor::writeValidationAVSFile(string filename)
{
    UBLOG(logINFO, "D3Q27Interactor::writeValidationAVSFile(" << filename << ") - start ");
    ofstream out(filename.c_str(), ios::out);
    if (!out)
        throw UbException(UB_EXARGS, "couldn't open file " + filename);

    int numpoints, numlines;
    vector<pair<GbPoint3D, GbPoint3D>> qsLineSet = this->getQsLineSet();
    numlines                                     = (unsigned)qsLineSet.size();
    numpoints                                    = numlines * 2;

    out << "# UCD-File created by D3Q27Interactor\n";
    out << numpoints << " " << numlines << " 0 0 0 " << endl;
    int nr = 1;
    for (int i = 0; i < numlines; i++) {
        out << nr++ << " " << qsLineSet[i].first.getX1Coordinate() << " " << qsLineSet[i].first.getX2Coordinate() << " "
            << qsLineSet[i].first.getX3Coordinate() << " \n";
        out << nr++ << " " << qsLineSet[i].second.getX1Coordinate() << " " << qsLineSet[i].second.getX2Coordinate()
            << " " << qsLineSet[i].second.getX3Coordinate() << " \n";
    }
    nr = 1;
    for (int i = 0; i < numlines; i++) {
        int el = nr + 1;
        out << i + 1 << " " << 2 << " line " << nr << " " << el << " " << endl;
        nr = el + 1;
    }
    UBLOG(logINFO, "D3Q27Interactor::writeValidationAVSFile(" << filename << ") - end");
}
