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
//! \addtogroup cpu_Interactors Interactors
//! \ingroup cpu_core core
//! \{
//! \author Sören Freudiger
//! \author Sebastian Geller
//! \author Ehsan Kian Far
//! \author Konstantin Kutscher
//=======================================================================================

#include "D3Q27TriFaceMeshInteractor.h"
#include <basics/utilities/UbLogger.h>
#include <basics/utilities/UbMath.h>

#include <basics/writer/WbWriterVtkXmlASCII.h>
#include <basics/writer/WbWriterVtkASCII.h>
#include <basics/writer/WbWriterVtkBinary.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>

#include <basics/Timer/Timer.h>

#include "BCArray3D.h"
#include "BCSet.h"
#include "Block3D.h"
#include "BoundaryConditions.h"
#include "Grid3D.h"
#include "LBMKernel.h"
#include "VelocityBC.h"
#include <geometry3d/GbCuboid3D.h>
#include <geometry3d/GbHalfSpace3D.h>
#include <geometry3d/GbMeshTools3D.h>
#include <geometry3d/GbSystem3D.h>

#include <geometry3d/GbTriFaceMesh3D.h>

//#include <omp.h>

#include "CoordinateTransformation3D.h"
#include <stack>

using namespace std;

D3Q27TriFaceMeshInteractor::D3Q27TriFaceMeshInteractor() : D3Q27Interactor() {  }
//////////////////////////////////////////////////////////////////////////
D3Q27TriFaceMeshInteractor::D3Q27TriFaceMeshInteractor(SPtr<Grid3D> /*grid*/, std::string /*name*/)
{
    
}
//////////////////////////////////////////////////////////////////////////
D3Q27TriFaceMeshInteractor::D3Q27TriFaceMeshInteractor(SPtr<GbTriFaceMesh3D> triFaceMesh, SPtr<Grid3D> grid,
                                                       SPtr<BC> BC, int type)
    : D3Q27Interactor(triFaceMesh, grid, BC, type)
{
    
}
//////////////////////////////////////////////////////////////////////////
D3Q27TriFaceMeshInteractor::D3Q27TriFaceMeshInteractor(SPtr<GbTriFaceMesh3D> triFaceMesh, SPtr<Grid3D> grid,
                                                       SPtr<BC> BC, int type, Interactor3D::Accuracy a)
    : D3Q27Interactor(triFaceMesh, grid, BC, type, a)
{
   
}
//////////////////////////////////////////////////////////////////////////
D3Q27TriFaceMeshInteractor::~D3Q27TriFaceMeshInteractor() = default;
//////////////////////////////////////////////////////////////////////////
void D3Q27TriFaceMeshInteractor::initInteractor(const real &timeStep)
{
    updateBlocks(); 
    setQs(timeStep);
}
//////////////////////////////////////////////////////////////////////////
bool D3Q27TriFaceMeshInteractor::setDifferencesToGbObject3D(const SPtr<Block3D> block)
{
    if (!block)
        return false;

    // UBLOG(logINFO, "D3Q27TriFaceMeshInteractor::setDifferencesToGbObject3D()");

    bcNodeIndicesMap[block] = set<std::vector<int>>();
    //   set< std::vector<int> >& transNodeIndices = bcNodeIndicesMap[block];
    solidNodeIndicesMap[block]         = set<UbTupleInt3>();
    set<UbTupleInt3> &solidNodeIndices = solidNodeIndicesMap[block];

    bool oneEntryGotBC = false; 
    SPtr<BoundaryConditions> bc;

    SPtr<ILBMKernel> kernel = block->getKernel();
    SPtr<BCArray3D> bcArray = kernel->getBCSet()->getBCArray();

    double internX1, internX2, internX3;

    int startIX1 = 0, startIX2 = 0, startIX3 = 0;
    int stopIX1 = (int)bcArray->getNX1(), stopIX2 = (int)bcArray->getNX2(), stopIX3 = (int)bcArray->getNX3();

    bool pointOnBoundary = false;

    for (int ix3 = startIX3; ix3 < stopIX3; ix3++) {
        for (int ix2 = startIX2; ix2 < stopIX2; ix2++) {
            for (int ix1 = startIX1; ix1 < stopIX1; ix1++) {
                Vector3D coords = grid.lock()->getNodeCoordinates(block, ix1, ix2, ix3);
                internX1        = coords[0];
                internX2        = coords[1];
                internX3        = coords[2];

                if (this->isSolid()) {
                    if (this->geoObject3D->isPointInGbObject3D(internX1, internX2, internX3)) {
                        if (bcArray->isFluid(ix1, ix2, ix3)) {
                            solidNodeIndices.insert(UbTupleInt3(ix1, ix2, ix3));
                            bcArray->setSolid(ix1, ix2, ix3);
                        }
                    }
                } else if (this->isInverseSolid()) {
                    // with inverse solid all nodes are OUTSIDE and on the boundary SOLID
                    if (!this->geoObject3D->isPointInGbObject3D(internX1, internX2, internX3, pointOnBoundary) ||
                        pointOnBoundary == true) {
                        if (bcArray->isFluid(ix1, ix2, ix3)) {
                            solidNodeIndices.insert(UbTupleInt3(ix1, ix2, ix3));
                            bcArray->setSolid(ix1, ix2, ix3);
                        }
                    }
                }
            }
        }
    }

    return oneEntryGotBC;
}
//////////////////////////////////////////////////////////////////////////
// E.F. /4/16/2013
void D3Q27TriFaceMeshInteractor::setQs(const real &timeStep)
{
    using namespace vf::lbm::dir;
    using namespace vf::basics::constant;

    UBLOGML(logDEBUG1, "\nLBMTriFaceMeshInteractor - setQs start ");
    if (!this->grid.lock())
        throw UbException(UB_EXARGS, "ups, no grid.lock()!!");

    if (this->reinitWithStoredQsFlag && !bcNodeIndicesAndQsMap.empty()) {
        this->reinitWithStoredQs(timeStep);
        return;
    }

    GbTriFaceMesh3D *mesh = dynamic_cast<GbTriFaceMesh3D *>(this->geoObject3D.get());

    //////////////////////////////////////////////////////////////////////////
    // init bcs
    //////////////////////////////////////////////////////////////////////////
    int nofAdapter = (int)this->BCs.size();
    if (nofAdapter == 0)
        std::cout
            << "WARNING - D3Q27TriFaceMeshInteractor::initInteractor Warning - no nodeAdapter available for " /*<<this->getName()*/
            << std::endl;
    bool needTimeDependence = false;
    for (int pos = 0; pos < nofAdapter; ++pos) {
        this->BCs[pos]->init(this, timeStep);
        if (this->BCs[pos]->isTimeDependent())
            needTimeDependence = true;
    }
    if (needTimeDependence)
        this->setTimeDependent();
    else
        this->unsetTimeDependent();

    //////////////////////////////////////////////////////////////////////////
    // grid.lock() info
    //////////////////////////////////////////////////////////////////////////
    int coarsestInitLevel = grid.lock()->getCoarsestInitializedLevel();
    int finestInitLevel   = grid.lock()->getFinestInitializedLevel();

    UbTupleInt3 blocknx = grid.lock()->getBlockNX();
    int blocknx1        = val<1>(blocknx); // applies to all levels
    int blocknx2        = val<2>(blocknx); // applies to all levels
    int blocknx3        = val<3>(blocknx); // applies to all levels

    // grobe Blocklaengen
    SPtr<CoordinateTransformation3D> trafo = grid.lock()->getCoordinateTransformator();
    double cblockDeltaX1, cblockDeltaX2, cblockDeltaX3, delta;
    cblockDeltaX1 = cblockDeltaX2 = cblockDeltaX3 = delta = c1o1 / (double)(1 << coarsestInitLevel);
    if (trafo) {
        cblockDeltaX1 = trafo->getX1CoordinateScaling() * delta;
        cblockDeltaX2 = trafo->getX2CoordinateScaling() * delta;
        cblockDeltaX3 = trafo->getX3CoordinateScaling() * delta;
    }
    // Level-specific block lengths and node spacing
    std::vector<std::vector<double>> nodeDeltaToNeigh(finestInitLevel + 1);
    std::vector<float> deltaMinX1(finestInitLevel + 1), deltaMinX2(finestInitLevel + 1),
        deltaMinX3(finestInitLevel + 1);
    std::vector<float> deltaMaxX1(finestInitLevel + 1), deltaMaxX2(finestInitLevel + 1),
        deltaMaxX3(finestInitLevel + 1);

    // In the Boltzman context, dx1==dx2==dx3 must be!!
    assert(UbMath::equal(cblockDeltaX1 / (double)blocknx1, cblockDeltaX2 / (double)blocknx2));
    assert(UbMath::equal(cblockDeltaX1 / (double)blocknx1, cblockDeltaX3 / (double)blocknx3));

    for (int level = coarsestInitLevel; level <= finestInitLevel; level++) {
        real nodeDeltaX1 = cblockDeltaX1 / (double)(blocknx1 * (1 << (level - coarsestInitLevel)));
        real nodeDeltaX2 = cblockDeltaX2 / (double)(blocknx2 * (1 << (level - coarsestInitLevel)));
        real nodeDeltaX3 = cblockDeltaX3 / (double)(blocknx3 * (1 << (level - coarsestInitLevel)));

        std::vector<real> distNeigh(D3Q27System::FENDDIR + 1, 0.0);
        D3Q27System::calcDistanceToNeighbors(distNeigh, nodeDeltaX1, nodeDeltaX2, nodeDeltaX3);

        nodeDeltaToNeigh[level].resize(D3Q27System::ENDDIR + 1, 0.0);
        for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
            nodeDeltaToNeigh[level][fdir] = distNeigh[fdir];
        }

        // TODO: set 5.0 as variable parameter in constructor, default 2.0
        deltaMinX1[level] = (float)(c5o1 * nodeDeltaX1); // no minus down there -deltaMin
        deltaMinX2[level] = (float)(c5o1 * nodeDeltaX2);
        deltaMinX3[level] = (float)(c5o1 * nodeDeltaX3);
        deltaMaxX1[level] = (float)(c5o1 * nodeDeltaX1);
        deltaMaxX2[level] = (float)(c5o1 * nodeDeltaX2);
        deltaMaxX3[level] = (float)(c5o1 * nodeDeltaX3);
    }

    //////////////////////////////////////////////////////////////////////////
    // TRIANGLES: q-determination
    //////////////////////////////////////////////////////////////////////////

    // initialize necessary variables (including blockDeltas of the coarse level)
    float triPoints[3][3];
    float vx1 = 0.0, vx2 = 0.0, vx3 = 0.0;
    std::vector<GbTriFaceMesh3D::TriFace> &triangles = *mesh->getTriangles();
    std::vector<GbTriFaceMesh3D::Vertex> &nodes      = *mesh->getNodes();
    std::map<SPtr<Block3D>, std::set<UbTupleInt3>> tmpSolidNodesFromOtherInteractors;

    int onePercent = UbMath::integerRounding(triangles.size() * c1o100);
    if (onePercent == 0)
        onePercent = 1;
    vf::basics::Timer setQTimer;
    setQTimer.start();
    UBLOG(logDEBUG3, " - setQs for " << (int)triangles.size() << " triangles");

    
    float blockMinX[3], blockMaxX[3], boxCenter[3], halfBoxSize[3];

    for (size_t t = 0; t < triangles.size(); t++) {
        //////////////////////////////////////////////////////////////////////////
        // Generate halfspace for the triangle and determine the min/max of the triangle
        //////////////////////////////////////////////////////////////////////////
        GbTriFaceMesh3D::TriFace &triangle = triangles[t];

        GbTriFaceMesh3D::Vertex &v1 = nodes[triangle.v1];
        GbTriFaceMesh3D::Vertex &v2 = nodes[triangle.v2];
        GbTriFaceMesh3D::Vertex &v3 = nodes[triangle.v3];

        if (this->isInverseSolid()) {
            triangle.nx *= (-1);
            triangle.ny *= (-1);
            triangle.nz *= (-1);
        }
        GbHalfSpace3D halfSpace(v1.x, v1.y, v1.z, triangle.nx, triangle.ny, triangle.nz);

        //////////////////////////////////////////////////////////////////////////
        // for GbMeshTools3D::triBoxOverlap
        //////////////////////////////////////////////////////////////////////////
        triPoints[0][0] = v1.x;
        triPoints[0][1] = v1.y;
        triPoints[0][2] = v1.z;
        triPoints[1][0] = v2.x;
        triPoints[1][1] = v2.y;
        triPoints[1][2] = v2.z;
        triPoints[2][0] = v3.x;
        triPoints[2][1] = v3.y;
        triPoints[2][2] = v3.z;

        double minX1 = triangle.getMinX(nodes);
        double maxX1 = triangle.getMaxX(nodes);
        double minX2 = triangle.getMinY(nodes);
        double maxX2 = triangle.getMaxY(nodes);
        double minX3 = triangle.getMinZ(nodes);
        double maxX3 = triangle.getMaxZ(nodes);

        //////////////////////////////////////////////////////////////////////////
        // Loop through all levels
        //////////////////////////////////////////////////////////////////////////
        double e1x1, e1x2, e1x3, e2x1, e2x2, e2x3, px1, px2, px3, a, f, sx1, sx2, sx3, u, qx1, qx2, qx3, v;
        bool gotQs = false;
        SPtr<BoundaryConditions> bc;

        for (int level = coarsestInitLevel; level <= finestInitLevel; level++) {
            //////////////////////////////////////////////////////////////////////////
            // Determine the level-specific BoundCube of the triangle and obtain the associated blocks
            //////////////////////////////////////////////////////////////////////////
            double boundCubeTriangleMinX1 = minX1 - deltaMinX1[level];
            double boundCubeTriangleMaxX1 = maxX1 + deltaMaxX1[level];
            double boundCubeTriangleMinX2 = minX2 - deltaMinX2[level];
            double boundCubeTriangleMaxX2 = maxX2 + deltaMaxX2[level];
            double boundCubeTriangleMinX3 = minX3 - deltaMinX3[level];
            double boundCubeTriangleMaxX3 = maxX3 + deltaMaxX3[level];

            GbCuboid3D boundingCubeTriangle(boundCubeTriangleMinX1, boundCubeTriangleMinX2, boundCubeTriangleMinX3,
                                            boundCubeTriangleMaxX1, boundCubeTriangleMaxX2, boundCubeTriangleMaxX3);

            std::vector<SPtr<Block3D>> triBlocks;
            grid.lock()->getBlocksByCuboid(level, boundCubeTriangleMinX1, boundCubeTriangleMinX2,
                                           boundCubeTriangleMinX3, boundCubeTriangleMaxX1, boundCubeTriangleMaxX2,
                                           boundCubeTriangleMaxX3, triBlocks);

            //////////////////////////////////////////////////////////////////////////
            // Loop over blocks of the level that contain the triangle
            //////////////////////////////////////////////////////////////////////////
            for (std::size_t b = 0; b < triBlocks.size(); b++) {
                SPtr<Block3D> block = triBlocks[b];

                ////////////////////////////////////////////////////////////////////////////
                //// Block triangle/test
                ////////////////////////////////////////////////////////////////////////////
                UbTupleDouble3 coords = grid.lock()->getBlockWorldCoordinates(block);
                UbTupleDouble3 deltas = grid.lock()->getBlockLengths(block);

                blockMinX[0] = (float)(val<1>(coords) - deltaMinX1[level]);
                blockMinX[1] = (float)(val<2>(coords) - deltaMinX2[level]);
                blockMinX[2] = (float)(val<3>(coords) - deltaMinX3[level]);

                blockMaxX[0] = (float)(val<1>(coords) + val<1>(deltas) + deltaMaxX1[level]);
                blockMaxX[1] = (float)(val<2>(coords) + val<2>(deltas) + deltaMaxX2[level]);
                blockMaxX[2] = (float)(val<3>(coords) + val<3>(deltas) + deltaMaxX3[level]);

                boxCenter[0] = (float)(c1o2 * (blockMaxX[0] + blockMinX[0]));
                boxCenter[1] = (float)(c1o2 * (blockMaxX[1] + blockMinX[1]));
                boxCenter[2] = (float)(c1o2 * (blockMaxX[2] + blockMinX[2]));

                halfBoxSize[0] = (float)(c1o2 * (blockMaxX[0] - blockMinX[0]));
                halfBoxSize[1] = (float)(c1o2 * (blockMaxX[1] - blockMinX[1]));
                halfBoxSize[2] = (float)(c1o2 * (blockMaxX[2] - blockMinX[2]));

                // if triangle "enlarged cube" does not intersect/touch -> no BC possible -> continue
                if (!GbMeshTools3D::triBoxOverlap(boxCenter, halfBoxSize, triPoints)) {
                    continue;
                }

                //////////////////////////////////////////////////////////////////////////
                // Examination of the individual nodes
                //////////////////////////////////////////////////////////////////////////
                bool blockGotBCs = false;

                SPtr<ILBMKernel> kernel  = block->getKernel();
                SPtr<BCArray3D> bcMatrix = kernel->getBCSet()->getBCArray();

                int indexMinX1 = 0;
                int indexMinX2 = 0;
                int indexMinX3 = 0;

                int indexMaxX1 = (int)bcMatrix->getNX1();
                int indexMaxX2 = (int)bcMatrix->getNX2();
                int indexMaxX3 = (int)bcMatrix->getNX3();

                std::set<std::vector<int>> &bcNodeIndices = this->bcNodeIndicesMap[block];
                double q, distance;

                double &nodeDx1 = nodeDeltaToNeigh[level][dP00];
                double &nodeDx2 = nodeDeltaToNeigh[level][d0P0];
                double &nodeDx3 = nodeDeltaToNeigh[level][d00P];

                // fuer OBB-Test
                double qEinflussDelta = 1.1 * sqrt(nodeDx1 * nodeDx1 + nodeDx2 * nodeDx2 + nodeDx3 * nodeDx3);

                for (int ix3 = indexMinX3; ix3 < indexMaxX3; ix3++) {
                    for (int ix2 = indexMinX2; ix2 < indexMaxX2; ix2++) {
                        for (int ix1 = indexMinX1; ix1 < indexMaxX1; ix1++) {
                            Vector3D pointplane1 = grid.lock()->getNodeCoordinates(block, ix1, ix2, ix3);
                            double internX1      = pointplane1[0];
                            double internX2      = pointplane1[1];
                            double internX3      = pointplane1[2];

                            if (bcMatrix->isSolid(ix1, ix2, ix3) || bcMatrix->isUndefined(ix1, ix2, ix3)) {
                                continue;
                            }

                            //////////////////////////////////////////////////////////////////////////
                            // Point in AABB of triangle?
                            //////////////////////////////////////////////////////////////////////////
                            // Ehsan changed
                            bool pointIsOnBoundary = true;
                            if (!boundingCubeTriangle.isPointInGbObject3D(internX1, internX2, internX3,
                                                                          pointIsOnBoundary)) {
                                continue;
                            }
                            // std::cout<<"internX3  "<<internX3<<"  internX2"<<internX2<<" internX1 "<<internX1<<"\n";
                            //////////////////////////////////////////////////////////////////////////
                            // Half-plane tests
                            //////////////////////////////////////////////////////////////////////////
                            distance = halfSpace.getDistance(internX1, internX2, internX3);
                            // Point in half plane? (no, if distance<0)
                            if (useHalfSpace &&
                                UbMath::less(distance, 0.0)) //== !halfSpace.ptInside(internX1,internX2,internX3) )
                            {
                                continue;
                            }

                            // CheapOBB test: if distance > qInfluenceDelta -> no q
                            if (UbMath::greater(fabs(distance), qEinflussDelta)) {
                                continue;
                            }

                            /////////////////////////////////////////////////////////////////////////////
                            // Ray tracing for discrete Boltzmann directions
                            /////////////////////////////////////////////////////////////////////////////
                            gotQs = false;
                            bc    = SPtr<BoundaryConditions>();

                            // RAYTRACING - diskrete LB-dir zu Dreick
                            // e1 = v1 - v0
                            e1x1 = v2.x - v1.x;
                            e1x2 = v2.y - v1.y;
                            e1x3 = v2.z - v1.z;

                            // e2 = v2 - v0
                            e2x1 = v3.x - v1.x;
                            e2x2 = v3.y - v1.y;
                            e2x3 = v3.z - v1.z;

                            // s = o - v0
                            sx1 = internX1 - v1.x;
                            sx2 = internX2 - v1.y;
                            sx3 = internX3 - v1.z;

                            for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
                                // p = d x e2
                                px1 = this->rayX2[fdir] * e2x3 - this->rayX3[fdir] * e2x2;
                                px2 = this->rayX3[fdir] * e2x1 - this->rayX1[fdir] * e2x3;
                                px3 = this->rayX1[fdir] * e2x2 - this->rayX2[fdir] * e2x1;

                                // a = e1 dot p
                                a = e1x1 * px1 + e1x2 * px2 + e1x3 * px3;
                                if (fabs(a) < 1.E-10)
                                    continue;
                                f = c1o1 / a;

                                // u = f * ( s dot p)
                                u = f * (sx1 * px1 + sx2 * px2 + sx3 * px3);
                                if (u < -1.E-10 || u > c1o1 + 1.E-10)
                                    continue;

                                // q = s x e1
                                qx1 = sx2 * e1x3 - sx3 * e1x2;
                                qx2 = sx3 * e1x1 - sx1 * e1x3;
                                qx3 = sx1 * e1x2 - sx2 * e1x1;

                                // v = f*(e2 dot q)
                                v = f * (this->rayX1[fdir] * qx1 + this->rayX2[fdir] * qx2 + this->rayX3[fdir] * qx3);
                                if (v < -1.E-10 || (u + v) > c1o1 + 1.E-10)
                                    continue;

                                // t = f * (e2 dot q)
                                q = f * (e2x1 * qx1 + e2x2 * qx2 + e2x3 * qx3);
                                q /= nodeDeltaToNeigh[level][fdir];
                                /////ehsan q/////////////////////////////////////////////////////////////////////
                                double det = triangle.nx * this->rayX1[fdir] + triangle.ny * this->rayX2[fdir] +
                                             triangle.nz * this->rayX3[fdir];

                                if (det > -1.E-10)
                                    continue;
                                double d  = triangle.nx * v1.x + triangle.ny * v1.y + triangle.nz * v1.z;
                                double x1 = -((-d * this->rayX1[fdir] - triangle.ny * this->rayX2[fdir] * internX1 -
                                               triangle.nz * this->rayX3[fdir] * internX1 +
                                               triangle.ny * this->rayX1[fdir] * internX2 +
                                               triangle.nz * this->rayX1[fdir] * internX3)) /
                                            det;
                                double y1 = -((-d * this->rayX2[fdir] + triangle.nx * this->rayX2[fdir] * internX1 -
                                               triangle.nx * this->rayX1[fdir] * internX2 -
                                               triangle.nz * this->rayX3[fdir] * internX2 +
                                               triangle.nz * this->rayX2[fdir] * internX3)) /
                                            det;
                                double z1 = -((-d * this->rayX3[fdir] + triangle.nx * this->rayX3[fdir] * internX1 +
                                               triangle.ny * this->rayX3[fdir] * internX2 -
                                               triangle.nx * this->rayX1[fdir] * internX3 -
                                               triangle.ny * this->rayX2[fdir] * internX3)) /
                                            det;
                                double q_ehsan =
                                    sqrt((x1 - internX1) * (x1 - internX1) + (y1 - internX2) * (y1 - internX2) +
                                         (z1 - internX3) * (z1 - internX3));
                                q_ehsan /= nodeDeltaToNeigh[level][fdir];
                                q = q_ehsan;
                                if (UbMath::greater(q, c1o1) || UbMath::lessEqual(q, 0.0))
                                    continue;

                                //Check found q for validity
                                if (UbMath::zero(q)) {
                                    // new (05/18/2010)
                                    // It can happen that with thin-walled geos there are points on a triangle
                                    // lie, get qs that go through the geo. these points will come later
                                    // but no longer tested for solidity because they didn't get a BC
                                    //--> there is at least a q==0.0 for one of the triangles -> put it solid there
                                    this->solidNodeIndicesMap[block].insert(UbTupleInt3(ix1, ix2, ix3));
                                    bcMatrix->setSolid(ix1, ix2, ix3);
                                    continue;
                                }

                                if (UbMath::inClosedInterval(q, c1o1, c1o1))
                                    q = c1o1;
                                if (UbMath::greater(q, 0.0) && UbMath::lessEqual(q, c1o1)) {
                                    gotQs = blockGotBCs = true;

                                    bc = bcMatrix->getBC(ix1, ix2, ix3);

                                    // SG 26.08.2010 if(!bc && !bcMatrix->isSolid())
                                    if (!bc) {
                                        bc = SPtr<BoundaryConditions>(new BoundaryConditions);
                                        ;
                                        bcMatrix->setBC(ix1, ix2, ix3, bc);
                                    } else if (UbMath::less(bc->getQ(fdir), q) &&
                                               UbMath::equal(-999.0, q)) // schon ein kuerzeres q voehanden?
                                    {
                                        // new:: May 18, 2010
                                        // to avoid incorrect qs that might go through the "wall".
                                        // only reset q if new q is smaller than existing one!
                                        // Also: especially on corners with two BC geos, that's all
                                        // more details valid
                                        continue;
                                    }

                                    bc->setBoundaryVelocityX1(vx1);
                                    bc->setBoundaryVelocityX2(vx2);
                                    bc->setBoundaryVelocityX3(vx3);

                                    for (int index = (int)this->BCs.size() - 1; index >= 0; --index)
                                        this->BCs[index]->adaptBCForDirection(*this, bc, internX1, internX2,
                                                                                     internX3, q, fdir);

                                    // for accelerated rereading
                                    if (this->reinitWithStoredQsFlag) {
                                        bcNodeIndicesAndQsMap[block][UbTupleInt3(ix1, ix2, ix3)].resize(
                                            D3Q27System::FENDDIR + 1 + 3, -1.0f);
                                        bcNodeIndicesAndQsMap[block][UbTupleInt3(ix1, ix2, ix3)][fdir] = float(q);
                                        bcNodeIndicesAndQsMap[block][UbTupleInt3(ix1, ix2, ix3)]
                                                             [D3Q27System::FENDDIR + 1 + 0] = float(internX1);
                                        bcNodeIndicesAndQsMap[block][UbTupleInt3(ix1, ix2, ix3)]
                                                             [D3Q27System::FENDDIR + 1 + 1] = float(internX2);
                                        bcNodeIndicesAndQsMap[block][UbTupleInt3(ix1, ix2, ix3)]
                                                             [D3Q27System::FENDDIR + 1 + 2] = float(internX3);
                                    }
                                }
                            }

                            if (gotQs) {
                                std::vector<int> p(3);
                                p[0] = ix1;
                                p[1] = ix2;
                                p[2] = ix3;
                                bcNodeIndices.insert(p);

                                for (int index = (int)this->BCs.size() - 1; index >= 0; --index)
                                    this->BCs[index]->adaptBC(*this, bc, internX1, internX2, internX3);
                            }
                        }
                    }
                }
            }
            // Unfortunately, dynamic points of the GbCuboids have to be deleted manually :-(
            boundingCubeTriangle.finalize();
        }
    }
    UBLOGML(logDEBUG1, "\nLBMTriFaceMeshInteractor - setQs end ");
}
//////////////////////////////////////////////////////////////////////////
// Procedure
// A – Determination of the q's
// 1. The blocks of the block grid are determined for each bounding cube of a triangle in the network
// 2. using a triangle/block intersection test, further irrelevant blocks are sorted out
// (for long triangles that are crooked in space and where the bounding cube is suboptimal)
// 3. Each node of these blocks is tested against the bound cube of the triangle
// 4. Nodes that are within the cube but “within” the network are sorted out using a half-plane test
// 5. For the remaining nodes, the q determination is carried out using efficient ray tracing algorithms
// for the discrete Boltzmann directions
// B – Setting the inactive blocks and solid nodes
// all blocks of the bounding cube of the network that received at least one BC were marked in A
// 1. for unmarked blocks, ONE pointInObject (triangular network) test is sufficient for the entire block if successful
// mark “not active”.
// 2. a recursive filling algorithm is carried out for marked blocks
//////////////////////////////////////////////////////////////////////////
void D3Q27TriFaceMeshInteractor::refineBlockGridToLevel(int level, real startDistance, real stopDistance)
{
    using namespace vf::basics::constant;

    UBLOG(logDEBUG1, "D3Q27TriFaceMeshInteractor::refineBlockGridToLevel - start");

    // ToDo: maybe check whether you can install a half-space check for StopDistance
    // or whether you're faster if you don't do a half-space test at all...
    if (!grid.lock())
        throw UbException(UB_EXARGS, "Grid isn't exist!");
    if (UbMath::greater(startDistance, 0.0))
        throw UbException(UB_EXARGS, "startDistance>0.0 not supported by this interactor");
    if (UbMath::less(stopDistance, 0.0))
        throw UbException(UB_EXARGS, "stopDistance<0.0  not supported by this interactor");

    SPtr<Grid3D> bgrid    = this->grid.lock();
    GbTriFaceMesh3D &mesh = dynamic_cast<GbTriFaceMesh3D &>(*this->geoObject3D.get());

    int coarsestLevel = bgrid->getCoarsestInitializedLevel();

    std::vector<GbTriFaceMesh3D::TriFace> &triangles = *mesh.getTriangles();
    std::vector<GbTriFaceMesh3D::Vertex> &nodes      = *mesh.getNodes();

    double minX1, minX2, minX3, maxX1, maxX2, maxX3;
    float blockMinX[3], blockMaxX[3], boxCenter[3], halfBoxSize[3];
    float triPoints[3][3];

    size_t nofTriangles = (int)triangles.size();

    //#pragma omp parallel
    //#pragma omp for
    for (size_t i = 0; i < nofTriangles; i++)
    // for(int i=0; i<nofTriangles; i++)
    {
        //#pragma omp master
        //{
        //    printf_s("num_threads=%d\n", omp_get_num_threads( ));
        //}

        GbTriFaceMesh3D::TriFace &triangle = triangles[i];

        GbTriFaceMesh3D::Vertex &v1 = nodes[triangle.v1];
        GbTriFaceMesh3D::Vertex &v2 = nodes[triangle.v2];
        GbTriFaceMesh3D::Vertex &v3 = nodes[triangle.v3];

        // dreick muss normal besitzen!
        assert(!UbMath::zero(triangle.nx) || !UbMath::zero(triangle.ny) || !UbMath::zero(triangle.nz));
        // Normale muss normiert sein!
        assert((fabs(std::sqrt(triangle.nx * triangle.nx + triangle.ny * triangle.ny + triangle.nz * triangle.nz)) -
                1.0f) < 1.0E-6);

        // Move halfspace around startDistance to normal, otherwise we will do it later
        // blocks to be tested on the back of the triangle not tested!!!
        GbHalfSpace3D halfSpace(
            v1.x + startDistance * triangle.nx, v1.y + startDistance * triangle.ny, v1.z + startDistance * triangle.nz,
            v2.x + startDistance * triangle.nx, v2.y + startDistance * triangle.ny, v2.z + startDistance * triangle.nz,
            v3.x + startDistance * triangle.nx, v3.y + startDistance * triangle.ny, v3.z + startDistance * triangle.nz);

        // Expand the bounding box with relevant dx -> to determine the blocks to be tested
        if (triangle.nx > 1.0E-8) {
            minX1 = triangle.getMinX(nodes) + c21o20 * triangle.nx * startDistance;
            maxX1 = triangle.getMaxX(nodes) + c21o20 * triangle.nx * stopDistance;
        } else {
            minX1 = triangle.getMinX(nodes) + c21o20 * triangle.nx * stopDistance;
            maxX1 = triangle.getMaxX(nodes) + c21o20 * triangle.nx * startDistance;
        }

        if (triangle.ny > 1.0E-8) {
            minX2 = triangle.getMinY(nodes) + c21o20 * triangle.ny * startDistance;
            maxX2 = triangle.getMaxY(nodes) + c21o20 * triangle.ny * stopDistance;
        } else {
            minX2 = triangle.getMinY(nodes) + c21o20 * triangle.ny * stopDistance;
            maxX2 = triangle.getMaxY(nodes) + c21o20 * triangle.ny * startDistance;
        }

        if (triangle.nz > 1.0E-8) {
            minX3 = triangle.getMinZ(nodes) + c21o20 * triangle.nz * startDistance;
            maxX3 = triangle.getMaxZ(nodes) + c21o20 * triangle.nz * stopDistance;
        } else {
            minX3 = triangle.getMinZ(nodes) + c21o20 * triangle.nz * stopDistance;
            maxX3 = triangle.getMaxZ(nodes) + c21o20 * triangle.nz * startDistance;
        }

        int flag = 0;
        // Get all the blocks that intersect extended BB level by level
        // and edit
        for (int l = coarsestLevel; l < level; l++) {
            std::vector<SPtr<Block3D>> consideredBlocks;
            bgrid->getBlocksByCuboid(l, minX1, minX2, minX3, maxX1, maxX2, maxX3, consideredBlocks);
            double x1a, x2a, x3a, x1b, x2b, x3b;

            for (size_t b = 0; b < consideredBlocks.size(); b++) {
                SPtr<Block3D> block = consideredBlocks[b];
                if (block->getLevel() >= level)
                    continue;

                // Determine the start coordinates of the block
                UbTupleDouble3 coords = bgrid->getBlockWorldCoordinates(block);
                UbTupleDouble3 deltas = bgrid->getBlockLengths(block);

                // Check whether the block is completely in the half space
                x1a = val<1>(coords);
                x1b = val<1>(coords) + val<1>(deltas);
                x2a = val<2>(coords);
                x2b = val<2>(coords) + val<2>(deltas);
                x3a = val<3>(coords);
                x3b = val<3>(coords) + val<3>(deltas);

                flag = 0;
                if (!halfSpace.ptInside(x1a, x2a, x3a))
                    flag |= (1 << 0); // 1
                if (!halfSpace.ptInside(x1b, x2a, x3a))
                    flag |= (1 << 1); // 2
                if (!halfSpace.ptInside(x1b, x2b, x3a))
                    flag |= (1 << 2); // 4
                if (!halfSpace.ptInside(x1a, x2b, x3a))
                    flag |= (1 << 3); // 8
                if (!halfSpace.ptInside(x1a, x2a, x3b))
                    flag |= (1 << 4); // 16
                if (!halfSpace.ptInside(x1b, x2a, x3b))
                    flag |= (1 << 5); // 32
                if (!halfSpace.ptInside(x1b, x2b, x3b))
                    flag |= (1 << 6); // 64
                if (!halfSpace.ptInside(x1a, x2b, x3b))
                    flag |= (1 << 7); // 128

                if (true && flag != 255) {
                    // determine block side (scalar product triangle-normal, vector (midTri->midCub) )
                    // depending on this, start or stopdistance must be used for the relevant block
                    // block is on the pos side -> stopdistance otherwise startdistance
                    double skalarprod = triangle.nx * (c1o2 * (x1a + x1b) - triangle.getX1Centroid(nodes)) +
                                        triangle.ny * (c1o2 * (x2a + x2b) - triangle.getX2Centroid(nodes)) +
                                        triangle.nz * (c1o2 * (x3a + x3b) - triangle.getX3Centroid(nodes));

                    double blockdelta = c21o20 * stopDistance;
                    if (skalarprod < 1.E-8)
                        blockdelta = -c21o20 * startDistance; // startDistance<0!!
                    else if (fabs(skalarprod) < 1.E-8)
                        blockdelta = c21o20 * UbMath::max(-startDistance, stopDistance);

                    // adjust block
                    blockMinX[0] = (float)(val<1>(coords) - blockdelta);
                    blockMinX[1] = (float)(val<2>(coords) - blockdelta);
                    blockMinX[2] = (float)(val<3>(coords) - blockdelta);

                    blockMaxX[0] = (float)(val<1>(coords) + val<1>(deltas) + blockdelta);
                    blockMaxX[1] = (float)(val<2>(coords) + val<2>(deltas) + blockdelta);
                    blockMaxX[2] = (float)(val<3>(coords) + val<3>(deltas) + blockdelta);

                    boxCenter[0] = (float)(c1o2 * (blockMaxX[0] + blockMinX[0]));
                    boxCenter[1] = (float)(c1o2 * (blockMaxX[1] + blockMinX[1]));
                    boxCenter[2] = (float)(c1o2 * (blockMaxX[2] + blockMinX[2]));

                    halfBoxSize[0] = (float)(c1o2 * (blockMaxX[0] - blockMinX[0]));
                    halfBoxSize[1] = (float)(c1o2 * (blockMaxX[1] - blockMinX[1]));
                    halfBoxSize[2] = (float)(c1o2 * (blockMaxX[2] - blockMinX[2]));

                    GbTriFaceMesh3D::Vertex &v1_ = nodes[triangle.v1];
                    GbTriFaceMesh3D::Vertex &v2_ = nodes[triangle.v2];
                    GbTriFaceMesh3D::Vertex &v3_ = nodes[triangle.v3];

                    triPoints[0][0] = v1_.x;
                    triPoints[0][1] = v1_.y;
                    triPoints[0][2] = v1_.z;
                    triPoints[1][0] = v2_.x;
                    triPoints[1][1] = v2_.y;
                    triPoints[1][2] = v2_.z;
                    triPoints[2][0] = v3_.x;
                    triPoints[2][1] = v3_.y;
                    triPoints[2][2] = v3_.z;

                    // if block triangle cuts, then it needs to be refined
                    if (GbMeshTools3D::triBoxOverlap(boxCenter, halfBoxSize, triPoints)) {
                        bgrid->expandBlock(block->getX1(), block->getX2(), block->getX3(), block->getLevel());
                    }
                }
            }
        }
    }
    UBLOG(logDEBUG1, " - refine done");
}
//////////////////////////////////////////////////////////////////////////
void D3Q27TriFaceMeshInteractor::updateMovedGeometry(const real &timeStep) {}
////////////////////////////////////////////////////////////////////////////
void D3Q27TriFaceMeshInteractor::recursiveGridFill(CbArray3D<FLAGS> &flagfield, const short &xs, const short &ys,
                                                   const short &zs, const FLAGS &type)
{
    // Algorithmus zum Füllen eines Polyeders, ausgehend vom Saatpunkt xs,ys,zs

    // Saatknoten einfärben
    if (flagfield(xs, ys, zs) == UNDEF_FLAG) {
        flagfield(xs, ys, zs) = type;

        if (flagfield.indicesInRange(xs + 1, ys, zs))
            this->recursiveGridFill(flagfield, xs + 1, ys, zs, type);
        if (flagfield.indicesInRange(xs, ys + 1, zs))
            this->recursiveGridFill(flagfield, xs, ys + 1, zs, type);
        if (flagfield.indicesInRange(xs, ys, zs + 1))
            this->recursiveGridFill(flagfield, xs, ys, zs + 1, type);
        if (flagfield.indicesInRange(xs - 1, ys, zs))
            this->recursiveGridFill(flagfield, xs - 1, ys, zs, type);
        if (flagfield.indicesInRange(xs, ys - 1, zs))
            this->recursiveGridFill(flagfield, xs, ys - 1, zs, type);
        if (flagfield.indicesInRange(xs, ys, zs - 1))
            this->recursiveGridFill(flagfield, xs, ys, zs - 1, type);
    }
}
//////////////////////////////////////////////////////////////////////////
void D3Q27TriFaceMeshInteractor::iterativeGridFill(CbArray3D<FLAGS> &flagfield, const short &xs, const short &ys,
                                                   const short &zs, const FLAGS &type)
{
    std::stack<UbTupleInt3> stck;
    stck.push(UbTupleInt3(xs, ys, zs));

    int x, y, z;

    while (!stck.empty()) {
        x = val<1>(stck.top());
        y = val<2>(stck.top());
        z = val<3>(stck.top());
        stck.pop();

        FLAGS &flagType = flagfield(x, y, z);

        if (flagType == UNDEF_FLAG) {
            flagType = type;

            if (flagfield.indicesInRange(x + 1, y, z))
                stck.push(UbTupleInt3(x + 1, y, z));
            if (flagfield.indicesInRange(x, y + 1, z))
                stck.push(UbTupleInt3(x, y + 1, z));
            if (flagfield.indicesInRange(x, y, z + 1))
                stck.push(UbTupleInt3(x, y, z + 1));
            if (flagfield.indicesInRange(x - 1, y, z))
                stck.push(UbTupleInt3(x - 1, y, z));
            if (flagfield.indicesInRange(x, y - 1, z))
                stck.push(UbTupleInt3(x, y - 1, z));
            if (flagfield.indicesInRange(x, y, z - 1))
                stck.push(UbTupleInt3(x, y, z - 1));
        }
    }
    return;
}
//////////////////////////////////////////////////////////////////////////
string D3Q27TriFaceMeshInteractor::toString()
{
    stringstream ss;
    ss << "D3Q27TriFaceMeshInteractor[label=D3Q27TriFaceMeshInteractor";
    if (this->isSolid())
        ss << ", solid";
    if (this->isInverseSolid())
        ss << ", inversesolid";
    if (this->isTimeDependent())
        ss << ", timedependent";
    if (geoObject3D != NULL)
        ss << ", D3Q27TriFaceMeshInteractor: " << geoObject3D->toString();
    ss << "]";

    return ss.str();
}
//////////////////////////////////////////////////////////////////////////
void D3Q27TriFaceMeshInteractor::reinitWithStoredQs(const real & /*timeStep*/)
{
    using namespace vf::basics::constant;

    // alle solid Bloecke wieder solid setzen
    std::vector<SPtr<Block3D>> &solidBlocks = this->getSolidBlockSet();
    for (size_t i = 0; i < solidBlocks.size(); i++) {
        solidBlocks[i]->setActive(false); //<- quick n dirty
    }

    // alle solid-nodes wieder solid setzen (solids die quasi in den TransBloecken liegen)
    std::map<SPtr<Block3D>, std::set<UbTupleInt3>>::iterator it1;
    for (it1 = this->solidNodeIndicesMap.begin(); it1 != this->solidNodeIndicesMap.end(); ++it1) {
        SPtr<Block3D> block = it1->first;

        SPtr<ILBMKernel> kernel           = block->getKernel();
        SPtr<BCArray3D> bcMatrix          = kernel->getBCSet()->getBCArray();
        std::set<UbTupleInt3> &indicesSet = it1->second;

        for (std::set<UbTupleInt3>::iterator setIt = indicesSet.begin(); setIt != indicesSet.end(); ++setIt) {
            bcMatrix->setSolid(val<1>(*setIt), val<2>(*setIt), val<3>(*setIt));
        }
    }

    // BCS WIEDERHERSTELLEN
    std::map<SPtr<Block3D>, std::map<UbTupleInt3, std::vector<float>>>::iterator it;
    for (it = bcNodeIndicesAndQsMap.begin(); it != bcNodeIndicesAndQsMap.end(); ++it) {
        SPtr<Block3D> block      = it->first;
        SPtr<ILBMKernel> kernel  = block->getKernel();
        SPtr<BCArray3D> bcMatrix = kernel->getBCSet()->getBCArray();

        std::map<UbTupleInt3, std::vector<float>>::iterator it2;
        for (it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
            const UbTupleInt3 &pos = it2->first;
            std::vector<float> qs  = it2->second;

            // SG_27.08.2010
            if (bcMatrix->isSolid(val<1>(pos), val<2>(pos), val<3>(pos)))
                continue;

            SPtr<BoundaryConditions> bc = bcMatrix->getBC(val<1>(pos), val<2>(pos), val<3>(pos));
            if (!bc) {
                bc = SPtr<BoundaryConditions>(new BoundaryConditions);
                bcMatrix->setBC(val<1>(pos), val<2>(pos), val<3>(pos), bc);
            }

            double x1w = qs[D3Q27System::FENDDIR + 1 + 0];
            double x2w = qs[D3Q27System::FENDDIR + 1 + 1];
            double x3w = qs[D3Q27System::FENDDIR + 1 + 2];

            // TODO: HACK DOES NOT BELONG HERE!!! - begin
            // it is a static object and the propeller is there
            // there are difficulties at the wing tips, that can come from
            // that there are too many bc flags and at speed this results in quite large values
            bc->setBoundaryVelocityX1(c0o1);
            bc->setBoundaryVelocityX2(c0o1);
            bc->setBoundaryVelocityX3(c0o1);
            // TODO: HACK DOES NOT BELONG HERE!!! - end

            bool gotQs = false;
            for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
                if (UbMath::greater(qs[fdir], -1.0) && UbMath::less(qs[fdir], bc->getQ(fdir))) {
                    gotQs = true;
                    for (size_t index = 0; index < this->BCs.size(); index++)
                        this->BCs[index]->adaptBCForDirection(*this, bc, x1w, x2w, x3w, qs[fdir], fdir);
                }
            }

            if (gotQs)
                for (size_t index = 0; index < this->BCs.size(); index++)
                    this->BCs[index]->adaptBC(*this, bc, x1w, x2w, x3w);
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void D3Q27TriFaceMeshInteractor::updateInteractor(const real &timestep)
{
    D3Q27Interactor::updateInteractor(timestep);
}

//! \}
