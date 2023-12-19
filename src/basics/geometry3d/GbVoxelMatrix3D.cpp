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
//! \addtogroup geometry3d
//! \ingroup basics
//! \{
//! \author Konstantin Kutscher, Soeren Textor, Sebastian Geller
//=======================================================================================
#include <geometry3d/GbVoxelMatrix3D.h>

#include <basics/utilities/UbFileInputASCII.h>
#include <basics/utilities/UbFileOutputASCII.h>
#include <basics/utilities/UbMath.h>
#include <geometry3d/CoordinateTransformation3D.h>
#include <geometry3d/GbTriangle3D.h>

#include <basics/utilities/UbSystem.h>
#include "basics/constants/NumericConstants.h"

using namespace std;

const float GbVoxelMatrix3D::SOLID = 1.0f;
const float GbVoxelMatrix3D::FLUID = 0.0f;

/*=======================================================*/
// Konstruktor
GbVoxelMatrix3D::GbVoxelMatrix3D(int nx1, int nx2, int nx3, float initVal, double lowerThreshold, double upperThreshold)
    : GbObject3D(), lowerThreshold(lowerThreshold), upperThreshold(upperThreshold), nodesX1(nx1), nodesX2(nx2),
      nodesX3(nx3), voxelMatrix(Matrix3D(nx1, nx2, nx3, initVal))
{
    this->setName("VoxelMatrix3D");
}
/*=======================================================*/
GbVoxelMatrix3D::GbVoxelMatrix3D() : GbObject3D() { this->setName("VoxelMatrix3D"); }
/*=======================================================*/
GbVoxelMatrix3D::GbVoxelMatrix3D(const Matrix3D &voxelMatrix, double lowerThreshold, double upperThreshold)
    : GbObject3D(), nodesX1((int)voxelMatrix.getNX1()), nodesX2((int)voxelMatrix.getNX2()),
      nodesX3((int)voxelMatrix.getNX3()), lowerThreshold(lowerThreshold), upperThreshold(upperThreshold),
      voxelMatrix(voxelMatrix)
{
    this->setName("VoxelMatrix3D");
}
/*=======================================================*/
GbVoxelMatrix3D *GbVoxelMatrix3D::clone()
{
    GbVoxelMatrix3D *vm = new GbVoxelMatrix3D(voxelMatrix, lowerThreshold, upperThreshold);
    vm->setVoxelMatrixMininum(minX1, minX2, minX3);
    vm->setVoxelMatrixDelta(deltaX1, deltaX2, deltaX3);
    return vm;
}
/*=======================================================*/
void GbVoxelMatrix3D::setCenterCoordinates(const double &x1, const double &x2, const double &x3)
{
    this->translate(x1 - getX1Centroid(), x2 - getX2Centroid(), x3 - getX3Centroid());
}
/*=======================================================*/
void GbVoxelMatrix3D::translate(const double &tx1, const double &tx2, const double &tx3)
{
    this->minX1 += tx1;
    this->minX2 += tx2;
    this->minX3 += tx3;
    this->notifyObserversObjectChanged();
}
/*=======================================================*/
void GbVoxelMatrix3D::setClosedVoidSpaceToSolid()
{
    voxelMatrixTemp = Matrix3D(nodesX1, nodesX2, nodesX3, SOLID);
    flagMatrix      = CbArray3D<char>(nodesX1, nodesX2, nodesX3, 0);

    for (int x3 = 0; x3 < nodesX3; x3++) {
        for (int x2 = 0; x2 < nodesX2; x2++) {
            for (int x1 = 0; x1 < nodesX1; x1++) {
                if (voxelMatrix(x1, x2, x3) == FLUID) {
                    UBLOG(logINFO, "setClosedVoidSpaceToSolid:start");
                    x1Nbr.push_back(x1);
                    x2Nbr.push_back(x2);
                    x3Nbr.push_back(x3);
                    int size = (int)x1Nbr.size();
                    while (size > 0) {
                        for (int i = 0; i < size; i++) {
                            findFluidNeighbor(x1Nbr[i], x2Nbr[i], x3Nbr[i]);
                        }

                        swap(x1Nbr, x1NbrTemp);
                        swap(x2Nbr, x2NbrTemp);
                        swap(x3Nbr, x3NbrTemp);

                        x1NbrTemp.clear();
                        x2NbrTemp.clear();
                        x3NbrTemp.clear();
                        size = (int)x1Nbr.size();
                    }
                    UBLOG(logINFO, "setClosedVoidSpaceToSolid:end");
                    voxelMatrix = voxelMatrixTemp;
                    return;
                }
            }
        }
    }
}
/*=======================================================*/
void GbVoxelMatrix3D::findFluidNeighbor(int x1, int x2, int x3)
{
    for (int k3 = -1; k3 <= 1; k3++) {
        for (int k2 = -1; k2 <= 1; k2++) {
            for (int k1 = -1; k1 <= 1; k1++) {
                int j1 = x1 + k1;
                int j2 = x2 + k2;
                int j3 = x3 + k3;
                if (j1 >= 0 && j1 < nodesX1 && j2 >= 0 && j2 < nodesX2 && j3 >= 0 && j3 < nodesX3) {
                    if (voxelMatrix(j1, j2, j3) == FLUID) {
                        if (flagMatrix(j1, j2, j3) == 0) {
                            voxelMatrixTemp(j1, j2, j3) = FLUID;
                            flagMatrix(j1, j2, j3)      = 1;
                            x1NbrTemp.push_back(j1);
                            x2NbrTemp.push_back(j2);
                            x3NbrTemp.push_back(j3);
                        }
                    }
                }
            }
        }
    }
}
/*=======================================================*/
void GbVoxelMatrix3D::calculateNumberOfSolidAndFluid()
{
    numberOfSolid = 0;

    for (int x3 = 0; x3 < nodesX3; x3++)
        for (int x2 = 0; x2 < nodesX2; x2++)
            for (int x1 = 0; x1 < nodesX1; x1++) {
                if (voxelMatrix(x1, x2, x3) == GbVoxelMatrix3D::SOLID) {
                    numberOfSolid++;
                }
            }

    numberOfFluid = (long)nodesX1 * (long)nodesX2 * (long)nodesX3 - numberOfSolid;
}
/*=======================================================*/
long GbVoxelMatrix3D::getNumberOfSolid() { return numberOfSolid; }
/*=======================================================*/
long GbVoxelMatrix3D::getNumberOfFluid() { return numberOfFluid; }
/*=======================================================*/
double GbVoxelMatrix3D::getIntersectionRaytraceFactor(const double &x1, const double &x2, const double &x3,
                                                      const double &rx1, const double &rx2, const double &rx3)
{
    if (!((UbMath::equal(rx1, 0.0) || UbMath::equal(fabs(rx1), 1.0) ||
           UbMath::equal(fabs(rx1), vf::basics::constant::c1oSqrt2) || UbMath::equal(fabs(rx1), vf::basics::constant::c1oSqrt3)) &&
          (UbMath::equal(rx2, 0.0) || UbMath::equal(fabs(rx2), 1.0) ||
           UbMath::equal(fabs(rx2), vf::basics::constant::c1oSqrt2) || UbMath::equal(fabs(rx2), vf::basics::constant::c1oSqrt3)) &&
          (UbMath::equal(rx3, 0.0) || UbMath::equal(fabs(rx3), 1.0) ||
           UbMath::equal(fabs(rx3), vf::basics::constant::c1oSqrt2) || UbMath::equal(fabs(rx3), vf::basics::constant::c1oSqrt3)))) {
        throw UbException(UB_EXARGS, "nur fuer diskrete Boltzmannrichungen implementiert!!!");
    }

    // nachbarindex ermitteln
    int ndx1 = 0, ndx2 = 0, ndx3 = 0;
    if (UbMath::greater(rx1, 0.0))
        ndx1 = 1;
    else if (UbMath::less(rx1, 0.0))
        ndx1 = -1;
    if (UbMath::greater(rx2, 0.0))
        ndx2 = 1;
    else if (UbMath::less(rx2, 0.0))
        ndx2 = -1;
    if (UbMath::greater(rx3, 0.0))
        ndx3 = 1;
    else if (UbMath::less(rx3, 0.0))
        ndx3 = -1;

    int nix1 = UbMath::integerRounding((x1 - minX1) / deltaX1) + ndx1;
    int nix2 = UbMath::integerRounding((x2 - minX2) / deltaX2) + ndx2;
    int nix3 = UbMath::integerRounding((x3 - minX3) / deltaX3) + ndx3;

    // test ob nachbar solid
    if (nix1 >= 0 && nix2 >= 0 && nix3 >= 0 && nix1 < (int)voxelMatrix.getNX1() && nix2 < (int)voxelMatrix.getNX2() &&
        nix3 < (int)voxelMatrix.getNX3()) {
        if (UbMath::equal(voxelMatrix(nix1, nix2, nix3), SOLID)) {
            // return halber abstand der beiden knoten
            return 0.5 * sqrt((ndx1 * deltaX1) * (ndx1 * deltaX1) + (ndx2 * deltaX2) * (ndx2 * deltaX2) +
                              (ndx3 * deltaX3) * (ndx3 * deltaX3));
        }
    }

    return 0.0;
}
/*=======================================================*/
bool GbVoxelMatrix3D::isPointInGbObject3D(const double &x1p, const double &x2p, const double &x3p)
{
    int ix1 = UbMath::integerRounding((x1p - minX1) / deltaX1);
    int ix2 = UbMath::integerRounding((x2p - minX2) / deltaX2);
    int ix3 = UbMath::integerRounding((x3p - minX3) / deltaX3);

    if (ix1 >= 0 && ix2 >= 0 && ix3 >= 0 && ix1 < (int)voxelMatrix.getNX1() && ix2 < (int)voxelMatrix.getNX2() &&
        ix3 < (int)voxelMatrix.getNX3()) {
        if (UbMath::equal(voxelMatrix(ix1, ix2, ix3), SOLID))
            return true;
    }
    return false;
}
/*=======================================================*/
bool GbVoxelMatrix3D::isPointInGbObject3D(const double &x1p, const double &x2p, const double &x3p,
                                          bool &pointIsOnBoundary)
{
    pointIsOnBoundary = false;

    return isPointInGbObject3D(x1p, x2p, x3p);
}
/*=======================================================*/
bool GbVoxelMatrix3D::isCellInsideGbObject3D(const double & /*x1p1*/, const double & /*x2p1*/, const double & /*x3p1*/,
                                             const double & /*x1p2*/, const double & /*x2p2*/, const double & /*x3p2*/)
{
    return false;
    // FIXME: unreachable cide
    ////dass haengt von der Konfigration ab, aber meist ist der Block groesser wie etliche Poren ...

    //   //indizes ermitteln
    // int startix1 = (int)std::floor((x1p1-minX1)/deltaX1+1E-13);
    // int startix2 = (int)std::floor((x2p1-minX2)/deltaX2+1E-13);
    // int startix3 = (int)std::floor((x3p1-minX3)/deltaX3+1E-13);

    // if (startix1<0) return false;
    // if (startix2<0) return false;
    // if (startix3<0) return false;

    // int maxiX1 = (int)voxelMatrix.getNX1()-1;
    // int maxiX2 = (int)voxelMatrix.getNX2()-1;
    // int maxiX3 = (int)voxelMatrix.getNX3()-1;

    // int endix1 = (int)std::ceil((x1p2-minX1)/deltaX1-1E-13);
    // int endix2 = (int)std::ceil((x2p2-minX2)/deltaX2-1E-13);
    // int endix3 = (int)std::ceil((x3p2-minX3)/deltaX3-1E-13);

    // if (endix1>maxiX1) return false;
    // if (endix2>maxiX2) return false;
    // if (endix3>maxiX3) return false;

    // for (int ix3 = startix3; ix3<=endix3; ix3++)
    //   for (int ix2 = startix2; ix2<=endix2; ix2++)
    //      for (int ix1 = startix1; ix1<=endix1; ix1++)
    //         if (UbMath::equal(voxelMatrix(ix1, ix2, ix3), FLUID))
    //            return false;
    // return true;
}
/*=======================================================*/
bool GbVoxelMatrix3D::isCellCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a,
                                              const double &x1b, const double &x2b, const double &x3b)
// Merksatz: cell oder deren Volumen schneidet oder beinhaltet komplette oder Teile der CuboidUmrandung
// returns true:
//  - cell cuts  GbVoxelMatrix3D
//  - cell boxes GbVoxelMatrix3D
// returns false:
//  - cell completely inside GbVoxelMatrix3D
//  - cell und cuboid3D haben kein gemeinsames Volumen
{
    // erstmal die dumm Loesung
    if (!(this->isCellInsideGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b)) &&
        this->isCellInsideOrCuttingGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b)) {
        return true;
    }

    return false;
}
/*=======================================================*/
bool GbVoxelMatrix3D::isCellInsideOrCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a,
                                                      const double &x1b, const double &x2b, const double &x3b)
// returns true:
//  - cell completely inside cuboid3D ( = cuboid3D boxes cell)
//  - cell cuts  cuboid3D
//  - cell boxes cuboid3D
// returns false:
//  - cell und cuboid3D haben kein gemeinsames Volumen
{
    // simpler check, da unser GbCuboid3D ein AABB is:
    //  anfA        midA         endA             anfB    midB    endB
    //   |            x<-- dxA -->|                 |<-dxB->x       |
    //                |<----------------- T --------------->|
    // ist |T| <= dxA + dxB -> overlap!

    if (UbMath::lessEqual(std::fabs(this->getX1Centroid() - 0.5 * (x1b + x1a) /*Tx1*/),
                          0.5 * (this->getLengthX1() + std::fabs(x1b - x1a) /*dx1A+dx1B*/))

        && UbMath::lessEqual(std::fabs(this->getX2Centroid() - 0.5 * (x2b + x2a) /*Tx2*/),
                             0.5 * (this->getLengthX2() + std::fabs(x2b - x2a) /*dx2A+dx2B*/))

        && UbMath::lessEqual(std::fabs(this->getX3Centroid() - 0.5 * (x3b + x3a) /*Tx3*/),
                             0.5 * (this->getLengthX3() + std::fabs(x3b - x3a) /*dx3A+dx3B*/))) {
        return true;
    }

    return false;
}
/*=======================================================*/
vector<GbTriangle3D *> GbVoxelMatrix3D::getSurfaceTriangleSet()
{
    vector<GbTriangle3D *> triangles;
    cerr
        << "vector<GbTriangle3D*> GbVoxelMatrix3D::getSurfaceTriangleSet() - benoetigt MARCHING_CUBE paket aus 3rdParty"
        << endl;
    return triangles;
}
/*=======================================================*/
void GbVoxelMatrix3D::addSurfaceTriangleSet(vector<UbTupleFloat3> & /*nodes*/, vector<UbTupleInt3> & /*triangles*/)
{
    UBLOG(logINFO, " GbVoxelMatrix3D addSurfaceTriangleSet start")
    if (!this->addSurfaceTriangleSetFlag) {
        UBLOG(logINFO, " GbVoxelMatrix3D addSurfaceTriangleSet end without TriangleSetCreation")
        return;
    }

    cerr << "void GbVoxelMatrix3D.addSurfaceTriangleSet  - benoetigt MARCHING_CUBE paket aus 3rdParty" << endl;

    UBLOG(logINFO, " GbVoxelMatrix3D addSurfaceTriangleSet end")
}
/*=======================================================*/
string GbVoxelMatrix3D::toString() { return "GbVoxelMatrix3D"; }
/*=======================================================*/
void GbVoxelMatrix3D::readMatrixFromVtiASCIIFile(std::string filename)

{
    // UBLOG(logINFO,"  - create GbVoxelMatrix3D");
    UbFileInputASCII in(filename);
    // ifstream in(filename.c_str(), ios::binary);
    if (!in)
        throw UbException(UB_EXARGS, "could not open file " + filename);
    in.readLine();
    in.readLine();
    in.readLine();
    in.readLine();
    in.readLine();

    voxelMatrix = Matrix3D(nodesX1, nodesX2, nodesX3, GbVoxelMatrix3D::FLUID);

    // UBLOG(logINFO,"  - init values");
    int val;
    for (int x3 = 0; x3 < nodesX3; x3++)
        for (int x2 = 0; x2 < nodesX2; x2++)
            for (int x1 = 0; x1 < nodesX1; x1++) {
                val = in.readInteger();
                // if( !UbMath::equal(val, 0.0f) )
                // if( UbMath::greater(val, threshold) )
                if ((double)val >= lowerThreshold && (double)val <= upperThreshold) {
                    (voxelMatrix)(x1, x2, x3) = GbVoxelMatrix3D::SOLID;
                }
            }
    // UBLOG(logINFO,"  - create GbVoxelMatrix3D done");
}

//////////////////////////////////////////////////////////////////////////
void GbVoxelMatrix3D::rotate90aroundX(double /*cX1*/, double cX2, double cX3)
{
    //   double tempMinPunktX1 = minX1-cX1;
    double tempMinPunktX2 = minX2 - cX2;
    double tempMinPunktX3 = getX3Maximum() - cX3;

    //   double tempMinPunktX1tf = tempMinPunktX1;
    double tempMinPunktX2tf = -tempMinPunktX3;
    double tempMinPunktX3tf = tempMinPunktX2;

    //   double minX1_temp = tempMinPunktX1tf+cX1;
    double minX2_temp = tempMinPunktX2tf + cX2;
    double minX3_temp = tempMinPunktX3tf + cX3;

    minX2 = minX2_temp;
    minX3 = minX3_temp;

    int nx1 = (int)voxelMatrix.getNX1();
    int nx2 = (int)voxelMatrix.getNX2();
    int nx3 = (int)voxelMatrix.getNX3();

    int nx1_new = nx1;
    int nx2_new = nx3;
    int nx3_new = nx2;

    double delta_temp = deltaX2;
    deltaX2           = deltaX3;
    deltaX3           = delta_temp;

    GbVoxelMatrix3D::Matrix3D voxelMatrix_temp(nx1_new, nx2_new, nx3_new);

    for (int x3 = 0; x3 < nx3; x3++) {
        for (int x2 = 0; x2 < nx2; x2++) {
            for (int x1 = 0; x1 < nx1; x1++) {
                voxelMatrix_temp(x1, nx3 - x3 - 1, x2) = this->voxelMatrix(x1, x2, x3);
            }
        }
    }
    std::swap(this->voxelMatrix, voxelMatrix_temp);
}
//////////////////////////////////////////////////////////////////////////
void GbVoxelMatrix3D::rotate90aroundX()
{
    double cX1 = this->getX1Centroid();
    double cX2 = this->getX2Centroid();
    double cX3 = this->getX3Centroid();

    rotate90aroundX(cX1, cX2, cX3);
}
//////////////////////////////////////////////////////////////////////////
void GbVoxelMatrix3D::rotate90aroundY(double cX1, double /*cX2*/, double cX3)
{
    double tempMinPunktX1 = getX1Maximum() - cX1;
    //   double tempMinPunktX2 = minX2-cX2;
    double tempMinPunktX3 = minX3 - cX3;

    double tempMinPunktX1tf = tempMinPunktX3;
    //   double tempMinPunktX2tf = tempMinPunktX2;
    double tempMinPunktX3tf = -tempMinPunktX1;

    double minX1_temp = tempMinPunktX1tf + cX1;
    //   double minX2_temp = tempMinPunktX2tf+cX2;
    double minX3_temp = tempMinPunktX3tf + cX3;

    minX1 = minX1_temp;
    minX3 = minX3_temp;

    int nx1 = (int)voxelMatrix.getNX1();
    int nx2 = (int)voxelMatrix.getNX2();
    int nx3 = (int)voxelMatrix.getNX3();

    int nx1_new = nx3;
    int nx2_new = nx2;
    int nx3_new = nx1;

    double delta_temp = deltaX1;
    deltaX1           = deltaX3;
    deltaX3           = delta_temp;

    GbVoxelMatrix3D::Matrix3D voxelMatrix_temp(nx1_new, nx2_new, nx3_new);

    for (int x3 = 0; x3 < nx3; x3++) {
        for (int x2 = 0; x2 < nx2; x2++) {
            for (int x1 = 0; x1 < nx1; x1++) {
                voxelMatrix_temp(x3, x2, nx1 - x1 - 1) = this->voxelMatrix(x1, x2, x3);
            }
        }
    }
    std::swap(this->voxelMatrix, voxelMatrix_temp);
}
//////////////////////////////////////////////////////////////////////////
void GbVoxelMatrix3D::rotate90aroundY()
{
    double cX1 = this->getX1Centroid();
    double cX2 = this->getX2Centroid();
    double cX3 = this->getX3Centroid();

    rotate90aroundY(cX1, cX2, cX3);
}
//////////////////////////////////////////////////////////////////////////
void GbVoxelMatrix3D::rotate90aroundZ(double cX1, double cX2, double /*cX3*/)
{
    double tempMinPunktX1 = minX1 - cX1;
    double tempMinPunktX2 = getX2Maximum() - cX2;
    //   double tempMinPunktX3 = minX3-cX3;

    double tempMinPunktX1tf = -tempMinPunktX2;
    double tempMinPunktX2tf = tempMinPunktX1;
    //   double tempMinPunktX3tf = tempMinPunktX3;

    double minX1_temp = tempMinPunktX1tf + cX1;
    double minX2_temp = tempMinPunktX2tf + cX2;
    //   double minX3_temp = tempMinPunktX3tf+cX3;

    minX1 = minX1_temp;
    minX2 = minX2_temp;

    int nx1 = (int)voxelMatrix.getNX1();
    int nx2 = (int)voxelMatrix.getNX2();
    int nx3 = (int)voxelMatrix.getNX3();

    int nx1_new = nx2;
    int nx2_new = nx1;
    int nx3_new = nx3;

    double delta_temp = deltaX1;
    deltaX1           = deltaX2;
    deltaX2           = delta_temp;

    GbVoxelMatrix3D::Matrix3D voxelMatrix_temp(nx1_new, nx2_new, nx3_new);

    for (int x3 = 0; x3 < nx3; x3++) {
        for (int x2 = 0; x2 < nx2; x2++) {
            for (int x1 = 0; x1 < nx1; x1++) {
                voxelMatrix_temp(nx2 - x2 - 1, x1, x3) = this->voxelMatrix(x1, x2, x3);
            }
        }
    }
    std::swap(this->voxelMatrix, voxelMatrix_temp);
}
//////////////////////////////////////////////////////////////////////////
void GbVoxelMatrix3D::rotate90aroundZ()
{
    double cX1 = this->getX1Centroid();
    double cX2 = this->getX2Centroid();
    double cX3 = this->getX3Centroid();

    rotate90aroundZ(cX1, cX2, cX3);
}
//////////////////////////////////////////////////////////////////////////
void GbVoxelMatrix3D::mirrorX()
{
    int nx1 = (int)voxelMatrix.getNX1();
    int nx2 = (int)voxelMatrix.getNX2();
    int nx3 = (int)voxelMatrix.getNX3();

    GbVoxelMatrix3D::Matrix3D voxelMatrix_temp(nx1, nx2, nx3);

    for (int x3 = 0; x3 < nx3; x3++) {
        for (int x2 = 0; x2 < nx2; x2++) {
            for (int x1 = 0; x1 < nx1; x1++) {
                voxelMatrix_temp(nx1 - x1 - 1, x2, x3) = this->voxelMatrix(x1, x2, x3);
            }
        }
    }
    std::swap(this->voxelMatrix, voxelMatrix_temp);
}
//////////////////////////////////////////////////////////////////////////
void GbVoxelMatrix3D::mirrorY()
{
    int nx1 = (int)voxelMatrix.getNX1();
    int nx2 = (int)voxelMatrix.getNX2();
    int nx3 = (int)voxelMatrix.getNX3();

    GbVoxelMatrix3D::Matrix3D voxelMatrix_temp(nx1, nx2, nx3);

    for (int x3 = 0; x3 < nx3; x3++) {
        for (int x2 = 0; x2 < nx2; x2++) {
            for (int x1 = 0; x1 < nx1; x1++) {
                voxelMatrix_temp(x1, nx2 - x2 - 1, x3) = this->voxelMatrix(x1, x2, x3);
            }
        }
    }
    std::swap(this->voxelMatrix, voxelMatrix_temp);
}
//////////////////////////////////////////////////////////////////////////
void GbVoxelMatrix3D::mirrorZ()
{
    int nx1 = (int)voxelMatrix.getNX1();
    int nx2 = (int)voxelMatrix.getNX2();
    int nx3 = (int)voxelMatrix.getNX3();

    GbVoxelMatrix3D::Matrix3D voxelMatrix_temp(nx1, nx2, nx3);

    for (int x3 = 0; x3 < nx3; x3++) {
        for (int x2 = 0; x2 < nx2; x2++) {
            for (int x1 = 0; x1 < nx1; x1++) {
                voxelMatrix_temp(x1, x2, nx3 - x3 - 1) = this->voxelMatrix(x1, x2, x3);
            }
        }
    }
    std::swap(this->voxelMatrix, voxelMatrix_temp);
}
//////////////////////////////////////////////////////////////////////////
void GbVoxelMatrix3D::rotateAroundY(double theta)
{
    int nx1 = (int)voxelMatrix.getNX1();
    int nx2 = (int)voxelMatrix.getNX2();
    int nx3 = (int)voxelMatrix.getNX3();

    GbVoxelMatrix3D::Matrix3D voxelMatrix_temp(nx1, nx2, nx3, FLUID);

    for (int x3 = 0; x3 < nx3; x3++) {
        for (int x2 = 0; x2 < nx2; x2++) {
            for (int x1 = 0; x1 < nx1; x1++) {
                if (voxelMatrix(x1, x2, x3) == SOLID) {
                    double rcX1 = minX1 + deltaX1 * x1;
                    double rcX3 = minX3 + deltaX3 * x3;

                    double nrcX1 = cos(theta) * rcX1 + sin(theta) * rcX3;
                    double nrcX3 = -sin(theta) * rcX1 + cos(theta) * rcX3;

                    int newX1 = UbMath::integerRounding((nrcX1 - minX1) / deltaX1);
                    int newX2 = x2;
                    int newX3 = UbMath::integerRounding((nrcX3 - minX3) / deltaX3);

                    if (newX1 > 0 && newX3 > 0 && newX1 < nx1 && newX3 < nx3) {
                        voxelMatrix_temp(newX1, newX2, newX3) = voxelMatrix(x1, x2, x3);
                    }

                    // int ix1, ix2, ix3;
                    // double ixx1 = (abs(nrcX1-minX1)/deltaX1);
                    // ix2 = x2;
                    // double ixx3 = (abs(nrcX3-minX3)/deltaX3);
                    // if (ixx1-(int)ixx1>.9999999999) ix1 = (int)ixx1+1; else ix1 = (int)ixx1;
                    ////if (ixx2-(int)ixx2>.9999999999) ix2 = (int)ixx2+1; else ix2 = (int)ixx2;
                    // if (ixx3-(int)ixx3>.9999999999) ix3 = (int)ixx3+1; else ix3 = (int)ixx3;

                    // if (ix1>=0 && ix3>=0 && ix1<nx1 && ix3<nx3)
                    //{
                    //   voxelMatrix_temp(ix1, ix2, ix3) = voxelMatrix(x1, x2, x3);
                    //}
                }
            }
        }
    }
    std::swap(voxelMatrix, voxelMatrix_temp);

    for (int x3 = 0; x3 < nx3; x3++)
        for (int x2 = 0; x2 < nx2; x2++)
            for (int x1 = 0; x1 < nx1; x1++) {
                int count = 0;
                for (int k3 = -1; k3 <= 1; k3++)
                    for (int k1 = -1; k1 <= 1; k1++) {
                        int j1 = x1 + k1;
                        int j3 = x3 + k3;
                        if (j1 >= 0 && j3 >= 0 && j1 < nx1 && j3 < nx3) {
                            if (voxelMatrix(j1, x2, j3) == SOLID) {
                                count++;
                            }
                        }
                    }
                if (count == 8) {
                    voxelMatrix(x1, x2, x3) = SOLID;
                }
            }
}
//////////////////////////////////////////////////////////////////////////
void GbVoxelMatrix3D::writeToLegacyVTKASCII(const std::string &fileName)
{
    string fn = fileName + ".ascii.vtk";

    FILE *file;
    file = fopen(fn.c_str(), "w");

    if (file == NULL) {
        std::string pathf = UbSystem::getPathFromString(fn);
        if (fn.size() > 0) {
            UbSystem::makeDirectory(pathf);
            file = fopen(fn.c_str(), "w");
        }
        if (file == NULL)
            throw UbException(UB_EXARGS, "can not open " + fn);
    }

    if (file == NULL)
        throw UbException(UB_EXARGS, "can not open " + fn);

    int nx1 = (int)voxelMatrix.getNX1();
    int nx2 = (int)voxelMatrix.getNX2();
    int nx3 = (int)voxelMatrix.getNX3();

    int nn = nx1 * nx2 * nx3;

    fprintf(file, "# vtk DataFile Version 2.0\n");
    fprintf(file, "vtk output\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET STRUCTURED_POINTS\n");
    fprintf(file, "DIMENSIONS %d %d %d\n", nx1, nx2, nx3);
    fprintf(file, "ORIGIN %g %g %g\n", minX1, minX2, minX3);
    fprintf(file, "SPACING %g %g %g\n", deltaX1, deltaX2, deltaX3);
    fprintf(file, "POINT_DATA %d\n", nn);
    fprintf(file, "SCALARS Geo float\n");
    fprintf(file, "LOOKUP_TABLE default\n");

    for (int k = 0; k < nx3; k++) {
        for (int j = 0; j < nx2; j++) {
            for (int i = 0; i < nx1; i++) {
                fprintf(file, "%g ", voxelMatrix(i, j, k));
            }
        }
    }

    fprintf(file, "\n");

    fclose(file);
}
//////////////////////////////////////////////////////////////////////////
void GbVoxelMatrix3D::writeToLegacyVTKBinary(const std::string &fileName)
{
    string fn = fileName + ".binary.vtk";

    FILE *file;
    file = fopen(fn.c_str(), "w");

    if (file == NULL) {
        std::string pathf = UbSystem::getPathFromString(fn);
        if (fn.size() > 0) {
            UbSystem::makeDirectory(pathf);
            file = fopen(fn.c_str(), "w");
        }
        if (file == NULL)
            throw UbException(UB_EXARGS, "can not open " + fn);
    }

    int nx1 = (int)voxelMatrix.getNX1();
    int nx2 = (int)voxelMatrix.getNX2();
    int nx3 = (int)voxelMatrix.getNX3();

    int nn = nx1 * nx2 * nx3;

    char LF = 0x0A;

    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "vtk output\n");
    fprintf(file, "BINARY\n");
    fprintf(file, "DATASET STRUCTURED_POINTS\n");
    fprintf(file, "DIMENSIONS %d %d %d\n", nx1, nx2, nx3);
    fprintf(file, "ORIGIN %g %g %g\n", minX1, minX2, minX3);
    fprintf(file, "SPACING %g %g %g\n", deltaX1, deltaX2, deltaX3);
    fprintf(file, "POINT_DATA %d\n", nn);
    fprintf(file, "SCALARS Geo float\n");
    fprintf(file, "LOOKUP_TABLE default");
    fclose(file);

    GbVoxelMatrix3D::Matrix3D voxelMatrix_temp(nx1, nx2, nx3);

    if (UbSystem::isLittleEndian()) {
        for (int x3 = 0; x3 < nx3; x3++) {
            for (int x2 = 0; x2 < nx2; x2++) {
                for (int x1 = 0; x1 < nx1; x1++) {
                    float tmp = this->voxelMatrix(x1, x2, x3);
                    UbSystem::swapByteOrder((unsigned char *)(&(tmp)), sizeof(float));
                    voxelMatrix_temp(x1, x2, x3) = tmp;
                }
            }
        }
    }

    file = fopen(fn.c_str(), "ab");

    fwrite(&LF, sizeof(char), 1, file);
    fwrite(voxelMatrix_temp.getStartAdressOfSortedArray(0, 0, 0), sizeof(float),
           voxelMatrix_temp.getDataVector().size(), file);
    fwrite(&LF, sizeof(char), 1, file);
    fclose(file);
}
//////////////////////////////////////////////////////////////////////////
void GbVoxelMatrix3D::writeToVTKImageDataASCII(const std::string &fileName)
{
    int nx1 = (int)voxelMatrix.getNX1();
    int nx2 = (int)voxelMatrix.getNX2();
    int nx3 = (int)voxelMatrix.getNX3();

    string fn = fileName + ".ascii.vti";

    FILE *file;
    file = fopen(fn.c_str(), "w");

    if (file == NULL) {
        std::string pathf = UbSystem::getPathFromString(fn);
        if (fn.size() > 0) {
            UbSystem::makeDirectory(pathf);
            file = fopen(fn.c_str(), "w");
        }
        if (file == NULL)
            throw UbException(UB_EXARGS, "can not open " + fn);
    }

    fprintf(file,
            "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"); // paraview
                                                                                                                  // 4.1
    // fprintf(file,"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"); //paraview 3.1
    fprintf(file, "  <ImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"%g %g %g\" Spacing=\"%g %g %g\">\n", 0,
            nx1 - 1, 0, nx2 - 1, 0, nx3 - 1, minX1, minX2, minX3, deltaX1, deltaX2, deltaX3);
    fprintf(file, "  <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, nx1 - 1, 0, nx2 - 1, 0, nx3 - 1);
    fprintf(file, "    <PointData Scalars=\"VoxelMatrix\">\n");
    fprintf(file, "      <DataArray type=\"Float32\" Name=\"VoxelMatrix\" format=\"ascii\" RangeMin=\"0\" "
                  "RangeMax=\"1\">\n        ");

    for (int k = 0; k < nx3; k++) {
        for (int j = 0; j < nx2; j++) {
            for (int i = 0; i < nx1; i++) {
                fprintf(file, "%g ", voxelMatrix(i, j, k));
            }
        }
    }

    fprintf(file, "\n      </DataArray>\n");
    fprintf(file, "    </PointData>\n");
    fprintf(file, "    <CellData>\n");
    fprintf(file, "    </CellData>\n");
    fprintf(file, "  </Piece>\n");
    fprintf(file, "  </ImageData>\n");
    fprintf(file, "</VTKFile>\n");

    fclose(file);
}
//////////////////////////////////////////////////////////////////////////
void GbVoxelMatrix3D::writeToVTKImageDataAppended(const std::string &fileName)
{
    int nx1 = (int)voxelMatrix.getNX1();
    int nx2 = (int)voxelMatrix.getNX2();
    int nx3 = (int)voxelMatrix.getNX3();

    string fn = fileName + ".appended.vti";

    FILE *file;
    file = fopen(fn.c_str(), "w");

    if (file == NULL) {
        std::string pathf = UbSystem::getPathFromString(fn);
        if (fn.size() > 0) {
            UbSystem::makeDirectory(pathf);
            file = fopen(fn.c_str(), "w");
        }
        if (file == NULL)
            throw UbException(UB_EXARGS, "can not open " + fn);
    }

    fprintf(file,
            "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"); // paraview
                                                                                                                  // 4.1
    fprintf(file, "  <ImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"%g %g %g\" Spacing=\"%g %g %g\">\n", 0,
            nx1 - 1, 0, nx2 - 1, 0, nx3 - 1, minX1, minX2, minX3, deltaX1, deltaX2, deltaX3);
    fprintf(file, "  <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, nx1 - 1, 0, nx2 - 1, 0, nx3 - 1);
    fprintf(file, "    <PointData Scalars=\"VoxelMatrix\">\n");
    fprintf(file, "      <DataArray type=\"Float32\" Name=\"VoxelMatrix\" format=\"appended\" RangeMin=\"0\" "
                  "RangeMax=\"1\" offset=\"0\" />\n");
    fprintf(file, "    </PointData>\n");
    fprintf(file, "    <CellData>\n");
    fprintf(file, "    </CellData>\n");
    fprintf(file, "  </Piece>\n");
    fprintf(file, "  </ImageData>\n");
    fprintf(file, "  <AppendedData encoding=\"raw\">\n");
    fprintf(file, "   _");
    fclose(file);

    file                    = fopen(fn.c_str(), "ab");
    unsigned long long size = (unsigned long long)voxelMatrix.getDataVector().size() * sizeof(float);
    fwrite(&size, sizeof(unsigned long long), 1, file);
    fwrite(voxelMatrix.getStartAdressOfSortedArray(0, 0, 0), sizeof(float), voxelMatrix.getDataVector().size(), file);
    fclose(file);

    file = fopen(fn.c_str(), "a");
    fprintf(file, "\n");
    fprintf(file, "  </AppendedData>\n");
    fprintf(file, "</VTKFile>\n");
    fclose(file);
}

//! \}
