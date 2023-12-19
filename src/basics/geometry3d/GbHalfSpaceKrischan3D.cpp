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
#include <geometry3d/GbHalfSpaceKrischan3D.h>

using namespace std;

/*==========================================================*/
GbHalfSpaceKrischan3D::GbHalfSpaceKrischan3D(GbTriangle3D *triangle)
{
    GbPoint3D *PointA = triangle->getPoint1();
    GbPoint3D *PointB = triangle->getPoint2();
    GbPoint3D *PointC = triangle->getPoint3();

    GbVector3D A(PointA->x1, PointA->x2, PointA->x3);
    GbVector3D BA(PointB->x1 - PointA->x1, PointB->x2 - PointA->x2, PointB->x3 - PointA->x3);
    GbVector3D CA(PointC->x1 - PointA->x1, PointC->x2 - PointA->x2, PointC->x3 - PointA->x3);
    GbVector3D BACA = BA.Cross(CA);
    // this->Normal = PointB->subtract(PointA)->cross(PointC->subtract(PointA))->normalize();
    BACA.Normalize();
    this->Normal = BACA;
    this->d      = this->Normal.Dot(A);
}
/*==========================================================*/
GbHalfSpaceKrischan3D::GbHalfSpaceKrischan3D(GbPoint3D *PointA, GbPoint3D *PointB, GbPoint3D *PointC)
{
    GbVector3D A(PointA->x1, PointA->x2, PointA->x3);
    GbVector3D BA(PointB->x1 - PointA->x1, PointB->x2 - PointA->x2, PointB->x3 - PointA->x3);
    GbVector3D CA(PointC->x1 - PointA->x1, PointC->x2 - PointA->x2, PointC->x3 - PointA->x3);
    GbVector3D BACA = BA.Cross(CA);
    // this->Normal = PointB->subtract(PointA)->cross(PointC->subtract(PointA))->normalize();
    BACA.Normalize();
    this->Normal = BACA;
    this->d      = this->Normal.Dot(A);
}
/*==========================================================*/
GbHalfSpaceKrischan3D::GbHalfSpaceKrischan3D(GbPoint3D *PointA, GbPoint3D *PointB)
{
    GbVector3D A(PointA->x1, PointA->x2, PointA->x3);
    GbVector3D B(PointB->x1, PointB->x2, PointB->x3);
    GbVector3D K(0.0, 0.0, 0.99); // the vector from PointA - third point

    GbVector3D PointBA  = B - A;
    GbVector3D PointBAK = PointBA.Cross(K);
    PointBAK.Normalize();
    this->Normal = PointBAK;
    this->d      = this->Normal.Dot(A);
}
/*==========================================================*/
GbHalfSpaceKrischan3D::GbHalfSpaceKrischan3D(const double &p1x, const double &p1y, const double &p1z, const double &p2x,
                                             const double &p2y, const double &p2z, const double &p3x, const double &p3y,
                                             const double &p3z)
{
    GbVector3D A(p1x, p1y, p1z);
    GbVector3D BA(p2x - p1x, p2y - p1y, p2z - p1z);
    GbVector3D CA(p3x - p1x, p3y - p1y, p3z - p1z);
    GbVector3D BACA = BA.Cross(CA);

    BACA.Normalize();
    this->Normal = BACA;
    this->d      = this->Normal.Dot(A);
}
/*==========================================================*/
GbHalfSpaceKrischan3D::GbHalfSpaceKrischan3D(double nx, double ny, double nz, double dist)
{
    this->Normal = GbVector3D(nx, ny, nz);
    this->Normal.Normalize();
    this->d = dist;
}
/*==========================================================*/
double GbHalfSpaceKrischan3D::getCellVolumeInsideGbObject3D(const double &x1a, const double &x2a, const double &x3a,
                                                            const double &x1b, const double &x2b, const double &x3b)
{

    double x1 = x1b - x1a;
    double x2 = x2b - x2a;
    double x3 = x3b - x3a;

    if (this->isCellInsideGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b))
        return 1.0 * x1 * x2 * x3;
    if (!(this->isCellCuttingGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b)))
        return 0.0;

    double alpha = 0.0;
    double internX1, internX2, internX3;

    for (int x1vers = 0; x1vers < 2; x1vers++) {
        for (int x2vers = 0; x2vers < 2; x2vers++) {
            for (int x3vers = 0; x3vers < 2; x3vers++) {
                internX1 = x1a + (x1b - x1a) * x1vers;
                internX2 = x2a + (x2b - x2a) * x2vers;
                internX3 = x3a + (x3b - x3a) * x3vers;

                // if point is INSIDE the halfspace, distance is smaller than zero
                // --> loop determines the minimum alpha...i.e. the alpha with maximum absolute value for all points
                // INSIDE the halfspace
                if (UbMath::lessEqual(this->getDistance(internX1, internX2, internX3), alpha))
                    alpha = this->getDistance(internX1, internX2, internX3);
                // cout<<zelltyp<<" "<<kugel->getDistance(internX1,internX2,internX3)<<" "<<alpha<<endl;
            } // end first for
        }     // end second for
    }         // end third for

    // PLIC needs alphas > 0.0
    alpha = (-1) * alpha;

    double n[3];
    n[0] = this->Normal[0];
    n[1] = this->Normal[1];
    n[2] = this->Normal[2];

    // cout << "Koordinaten:  "<<x1<<" "<<x2<<" "<<x3<<endl;
    // cout << "Deltas:       "<<deltaX1<<" "<<deltaX2<<" "<<deltaX3<<endl;
    // cout << "Halbe Zelle:  "<<halfcelldelta<<endl;

    // cout<<"Centroid:  "<<kugel->getX1Centroid()<<" "<<kugel->getX2Centroid()<<" "<<kugel->getX3Centroid()<<endl;

    // cout<<"Normals: "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;

    double normLength;
    normLength = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    n[0] /= normLength;
    n[1] /= normLength;
    n[2] /= normLength;

    if (UbMath::less(n[0], 0.0))
        n[0] = -n[0];
    if (UbMath::less(n[1], 0.0))
        n[1] = -n[1];
    if (UbMath::less(n[2], 0.0))
        n[2] = -n[2];

    // cout<<"Normals: "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;

    double dummy;
    if (UbMath::greater(n[0], n[1])) {
        dummy = n[1];
        n[1]  = n[0];
        n[0]  = dummy;
    }
    if (UbMath::greater(n[1], n[2])) {
        dummy = n[2];
        n[2]  = n[1];
        n[1]  = dummy;
    }
    if (UbMath::greater(n[0], n[1])) {
        dummy = n[1];
        n[1]  = n[0];
        n[0]  = dummy;
    }

    // cout<<"Normals: "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;

    double n1, n2, n3;
    n1 = n[0];
    n2 = n[1];
    n3 = n[2];

    double preresult = 0.0, result = 0.0;

    // 1D Check
    if (UbMath::lessEqual(n1, 0.00001) && UbMath::lessEqual(n2, 0.00001)) {
        result = alpha * x1 * x2;
    }
    // 2D Check
    else if (UbMath::lessEqual(n1, 0.00001)) {
        preresult = (2 * n2 * n3);
        result    = (alpha * alpha) / preresult;

        if (UbMath::greater(alpha, n2 * x2)) {
            result += -(alpha - n2 * x2) * (alpha - n2 * x2) / preresult;
        }
        if (UbMath::greater(alpha, n3 * x3)) {
            result += -(alpha - n3 * x3) * (alpha - n3 * x3) / preresult;
        }
        if (UbMath::greater(alpha, n2 * x2 + n3 * x3)) {
            result += (alpha - n2 * x2 - n3 * x3) * (alpha - n2 * x2 - n3 * x3) / preresult;
        }

        // tiefenrichtung mit einmultiplizieren...
        result *= x1;
    }
    // 3D Check
    else {
        preresult = 6 * n1 * n2 * n3;

        result = alpha * alpha * alpha / preresult;

        if (UbMath::greater(alpha, n1 * x1)) {
            result += -((alpha - n1 * x1) * (alpha - n1 * x1) * (alpha - n1 * x1)) / preresult;
        }
        if (UbMath::greater(alpha, n2 * x2)) {
            result += -((alpha - n2 * x2) * (alpha - n2 * x2) * (alpha - n2 * x2)) / preresult;
        }
        if (UbMath::greater(alpha, n3 * x3)) {
            result += -((alpha - n3 * x3) * (alpha - n3 * x3) * (alpha - n3 * x3)) / preresult;
        }
        if (UbMath::greater(alpha, (n1 * x1 + n2 * x2))) {
            result += ((alpha - (n1 * x1 + n2 * x2)) * (alpha - (n1 * x1 + n2 * x2)) * (alpha - (n1 * x1 + n2 * x2))) /
                      preresult;
        }
        if (UbMath::greater(alpha, (n1 * x1 + n3 * x3))) {
            result += ((alpha - (n1 * x1 + n3 * x3)) * (alpha - (n1 * x1 + n3 * x3)) * (alpha - (n1 * x1 + n3 * x3))) /
                      preresult;
        }
        if (UbMath::greater(alpha, (n2 * x2 + n3 * x3))) {
            result += ((alpha - (n2 * x2 + n3 * x3)) * (alpha - (n2 * x2 + n3 * x3)) * (alpha - (n2 * x2 + n3 * x3))) /
                      preresult;
        }

        // NEW
        if (UbMath::greater(alpha, (n1 * x1 + n2 * x2 + n3 * x3))) {
            result += -((alpha - (n1 * x1 + n2 * x2 + n3 * x3)) * (alpha - (n1 * x1 + n2 * x2 + n3 * x3)) *
                        (alpha - (n1 * x1 + n2 * x2 + n3 * x3))) /
                      preresult;
        }
    }

    if (!UbMath::inClosedInterval(result / (x1 * x2 * x3), -0.01, 1.01)) {
        stringstream errMsg;

        errMsg << "Danger...Fuellstand " << result << " nicht im Interfall [0.0..1.0]" << endl;
        errMsg << "NormVec: " << n1 << " " << n2 << " " << n3 << endl;
        errMsg << "Cell:    " << x1 << " " << x2 << " " << x3 << endl;
        errMsg << "Alpha:   " << alpha << endl;

        throw UbException(UB_EXARGS, errMsg.str());
    }

    return result;

    // double eps=0.0;
    // if( UbMath::equal(n1,0.0) && UbMath::equal(n2,0.0) )
    //{
    //   eps = alpha/n3;
    //}
    // else if( UbMath::equal(n1,0.0) )
    //{
    //   double dim1,dim2;
    //   dim1 = alpha/n2;
    //   dim2 = alpha/n3;

    //   eps = 0.5*dim1*dim2;
    //   if( UbMath::greater(dim1,1.0) )   eps -= 0.5*(dim1-1.0)*dim2/dim1*(dim1-1.0);
    //   if( UbMath::greater(dim2,1.0) )   eps -= 0.5*(dim2-1.0)*dim1/dim2*(dim2-1.0);
    //}
    // else
    //{
    //   eps = alpha*alpha*alpha;
    //   if( UbMath::greater(alpha,n1) )
    //      eps -= (alpha-n1)*(alpha-n1)*(alpha-n1);
    //   if( UbMath::greater(alpha,n2) )
    //      eps -= (alpha-n2)*(alpha-n2)*(alpha-n2);
    //   if( UbMath::greater(alpha,n3) )
    //      eps -= (alpha-n3)*(alpha-n3)*(alpha-n3);

    //   if( UbMath::greater(alpha,n1+n2) )
    //      eps += (alpha-n1-n2)*(alpha-n1-n2)*(alpha-n1-n2);
    //   if( UbMath::greater(alpha,n1+n3) )
    //      eps += (alpha-n1-n3)*(alpha-n1-n3)*(alpha-n1-n3);
    //   if( UbMath::greater(alpha,n2+n3) )
    //      eps += (alpha-n2-n3)*(alpha-n2-n3)*(alpha-n2-n3);

    //   //attention: use without delta_i
    //   eps = eps / (6*n[0]*n[1]*n[2]);

    //   eps = eps / (deltaX1*deltaX2*deltaX3);
    //}

    // return(eps) ;
    // cout << "alpha ist " << alpha << endl;
    // cout << "fillLevel ist " << eps << endl;
}
/*==========================================================*/

//! \}
