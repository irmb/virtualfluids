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
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#include <GbCuboid3D.h>
#include <GbSystem3D.h>
#include <GbTriangle3D.h>

#include <basics/utilities/UbMath.h>

using namespace std;

/*=======================================================*/
// Konstruktor
GbCuboid3D::GbCuboid3D() : GbObject3D()
{
    this->setName("cuboid");
    this->p1 = new GbPoint3D(0.0, 0.0, 0.0);
    this->p2 = new GbPoint3D(0.0, 0.0, 0.0);
    this->p1->addObserver(this);
    this->p2->addObserver(this);
}
/*=======================================================*/
GbCuboid3D::GbCuboid3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b, const double &x2b,
                       const double &x3b)
    : GbObject3D()
{
    this->setName("cuboid");
    this->p1 = new GbPoint3D(x1a, x2a, x3a);
    this->p1->addObserver(this);
    this->p2 = new GbPoint3D(x1b, x2b, x3b);
    this->p2->addObserver(this);
}
/*=======================================================*/
GbCuboid3D::GbCuboid3D(GbPoint3D *p1, GbPoint3D *p2) : GbObject3D()
{
    this->setName("cuboid");
    if (!p1 || !p2)
        throw UbException(UB_EXARGS, "one point ==NULL");
    this->p1 = p1;
    this->p1->addObserver(this);
    this->p2 = p2;
    this->p2->addObserver(this);
}
/*=======================================================*/
GbCuboid3D::GbCuboid3D(GbCuboid3D *cuboid) : GbObject3D()
{
    this->setName("cuboid");
    if (!cuboid->getPoint1() || !cuboid->getPoint2())
        throw UbException(UB_EXARGS, "cuboid ==NULL");
    this->p1 = cuboid->getPoint1()->clone();
    this->p1->addObserver(this);
    this->p2 = cuboid->getPoint2()->clone();
    this->p2->addObserver(this);
}
/*=======================================================*/
// Destruktor
GbCuboid3D::~GbCuboid3D()
{
    // cout<<"~GbCuboid3D()"<<endl;
    if (this->p1)
        this->p1->removeObserver(this);
    if (this->p2)
        this->p2->removeObserver(this);
}
/*=======================================================*/
void GbCuboid3D::finalize()
{
    if (this->p1) {
        this->p1->removeObserver(this);
        this->p1->finalize();
        delete this->p1;
        this->p1 = NULL;
    }
    if (this->p2) {
        this->p2->removeObserver(this);
        this->p2->finalize();
        delete this->p2;
        this->p2 = NULL;
    }
}
/*=======================================================*/
void GbCuboid3D::setPoint1(GbPoint3D *point1)
{
    if (this->p1)
        this->p1->removeObserver(this);
    this->p1 = point1;
    this->p1->addObserver(this);

    this->notifyObserversObjectChanged();
}
/*=======================================================*/
void GbCuboid3D::setPoint2(GbPoint3D *point2)
{
    if (this->p2)
        this->p2->removeObserver(this);
    this->p2 = point2;
    this->p2->addObserver(this);

    this->notifyObserversObjectChanged();
}
/*=======================================================*/
void GbCuboid3D::setPoints(GbPoint3D *point1, GbPoint3D *point2)
{
    if (this->p1)
        this->p1->removeObserver(this);
    if (this->p2)
        this->p2->removeObserver(this);

    this->p1 = point1;
    this->p2 = point2;

    this->p1->addObserver(this);
    this->p2->addObserver(this);

    this->notifyObserversObjectChanged();
}
/*=======================================================*/
void GbCuboid3D::setCenterCoordinates(const double &x1, const double &x2, const double &x3)
{
    this->translate(x1 - getX1Centroid(), x2 - getX2Centroid(), x3 - getX3Centroid());
}
/*=======================================================*/
double GbCuboid3D::getX1Centroid() { return (0.5 * (p1->x1 + p2->x1)); }
/*=======================================================*/
double GbCuboid3D::getX1Minimum() { return (this->p1->x1 < this->p2->x1 ? this->p1->x1 : this->p2->x1); }
/*=======================================================*/
double GbCuboid3D::getX1Maximum() { return (this->p1->x1 > this->p2->x1 ? this->p1->x1 : this->p2->x1); }
/*=======================================================*/
double GbCuboid3D::getX2Centroid() { return (0.5 * (p1->x2 + p2->x2)); }
/*=======================================================*/
double GbCuboid3D::getX2Minimum() { return (this->p1->x2 < this->p2->x2 ? this->p1->x2 : this->p2->x2); }
/*=======================================================*/
double GbCuboid3D::getX2Maximum() { return (this->p1->x2 > this->p2->x2 ? this->p1->x2 : this->p2->x2); }
/*=======================================================*/
double GbCuboid3D::getX3Centroid() { return (0.5 * (p1->x3 + p2->x3)); }
/*=======================================================*/
double GbCuboid3D::getX3Minimum() { return (this->p1->x3 < this->p2->x3 ? this->p1->x3 : this->p2->x3); }
/*=======================================================*/
double GbCuboid3D::getX3Maximum() { return (this->p1->x3 > this->p2->x3 ? this->p1->x3 : this->p2->x3); }
/*=======================================================*/
double GbCuboid3D::getLengthX1() { return (this->getX1Maximum() - this->getX1Minimum()); }
/*=======================================================*/
double GbCuboid3D::getLengthX2() { return (this->getX2Maximum() - this->getX2Minimum()); }
/*=======================================================*/
double GbCuboid3D::getLengthX3() { return (this->getX3Maximum() - this->getX3Minimum()); }
/*=======================================================*/
bool GbCuboid3D::isPointInGbObject3D(const double &x1p, const double &x2p, const double &x3p)
{
    // true, wenn 'in Object' oder 'auf Boundary'!
    if (UbMath::less(x1p, this->getX1Minimum()))
        return false;
    else if (UbMath::less(x2p, this->getX2Minimum()))
        return false;
    else if (UbMath::less(x3p, this->getX3Minimum()))
        return false;
    else if (UbMath::greater(x1p, this->getX1Maximum()))
        return false;
    else if (UbMath::greater(x2p, this->getX2Maximum()))
        return false;
    else if (UbMath::greater(x3p, this->getX3Maximum()))
        return false;

    return true;
}
/*=======================================================*/
bool GbCuboid3D::isPointInGbObject3D(const double &x1p, const double &x2p, const double &x3p, bool &pointIsOnBoundary)
{
    pointIsOnBoundary = false;

    // true, wenn 'in Object' oder 'auf Boundary'!
    if (UbMath::less(x1p, this->getX1Minimum()))
        return false;
    else if (UbMath::less(x2p, this->getX2Minimum()))
        return false;
    else if (UbMath::less(x3p, this->getX3Minimum()))
        return false;
    else if (UbMath::greater(x1p, this->getX1Maximum()))
        return false;
    else if (UbMath::greater(x2p, this->getX2Maximum()))
        return false;
    else if (UbMath::greater(x3p, this->getX3Maximum()))
        return false;

    if (UbMath::equal(x1p, this->getX1Minimum()))
        pointIsOnBoundary = true;
    else if (UbMath::equal(x2p, this->getX2Minimum()))
        pointIsOnBoundary = true;
    else if (UbMath::equal(x3p, this->getX3Minimum()))
        pointIsOnBoundary = true;
    else if (UbMath::equal(x1p, this->getX1Maximum()))
        pointIsOnBoundary = true;
    else if (UbMath::equal(x2p, this->getX2Maximum()))
        pointIsOnBoundary = true;
    else if (UbMath::equal(x3p, this->getX3Maximum()))
        pointIsOnBoundary = true;

    return true;
}
/*=======================================================*/
bool GbCuboid3D::isCellInsideGbObject3D(const double &x1p1, const double &x2p1, const double &x3p1, const double &x1p2,
                                        const double &x2p2, const double &x3p2)
{
    if (UbMath::less(x1p1, this->getX1Minimum()))
        return false;
    else if (UbMath::less(x2p1, this->getX2Minimum()))
        return false;
    else if (UbMath::less(x3p1, this->getX3Minimum()))
        return false;
    else if (UbMath::greater(x1p2, this->getX1Maximum()))
        return false;
    else if (UbMath::greater(x2p2, this->getX2Maximum()))
        return false;
    else if (UbMath::greater(x3p2, this->getX3Maximum()))
        return false;

    return true;
}
/*=======================================================*/
bool GbCuboid3D::isCellCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                         const double &x2b, const double &x3b)
// Merksatz: cell oder deren Volumen schneidet oder beinhaltet komplette oder Teile der CuboidUmrandung
// returns true:
//  - cell cuts  cuboid3D
//  - cell boxes cuboid3D
// returns false:
//  - cell completely inside cuboid3D ( = cuboid3D boxes cell)
//  - cell und cuboid3D haben kein gemeinsames Volumen
{
    // erstmal die dumm Loesung
    if (!this->isCellInsideGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b) &&
        this->isCellInsideOrCuttingGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b)) {
        return true;
    }

    return false;

    // GbCuboid3D* cube = GbSystem3D::clipRectangle3D(*this->p1, *this->p2, x1a,x2a,x3a,x1b,x2b,x3b);
    // if(cube)
    //{
    //   cube->finalize();
    //   delete cube;
    //   return true;
    //}

    // return false;
}
/*=======================================================*/
bool GbCuboid3D::isCellInsideOrCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a,
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

    // if(   this->isCellInsideGbObject3D(x1a,x2a,x3a,x1b,x2b,x3b)
    //    || this->isCellCuttingGbObject3D(x1a,x2a,x3a,x1b,x2b,x3b) ) return true;
    //
    // return false;
}
/*=======================================================*/
vector<GbTriangle3D *> GbCuboid3D::getSurfaceTriangleSet()
{
    vector<GbTriangle3D *> triangles;
    GbPoint3D p1(getX1Minimum(), getX2Minimum(), getX3Minimum());
    GbPoint3D p2(getX1Maximum(), getX2Minimum(), getX3Minimum());
    GbPoint3D p3(getX1Maximum(), getX2Maximum(), getX3Minimum());
    GbPoint3D p4(getX1Minimum(), getX2Maximum(), getX3Minimum());
    GbPoint3D p5(getX1Minimum(), getX2Minimum(), getX3Maximum());
    GbPoint3D p6(getX1Maximum(), getX2Minimum(), getX3Maximum());
    GbPoint3D p7(getX1Maximum(), getX2Maximum(), getX3Maximum());
    GbPoint3D p8(getX1Minimum(), getX2Maximum(), getX3Maximum());

    GbPoint3D pUnten(getX1Centroid(), getX2Centroid(), getX3Minimum());
    GbPoint3D pOben(getX1Centroid(), getX2Centroid(), getX3Maximum());
    GbPoint3D pLinks(getX1Minimum(), getX2Centroid(), getX3Centroid());
    GbPoint3D pRechts(getX1Maximum(), getX2Centroid(), getX3Centroid());
    GbPoint3D pVorne(getX1Centroid(), getX2Minimum(), getX3Centroid());
    GbPoint3D pHinten(getX1Centroid(), getX2Maximum(), getX3Centroid());

    //"unten"
    triangles.push_back(new GbTriangle3D(new GbPoint3D(p2), new GbPoint3D(pUnten), new GbPoint3D(p3)));
    triangles.push_back(new GbTriangle3D(new GbPoint3D(p2), new GbPoint3D(p1), new GbPoint3D(pUnten)));
    triangles.push_back(new GbTriangle3D(new GbPoint3D(p3), new GbPoint3D(pUnten), new GbPoint3D(p4)));
    triangles.push_back(new GbTriangle3D(new GbPoint3D(p1), new GbPoint3D(p4), new GbPoint3D(pUnten)));
    //"oben"
    triangles.push_back(new GbTriangle3D(new GbPoint3D(p5), new GbPoint3D(p6), new GbPoint3D(pOben)));
    triangles.push_back(new GbTriangle3D(new GbPoint3D(p6), new GbPoint3D(p7), new GbPoint3D(pOben)));
    triangles.push_back(new GbTriangle3D(new GbPoint3D(p5), new GbPoint3D(pOben), new GbPoint3D(p8)));
    triangles.push_back(new GbTriangle3D(new GbPoint3D(pOben), new GbPoint3D(p7), new GbPoint3D(p8)));
    //"links"
    triangles.push_back(new GbTriangle3D(new GbPoint3D(p4), new GbPoint3D(p1), new GbPoint3D(pLinks)));
    triangles.push_back(new GbTriangle3D(new GbPoint3D(p4), new GbPoint3D(pLinks), new GbPoint3D(p8)));
    triangles.push_back(new GbTriangle3D(new GbPoint3D(p8), new GbPoint3D(pLinks), new GbPoint3D(p5)));
    triangles.push_back(new GbTriangle3D(new GbPoint3D(pLinks), new GbPoint3D(p1), new GbPoint3D(p5)));
    //"rechts"
    triangles.push_back(new GbTriangle3D(new GbPoint3D(p2), new GbPoint3D(p3), new GbPoint3D(pRechts)));
    triangles.push_back(new GbTriangle3D(new GbPoint3D(pRechts), new GbPoint3D(p3), new GbPoint3D(p7)));
    triangles.push_back(new GbTriangle3D(new GbPoint3D(p2), new GbPoint3D(pRechts), new GbPoint3D(p6)));
    triangles.push_back(new GbTriangle3D(new GbPoint3D(pRechts), new GbPoint3D(p7), new GbPoint3D(p6)));
    //"hinten"
    triangles.push_back(new GbTriangle3D(new GbPoint3D(p3), new GbPoint3D(p4), new GbPoint3D(pHinten)));
    triangles.push_back(new GbTriangle3D(new GbPoint3D(p3), new GbPoint3D(pHinten), new GbPoint3D(p7)));
    triangles.push_back(new GbTriangle3D(new GbPoint3D(p7), new GbPoint3D(pHinten), new GbPoint3D(p8)));
    triangles.push_back(new GbTriangle3D(new GbPoint3D(pHinten), new GbPoint3D(p4), new GbPoint3D(p8)));
    //"vorne"
    triangles.push_back(new GbTriangle3D(new GbPoint3D(p1), new GbPoint3D(p2), new GbPoint3D(pVorne)));
    triangles.push_back(new GbTriangle3D(new GbPoint3D(pVorne), new GbPoint3D(p2), new GbPoint3D(p6)));
    triangles.push_back(new GbTriangle3D(new GbPoint3D(p1), new GbPoint3D(pVorne), new GbPoint3D(p5)));
    triangles.push_back(new GbTriangle3D(new GbPoint3D(pVorne), new GbPoint3D(p6), new GbPoint3D(p5)));
    return triangles;
}
/*=======================================================*/
void GbCuboid3D::addSurfaceTriangleSet(vector<UbTupleFloat3> &nodes, vector<UbTupleInt3> &triangles)
{
    /*0*/ nodes.push_back(makeUbTuple((float)getX1Minimum(), (float)getX2Minimum(), (float)getX3Minimum()));
    /*1*/ nodes.push_back(makeUbTuple((float)getX1Maximum(), (float)getX2Minimum(), (float)getX3Minimum()));
    /*2*/ nodes.push_back(makeUbTuple((float)getX1Maximum(), (float)getX2Maximum(), (float)getX3Minimum()));
    /*3*/ nodes.push_back(makeUbTuple((float)getX1Minimum(), (float)getX2Maximum(), (float)getX3Minimum()));

    /*4*/ nodes.push_back(makeUbTuple((float)getX1Minimum(), (float)getX2Minimum(), (float)getX3Maximum()));
    /*5*/ nodes.push_back(makeUbTuple((float)getX1Maximum(), (float)getX2Minimum(), (float)getX3Maximum()));
    /*6*/ nodes.push_back(makeUbTuple((float)getX1Maximum(), (float)getX2Maximum(), (float)getX3Maximum()));
    /*7*/ nodes.push_back(makeUbTuple((float)getX1Minimum(), (float)getX2Maximum(), (float)getX3Maximum()));

    //"unten"
    triangles.push_back(makeUbTuple(0, 1, 2));
    triangles.push_back(makeUbTuple(0, 2, 3));
    //"oben"
    triangles.push_back(makeUbTuple(4, 5, 6));
    triangles.push_back(makeUbTuple(4, 6, 7));
    //"links"
    triangles.push_back(makeUbTuple(0, 3, 7));
    triangles.push_back(makeUbTuple(0, 7, 4));
    //"rechts"
    triangles.push_back(makeUbTuple(1, 2, 6));
    triangles.push_back(makeUbTuple(1, 6, 5));
    //"hinten"
    triangles.push_back(makeUbTuple(3, 2, 7));
    triangles.push_back(makeUbTuple(2, 7, 6));
    //"vorne"
    triangles.push_back(makeUbTuple(0, 1, 5));
    triangles.push_back(makeUbTuple(0, 5, 4));
}
/*=======================================================*/
string GbCuboid3D::toString()
{
    stringstream ss;
    ss << "GbCuboid3D[";
    ss << "p1=" << this->p1->toString();
    ss << ", p2=" << this->p2->toString();
    ss << "]";
    return ss.str();
}
/*=======================================================*/
GbPoint3D *GbCuboid3D::calculateInterSectionPoint3D(GbPoint3D & /*point1*/, GbPoint3D & /*point2*/)
{
    throw UbException(UB_EXARGS, "not correct implemented");
}
/*=======================================================*/
GbLine3D *GbCuboid3D::createClippedLine3D(GbPoint3D &point1, GbPoint3D &point2)
{
    return GbSystem3D::createClipLine3D(point1, point2, p1->getX1Coordinate(), p1->getX2Coordinate(),
                                        p1->getX3Coordinate(), p2->getX1Coordinate(), p2->getX2Coordinate(),
                                        p2->getX3Coordinate());
}
/*==========================================================*/
void GbCuboid3D::objectChanged(UbObservable *changedObject)
{
    GbPoint3D *point = dynamic_cast<GbPoint3D *>(changedObject);
    if (!point || (this->p1 != point && this->p2 != point))
        return;

    this->notifyObserversObjectChanged();
}
/*==========================================================*/
void GbCuboid3D::objectWillBeDeleted(UbObservable *objectForDeletion)
{
    if (this->p1) {
        UbObservable *observedObj = dynamic_cast<UbObservable *>(this->p1);
        if (objectForDeletion == observedObj) {
            this->p1 = NULL;
        }
    }
    if (this->p2) {
        UbObservable *observedObj = dynamic_cast<UbObservable *>(this->p2);
        if (objectForDeletion == observedObj) {
            this->p2 = NULL;
        }
    }
    // ACHTUNG: eigentlich muessten in allen methoden von GbLine if abfragen fuer NULL pointer hin... toDo
}
/*=======================================================*/
void GbCuboid3D::translate(const double &tx1, const double &tx2, const double &tx3)
{
    this->p1->translate(tx1, tx2, tx3);
    this->p2->translate(tx1, tx2, tx3);
    this->notifyObserversObjectChanged();
}
/*=======================================================*/
void GbCuboid3D::scale(const double &sx1, const double &sx2, const double &sx3)
{
    double lenX1 = this->getLengthX1();
    double lenX2 = this->getLengthX2();
    double lenX3 = this->getLengthX3();

    double deltaX1 = lenX1 * sx1 - lenX1;
    double deltaX2 = lenX2 * sx2 - lenX2;
    double deltaX3 = lenX3 * sx3 - lenX3;

    double p1X1 = this->p1->getX1Coordinate();
    double p1X2 = this->p1->getX2Coordinate();
    double p1X3 = this->p1->getX3Coordinate();

    double p2X1 = this->p2->getX1Coordinate();
    double p2X2 = this->p2->getX2Coordinate();
    double p2X3 = this->p2->getX3Coordinate();

    this->p1->setCoordinates(p1X1 - 0.5 * deltaX1, p1X2 - 0.5 * deltaX2, p1X3 - 0.5 * deltaX3);

    this->p2->setCoordinates(p2X1 + 0.5 * deltaX1, p2X2 + 0.5 * deltaX2, p2X3 + 0.5 * deltaX3);
}
/*==========================================================*/
double GbCuboid3D::getCellVolumeInsideGbObject3D(const double &x1a, const double &x2a, const double &x3a,
                                                 const double &x1b, const double &x2b, const double &x3b)
{
    if (this->isCellInsideGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b))
        return 1.0 * (x1b - x1a) * (x2b - x2a) * (x3b - x3a);
    if (!(this->isCellCuttingGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b)))
        return 0.0;

    GbCuboid3D *cube = GbSystem3D::clipRectangle3D(*this->p1, *this->p2, x1a, x2a, x3a, x1b, x2b, x3b);

    if (cube) {
        double eps;
        eps = (cube->getLengthX1()) * (cube->getLengthX2()) * (cube->getLengthX3());
        cube->finalize();
        delete cube;
        return eps;
    }
    return 0.0;
}
/*==========================================================*/
double GbCuboid3D::getIntersectionRaytraceFactor(const double &x1, const double &x2, const double &x3,
                                                 const double &rx1, const double &rx2, const double &rx3)
{
    double minB[3]   = { this->getX1Minimum(), this->getX2Minimum(), this->getX3Minimum() };
    double maxB[3]   = { this->getX1Maximum(), this->getX2Maximum(), this->getX3Maximum() };
    double origin[3] = { x1, x2, x3 };    // point
    double dir[3]    = { rx1, rx2, rx3 }; // ray

    bool inside = true;
    char quadrant[3];
    int whichPlane;
    double maxT[3];
    double candidatePlane[3];

    /* Find candidate planes; this loop can be avoided if
    rays cast all from the eye(assume perpsective view) */
    for (int i = 0; i < 3; i++) {
        if (origin[i] < minB[i]) {
            quadrant[i]       = 1 /*LEFT*/;
            candidatePlane[i] = minB[i];
            inside            = false;
        } else if (origin[i] > maxB[i]) {
            quadrant[i]       = 0 /*RIGHT*/;
            candidatePlane[i] = maxB[i];
            inside            = false;
        } else {
            quadrant[i] = 2 /*MIDDLE*/;
        }
    }
    /* Ray origin inside bounding box */
    if (inside) {
        // throw UbException(UB_EXARGS,"not done");
        return 0.0;
    }

    /* Calculate T distances to candidate planes */
    for (int i = 0; i < 3; i++) {
        if (quadrant[i] != 2 /*MIDDLE*/ && fabs(dir[i]) > 1.E-10) {
            maxT[i] = (candidatePlane[i] - origin[i]) / dir[i];
        } else
            maxT[i] = -1.0;
    }

    /* Get largest of the maxT's for final choice of intersection */
    whichPlane = 0;
    for (int i = 1; i < 3; i++)
        if (maxT[whichPlane] < maxT[i])
            whichPlane = i;

    /* Check final candidate actually inside box */
    if (maxT[whichPlane] < -1.E-10)
        return -1.0;
    double dummy;
    for (int i = 0; i < 3; i++) {
        if (whichPlane != i) {
            dummy = origin[i] + maxT[whichPlane] * dir[i];
            if (dummy < minB[i] || dummy > maxB[i])
                return -1.0;
        }
    }

    return maxT[whichPlane]; /* ray hits box */
}

//! \}
