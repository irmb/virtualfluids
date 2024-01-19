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
#include <basics/utilities/UbInfinity.h>
#include <geometry3d/GbCylinder3D.h>
#include <geometry3d/GbPoint3D.h>
#include <geometry3d/GbSystem3D.h>
#include <geometry3d/GbTriangle3D.h>

using namespace std;

// Konstruktor
/*==========================================================*/
GbCylinder3D::GbCylinder3D()

{
    this->setName("cylinder");
    GbPoint3D *p1 = new GbPoint3D();
    GbPoint3D *p2 = new GbPoint3D();
    mLine         = new GbLine3D(p1, p2);
    this->mLine->addObserver(this);
    mRad         = 0.0;
    cylinderType = GbCylinder3D::NOTPARALLELTOAXIS;
    this->mLine->addObserver(this);
    this->calculateValues();
}
/*=======================================================*/
GbCylinder3D::GbCylinder3D(GbCylinder3D *cylinder)
{
    this->setName("cylinder");
    mRad         = cylinder->getRadius();
    cylinderType = cylinder->cylinderType;
    mLine        = cylinder->getLine()->clone();

    this->mLine->addObserver(this);
    this->calculateValues();
}
/*==========================================================*/
GbCylinder3D::GbCylinder3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                           const double &x2b, const double &x3b, const double &rad)
{
    this->setName("cylinder");
    mLine = new GbLine3D;
    // Min/Max, damit gewaehrleistet ist, dass Startpunkt immer der "Achs-Minimale" ist
    // Anm.: bin nich tsicher ob weiter unten irgendwelche Algos drauf beruhen...
    //      geht nat nur solange, zylinder achs-parallel, aber das ist erzeit so!!!
    mLine->setPoints(new GbPoint3D(min(x1a, x1b), min(x2a, x2b), min(x3a, x3b)),
                     new GbPoint3D(max(x1a, x1b), max(x2a, x2b), max(x3a, x3b)));
    // mLine->setPoints( new GbPoint3D(x1a,x2a,x3a),new GbPoint3D(x1b, x2b ,x3b ));
    this->mLine->addObserver(this);
    mRad = fabs(rad);

    this->calculateValues();
}
/*==========================================================*/
GbCylinder3D::GbCylinder3D(GbPoint3D *p1, GbPoint3D *p2, const double &rad)
{
    this->setName("cylinder");
    mRad = rad;

    mLine = new GbLine3D(p1, p2);
    this->mLine->addObserver(this);
    this->calculateValues();
}
/*==========================================================*/
GbCylinder3D::GbCylinder3D(GbLine3D *line, const double &rad)
{
    this->setName("cylinder");
    mRad = rad;

    this->mLine = line;
    this->mLine->addObserver(this);

    this->calculateValues();
}
/*==========================================================*/
// Destruktor
GbCylinder3D::~GbCylinder3D()
{
    if (mLine)
        this->mLine->removeObserver(this);
    mLine = NULL;
}
/*=======================================================*/
void GbCylinder3D::calculateValues()
{
    double x1a = mLine->getPoint1()->x1;
    double x1b = mLine->getPoint2()->x1;
    double x2a = mLine->getPoint1()->x2;
    double x2b = mLine->getPoint2()->x2;
    double x3a = mLine->getPoint1()->x3;
    double x3b = mLine->getPoint2()->x3;

    if (x1a != x1b && x2a == x2b && x3a == x3b)
        this->cylinderType = X1PARALLEL;
    else if (x2a != x2b && x1a == x1b && x3a == x3b)
        this->cylinderType = X2PARALLEL;
    else if (x3a != x3b && x1a == x1b && x2a == x2b)
        this->cylinderType = X3PARALLEL;
    // nach dem serialisieren ruft er den Standardkonstruktor auf wo alles 0 ist und bricht sonst hier ab
    else if (x3a == x3b && x1a == x1b && x2a == x2b)
        this->cylinderType = X1PARALLEL;
    else
        this->cylinderType = NOTPARALLELTOAXIS;

    if ((this->cylinderType & NOTPARALLELTOAXIS) == NOTPARALLELTOAXIS)
        throw UbException(UB_EXARGS,
                          "derzeit nur zu Achsen orthogonale Zylinder erlaubt... isPointInObject3D funzt sonst ned");

    if (this->isParallelToX1Axis()) {
        minX1 = mLine->getX1Minimum();
        maxX1 = mLine->getX1Maximum();
        minX2 = mLine->getX2Centroid() - mRad;
        maxX2 = mLine->getX2Centroid() + mRad;
        minX3 = mLine->getX3Centroid() - mRad;
        maxX3 = mLine->getX3Centroid() + mRad;
    } else if (this->isParallelToX2Axis()) {
        minX1 = mLine->getX1Centroid() - mRad;
        maxX1 = mLine->getX1Centroid() + mRad;
        minX2 = mLine->getX2Minimum();
        maxX2 = mLine->getX2Maximum();
        minX3 = mLine->getX3Centroid() - mRad;
        maxX3 = mLine->getX3Centroid() + mRad;
    } else if (this->isParallelToX3Axis()) {
        minX1 = mLine->getX1Centroid() - mRad;
        maxX1 = mLine->getX1Centroid() + mRad;
        minX2 = mLine->getX2Centroid() - mRad;
        maxX2 = mLine->getX2Centroid() + mRad;
        minX3 = mLine->getX3Minimum();
        maxX3 = mLine->getX3Maximum();
    }

    centerX1 = mLine->getX1Centroid();
    centerX2 = mLine->getX2Centroid();
    centerX3 = mLine->getX3Centroid();
}

/*=======================================================*/
void GbCylinder3D::finalize()
{
    if (this->mLine) {
        mLine->finalize();
        delete mLine;
        mLine = NULL;
    }
}
/*=======================================================*/
double GbCylinder3D::getHeight()
{
    if (mLine)
        return mLine->getLength();

    return 0.0;
}
/*=======================================================*/
GbPoint3D *GbCylinder3D::getPoint1()
{
    if (this->mLine)
        return this->mLine->getPoint1();
    return NULL;
}
/*=======================================================*/
GbPoint3D *GbCylinder3D::getPoint2()
{
    if (this->mLine)
        return this->mLine->getPoint2();
    return NULL;
}
/*=======================================================*/
void GbCylinder3D::setRadius(const double &radius)
{
    this->mRad = std::fabs(radius);
    this->notifyObserversObjectChanged();
}
/*=======================================================*/
void GbCylinder3D::setLine(GbLine3D *line)
{
    if (this->mLine)
        this->mLine->removeObserver(this);
    this->mLine = line;
    this->mLine->addObserver(this);
    this->calculateValues();

    this->notifyObserversObjectChanged();
}
/*=======================================================*/
void GbCylinder3D::setPoint1(const double &x1, const double &x2, const double &x3)
{
    if (!mLine->getPoint1())
        throw UbException(UB_EXARGS, "line has no point1");
    mLine->getPoint1()->setCoordinates(x1, x2, x3);
    this->calculateValues();

    // this->notifyObserversObjectChanged(); //wird automatisch aufgerufen, da der point (this) benachrichtigt...
}
/*=======================================================*/
void GbCylinder3D::setPoint2(const double &x1, const double &x2, const double &x3)
{
    if (!mLine->getPoint2())
        throw UbException(UB_EXARGS, "line has no point2");
    mLine->getPoint2()->setCoordinates(x1, x2, x3);
    this->calculateValues();

    // this->notifyObserversObjectChanged(); //wird automatisch aufgerufen, da der point (this) benachrichtigt...
}
/*==========================================================*/
bool GbCylinder3D::isPointInGbObject3D(const double &x1p, const double &x2p, const double &x3p)
{
    // true, wenn 'in Object' oder 'auf Boundary'!
    if (this->isParallelToX1Axis() && (UbMath::less(x1p, minX1) || UbMath::greater(x1p, maxX1)))
        return false;
    else if (this->isParallelToX2Axis() && (UbMath::less(x2p, minX2) || UbMath::greater(x2p, maxX2)))
        return false;
    else if (this->isParallelToX3Axis() && (UbMath::less(x3p, minX3) || UbMath::greater(x3p, maxX3)))
        return false;
    else if (this->isNotParallelToAxis())
        throw UbException(UB_EXARGS,
                          "derzeit nur zu Achsen orthogonale Zylinder erlaubt... isPointInObject3D funzt sonst ned");

    return UbMath::lessEqual(fabs(mLine->getDistance(x1p, x2p, x3p)), fabs(mRad));
}
/*==========================================================*/
bool GbCylinder3D::isPointInGbObject3D(const double &x1p, const double &x2p, const double &x3p, bool &pointIsOnBoundary)
{
    // funzt derzeit nur bei achsparallelen cylindern
    pointIsOnBoundary = false;

    if (this->isParallelToX1Axis() && (UbMath::less(x1p, minX1) || UbMath::greater(x1p, maxX1)))
        return false;
    else if (this->isParallelToX2Axis() && (UbMath::less(x2p, minX2) || UbMath::greater(x2p, maxX2)))
        return false;
    else if (this->isParallelToX3Axis() && (UbMath::less(x3p, minX3) || UbMath::greater(x3p, maxX3)))
        return false;
    else if (this->isNotParallelToAxis())
        throw UbException(UB_EXARGS,
                          "derzeit nur zu Achsen orthogonale Zylinder erlaubt... isPointInObject3D funzt sonst ned");

    // true, wenn 'in Object' oder 'auf Boundary'!

    double dis = mLine->getDistance(x1p, x2p, x3p);

    if (UbMath::equal(dis, mRad))
        pointIsOnBoundary = true;

    if (this->isParallelToX1Axis() && (UbMath::equal(x1p, minX1) || UbMath::equal(x1p, maxX1)))
        pointIsOnBoundary = true;
    else if (this->isParallelToX2Axis() && (UbMath::equal(x2p, minX2) || UbMath::equal(x2p, maxX2)))
        pointIsOnBoundary = true;
    else if (this->isParallelToX3Axis() && (UbMath::equal(x3p, minX3) || UbMath::equal(x3p, maxX3)))
        pointIsOnBoundary = true;

    return UbMath::lessEqual(dis, mRad);
}
/*==========================================================*/
string GbCylinder3D::toString()
{
    stringstream ss;
    ss << "GbCylinder3D[";
    ss << "line=" << this->mLine->toString();
    ss << ", r=" << this->mRad;
    ss << "]";
    return (ss.str());
}
/*=======================================================*/
bool GbCylinder3D::isCellInsideGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                          const double &x2b, const double &x3b)
{
    if (this->isPointInGbObject3D(x1a, x2a, x3a) && this->isPointInGbObject3D(x1b, x2a, x3a) &&
        this->isPointInGbObject3D(x1b, x2b, x3a) && this->isPointInGbObject3D(x1a, x2b, x3a) &&
        this->isPointInGbObject3D(x1a, x2a, x3b) && this->isPointInGbObject3D(x1b, x2a, x3b) &&
        this->isPointInGbObject3D(x1b, x2b, x3b) && this->isPointInGbObject3D(x1a, x2b, x3b)) {
        return true;
    }
    return false;
}
/*==========================================================*/
bool GbCylinder3D::isCellCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                           const double &x2b, const double &x3b)
// Merksatz: cell oder deren Volumen schneidet oder beinhaltet komplette oder Teile der CuboidUmrandung
// returns true:
//  - cell cuts  cylinder3D
//  - cell boxes cylinder3D
// returns false:
//  - cell completely inside cylinder3D ( = cylinder3D boxes cell)
//  - cell und cylinder3D haben kein gemeinsames Volumen
{
    // erstmal wieder die dumm Loesung
    if (this->isCellInsideOrCuttingGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b) &&
        !this->isCellInsideGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b)) {
        return true;
    }

    return false;
}
/*==========================================================*/
bool GbCylinder3D::isCellInsideOrCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a,
                                                   const double &x1b, const double &x2b, const double &x3b)
// returns true:
//  - cell completely inside cylinder3D ( = cylinder3D boxes cell)
//  - cell cuts  cylinder3D
//  - cell boxes cylinder3D
// returns false:
//  - cell und cylinder3D haben kein gemeinsames Volumen
{
    double dmin = 0.0;

    if (this->isParallelToX1Axis()) {
        // check liegt Cell komplett !x1-ausserhalb"?
        if (UbMath::less(x1a, minX1) && UbMath::less(x1b, minX1))
            return false;
        if (UbMath::greater(x1a, maxX1) && UbMath::greater(x1b, maxX1))
            return false;

        // mittelpunkt kreis-querschnitt
        double &midX2 = mLine->getPoint1()->x2;
        double &midX3 = mLine->getPoint1()->x3;
        if (UbMath::less(midX2, x2a))
            dmin += std::pow(midX2 - x2a, 2.0);
        else if (UbMath::greater(midX2, x2b))
            dmin += std::pow(midX2 - x2b, 2.0);
        if (UbMath::less(midX3, x3a))
            dmin += std::pow(midX3 - x3a, 2.0);
        else if (UbMath::greater(midX3, x3b))
            dmin += std::pow(midX3 - x3b, 2.0);
        if (UbMath::lessEqual(dmin, mRad * mRad))
            return true;

        return false;
    } else if (this->isParallelToX2Axis()) {
        // check liegt Cell komplett !x2-ausserhalb"?
        if (UbMath::less(x2a, minX2) && UbMath::less(x2b, minX2))
            return false;
        if (UbMath::greater(x2a, maxX2) && UbMath::greater(x2b, maxX2))
            return false;

        // mittelpunkt kreis-querschnitt
        double &midX1 = mLine->getPoint1()->x1;
        double &midX3 = mLine->getPoint1()->x3;
        if (UbMath::less(midX1, x1a))
            dmin += std::pow(midX1 - x1a, 2.0);
        else if (UbMath::greater(midX1, x1b))
            dmin += std::pow(midX1 - x1b, 2.0);
        if (UbMath::less(midX3, x3a))
            dmin += std::pow(midX3 - x3a, 2.0);
        else if (UbMath::greater(midX3, x3b))
            dmin += std::pow(midX3 - x3b, 2.0);
        if (UbMath::lessEqual(dmin, mRad * mRad))
            return true;

    } else if (this->isParallelToX3Axis()) {
        // check liegt Cell komplett !x3-ausserhalb"?
        if (UbMath::less(x3a, minX3) && UbMath::less(x3b, minX3))
            return false;
        if (UbMath::greater(x3a, maxX3) && UbMath::greater(x3b, maxX3))
            return false;

        // mittelpunkt kreis-querschnitt
        double &midX1 = mLine->getPoint1()->x1;
        double &midX2 = mLine->getPoint1()->x2;
        if (UbMath::less(midX1, x1a))
            dmin += std::pow(midX1 - x1a, 2.0);
        else if (UbMath::greater(midX1, x1b))
            dmin += std::pow(midX1 - x1b, 2.0);
        if (UbMath::less(midX2, x2a))
            dmin += std::pow(midX2 - x2a, 2.0);
        else if (UbMath::greater(midX2, x2b))
            dmin += std::pow(midX2 - x2b, 2.0);
        if (UbMath::lessEqual(dmin, mRad * mRad))
            return true;
    }

    return false;
}
/*==========================================================*/
GbLine3D *GbCylinder3D::createClippedLine3D(GbPoint3D &point1, GbPoint3D &point2)
{
    // liefert immer "innere" linie, also der teil, der vom Zylinder "abgeschnitten" wurde!
    // funktioniert derzeit nur mit achsenparallelen Zylindern!
    vector<GbPoint3D *> schnittpunkte;

    double xa, ya, za, xb, yb, zb, xm, ym, zStart, zEnd, t1, t2;

    if (this->isParallelToX1Axis()) {
        xa     = point1.getX2Coordinate();
        ya     = point1.getX3Coordinate();
        za     = point1.getX1Coordinate();
        xb     = point2.getX2Coordinate();
        yb     = point2.getX3Coordinate();
        zb     = point2.getX1Coordinate();
        xm     = mLine->getPoint1()->getX2Coordinate();
        ym     = mLine->getPoint1()->getX3Coordinate();
        zStart = mLine->getPoint1()->getX1Coordinate();
        zEnd   = mLine->getPoint2()->getX1Coordinate();
    } else if (this->isParallelToX2Axis()) {
        xa     = point1.getX1Coordinate();
        ya     = point1.getX3Coordinate();
        za     = point1.getX2Coordinate();
        xb     = point2.getX1Coordinate();
        yb     = point2.getX3Coordinate();
        zb     = point2.getX2Coordinate();
        xm     = mLine->getPoint1()->getX1Coordinate();
        ym     = mLine->getPoint1()->getX3Coordinate();
        zStart = mLine->getPoint1()->getX2Coordinate();
        zEnd   = mLine->getPoint2()->getX2Coordinate();
    } else if (this->isParallelToX3Axis()) {
        xa     = point1.getX1Coordinate();
        ya     = point1.getX2Coordinate();
        za     = point1.getX3Coordinate();
        xb     = point2.getX1Coordinate();
        yb     = point2.getX2Coordinate();
        zb     = point2.getX3Coordinate();
        xm     = mLine->getPoint1()->getX1Coordinate();
        ym     = mLine->getPoint1()->getX2Coordinate();
        zStart = mLine->getPoint1()->getX3Coordinate();
        zEnd   = mLine->getPoint2()->getX3Coordinate();
    } else
        throw UbException(UB_EXARGS, "funktioniert derzeit nur mit achsenparallelen Zylindern");

    // Bestimmung des Schnittpunktes mit unendlich ausgedehntem Zylinder
    double r   = mRad;
    double r2  = r * r;
    double xa2 = xa * xa;
    double xb2 = xb * xb;
    double ya2 = ya * ya;
    double yb2 = yb * yb;
    double xm2 = xm * xm;
    double ym2 = ym * ym;

    double wurzel = 2.0 * xa * xm * yb2 + 2.0 * ya * ym * xb2 - 2.0 * xa * xb * r2 + 2.0 * xa * xb * ym2 -
                    2.0 * ya * yb * r2 + 2.0 * xa2 * yb * ym + 2.0 * xa * xm * ya * ym - 2.0 * xa * xm * yb * ym -
                    2.0 * ya * ym * xb * xm + 2.0 * xb * xm * yb * ym + 2.0 * ya * yb * xa * xb -
                    2.0 * ya * yb * xa * xm - 2.0 * ya * yb * xb * xm - 2.0 * xa * xb * ya * ym -
                    2.0 * xa * xb * yb * ym + 2.0 * xb * xm * ya2 + 2.0 * ya * yb * xm2 - xa2 * yb2 - xb2 * ya2 +
                    xa2 * r2 - xa2 * ym2 + xb2 * r2 - xb2 * ym2 + ya2 * r2 - ya2 * xm2 + yb2 * r2 - yb2 * xm2;
    double nenner  = -2.0 * (ya * yb + xa * xb) + xa2 + xb2 + ya2 + yb2;
    double zaehler = 2.0 * (-xa * xm + xb * xm - ya * ym + yb * ym) + xa2 - xb2 + ya2 - yb2;

    if (UbMath::greaterEqual(wurzel, 0.0) && !UbMath::zero(nenner)) // fabs(nenner)>1.E-13)
    {
        t1 = (zaehler + 2.0 * sqrt(wurzel)) / nenner;
        t2 = (zaehler - 2.0 * sqrt(wurzel)) / nenner;

        if (UbMath::inClosedInterval(t1, -1.0, 1.0)) // Schnittpunkt innerhalb der Strecke
        {
            double x = xa * (0.5 - 0.5 * t1) + xb * (0.5 + 0.5 * t1);
            double y = ya * (0.5 - 0.5 * t1) + yb * (0.5 + 0.5 * t1);
            double z = za * (0.5 - 0.5 * t1) + zb * (0.5 + 0.5 * t1);

            if (UbMath::inClosedInterval(z, zStart, zEnd)) // zWert muss sich innerhal der cylinderlaenge befinden
            {
                if (this->isParallelToX1Axis())
                    schnittpunkte.push_back(new GbPoint3D(z, x, y));
                else if (this->isParallelToX2Axis())
                    schnittpunkte.push_back(new GbPoint3D(x, z, y));
                else if (this->isParallelToX3Axis())
                    schnittpunkte.push_back(new GbPoint3D(x, y, z));
            }
        }
        if (fabs(t2 - t1) > 1.E-13 && UbMath::inClosedInterval(t2, -1.0, 1.0)) // Schnittpunkt innerhalb der Strecke
        {
            double x = xa * (0.5 - 0.5 * t2) + xb * (0.5 + 0.5 * t2);
            double y = ya * (0.5 - 0.5 * t2) + yb * (0.5 + 0.5 * t2);
            double z = za * (0.5 - 0.5 * t2) + zb * (0.5 + 0.5 * t2);

            if (UbMath::inClosedInterval(z, zStart, zEnd)) // zWert muss sich innerhal der cylinderlaenge befinden
            {
                if (this->isParallelToX1Axis())
                    schnittpunkte.push_back(new GbPoint3D(z, x, y));
                else if (this->isParallelToX2Axis())
                    schnittpunkte.push_back(new GbPoint3D(x, z, y));
                else if (this->isParallelToX3Axis())
                    schnittpunkte.push_back(new GbPoint3D(x, y, z));
            }
        }
    }
    // wenn nenner==0 -> Strecke parallel zu Zylinder! Es muss noch auf Schnittpunkt mit "Deckeln" geprueft werden

    // Schnittpunkt mit Seitenflaechen bestimmen
    // hierzu wird der schnittpunkt der gegebnen strecke mit den seitenflaechenberechnet
    // als erstes "schaut man seitlich auf den Zylinder" --> kreisflaechen wird als strecke darsgestellt
    // mit diesen "strecken" berechnet man Schnittpunkte.
    // anschliessend wird geprueft, ob der berechnete Schnittpunkt ueberhaupt im kreis liegt
    // falls ja --> Schnittpunkt vorhanden

    double x1a, y1a, z1a, x1b, y1b, z1b, // uebergebene Strecke
        x2a, y2a, x2b, y2b,              // erste "Kreisstrecke"
        x3a, y3a, x3b, y3b,              // zweite "Kreisstrecke"
        y2m, /*z2m,*/ y3m, z3m;
    double nenner1ab;

    if (this->isParallelToX1Axis()) {
        x1a = point1.getX1Coordinate();
        y1a = point1.getX2Coordinate();
        z1a = point1.getX3Coordinate();
        x1b = point2.getX1Coordinate();
        y1b = point2.getX2Coordinate();
        z1b = point2.getX3Coordinate();

        x2a = mLine->getPoint1()->getX1Coordinate();
        y2m = mLine->getPoint1()->getX2Coordinate();
        //      z2m=mLine->getPoint1()->getX3Coordinate();
        y2a = y2m + mRad;
        x2b = mLine->getPoint1()->getX1Coordinate();
        y2b = y2m - mRad;

        x3a = mLine->getPoint2()->getX1Coordinate(); //
        y3m = mLine->getPoint2()->getX2Coordinate();
        z3m = mLine->getPoint2()->getX3Coordinate();
        y3a = y3m + mRad;
        x3b = mLine->getPoint2()->getX1Coordinate();
        y3b = y3m - mRad;
    } else if (this->isParallelToX2Axis()) {
        x1a = point1.getX2Coordinate();
        y1a = point1.getX3Coordinate();
        z1a = point1.getX1Coordinate();
        x1b = point2.getX2Coordinate();
        y1b = point2.getX3Coordinate();
        z1b = point2.getX1Coordinate();

        x2a = mLine->getPoint1()->getX2Coordinate();
        y2m = mLine->getPoint1()->getX3Coordinate();
        //      z2m=mLine->getPoint1()->getX1Coordinate();
        y2a = y2m + mRad;
        x2b = mLine->getPoint1()->getX2Coordinate();
        y2b = y2m - mRad;

        x3a = mLine->getPoint2()->getX2Coordinate(); //
        y3m = mLine->getPoint2()->getX3Coordinate();
        z3m = mLine->getPoint2()->getX1Coordinate();
        y3a = y3m + mRad;
        x3b = mLine->getPoint2()->getX2Coordinate();
        y3b = y3m - mRad;
    } else if (this->isParallelToX3Axis()) {
        x1a = point1.getX3Coordinate();
        y1a = point1.getX2Coordinate();
        z1a = point1.getX1Coordinate();
        x1b = point2.getX3Coordinate();
        y1b = point2.getX2Coordinate();
        z1b = point2.getX1Coordinate();

        x2a = mLine->getPoint1()->getX3Coordinate();
        y2m = mLine->getPoint1()->getX2Coordinate();
        //      z2m=mLine->getPoint1()->getX1Coordinate();
        y2a = y2m + mRad;
        x2b = mLine->getPoint1()->getX3Coordinate();
        y2b = y2m - mRad;

        x3a = mLine->getPoint2()->getX3Coordinate(); //
        y3m = mLine->getPoint2()->getX2Coordinate();
        z3m = mLine->getPoint2()->getX1Coordinate();
        y3a = y3m + mRad;
        x3b = mLine->getPoint2()->getX3Coordinate();
        y3b = y3m - mRad;
    } else
        throw UbException(UB_EXARGS, "funktioniert derzeit nur mit achsenparallelen Zylindern");

    nenner1ab = -y1a * x2a + y1a * x2b + y1b * x2a - y1b * x2b + x1a * y2a - x1a * y2b - x1b * y2a + x1b * y2b;
    // double nenner2 = x1a*y2a-x1a*y2b-x1b*y2a+x1b*y2b-y1a*x2a+y1a*x2b+y1b*x2a-y1b*x2b;
    if (fabs(nenner1ab) > 1.E-13) // andernfalls sind die beiden Strecken parallel
    {
        // tStrecke ist fuer gegebene Strecke!
        double t1ab = (-y1a * x2a + y1a * x2b - 2.0 * y2a * x2b + x1a * y2a - x1a * y2b - x1b * y2b + 2.0 * y2b * x2a +
                       x1b * y2a - y1b * x2a + y1b * x2b) /
                      nenner1ab;
        // double tStrecke =
        // -(-x1a*y2a+x1a*y2b+2.0*y2a*x2b+y1a*x2a-2.0*x2a*y2b-y1a*x2b+y1b*x2a-y1b*x2b-x1b*y2a+x1b*y2b)/nenner2; wenn -1
        // <= t2 <= +1 -> SP mit strecke
        if (UbMath::inClosedInterval(t1ab, -1.0, 1.0)) // Schnittpunkt innerhalb der Strecke
        {
            double x, y, z, abstand_ist;
            if (this->isParallelToX1Axis()) {
                x           = x1a * (0.5 - 0.5 * t1ab) + x1b * (0.5 + 0.5 * t1ab);
                y           = y1a * (0.5 - 0.5 * t1ab) + y1b * (0.5 + 0.5 * t1ab);
                z           = z1a * (0.5 - 0.5 * t1ab) + z1b * (0.5 + 0.5 * t1ab);
                abstand_ist = sqrt((y3m - y) * (y3m - y) + (z3m - z) * (z3m - z));
            } else if (this->isParallelToX2Axis()) {
                y           = x1a * (0.5 - 0.5 * t1ab) + x1b * (0.5 + 0.5 * t1ab);
                z           = y1a * (0.5 - 0.5 * t1ab) + y1b * (0.5 + 0.5 * t1ab);
                x           = z1a * (0.5 - 0.5 * t1ab) + z1b * (0.5 + 0.5 * t1ab);
                abstand_ist = sqrt((y3m - z) * (y3m - z) + (z3m - x) * (z3m - x));
            } else if (this->isParallelToX3Axis()) {
                z           = x1a * (0.5 - 0.5 * t1ab) + x1b * (0.5 + 0.5 * t1ab);
                y           = y1a * (0.5 - 0.5 * t1ab) + y1b * (0.5 + 0.5 * t1ab);
                x           = z1a * (0.5 - 0.5 * t1ab) + z1b * (0.5 + 0.5 * t1ab);
                abstand_ist = sqrt((y3m - y) * (y3m - y) + (z3m - x) * (z3m - x));
            } else
                throw UbException(UB_EXARGS, "funktioniert derzeit nur mit achsenparallelen Zylindern");

            // pruefen, ob Punkt Element von Kreisflaeche
            // double abstand_ist=sqrt((y2m-y)*(y2m-y)+(z2m-z)*(z2m-z));
            if (UbMath::lessEqual(abstand_ist, mRad)) // Punkt ist Schnittpunkt
            {
                bool exists = false;
                for (int pos = 0; pos < (int)schnittpunkte.size(); ++pos) {
                    if (fabs(schnittpunkte[pos]->getX1Coordinate() - x) < 1.E-13 &&
                        fabs(schnittpunkte[pos]->getX2Coordinate() - y) < 1.E-13 &&
                        fabs(schnittpunkte[pos]->getX3Coordinate() - z) < 1.E-13)
                        exists = true;
                }

                if (!exists)
                    schnittpunkte.push_back(new GbPoint3D(x, y, z));
            }
        }
    }

    nenner1ab = -y1a * x3a + y1a * x3b + y1b * x3a - y1b * x3b + x1a * y3a - x1a * y3b - x1b * y3a + x1b * y3b;

    if (fabs(nenner1ab) > 1.E-13) // andernfalls sind die beiden Strecken parallel
    {
        // tStrecke ist fuer gegebene Strecke!
        double t1ab = (-y1a * x3a + y1a * x3b - x1b * y3b - 2.0 * y3a * x3b - x1a * y3b + 2.0 * y3b * x3a + x1a * y3a +
                       x1b * y3a - y1b * x3a + y1b * x3b) /
                      nenner1ab;

        if (UbMath::inClosedInterval(t1ab, -1.0, 1.0)) // Schnittpunkt innerhalb der Strecke
        {
            double x, y, z, abstand_ist;
            if (this->isParallelToX1Axis()) {
                x           = x1a * (0.5 - 0.5 * t1ab) + x1b * (0.5 + 0.5 * t1ab);
                y           = y1a * (0.5 - 0.5 * t1ab) + y1b * (0.5 + 0.5 * t1ab);
                z           = z1a * (0.5 - 0.5 * t1ab) + z1b * (0.5 + 0.5 * t1ab);
                abstand_ist = sqrt((y3m - y) * (y3m - y) + (z3m - z) * (z3m - z));
            } else if (this->isParallelToX2Axis()) {
                y           = x1a * (0.5 - 0.5 * t1ab) + x1b * (0.5 + 0.5 * t1ab);
                z           = y1a * (0.5 - 0.5 * t1ab) + y1b * (0.5 + 0.5 * t1ab);
                x           = z1a * (0.5 - 0.5 * t1ab) + z1b * (0.5 + 0.5 * t1ab);
                abstand_ist = sqrt((y3m - z) * (y3m - z) + (z3m - x) * (z3m - x));
            } else if (this->isParallelToX3Axis()) {
                z           = x1a * (0.5 - 0.5 * t1ab) + x1b * (0.5 + 0.5 * t1ab);
                y           = y1a * (0.5 - 0.5 * t1ab) + y1b * (0.5 + 0.5 * t1ab);
                x           = z1a * (0.5 - 0.5 * t1ab) + z1b * (0.5 + 0.5 * t1ab);
                abstand_ist = sqrt((y3m - y) * (y3m - y) + (z3m - x) * (z3m - x));
            } else
                throw UbException(UB_EXARGS, "cylinder must be parallel to one axis");

            // pruefen, ob Punkt Element von Kreisflaeche
            // double abstand_ist=sqrt((y2m-y)*(y2m-y)+(z2m-z)*(z2m-z));

            if (UbMath::lessEqual(abstand_ist, mRad)) // Punkt ist Schnittpunkt
            {
                bool exists = false;
                for (int pos = 0; pos < (int)schnittpunkte.size(); ++pos) {
                    if (fabs(schnittpunkte[pos]->getX1Coordinate() - x) < 1.E-13 &&
                        fabs(schnittpunkte[pos]->getX2Coordinate() - y) < 1.E-13 &&
                        fabs(schnittpunkte[pos]->getX3Coordinate() - z) < 1.E-13)
                        exists = true;
                }

                if (!exists)
                    schnittpunkte.push_back(new GbPoint3D(x, y, z));
            }
        }
    }

    int nofSchnittpunkte = (int)schnittpunkte.size();
    if (nofSchnittpunkte == 0)
        return NULL;
    else if (nofSchnittpunkte > 2)
        throw UbException(UB_EXARGS, "more than three intersection points - not possible");
    else if (nofSchnittpunkte == 2)
        return new GbLine3D(schnittpunkte[0], schnittpunkte[1]);
    else if (nofSchnittpunkte == 1) {
        if (this->isPointInGbObject3D(&point1))
            return new GbLine3D(schnittpunkte[0], new GbPoint3D(point1));
        else if (this->isPointInGbObject3D(&point2))
            return new GbLine3D(schnittpunkte[0], new GbPoint3D(point2));
        else
            return new GbLine3D(
                schnittpunkte[0],
                new GbPoint3D(*(schnittpunkte[0]))); // strecke beruehrt clippedLine reduziert sich auf einen Punkt!!!
    }

    return NULL;
}
/*==========================================================*/
vector<GbTriangle3D *> GbCylinder3D::getSurfaceTriangleSet()
{
    double x1ma, x1mb, x2m, x3m;
    if (this->isParallelToX1Axis()) {
        x1ma = this->getX1Minimum();
        x1mb = this->getX1Maximum();
        x2m  = this->getX2Centroid();
        x3m  = this->getX3Centroid();
    } else if (this->isParallelToX2Axis()) {
        x1ma = this->getX2Minimum();
        x1mb = this->getX2Maximum();
        x2m  = this->getX1Centroid();
        x3m  = this->getX3Centroid();
    } else if (this->isParallelToX3Axis()) {
        x1ma = this->getX3Minimum();
        x1mb = this->getX3Maximum();
        x2m  = this->getX2Centroid();
        x3m  = this->getX1Centroid();
    } else
        throw UbException(UB_EXARGS, "cylinder not axis prallel");

    vector<GbTriangle3D *> triangles;

    int segmentsCircle = 20;
    double deltaPhi    = UbMath::PI / (double)segmentsCircle;

    double phiX1a, phiX1b;
    double x1a, x2a, x3a, x1b, x2b, x3b, x1c, x2c, x3c, x1d, x2d, x3d;

    double dXCylinder    = fabs((x1mb - x1ma)) / (double)segmentsCircle;
    int segmentsCylinder = (int)(fabs(x1mb - x1ma) / dXCylinder);
    for (int segCyl = 0; segCyl < segmentsCylinder; segCyl++) {
        x1a = x1d = x1ma + segCyl * dXCylinder;
        x1b = x1c = x1a + dXCylinder;

        for (phiX1a = 2.0 * UbMath::PI; phiX1a > 0; phiX1a -= deltaPhi) {
            phiX1b = phiX1a + deltaPhi;

            x2a = x2m + mRad * std::sin(phiX1a);
            x3a = x3m + mRad * std::cos(phiX1a);
            x2b = x2m + mRad * std::sin(phiX1b);
            x3b = x3m + mRad * std::cos(phiX1b);

            if (this->isParallelToX1Axis()) {
                triangles.push_back(new GbTriangle3D(new GbPoint3D(x1b, x2b, x3b), new GbPoint3D(x1b, x2a, x3a),
                                                     new GbPoint3D(x1a, x2a, x3a)));
                triangles.push_back(new GbTriangle3D(new GbPoint3D(x1a, x2a, x3a), new GbPoint3D(x1a, x2b, x3b),
                                                     new GbPoint3D(x1b, x2b, x3b)));
            } else if (this->isParallelToX2Axis()) {
                triangles.push_back(new GbTriangle3D(new GbPoint3D(x2b, x1b, x3b), new GbPoint3D(x2a, x1b, x3a),
                                                     new GbPoint3D(x2a, x1a, x3a)));
                triangles.push_back(new GbTriangle3D(new GbPoint3D(x2a, x1a, x3a), new GbPoint3D(x2b, x1a, x3b),
                                                     new GbPoint3D(x2b, x1b, x3b)));
            } else if (this->isParallelToX3Axis()) {
                triangles.push_back(new GbTriangle3D(new GbPoint3D(x3b, x2b, x1b), new GbPoint3D(x3a, x2a, x1b),
                                                     new GbPoint3D(x3a, x2a, x1a)));
                triangles.push_back(new GbTriangle3D(new GbPoint3D(x3a, x2a, x1a), new GbPoint3D(x3b, x2b, x1a),
                                                     new GbPoint3D(x3b, x2b, x1b)));
            }
        }
    }

    int segmentsSide = (int)(mRad / dXCylinder);
    double radius0, radius1;
    for (int segCyl = 0; segCyl < segmentsSide; segCyl++) {
        radius0 = segCyl * dXCylinder;
        radius1 = radius0 + dXCylinder;
        if (segCyl == segmentsSide - 1)
            radius1 = mRad;

        for (phiX1a = 2.0 * UbMath::PI; phiX1a > 0; phiX1a -= deltaPhi) {
            phiX1b = phiX1a + deltaPhi;

            x2a = x2m + radius0 * std::sin(phiX1a);
            x3a = x3m + radius0 * std::cos(phiX1a);
            x2b = x2m + radius0 * std::sin(phiX1b);
            x3b = x3m + radius0 * std::cos(phiX1b);
            x2c = x2m + radius1 * std::sin(phiX1b);
            x3c = x3m + radius1 * std::cos(phiX1b);
            x2d = x2m + radius1 * std::sin(phiX1a);
            x3d = x3m + radius1 * std::cos(phiX1a);

            if (this->isParallelToX1Axis()) {
                triangles.push_back(new GbTriangle3D(new GbPoint3D(x1ma, x2a, x3a), new GbPoint3D(x1ma, x2b, x3b),
                                                     new GbPoint3D(x1ma, x2c, x3c)));
                triangles.push_back(new GbTriangle3D(new GbPoint3D(x1ma, x2c, x3c), new GbPoint3D(x1ma, x2d, x3d),
                                                     new GbPoint3D(x1ma, x2a, x3a)));
                triangles.push_back(new GbTriangle3D(new GbPoint3D(x1mb, x2c, x3c), new GbPoint3D(x1mb, x2b, x3b),
                                                     new GbPoint3D(x1mb, x2a, x3a)));
                triangles.push_back(new GbTriangle3D(new GbPoint3D(x1mb, x2a, x3a), new GbPoint3D(x1mb, x2d, x3d),
                                                     new GbPoint3D(x1mb, x2c, x3c)));
            } else if (this->isParallelToX2Axis()) {
                triangles.push_back(new GbTriangle3D(new GbPoint3D(x2a, x1ma, x3a), new GbPoint3D(x2b, x1ma, x3b),
                                                     new GbPoint3D(x2c, x1ma, x3c)));
                triangles.push_back(new GbTriangle3D(new GbPoint3D(x2c, x1ma, x3c), new GbPoint3D(x2d, x1ma, x3d),
                                                     new GbPoint3D(x2a, x1ma, x3a)));
                triangles.push_back(new GbTriangle3D(new GbPoint3D(x2c, x1mb, x3c), new GbPoint3D(x2b, x1mb, x3b),
                                                     new GbPoint3D(x2a, x1mb, x3a)));
                triangles.push_back(new GbTriangle3D(new GbPoint3D(x2a, x1mb, x3a), new GbPoint3D(x2d, x1mb, x3d),
                                                     new GbPoint3D(x2c, x1mb, x3c)));
            } else if (this->isParallelToX3Axis()) {
                triangles.push_back(new GbTriangle3D(new GbPoint3D(x3a, x2a, x1ma), new GbPoint3D(x3b, x2b, x1ma),
                                                     new GbPoint3D(x3c, x2c, x1ma)));
                triangles.push_back(new GbTriangle3D(new GbPoint3D(x3c, x2c, x1ma), new GbPoint3D(x3d, x2d, x1ma),
                                                     new GbPoint3D(x3a, x2a, x1ma)));
                triangles.push_back(new GbTriangle3D(new GbPoint3D(x3c, x2c, x1mb), new GbPoint3D(x3b, x2b, x1mb),
                                                     new GbPoint3D(x3a, x2a, x1mb)));
                triangles.push_back(new GbTriangle3D(new GbPoint3D(x3a, x2a, x1mb), new GbPoint3D(x3d, x2d, x1mb),
                                                     new GbPoint3D(x3c, x2c, x1mb)));
            }
        }
    }

    return triangles;
}
/*==========================================================*/
void GbCylinder3D::addSurfaceTriangleSet(vector<UbTupleFloat3> &nodes, vector<UbTupleInt3> &triangles)
{
    float x1ma, x1mb, x2m, x3m;
    if (this->isParallelToX1Axis()) {
        x1ma = (float)this->getX1Minimum();
        x1mb = (float)this->getX1Maximum();
        x2m  = (float)this->getX2Centroid();
        x3m  = (float)this->getX3Centroid();
    } else if (this->isParallelToX2Axis()) {
        x1ma = (float)this->getX2Minimum();
        x1mb = (float)this->getX2Maximum();
        x2m  = (float)this->getX1Centroid();
        x3m  = (float)this->getX3Centroid();
    } else if (this->isParallelToX3Axis()) {
        x1ma = (float)this->getX3Minimum();
        x1mb = (float)this->getX3Maximum();
        x2m  = (float)this->getX2Centroid();
        x3m  = (float)this->getX1Centroid();
    } else
        throw UbException(UB_EXARGS, "cylinder not axis prallel");

    int segmentsCircle = 20;
    double deltaPhi    = UbMath::PI / (double)segmentsCircle;

    double phiX1a, phiX1b;
    float x1a, x2a, x3a, x1b, x2b, x3b, x1c, x2c, x3c, x1d, x2d, x3d;

    double dXCylinder    = fabs((x1mb - x1ma)) / (double)segmentsCircle;
    int segmentsCylinder = (int)(fabs(x1mb - x1ma) / dXCylinder);
    int nodenr           = 0;
    for (int segCyl = 0; segCyl < segmentsCylinder; segCyl++) {
        x1a = x1d = (float)(x1ma + segCyl * dXCylinder);
        x1b = x1c = (float)(x1a + dXCylinder);

        for (phiX1a = 2.0 * UbMath::PI; phiX1a > 0; phiX1a -= deltaPhi) {
            phiX1b = phiX1a + deltaPhi;

            x2a = (float)(x2m + mRad * std::sin(phiX1a));
            x3a = (float)(x3m + mRad * std::cos(phiX1a));
            x2b = (float)(x2m + mRad * std::sin(phiX1b));
            x3b = (float)(x3m + mRad * std::cos(phiX1b));

            if (this->isParallelToX1Axis()) {
                nodes.push_back(makeUbTuple(x1b, x2b, x3b));
                nodes.push_back(makeUbTuple(x1b, x2a, x3a));
                nodes.push_back(makeUbTuple(x1a, x2a, x3a));
                nodes.push_back(makeUbTuple(x1a, x2a, x3a));
                nodes.push_back(makeUbTuple(x1a, x2b, x3b));
                nodes.push_back(makeUbTuple(x1b, x2b, x3b));
            } else if (this->isParallelToX2Axis()) {
                nodes.push_back(makeUbTuple(x2b, x1b, x3b));
                nodes.push_back(makeUbTuple(x2a, x1b, x3a));
                nodes.push_back(makeUbTuple(x2a, x1a, x3a));
                nodes.push_back(makeUbTuple(x2a, x1a, x3a));
                nodes.push_back(makeUbTuple(x2b, x1a, x3b));
                nodes.push_back(makeUbTuple(x2b, x1b, x3b));
            } else if (this->isParallelToX3Axis()) {
                nodes.push_back(makeUbTuple(x3b, x2b, x1b));
                nodes.push_back(makeUbTuple(x3a, x2a, x1b));
                nodes.push_back(makeUbTuple(x3a, x2a, x1a));
                nodes.push_back(makeUbTuple(x3a, x2a, x1a));
                nodes.push_back(makeUbTuple(x3b, x2b, x1a));
                nodes.push_back(makeUbTuple(x3b, x2b, x1b));
            }
            triangles.push_back(makeUbTuple(nodenr, nodenr + 1, nodenr + 2));
            nodenr += 3;
            triangles.push_back(makeUbTuple(nodenr, nodenr + 1, nodenr + 2));
            nodenr += 3;
        }
    }

    int segmentsSide = (int)(mRad / dXCylinder);
    double radius0, radius1;
    for (int segCyl = 0; segCyl < segmentsSide; segCyl++) {
        radius0 = segCyl * dXCylinder;
        radius1 = radius0 + dXCylinder;
        if (segCyl == segmentsSide - 1)
            radius1 = mRad;

        for (phiX1a = 2.0 * UbMath::PI; phiX1a > 0; phiX1a -= deltaPhi) {
            phiX1b = phiX1a + deltaPhi;

            x2a = x2m + (float)(radius0 * std::sin(phiX1a));
            x3a = x3m + (float)(radius0 * std::cos(phiX1a));
            x2b = x2m + (float)(radius0 * std::sin(phiX1b));
            x3b = x3m + (float)(radius0 * std::cos(phiX1b));
            x2c = x2m + (float)(radius1 * std::sin(phiX1b));
            x3c = x3m + (float)(radius1 * std::cos(phiX1b));
            x2d = x2m + (float)(radius1 * std::sin(phiX1a));
            x3d = x3m + (float)(radius1 * std::cos(phiX1a));

            if (this->isParallelToX1Axis()) {
                nodes.push_back(makeUbTuple(x1ma, x2a, x3a));
                nodes.push_back(makeUbTuple(x1ma, x2b, x3b));
                nodes.push_back(makeUbTuple(x1ma, x2c, x3c));
                nodes.push_back(makeUbTuple(x1ma, x2c, x3c));
                nodes.push_back(makeUbTuple(x1ma, x2d, x3d));
                nodes.push_back(makeUbTuple(x1ma, x2a, x3a));
                nodes.push_back(makeUbTuple(x1mb, x2c, x3c));
                nodes.push_back(makeUbTuple(x1mb, x2b, x3b));
                nodes.push_back(makeUbTuple(x1mb, x2a, x3a));
                nodes.push_back(makeUbTuple(x1mb, x2a, x3a));
                nodes.push_back(makeUbTuple(x1mb, x2d, x3d));
                nodes.push_back(makeUbTuple(x1mb, x2c, x3c));
            } else if (this->isParallelToX2Axis()) {
                nodes.push_back(makeUbTuple(x2a, x1ma, x3a));
                nodes.push_back(makeUbTuple(x2b, x1ma, x3b));
                nodes.push_back(makeUbTuple(x2c, x1ma, x3c));
                nodes.push_back(makeUbTuple(x2c, x1ma, x3c));
                nodes.push_back(makeUbTuple(x2d, x1ma, x3d));
                nodes.push_back(makeUbTuple(x2a, x1ma, x3a));
                nodes.push_back(makeUbTuple(x2c, x1mb, x3c));
                nodes.push_back(makeUbTuple(x2b, x1mb, x3b));
                nodes.push_back(makeUbTuple(x2a, x1mb, x3a));
                nodes.push_back(makeUbTuple(x2a, x1mb, x3a));
                nodes.push_back(makeUbTuple(x2d, x1mb, x3d));
                nodes.push_back(makeUbTuple(x2c, x1mb, x3c));
            } else if (this->isParallelToX3Axis()) {
                nodes.push_back(makeUbTuple(x3a, x2a, x1ma));
                nodes.push_back(makeUbTuple(x3b, x2b, x1ma));
                nodes.push_back(makeUbTuple(x3c, x2c, x1ma));
                nodes.push_back(makeUbTuple(x3c, x2c, x1ma));
                nodes.push_back(makeUbTuple(x3d, x2d, x1ma));
                nodes.push_back(makeUbTuple(x3a, x2a, x1ma));
                nodes.push_back(makeUbTuple(x3c, x2c, x1mb));
                nodes.push_back(makeUbTuple(x3b, x2b, x1mb));
                nodes.push_back(makeUbTuple(x3a, x2a, x1mb));
                nodes.push_back(makeUbTuple(x3a, x2a, x1mb));
                nodes.push_back(makeUbTuple(x3d, x2d, x1mb));
                nodes.push_back(makeUbTuple(x3c, x2c, x1mb));
            }

            triangles.push_back(makeUbTuple(nodenr, nodenr + 1, nodenr + 2));
            nodenr += 3;
            triangles.push_back(makeUbTuple(nodenr, nodenr + 1, nodenr + 2));
            nodenr += 3;
            triangles.push_back(makeUbTuple(nodenr, nodenr + 1, nodenr + 2));
            nodenr += 3;
            triangles.push_back(makeUbTuple(nodenr, nodenr + 1, nodenr + 2));
            nodenr += 3;
        }
    }
}
/*==========================================================*/
void GbCylinder3D::addSurfaceTriangleSetSegments(vector<UbTupleFloat3> &nodes, vector<UbTupleInt3> &triangles,
                                                 int segmentsRound, int segmentsHeight)
{
    float x1ma, x1mb, x2m, x3m;
    if (this->isParallelToX1Axis()) {
        x1ma = (float)this->getX1Minimum();
        x1mb = (float)this->getX1Maximum();
        x2m  = (float)this->getX2Centroid();
        x3m  = (float)this->getX3Centroid();
    } else if (this->isParallelToX2Axis()) {
        x1ma = (float)this->getX2Minimum();
        x1mb = (float)this->getX2Maximum();
        x2m  = (float)this->getX1Centroid();
        x3m  = (float)this->getX3Centroid();
    } else if (this->isParallelToX3Axis()) {
        x1ma = (float)this->getX3Minimum();
        x1mb = (float)this->getX3Maximum();
        x2m  = (float)this->getX2Centroid();
        x3m  = (float)this->getX1Centroid();
    } else
        throw UbException(UB_EXARGS, "cylinder not axis prallel");

    int segmentsCircle = segmentsRound;
    double deltaPhi    = UbMath::PI / (double)segmentsCircle;

    double phiX1a, phiX1b;
    float x1a, x2a, x3a, x1b, x2b, x3b, x1c, x2c, x3c, x1d, x2d, x3d;

    double dXCylinder = fabs((x1mb - x1ma)) / (double)segmentsHeight; // hier evtl. segmentsheight
    int segmentsCylinder = (int)(fabs(x1mb - x1ma) / dXCylinder);
    int nodenr           = 0;
    for (int segCyl = 0; segCyl < segmentsCylinder; segCyl++) {
        x1a = x1d = (float)(x1ma + segCyl * dXCylinder);
        x1b = x1c = (float)(x1a + dXCylinder);

        // for(phiX1a=2.0*UbMath::PI; phiX1a>0.0; phiX1a-=deltaPhi)
        for (phiX1a = 0.0; phiX1a < 2.0 * UbMath::PI - 0.5 * deltaPhi; phiX1a += deltaPhi) {
            phiX1b = phiX1a + deltaPhi;

            x2a = (float)(x2m + mRad * std::sin(phiX1a));
            x3a = (float)(x3m + mRad * std::cos(phiX1a));
            x2b = (float)(x2m + mRad * std::sin(phiX1b));
            x3b = (float)(x3m + mRad * std::cos(phiX1b));

            if (this->isParallelToX1Axis()) {
                nodes.push_back(makeUbTuple(x1b, x2b, x3b));
                nodes.push_back(makeUbTuple(x1b, x2a, x3a));
                nodes.push_back(makeUbTuple(x1a, x2a, x3a));
                nodes.push_back(makeUbTuple(x1a, x2a, x3a));
                nodes.push_back(makeUbTuple(x1a, x2b, x3b));
                nodes.push_back(makeUbTuple(x1b, x2b, x3b));
            } else if (this->isParallelToX2Axis()) {
                nodes.push_back(makeUbTuple(x2b, x1b, x3b));
                nodes.push_back(makeUbTuple(x2a, x1b, x3a));
                nodes.push_back(makeUbTuple(x2a, x1a, x3a));
                nodes.push_back(makeUbTuple(x2a, x1a, x3a));
                nodes.push_back(makeUbTuple(x2b, x1a, x3b));
                nodes.push_back(makeUbTuple(x2b, x1b, x3b));
            } else if (this->isParallelToX3Axis()) {
                nodes.push_back(makeUbTuple(x3b, x2b, x1b));
                nodes.push_back(makeUbTuple(x3a, x2a, x1b));
                nodes.push_back(makeUbTuple(x3a, x2a, x1a));
                nodes.push_back(makeUbTuple(x3a, x2a, x1a));
                nodes.push_back(makeUbTuple(x3b, x2b, x1a));
                nodes.push_back(makeUbTuple(x3b, x2b, x1b));
            }
            triangles.push_back(makeUbTuple(nodenr, nodenr + 1, nodenr + 2));
            nodenr += 3;
            triangles.push_back(makeUbTuple(nodenr, nodenr + 1, nodenr + 2));
            nodenr += 3;
        }
    }

    int segmentsSide = (int)(mRad / dXCylinder);
    double radius0, radius1;
    for (int segCyl = 0; segCyl < segmentsSide; segCyl++) {
        radius0 = segCyl * dXCylinder;
        radius1 = radius0 + dXCylinder;
        if (segCyl == segmentsSide - 1)
            radius1 = mRad;

        // for(phiX1a=2.0*UbMath::PI; phiX1a>0.0; phiX1a-=deltaPhi)
        for (phiX1a = 0.0; phiX1a < 2.0 * UbMath::PI - 0.5 * deltaPhi; phiX1a += deltaPhi) {
            phiX1b = phiX1a + deltaPhi;

            x2a = x2m + (float)(radius0 * std::sin(phiX1a));
            x3a = x3m + (float)(radius0 * std::cos(phiX1a));
            x2b = x2m + (float)(radius0 * std::sin(phiX1b));
            x3b = x3m + (float)(radius0 * std::cos(phiX1b));
            x2c = x2m + (float)(radius1 * std::sin(phiX1b));
            x3c = x3m + (float)(radius1 * std::cos(phiX1b));
            x2d = x2m + (float)(radius1 * std::sin(phiX1a));
            x3d = x3m + (float)(radius1 * std::cos(phiX1a));

            if (this->isParallelToX1Axis()) {
                nodes.push_back(makeUbTuple(x1ma, x2a, x3a));
                nodes.push_back(makeUbTuple(x1ma, x2b, x3b));
                nodes.push_back(makeUbTuple(x1ma, x2c, x3c));
                nodes.push_back(makeUbTuple(x1ma, x2c, x3c));
                nodes.push_back(makeUbTuple(x1ma, x2d, x3d));
                nodes.push_back(makeUbTuple(x1ma, x2a, x3a));
                nodes.push_back(makeUbTuple(x1mb, x2c, x3c));
                nodes.push_back(makeUbTuple(x1mb, x2b, x3b));
                nodes.push_back(makeUbTuple(x1mb, x2a, x3a));
                nodes.push_back(makeUbTuple(x1mb, x2a, x3a));
                nodes.push_back(makeUbTuple(x1mb, x2d, x3d));
                nodes.push_back(makeUbTuple(x1mb, x2c, x3c));
            } else if (this->isParallelToX2Axis()) {
                nodes.push_back(makeUbTuple(x2a, x1ma, x3a));
                nodes.push_back(makeUbTuple(x2b, x1ma, x3b));
                nodes.push_back(makeUbTuple(x2c, x1ma, x3c));
                nodes.push_back(makeUbTuple(x2c, x1ma, x3c));
                nodes.push_back(makeUbTuple(x2d, x1ma, x3d));
                nodes.push_back(makeUbTuple(x2a, x1ma, x3a));
                nodes.push_back(makeUbTuple(x2c, x1mb, x3c));
                nodes.push_back(makeUbTuple(x2b, x1mb, x3b));
                nodes.push_back(makeUbTuple(x2a, x1mb, x3a));
                nodes.push_back(makeUbTuple(x2a, x1mb, x3a));
                nodes.push_back(makeUbTuple(x2d, x1mb, x3d));
                nodes.push_back(makeUbTuple(x2c, x1mb, x3c));
            } else if (this->isParallelToX3Axis()) {
                nodes.push_back(makeUbTuple(x3a, x2a, x1ma));
                nodes.push_back(makeUbTuple(x3b, x2b, x1ma));
                nodes.push_back(makeUbTuple(x3c, x2c, x1ma));
                nodes.push_back(makeUbTuple(x3c, x2c, x1ma));
                nodes.push_back(makeUbTuple(x3d, x2d, x1ma));
                nodes.push_back(makeUbTuple(x3a, x2a, x1ma));
                nodes.push_back(makeUbTuple(x3c, x2c, x1mb));
                nodes.push_back(makeUbTuple(x3b, x2b, x1mb));
                nodes.push_back(makeUbTuple(x3a, x2a, x1mb));
                nodes.push_back(makeUbTuple(x3a, x2a, x1mb));
                nodes.push_back(makeUbTuple(x3d, x2d, x1mb));
                nodes.push_back(makeUbTuple(x3c, x2c, x1mb));
            }

            triangles.push_back(makeUbTuple(nodenr, nodenr + 1, nodenr + 2));
            nodenr += 3;
            triangles.push_back(makeUbTuple(nodenr, nodenr + 1, nodenr + 2));
            nodenr += 3;
            triangles.push_back(makeUbTuple(nodenr, nodenr + 1, nodenr + 2));
            nodenr += 3;
            triangles.push_back(makeUbTuple(nodenr, nodenr + 1, nodenr + 2));
            nodenr += 3;
        }
    }
}

/*==========================================================*/
void GbCylinder3D::objectChanged(UbObservable *changedObject)
{
    GbLine3D *line = dynamic_cast<GbLine3D *>(changedObject);
    if (!line || this->mLine != line)
        return;

    this->notifyObserversObjectChanged();
}
/*==========================================================*/
void GbCylinder3D::objectWillBeDeleted(UbObservable *objectForDeletion)
{
    if (this->mLine) {
        UbObservable *observedObj = dynamic_cast<UbObservable *>(this->mLine);
        if (objectForDeletion == observedObj) {
            this->mLine = NULL;
        }
    }
}
/*=======================================================*/
void GbCylinder3D::scale(const double &sx1, const double &sx2, const double &sx3)
{
    if (this->isParallelToX1Axis()) {
        if (!UbMath::equal(sx2, sx3))
            throw UbException(UB_EXARGS, "|| to x1 -> different scaling sx2 and sx3 not possible");
        this->mRad *= sx2;
    } else if (this->isParallelToX2Axis()) {
        if (!UbMath::equal(sx1, sx3))
            throw UbException(UB_EXARGS, "|| to x2 -> different scaling sx1 and sx3 not possible");
        this->mRad *= sx1;
    } else if (this->isParallelToX3Axis()) {
        if (!UbMath::equal(sx1, sx2))
            throw UbException(UB_EXARGS, "|| to x3 -> different scaling sx1 and sx2 not possible");
        this->mRad *= sx1;
    } else
        throw UbException(UB_EXARGS, "unknown direction");

    this->mLine->scale(sx1, sx2, sx3);
    // notify observer wird automatisch aufgerufen
}

/*==========================================================*/
double GbCylinder3D::getIntersectionRaytraceFactor(const double &x1, const double &x2, const double &x3,
                                                   const double &rx1, const double &rx2, const double &rx3)
{
    /*
    Distance D of the intersection between a Ray((ox1,ox2,ox3),(dx1,dx2,dx3)) and a Plane P: ax+by+cz+d=0
    dc = a*dx1 + b*dx2 + c*dx3
    dw = a*ox1 + b*ox2 + c*ox3 + d
    D =   - dw / dc
    */
    double px1, px2, px3;
    double d = Ub::inf; // Distance to Min or Max Plane of the Zylinder
                        // final distance should be less that d

    if (this->isParallelToX1Axis()) {
        if (UbMath::equal(x1, minX1) && UbMath::negative(rx1))
            return -1.0;
        else if (UbMath::equal(x1, maxX1) && UbMath::positive(rx1))
            return -1.0;

        // falls die Linie nicht parallel zu den Seitenflaechen ist
        if (x1 < minX1 || x1 > maxX1) // nur fuer punkte links und rechts des cylinders
        {
            px1 = (x1 < minX1 ? minX1 : maxX1);
            // falls die Linie nicht parallel zu den Seitenflaechen ist
            if (!UbMath::zero(rx1)) {
                // Plane a= 0, b= 1, c=0 d= -1*px2
                d   = -1.0 * (x1 - px1) / rx1;
                px2 = x2 + d * rx2;
                px3 = x3 + d * rx3;

                if (UbMath::greater(mLine->getDistance(px1, px2, px3), mRad)) {
                    if (x1 < minX1 && rx1 > 0.0)
                        d = Ub::inf; // punkt liegt "links" vom cylinder und strahl hat evtl weiteren SP auf oberflaeche
                    else if (x1 > maxX1 && rx1 < 0.0)
                        d = Ub::inf;
                    else
                        return -1.0;
                } else
                    return d;
            } else
                return -1.0;
        } else {
            if (UbMath::negative(rx1))
                d = -1.0 * (x1 - minX1) / rx1;
            else if (UbMath::positive(rx1))
                d = -1.0 * (x1 - maxX1) / rx1;
        }
    } else if (this->isParallelToX2Axis()) {
        if (UbMath::equal(x2, minX2) && UbMath::negative(rx2))
            return -1;
        else if (UbMath::equal(x2, maxX2) && UbMath::positive(rx2))
            return -1;

        if (minX2 > x2 || x2 > maxX2) {
            px2 = (x2 < minX2 ? minX2 : maxX2);
            // falls die Linie nicht parallel zu den Seitenflaechen ist
            if (!UbMath::zero(rx2)) {
                // Plane a= 0, b= 1, c=0 d= -1*px2
                d   = -1 * (x2 - px2) / rx2;
                px1 = x1 + d * rx1;
                px3 = x3 + d * rx3;

                if (UbMath::greater(mLine->getDistance(px1, px2, px3), mRad)) {
                    if (x2 < minX2 && rx2 > 0.0)
                        d = Ub::inf; // punkt liegt "links oberhalb" vom cylinder und strahl mit pos x1 hat evtl
                                     // weiteren SP auf oberflaeche
                    else if (x2 > maxX2 && rx2 < 0.0)
                        d = Ub::inf;
                    else
                        return -1.0;
                } else
                    return d;
            } else
                return -1.0;
        } else {
            if (UbMath::negative(rx2))
                d = -1.0 * (x2 - minX2) / rx2;
            else if (UbMath::positive(rx2))
                d = -1.0 * (x2 - maxX2) / rx2;
        }
    } else if (this->isParallelToX3Axis()) {
        if (UbMath::equal(x3, minX3) && UbMath::negative(rx3))
            return -1.0;
        else if (UbMath::equal(x3, maxX3) && UbMath::positive(rx3))
            return -1.0;

        if (minX3 > x3 || x3 > maxX3) {
            px3 = (x3 < minX3 ? minX3 : maxX3);
            // falls die Linie nicht parallel zu den Seitenflaechen ist
            if (!UbMath::zero(rx3)) {
                // Plane a= 0, b= 0, c=1 d= -1*px3
                d   = -1.0 * (x3 - px3) / rx3;
                px2 = x2 + d * rx2;
                px1 = x1 + d * rx1;
                if (UbMath::greater(mLine->getDistance(px1, px2, px3), mRad)) {
                    if (x3 < minX3 && rx3 > 0.0)
                        d = Ub::inf;
                    else if (x3 > maxX3 && rx3 < 0.0)
                        d = Ub::inf;
                    else
                        return -1.0;
                } else
                    return d;
            } else
                return -1.0;
        } else {
            if (UbMath::negative(rx3))
                d = -1.0 * (x3 - minX3) / rx3;
            else if (UbMath::positive(rx3))
                d = -1.0 * (x3 - maxX3) / rx3;
        }
    } else
        throw UbException(UB_EXARGS, "funzt nur bei achsen parallelem cylinder");
    //////////////////////////////////////////////////////////////////////////
    // Q berechnen fuer Infinity Zylinder
    double axisX1 = mLine->getPoint2()->x1 - mLine->getPoint1()->x1; /* Axis of the cylinder   */
    double axisX2 = mLine->getPoint2()->x2 - mLine->getPoint1()->x2; /* mit p1 als base of cylinder */
    double axisX3 = mLine->getPoint2()->x3 - mLine->getPoint1()->x3;

    // double dirlen = mLine->getLength();
    // double abs, t, s;

    double RCx1 = x1 - mLine->getPoint1()->x1;
    double RCx2 = x2 - mLine->getPoint1()->x2;
    double RCx3 = x3 - mLine->getPoint1()->x3;

    // n = ray x axis
    double nx1     = rx2 * axisX3 - rx3 * axisX2;
    double nx2     = rx3 * axisX1 - rx1 * axisX3;
    double nx3     = rx1 * axisX2 - rx2 * axisX1;
    double nLength = nx1 * nx1 + nx2 * nx2 + nx3 * nx3;

    double abs;
    if (UbMath::zero(nLength)) { /* ray parallel to cyl  */
        // abs = RC dot axis
        double tmpabs = RCx1 * axisX1 + RCx2 * axisX2 + RCx3 * axisX3;
        double dx1    = RCx1 - tmpabs * axisX1;
        double dx2    = RCx2 - tmpabs * axisX2;
        double dx3    = RCx3 - tmpabs * axisX3;
        if (UbMath::greater(dx1 * dx1 + dx2 * dx2 + dx3 * dx3, mRad * mRad))
            return -1.0;
    }

    // normalize "n"
    nLength           = std::sqrt(nLength);
    double invnLength = 1.0 / nLength;
    nx1 *= invnLength;
    nx2 *= invnLength;
    nx3 *= invnLength;

    // shortest distance  = fabs( RC dot n )
    abs = std::fabs(RCx1 * nx1 + RCx2 * nx2 + RCx3 * nx3);

    if (UbMath::lessEqual(abs, mRad)) { /* if ray hits cylinder */
        // Ox1 = RC x axis
        double Ox1 = RCx2 * axisX3 - RCx3 * axisX2;
        double Ox2 = RCx3 * axisX1 - RCx1 * axisX3;
        double Ox3 = RCx1 * axisX2 - RCx2 * axisX1;
        // t = - O dot n / nLength;
        double t = -(Ox1 * nx1 + Ox2 * nx2 + Ox3 * nx3) / nLength;

        // O = n x axis;
        Ox1 = nx2 * axisX3 - nx3 * axisX2;
        Ox2 = nx3 * axisX1 - nx1 * axisX3;
        Ox3 = nx1 * axisX2 - nx2 * axisX1;

        // normalize O
        invnLength = 1.0 / std::sqrt(Ox1 * Ox1 + Ox2 * Ox2 + Ox3 * Ox3);
        Ox1 *= invnLength;
        Ox2 *= invnLength;
        Ox3 *= invnLength;

        double s = std::fabs(sqrt(mRad * mRad - abs * abs) / (rx1 * Ox1 + rx2 * Ox2 + rx3 * Ox3));

        // Wert a) t-s: entering distance
        //     b) t+s: exiting  distance
        //
        // -> we only consider factors in ray-dir -> means positive values!
        //    (s is always positive)

        if (t > s) {
            return UbMath::min(t - s, d);
        } else if ((t + s) > 0) {
            return UbMath::min(t + s, d);
        }
    }

    return -1.0;
}
/*==========================================================*/

//! \}
