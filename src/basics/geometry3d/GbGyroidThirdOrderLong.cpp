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
//! \file GbGyroidThirdOrderLong.cpp
//! \ingroup geometry3d
//! \author Hussein Alihussein
//=======================================================================================

#include <GbGyroidThirdOrderLong.h>

#ifdef BUILD_USE_BOOST

#include <basics/utilities/UbMath.h>

#include <geometry3d/GbSystem3D.h>
#include <geometry3d/GbTriangle3D.h>

#include <boost/math/tools/roots.hpp>


using namespace std;
using boost::math::tools::bisect;

// Konstruktor
GbGyroidThirdOrderLong::GbGyroidThirdOrderLong() //: GbObject3D()
{

}
/*=======================================================*/
// Konstruktor
GbGyroidThirdOrderLong::GbGyroidThirdOrderLong(const double& x1a, const double& x2a, const double& x3a, const double& x1b, const double& x2b, const double& x3b, const double& edgeLength, const double& dx, const double& thickness) :GbObject3D()
{
    this->p1 = new GbPoint3D(x1a, x2a, x3a);
    this->p2 = new GbPoint3D(x1b, x2b, x3b);
    this->p1->addObserver(this);
    this->p2->addObserver(this);

    this->p3 = new GbPoint3D(x1a, x2a, x3a);
    this->p4 = new GbPoint3D(x1b, x2b, x3b);
    this->p3->addObserver(this);
    this->p4->addObserver(this);

    this->edgeLength = edgeLength;
    this->dx = dx;
    this->thickness = thickness;
}
/*=======================================================*/
GbGyroidThirdOrderLong::GbGyroidThirdOrderLong(GbGyroidThirdOrderLong * imp)
{
}
/*=======================================================*/
// Destruktor
GbGyroidThirdOrderLong::~GbGyroidThirdOrderLong()
{
    if (this->p1)
        this->p1->removeObserver(this);
    if (this->p2)
        this->p2->removeObserver(this);
    if (this->p3)
        this->p3->removeObserver(this);
    if (this->p4)
        this->p4->removeObserver(this);
}
/*=======================================================*/
struct TerminationCondition {
    bool operator() (double min, double max) {
        return abs(min - max) <= 10e-10;
    }
};
/*============================================M-===========*/
struct FunctionToApproximate {
    double x, y, z;
    double dir1, dir2, dir3, L;
    double operator() (double q) {
        return sin(2.*M_PI / L*(x + q*dir1))*cos(2.*M_PI / L*(y + q*dir2)) + sin(2.*M_PI / L*(y + q*dir2))*cos(2.*M_PI / L*(z + q*dir3)) + sin(2.*M_PI / L*(z + q*dir3))*cos(2.*M_PI / L*(x + q*dir1));
    }
};
/*=======================================================*/
struct FunctionGyroidThirdOrder {
    double x, y, z;
    double dir1, dir2, dir3, L;
    double h;
    
    double t17, t3, t2, t18, t20, t8, t13, t5, t9, t6, t11, t14;
    double f300, f210, f201, f120, f102, f030, f021, f012, f003, f200, f110, f101, f020, f011, f002, f100, f010, f001, f000;

    double repeatedTerm, repeatedTermRoot;
    double T1, T2, T3, T4, T5, T6, T7, T8, T9, Gyroidh;

    double operator() (double q) {
    //sins and cosines combinations 
     t2  = sin((2. * M_PI*(x+q*dir1)) / L)*sin((2. * M_PI*(y+q*dir2)) / L);
     t3  = sin((2. * M_PI*(x+q*dir1)) / L)*sin((2. * M_PI*(z+q*dir3)) / L);
     t5  = cos((2. * M_PI*(y+q*dir2)) / L)*sin((2. * M_PI*(x+q*dir1)) / L);
     t6  = cos((2. * M_PI*(z+q*dir3)) / L)*sin((2. * M_PI*(x+q*dir1)) / L);
     t8  = sin((2. * M_PI*(y+q*dir2)) / L)*sin((2. * M_PI*(z+q*dir3)) / L);
     t9  = cos((2. * M_PI*(x+q*dir1)) / L)*sin((2. * M_PI*(y+q*dir2)) / L);
     t11 = cos((2. * M_PI*(z+q*dir3)) / L)*sin((2. * M_PI*(y+q*dir2)) / L);
     t13 = cos((2. * M_PI*(x+q*dir1)) / L)*sin((2. * M_PI*(z+q*dir3)) / L);
     t14 = cos((2. * M_PI*(y+q*dir2)) / L)*sin((2. * M_PI*(z+q*dir3)) / L);
     t17 = cos((2. * M_PI*(x+q*dir1)) / L)*cos((2. * M_PI*(y+q*dir2)) / L);
     t18 = cos((2. * M_PI*(x+q*dir1)) / L)*cos((2. * M_PI*(z+q*dir3)) / L);
     t20 = cos((2. * M_PI*(y+q*dir2)) / L)*cos((2. * M_PI*(z+q*dir3)) / L);

    //Gyroid third order derivatives
     f300 = (8. * pow(M_PI, 3.)*(-t17 + t3)) / pow(L, 3.);
     f210 = (8. * pow(M_PI, 3.)*t2) / pow(L, 3.);
     f201 = (-8. * pow(M_PI, 3.)*t18) / pow(L, 3.);
     f120 = (-8. * pow(M_PI, 3.)*t17) / pow(L, 3.);
     f102 = (8. * pow(M_PI, 3.)*t3) / pow(L, 3.);
     f030 = (8. * pow(M_PI, 3.)*(t2 - t20)) / pow(L, 3.);
     f021 = (8. * pow(M_PI, 3.)*t8) / pow(L, 3.);
     f012 = (-8. * pow(M_PI, 3.)*t20) / pow(L, 3.);
     f003 = (8. * pow(M_PI, 3.)*(-t18 + t8)) / pow(L, 3.);

    //Gyroid second order derivatives
     f200 = (-4. * pow(M_PI, 2.)*(t13 + t5)) / pow(L, 2.);
     f110 = (-4. * pow(M_PI, 2.)*t9) / pow(L, 2.);
     f101 = (-4. * pow(M_PI, 2.)*t6) / pow(L, 2.);
     f020 = (-4. * pow(M_PI, 2.)*(t11 + t5)) / pow(L, 2.);
     f011 = (-4. * pow(M_PI, 2.)*t14) / pow(L, 2.);
     f002 = (-4. * pow(M_PI, 2.)*(t11 + t13)) / pow(L, 2.);

    //Gyroid first order derivatives
     f100 = (2. * M_PI*(t17 - t3)) / L;
     f010 = (2. * M_PI*(-t2 + t20)) / L;
     f001 = (2. * M_PI*(t18 - t8)) / L;

    //Gyroid 
     f000 = t11 + t13 + t5;

     repeatedTerm = f100*f100 + f010*f010 + f001*f001;
     repeatedTermRoot = sqrt(repeatedTerm);

     T1 = f001*f002 + f010*f011 + f100*f101;
     T2 = f001*f011 + f010*f020 + f100*f110;
     T3 = f001*f101 + f010*f110 + f100*f200;
     T4 = f002*f011 + f001*f012 + f011*f020 + f010*f021 + f101*f110;
     T5 = f002*f101 + f001*f102 + f011*f110 + f101*f200 + f100*f201;
     T6 = f011*f101 + f020*f110 + f010*f120 + f110*f200 + f100*f210;
     T7 = f001*f002*h + f010*f011*h + f100*f101*h;
     T8 = f001*f011*h + f010*f020*h + f100*f110*h;
     T9 = f001*f101*h + f010*f110*h + f100*f200*h;


     Gyroidh = 2 * h*sqrt(pow(f001 - (T1*h) / (2.*repeatedTermRoot), 2) + pow(f010 - (T2*h) / (2.*repeatedTermRoot), 2) + pow(f100 - (T3*h) / (2.*repeatedTermRoot), 2))
        - (3 * h*sqrt(pow(f001 - (T1*h) / (3.*repeatedTermRoot), 2) + pow(f010 - (T2*h) / (3.*repeatedTermRoot), 2) + pow(f100 - (T3*h) / (3.*repeatedTermRoot), 2))) / 2.
        - (3 * h*sqrt(pow(f001 - (T1*h) / (3.*repeatedTermRoot) + (h*((T7 - 3 * f001*repeatedTermRoot)*
        (4 * pow(T1, 2)*h - 4 * (pow(f002, 2) + f001*f003 + pow(f011, 2) + f010*f012 + pow(f101, 2) + f100*f102)*h*repeatedTerm + 12 * f002*pow(repeatedTerm, 1.5)) +
            (T8 - 3 * f010*repeatedTermRoot)*(4 * T1*T2*h - 4 * (T4)*h*repeatedTerm + 12 * f011*pow(repeatedTerm, 1.5)) +
            (T9 - 3 * f100*repeatedTermRoot)*(4 * T1*T3*h - 4 * (T5)*h*repeatedTerm + 12 * f101*pow(repeatedTerm, 1.5)))) /
            (108.*sqrt(pow(f001 - (T1*h) / (3.*repeatedTermRoot), 2) + pow(f010 - (T2*h) / (3.*repeatedTermRoot), 2) + pow(f100 - (T3*h) / (3.*repeatedTermRoot), 2))*
                pow(repeatedTerm, 2)), 2) + pow(f010 - (T2*h) / (3.*repeatedTermRoot) +
                (h*((T7 - 3 * f001*repeatedTermRoot)*(4 * T1*T2*h - 4 * (T4)*h*repeatedTerm + 12 * f011*pow(repeatedTerm, 1.5)) +
                    (T8 - 3 * f010*repeatedTermRoot)*(4 * pow(T2, 2)*h - 4 * (pow(f011, 2) + pow(f020, 2) + f001*f021 + f010*f030 + pow(f110, 2) + f100*f120)*h*repeatedTerm + 12 * f020*pow(repeatedTerm, 1.5)) +
                    (T9 - 3 * f100*repeatedTermRoot)*(4 * T2*T3*h - 4 * (T6)*h*repeatedTerm + 12 * f110*pow(repeatedTerm, 1.5)))) /
                    (108.*sqrt(pow(f001 - (T1*h) / (3.*repeatedTermRoot), 2) + pow(f010 - (T2*h) / (3.*repeatedTermRoot), 2) + pow(f100 - (T3*h) / (3.*repeatedTermRoot), 2))*
                        pow(repeatedTerm, 2)), 2) + pow(f100 - (T3*h) / (3.*repeatedTermRoot) +
                        (h*((T7 - 3 * f001*repeatedTermRoot)*(4 * T1*T3*h - 4 * (T5)*h*repeatedTerm + 12 * f101*pow(repeatedTerm, 1.5)) +
                            (T8 - 3 * f010*repeatedTermRoot)*(4 * T2*T3*h - 4 * (T6)*h*repeatedTerm + 12 * f110*pow(repeatedTerm, 1.5)) +
                            (T9 - 3 * f100*repeatedTermRoot)*(4 * pow(T3, 2)*h - 4 * (pow(f101, 2) + pow(f110, 2) + pow(f200, 2) + f001*f201 + f010*f210 + f100*f300)*h*repeatedTerm + 12 * f200*pow(repeatedTerm, 1.5)))) /
                            (108.*sqrt(pow(f001 - (T1*h) / (3.*repeatedTermRoot), 2) + pow(f010 - (T2*h) / (3.*repeatedTermRoot), 2) + pow(f100 - (T3*h) / (3.*repeatedTermRoot), 2))*
                                pow(repeatedTerm, 2)), 2))) / 2. + f000;
    
        return Gyroidh;
    }
};
/*==========================================================*/
bool GbGyroidThirdOrderLong::isPointInGbObject3D(const double& x1, const double& x2, const double& x3)
{
    //double f = sin(2.*M_PI*x1/edgeLength)*cos(2.*M_PI*x2 / edgeLength) + sin(2.*M_PI*x2 / edgeLength)*cos(2.*M_PI*x3 / edgeLength) + sin(2.*M_PI*x3 / edgeLength)*cos(2.*M_PI*x1 / edgeLength);
    //evaluateImplicitFunction(x1,x2,x3, 0., 0., 0.)
    double f1 = evaluateImplicitFunction(x1, x2, x3, 1.);
    double f2 = evaluateImplicitFunction(x1, x2, x3, -1.);
    return UbMath::lessEqual(f1,0.) && UbMath::greaterEqual(f2,0.);
}

/*==========================================================*/
double GbGyroidThirdOrderLong::getIntersectionRaytraceFactor(const double& x1, const double& x2, const double& x3, const double& rx1, const double& rx2, const double& rx3)
{
    double from = 0;  // The solution must lie in the interval [from, to], additionally f(from) <= 0 && f(to) >= 0
    double to = dx*sqrt(rx1*rx1+ rx2*rx2+ rx3*rx3);
    FunctionGyroidThirdOrder f;
    f.x =x1 ;
    f.y =x2 ;
    f.z =x3 ;
    f.dir1 = rx1;
    f.dir2 = rx2;
    f.dir3 = rx3;
    f.L = edgeLength;
    f.h = thickness;
    if (f(from)*f(to)<0)
        {
        std::pair<double, double> result = bisect(f, from, to, TerminationCondition());
        double root = (result.first + result.second) / 2;
        return root;
        }
    f.h = -thickness;
    if (f(from)*f(to) < 0)
    {
        std::pair<double, double> result = bisect(f, from, to, TerminationCondition());
        double root = (result.first + result.second) / 2;
        return root;
    }
    else
    {
        return 999;
    }
    
}
/*=======================================================*/
double GbGyroidThirdOrderLong::evaluateImplicitFunction(const double& x1, const double& x2, const double& x3, const double& position)
{
    double to = 0.;
    FunctionGyroidThirdOrder f;
    f.x = x1;
    f.y = x2;
    f.z = x3;
    f.dir1 = 0.;
    f.dir2 = 0.;
    f.dir3 = 0.;
    f.L = edgeLength;
    f.h = position*thickness;
    return f(to);
}
/*=======================================================*/
double GbGyroidThirdOrderLong::getX1Centroid()
{
    return (0.5*(p1->x1 + p2->x1));
}
/*=======================================================*/
double GbGyroidThirdOrderLong::getX1Minimum()
{
    return (this->p1->x1 < this->p2->x1 ? this->p1->x1 : this->p2->x1);
}
/*=======================================================*/
double GbGyroidThirdOrderLong::getX1Maximum()
{
    return (this->p1->x1 > this->p2->x1 ? this->p1->x1 : this->p2->x1);
}
/*=======================================================*/
double GbGyroidThirdOrderLong::getX2Centroid()
{
    return (0.5*(p1->x2 + p2->x2));
}
/*=======================================================*/
double GbGyroidThirdOrderLong::getX2Minimum()
{
    return (this->p1->x2 < this->p2->x2 ? this->p1->x2 : this->p2->x2);
}
/*=======================================================*/
double GbGyroidThirdOrderLong::getX2Maximum()
{
    return (this->p1->x2 > this->p2->x2 ? this->p1->x2 : this->p2->x2);
}
/*=======================================================*/
double GbGyroidThirdOrderLong::getX3Centroid()
{
    return (0.5*(p1->x3 + p2->x3));
}
/*=======================================================*/
double GbGyroidThirdOrderLong::getX3Minimum()
{
    return (this->p1->x3 < this->p2->x3 ? this->p1->x3 : this->p2->x3);
}
/*=======================================================*/
double GbGyroidThirdOrderLong::getX3Maximum()
{
    return (this->p1->x3 > this->p2->x3 ? this->p1->x3 : this->p2->x3);
}
/*=======================================================*/
bool GbGyroidThirdOrderLong::isCellInsideGbObject3D(const double& x1a, const double& x2a, const double& x3a, const double& x1b, const double& x2b, const double& x3b)
{
    if (this->isPointInGbObject3D(x1a, x2a, x3a)
        && this->isPointInGbObject3D(x1b, x2a, x3a)
        && this->isPointInGbObject3D(x1b, x2b, x3a)
        && this->isPointInGbObject3D(x1a, x2b, x3a)
        && this->isPointInGbObject3D(x1a, x2a, x3b)
        && this->isPointInGbObject3D(x1b, x2a, x3b)
        && this->isPointInGbObject3D(x1b, x2b, x3b)
        && this->isPointInGbObject3D(x1a, x2b, x3b))
    {
        return true;
    }
    return false;
}
/*=======================================================*/
bool GbGyroidThirdOrderLong::isCellInsideOrCuttingGbObject3D(const double& x1a, const double& x2a, const double& x3a, const double& x1b, const double& x2b, const double& x3b)
{
    if ((this->isPointInGbObject3D(x1a, x2a, x3a) == false)
        && (this->isPointInGbObject3D(x1b, x2a, x3a) == false)
        && (this->isPointInGbObject3D(x1b, x2b, x3a) == false)
        && (this->isPointInGbObject3D(x1a, x2b, x3a) == false)
        && (this->isPointInGbObject3D(x1a, x2a, x3b) == false)
        && (this->isPointInGbObject3D(x1b, x2a, x3b) == false)
        && (this->isPointInGbObject3D(x1b, x2b, x3b) == false)
        && (this->isPointInGbObject3D(x1a, x2b, x3b) == false))
    {
        return false;
    }
    return true;
}
/*=======================================================*/
bool GbGyroidThirdOrderLong::isCellCuttingGbObject3D(const double& x1a, const double& x2a, const double& x3a, const double& x1b, const double& x2b, const double& x3b)
{
    if (!this->isCellInsideGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b)
        && this->isCellInsideOrCuttingGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b))
    {
        return true;
    }
    return false;
}
/*==========================================================*/
void GbGyroidThirdOrderLong::objectChanged(UbObservable *changedObject)
{
    GbPoint3D *point = dynamic_cast<GbPoint3D *>(changedObject);
    if (!point || (this->p1 != point && this->p2 != point && this->p3 != point && this->p4 != point))
        return;

    this->notifyObserversObjectChanged();
}
/*==========================================================*/
void GbGyroidThirdOrderLong::objectWillBeDeleted(UbObservable *objectForDeletion)
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
    if (this->p3) {
        UbObservable *observedObj = dynamic_cast<UbObservable *>(this->p3);
        if (objectForDeletion == observedObj) {
            this->p3 = NULL;
        }
    }
    if (this->p4) {
        UbObservable *observedObj = dynamic_cast<UbObservable *>(this->p4);
        if (objectForDeletion == observedObj) {
            this->p4 = NULL;
        }
    }
}

#endif
