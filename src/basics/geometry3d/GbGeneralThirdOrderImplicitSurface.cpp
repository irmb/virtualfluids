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
//! \file GbGeneralThirdOrderImplicitSurface.cpp
//! \ingroup geometry3d
//! \author Hussein Alihussein
//=======================================================================================

#include <GbGeneralThirdOrderImplicitSurface.h>

#ifdef BUILD_USE_BOOST

#include <basics/utilities/UbMath.h>

#include <geometry3d/GbSystem3D.h>
#include <geometry3d/GbTriangle3D.h>

#include <boost/math/tools/roots.hpp>


using namespace std;
using boost::math::tools::bisect;

/*=======================================================*/
// ObObjectCreator* GbGeneralThirdOrderImplicitSurface::getCreator()
// {
// 	 GbObject3DCreator instance;
// 	return &instance;
// }
/*=======================================================*/
// Konstruktor
GbGeneralThirdOrderImplicitSurface::GbGeneralThirdOrderImplicitSurface() //: GbObject3D()
{

}
/*=======================================================*/
// Konstruktor
GbGeneralThirdOrderImplicitSurface::GbGeneralThirdOrderImplicitSurface(const double& x1a, const double& x2a, const double& x3a, const double& x1b, const double& x2b, const double& x3b, const double& edgeLength, const double& dx, const double& thickness) :GbObject3D()
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

GbGeneralThirdOrderImplicitSurface::GbGeneralThirdOrderImplicitSurface(GbGeneralThirdOrderImplicitSurface * imp)
{
}
/*=======================================================*/
// Destruktor
GbGeneralThirdOrderImplicitSurface::~GbGeneralThirdOrderImplicitSurface()
{
   
}
/*=======================================================*/
struct TerminationCondition {
	bool operator() (double min, double max) {
		return abs(min - max) <= 10e-10;
	}
};
/*=======================================================*/
struct FunctionGyroidThirdOrder {
	double x, y, z;
	double dir1, dir2, dir3, L;
	double h;
	

	double operator() (double q) {
	double f111,f300, f210, f201, f120, f102, f030, f021, f012, f003, f200, f110, f101, f020, f011, f002, f100, f010, f001, f;

	double repeatedTerm1,repeatedTerm2,repeatedTerm3,repeatedTerm4,repeatedTerm5,repeatedTerm6,repeatedTerm7,repeatedTerm8,repeatedTerm9,repeatedTerm10,repeatedTerm11,repeatedTerm12,repeatedTerm13,repeatedTerm14,repeatedTerm15,repeatedTerm16,repeatedTerm17,distancedImplicit;

	double t2  = sin((2. * M_PI*(x+q*dir1)) / L)*sin((2. * M_PI*(y+q*dir2)) / L);
	double t3  = sin((2. * M_PI*(x+q*dir1)) / L)*sin((2. * M_PI*(z+q*dir3)) / L);
	double t5  = cos((2. * M_PI*(y+q*dir2)) / L)*sin((2. * M_PI*(x+q*dir1)) / L);
	double t6  = cos((2. * M_PI*(z+q*dir3)) / L)*sin((2. * M_PI*(x+q*dir1)) / L);
	double t8  = sin((2. * M_PI*(y+q*dir2)) / L)*sin((2. * M_PI*(z+q*dir3)) / L);
	double t9  = cos((2. * M_PI*(x+q*dir1)) / L)*sin((2. * M_PI*(y+q*dir2)) / L);
	double t11 = cos((2. * M_PI*(z+q*dir3)) / L)*sin((2. * M_PI*(y+q*dir2)) / L);
	double t13 = cos((2. * M_PI*(x+q*dir1)) / L)*sin((2. * M_PI*(z+q*dir3)) / L);
	double t14 = cos((2. * M_PI*(y+q*dir2)) / L)*sin((2. * M_PI*(z+q*dir3)) / L);
	double t17 = cos((2. * M_PI*(x+q*dir1)) / L)*cos((2. * M_PI*(y+q*dir2)) / L);
	double t18 = cos((2. * M_PI*(x+q*dir1)) / L)*cos((2. * M_PI*(z+q*dir3)) / L);
	double t20 = cos((2. * M_PI*(y+q*dir2)) / L)*cos((2. * M_PI*(z+q*dir3)) / L);
	
	
	//Implicit function third order derivatives
	 double L_cubed = pow(L, 3.);
	 double PI_cubed =  pow(M_PI, 3.);
	 f300 = (8. * PI_cubed*(-t17 + t3)) / L_cubed;
	 f210 = (8. * PI_cubed*t2) / L_cubed;
	 f201 = (-8. * PI_cubed*t18) / L_cubed;
	 f120 = (-8. * PI_cubed*t17) / L_cubed;
	 f102 = (8. * PI_cubed*t3) / L_cubed;
	 f030 = (8. * PI_cubed*(t2 - t20)) / L_cubed;
	 f021 = (8. * PI_cubed*t8) / L_cubed;
	 f012 = (-8. * PI_cubed*t20) / L_cubed;
	 f003 = (8. * PI_cubed*(-t18 + t8)) / L_cubed;
	 f111 = 0.;
	//Implicit function second order derivatives		
	 double L_squared = pow(L, 2.);
	 double PI_squared =  pow(M_PI, 2.);
	 f200 = (-4. * PI_squared*(t13 + t5)) / L_squared;
	 f110 = (-4. * PI_squared*t9) / L_squared;
	 f101 = (-4. * PI_squared*t6) / L_squared;
	 f020 = (-4. * PI_squared*(t11 + t5)) / L_squared;
	 f011 = (-4. * PI_squared*t14) / L_squared;
	 f002 = (-4. * PI_squared*(t11 + t13)) / L_squared;

	//Implicit function first order derivatives
	 f100 = (2. * M_PI*(t17 - t3)) / L;
	 f010 = (2. * M_PI*(-t2 + t20)) / L;
	 f001 = (2. * M_PI*(t18 - t8)) / L;

	//Implicit function 
	 f = t11 + t13 + t5;

repeatedTerm1 =sqrt(pow(f001,2) + pow(f010,2) + pow(f100,2));
repeatedTerm2 =pow(pow(f001,2) + pow(f010,2) + pow(f100,2),1.5);

repeatedTerm3 = 2.*f001*f002 + 2.*f010*f011 + 2.*f100*f101;
repeatedTerm4 = 2.*f001*f011 + 2.*f010*f020 + 2.*f100*f110;
repeatedTerm5 = 2.*f001*f101 + 2.*f010*f110 + 2.*f100*f200;

repeatedTerm6 = (0.16666666666666666666666666666667*(2.*f002*f011 + 2.*f001*f012 + 2.*f011*f020 + 2.*f010*f021 + 2.*f101*f110 + 2.*f100*f111)*h)/repeatedTerm1;
repeatedTerm7 = (0.16666666666666666666666666666667*(2.*f002*f101 + 2.*f001*f102 + 2.*f011*f110 + 2.*f010*f111 + 2.*f101*f200 + 2.*f100*f201)*h)/repeatedTerm1;
repeatedTerm8 = (0.16666666666666666666666666666667*(2.*f011*f101 + 2.*f020*f110 + 2.*f001*f111 + 2.*f010*f120 + 2.*f110*f200 + 2.*f100*f210)*h)/repeatedTerm1;

repeatedTerm9  = (0.16666666666666666666666666666667*(repeatedTerm3)*h)/repeatedTerm1;
repeatedTerm10 = (0.16666666666666666666666666666667*(repeatedTerm4)*h)/repeatedTerm1;
repeatedTerm11 = (0.16666666666666666666666666666667*(repeatedTerm5)*h)/repeatedTerm1;



repeatedTerm12 = (0.083333333333333333333333333333333*(repeatedTerm3)*(repeatedTerm4)*h);
repeatedTerm13 = (0.083333333333333333333333333333333*(repeatedTerm3)*(repeatedTerm5)*h);
repeatedTerm14 = (0.083333333333333333333333333333333*(repeatedTerm4)*(repeatedTerm5)*h);

repeatedTerm15 = sqrt(pow(f001 - repeatedTerm9,2) + pow(f010 - repeatedTerm10,2) + pow(f100 - repeatedTerm11,2));

repeatedTerm16 = (f011 + repeatedTerm12/repeatedTerm2 - repeatedTerm6);
repeatedTerm17  = (f101 + repeatedTerm13/repeatedTerm2 - repeatedTerm7);


distancedImplicit = 0.5*(f - 1.*repeatedTerm1*h) - 4.*(f - 0.5*repeatedTerm1*h - 0.5*h*sqrt(pow(f001 - (0.25*(repeatedTerm3)*h)/repeatedTerm1,2) + pow(f010 - (0.25*(repeatedTerm4)*h)/repeatedTerm1,2) + pow(f100 - (0.25*(repeatedTerm5)*h)/repeatedTerm1,2))) + 
4.5*(f - 0.33333333333333333333333333333333*repeatedTerm1*h - 0.33333333333333333333333333333333*h*repeatedTerm15 - 0.33333333333333333333333333333333*h*sqrt(pow(f001 - repeatedTerm9 - (0.16666666666666666666666666666667*h*(2.*(f001 - repeatedTerm9)*(f002 + (0.083333333333333333333333333333333*pow(repeatedTerm3,2)*h)/repeatedTerm2 - (0.16666666666666666666666666666667*(2.*pow(f002,2) + 2.*f001*f003 + 2.*pow(f011,2) + 2.*f010*f012 + 2.*pow(f101,2) + 2.*f100*f102)*h)/repeatedTerm1) + 2.*(f010 - repeatedTerm10)*repeatedTerm16 + 2.*(f100 - repeatedTerm11)*repeatedTerm17))/repeatedTerm15,2) + pow(f010 - repeatedTerm10 - (0.16666666666666666666666666666667*h*(2.*(f001 - repeatedTerm9)*repeatedTerm16 + 2.*(f010 - repeatedTerm10)*(f020 + (0.083333333333333333333333333333333*pow(repeatedTerm4,2)*h)/repeatedTerm2 - (0.16666666666666666666666666666667*(2.*pow(f011,2) + 2.*pow(f020,2) + 2.*f001*f021 + 2.*f010*f030 + 2.*pow(f110,2) + 2.*f100*f120)*h)/repeatedTerm1) + 2.*(f100 - repeatedTerm11)*(f110 + repeatedTerm14/repeatedTerm2 - repeatedTerm8)))/repeatedTerm15,2) + pow(f100 - repeatedTerm11 - (0.16666666666666666666666666666667*h*(2.*(f001 - repeatedTerm9)*repeatedTerm17 + 2.*(f010 - repeatedTerm10)*(f110 + repeatedTerm14/repeatedTerm2 - repeatedTerm8) + 2.*(f100 - repeatedTerm11)*(f200 + (0.083333333333333333333333333333333*pow(repeatedTerm5,2)*h)/repeatedTerm2 - (0.16666666666666666666666666666667*(2.*pow(f101,2) + 2.*pow(f110,2) + 2.*pow(f200,2) + 2.*f001*f201 + 2.*f010*f210 + 2.*f100*f300)*h)/repeatedTerm1)))/repeatedTerm15,2)));
	
		return distancedImplicit;
	}
};
/*==========================================================*/
bool GbGeneralThirdOrderImplicitSurface::isPointInGbObject3D(const double& x1, const double& x2, const double& x3)
{
	//double f = sin(2.*M_PI*x1/edgeLength)*cos(2.*M_PI*x2 / edgeLength) + sin(2.*M_PI*x2 / edgeLength)*cos(2.*M_PI*x3 / edgeLength) + sin(2.*M_PI*x3 / edgeLength)*cos(2.*M_PI*x1 / edgeLength);
	//evaluateImplicitFunction(x1,x2,x3, 0., 0., 0.)
	double f1 = evaluateImplicitFunction(x1, x2, x3, 1.);
	double f2 = evaluateImplicitFunction(x1, x2, x3, -1.);
	// 	if (f < 10.0E-15 && f > -10.0E-15)
		//if (fabs(f) <= 10e-15)
	 //if (f <= 0)
	 if (UbMath::less(x1, this->getX1Minimum()))
        return false;
    else if (UbMath::less(x2, this->getX2Minimum()))
        return false;
    else if (UbMath::less(x3, this->getX3Minimum()))
        return false;
    else if (UbMath::greater(x1, this->getX1Maximum()))
        return false;
    else if (UbMath::greater(x2, this->getX2Maximum()))
        return false;
    else if (UbMath::greater(x3, this->getX3Maximum()))
        return false;
	 if (UbMath::lessEqual(f1,0.) && UbMath::greaterEqual(f2,0.) )
		{
			return true;
		}
		else
		{
			return false;
		} 

}

/*==========================================================*/
double GbGeneralThirdOrderImplicitSurface::getIntersectionRaytraceFactor(const double& x1, const double& x2, const double& x3, const double& rx1, const double& rx2, const double& rx3)
{
	double from = 0;  // The solution must lie in the interval [from, to], additionally f(from) <= 0 && f(to) >= 0
	double to = dx*sqrt(rx1*rx1+ rx2*rx2+ rx3*rx3);
	FunctionGyroidThirdOrder f;
	//FunctionToApproximate f;
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
double GbGeneralThirdOrderImplicitSurface::evaluateImplicitFunction(const double& x1, const double& x2, const double& x3, const double& position)
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
double GbGeneralThirdOrderImplicitSurface::getX1Centroid()
{
	return (0.5*(p1->x1 + p2->x1));
}
/*=======================================================*/
double GbGeneralThirdOrderImplicitSurface::getX1Minimum()
{
	return (this->p1->x1 < this->p2->x1 ? this->p1->x1 : this->p2->x1);
}
/*=======================================================*/
double GbGeneralThirdOrderImplicitSurface::getX1Maximum()
{
	return (this->p1->x1 > this->p2->x1 ? this->p1->x1 : this->p2->x1);
}
/*=======================================================*/
double GbGeneralThirdOrderImplicitSurface::getX2Centroid()
{
	return (0.5*(p1->x2 + p2->x2));
}
/*=======================================================*/
double GbGeneralThirdOrderImplicitSurface::getX2Minimum()
{
	return (this->p1->x2 < this->p2->x2 ? this->p1->x2 : this->p2->x2);
}
/*=======================================================*/
double GbGeneralThirdOrderImplicitSurface::getX2Maximum()
{
	return (this->p1->x2 > this->p2->x2 ? this->p1->x2 : this->p2->x2);
}
/*=======================================================*/
double GbGeneralThirdOrderImplicitSurface::getX3Centroid()
{
	return (0.5*(p1->x3 + p2->x3));
}
/*=======================================================*/
double GbGeneralThirdOrderImplicitSurface::getX3Minimum()
{
	return (this->p1->x3 < this->p2->x3 ? this->p1->x3 : this->p2->x3);
}
/*=======================================================*/
double GbGeneralThirdOrderImplicitSurface::getX3Maximum()
{
	return (this->p1->x3 > this->p2->x3 ? this->p1->x3 : this->p2->x3);
}
/*=======================================================*/
bool GbGeneralThirdOrderImplicitSurface::isCellInsideGbObject3D(const double& x1a, const double& x2a, const double& x3a, const double& x1b, const double& x2b, const double& x3b)
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
bool GbGeneralThirdOrderImplicitSurface::isCellInsideOrCuttingGbObject3D(const double& x1a, const double& x2a, const double& x3a, const double& x1b, const double& x2b, const double& x3b)
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
bool GbGeneralThirdOrderImplicitSurface::isCellCuttingGbObject3D(const double& x1a, const double& x2a, const double& x3a, const double& x1b, const double& x2b, const double& x3b)
{
	if (!this->isCellInsideGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b)
		&& this->isCellInsideOrCuttingGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b))
	{
		return true;
	}
	return false;
}
/*=======================================================*/
void GbGeneralThirdOrderImplicitSurface::addSurfaceTriangleSet(vector<UbTupleFloat3>& nodes, vector<UbTupleInt3>& triangles)
{
	/*0*/nodes.push_back(makeUbTuple((float)getX1Minimum(), (float)getX2Minimum(), (float)getX3Minimum()));
	/*1*/nodes.push_back(makeUbTuple((float)getX1Maximum(), (float)getX2Minimum(), (float)getX3Minimum()));
	/*2*/nodes.push_back(makeUbTuple((float)getX1Maximum(), (float)getX2Maximum(), (float)getX3Minimum()));
	/*3.*/nodes.push_back(makeUbTuple((float)getX1Minimum(), (float)getX2Maximum(), (float)getX3Minimum()));

	/*4*/nodes.push_back(makeUbTuple((float)getX1Minimum(), (float)getX2Minimum(), (float)getX3Maximum()));
	/*5*/nodes.push_back(makeUbTuple((float)getX1Maximum(), (float)getX2Minimum(), (float)getX3Maximum()));
	/*6*/nodes.push_back(makeUbTuple((float)getX1Maximum(), (float)getX2Maximum(), (float)getX3Maximum()));
	/*7*/nodes.push_back(makeUbTuple((float)getX1Minimum(), (float)getX2Maximum(), (float)getX3Maximum()));

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
/*==========================================================*/
void GbGeneralThirdOrderImplicitSurface::objectChanged(UbObservable *changedObject)
{
    GbPoint3D *point = dynamic_cast<GbPoint3D *>(changedObject);
    if (!point || (this->p1 != point && this->p2 != point && this->p3 != point && this->p4 != point))
        return;

    this->notifyObserversObjectChanged();
}
/*==========================================================*/
void GbGeneralThirdOrderImplicitSurface::objectWillBeDeleted(UbObservable *objectForDeletion)
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
    // ACHTUNG: eigentlich muessten in allen methoden von GbLine if abfragen fuer NULL pointer hin... toDo
}

#endif
