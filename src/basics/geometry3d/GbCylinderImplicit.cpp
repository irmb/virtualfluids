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
//! \author Hussein Alihussein
//=======================================================================================

#include <GbCylinderImplicit.h>

#ifdef VF_BOOST

#include <basics/utilities/UbMath.h>

#include <geometry3d/GbSystem3D.h>
#include <geometry3d/GbTriangle3D.h>

#include <boost/math/tools/roots.hpp>


using namespace std;
using boost::math::tools::bisect;


/*=======================================================*/
GbCylinderImplicit::GbCylinderImplicit() 
{

}
/*=======================================================*/
// Konstruktor
GbCylinderImplicit::GbCylinderImplicit(const double& x1a, const double& x2a, const double& x3a, const double& x1b, const double& x2b, const double& x3b ,const double& dx,  const double& OuterRadius, const double& InnerRadius) :GbObject3D()
{
    this->p1 = new GbPoint3D(x1a, x2a, x3a);
    this->p2 = new GbPoint3D(x1b, x2b, x3b);
    this->p1->addObserver(this);
    this->p2->addObserver(this);
    
    //this->p3 = new GbPoint3D(x1Lbound, x2Lbound, x3Lbound);
    //this->p4 = new GbPoint3D(x1Ubound, x2Ubound, x3Ubound);

    this->dx = dx;
    this->OuterRadius=OuterRadius;
    this->InnerRadius=InnerRadius;


    this->f.x0 = x1a;
    this->f.y0 = x2a;
    this->f.z0 = x3a;
    this->f.x1 = x1b;
    this->f.y1 = x2b;
    this->f.z1 = x3b;
    this->f.r=InnerRadius;
    this->f.R=OuterRadius;

    this-> f.init();
    this->cylinderLengthSquared = pow(x1a - x1b,2) + pow(x2a - x2b,2) + pow(x3a - x3b,2);
    //base-circles local coordinate system
    this->e1x = (pow(x2a - x2b,2) + pow(x3a - x3b,2) + x1b*(x2a - x2b + x3a - x3b) + x1a*(-x2a + x2b - x3a + x3b))/(this->cylinderLengthSquared*sqrt(2 + (2*(-((x2a - x2b)*(x3a - x3b)) + x1b*(x2a - x2b + x3a - x3b) + x1a*(-x2a + x2b - x3a + x3b)))/
       this->cylinderLengthSquared));



    this->e1y= (pow(x1a,2) + pow(x1b,2) + x1b*(x2a - x2b) + x1a*(-2*x1b - x2a + x2b) + (x3a - x3b)*(-x2a + x2b + x3a - x3b))/(this->cylinderLengthSquared*sqrt(2 + (2*(-((x2a - x2b)*(x3a - x3b)) + x1b*(x2a - x2b + x3a - x3b) + x1a*(-x2a + x2b - x3a + x3b)))/
       this->cylinderLengthSquared));



    this->e1z = (pow(x1a,2) + pow(x1b,2) + x1b*(x3a - x3b) + x1a*(-2*x1b - x3a + x3b) + (x2a - x2b)*(x2a - x2b - x3a + x3b))/(this->cylinderLengthSquared*sqrt(2 + (2*(-((x2a - x2b)*(x3a - x3b)) + x1b*(x2a - x2b + x3a - x3b) + x1a*(-x2a + x2b - x3a + x3b)))/
       this->cylinderLengthSquared));



    this->e2x = (-x2a + x2b + x3a - x3b)/(sqrt(2)*sqrt(this->cylinderLengthSquared)*
     sqrt(1 + (-((x2a - x2b)*(x3a - x3b)) + x1b*(x2a - x2b + x3a - x3b) + x1a*(-x2a + x2b - x3a + x3b))/(this->cylinderLengthSquared)));



    this->e2y= (x1a - x1b - x3a + x3b)/(sqrt(2)*sqrt(this->cylinderLengthSquared)*
     sqrt(1 + (-((x2a - x2b)*(x3a - x3b)) + x1b*(x2a - x2b + x3a - x3b) + x1a*(-x2a + x2b - x3a + x3b))/(this->cylinderLengthSquared)));




    this->e2z = (-x1a + x1b + x2a - x2b)/(sqrt(2)*sqrt(this->cylinderLengthSquared)*
     sqrt(1 + (-((x2a - x2b)*(x3a - x3b)) + x1b*(x2a - x2b + x3a - x3b) + x1a*(-x2a + x2b - x3a + x3b))/(this->cylinderLengthSquared)));
  



    
}
/*=======================================================*/
GbCylinderImplicit::GbCylinderImplicit(GbCylinderImplicit * imp)
{
}
/*=======================================================*/
// Destruktor
GbCylinderImplicit::~GbCylinderImplicit()
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
// struct TerminationCondition {
//     bool operator() (double min, double max) {
//         return abs(min - max) <= 10e-10;
//     }
// };
// /*=======================================================*/
// struct FunctionOfImplicitCylinder {
//     double x, y, z;
//     double x0, y0, z0;
//     double x1, y1, z1;   
//     double r;
    
//     double l = sqrt(pow(x0 - x1,2) + pow(y0 - y1,2) + pow(z0 - z1,2));

//     double f1 = pow(x - x0,2) + pow(y - y0,2) + pow(z - z0,2) - 
//     (pow(r,2) + pow((x - x0)*(x0 - x1) + (y - y0)*(y0 - y1) + (z - z0)*(z0 - z1),2)/(l*l));
//     double f2= -(
//    sqrt(pow(x0,2) + pow(y0,2) + pow(z0,2)) + (x*(-x0 + x1) + y*(-y0 + y1) + z*(-z0 + z1))/l );

//     double f3=(x*(-x0 + x1) + y*(-y0 + y1) + z*(-z0 + z1))/l - sqrt(pow(x1,2) + pow(y1,2) + pow(z1,2));




//     bool operator() () {
// //    return pow(x - x0,2) + pow(y - y0,2) + pow(z - z0,2) <= 
// //     pow(r,2) + pow((x - x0)*(x0 - x1) + (y - y0)*(y0 - y1) + (z - z0)*(z0 - z1),2)/(l*l) && 
// //    sqrt(pow(x0,2) + pow(y0,2) + pow(z0,2)) + (x*(-x0 + x1) + y*(-y0 + y1) + z*(-z0 + z1))/l >= 0 && 
// //    (x*(-x0 + x1) + y*(-y0 + y1) + z*(-z0 + z1))/l <= sqrt(pow(x1,2) + pow(y1,2) + pow(z1,2));
// return f1<=0 && f2<=0 && f3<=0;
    
//     }
// };
/*==========================================================*/
bool GbCylinderImplicit::isPointInGbObject3D(const double& x1, const double& x2, const double& x3)
{
    this->f.x = x1;
    this->f.y = x2;
    this->f.z = x3;
    
    double t1=this->f.f1(0);
    //double t2=this->f.f2(0);
    //double t3=this->f.f3(0);

    return t1<=0 ;
    //return t1<=0 && t2<=0 && t3<=0;
}
/*==========================================================*/
bool GbCylinderImplicit::isPointInGbObject3D(const double& x1, const double& x2, const double& x3, bool &pointIsOnBoundary)
{
    this->f.x = x1;
    this->f.y = x2;
    this->f.z = x3;
    
    double t1=this->f.f1(0);
    //double t2=this->f.f2(0);
    //double t3=this->f.f3(0);


    //pointIsOnBoundary = ub_math::equal(t1,0.) && ub_math::equal(t2,0.) && ub_math::equal(t3,0.);  
    pointIsOnBoundary = ub_math::equal(t1,0.);  
 

    return t1<=0 ;
    //return t1<=0 && t2<=0 && t3<=0;

}
/*==========================================================*/
double GbCylinderImplicit::getIntersectionRaytraceFactor(const double& x1, const double& x2, const double& x3, const double& rx1, const double& rx2, const double& rx3)
{
    //return 1.;
    double from = 0;  // The solution must lie in the interval [from, to], additionally f(from) <= 0 && f(to) >= 0
    double to = dx/* *sqrt(rx1*rx1+ rx2*rx2+ rx3*rx3) */;
    //FunctionOfImplicitCylinder g=this->f;
    f.x =x1 ;
    f.y =x2 ;
    f.z =x3 ;

    f.dir1 = rx1/sqrt(rx1*rx1+ rx2*rx2+ rx3*rx3);
    f.dir2 = rx2/sqrt(rx1*rx1+ rx2*rx2+ rx3*rx3);
    f.dir3 = rx3/sqrt(rx1*rx1+ rx2*rx2+ rx3*rx3);
    
    double f1= f.f1(0);
    //double f2=f.f2(0);
    //double f3=f.f3(0);
    
   
   double g1= f.f1(to);
   //double g2= f.f2(to);
   //double g3= f.f3(to);

//int ss=std::signbit(f1);

//bool test1 = f1<0 && f2<0 && f3<0;
//bool test2 = g1<0 && g2<0 && g3<0;
//bool test3 = f1>0 && f2>0 && f3>0;
//bool test4 = g1>0 && g2>0 && g3>0;


//if(test1 || test2 || test3 || test4)
//{    
 

    
double p1=f1*g1;
//double p2=f2*g2;
//double p3=f3*g3;
if (p1>0) {
    return 999.;
}
else {
std::pair<double, double> result1 = bisect(f.f1, from, to, TerminationCondition());
           double root1 = (result1.first + result1.second) / 2;
            return root1;
}


// if (p1<0) {
//         if (p2 < 0) {
//             //std::cout << "f1&f2" << std::endl;
//             std::pair<double, double> result1 = bisect(f.f1, from, to, TerminationCondition());
//             double root1 = (result1.first + result1.second) / 2;
//             std::pair<double, double> result2 = bisect(f.f2, from, to, TerminationCondition());
//             double root2 = (result2.first + result2.second) / 2;
//             return UbMath::max(root1,root2);
//         } else if (p3 < 0) {
//             //std::cout << "f1&f3" << std::endl;
//             std::pair<double, double> result1 = bisect(f.f1, from, to, TerminationCondition());
//             double root1 = (result1.first + result1.second) / 2;
//             std::pair<double, double> result2 = bisect(f.f3, from, to, TerminationCondition());
//             double root2 = (result2.first + result2.second) / 2;
//             return UbMath::max(root1,root2);
//         } else {
//             //std::cout << "f1" << std::endl;
//             std::pair<double, double> result1 = bisect(f.f1, from, to, TerminationCondition());
//             double root1 = (result1.first + result1.second) / 2;
//             return root1;
//         }
//     } else  {
//         if (p3 < 0) {
//             //std::cout << "f3" << std::endl;
//             //std::pair<double, double> result1 = bisect(f.f2, from, to, TerminationCondition());
//             //double root1 = (result1.first + result1.second) / 2;
//             std::pair<double, double> result2 = bisect(f.f3, from, to, TerminationCondition());
//             double root2 = (result2.first + result2.second) / 2;
//             return root2;
//         } else if (p2 < 0)  {
//             //std::cout << "f2" << std::endl;
//             std::pair<double, double> result1 = bisect(f.f2, from, to, TerminationCondition());
//             double root1 = (result1.first + result1.second) / 2;
//             return root1;
//         } else {
//         //std::cout << "1" << std::endl;
//         return 999.;
//     }
//     }     
// }
// else
//  {
//     return 999.;
//  }
}
/*=======================================================*/
// double GbCylinderImplicit::evaluateImplicitFunction(const double& x1, const double& x2, const double& x3)
// {
//     FunctionOfImplicitCylinder f;
//     f.x = x1;
//     f.y = x2;
//     f.z = x3;

//     // f.x0 = p1->x1;
//     // f.y0 = p1->x2;
//     // f.z0 = p1->x3;
//     // f.x1 = p2->x1;
//     // f.y1 = p2->x2;
//     // f.z1 = p2->x3;

//     return f();
// }
/*=======================================================*/
double GbCylinderImplicit::getX1Centroid()
{
    return (0.5*(p1->x1 + p2->x1));
}
/*=======================================================*/
double GbCylinderImplicit::getX1Minimum()
{
    //return (this->p1->x1 < this->p2->x1 ? this->p1->x1 : this->p2->x1);
    //return this->p3->x1;
    return ub_math::min(-(sqrt(pow(e1x,2) + pow(e2x,2))*this->OuterRadius) + p1->x1,
                        -(sqrt(pow(e1x,2) + pow(e2x,2))*this->OuterRadius) + p2->x1);
}
/*=======================================================*/
double GbCylinderImplicit::getX1Maximum()
{
    //return (this->p1->x1 > this->p2->x1 ? this->p1->x1 : this->p2->x1);
    //return this->p4->x1;
    return ub_math::max(sqrt(pow(e1x,2) + pow(e2x,2))*this->OuterRadius + p1->x1,
                        sqrt(pow(e1x,2) + pow(e2x,2))*this->OuterRadius + p2->x1);

}
/*=======================================================*/
double GbCylinderImplicit::getX2Centroid()
{
    return (0.5*(p1->x2 + p2->x2));
}
/*=======================================================*/
double GbCylinderImplicit::getX2Minimum()
{
    //return (this->p1->x2 < this->p2->x2 ? this->p1->x2 : this->p2->x2);
    //return this->p3->x2;
    return ub_math::min(-(sqrt(pow(e1y,2) + pow(e2y,2))*this->OuterRadius) + p1->x2,
                        -(sqrt(pow(e1y,2) + pow(e2y,2))*this->OuterRadius) + p2->x2);

}
/*=======================================================*/
double GbCylinderImplicit::getX2Maximum()
{
    //return (this->p1->x2 > this->p2->x2 ? this->p1->x2 : this->p2->x2);
    //return this->p4->x2;
    return ub_math::max(sqrt(pow(e1y,2) + pow(e2y,2))*this->OuterRadius + p1->x2,
                        sqrt(pow(e1y,2) + pow(e2y,2))*this->OuterRadius + p2->x2);

}
/*=======================================================*/
double GbCylinderImplicit::getX3Centroid()
{
    return (0.5*(p1->x3 + p2->x3));
}
/*=======================================================*/
double GbCylinderImplicit::getX3Minimum()
{
    //return (this->p1->x3 < this->p2->x3 ? this->p1->x3 : this->p2->x3);
    //return this->p3->x3;
    return ub_math::min(-(sqrt(pow(e1z,2) + pow(e2z,2))*this->OuterRadius) + p1->x3,
                        -(sqrt(pow(e1z,2) + pow(e2z,2))*this->OuterRadius) + p2->x3);


}
/*=======================================================*/
double GbCylinderImplicit::getX3Maximum()
{
    //return (this->p1->x3 > this->p2->x3 ? this->p1->x3 : this->p2->x3);
    //return this->p4->x3;
    return ub_math::max(sqrt(pow(e1z,2) + pow(e2z,2))*this->OuterRadius + p1->x3,
                        sqrt(pow(e1z,2) + pow(e2z,2))*this->OuterRadius + p2->x3);
}
/*=======================================================*/
bool GbCylinderImplicit::isCellInsideGbObject3D(const double& x1a, const double& x2a, const double& x3a, const double& x1b, const double& x2b, const double& x3b)
{
    if (this->isPointInGbObject3D(x1a, x2a, x3a)
        && this->isPointInGbObject3D(x1b, x2a, x3a)
        && this->isPointInGbObject3D(x1b, x2b, x3a)
        && this->isPointInGbObject3D(x1a, x2b, x3a)
        && this->isPointInGbObject3D(x1a, x2a, x3b)
        && this->isPointInGbObject3D(x1b, x2a, x3b)
        && this->isPointInGbObject3D(x1b, x2b, x3b)
        && this->isPointInGbObject3D(x1a, x2b, x3b)
        )
    {
        return true;
    }
    return false;
}
/*=======================================================*/
bool GbCylinderImplicit::isCellInsideOrCuttingGbObject3D(const double& x1a, const double& x2a, const double& x3a, const double& x1b, const double& x2b, const double& x3b)
{
    if ((this->isPointInGbObject3D(x1a, x2a, x3a) == false)
        && (this->isPointInGbObject3D(x1b, x2a, x3a) == false)
        && (this->isPointInGbObject3D(x1b, x2b, x3a) == false)
        && (this->isPointInGbObject3D(x1a, x2b, x3a) == false)
        && (this->isPointInGbObject3D(x1a, x2a, x3b) == false)
        && (this->isPointInGbObject3D(x1b, x2a, x3b) == false)
        && (this->isPointInGbObject3D(x1b, x2b, x3b) == false)
        && (this->isPointInGbObject3D(x1a, x2b, x3b) == false)
        )
    {
        return false;
    }
    return true;
}
/*=======================================================*/
bool GbCylinderImplicit::isCellCuttingGbObject3D(const double& x1a, const double& x2a, const double& x3a, const double& x1b, const double& x2b, const double& x3b)
{
    if (!this->isCellInsideGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b)
        && this->isCellInsideOrCuttingGbObject3D(x1a, x2a, x3a, x1b, x2b, x3b))
    {
        return true;
    }
    return false;
}
/*==========================================================*/
void GbCylinderImplicit::objectChanged(UbObservable *changedObject)
{
    GbPoint3D *point = dynamic_cast<GbPoint3D *>(changedObject);
    if (!point || (this->p1 != point && this->p2 != point && this->p3 != point && this->p4 != point))
        return;

    this->notifyObserversObjectChanged();
}
/*==========================================================*/
void GbCylinderImplicit::objectWillBeDeleted(UbObservable *objectForDeletion)
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

//! \}
