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
//! \file GbMeshTools3D.h
//! \ingroup geometry3D
//! \author Konstantin Kutscher, Soeren Textor, Sebastian Geller
//=======================================================================================
#ifndef GBMESHTOOLS3D_H
#define GBMESHTOOLS3D_H

#include <iostream>
#include <sstream>
#include <vector>

#include <basics/utilities/UbMath.h>

namespace GbMeshTools3D
{
inline int planeBoxOverlap(float normal[3], float vert[3], float maxbox[3]) // -NJMP-
{
    int q;
    float vmin[3], vmax[3], v;

    for (q = 0; q <= 2; q++) {
        v = vert[q]; // -NJMP-
        if (normal[q] > 0.0f) {
            vmin[q] = -maxbox[q] - v; // -NJMP-
            vmax[q] = maxbox[q] - v;  // -NJMP-
        } else {
            vmin[q] = maxbox[q] - v;  // -NJMP-
            vmax[q] = -maxbox[q] - v; // -NJMP-
        }
    }
    if ((normal[0] * vmin[0] + normal[1] * vmin[1] + normal[2] * vmin[2]) > 0.0f)
        return 0; // -NJMP-
    if ((normal[0] * vmax[0] + normal[1] * vmax[1] + normal[2] * vmax[2]) >= 0.0f)
        return 1; // -NJMP-
    return 0;
}

// Testet auf schnittpunkt Box <-> Dreieck
// boxcenter   = Mittelpunkt der Box
// boxhalfsize = Halbe Laenge/Hoehe/Breite der Box
// triverts    = Koordinaten der Deickspunkte
inline int triBoxOverlap(float boxcenter[3], float boxhalfsize[3], float triverts[3][3])
{
    /*    use separating axis theorem to test overlap between triangle and box */
    /*    need to test for overlap in these directions: */
    /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
    /*       we do not even need to test these) */
    /*    2) normal of the triangle */
    /*    3) crossproduct(edge from tri, {x,y,z}-directin) */
    /*       this gives 3x3=9 more tests */

    float v0[3], v1[3], v2[3];

    //   float axis[3];

    float min, max, p0, p1, p2, rad, fex, fey, fez; // -NJMP- "d" local variable removed
    float normal[3], e0[3], e1[3], e2[3];

    /* This is the fastest branch on Sun */
    /* move everything so that the boxcenter is in (0,0,0) */
    // SUB(v0,triverts[0],boxcenter);
    //#define SUB(dest,v1,v2)
    v0[0] = triverts[0][0] - boxcenter[0];
    v0[1] = triverts[0][1] - boxcenter[1];
    v0[2] = triverts[0][2] - boxcenter[2];

    // SUB(v1,triverts[1],boxcenter);
    //#define SUB(dest,v1,v2)
    v1[0] = triverts[1][0] - boxcenter[0];
    v1[1] = triverts[1][1] - boxcenter[1];
    v1[2] = triverts[1][2] - boxcenter[2];

    // SUB(v2,triverts[2],boxcenter);
    //#define SUB(dest,v1,v2)
    v2[0] = triverts[2][0] - boxcenter[0];
    v2[1] = triverts[2][1] - boxcenter[1];
    v2[2] = triverts[2][2] - boxcenter[2];

    /* compute triangle edges */
    // SUB(e0,v1,v0);      /* tri edge 0 */
    //#define SUB(dest,v1,v2)
    e0[0] = v1[0] - v0[0];
    e0[1] = v1[1] - v0[1];
    e0[2] = v1[2] - v0[2];

    // SUB(e1,v2,v1);      /* tri edge 1 */
    //#define SUB(dest,v1,v2)
    e1[0] = v2[0] - v1[0];
    e1[1] = v2[1] - v1[1];
    e1[2] = v2[2] - v1[2];

    // SUB(e2,v0,v2);      /* tri edge 2 */
    //#define SUB(dest,v1,v2)
    e2[0] = v0[0] - v2[0];
    e2[1] = v0[1] - v2[1];
    e2[2] = v0[2] - v2[2];

    /* Bullet 3:  */
    /*  test the 9 tests first (this was faster) */
    fex = fabsf(e0[0]);
    fey = fabsf(e0[1]);
    fez = fabsf(e0[2]);

    // AXISTEST_X01(e0[2], e0[1], fez, fey);
    //#define AXISTEST_X01(a, b, fa, fb)
    p0 = e0[2] * v0[1] - e0[1] * v0[2];
    p2 = e0[2] * v2[1] - e0[1] * v2[2];
    if (p0 < p2) {
        min = p0;
        max = p2;
    } else {
        min = p2;
        max = p0;
    }
    rad = fez * boxhalfsize[1] + fey * boxhalfsize[2];
    if (min > rad || max < -rad)
        return 0;

    // AXISTEST_Y02(e0[2], e0[0], fez, fex);
    //#define AXISTEST_Y02(a, b, fa, fb)
    p0 = -e0[2] * v0[0] + e0[0] * v0[2];
    p2 = -e0[2] * v2[0] + e0[0] * v2[2];
    if (p0 < p2) {
        min = p0;
        max = p2;
    } else {
        min = p2;
        max = p0;
    }
    rad = fez * boxhalfsize[0] + fex * boxhalfsize[2];
    if (min > rad || max < -rad)
        return 0;

    // AXISTEST_Z12(e0[1], e0[0], fey, fex);
    //#define AXISTEST_Z12(a, b, fa, fb)
    p1 = e0[1] * v1[0] - e0[0] * v1[1];
    p2 = e0[1] * v2[0] - e0[0] * v2[1];
    if (p2 < p1) {
        min = p2;
        max = p1;
    } else {
        min = p1;
        max = p2;
    }
    rad = fey * boxhalfsize[0] + fex * boxhalfsize[1];
    if (min > rad || max < -rad)
        return 0;

    fex = fabsf(e1[0]);
    fey = fabsf(e1[1]);
    fez = fabsf(e1[2]);

    // AXISTEST_X01(e1[2], e1[1], fez, fey);
    //#define AXISTEST_X01(a, b, fa, fb)
    p0 = e1[2] * v0[1] - e1[1] * v0[2];
    p2 = e1[2] * v2[1] - e1[1] * v2[2];
    if (p0 < p2) {
        min = p0;
        max = p2;
    } else {
        min = p2;
        max = p0;
    }
    rad = fez * boxhalfsize[1] + fey * boxhalfsize[2];
    if (min > rad || max < -rad)
        return 0;

    // AXISTEST_Y02(e1[2], e1[0], fez, fex);
    //#define AXISTEST_Y02(a, b, fa, fb)
    p0 = -e1[2] * v0[0] + e1[0] * v0[2];
    p2 = -e1[2] * v2[0] + e1[0] * v2[2];
    if (p0 < p2) {
        min = p0;
        max = p2;
    } else {
        min = p2;
        max = p0;
    }
    rad = fez * boxhalfsize[0] + fex * boxhalfsize[2];
    if (min > rad || max < -rad)
        return 0;

    // AXISTEST_Z0(e1[1], e1[0], fey, fex);
    //#define AXISTEST_Z0(a, b, fa, fb)
    p0 = e1[1] * v0[0] - e1[0] * v0[1];
    p1 = e1[1] * v1[0] - e1[0] * v1[1];
    if (p0 < p1) {
        min = p0;
        max = p1;
    } else {
        min = p1;
        max = p0;
    }
    rad = fey * boxhalfsize[0] + fex * boxhalfsize[2];
    if (min > rad || max < -rad)
        return 0;

    fex = fabsf(e2[0]);
    fey = fabsf(e2[1]);
    fez = fabsf(e2[2]);
    // AXISTEST_X2(e2[2], e2[1], fez, fey);
    //#define AXISTEST_X2(a, b, fa, fb)
    p0 = e2[2] * v0[1] - e2[1] * v0[2];
    p1 = e2[2] * v1[1] - e2[1] * v1[2];
    if (p0 < p1) {
        min = p0;
        max = p1;
    } else {
        min = p1;
        max = p0;
    }
    rad = fez * boxhalfsize[1] + fey * boxhalfsize[2];
    if (min > rad || max < -rad)
        return 0;

    // AXISTEST_Y1(e2[2], e2[0], fez, fex);
    //#define AXISTEST_Y1(a, b, fa, fb)
    p0 = -e2[2] * v0[0] + e2[0] * v0[2];
    p1 = -e2[2] * v1[0] + e2[0] * v1[2];
    if (p0 < p1) {
        min = p0;
        max = p1;
    } else {
        min = p1;
        max = p0;
    }
    rad = fez * boxhalfsize[0] + fex * boxhalfsize[2];
    if (min > rad || max < -rad)
        return 0;

    // AXISTEST_Z12(e2[1], e2[0], fey, fex);
    //#define AXISTEST_Z12(a, b, fa, fb)
    p1 = e2[1] * v1[0] - e2[0] * v1[1];
    p2 = e2[1] * v2[0] - e2[0] * v2[1];
    if (p2 < p1) {
        min = p2;
        max = p1;
    } else {
        min = p1;
        max = p2;
    }
    rad = fey * boxhalfsize[0] + fex * boxhalfsize[1];
    if (min > rad || max < -rad)
        return 0;

    /* Bullet 1: */
    /*  first test overlap in the {x,y,z}-directions */
    /*  find min, max of the triangle each direction, and test for overlap in */
    /*  that direction -- this is equivalent to testing a minimal AABB around */
    /*  the triangle against the AABB */
    /* test in X-direction */
    // FINDMINMAX(v0[0],v1[0],v2[0],min,max);
    min = (float)UbMath::min(v0[0], v1[0], v2[0]);
    max = (float)UbMath::max(v0[0], v1[0], v2[0]);
    if (min > boxhalfsize[0] || max < -boxhalfsize[0])
        return 0;

    /* test in Y-direction */
    // FINDMINMAX(v0[1],v1[1],v2[1],min,max);
    min = (float)UbMath::min(v0[1], v1[1], v2[1]);
    max = (float)UbMath::max(v0[1], v1[1], v2[1]);
    if (min > boxhalfsize[1] || max < -boxhalfsize[1])
        return 0;

    /* test in Z-direction */
    // FINDMINMAX(v0[2],v1[2],v2[2],min,max);
    min = (float)UbMath::min(v0[2], v1[2], v2[2]);
    max = (float)UbMath::max(v0[2], v1[2], v2[2]);

    if (min > boxhalfsize[2] || max < -boxhalfsize[2])
        return 0;

    /* Bullet 2: */
    /*  test if the box intersects the plane of the triangle */
    /*  compute plane equation of triangle: normal*x+d=0 */
    // CROSS(normal,e0,e1);
    //#define CROSS(dest,v1,v2)
    normal[0] = e0[1] * e1[2] - e0[2] * e1[1];
    normal[1] = e0[2] * e1[0] - e0[0] * e1[2];
    normal[2] = e0[0] * e1[1] - e0[1] * e1[0];

    // -NJMP- (line removed here)
    if (!planeBoxOverlap(normal, v0, boxhalfsize))
        return 0; // -NJMP-
    return 1;     /* box and triangle overlaps */
}

} // namespace GbMeshTools3D

#endif

// original - NICHT LOESCHEN - von kroete lynn

//#define X 0
//#define Y 1
//#define Z 2
//
//#define CROSS(dest,v1,v2)
//   dest[0]=v1[1]*v2[2]-v1[2]*v2[1];
//   dest[1]=v1[2]*v2[0]-v1[0]*v2[2];
//   dest[2]=v1[0]*v2[1]-v1[1]*v2[0];
//
//#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
//
//#define SUB(dest,v1,v2)
//   dest[0]=v1[0]-v2[0];
//   dest[1]=v1[1]-v2[1];
//   dest[2]=v1[2]-v2[2];
//
//#define FINDMINMAX(x0,x1,x2,min,max)
//   min = max = x0;
//   if(x1<min) min=x1;
//   if(x1>max) max=x1;
//   if(x2<min) min=x2;
//   if(x2>max) max=x2;
//
// int planeBoxOverlap(float normal[3], float vert[3], float maxbox[3])	// -NJMP-
//{
//   int q;
//   float vmin[3],vmax[3],v;
//
//   for(q=X;q<=Z;q++)
//   {
//      v=vert[q];					// -NJMP-
//      if(normal[q]>0.0f)
//      {
//         vmin[q]=-maxbox[q] - v;	// -NJMP-
//         vmax[q]= maxbox[q] - v;	// -NJMP-
//      }
//      else
//      {
//         vmin[q]= maxbox[q] - v;	// -NJMP-
//         vmax[q]=-maxbox[q] - v;	// -NJMP-
//      }
//   }
//   if(DOT(normal,vmin)>0.0f) return 0;	// -NJMP-
//   if(DOT(normal,vmax)>=0.0f) return 1;	// -NJMP-
//   return 0;
//}
//
///*======================== X-tests ========================
//
//#define AXISTEST_X01(a, b, fa, fb)
//   p0 = a*v0[Y] - b*v0[Z];
//   p2 = a*v2[Y] - b*v2[Z];
//      if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;}
//   rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];
//   if(min>rad || max<-rad) return 0;
//
//#define AXISTEST_X2(a, b, fa, fb)
//   p0 = a*v0[Y] - b*v0[Z];
//   p1 = a*v1[Y] - b*v1[Z];
//      if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;}
//   rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];
//   if(min>rad || max<-rad) return 0;
//
//
//
///*======================== Y-tests ========================
//
//#define AXISTEST_Y02(a, b, fa, fb)
//   p0 = -a*v0[X] + b*v0[Z];
//   p2 = -a*v2[X] + b*v2[Z];
//   if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;}
//   rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];
//   if(min>rad || max<-rad) return 0;
//
//#define AXISTEST_Y1(a, b, fa, fb)
//   p0 = -a*v0[X] + b*v0[Z];
//   p1 = -a*v1[X] + b*v1[Z];
//      if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;}
//   rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];
//   if(min>rad || max<-rad) return 0;
//
//======================== Z-tests ========================
//
//#define AXISTEST_Z12(a, b, fa, fb)
//   p1 = a*v1[X] - b*v1[Y];
//   p2 = a*v2[X] - b*v2[Y];
//      if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;}
//   rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];
//   if(min>rad || max<-rad) return 0;
//
//#define AXISTEST_Z0(a, b, fa, fb)
//   p0 = a*v0[X] - b*v0[Y];
//   p1 = a*v1[X] - b*v1[Y];
//      if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;}
//   rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];
//   if(min>rad || max<-rad) return 0;
//
// int triBoxOverlap(float boxcenter[3],float boxhalfsize[3],float triverts[3][3])
//{
//      //use separating axis theorem to test overlap between triangle and box
//      //need to test for overlap in these directions:
//      //1) the {x,y,z}-directions (actually, since we use the AABB of the triangle
//      //   we do not even need to test these)
//      //2) normal of the triangle
//      //3) crossproduct(edge from tri, {x,y,z}-directin)
//      //   this gives 3x3=9 more tests
//
//   float v0[3],v1[3],v2[3];
//
//   //   float axis[3];
//
//   float min,max,p0,p1,p2,rad,fex,fey,fez;		// -NJMP- "d" local variable removed
//   float normal[3],e0[3],e1[3],e2[3];
//
//   /* This is the fastest branch on Sun */
//   /* move everything so that the boxcenter is in (0,0,0) */
//   SUB(v0,triverts[0],boxcenter);
//   SUB(v1,triverts[1],boxcenter);
//   SUB(v2,triverts[2],boxcenter);
//
//   /* compute triangle edges */
//   SUB(e0,v1,v0);      /* tri edge 0 */
//   SUB(e1,v2,v1);      /* tri edge 1 */
//   SUB(e2,v0,v2);      /* tri edge 2 */
//
//   /* Bullet 3:  */
//   /*  test the 9 tests first (this was faster) */
//   fex = fabsf(e0[X]);
//   fey = fabsf(e0[Y]);
//   fez = fabsf(e0[Z]);
//
//   AXISTEST_X01(e0[Z], e0[Y], fez, fey);
//   AXISTEST_Y02(e0[Z], e0[X], fez, fex);
//   AXISTEST_Z12(e0[Y], e0[X], fey, fex);
//   fex = fabsf(e1[X]);
//   fey = fabsf(e1[Y]);
//   fez = fabsf(e1[Z]);
//
//   AXISTEST_X01(e1[Z], e1[Y], fez, fey);
//   AXISTEST_Y02(e1[Z], e1[X], fez, fex);
//   AXISTEST_Z0(e1[Y], e1[X], fey, fex);
//
//   fex = fabsf(e2[X]);
//   fey = fabsf(e2[Y]);
//   fez = fabsf(e2[Z]);
//   AXISTEST_X2(e2[Z], e2[Y], fez, fey);
//   AXISTEST_Y1(e2[Z], e2[X], fez, fex);
//   AXISTEST_Z12(e2[Y], e2[X], fey, fex);
//
//   /* Bullet 1: */
//   /*  first test overlap in the {x,y,z}-directions */
//   /*  find min, max of the triangle each direction, and test for overlap in */
//   /*  that direction -- this is equivalent to testing a minimal AABB around */
//   /*  the triangle against the AABB */
//   /* test in X-direction */
//   FINDMINMAX(v0[X],v1[X],v2[X],min,max);
//   if(min>boxhalfsize[X] || max<-boxhalfsize[X]) return 0;
//
//   /* test in Y-direction */
//   FINDMINMAX(v0[Y],v1[Y],v2[Y],min,max);
//   if(min>boxhalfsize[Y] || max<-boxhalfsize[Y]) return 0;
//
//   /* test in Z-direction */
//   FINDMINMAX(v0[Z],v1[Z],v2[Z],min,max);
//   if(min>boxhalfsize[Z] || max<-boxhalfsize[Z]) return 0;
//
//   /* Bullet 2: */
//   /*  test if the box intersects the plane of the triangle */
//   /*  compute plane equation of triangle: normal*x+d=0 */
//   CROSS(normal,e0,e1);
//
//   // -NJMP- (line removed here)
//   if(!planeBoxOverlap(normal,v0,boxhalfsize)) return 0;	// -NJMP-
//   return 1;   /* box and triangle overlaps */
//}
