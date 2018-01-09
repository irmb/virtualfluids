//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef VECTOR_3D_H
#define VECTOR_3D_H

#include <string>

/*=========================================================================*/
/*  Vector3D                                                             */
/*                                                                         */
/**
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.4 - 04.10.07
*/ 

/*
usage: ...
*/
#include "VirtualFluidsBasics_EXPORT.h"

   class VirtualFluidsBasics_EXPORT Vector3D
   {
   public:
      // construction
      Vector3D(); 
      Vector3D(const double& fX1, const double& fX2, const double& fX3);
      Vector3D(const Vector3D& rkV);

      std::string toString() const;

      // coordinate access
      operator const double*() const;
      operator double*();
      double   operator[](const int& i) const;
      double&  operator[](const int& i);
      double   X1() const;
      double&  X1();
      double   X2() const;
      double&  X2();                                    
      double   X3() const;
      double&  X3();

      // assignment
      Vector3D& operator=(const Vector3D& rkV);

      // comparison
      bool operator==(const Vector3D& rkV) const;
      bool operator!=(const Vector3D& rkV) const;
      bool operator< (const Vector3D& rkV) const;
      bool operator<=(const Vector3D& rkV) const;
      bool operator> (const Vector3D& rkV) const;
      bool operator>=(const Vector3D& rkV) const;

      // arithmetic operations
      Vector3D operator+(const Vector3D& rkV) const;
      Vector3D operator-(const Vector3D& rkV) const;
      Vector3D operator*(const double& fScalar) const;
      Vector3D operator/(const double& fScalar) const;
      Vector3D operator-() const;

      // arithmetic updates
      Vector3D& operator+= (const Vector3D& rkV);
      Vector3D& operator-= (const Vector3D& rkV);
      Vector3D& operator*= (const double& fScalar);
      Vector3D& operator/= (const double& fScalar);

      Vector3D Add(Vector3D& vector);
      Vector3D Subtract(Vector3D& vector);
      Vector3D Scale(const double& x);

      // vector operations
      double Length () const;
      double SquaredLength () const;
      double Dot (const Vector3D& rkV) const;
      double Normalize ();

      // The cross products are computed using the right-handed rule.  Be aware
      // that some graphics APIs use a left-handed rule.  If you have to compute
      // a cross product with these functions and send the result to the API
      // that expects left-handed, you will need to change sign on the vector
      // (replace each component value c by -c).
      Vector3D Cross (const Vector3D& rkV) const;
      Vector3D UnitCross (const Vector3D& rkV) const;

      // Compute the barycentric coordinates of the point with respect to the
      // tetrahedron <V0,V1,V2,V3>, P = b0*V0 + b1*V1 + b2*V2 + b3*V3, where
      // b0 + b1 + b2 + b3 = 1.
      void GetBarycentrics (const Vector3D& rkV0, const Vector3D& rkV1, const Vector3D& rkV2, const Vector3D& rkV3, double afBary[4]) const;

      // Gram-Schmidt orthonormalization.  Take linearly independent vectors
      // U, V, and W and compute an orthonormal set (unit length, mutually
      // perpendicular).
      static void Orthonormalize (Vector3D& rkU, Vector3D& rkV, Vector3D& rkW);
      static void Orthonormalize (Vector3D* akV);

      // Input W must be initialized to a nonzero vector, output is {U,V,W},
      // an orthonormal basis.  A hint is provided about whether or not W
      // is already unit length.
      static void GenerateOrthonormalBasis (Vector3D& rkU, Vector3D& rkV, Vector3D& rkW, bool bUnitLengthW);

      // special vectors
      static const Vector3D ZERO;
      static const Vector3D UNIT_X1;
      static const Vector3D UNIT_X2;
      static const Vector3D UNIT_X3;

   #ifdef CAB_RCF
         template<class Archive>
         void serialize(Archive & ar, const unsigned int version)
         {
            ar & m_afTuple;
         }
   #endif //CAB_RCF

   protected:
      // support for comparisons
      int CompareArrays (const Vector3D& rkV) const;

      double m_afTuple[3];
   };
   
   //globaler multiplaktor mit skalar
   Vector3D operator*(const double& fScalar, const Vector3D& rkV);
   std::ostream& operator<<(std::ostream& os, const Vector3D& rkV);

  
#endif
