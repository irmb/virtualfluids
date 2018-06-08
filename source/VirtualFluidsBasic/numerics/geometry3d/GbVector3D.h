//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef GBVECTOR3D_H
#define GBVECTOR3D_H
                                                                   
#include <cfloat>                               
#include <cassert> 
#include <string>

#ifdef CAB_RCF
   #include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif //CAB_RCF

#include <PointerDefinitions.h>

class GbPoint3D;

class GbVector3D 
{
public:
   // construction
    GbVector3D (); 
    GbVector3D (const double& fX1, const double& fX2, const double& fX3);
    GbVector3D (const GbVector3D& rkV);
    GbVector3D (const GbPoint3D& rkV);

    std::string toString();

    // coordinate access
    operator const double* () const;
    operator double* ();
    double operator[] (int i) const;
    double& operator[] (int i);
    double X1 () const;
    double& X1 ();
    double X2 () const;
    double& X2 ();                                    
    double X3 () const;
    double& X3 ();

    // assignment
    GbVector3D& operator= (const GbVector3D& rkV);

    // comparison
    bool operator== (const GbVector3D& rkV) const;
    bool operator!= (const GbVector3D& rkV) const;
    bool operator<  (const GbVector3D& rkV) const;
    bool operator<= (const GbVector3D& rkV) const;
    bool operator>  (const GbVector3D& rkV) const;
    bool operator>= (const GbVector3D& rkV) const;

    // arithmetic operations
    GbVector3D operator+ (const GbVector3D& rkV) const;
    GbVector3D operator- (const GbVector3D& rkV) const;
    GbVector3D operator* (const double& fScalar) const;
    GbVector3D operator/ (const double& fScalar) const;
    GbVector3D operator- () const;

    // arithmetic updates
    GbVector3D& operator+= (const GbVector3D& rkV);
    GbVector3D& operator-= (const GbVector3D& rkV);
    GbVector3D& operator*= (const double& fScalar);
    GbVector3D& operator/= (const double& fScalar);

    GbVector3D Add(GbVector3D& vector);
    GbVector3D Subtract(GbVector3D& vector);
    GbVector3D Scale(const double& x);

    // vector operations
    double Length () const;
    double SquaredLength () const;
    double Dot (const GbVector3D& rkV) const;
    double Normalize ();

    // The cross products are computed using the right-handed rule.  Be aware
    // that some graphics APIs use a left-handed rule.  If you have to compute
    // a cross product with these functions and send the result to the API
    // that expects left-handed, you will need to change sign on the vector
    // (replace each component value c by -c).
    GbVector3D Cross (const GbVector3D& rkV) const;
    GbVector3D UnitCross (const GbVector3D& rkV) const;

    // Compute the barycentric coordinates of the point with respect to the
    // tetrahedron <V0,V1,V2,V3>, P = b0*V0 + b1*V1 + b2*V2 + b3*V3, where
    // b0 + b1 + b2 + b3 = 1.
    void GetBarycentrics (const GbVector3D& rkV0,
                          const GbVector3D& rkV1, const GbVector3D& rkV2,
                          const GbVector3D& rkV3, double afBary[4]) const;

    // Gram-Schmidt orthonormalization.  Take linearly independent vectors
    // U, V, and W and compute an orthonormal set (unit length, mutually
    // perpendicular).
    static void Orthonormalize (GbVector3D& rkU, GbVector3D& rkV, GbVector3D& rkW);
    static void Orthonormalize (GbVector3D* akV);

    // Input W must be initialized to a nonzero vector, output is {U,V,W},
    // an orthonormal basis.  A hint is provided about whether or not W
    // is already unit length.
    static void GenerateOrthonormalBasis (GbVector3D& rkU, GbVector3D& rkV,
                                          GbVector3D& rkW, bool bUnitLengthW);

    // special vectors
    static const GbVector3D ZERO;
    static const GbVector3D UNIT_X1;
    static const GbVector3D UNIT_X2;
    static const GbVector3D UNIT_X3;

#ifdef CAB_RCF
    template<class Archive>
    void SF_SERIALIZE(Archive & ar)
    {
       ar & m_afTuple;
    }
#endif //CAB_RCF
private:
    // support for comparisons
    int CompareArrays (const GbVector3D& rkV) const;

    double m_afTuple[3];
};

GbVector3D operator* (const double& fScalar, const GbVector3D& rkV);

#ifdef RCF_USE_SF_SERIALIZATION
   UB_AUTO_RUN_NAMED(   SF::registerType<GbVector3D  >("GbVector3D  "), SF_GbVector3D     );
#endif //RCF_USE_SF_SERIALIZATION

#endif //GBVECTOR3D_H
