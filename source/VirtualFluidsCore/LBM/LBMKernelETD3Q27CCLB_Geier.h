//Cascaded Cumulant LBM

#ifndef LBMKernelETD3Q27CCLB_Geier_H
#define LBMKernelETD3Q27CCLB_Geier_H

#include "LBMKernelETD3Q27.h"
#include "D3Q27ETBCProcessor.h"
#include "D3Q27System.h"
#include <boost/serialization/export.hpp>
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

class LBMKernelETD3Q27CCLB_Geier;
typedef boost::shared_ptr<LBMKernelETD3Q27CCLB_Geier> LBMKernelETD3Q27CCLB_GeierPtr;


class LBMKernelETD3Q27CCLB_Geier :  public LBMKernelETD3Q27
{
public:
   LBMKernelETD3Q27CCLB_Geier();
   LBMKernelETD3Q27CCLB_Geier(int nx1, int nx2, int nx3, int option);
   ~LBMKernelETD3Q27CCLB_Geier(void);
   void calculate();
   LBMKernel3DPtr clone();
   double getCallculationTime();

protected:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<LBMKernelETD3Q27>(*this);
      ar & OxyyMxzz; 
      ar & option;
   }

   void collideAll();  
   void collideAll2(); 
   void collideAll3();
   void init();
   LBMReal f[D3Q27System::ENDF+1];
   LBMReal rho, vx, vy, vz, durchrho;
   LBMReal w1,w2,w3,w4,w5,w6,w7,w8,w9,w10;
   LBMReal s9,c1o27,c2o3;

   LBMReal M_zXX, M_zYY, M_zZZ, M_zXY,    M_zXZ,  M_zYZ,
      M_zXXY,    M_zXYY,    M_zXXZ,    M_zXZZ,   M_zYYZ,  M_zYZZ,  M_zXYZ,
      M_zXXYY,   M_zXXZZ,   M_zYYZZ,   M_zXXYZ,  M_zXYYZ,  M_zXYZZ,
      M_zXXYYZ,  M_zXXYZZ,  M_zXYYZZ,  M_zXXYYZZ;
   LBMReal mu200, mu020, mu002,mu110, mu101,mu011,mu210,mu120,mu102,mu111, 
      mu201,mu021,mu012,mu220,mu121,mu202,mu211,mu112,mu022,mu221,mu122,mu212,mu222,mu000,mu100,mu010,mu001;
   LBMReal vx_sq, vy_sq, vz_sq, vx_vy, vx_vz, vy_vz, vx_vy_vz;
   LBMReal MXXpMYYpMZZ,MXXmMYY, MXXmMZZ,
      MXXYpMYZZ,MXXZpMYYZ,MXYYpMXZZ,  MXXYmMYZZ,MXXZmMYYZ,MXYYmMXZZ,
      MXXYYppp,MXXYYpm2p, MXXYYppm2;
   UbTimer timer;

   LBMReal OxyyMxzz;
   int option;

   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr   zeroDistributions;

   mu::value_type muX1,muX2,muX3;
   mu::value_type muDeltaT;
   mu::value_type muNue;
   LBMReal forcingX1;
   LBMReal forcingX2;
   LBMReal forcingX3;
};

#endif
