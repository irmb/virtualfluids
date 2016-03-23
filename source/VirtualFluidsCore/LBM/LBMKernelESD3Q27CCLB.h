//Cascaded Cumulant LBM

#ifndef LBMKernelESD3Q27CCLB_H
#define LBMKernelESD3Q27CCLB_H

#include "LBMKernelETD3Q27.h"
#include "D3Q27ETBCProcessor.h"
#include "D3Q27System.h"
#include <boost/serialization/export.hpp>
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"
#include "Grid3D.h"


class LBMKernelESD3Q27CCLB;
typedef boost::shared_ptr<LBMKernelESD3Q27CCLB> LBMKernelESD3Q27CCLBPtr;


class LBMKernelESD3Q27CCLB :  public LBMKernelETD3Q27
{
public:
   LBMKernelESD3Q27CCLB();
   LBMKernelESD3Q27CCLB(int nx1, int nx2, int nx3, Grid3DPtr grid);
   ~LBMKernelESD3Q27CCLB(void);
   void calculate();
   LBMKernel3DPtr clone();
   double getCallculationTime();
   void initNeighbours();

protected:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<LBMKernelETD3Q27>(*this);
   }

   void collideAll();  

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

   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr   zeroDistributions;

   mu::value_type muX1,muX2,muX3;
   mu::value_type muDeltaT;
   mu::value_type muNue;
   LBMReal forcingX1;
   LBMReal forcingX2;
   LBMReal forcingX3;
   Grid3DPtr grid;
   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr blockN100L, blockN010L, blockN001L, blockN110L, blockN011L, blockN101L, blockN111L;
   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr blockN100NL, blockN010NL, blockN001NL, blockN110NL, blockN011NL, blockN101NL, blockN111NL;
   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr PxL, PxyL, PxyzL, PyL, PyzL, PxzL, PzL;
   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr PxNL, PxyNL, PxyzNL, PyNL, PyzNL, PxzNL, PzNL;
   int blockNX1, blockNX2, blockNX3;
};

#endif
