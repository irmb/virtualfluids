//Cascaded Cumulant LBM

#ifndef LBMKernelETD3Q27CCLBEX2_H
#define LBMKernelETD3Q27CCLBEX2_H

#include "LBMKernelETD3Q27.h"
#include "D3Q27ETBCProcessor.h"
#include "D3Q27System.h"
#include <boost/serialization/export.hpp>
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"
#include "Grid3D.h"
#include "D3Q27EsoTwist3DSoA.h"

class LBMKernelETD3Q27CCLBex2;
typedef boost::shared_ptr<LBMKernelETD3Q27CCLBex2> LBMKernelETD3Q27CCLBex2Ptr;


class LBMKernelETD3Q27CCLBex2 :  public LBMKernelETD3Q27
{
public:
   LBMKernelETD3Q27CCLBex2();
   LBMKernelETD3Q27CCLBex2(int nx1, int nx2, int nx3, int option, Grid3DPtr grid);
   ~LBMKernelETD3Q27CCLBex2(void);
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
   }

   void collideAll();  
   void init();
   void setNewDistributions(int x1, int x2, int x3);
   LBMReal f[D3Q27System::ENDF+1];
   LBMReal rho, vx, vy, vz, durchrho;
   LBMReal w1,w2,w3,w4,w5,w6,w7,w8,w9,w10;
   LBMReal s9,c1o27,c2o3;

   LBMReal mu200, mu020, mu002,mu110, mu101,mu011,mu210,mu120,mu102,mu111, 
      mu201,mu021,mu012,mu220,mu121,mu202,mu211,mu112,mu022,mu221,mu122,mu212,mu222,mu000,mu100,mu010,mu001;
   LBMReal vx_sq, vy_sq, vz_sq, vx_vy, vx_vz, vy_vz, vx_vy_vz;

   UbTimer timer;

   LBMReal OxyyMxzz;
   int option;

   Distributions d;

   mu::value_type muX1,muX2,muX3;
   mu::value_type muDeltaT;
   mu::value_type muNue;
   LBMReal forcingX1;
   LBMReal forcingX2;
   LBMReal forcingX3;
   Grid3DPtr grid;
   int lX1, lX2, lX3;
   int bX1;
   int bX2;
   int bX3;
   int level;
   UbTupleInt3 blockNX;
};

#endif
