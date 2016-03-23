#ifndef LBMKernelETD3Q27CCLBWithSpongeLayer_H
#define LBMKernelETD3Q27CCLBWithSpongeLayer_H

#include "LBMKernelETD3Q27CCLB.h"
#include "D3Q27ETBCProcessor.h"
#include "D3Q27System.h"
#include <boost/serialization/export.hpp>
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

class LBMKernelETD3Q27CCLBWithSpongeLayer;
typedef boost::shared_ptr<LBMKernelETD3Q27CCLBWithSpongeLayer> LBMKernelETD3Q27CCLBWithSpongeLayerPtr;

//! \brief   Cascaded Cumulant LBM kernel. 
//! \details CFD solver with sponge layer that use Cascaded Cumulant Lattice Boltzmann method for D3Q27 model <br>
//! variable spongeFactor is depending on funktion in muSpongeLayer and varies between 1 (do nothing) und 0.5 (collFactor about 1);
//! Initialization in test case (example): <br>
//! \code{.cpp}
//! int sizeSP=8; //width of sponge layer in blocks 
//! mu::Parser spongeLayer;
//! spongeLayer.SetExpr("x1>=(sizeX-sizeSP)/dt ? (sizeX-(x1+1))/sizeSP/2.0 + 0.5 : 1.0");
//! spongeLayer.DefineConst("sizeX", nx[0]*blocknx[0]); // width of grid for X in coarse nodes
//! spongeLayer.DefineConst("sizeSP", sizeSP*blocknx[0]); // width of sponge layer in coarse nodes
//! kernel->setWithSpongeLayer(true);
//! kernel->setSpongeLayer(spongeLayer);
//! \endcode
//! \author  K. Kucher, M. Geier
class LBMKernelETD3Q27CCLBWithSpongeLayer :  public LBMKernelETD3Q27CCLB
{
public:
   LBMKernelETD3Q27CCLBWithSpongeLayer();
   //! Constructor
   //! \param nx1 number of nodes in x dimension
   //! \param nx2 number of nodes in y dimension
   //! \param nx3 number of nodes in z dimension
   //! \param p   set relaxation parameter: NORMAL is OxyyMxzz = 1.0 and MAGIC is OxyyMxzz = 2.0 +(-collFactor)
   LBMKernelETD3Q27CCLBWithSpongeLayer(int nx1, int nx2, int nx3, Parameter p);
   virtual ~LBMKernelETD3Q27CCLBWithSpongeLayer(void);
   LBMKernel3DPtr clone();
   void calculate();
protected:
   void init();
   LBMReal OxyyMxzz;

   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<LBMKernelETD3Q27CCLB>(*this);
      ar & OxyyMxzz;
   }

   void collideAll();  
};

#endif
