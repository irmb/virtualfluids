#ifndef IncompressibleCumulantWithSpongeLayerLBMKernel_H
#define IncompressibleCumulantWithSpongeLayerLBMKernel_H

#include "IncompressibleCumulantLBMKernel.h"
#include "BCProcessor.h"
#include "D3Q27System.h"
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

class IncompressibleCumulantWithSpongeLayerLBMKernel;
typedef SPtr<IncompressibleCumulantWithSpongeLayerLBMKernel> LBMKernelETD3Q27CCLBWithSpongeLayerPtr;

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
//! \author  K. Kucher, M. Geier, A. Karanchuk
class IncompressibleCumulantWithSpongeLayerLBMKernel :  public IncompressibleCumulantLBMKernel
{
public:
   IncompressibleCumulantWithSpongeLayerLBMKernel();
   //! Constructor
   //! \param nx1 number of nodes in x dimension
   //! \param nx2 number of nodes in y dimension
   //! \param nx3 number of nodes in z dimension
   //! \param p   set relaxation parameter: NORMAL is OxyyMxzz = 1.0 and MAGIC is OxyyMxzz = 2.0 +(-collFactor)
   IncompressibleCumulantWithSpongeLayerLBMKernel(int nx1, int nx2, int nx3, Parameter p);
   virtual ~IncompressibleCumulantWithSpongeLayerLBMKernel(void);
   SPtr<LBMKernel> clone();
   void calculate();
   void initRelaxFactor(int vdir, double vL1, double vdx, double vSP);
   //! \param vdir where the sponge layer is placed
   //! \param vL1 length of simulation domain
   //! \param vdx subgrid space 
   //! \param vSP length of sponge layer
   void setRelaxFactorParam(int vdir, double vL1, double vdx, double vSP);
protected:
   void init();
  LBMReal OxyyMxzz;
  int direction;
  double L1;
  double dx;
  double SP;

   void collideAll();  
};

#endif
