#ifndef IncompressibleCumulantWithSpongeLayerLBMKernel_H
#define IncompressibleCumulantWithSpongeLayerLBMKernel_H

#include "IncompressibleCumulantLBMKernel.h"

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
   virtual ~IncompressibleCumulantWithSpongeLayerLBMKernel(void);
   SPtr<LBMKernel> clone();
   void calculate(int step);
   void initRelaxFactor(int vdir, double vL1, double vdx, double vSP);
   //! \param vdir where the sponge layer is placed
   //! \param vL1 length of simulation domain
   //! \param vdx subgrid space 
   //! \param vSP length of sponge layer
   void setRelaxFactorParam(int vdir, double vL1, double vdx, double vSP);
protected:
  void initDataSet();
  LBMReal OxyyMxzz;
  int direction;
  double L1;
  double dx;
  double SP;
};

#endif
