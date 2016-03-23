//Cascaded Cumulant LBM

#ifndef LBMKernelETD3Q27CCLB_H
#define LBMKernelETD3Q27CCLB_H

#include "LBMKernelETD3Q27.h"
#include "D3Q27ETBCProcessor.h"
#include "D3Q27System.h"
#include <boost/serialization/export.hpp>
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

class LBMKernelETD3Q27CCLB;
typedef boost::shared_ptr<LBMKernelETD3Q27CCLB> LBMKernelETD3Q27CCLBPtr;

//! \brief   Cascaded Cumulant LBM kernel. 
//! \details CFD solver that use Cascaded Cumulant Lattice Boltzmann method for D3Q27 model
//! \author  K. Kucher, M. Geier
class LBMKernelETD3Q27CCLB :  public LBMKernelETD3Q27
{
public:
   //! This option set relaxation parameter: NORMAL  
   enum Parameter{NORMAL, MAGIC};
public:
   LBMKernelETD3Q27CCLB();
   //! Constructor
   //! \param nx1 number of nodes in x dimension
   //! \param nx2 number of nodes in y dimension
   //! \param nx3 number of nodes in z dimension
   //! \param p   set relaxation parameter: NORMAL is OxyyMxzz = 1.0 and MAGIC is OxyyMxzz = 2.0 +(-collFactor)
   LBMKernelETD3Q27CCLB(int nx1, int nx2, int nx3, Parameter p);
   virtual ~LBMKernelETD3Q27CCLB(void);
   virtual void calculate();
   virtual LBMKernel3DPtr clone();
   double getCallculationTime();

protected:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<LBMKernelETD3Q27>(*this);
      ar & OxyyMxzz; 
      ar & parameter;
   }

   virtual void collideAll();  
   virtual void init();
   LBMReal f[D3Q27System::ENDF+1];

   UbTimer timer;

   LBMReal OxyyMxzz;
   Parameter parameter;

   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
   CbArray4D<LBMReal,IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr   zeroDistributions;

   mu::value_type muX1,muX2,muX3;
   mu::value_type muDeltaT;
   mu::value_type muNu;
   LBMReal forcingX1;
   LBMReal forcingX2;
   LBMReal forcingX3;
};

#endif
