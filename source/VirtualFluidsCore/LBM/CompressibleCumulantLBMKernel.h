#ifndef CompressibleCumulantLBMKernel_h__
#define CompressibleCumulantLBMKernel_h__

#include "LBMKernel.h"
#include "BCProcessor.h"
#include "D3Q27System.h"
#include <boost/serialization/export.hpp>
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

class CompressibleCumulantLBMKernel;
typedef boost::shared_ptr<CompressibleCumulantLBMKernel> CompressibleCumulantLBMKernelPtr;

//! \brief   compressible cumulant LBM kernel. 
//! \details CFD solver that use Cascaded Cumulant Lattice Boltzmann method for D3Q27 model
//! \author  K. Kutscher, M. Geier
class CompressibleCumulantLBMKernel :  public LBMKernel
{
public:
   //! This option set relaxation parameter: NORMAL  
   enum Parameter{NORMAL, MAGIC};
public:
   CompressibleCumulantLBMKernel();
   //! Constructor
   //! \param nx1 number of nodes in x dimension
   //! \param nx2 number of nodes in y dimension
   //! \param nx3 number of nodes in z dimension
   //! \param p   set relaxation parameter: NORMAL is OxyyMxzz = 1.0 and MAGIC is OxyyMxzz = 2.0 +(-collFactor)
   CompressibleCumulantLBMKernel(int nx1, int nx2, int nx3, Parameter p);
   virtual ~CompressibleCumulantLBMKernel(void);
   virtual void calculate();
   virtual LBMKernelPtr clone();
   double getCallculationTime();

protected:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<LBMKernel>(*this);
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
#endif // CompressibleCumulantLBMKernel_h__


