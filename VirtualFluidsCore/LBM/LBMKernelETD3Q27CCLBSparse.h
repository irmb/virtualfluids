//Cascaded Cumulant LBM

#ifndef LBMKernelETD3Q27CCLBSparse_H
#define LBMKernelETD3Q27CCLBSparse_H

#include "LBMKernelETD3Q27.h"
#include "D3Q27ETBCProcessor.h"
#include "D3Q27System.h"
#include <boost/serialization/export.hpp>
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"
//#include <SparseMatrix3D.h>
//#include <SparseMatrix4D.h>


class LBMKernelETD3Q27CCLBSparse;
typedef boost::shared_ptr<LBMKernelETD3Q27CCLBSparse> LBMKernelETD3Q27CCLBSparsePtr;

//! \brief   Cascaded Cumulant LBM kernel. 
//! \details CFD solver that use Cascaded Cumulant Lattice Boltzmann method for D3Q27 model
//! \author  K. Kucher, M. Geier
class LBMKernelETD3Q27CCLBSparse :  public LBMKernelETD3Q27
{
public:
   //! This option set relaxation parameter: NORMAL  
   enum Parameter{NORMAL, MAGIC};
public:
   LBMKernelETD3Q27CCLBSparse();
   //! Constructor
   //! \param nx1 number of nodes in x dimension
   //! \param nx2 number of nodes in y dimension
   //! \param nx3 number of nodes in z dimension
   //! \param p   set relaxation parameter: NORMAL is OxyyMxzz = 1.0 and MAGIC is OxyyMxzz = 2.0 +(-collFactor)
   LBMKernelETD3Q27CCLBSparse(int nx1, int nx2, int nx3, Parameter p);
   virtual ~LBMKernelETD3Q27CCLBSparse(void);
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
      ar & parameter;
   }

   void collideAll();  
   void init();
   LBMReal f[D3Q27System::ENDF+1];

   UbTimer timer;

   LBMReal OxyyMxzz;
   Parameter parameter;
   
   mu::value_type muX1,muX2,muX3;
   mu::value_type muDeltaT;
   mu::value_type muNu;
   LBMReal forcingX1;
   LBMReal forcingX2;
   LBMReal forcingX3;
};

#endif
