#ifndef NoSlipBCAdapter_H
#define NoSlipBCAdapter_H

#include "BCAdapter.h"

/*=========================================================================*/
/*  D3Q27NoSlipBCAdapter                                                   */
/*                                                                         */
/**
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 06.09.06
*/ 

/*
usage: ...
*/

class NoSlipBCAdapter : public BCAdapter
{
public:
   NoSlipBCAdapter()
    : BCAdapter()
   {
   }
   NoSlipBCAdapter(const short& secondaryBcOption)
      : BCAdapter(secondaryBcOption)
   {
   }

   void init(const D3Q27Interactor* const& interactor, const double& time=0) {}
   void update(const D3Q27Interactor* const& interactor, const double& time=0) {}

   void adaptBCForDirection( const D3Q27Interactor& interactor, BoundaryConditionsPtr bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& q, const int& fdirection, const double& time=0 )
   {
      bc->setNoSlipBoundaryFlag(D3Q27System::INVDIR[fdirection],secondaryBcOption);
      bc->setQ((float)q,fdirection);
   }
   void adaptBC( const D3Q27Interactor& interactor, BoundaryConditionsPtr bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& time=0 ) 
   {
      bc->setBcAlgorithmType(algorithmType);
   }

private:

};
#endif //NoSlipBCAdapter_H
