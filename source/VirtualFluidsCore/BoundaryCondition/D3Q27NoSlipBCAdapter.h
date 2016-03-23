#ifndef D3Q27NOSLIPBCADAPTER_H
#define D3Q27NOSLIPBCADAPTER_H

#include "D3Q27BoundaryConditionAdapter.h"

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

class D3Q27NoSlipBCAdapter : public D3Q27BoundaryConditionAdapter
{
public:
   D3Q27NoSlipBCAdapter()
    : D3Q27BoundaryConditionAdapter()
   {
   }
   D3Q27NoSlipBCAdapter(const short& secondaryBcOption)
      : D3Q27BoundaryConditionAdapter(secondaryBcOption)
   {
   }

   void init(const D3Q27Interactor* const& interactor, const double& time=0) {}
   void update(const D3Q27Interactor* const& interactor, const double& time=0) {}

   void adaptBCForDirection( const D3Q27Interactor& interactor, D3Q27BoundaryConditionPtr bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& q, const int& fdirection, const double& time=0 )
   {
      bc->setNoSlipBoundaryFlag(D3Q27System::INVDIR[fdirection],secondaryBcOption);
      bc->setQ((float)q,fdirection);
   }
   void adaptBC( const D3Q27Interactor& interactor, D3Q27BoundaryConditionPtr bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& time=0 ) 
   {

   }

private:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<D3Q27BoundaryConditionAdapter>(*this);
   }
};
#endif //D3Q27NOSLIPBCADAPTER_H
