//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef SlipBCAdapter_H
#define SlipBCAdapter_H

#ifdef CAB_RCF
#include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif

#include "BCAdapter.h"

/*=======================================================*/
// D3Q27SlipBCAdapterCreator
// class D3Q27SlipBCAdapterCreator : public ObObjectCreator
//{
// public:
//   static D3Q27SlipBCAdapterCreator* getInstance()
//   {
//      static D3Q27SlipBCAdapterCreator instance;
//      return &instance;
//   }
//
//   ObObject* createObObject();
//
//   std::string getTypeID() { return "D3Q27SlipBCAdapter";}
//   std::string toString()  { return "D3Q27SlipBCAdapterCreator"; }
//
// private:
//   D3Q27SlipBCAdapterCreator( const D3Q27SlipBCAdapterCreator& );                  //no copy allowed
//   const D3Q27SlipBCAdapterCreator& operator=( const D3Q27SlipBCAdapterCreator& ); //no copy allowed
//   D3Q27SlipBCAdapterCreator() : ObObjectCreator() {}
//};
//
//#ifndef SWIG
// UB_AUTO_RUN_NAMED(
// D3Q27BCAdapterFactory::getInstance()->addObObjectCreator(D3Q27SlipBCAdapterCreator::getInstance()),
// CAB_D3Q27SlipBCAdapterCreator); #endif

/*=========================================================================*/
/*  D3Q27SlipBCAdapter                                                     */
/*                                                                         */
/**
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 06.09.06
*/

/*
usage: ...
*/

class SlipBCAdapter : public BCAdapter
{
public:
    SlipBCAdapter() : BCAdapter() {}
    SlipBCAdapter(const short &secondaryBcOption) : BCAdapter(secondaryBcOption) {}

    //------------- implements D3Q27BoundaryConditionAdapter ----- start

    void init(const D3Q27Interactor *const &interactor, const double &timestep = 0) override {}
    void update(const D3Q27Interactor *const &interactor, const double &timestep = 0) override {}

    void adaptBCForDirection(const D3Q27Interactor & /*interactor*/, SPtr<BoundaryConditions> bc,
                             const double & /*worldX1*/, const double & /*worldX2*/, const double & /*worldX3*/,
                             const double &q, const int &fdirection, const double & /*time*/ = 0) override
    {
        bc->setSlipBoundaryFlag(D3Q27System::INVDIR[fdirection], secondaryBcOption);
        bc->setQ((float)q, fdirection);
    }
    void adaptBC(const D3Q27Interactor &interactor, SPtr<BoundaryConditions> bc, const double &worldX1,
                 const double &worldX2, const double &worldX3, const double &time = 0) override;

    //------------- implements D3Q27BoundaryConditionAdapter ----- end

private:
};

#endif
