//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef D3Q27BOUNDARYCONDITIONADAPTER_H
#define D3Q27BOUNDARYCONDITIONADAPTER_H

#ifdef CAB_RCF
   #include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include <boost/shared_ptr.hpp>

class D3Q27BoundaryConditionAdapter;
typedef boost::shared_ptr<D3Q27BoundaryConditionAdapter> D3Q27BoundaryConditionAdapterPtr;

#include "D3Q27BoundaryCondition.h"
#include "basics/objects/ObObject.h"
#include "basics/objects/ObObjectCreator.h"
#include "basics/utilities/UbFileOutput.h"
#include "basics/utilities/UbFileInput.h"
#include "basics/utilities/UbAutoRun.hpp"


/*=========================================================================*/
/*  D3Q27BoundaryConditionAdapter                                          */
/*                                                                         */
/**
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 06.09.06
*/ 

/*
usage: ...
*/

class D3Q27Interactor;

class D3Q27BoundaryConditionAdapter
{
public:
   D3Q27BoundaryConditionAdapter() 
      :  secondaryBcOption(0)
       , type(0)
   {
   }
   D3Q27BoundaryConditionAdapter(const short& secondaryBcOption) 
      :  secondaryBcOption(secondaryBcOption) 
       , type(0)
   {
   }
   virtual ~D3Q27BoundaryConditionAdapter() {}

   //methods
   bool isTimeDependent() { return((this->type & TIMEDEPENDENT) ==  TIMEDEPENDENT); }

   virtual short getSecondaryBcOption() { return this->secondaryBcOption; }
   virtual void  setSecondaryBcOption(const short& val) { this->secondaryBcOption=val; }

   virtual void init(const D3Q27Interactor* const& interactor, const double& time=0) = 0;
   virtual void update(const D3Q27Interactor* const& interactor, const double& time=0) = 0;

   virtual void adaptBC( const D3Q27Interactor& interactor, D3Q27BoundaryConditionPtr bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& time=0 ) = 0;
   virtual void adaptBCForDirection( const D3Q27Interactor& interactor, D3Q27BoundaryConditionPtr bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& q, const int& fdirection, const double& time=0 ) = 0;

protected:
   short secondaryBcOption;

   char  type;

   static const char   TIMEDEPENDENT = 1<<0;//'1';
   static const char   TIMEPERIODIC  = 1<<1;//'2';

private:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & secondaryBcOption;
   }
};


#endif //D3Q27BOUNDARYCONDITIONADAPTER_H
