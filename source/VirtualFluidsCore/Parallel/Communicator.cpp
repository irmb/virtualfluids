#include "Communicator.h"
#include <basics/utilities/UbException.h>

SPtr<Communicator> Communicator::instance = SPtr<Communicator>();
//////////////////////////////////////////////////////////////////////////
SPtr<Communicator> Communicator::getInstance()
{
   if( !instance )
      UB_THROW(UbException(UB_EXARGS,"Communicator isn't initialized correctly! You can not create a new instance of abstract Communicator class!"));
   return instance;
}

