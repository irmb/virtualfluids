#include "NestedTimeStep.h"

#include <iostream>

#include "Core/RealConstants.h"

#include "BoundaryConditions/BoundaryCondition.h"
#include "Communication/Communicator.h"
#include "CellUpdate/CellUpdate.h"
#include "FluxComputation/FluxComputation.h"
#include "Interface/Interface.h"

void TimeStepping::nestedTimeStep( SPtr<DataBase> dataBase, 
                                   Parameters parameters,
                                   uint level )
{
    //////////////////////////////////////////////////////////////////////////

    if( level != 0 ) parameters.dt /= two;
    if( level != 0 ) parameters.dx /= two;

    //////////////////////////////////////////////////////////////////////////

    if( level != dataBase->numberOfLevels - 1 )
    {
        Interface::runFineToCoarse( dataBase, level );
    }

    //////////////////////////////////////////////////////////////////////////

    for( SPtr<BoundaryCondition> bc : dataBase->boundaryConditions )
    {
        bc->runBoundaryConditionKernel( dataBase, parameters, level );
    }

    //////////////////////////////////////////////////////////////////////////
    
    if( !dataBase->communicators.empty() )
    {
        //for( uint direction = 0; direction < 6; direction++ )
        //{
        //    if( dataBase->communicators[level][direction] != nullptr )
        //    {
        //        dataBase->communicators[level][direction]->exchangeData(dataBase);
        //    }
        //}
        for( uint direction = 0; direction < 6; direction++ )
        {
            if( dataBase->communicators[level][direction] != nullptr )
            {
                dataBase->communicators[level][direction]->sendData(dataBase);
            }
        }
        for( uint direction = 0; direction < 6; direction++ )
        {
            if( dataBase->communicators[level][direction] != nullptr )
            {
                dataBase->communicators[level][direction]->recvData(dataBase);
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////

    if( level != dataBase->numberOfLevels - 1 )
    {
        Interface::runCoarseToFine( dataBase, level );

        nestedTimeStep( dataBase, parameters, level + 1 );
        nestedTimeStep( dataBase, parameters, level + 1 );
    }

    //////////////////////////////////////////////////////////////////////////

    FluxComputation::run( dataBase, parameters, level );

    //////////////////////////////////////////////////////////////////////////

    CellUpdate::run( dataBase, parameters, level );

    //////////////////////////////////////////////////////////////////////////
}

