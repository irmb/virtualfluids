#include "NestedTimeStep.h"

#include <iostream>

#include "Core/RealConstants.h"

#include "BoundaryConditions/BoundaryCondition.h"
#include "Communication/Communicator.h"
#include "CellUpdate/CellUpdate.h"
#include "FluxComputation/FluxComputation.h"
#include "Interface/Interface.h"
#include "Initializer/Initializer.h"
#include "CudaUtility/CudaUtility.h"

void TimeStepping::nestedTimeStep( SPtr<DataBase> dataBase, 
                                   Parameters parameters,
                                   uint level )
{
    //////////////////////////////////////////////////////////////////////////

    if( level != 0 ) parameters.dt /= c2o1;
    if( level != 0 ) parameters.dx /= c2o1;

    //////////////////////////////////////////////////////////////////////////

    //set different viscosity on specific levels
    //if( level >= 3 ) parameters.mu = 1.0e-3;

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

    FluxComputation::run( dataBase, parameters, level );
    
    //////////////////////////////////////////////////////////////////////////
    
    if( !dataBase->communicators.empty() )
    {
        //////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////
        //
        // TODO: Check if this actually works in 3D!!!
        //
        //////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////
        //for( uint direction = 0; direction < 6; direction++ )
        //{
        //    if( dataBase->communicators[level][direction] != nullptr )
        //    {
        //        dataBase->communicators[level][direction]->exchangeData(dataBase);
        //    }
        //}
        //////////////////////////////////////////////////////////////////////////
        //for( uint direction = 0; direction < 6; direction++ )
        //{
        //    if( dataBase->communicators[level][direction] != nullptr )
        //    {
        //        dataBase->communicators[level][direction]->sendData(dataBase);
        //    }
        //}
        //for( uint direction = 0; direction < 6; direction++ )
        //{
        //    if( dataBase->communicators[level][direction] != nullptr )
        //    {
        //        dataBase->communicators[level][direction]->recvData(dataBase);
        //    }
        //}
        //////////////////////////////////////////////////////////////////////////
        // X
        //////////////////////////////////////////////////////////////////////////
        if( dataBase->communicators[level][0] != nullptr ) dataBase->communicators[level][0]->sendData(dataBase);
        if( dataBase->communicators[level][1] != nullptr ) dataBase->communicators[level][1]->sendData(dataBase);

        if( dataBase->communicators[level][0] != nullptr ) dataBase->communicators[level][0]->recvData(dataBase);
        if( dataBase->communicators[level][1] != nullptr ) dataBase->communicators[level][1]->recvData(dataBase);
        //////////////////////////////////////////////////////////////////////////
        // Y
        //////////////////////////////////////////////////////////////////////////
        if( dataBase->communicators[level][2] != nullptr ) dataBase->communicators[level][2]->sendData(dataBase);
        if( dataBase->communicators[level][3] != nullptr ) dataBase->communicators[level][3]->sendData(dataBase);
        
        if( dataBase->communicators[level][2] != nullptr ) dataBase->communicators[level][2]->recvData(dataBase);
        if( dataBase->communicators[level][3] != nullptr ) dataBase->communicators[level][3]->recvData(dataBase);
        //////////////////////////////////////////////////////////////////////////
        // Z
        //////////////////////////////////////////////////////////////////////////
        if( dataBase->communicators[level][4] != nullptr ) dataBase->communicators[level][4]->sendData(dataBase);
        if( dataBase->communicators[level][5] != nullptr ) dataBase->communicators[level][5]->sendData(dataBase);
        
        if( dataBase->communicators[level][4] != nullptr ) dataBase->communicators[level][4]->recvData(dataBase);
        if( dataBase->communicators[level][5] != nullptr ) dataBase->communicators[level][5]->recvData(dataBase);
    }

    //////////////////////////////////////////////////////////////////////////

    FluxComputation::run( dataBase, parameters, level, true );

    //////////////////////////////////////////////////////////////////////////

    CudaUtility::synchronizeCudaDevice();

    //////////////////////////////////////////////////////////////////////////

    if( level != dataBase->numberOfLevels - 1 )
    {
        Interface::runCoarseToFine( dataBase, level );

        nestedTimeStep( dataBase, parameters, level + 1 );
        nestedTimeStep( dataBase, parameters, level + 1 );
    }

    //////////////////////////////////////////////////////////////////////////

    CellUpdate::run( dataBase, parameters, level );

    //////////////////////////////////////////////////////////////////////////
}

