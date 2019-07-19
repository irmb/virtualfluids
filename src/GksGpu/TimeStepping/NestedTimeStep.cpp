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

    CudaUtility::synchronizeCudaDevice();

    //////////////////////////////////////////////////////////////////////////

    FluxComputation::run( dataBase, parameters, level ); // comment out to disable commhiding
    
    //////////////////////////////////////////////////////////////////////////
    
    if( !dataBase->communicators.empty() )
    {
        //////////////////////////////////////////////////////////////////////////
        // X
        //////////////////////////////////////////////////////////////////////////
        if( dataBase->communicators[level][0] != nullptr ) dataBase->communicators[level][0]->sendData(dataBase, Communicator::tagSendNegative);
        if( dataBase->communicators[level][1] != nullptr ) dataBase->communicators[level][1]->sendData(dataBase, Communicator::tagSendPositive);

        if( dataBase->communicators[level][0] != nullptr ) dataBase->communicators[level][0]->recvData(dataBase, Communicator::tagSendPositive);
        if( dataBase->communicators[level][1] != nullptr ) dataBase->communicators[level][1]->recvData(dataBase, Communicator::tagSendNegative);
        //////////////////////////////////////////////////////////////////////////
        // Y
        //////////////////////////////////////////////////////////////////////////
        if( dataBase->communicators[level][2] != nullptr ) dataBase->communicators[level][2]->sendData(dataBase, Communicator::tagSendNegative);
        if( dataBase->communicators[level][3] != nullptr ) dataBase->communicators[level][3]->sendData(dataBase, Communicator::tagSendPositive);
        
        if( dataBase->communicators[level][2] != nullptr ) dataBase->communicators[level][2]->recvData(dataBase, Communicator::tagSendPositive);
        if( dataBase->communicators[level][3] != nullptr ) dataBase->communicators[level][3]->recvData(dataBase, Communicator::tagSendNegative);
        //////////////////////////////////////////////////////////////////////////
        // Z
        //////////////////////////////////////////////////////////////////////////
        if( dataBase->communicators[level][4] != nullptr ) dataBase->communicators[level][4]->sendData(dataBase, Communicator::tagSendNegative);
        if( dataBase->communicators[level][5] != nullptr ) dataBase->communicators[level][5]->sendData(dataBase, Communicator::tagSendPositive);
        
        if( dataBase->communicators[level][4] != nullptr ) dataBase->communicators[level][4]->recvData(dataBase, Communicator::tagSendPositive);
        if( dataBase->communicators[level][5] != nullptr ) dataBase->communicators[level][5]->recvData(dataBase, Communicator::tagSendNegative);
    }

    //////////////////////////////////////////////////////////////////////////

    FluxComputation::run( dataBase, parameters, level, true ); // comment out to disable commhiding
    
    //CudaUtility::synchronizeCudaDevice();                   // comment in to disable commhiding
    //FluxComputation::run( dataBase, parameters, level );    // comment in to disable commhiding

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

