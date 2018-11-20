#include "NestedTimeStep.h"

#include "Core/RealConstants.h"

#include "CellUpdate/CellUpdate.h"

void TimeStepping::nestedTimeStep( SPtr<DataBase> dataBase, 
                                   Parameters parameters,
                                   uint level )
{
    if( level != 0 ) parameters.dt /= two;
    if( level != 0 ) parameters.dx /= two;

    //if( level != dataBase->numberOfLevels - 1 ){
    //    runFineToCoarseKernel( dataBase, level ); getLastCudaError();
    //}

    //for( std::shared_ptr<BoundaryCondition> bc : bcList ){
    //    bc->runBoundaryConditionKernel( dataBase->toStruct(), parameters, type, level ); getLastCudaError();
    //}

    if( level != dataBase->numberOfLevels - 1 ){

        //runCoarseToFineKernel( dataBase, type, level ); getLastCudaError();

        nestedTimeStep( dataBase, parameters, level + 1 );
        nestedTimeStep( dataBase, parameters, level + 1 );
    }

    //runFluxKernel( dataBase, parameters, type, level ); getLastCudaError();

    CellUpdate::updateCells( dataBase, parameters, level );
}

