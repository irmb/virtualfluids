#include "NestedTimeStep.h"

#include "Core/RealConstants.h"

#include "BoundaryConditions/BoundaryCondition.h"

#include "CellUpdate/CellUpdate.h"

#include "FluxComputation/FluxComputation.h"

#include "Interface/Interface.h"

void TimeStepping::nestedTimeStep( SPtr<DataBase> dataBase, 
                                   Parameters parameters,
                                   SPtr<Communicator> communicator,
                                   uint level )
{
    if( level != 0 ) parameters.dt /= two;
    if( level != 0 ) parameters.dx /= two;

    if( level != dataBase->numberOfLevels - 1 ){
        Interface::runFineToCoarse( dataBase, level );
    }

    for( SPtr<BoundaryCondition> bc : dataBase->boundaryConditions ){
        bc->runBoundaryConditionKernel( dataBase, parameters, level );
    }

    if( level == 0 && communicator != nullptr){
        communicator->exchangeData( dataBase );
    }

    if( level != dataBase->numberOfLevels - 1 ){
    
        Interface::runCoarseToFine( dataBase, level );

        nestedTimeStep( dataBase, parameters, communicator, level + 1 );
        nestedTimeStep( dataBase, parameters, communicator, level + 1 );
    }

    FluxComputation::run( dataBase, parameters, level );

    CellUpdate::run( dataBase, parameters, level );
}

