#include "Initializer.h"

#include <sstream>
#define _USE_MATH_DEFINES
#include <math.h>

#include "Core/PointerDefinitions.h"
#include "Core/RealConstants.h"

#include "DataBase/DataBaseStruct.h"

#include "Definitions/MemoryAccessPattern.h"
#include "Definitions/PassiveScalar.h"

void Initializer::interpret(SPtr<DataBase> dataBase, std::function<ConservedVariables(Vec3)> initialCondition)
{
    for( uint cellIdx = 0; cellIdx < dataBase->numberOfCells; cellIdx++ ){

        Vec3 cellCenter = dataBase->getCellCenter( cellIdx );

        ConservedVariables cellCons = initialCondition(cellCenter);

        dataBase->dataHost[ RHO__(cellIdx, dataBase->numberOfCells) ] = cellCons.rho ;
        dataBase->dataHost[ RHO_U(cellIdx, dataBase->numberOfCells) ] = cellCons.rhoU;
        dataBase->dataHost[ RHO_V(cellIdx, dataBase->numberOfCells) ] = cellCons.rhoV;
        dataBase->dataHost[ RHO_W(cellIdx, dataBase->numberOfCells) ] = cellCons.rhoW;
        dataBase->dataHost[ RHO_E(cellIdx, dataBase->numberOfCells) ] = cellCons.rhoE;
    #ifdef USE_PASSIVE_SCALAR
	    dataBase->dataHost[ RHO_S_1(cellIdx, dataBase->numberOfCells) ] = cellCons.rhoS_1;
	    dataBase->dataHost[ RHO_S_2(cellIdx, dataBase->numberOfCells) ] = cellCons.rhoS_2;
    #endif // USE_PASSIVE_SCALAR
    }

    return;
}
