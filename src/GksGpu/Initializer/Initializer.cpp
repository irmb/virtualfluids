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

        Vec3 cellCenter = getCellCenter(dataBase, cellIdx);

        ConservedVariables cellCons = initialCondition(cellCenter);

        dataBase->dataHost[ RHO__(cellIdx, dataBase->numberOfCells) ] = cellCons.rho ;
        dataBase->dataHost[ RHO_U(cellIdx, dataBase->numberOfCells) ] = cellCons.rhoU;
        dataBase->dataHost[ RHO_V(cellIdx, dataBase->numberOfCells) ] = cellCons.rhoV;
        dataBase->dataHost[ RHO_E(cellIdx, dataBase->numberOfCells) ] = cellCons.rhoE;
    #ifdef USE_PASSIVE_SCALAR
	    dataBase->dataHost[ RHO_S(cellIdx, dataBase->numberOfCells) ] = cellCons.rhoS;
    #endif // USE_PASSIVE_SCALAR
    }

    return;
}

Vec3 Initializer::getCellCenter(SPtr<DataBase> dataBase, uint cellIdx)
{
    Vec3 cellCenter;

    for( uint node = 0; node < 8; node++ )
    {
        cellCenter = cellCenter + dataBase->nodeCoordinates[ dataBase->cellToNode[ cellIdx ][ node ] ];
    }

    cellCenter.x /= eight;
    cellCenter.y /= eight;
    cellCenter.z /= eight;

    return cellCenter;
}
