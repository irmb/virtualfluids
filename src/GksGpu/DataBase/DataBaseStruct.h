#ifndef DataBaseStruct_H
#define DataBaseStruct_H

#include "Core/DataTypes.h"

#include <VirtualFluidsDefinitions.h>

struct VF_PUBLIC DataBaseStruct
{
    uint  numberOfCells;
    uint  numberOfFaces;
    
    uint  numberOfCoarseGhostCells;
    uint  numberOfFineGhostCells;

    uint* cellToCell;
    uint* faceToCell;

    uint* parentCell;

    uint* fineToCoarse;
    uint* coarseToFine;

    real* faceCenter;
    real* cellCenter;

    CellProperties* cellProperties;

    real* data;
    real* dataUpdate;

    real* massFlux;
};


#endif