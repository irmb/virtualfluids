#ifndef DataBaseStruct_H
#define DataBaseStruct_H

#include "Core/DataTypes.h"

#include <VirtualFluidsDefinitions.h>

namespace GksGpu{ 

struct VIRTUALFLUIDS_GPU_EXPORT DataBaseStruct
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

    char* faceOrientation;

    real*            data;
    realAccumulator* dataUpdate;

    real* massFlux;

    realAccumulator* diffusivity;

    int* crashCellIndex;
};

} // namespace GksGpu


#endif