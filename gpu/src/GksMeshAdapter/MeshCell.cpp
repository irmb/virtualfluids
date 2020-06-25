#include "MeshCell.h"

#include "GridGenerator/grid/NodeValues.h"

MeshCell::MeshCell(){

    level   = INVALID_INDEX;
    gridIdx = INVALID_INDEX;


    for( uint& index : this->cellToNode     ) index = INVALID_INDEX;
    for( uint& index : this->cellToEdgeNode ) index = INVALID_INDEX;
    for( uint& index : this->cellToFaceNode ) index = INVALID_INDEX;
    for( uint& index : this->cellToCell     ) index = INVALID_INDEX;
    for( uint& index : this->children       ) index = INVALID_INDEX;
    
    parent = INVALID_INDEX;

    for( bool& flag : this->faceExists    ) flag = false;

    isGhostCell = false;

    isWall = false;

    isFluxBC = false;

    isInsulated = false;

    isRecvCell = false;
}

bool MeshCell::isCoarseGhostCell()
{
    return this->type == FLUID_FCC;
}

bool MeshCell::isFineGhostCell()
{
    return this->type == FLUID_CFF;
}
