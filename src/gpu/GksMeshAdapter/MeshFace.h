#ifndef MESH_FACE_H
#define MESH_FACE_H

#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"
#include "Core/ArrayTypes.h"

struct VIRTUALFLUIDS_GPU_EXPORT MeshFace
{
    //////////////////////////////////////////////////////////////////////////

    //      o 2                 
    //     /|                   
    //    / |                   
    //   o 3|     n            
    //   | -+--------->        
    //   |  o 1                
    //   | /                    
    //   |/                     
    //   o 0                    
    //
    //

    uint_4 faceToNode;

    //////////////////////////////////////////////////////////////////////////

    uint posCell;
    uint negCell;

    uint posCellCoarse;
    uint negCellCoarse;

    //////////////////////////////////////////////////////////////////////////

    Vec3 faceCenter;

    //////////////////////////////////////////////////////////////////////////

    char orientation;

    uint level;

    MeshFace();
};


#endif