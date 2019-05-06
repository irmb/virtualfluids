#ifndef MESH_CELL_H
#define MESH_CELL_H

#include <array>

#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"
#include "Core/ArrayTypes.h"

struct VF_PUBLIC MeshCell{

    uint level;
    uint gridIdx;

    //////////////////////////////////////////////////////////////////////////

    uint_8 cellToNode;           uint_12 cellToEdgeNode;         uint_6 cellToFaceNode;

    // for sorting see LBM f numbering
    //                                   8                              
    //     4 o--------o 0            o--------o                    o--------o 
    //      /|       /|            7/|      4/|                   /|  4    /|  
    //     / |      / |            /3| 11   / |                  / |    2 / |  
    //  6 o--------o 2|           o--------o  |0                o--------o  |  
    //    |  |     |  |           |  |  10 |  |                 |1 |     | 0|  
    //    |5 o-----|--o 1        1|  o-----|--o                 |  o-----|--o  
    //    | /      | /            |5/     2| /6                 | / 3    | /   
    //    |/       |/             |/   9   |/                   |/    5  |/    
    //  7 o--------o 3            o--------o                    o--------o      
    //
    //  z | / y
    //    |/
    //    +---- x
    //     

    //////////////////////////////////////////////////////////////////////////

    uint_27 cellToCell;

    // for sorting see LBM f numbering

    //////////////////////////////////////////////////////////////////////////

    uint_8 children;

    // for sorting see cellToNode

    uint  parent;

    //////////////////////////////////////////////////////////////////////////

    Vec3   cellCenter;

    //////////////////////////////////////////////////////////////////////////

    // order: +x, -x, +y, -y, +z, -z (see cellToCell)
    bool_6 faceExists;

    bool isGhostCell;

    bool isWall;

    bool isFluxBC;

    bool isInsulated;

    char type;

    MeshCell();

    bool isCoarseGhostCell();

    bool isFineGhostCell();
};


#endif