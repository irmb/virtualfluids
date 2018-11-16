#ifndef GKS_MESH_ADAPTER_H
#define GKS_MESH_ADAPTER_H

#include <memory>
#include <array>
#include <vector>

#include "Core/DataTypes.h"
#include "Core/PointerDefinitions.h"

#include "MeshCell.h"
#include "MeshFace.h"


class MultipleGridBuilder;

class VF_PUBLIC GksMeshAdapter{

public:

    uint numberOfLevels;

    double dxCoarse;

    //////////////////////////////////////////////////////////////////////////

    std::vector<Vec3> nodes;

    //////////////////////////////////////////////////////////////////////////
    //
    //    C e l l    s o r t i n g :
    //
    //  | Level 0                   | Level 1                   | Level 2                   |
    //  | FluidCells   | GhostCells | FluidCells   | GhostCells | FluidCells   | GhostCells | 
    //
    //  GhostCells: not included in Cell update, i.e. BoundaryCells and FCC-Cells
    //  FluidCells: all other, including CFF-Cells
    //
    std::vector<MeshCell> cells;

    std::vector<uint> numberOfCellsPerLevel;
    std::vector<uint> numberOfBulkCellsPerLevel;
    std::vector<uint> startOfCellsPerLevel;

    //////////////////////////////////////////////////////////////////////////
    //
    //    F a c e    s o r t i n g :
    //
    //  | Level 0                              | Level 1                              | Level 2                              |
    //  | x-normal   | y-normal   | z-normal   | x-normal   | y-normal   | z-normal   | x-normal   | y-normal   | z-normal   |
    //
    std::vector<MeshFace> faces;

    std::vector<uint> numberOfFacesPerLevelXYZ;
    std::vector<uint> startOfFacesPerLevelXYZ;

    //std::vector<uint> solidFaces;

    //////////////////////////////////////////////////////////////////////////
    //
    //    F i n e T o C o a r s e    s o r t i n g :
    //
    //  | Coarse Cell Idx | Child Idcs | ...
    //
    std::vector<uint_9> fineToCoarse;

    std::vector<uint> numberOfFineToCoarsePerLevel;
    std::vector<uint> startOfFineToCoarsePerLevel;

    //////////////////////////////////////////////////////////////////////////
    //
    //    F i n e T o C o a r s e    s o r t i n g :
    //
    //  | Coarse Cell Idx | Coarse Neighbor Idcs | Child Idcs | ...
    //
    std::vector<uint_15> coarseToFine;

    std::vector<uint> numberOfCoarseToFinePerLevel;
    std::vector<uint> startOfCoarseToFinePerLevel;

    //////////////////////////////////////////////////////////////////////////
    // 
    // Connectivity from LBM grid to GKS Mesh
    //
    //    cellIdx = gridToMesh[ level ][ gridIdx ]
    //
    std::vector< std::vector<uint> > gridToMesh;

    //////////////////////////////////////////////////////////////////////////

    uint_8 cornerCells;

    //////////////////////////////////////////////////////////////////////////

public:

    void inputGrid( SPtr<MultipleGridBuilder> gridBuilder );

    void findQuadtreeConnectivity( SPtr<MultipleGridBuilder> gridBuilder );

    void findCellToCellConnectivity( SPtr<MultipleGridBuilder> gridBuilder );

    void countCells();

    void partitionCells();

    void refreshCellConnectivity(const std::vector<uint>& idxMap);

    void findCornerCells( SPtr<MultipleGridBuilder> gridBuilder, real z0 );

    void generateNodes( SPtr<MultipleGridBuilder> gridBuilder );

    void inputQs( SPtr<MultipleGridBuilder> gridBuilder );

    void morphNodes();

    void repairTriangularCells();

    void computeCellGeometry();

    void generateSolidGhostCells();

    void computeGhostCellGeometry();

    void generateFaces( SPtr<MultipleGridBuilder> gridBuilder );

    void generateSolidFaces( SPtr<MultipleGridBuilder> gridBuilder );

    void sortFaces();

    void countFaces();

    void findSolidFaces();

    void generateInterfaceConnectivity();

    void writeMeshVTK( std::string filename );

    void writeMeshFaceVTK( std::string filename );

    void writeMeshCellToCellVTK( std::string filename );

    void writeMeshFaceToCellVTK( std::string filename );

    //////////////////////////////////////////////////////////////////////////

    double getDx(uint level);
};


#endif