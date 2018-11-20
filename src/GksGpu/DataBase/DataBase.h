#ifndef DataBase_H
#define DataBase_H

#include <memory>
#include <string>
#include <vector>

#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"
#include "Core/ArrayTypes.h"

#include "VirtualFluidsDefinitions.h"

class  GksMeshAdapter;

class  DataBaseAllocator;
struct DataBase;
struct PerLevelCounts;
struct DataBaseStruct;

struct VF_PUBLIC DataBase : public std::enable_shared_from_this<DataBase>
{
    //////////////////////////////////////////////////////////////////////////
    // Management
    //////////////////////////////////////////////////////////////////////////

    std::shared_ptr<DataBaseAllocator> myAllocator;

    //////////////////////////////////////////////////////////////////////////
    // Sizes
    //////////////////////////////////////////////////////////////////////////

    uint numberOfNodes;

    uint numberOfCells;

    uint numberOfFaces;

    uint numberOfLevels;

    uint numberOfCoarseGhostCells;

    uint numberOfFineGhostCells;

    std::vector<PerLevelCounts> perLevelCount;

    //////////////////////////////////////////////////////////////////////////
    // Host only geometry and connectivity
    //////////////////////////////////////////////////////////////////////////

    std::vector<Vec3>   nodeCoordinates;

    std::vector<uint_8> cellToNode;
    std::vector<uint_4> faceToNode;

    //////////////////////////////////////////////////////////////////////////
    // Host/Device geometry and connectivity - READ ONLY
    //////////////////////////////////////////////////////////////////////////

    uint* cellToCell;     // 6

    uint* faceToCell;     // 2

    uint* parentCell;     // 1

    real* faceCenter;     // 3
    real* cellCenter;     // 3

    bool* faceIsWall;     // 1

    uint* fineToCoarse;   // 9
    uint* coarseToFine;   // 15

    //////////////////////////////////////////////////////////////////////////
    // Host/Device data - READ MODIFY
    //////////////////////////////////////////////////////////////////////////

    real* data;
    real* dataUpdate;

    real* massFlux;

    //////////////////////////////////////////////////////////////////////////
    // Host only data
    //////////////////////////////////////////////////////////////////////////

    std::vector<real> dataHost;

    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    // Methods
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    DataBase( std::string type );
    ~DataBase();

    //void setMesh( std::shared_ptr<MeshGeneratorQuadTree> mesh );

    void setMesh( GksMeshAdapter& adapter );

    void copyDataHostToDevice();

    void copyDataDeviceToHost();

    void copyDataDeviceToHost( real* dataHost );

    DataBaseStruct toStruct();

    //////////////////////////////////////////////////////////////////////////

    uint getCellLevel( uint cellIdx );
    uint getFaceLevel( uint faceIdx );

    bool isGhostCell( uint cellIdx );

    std::string getDeviceType();
};

struct PerLevelCounts
{
    uint numberOfCells;
    uint startOfCells;

    uint numberOfBulkCells;

    uint numberOfFaces;

    uint numberOfFacesX;
    uint startOfFacesX;

    uint numberOfFacesY;
    uint startOfFacesY;

    uint numberOfFacesZ;
    uint startOfFacesZ;

    uint numberOfCoarseToFine;
    uint startOfCoarseToFine;

    uint numberOfFineToCoarse;
    uint startOfFineToCoarse;
};

#endif