#ifndef DataBase_H
#define DataBase_H

#include <memory>
#include <string>
#include <vector>
#include <array>

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"
#include "Core/ArrayTypes.h"

#include "VirtualFluidsDefinitions.h"
#include "GksGpu_export.h"

#include "Definitions/AccumulatorDataType.h"

#include "CellProperties/CellProperties.cuh"

class  GksMeshAdapter;

namespace GksGpu { 

struct BoundaryCondition;
class  DataBaseAllocator;
struct DataBase;
struct PerLevelCounts;
struct DataBaseStruct;
struct Communicator;

struct GKSGPU_EXPORT DataBase : public std::enable_shared_from_this<DataBase>
{
    //////////////////////////////////////////////////////////////////////////
    // Management
    //////////////////////////////////////////////////////////////////////////

    SPtr<DataBaseAllocator> myAllocator;

    std::vector< SPtr<BoundaryCondition> > boundaryConditions;

    std::vector< std::array< SPtr< Communicator >, 6 > > communicators;

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

    std::vector<CellProperties> cellPropertiesHost;

    //////////////////////////////////////////////////////////////////////////
    // Host/Device geometry and connectivity - READ ONLY
    //////////////////////////////////////////////////////////////////////////

    uint* cellToCell;     // 6

    uint* faceToCell;     // 2

    uint* parentCell;     // 1

    real* faceCenter;     // 3
    real* cellCenter;     // 3

    CellProperties* cellProperties;     // 1 x byte

    char* faceOrientation;

    uint* fineToCoarse;   // 9
    uint* coarseToFine;   // 9

    //////////////////////////////////////////////////////////////////////////
    // Host/Device data - READ MODIFY
    //////////////////////////////////////////////////////////////////////////

    real*            data;
    realAccumulator* dataUpdate;

    real* massFlux;

    realAccumulator* diffusivity;

    int* crashCellIndex;

    //////////////////////////////////////////////////////////////////////////
    // Host only data
    //////////////////////////////////////////////////////////////////////////

    std::vector<real> dataHost;

    std::vector<real> diffusivityHost;

    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    // Methods
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    DataBase( std::string type );
    ~DataBase();

    //void setMesh( std::shared_ptr<MeshGeneratorQuadTree> mesh );

    void setMesh( GksMeshAdapter& adapter );

    void setCommunicators( GksMeshAdapter& adapter );

    void copyDataHostToDevice();

    void copyDataDeviceToHost();

    void copyDataDeviceToHost( real* dataHost );

    int getCrashCellIndex();

    DataBaseStruct toStruct();

    //////////////////////////////////////////////////////////////////////////

    uint getCellLevel( uint cellIdx );
    uint getFaceLevel( uint faceIdx );

    Vec3 getCellCenter( uint cellIdx );

    bool isGhostCell( uint cellIdx );

    std::string getDeviceType();
};

struct GKSGPU_EXPORT PerLevelCounts
{
    uint numberOfCells;
    uint startOfCells;

    uint numberOfBulkCells;

    uint numberOfFaces;

    uint numberOfInnerFaces;

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

} // namespace GksGpu

#endif