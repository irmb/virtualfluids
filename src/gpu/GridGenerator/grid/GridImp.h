#ifndef GRID_IMP_H
#define GRID_IMP_H

#include <array>

#include "Core/LbmOrGks.h"

#include "global.h"

#include "grid/distributions/Distribution.h"
#include "grid/Grid.h"
#include "grid/Cell.h"
#include "grid/Field.h" 

class TriangularMesh;
struct Vertex;
struct Triangle;
class GridStrategy;
class GridInterface;
class Object;
class BoundingBox;
class TriangularMeshDiscretizationStrategy;

#ifdef __GNUC__
    #ifndef __clang__
        #pragma push
        #pragma diag_suppress = 3156
    #endif
#endif

//GCC:  warning #3156-D: extern declaration of the entity DIRECTIONS is treated as a static definition
extern CONSTANT int DIRECTIONS[DIR_END_MAX][DIMENSION];

#ifdef __GNUC__
    #ifndef __clang__
        #pragma pop
    #endif
#endif

class GRIDGENERATOR_EXPORT GridImp : public enableSharedFromThis<GridImp>, public Grid
{
private:
    CUDA_HOST GridImp();
    CUDA_HOST GridImp(Object* object, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution d, uint level);

public:
    virtual HOSTDEVICE ~GridImp();
    static CUDA_HOST SPtr<GridImp> makeShared(Object* object, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution d, uint level);

private:
    CUDA_HOST void initalNumberOfNodesAndSize();
    HOSTDEVICE Cell getOddCellFromIndex(uint index) const;
    HOSTDEVICE bool isValidSolidStopper(uint index) const;
	HOSTDEVICE bool shouldBeBoundarySolidNode(uint index) const;
	HOSTDEVICE bool isValidEndOfGridStopper(uint index) const;
    HOSTDEVICE bool isValidEndOfGridBoundaryStopper(uint index) const;
    HOSTDEVICE bool isOutSideOfGrid(Cell &cell) const;
    HOSTDEVICE bool contains(Cell &cell, char type) const;
    HOSTDEVICE void setNodeTo(Cell &cell, char type);

    HOSTDEVICE bool nodeInPreviousCellIs(int index, char type) const;
    HOSTDEVICE bool nodeInCellIs(Cell& cell, char type) const override;

    HOSTDEVICE uint getXIndex(real x) const;
    HOSTDEVICE uint getYIndex(real y) const;
    HOSTDEVICE uint getZIndex(real z) const;

    uint level;

    real startX = 0.0, startY = 0.0, startZ = 0.0;
    real endX, endY, endZ;
    real delta = 1.0;

    bool xOddStart = false, yOddStart = false, zOddStart = false;

	uint nx, ny, nz;

	uint size;
    uint sparseSize;
    bool periodicityX = false, periodicityY = false, periodicityZ = false;

    Field field;
    Object* object;
    GridInterface* gridInterface;

    int *neighborIndexX, *neighborIndexY, *neighborIndexZ, *neighborIndexNegative;
    int *sparseIndices;

    std::vector<uint> fluidNodeIndices;
    std::vector<uint> fluidNodeIndicesBorder;

	uint *qIndices;     //maps from matrix index to qIndex
	real *qValues;
    uint *qPatches;

    bool innerRegionFromFinerGrid;

    uint numberOfLayers;

	SPtr<GridStrategy> gridStrategy;
    TriangularMeshDiscretizationStrategy* triangularMeshDiscretizationStrategy;

	uint numberOfSolidBoundaryNodes;

    bool enableFixRefinementIntoTheWall;

public:
    CUDA_HOST void inital(const SPtr<Grid> fineGrid, uint numberOfLayers) override;
    CUDA_HOST void setOddStart( bool xOddStart, bool yOddStart, bool zOddStart ) override;
    HOSTDEVICE void fixOddCell(uint index);

    CUDA_HOST void setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ) override;
    void setPeriodicityX(bool periodicity) override;
    void setPeriodicityY(bool periodicity) override;
    void setPeriodicityZ(bool periodicity) override;

    bool getPeriodicityX() override;
    bool getPeriodicityY() override;
    bool getPeriodicityZ() override;

    void setEnableFixRefinementIntoTheWall( bool enableFixRefinementIntoTheWall ) override;

    HOSTDEVICE void setCellTo(uint index, char type);
    HOSTDEVICE void setNonStopperOutOfGridCellTo(uint index, char type);

    HOSTDEVICE uint transCoordToIndex(const real &x, const real &y, const real &z) const override;
    HOSTDEVICE void transIndexToCoords(uint index, real &x, real &y, real &z) const override;

    CUDA_HOST virtual void findGridInterface(SPtr<Grid> grid, LbmOrGks lbmOrGks) override;

    HOSTDEVICE void repairGridInterfaceOnMultiGPU(SPtr<Grid> fineGrid) override;

    CUDA_HOST virtual void limitToSubDomain(SPtr<BoundingBox> subDomainBox, LbmOrGks lbmOrGks) override;

    CUDA_HOST void freeMemory() override;

    CUDA_HOST uint getLevel(real levelNull) const;
    CUDA_HOST uint getLevel() const;
    CUDA_HOST void setTriangularMeshDiscretizationStrategy(TriangularMeshDiscretizationStrategy* triangularMeshDiscretizationStrategy);
    CUDA_HOST TriangularMeshDiscretizationStrategy* getTriangularMeshDiscretizationStrategy();

	CUDA_HOST uint getNumberOfSolidBoundaryNodes() const override;
	CUDA_HOST void setNumberOfSolidBoundaryNodes(uint numberOfSolidBoundaryNodes) override;

	CUDA_HOST real getQValue(const uint index, const uint dir) const override;
	CUDA_HOST uint getQPatch(const uint index) const override;

    CUDA_HOST void setInnerRegionFromFinerGrid( bool innerRegionFromFinerGrid ) override;

    CUDA_HOST void setNumberOfLayers( uint numberOfLayers ) override;

public:
    Distribution distribution;

    HOSTDEVICE void initalNodeToOutOfGrid(uint index);

    HOSTDEVICE void findInnerNode(uint index);

    HOSTDEVICE void discretize(Object* object, char innerType, char outerType);

    bool isInside(const Cell& cell) const;

    HOSTDEVICE void setInnerBasedOnFinerGrid(const SPtr<Grid> fineGrid);
    
    HOSTDEVICE void addOverlap();
    HOSTDEVICE void setOverlapTmp( uint index );
    HOSTDEVICE void setOverlapFluid( uint index );

    HOSTDEVICE void fixRefinementIntoWall(uint xIndex, uint yIndex, uint zIndex, int dir);
    HOSTDEVICE void findStopperNode(uint index);
	HOSTDEVICE void findEndOfGridStopperNode(uint index);
	HOSTDEVICE void findSolidStopperNode(uint index);
	HOSTDEVICE void findBoundarySolidNode(uint index);

    HOSTDEVICE void findGridInterfaceCF(uint index, GridImp& finerGrid, LbmOrGks lbmOrGks);
    HOSTDEVICE void findGridInterfaceFC(uint index, GridImp& finerGrid);
    HOSTDEVICE void findOverlapStopper(uint index, GridImp& finerGrid);
    HOSTDEVICE void findInvalidBoundaryNodes(uint index);

    HOSTDEVICE void setNodeTo(uint index, char type);
    HOSTDEVICE bool isNode(uint index, char type) const;
    HOSTDEVICE bool nodeInNextCellIs(int index, char type) const;
    HOSTDEVICE bool hasAllNeighbors(uint index) const;
    HOSTDEVICE bool hasNeighborOfType(uint index, char type)const;
    HOSTDEVICE bool cellContainsOnly(Cell &cell, char type) const;
    HOSTDEVICE bool cellContainsOnly(Cell &cell, char typeA, char typeB) const;

    HOSTDEVICE const Object* getObject() const override;

    HOSTDEVICE Field getField() const;
    HOSTDEVICE char getFieldEntry(uint index) const override;
    HOSTDEVICE void setFieldEntry(uint matrixIndex, char type) override;


    HOSTDEVICE real getDelta() const override;
    HOSTDEVICE uint getSize() const override;
    HOSTDEVICE uint getSparseSize() const override;
    HOSTDEVICE int getSparseIndex(uint matrixIndex) const override;
    CUDA_HOST real* getDistribution() const override;
    CUDA_HOST int* getDirection() const override;
    CUDA_HOST int getStartDirection() const override;
    CUDA_HOST int getEndDirection() const override;

    HOSTDEVICE Vertex getMinimumOnNode(Vertex exact) const override;
    HOSTDEVICE Vertex getMaximumOnNode(Vertex exact) const override;

    HOSTDEVICE real getStartX() const override;
    HOSTDEVICE real getStartY() const override;
    HOSTDEVICE real getStartZ() const override;
    HOSTDEVICE real getEndX() const override;
    HOSTDEVICE real getEndY() const override;
    HOSTDEVICE real getEndZ() const override;
    HOSTDEVICE uint getNumberOfNodesX() const override;
    HOSTDEVICE uint getNumberOfNodesY() const override;
    HOSTDEVICE uint getNumberOfNodesZ() const override;
    CUDA_HOST void getNodeValues(real *xCoords, real *yCoords, real *zCoords, uint *neighborX, uint *neighborY, uint *neighborZ, uint *neighborNegative, uint *geo) const override;

    HOSTDEVICE uint getNumberOfNodesCF() const override;
    HOSTDEVICE uint getNumberOfNodesFC() const override;
    CUDA_HOST void getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf) const override;
    CUDA_HOST static void getGridInterface(uint* gridInterfaceList, const uint* oldGridInterfaceList, uint size);
    CUDA_HOST void getGridInterfaceIndicesFCCBorderBulk(uint *iCellFccBorder, uint &intFCBorderKfc, uint *&iCellFccBulk, uint &intFCBulkKfc, int level) const override;

    int* getNeighborsX() const override;
    int* getNeighborsY() const override;
    int* getNeighborsZ() const override;
    int* getNeighborsNegative() const override;

    CUDA_HOST uint* getCF_coarse() const override;
    CUDA_HOST uint* getCF_fine() const override;
    CUDA_HOST uint* getCF_offset() const override;


    CUDA_HOST uint* getFC_coarse() const override;
    CUDA_HOST uint* getFC_fine() const override;
    CUDA_HOST uint* getFC_offset() const override;

    SPtr<GridStrategy> getGridStrategy() const override;


    HOSTDEVICE void print() const;


public:
    CUDA_HOST virtual void findSparseIndices(SPtr<Grid> fineGrid) override;

    CUDA_HOST void updateSparseIndices();
    HOSTDEVICE void setNeighborIndices(uint index);
    HOSTDEVICE real getFirstFluidNode(real coords[3], int direction, real startCoord) const override;
    HOSTDEVICE real getLastFluidNode(real coords[3], int direction, real startCoord) const override;
private:
    HOSTDEVICE void setStopperNeighborCoords(uint index);
    HOSTDEVICE void getNeighborCoords(real &neighborX, real &neighborY, real &neighborZ, real x, real y, real z) const;
    HOSTDEVICE real getNeighborCoord(bool periodicity, real endCoord, real coords[3], int direction) const;
    HOSTDEVICE void getNegativeNeighborCoords(real &neighborX, real &neighborY, real &neighborZ, real x, real y, real z) const;
    HOSTDEVICE real getNegativeNeighborCoord(bool periodicity, real endCoord, real coords[3], int direction) const;
    

    HOSTDEVICE int getSparseIndex(const real &expectedX, const real &expectedY, const real &expectedZ) const;

    HOSTDEVICE static real getMinimumOnNodes(const real& minExact, const real& decimalStart, const real& delta);
    HOSTDEVICE static real getMaximumOnNodes(const real& maxExact, const real& decimalStart, const real& delta);

public:
    HOSTDEVICE BoundingBox getBoundingBoxOnNodes(Triangle &triangle) const;

    CUDA_HOST void mesh(Object* object) override;

    CUDA_HOST void mesh(TriangularMesh &geometry) override;
    HOSTDEVICE void mesh(Triangle &triangle);

    CUDA_HOST void closeNeedleCells() override;
    HOSTDEVICE bool closeCellIfNeedle(uint index);

    CUDA_HOST void closeNeedleCellsThinWall() override;
    HOSTDEVICE bool closeCellIfNeedleThinWall(uint index);

    CUDA_HOST void findQs(Object* object) override;
    CUDA_HOST void findQs(TriangularMesh &triangularMesh);
    HOSTDEVICE void findQs(Triangle &triangle);

    CUDA_HOST void findQsPrimitive(Object* object);
private:

    enum class qComputationStageType{
        FindSolidBoundaryNodes,
        ComputeQs
    } qComputationStage;

public:
    CUDA_HOST void enableFindSolidBoundaryNodes() override{ qComputationStage = qComputationStageType::FindSolidBoundaryNodes; }
    CUDA_HOST void enableComputeQs() override{ qComputationStage = qComputationStageType::ComputeQs; }

private:
    HOSTDEVICE void setDebugPoint(uint index, int pointValue);
	HOSTDEVICE void calculateQs(const Vertex &point, const Triangle &triangle) const;
	HOSTDEVICE void calculateQs(const uint index, const Vertex &point, const Triangle &triangle) const;
	CUDA_HOST void calculateQs(const uint index, const Vertex &point, Object* object) const;

    CUDA_HOST bool checkIfAtLeastOneValidQ(const uint index, const Vertex &point, const Triangle &triangle) const;

    CUDA_HOST bool checkIfAtLeastOneValidQ(const uint index, const Vertex &point, Object* object) const;

public:

    void findCommunicationIndices(int direction, SPtr<BoundingBox> subDomainBox, LbmOrGks lbmOrGks) override;
    void findCommunicationIndex( uint index, real coordinate, real limit, int direction );

    uint getNumberOfSendNodes(int direction) override;
    uint getNumberOfReceiveNodes(int direction) override;

    uint getSendIndex(int direction, uint index) override;
    uint getReceiveIndex(int direction, uint index) override;

    bool isSendNode(int index) const override;
    bool isReceiveNode(int index) const override;

    void repairCommunicationInices(int direction) override;

    void findFluidNodeIndices(bool splitDomain) override;
    void findFluidNodeIndicesBorder() override;

    uint getNumberOfFluidNodes() const override;
    CUDA_HOST void getFluidNodeIndices(uint *fluidNodeIndices) const override;

    uint getNumberOfFluidNodesBorder() const override;
    void getFluidNodeIndicesBorder(uint *fluidNodeIndicesBorder) const override;


public:

    struct CommunicationIndices
    {
        std::vector<uint> sendIndices;
        std::vector<uint> receiveIndices;
    };

    std::array<CommunicationIndices, 6> communicationIndices;


private:
    friend class GridGpuStrategy;
    friend class GridCpuStrategy;
};

#endif
