//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __         
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |        
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |        
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |        
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____    
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|   
//      \    \  |    |   ________________________________________________________________    
//       \    \ |    |  |  ______________________________________________________________|   
//        \    \|    |  |  |         __          __     __     __     ______      _______    
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)   
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______    
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/   
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file GridImp.h
//! \ingroup grid
//! \author Soeren Peters, Stephan Lenz, Martin Schoenherr
//=======================================================================================
#ifndef GRID_IMP_H
#define GRID_IMP_H

#include <array>
#include <vector>

#include "Core/LbmOrGks.h"

#include "gpu/GridGenerator/global.h"

#include "gpu/GridGenerator/grid/distributions/Distribution.h"
#include "gpu/GridGenerator/grid/Grid.h"
#include "gpu/GridGenerator/grid/Cell.h"
#include "gpu/GridGenerator/grid/Field.h" 

class TriangularMesh;
struct Vertex;
struct Triangle;
class GridInterface;
class Object;
class BoundingBox;
class TriangularMeshDiscretizationStrategy;


#ifdef __GNUC__
    #ifndef __clang__
        #pragma push
        #pragma nv_diag_suppress = 3156
    #endif
#endif

// GCC:  warning #3156-D: extern declaration of the entity DIRECTIONS is treated as a static definition
extern int DIRECTIONS[DIR_END_MAX][DIMENSION];

#ifdef __GNUC__
    #ifndef __clang__
        #pragma pop
    #endif
#endif

class GRIDGENERATOR_EXPORT GridImp : public enableSharedFromThis<GridImp>, public Grid
{
protected:
    GridImp() = default;
    GridImp(Object* object, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, Distribution d, uint level);

public:
    static SPtr<GridImp> makeShared(Object* object, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, std::string d3Qxx, uint level);
    virtual ~GridImp() = default;

private:
    void initalNumberOfNodesAndSize();
    Cell getOddCellFromIndex(uint index) const;
    bool isValidSolidStopper(uint index) const;
	bool shouldBeBoundarySolidNode(uint index) const;
	bool isValidEndOfGridStopper(uint index) const;
    bool isValidEndOfGridBoundaryStopper(uint index) const;
    bool isOutSideOfGrid(Cell &cell) const;
    bool contains(Cell &cell, char type) const;
    void setNodeTo(Cell &cell, char type);

    bool nodeInPreviousCellIs(int index, char type) const;
    bool nodeInCellIs(Cell& cell, char type) const override;


    uint getXIndex(real x) const;
    uint getYIndex(real y) const;
    uint getZIndex(real z) const;

    uint level;

    real startX = 0.0, startY = 0.0, startZ = 0.0;
    real endX, endY, endZ;
    real delta = 1.0;

    bool xOddStart = false, yOddStart = false, zOddStart = false;

	uint nx, ny, nz;

	uint size;
    uint sparseSize;
    bool periodicityX = false, periodicityY = false, periodicityZ = false;

    Object* object;
    GridInterface *gridInterface;

    int *sparseIndices;

    std::vector<uint> fluidNodeIndices;                 // run on CollisionTemplate::Default
    std::vector<uint> fluidNodeIndicesBorder;           // run on subdomain border nodes (CollisionTemplate::SubDomainBorder)
    std::vector<uint> fluidNodeIndicesMacroVars;        // run on CollisionTemplate::MacroVars
    std::vector<uint> fluidNodeIndicesApplyBodyForce;   // run on CollisionTemplate::ApplyBodyForce
    std::vector<uint> fluidNodeIndicesAllFeatures;      // run on CollisionTemplate::AllFeatures

	uint *qIndices;     //maps from matrix index to qIndex
	real *qValues;
    uint *qPatches;

    bool innerRegionFromFinerGrid;

    uint numberOfLayers;

    TriangularMeshDiscretizationStrategy *triangularMeshDiscretizationStrategy;

    uint numberOfSolidBoundaryNodes = 0;

    bool enableFixRefinementIntoTheWall;

    std::vector<SideType> bcAlreadySet;

protected:
    Field field;
    int *neighborIndexX, *neighborIndexY, *neighborIndexZ, *neighborIndexNegative;

public:
    void inital(const SPtr<Grid> fineGrid, uint numberOfLayers) override;
    void setOddStart(bool xOddStart, bool yOddStart, bool zOddStart) override;
    void fixOddCell(uint index);

    void setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ) override;
    void setPeriodicityX(bool periodicity) override;
    void setPeriodicityY(bool periodicity) override;
    void setPeriodicityZ(bool periodicity) override;

    bool getPeriodicityX() const override;
    bool getPeriodicityY() const override;
    bool getPeriodicityZ() const override;

    void setEnableFixRefinementIntoTheWall(bool enableFixRefinementIntoTheWall) override;

    void setCellTo(uint index, char type);
    void setNonStopperOutOfGridCellTo(uint index, char type);

    uint transCoordToIndex(const real &x, const real &y, const real &z) const override;
    void transIndexToCoords(uint index, real &x, real &y, real &z) const override;

    void findGridInterface(SPtr<Grid> grid, LbmOrGks lbmOrGks) override;

    void repairGridInterfaceOnMultiGPU(SPtr<Grid> fineGrid) override;

    void limitToSubDomain(SPtr<BoundingBox> subDomainBox, LbmOrGks lbmOrGks) override;

    void freeMemory() override;

    uint getLevel(real levelNull) const;
    uint getLevel() const;

    void setTriangularMeshDiscretizationStrategy(TriangularMeshDiscretizationStrategy *triangularMeshDiscretizationStrategy);
    TriangularMeshDiscretizationStrategy *getTriangularMeshDiscretizationStrategy();

    uint getNumberOfSolidBoundaryNodes() const override;
    void setNumberOfSolidBoundaryNodes(uint numberOfSolidBoundaryNodes) override;

    real getQValue(const uint index, const uint dir) const override;
    uint getQPatch(const uint index) const override;

    void setInnerRegionFromFinerGrid(bool innerRegionFromFinerGrid) override;

    void setNumberOfLayers(uint numberOfLayers) override;

    virtual std::vector<SideType> getBCAlreadySet();
    virtual void addBCalreadySet(SideType side);

public:
    Distribution distribution;

    void initalNodeToOutOfGrid(uint index);

    void findInnerNodes();
    void findInnerNode(uint index);

    void discretize(Object *object, char innerType, char outerType);

    bool isInside(const Cell &cell) const;

    void setInnerBasedOnFinerGrid(const SPtr<Grid> fineGrid);

    void addOverlap();
    void setOverlapTmp(uint index);
    void setOverlapFluid(uint index);

    void fixRefinementIntoWall(uint xIndex, uint yIndex, uint zIndex, int dir);
    void findStopperNode(uint index);
    void findEndOfGridStopperNode(uint index);
    void findSolidStopperNode(uint index);
    void findBoundarySolidNode(uint index);

    void findGridInterfaceCF(uint index, GridImp &finerGrid, LbmOrGks lbmOrGks);
    void findGridInterfaceFC(uint index, GridImp &finerGrid);
    void findOverlapStopper(uint index, GridImp &finerGrid);
    void findInvalidBoundaryNodes(uint index);

    void setNodeTo(uint index, char type);
    bool isNode(uint index, char type) const;
    bool nodeInNextCellIs(int index, char type) const;
    bool hasAllNeighbors(uint index) const;
    bool hasNeighborOfType(uint index, char type) const;
    bool nodeHasBC(uint index) const override;
    bool cellContainsOnly(Cell &cell, char type) const;
    bool cellContainsOnly(Cell &cell, char typeA, char typeB) const;

    const Object* getObject() const override;

    Field getField() const;
    char getFieldEntry(uint index) const override;
    void setFieldEntry(uint matrixIndex, char type) override;


    real getDelta() const override;
    uint getSize() const override;
    uint getSparseSize() const override;
    int getSparseIndex(uint matrixIndex) const override;
    real* getDistribution() const override;
    int* getDirection() const override;
    int getStartDirection() const override;
    int getEndDirection() const override;

    Vertex getMinimumOnNode(Vertex exact) const override;
    Vertex getMaximumOnNode(Vertex exact) const override;

    real getStartX() const override;
    real getStartY() const override;
    real getStartZ() const override;
    real getEndX() const override;
    real getEndY() const override;
    real getEndZ() const override;
    uint getNumberOfNodesX() const override;
    uint getNumberOfNodesY() const override;
    uint getNumberOfNodesZ() const override;
    void getNodeValues(real *xCoords, real *yCoords, real *zCoords, uint *neighborX, uint *neighborY, uint *neighborZ, uint *neighborNegative, uint *geo) const override;

    uint getNumberOfNodesCF() const override;
    uint getNumberOfNodesFC() const override;
    void getGridInterfaceIndices(uint *iCellCfc, uint *iCellCff, uint *iCellFcc, uint *iCellFcf) const override;

    static void getGridInterface(uint *gridInterfaceList, const uint *oldGridInterfaceList, uint size);

    bool isSparseIndexInFluidNodeIndicesBorder(uint &sparseIndex) const override;

    int *getNeighborsX() const override;
    int* getNeighborsY() const override;
    int* getNeighborsZ() const override;
    int* getNeighborsNegative() const override;

    uint *getCF_coarse() const override;
    uint *getCF_fine() const override;
    uint *getCF_offset() const override;

    uint *getFC_coarse() const override;
    uint *getFC_fine() const override;
    uint *getFC_offset() const override;

    void print() const;

public:
    virtual void findSparseIndices(SPtr<Grid> fineGrid) override;

    void findForGridInterfaceNewIndices(SPtr<GridImp> fineGrid);
    void updateSparseIndices();
    void setNeighborIndices(uint index);
    real getFirstFluidNode(real coords[3], int direction, real startCoord) const override;
    real getLastFluidNode(real coords[3], int direction, real startCoord) const override;
protected:
    virtual void setStopperNeighborCoords(uint index);
private:
    void getNeighborCoords(real &neighborX, real &neighborY, real &neighborZ, real x, real y, real z) const;
    real getNeighborCoord(bool periodicity, real endCoord, real coords[3], int direction) const;
    void getNegativeNeighborCoords(real &neighborX, real &neighborY, real &neighborZ, real x, real y, real z) const;
    real getNegativeNeighborCoord(bool periodicity, real endCoord, real coords[3], int direction) const;
    

    virtual int getSparseIndex(const real &expectedX, const real &expectedY, const real &expectedZ) const;

    static real getMinimumOnNodes(const real &minExact, const real &decimalStart, const real &delta);
    static real getMaximumOnNodes(const real &maxExact, const real &decimalStart, const real &delta);

public:
    BoundingBox getBoundingBoxOnNodes(Triangle &triangle) const;

    void mesh(Object *object) override;

    void mesh(TriangularMesh &geometry) override;
    void mesh(Triangle &triangle);

    void closeNeedleCells() override;
    bool closeCellIfNeedle(uint index);

    void closeNeedleCellsThinWall() override;
    bool closeCellIfNeedleThinWall(uint index);

    void findQs(Object *object) override;
    void findQs(TriangularMesh &triangularMesh);
    void findQs(Triangle &triangle);

    void findQsPrimitive(Object *object);

private:

    enum class qComputationStageType{
        FindSolidBoundaryNodes,
        ComputeQs
    } qComputationStage;

public:
    void enableFindSolidBoundaryNodes() override
    {
        qComputationStage = qComputationStageType::FindSolidBoundaryNodes;
    }
    void enableComputeQs() override { qComputationStage = qComputationStageType::ComputeQs; }

private:
    void setDebugPoint(uint index, int pointValue);
    void calculateQs(const Vertex &point, const Triangle &triangle) const;
    void calculateQs(const uint index, const Vertex &point, const Triangle &triangle) const;
    void calculateQs(const uint index, const Vertex &point, Object *object) const;

    bool checkIfAtLeastOneValidQ(const uint index, const Vertex &point, const Triangle &triangle) const;

    bool checkIfAtLeastOneValidQ(const uint index, const Vertex &point, Object *object) const;

    void allocateQs();

public:
    void findCommunicationIndices(int direction, SPtr<BoundingBox> subDomainBox, LbmOrGks lbmOrGks) override;
    void findCommunicationIndex(uint index, real coordinate, real limit, int direction);

    uint getNumberOfSendNodes(int direction) override;
    uint getNumberOfReceiveNodes(int direction) override;

    uint getSendIndex(int direction, uint index) override;
    uint getReceiveIndex(int direction, uint index) override;

    bool isSendNode(int index) const override;
    bool isReceiveNode(int index) const override;

    void repairCommunicationIndices(int direction) override;

    void findFluidNodeIndices(bool splitDomain) override;
    void findFluidNodeIndicesBorder() override;

    uint getNumberOfFluidNodes() const override;
    void getFluidNodeIndices(uint *fluidNodeIndices) const override;

    uint getNumberOfFluidNodesBorder() const override;
    void getFluidNodeIndicesBorder(uint *fluidNodeIndicesBorder) const override;

    void addFluidNodeIndicesMacroVars(std::vector<uint> _fluidNodeIndicesMacroVars) override;
    void addFluidNodeIndicesApplyBodyForce(std::vector<uint> _fluidNodeIndicesApplyBodyForce) override;
    void addFluidNodeIndicesAllFeatures(std::vector<uint> _fluidNodeIndicesAllFeatures) override;
    void sortFluidNodeIndicesMacroVars() override;
    void sortFluidNodeIndicesApplyBodyForce() override;
    void sortFluidNodeIndicesAllFeatures() override;

    uint getNumberOfFluidNodeIndicesMacroVars() const override;
    uint getNumberOfFluidNodeIndicesApplyBodyForce() const override;
    uint getNumberOfFluidNodeIndicesAllFeatures() const override; 
    void getFluidNodeIndicesMacroVars(uint *fluidNodeIndicesMacroVars) const override;
    void getFluidNodeIndicesApplyBodyForce(uint *fluidNodeIndicesApplyBodyForce) const override;
    void getFluidNodeIndicesAllFeatures(uint *fluidNodeIndicesAllFeatures) const override;

public:
    struct CommunicationIndices {
        std::vector<uint> sendIndices;
        std::vector<uint> receiveIndices;
    };

    std::array<CommunicationIndices, 6> communicationIndices;

};

#endif
