#ifndef BoundaryCondition_H
#define BoundaryCondition_H

#include <functional>

#include <memory>
#include <vector>

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "Parameters/Parameters.h"

class  GksMeshAdapter;
class  DataBaseAllocator;
struct DataBase;
struct BoundaryConditionStruct;

struct BoundaryConditionStruct
{
    uint  numberOfCells;

    uint* ghostCells;
    uint* domainCells;
    uint* secondCells;
};

struct VF_PUBLIC BoundaryCondition : virtual public BoundaryConditionStruct, public std::enable_shared_from_this<BoundaryCondition>
{
    SPtr<DataBaseAllocator> myAllocator;

    std::vector<uint> numberOfCellsPerLevel;
    std::vector<uint> startOfCellsPerLevel;

    BoundaryCondition( SPtr<DataBase> dataBase );

    ~BoundaryCondition();

    virtual void findBoundaryCells( GksMeshAdapter& adapter, bool allowGhostCells, std::function<bool(Vec3)> boundaryFinder);

    virtual bool isWall() = 0;

    virtual bool isFluxBC();

    virtual bool secondCellsNeeded();

    virtual void runBoundaryConditionKernel( const SPtr<DataBase> dataBase,
                                             const Parameters parameters, 
                                             const uint level ) = 0;

};

#endif
