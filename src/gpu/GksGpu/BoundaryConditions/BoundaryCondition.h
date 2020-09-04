#ifndef BoundaryCondition_H
#define BoundaryCondition_H

#include <functional>

#include <memory>
#include <vector>

#include "VirtualFluidsDefinitions.h"
#include "GksGpu_export.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "Parameters/Parameters.h"

class  GksMeshAdapter;

namespace GksGpu{

class  DataBaseAllocator;
struct DataBase;

struct BoundaryConditionStruct
{
    uint  numberOfCells;

    uint* ghostCells;
    uint* domainCells;
    uint* secondCells;
};

struct GKSGPU_EXPORT BoundaryCondition : virtual public BoundaryConditionStruct, public std::enable_shared_from_this<BoundaryCondition>
{
    SPtr<DataBaseAllocator> myAllocator;

    std::vector<uint> numberOfCellsPerLevel;
    std::vector<uint> startOfCellsPerLevel;

    BoundaryCondition( SPtr<DataBase> dataBase );

    ~BoundaryCondition();

    virtual void findBoundaryCells( GksMeshAdapter& adapter, bool allowGhostCells, std::function<bool(Vec3)> boundaryFinder);

    virtual bool isWall() = 0;

    virtual bool isFluxBC();
    
    virtual bool isInsulated();

    virtual bool secondCellsNeeded();

    virtual void runBoundaryConditionKernel( const SPtr<DataBase> dataBase,
                                             const Parameters parameters, 
                                             const uint level ) = 0;

    BoundaryConditionStruct toStruct()
    {
        BoundaryConditionStruct boundaryCondition;

        boundaryCondition.numberOfCells = this->numberOfCells;

        boundaryCondition.ghostCells      = this->ghostCells;
        boundaryCondition.domainCells     = this->domainCells;
        boundaryCondition.secondCells     = this->secondCells;

        return boundaryCondition;
    }

};

} // namespace GksGpu

#endif
