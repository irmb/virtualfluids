#ifndef Periodic_CUH
#define Periodic_CUH

#include <memory>

#include "VirtualFluidsDefinitions.h"
#include "GksGpu_export.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "BoundaryConditions/BoundaryCondition.h"

namespace GksGpu{

struct GKSGPU_EXPORT Periodic : public BoundaryCondition
{
    Periodic( SPtr<DataBase> dataBase ) : BoundaryCondition( dataBase ){}

    virtual bool isWall();

    virtual void findBoundaryCells( GksMeshAdapter& adapter, bool allowGhostCells, std::function<bool(Vec3)> boundaryFinder) override;

    virtual void runBoundaryConditionKernel(const SPtr<DataBase> dataBase,
                                            const Parameters parameters, 
                                            const uint level) override;

    BoundaryConditionStruct toStruct()
    {
        BoundaryConditionStruct boundaryCondition;

        boundaryCondition.numberOfCells = this->numberOfCells;

        boundaryCondition.ghostCells    = this->ghostCells;
        boundaryCondition.domainCells   = this->domainCells;
        boundaryCondition.secondCells   = this->secondCells;

        return boundaryCondition;
    }
};

} // namespace GksGpu

#endif
