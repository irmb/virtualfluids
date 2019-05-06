#ifndef AdiabaticWall_CUH
#define AdiabaticWall_CUH

#include <memory>

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "BoundaryConditions/BoundaryCondition.h"

struct AdiabaticWallStruct
{
    uint  numberOfCells;

    uint* ghostCells;
    uint* domainCells;
    uint* secondCells;

    Vec3 velocity;

    bool useSecondCells;
};

struct VF_PUBLIC AdiabaticWall : public BoundaryCondition //, public IsothermalWallStruct
{
    Vec3 velocity;

    bool useSecondCells;

    AdiabaticWall( SPtr<DataBase> dataBase, Vec3 velocity, bool useSecondCells );

    virtual bool isWall() override;

    virtual bool isInsulated() override;

    virtual bool secondCellsNeeded() override;

    virtual void runBoundaryConditionKernel(const SPtr<DataBase> dataBase,
                                            const Parameters parameters, 
                                            const uint level) override;

    AdiabaticWallStruct toStruct()
    {
        AdiabaticWallStruct boundaryCondition;

        boundaryCondition.numberOfCells = this->numberOfCells;

        boundaryCondition.ghostCells    = this->ghostCells;
        boundaryCondition.domainCells   = this->domainCells;
        boundaryCondition.secondCells   = this->secondCells;

        boundaryCondition.velocity      = this->velocity;

        boundaryCondition.useSecondCells  = this->useSecondCells;

        return boundaryCondition;
    }
};

#endif
