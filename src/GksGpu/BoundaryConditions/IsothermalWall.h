#ifndef IsothermalWall_CUH
#define IsothermalWall_CUH

#include <memory>

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "BoundaryConditions/BoundaryCondition.h"

//struct IsothermalWallStruct : virtual public BoundaryConditionStruct
//{
//    Vec3 velocity;
//    real lambda;
//    real S;
//};

struct IsothermalWallStruct
{
    uint  numberOfCells;

    uint* ghostCells;
    uint* domainCells;
    uint* secondCells;

    Vec3 velocity;
    real lambda;
    real S;
};

struct VF_PUBLIC IsothermalWall : public BoundaryCondition //, public IsothermalWallStruct
{
    Vec3 velocity;
    real lambda;
    real S;

    IsothermalWall( SPtr<DataBase> dataBase, Vec3 velocity, real lambda, real S );

    virtual bool isWall() override;

    virtual bool secondCellsNeeded() override;

    virtual void runBoundaryConditionKernel(const SPtr<DataBase> dataBase,
                                            const Parameters parameters, 
                                            const uint level) override;

    IsothermalWallStruct toStruct()
    {
        IsothermalWallStruct boundaryCondition;

        boundaryCondition.numberOfCells = this->numberOfCells;

        boundaryCondition.ghostCells    = this->ghostCells;
        boundaryCondition.domainCells   = this->domainCells;
        boundaryCondition.secondCells   = this->secondCells;

        boundaryCondition.velocity      = this->velocity;
        boundaryCondition.lambda        = this->lambda;
        boundaryCondition.S             = this->S;

        return boundaryCondition;
    }
};

#endif
