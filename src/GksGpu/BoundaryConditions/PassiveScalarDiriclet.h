#ifndef PassiveScalarDiriclet_CUH
#define PassiveScalarDiriclet_CUH

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

struct PassiveScalarDiricletStruct
{
    uint  numberOfCells;

    uint* ghostCells;
    uint* domainCells;
    uint* secondCells;

    real S_1;
    real S_2;
};

struct VF_PUBLIC PassiveScalarDiriclet : public BoundaryCondition
{
    real S_1;
    real S_2;

    PassiveScalarDiriclet( SPtr<DataBase> dataBase, real S_1, real S_2 );

    virtual bool isWall() override;

    virtual bool secondCellsNeeded() override;

    virtual void runBoundaryConditionKernel(const SPtr<DataBase> dataBase,
                                            const Parameters parameters, 
                                            const uint level) override;

    PassiveScalarDiricletStruct toStruct()
    {
        PassiveScalarDiricletStruct boundaryCondition;

        boundaryCondition.numberOfCells = this->numberOfCells;

        boundaryCondition.ghostCells    = this->ghostCells;
        boundaryCondition.domainCells   = this->domainCells;
        boundaryCondition.secondCells   = this->secondCells;

        boundaryCondition.S_1           = this->S_1;
        boundaryCondition.S_2           = this->S_2;

        return boundaryCondition;
    }
};

#endif
