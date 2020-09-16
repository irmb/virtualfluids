#ifndef IsothermalWall_CUH
#define IsothermalWall_CUH

#include <memory>


#include "GksGpu_export.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "BoundaryConditions/BoundaryCondition.h"

namespace GksGpu{

struct IsothermalWallStruct
{
    uint  numberOfCells;

    uint* ghostCells;
    uint* domainCells;
    uint* secondCells;

    Vec3 velocity;
    real lambda;
    real S_1;
    real S_2;

    bool useSecondCells;
};

struct GKSGPU_EXPORT IsothermalWall : public BoundaryCondition //, public IsothermalWallStruct
{
    Vec3 velocity;
    real lambda;
    real S_1;
    real S_2;

    bool useSecondCells;

    IsothermalWall( SPtr<DataBase> dataBase, Vec3 velocity, real lambda, bool useSecondCells, real S_1 = 0.0, real S_2 = 0.0 );

    virtual bool isWall() override;

    virtual bool secondCellsNeeded() override;

    virtual void runBoundaryConditionKernel(const SPtr<DataBase> dataBase,
                                            const Parameters parameters, 
                                            const uint level) override;

    IsothermalWallStruct toStruct()
    {
        IsothermalWallStruct boundaryCondition;

        boundaryCondition.numberOfCells = this->numberOfCells;

        boundaryCondition.ghostCells      = this->ghostCells;
        boundaryCondition.domainCells     = this->domainCells;
        boundaryCondition.secondCells     = this->secondCells;

        boundaryCondition.velocity        = this->velocity;
        boundaryCondition.lambda          = this->lambda;
        boundaryCondition.S_1             = this->S_1;
        boundaryCondition.S_2             = this->S_2;

        boundaryCondition.useSecondCells  = this->useSecondCells;

        return boundaryCondition;
    }
};

} // namespace GksGpu

#endif
