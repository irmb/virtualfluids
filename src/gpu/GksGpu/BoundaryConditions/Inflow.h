#ifndef Inflow_CUH
#define Inflow_CUH

#include <memory>

#include "VirtualFluidsDefinitions.h"
#include "GksGpu_export.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "BoundaryConditions/BoundaryCondition.h"

namespace GksGpu{

struct InflowStruct
{
    uint  numberOfCells;

    uint* ghostCells;
    uint* domainCells;
    uint* secondCells;

    Vec3 velocity;
    real lambda;
    real rho;
    real S_1;
    real S_2;

    real a0, a1, a2;
};

struct GKSGPU_EXPORT Inflow : public BoundaryCondition //, public IsothermalWallStruct
{
    Vec3 velocity;
    real lambda;
    real rho;
    real S_1;
    real S_2;

    real a0, a1, a2;

    Inflow( SPtr<DataBase> dataBase, Vec3 velocity, real lambda, real rho, real a0, real a1, real a2, real S_1 = 0.0, real S_2 = 0.0 );

    virtual bool isWall() override;

    virtual bool secondCellsNeeded() override;

    virtual void runBoundaryConditionKernel(const SPtr<DataBase> dataBase,
                                            const Parameters parameters, 
                                            const uint level) override;

    InflowStruct toStruct()
    {
        InflowStruct boundaryCondition;

        boundaryCondition.numberOfCells = this->numberOfCells;

        boundaryCondition.ghostCells      = this->ghostCells;
        boundaryCondition.domainCells     = this->domainCells;
        boundaryCondition.secondCells     = this->secondCells;

        boundaryCondition.velocity        = this->velocity;
        boundaryCondition.lambda          = this->lambda;
        boundaryCondition.rho             = this->rho;
        boundaryCondition.S_1             = this->S_1;
        boundaryCondition.S_2             = this->S_2;

        boundaryCondition.a0              = this->a0;
        boundaryCondition.a1              = this->a1;
        boundaryCondition.a2              = this->a2;

        return boundaryCondition;
    }
};

} // namespace GksGpu

#endif
