#ifndef Inflow_CUH
#define Inflow_CUH

#include <memory>

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "BoundaryConditions/BoundaryCondition.h"

struct InflowStruct
{
    uint  numberOfCells;

    uint* ghostCells;
    uint* domainCells;
    uint* secondCells;

    Vec3 velocity;
    real lambda;
    real rho;
    real S;

    real a0, a1, a2;
};

struct VF_PUBLIC Inflow : public BoundaryCondition //, public IsothermalWallStruct
{
    Vec3 velocity;
    real lambda;
    real rho;
    real S;

    real a0, a1, a2;

    Inflow( SPtr<DataBase> dataBase, Vec3 velocity, real lambda, real rho, real S, real a0, real a1, real a2 );

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
        boundaryCondition.S               = this->S;

        boundaryCondition.a0              = this->a0;
        boundaryCondition.a1              = this->a1;
        boundaryCondition.a2              = this->a2;

        return boundaryCondition;
    }
};

#endif
