#ifndef MassCompensation_CUH
#define MassCompensation_CUH

#include <memory>

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "FlowStateData/FlowStateData.cuh"

#include "BoundaryConditions/BoundaryCondition.h"

namespace GksGpu{

struct MassCompensationStruct
{
    uint  numberOfCells;

    uint* ghostCells;
    uint* domainCells;
    uint* secondCells;

    real rho;
    real velocity;
    real lambda;
};

struct VF_PUBLIC MassCompensation : public BoundaryCondition //, public IsothermalWallStruct
{
    real rho;
    real velocity;
    real lambda;

    MassCompensation( SPtr<DataBase> dataBase, real rho, real velocity, real lambda );

    virtual bool isWall() override;

    virtual bool isFluxBC() override;

    virtual bool secondCellsNeeded() override;

    virtual void runBoundaryConditionKernel(const SPtr<DataBase> dataBase,
                                            const Parameters parameters, 
                                            const uint level) override;

    MassCompensationStruct toStruct()
    {
        MassCompensationStruct boundaryCondition;

        boundaryCondition.numberOfCells = this->numberOfCells;

        boundaryCondition.ghostCells    = this->ghostCells;
        boundaryCondition.domainCells   = this->domainCells;
        boundaryCondition.secondCells   = this->secondCells;

        boundaryCondition.rho           = this->rho;
        boundaryCondition.velocity      = this->velocity;
        boundaryCondition.lambda        = this->lambda;

        return boundaryCondition;
    }
};

} // namespace GksGpu

#endif
