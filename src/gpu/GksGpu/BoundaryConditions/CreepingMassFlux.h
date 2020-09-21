#ifndef CreepingMassFlux_CUH
#define CreepingMassFlux_CUH

#include <memory>


#include "GksGpu_export.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "FlowStateData/FlowStateData.cuh"

#include "BoundaryConditions/BoundaryCondition.h"

namespace GksGpu{

struct CreepingMassFluxStruct
{
    uint  numberOfCells;

    uint* ghostCells;
    uint* domainCells;
    uint* secondCells;

    real rho;
    real velocity;
    real lambda;
};

struct GKSGPU_EXPORT CreepingMassFlux : public BoundaryCondition //, public IsothermalWallStruct
{
    real rho;
    real velocity;
    real lambda;

    CreepingMassFlux( SPtr<DataBase> dataBase, real rho, real velocity, real lambda );

    virtual bool isWall() override;

    virtual bool isFluxBC() override;

    virtual bool secondCellsNeeded() override;

    virtual void runBoundaryConditionKernel(const SPtr<DataBase> dataBase,
                                            const Parameters parameters, 
                                            const uint level) override;

    CreepingMassFluxStruct toStruct()
    {
        CreepingMassFluxStruct boundaryCondition;

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
