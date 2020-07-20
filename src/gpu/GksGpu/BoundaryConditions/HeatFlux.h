#ifndef HeatFlux_CUH
#define HeatFlux_CUH

#include <memory>

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "FlowStateData/FlowStateData.cuh"

#include "BoundaryConditions/BoundaryCondition.h"

namespace GksGpu{

struct HeatFluxStruct
{
    uint  numberOfCells;

    uint* ghostCells;
    uint* domainCells;
    uint* secondCells;

    real HRRPUA;
};

struct VIRTUALFLUIDS_GPU_EXPORT HeatFlux : public BoundaryCondition //, public IsothermalWallStruct
{
    real HRRPUA;

    HeatFlux( SPtr<DataBase> dataBase, real HRRPUA );

    virtual bool isWall() override;

    virtual bool isFluxBC() override;

    virtual bool secondCellsNeeded() override;

    virtual void runBoundaryConditionKernel(const SPtr<DataBase> dataBase,
                                            const Parameters parameters, 
                                            const uint level) override;

    HeatFluxStruct toStruct()
    {
        HeatFluxStruct boundaryCondition;

        boundaryCondition.numberOfCells = this->numberOfCells;

        boundaryCondition.ghostCells    = this->ghostCells;
        boundaryCondition.domainCells   = this->domainCells;
        boundaryCondition.secondCells   = this->secondCells;

        boundaryCondition.HRRPUA        = this->HRRPUA;

        return boundaryCondition;
    }
};

} // namespace GksGpu

#endif
