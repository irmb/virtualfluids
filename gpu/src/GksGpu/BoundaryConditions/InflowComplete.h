#ifndef InflowComplete_CUH
#define InflowComplete_CUH

#include <memory>

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "FlowStateData/FlowStateData.cuh"

#include "BoundaryConditions/BoundaryCondition.h"

namespace GksGpu{

struct InflowCompleteStruct
{
    uint  numberOfCells;

    uint* ghostCells;
    uint* domainCells;
    uint* secondCells;

    PrimitiveVariables prim;
};

struct VF_PUBLIC InflowComplete : public BoundaryCondition //, public IsothermalWallStruct
{
    PrimitiveVariables prim;

    InflowComplete( SPtr<DataBase> dataBase, PrimitiveVariables prim );

    virtual bool isWall() override;

    virtual bool isFluxBC() override;

    virtual bool secondCellsNeeded() override;

    virtual void runBoundaryConditionKernel(const SPtr<DataBase> dataBase,
                                            const Parameters parameters, 
                                            const uint level) override;

    InflowCompleteStruct toStruct()
    {
        InflowCompleteStruct boundaryCondition;

        boundaryCondition.numberOfCells = this->numberOfCells;

        boundaryCondition.ghostCells    = this->ghostCells;
        boundaryCondition.domainCells   = this->domainCells;
        boundaryCondition.secondCells   = this->secondCells;

        boundaryCondition.prim          = prim;

        return boundaryCondition;
    }
};

} // namespace GksGpu

#endif
