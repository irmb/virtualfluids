#ifndef Open_CUH
#define Open_CUH

#include <memory>


#include "GksGpu_export.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "FlowStateData/FlowStateData.cuh"

#include "BoundaryConditions/BoundaryCondition.h"

namespace GksGpu{

struct OpenStruct
{
    uint  numberOfCells;

    uint* ghostCells;
    uint* domainCells;
    uint* secondCells;

    PrimitiveVariables prim;

    real velocityLimiter;
};

struct GKSGPU_EXPORT Open : public BoundaryCondition //, public IsothermalWallStruct
{
    PrimitiveVariables prim;

    real velocityLimiter;

    Open( SPtr<DataBase> dataBase, PrimitiveVariables prim, real velocityLimiter );

    virtual bool isWall() override;

    virtual bool secondCellsNeeded() override;

    virtual void runBoundaryConditionKernel(const SPtr<DataBase> dataBase,
                                            const Parameters parameters, 
                                            const uint level) override;

    OpenStruct toStruct()
    {
        OpenStruct boundaryCondition;

        boundaryCondition.numberOfCells   = this->numberOfCells;

        boundaryCondition.ghostCells      = this->ghostCells;
        boundaryCondition.domainCells     = this->domainCells;
        boundaryCondition.secondCells     = this->secondCells;

        boundaryCondition.prim            = this->prim;

        boundaryCondition.velocityLimiter = this->velocityLimiter;

        return boundaryCondition;
    }
};

} // namespace GksGpu

#endif
