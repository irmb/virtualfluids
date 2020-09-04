#ifndef Extrapolation_CUH
#define Extrapolation_CUH

#include <memory>

#include "VirtualFluidsDefinitions.h"
#include "GksGpu_export.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "BoundaryConditions/BoundaryCondition.h"

namespace GksGpu{

struct ExtrapolationStruct
{
    uint  numberOfCells;

    uint* ghostCells;
    uint* domainCells;
    uint* secondCells;
};

struct GKSGPU_EXPORT Extrapolation : public BoundaryCondition //, public IsothermalWallStruct
{
    Extrapolation( SPtr<DataBase> dataBase );

    virtual bool isWall() override;

    virtual bool secondCellsNeeded() override;

    virtual void runBoundaryConditionKernel(const SPtr<DataBase> dataBase,
                                            const Parameters parameters, 
                                            const uint level) override;

    ExtrapolationStruct toStruct()
    {
        ExtrapolationStruct boundaryCondition;

        boundaryCondition.numberOfCells = this->numberOfCells;

        boundaryCondition.ghostCells    = this->ghostCells;
        boundaryCondition.domainCells   = this->domainCells;
        boundaryCondition.secondCells   = this->secondCells;

        return boundaryCondition;
    }
};

} // namespace GksGpu

#endif
