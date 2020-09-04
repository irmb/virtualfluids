#ifndef Symmetry_CUH
#define Symmetry_CUH

#include <memory>

#include "VirtualFluidsDefinitions.h"
#include "GksGpu_export.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "BoundaryConditions/BoundaryCondition.h"

namespace GksGpu{

struct SymmetryStruct
{
    uint  numberOfCells;

    uint* ghostCells;
    uint* domainCells;
    uint* secondCells;

    char direction;
};

struct GKSGPU_EXPORT Symmetry : public BoundaryCondition //, public IsothermalWallStruct
{
    char direction;

    Symmetry( SPtr<DataBase> dataBase, char direction );

    virtual bool isWall() override;

    virtual bool secondCellsNeeded() override;

    virtual void runBoundaryConditionKernel(const SPtr<DataBase> dataBase,
                                            const Parameters parameters, 
                                            const uint level) override;

    SymmetryStruct toStruct()
    {
        SymmetryStruct boundaryCondition;

        boundaryCondition.numberOfCells = this->numberOfCells;

        boundaryCondition.ghostCells    = this->ghostCells;
        boundaryCondition.domainCells   = this->domainCells;
        boundaryCondition.secondCells   = this->secondCells;

        boundaryCondition.direction     = this->direction;

        return boundaryCondition;
    }
};

} // namespace GksGpu

#endif
