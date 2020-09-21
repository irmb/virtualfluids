#ifndef SalinasVazquez_CUH
#define SalinasVazquez_CUH

#include <memory>


#include "GksGpu_export.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "BoundaryConditions/BoundaryCondition.h"

namespace GksGpu{

struct SalinasVazquezStruct
{
    uint  numberOfCells;

    uint* ghostCells;
    uint* domainCells;
    uint* secondCells;

    real lambdaMX;
    real lambdaPX;

    real a0, a1, a2, a3;

    bool useSecondCells;
};

struct GKSGPU_EXPORT SalinasVazquez : public BoundaryCondition //, public IsothermalWallStruct
{
    real lambdaMX;
    real lambdaPX;

    real a0, a1, a2, a3;

    bool useSecondCells;

    SalinasVazquez( SPtr<DataBase> dataBase, real lambdaMX, real lambdaPX, real a0, real a1, real a2, real a3, bool useSecondCells );

    virtual bool isWall() override;

    virtual bool secondCellsNeeded() override;

    virtual void runBoundaryConditionKernel(const SPtr<DataBase> dataBase,
                                            const Parameters parameters, 
                                            const uint level) override;

    SalinasVazquezStruct toStruct()
    {
        SalinasVazquezStruct boundaryCondition;

        boundaryCondition.numberOfCells = this->numberOfCells;

        boundaryCondition.ghostCells      = this->ghostCells;
        boundaryCondition.domainCells     = this->domainCells;
        boundaryCondition.secondCells     = this->secondCells;

        boundaryCondition.lambdaMX        = this->lambdaMX;
        boundaryCondition.lambdaPX        = this->lambdaPX;

        boundaryCondition.a0              = this->a0;
        boundaryCondition.a1              = this->a1;
        boundaryCondition.a2              = this->a2;
        boundaryCondition.a3              = this->a3;

        boundaryCondition.useSecondCells  = this->useSecondCells;

        return boundaryCondition;
    }
};

} // namespace GksGpu

#endif
