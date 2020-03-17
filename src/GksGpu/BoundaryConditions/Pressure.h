#ifndef Pressure_CUH
#define Pressure_CUH

#include <memory>

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "BoundaryConditions/BoundaryCondition.h"

namespace GksGpu{

//struct IsothermalWallStruct : virtual public BoundaryConditionStruct
//{
//    Vec3 velocity;
//    real lambda;
//    real S;
//};

struct PressureStruct
{
    uint  numberOfCells;

    uint* ghostCells;
    uint* domainCells;
    uint* secondCells;

    real p0;
};

struct VF_PUBLIC Pressure : public BoundaryCondition
{
    real p0;

    Pressure( SPtr<DataBase> dataBase, real p0 );

    virtual bool isWall() override;

    virtual bool secondCellsNeeded() override;

    virtual void runBoundaryConditionKernel(const SPtr<DataBase> dataBase,
                                            const Parameters parameters, 
                                            const uint level) override;

    PressureStruct toStruct()
    {
        PressureStruct boundaryCondition;

        boundaryCondition.numberOfCells = this->numberOfCells;

        boundaryCondition.ghostCells    = this->ghostCells;
        boundaryCondition.domainCells   = this->domainCells;
        boundaryCondition.secondCells   = this->secondCells;

        boundaryCondition.p0            = this->p0;

        return boundaryCondition;
    }
};

} // namespace GksGpu

#endif
