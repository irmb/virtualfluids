#ifndef ConcreteHeatFlux_CUH
#define ConcreteHeatFlux_CUH

#include <memory>

#include <thrust/device_vector.h>

#include "VirtualFluidsDefinitions.h"
#include "GksGpu_export.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "FlowStateData/FlowStateData.cuh"

#include "BoundaryConditions/BoundaryCondition.h"

namespace GksGpu{

struct ConcreteHeatFluxStruct
{
    uint  numberOfCells;

    uint* ghostCells;
    uint* domainCells;
    uint* secondCells;

    uint  numberOfPoints;

    real* temperatures;

    real* heatConductivity;

    real temperatureConductivity;
    real density;
    real specificHeatCapacity;

    real L;
    real ambientTemperature;
};

struct GKSGPU_EXPORT ConcreteHeatFlux : public BoundaryCondition //, public IsothermalWallStruct
{
    real* temperatures;

    uint numberOfPoints;

    real temperatureConductivity;
    real density;
    real specificHeatCapacity;

    real L;
    real ambientTemperature;

    std::vector<uint> ghostCellsHost ;
    std::vector<uint> domainCellsHost;
    std::vector<uint> secondCellsHost;

    std::vector<real> temperaturesHost;

    ~ConcreteHeatFlux();

    ConcreteHeatFlux( SPtr<DataBase> dataBase, uint numberOfPoints, real temperatureConductivity, real density, real specificHeatCapacity, real L, real ambientTemperature );

    void init();

    void download();

    virtual bool isWall() override;

    virtual bool isInsulated() override;

    virtual bool isFluxBC() override;

    virtual bool secondCellsNeeded() override;

    virtual void runBoundaryConditionKernel(const SPtr<DataBase> dataBase,
                                            const Parameters parameters, 
                                            const uint level) override;

    ConcreteHeatFluxStruct toStruct()
    {
        ConcreteHeatFluxStruct boundaryCondition;

        boundaryCondition.numberOfCells  = this->numberOfCells;

        boundaryCondition.ghostCells     = this->ghostCells;
        boundaryCondition.domainCells    = this->domainCells;
        boundaryCondition.secondCells    = this->secondCells;

        boundaryCondition.temperatures   = this->temperatures;
        boundaryCondition.numberOfPoints = this->numberOfPoints;

        boundaryCondition.temperatureConductivity = this->temperatureConductivity;
        boundaryCondition.density                 = this->density;
        boundaryCondition.specificHeatCapacity    = this->specificHeatCapacity;

        boundaryCondition.L                  = this->L;
        boundaryCondition.ambientTemperature = this->ambientTemperature;

        return boundaryCondition;
    }

    void writeVTKFile( SPtr<DataBase> dataBase, Parameters& parameters, std::string filename );
};

} // namespace GksGpu

#endif
