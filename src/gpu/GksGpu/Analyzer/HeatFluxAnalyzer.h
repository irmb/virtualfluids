#ifndef  HeatFluxAnalyzer_H
#define  HeatFluxAnalyzer_H

#include <vector>
#include <string>


#include "GksGpu_export.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"

#include "GksGpu/BoundaryConditions/BoundaryCondition.h"

#include "FlowStateData/FlowStateData.cuh"

#include "Parameters/Parameters.h"

namespace GksGpu {

struct DataBase;

class GKSGPU_EXPORT HeatFluxAnalyzer
{
private:

    SPtr<DataBase> dataBase;
    SPtr<GksGpu::BoundaryCondition> boundaryCondition;

    uint outputIter;

    uint analyzeIter;

    std::vector<real> heatFluxTimeSeries;

    real lambdaHot;
    real lambdaCold;

    real L;

public:

    HeatFluxAnalyzer( SPtr<DataBase> dataBase, SPtr<GksGpu::BoundaryCondition> boundaryCondition, uint analyzeIter, uint outputIter, real lambdaHot, real lambdaCold, real L );

    bool run( uint iter, Parameters parameters );

    void writeToFile( std::string filename );

};

} // namespace GksGpu

#endif
