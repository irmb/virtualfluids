#ifndef  EnstrophyAnalyzer_H
#define  EnstrophyAnalyzer_H

#include <vector>
#include <string>

#include "VirtualFluidsDefinitions.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"

#include "Parameters/Parameters.h"

#include "FlowStateData/FlowStateData.cuh"

namespace GksGpu {

struct DataBase;

class VIRTUALFLUIDS_GPU_EXPORT EnstrophyAnalyzer
{
private:

    SPtr<DataBase> dataBase;

    Parameters parameters;

    uint outputIter;

    uint analyzeIter;

    std::vector<real> enstrophyTimeSeries;

public:

    EnstrophyAnalyzer( SPtr<DataBase> dataBase, Parameters parameters, uint analyzeIter, uint outputIter );

    bool run( uint iter );

    void writeToFile( std::string filename );

};

} // namespace GksGpu

#endif
