#ifndef  KineticEngergyAnalyzer_H
#define  KineticEngergyAnalyzer_H

#include <vector>
#include <string>

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"

#include "FlowStateData/FlowStateData.cuh"

struct DataBase;

class VF_PUBLIC KineticEnergyAnalyzer
{
private:

    SPtr<DataBase> dataBase;

    uint outputIter;

    uint analyzeIter;

    std::vector<real> kineticEnergyTimeSeries;

public:

    KineticEnergyAnalyzer( SPtr<DataBase> dataBase, uint analyzeIter, uint outputIter );

    bool run( uint iter );

    void writeToFile( std::string filename );

};

#endif
