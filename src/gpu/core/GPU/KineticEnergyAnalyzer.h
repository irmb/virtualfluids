#ifndef  KineticEngergyAnalyzer_H
#define  KineticEngergyAnalyzer_H

#include <vector>
#include <string>



#include <basics/PointerDefinitions.h>
#include <basics/DataTypes.h>


class Parameter;

class KineticEnergyAnalyzer
{
private:

    SPtr<Parameter> para;

    uint analyzeIter;

    std::vector<real> kineticEnergyTimeSeries;

public:

    KineticEnergyAnalyzer( SPtr<Parameter> para, uint analyzeIter );

    bool run( uint iter );

    void writeToFile( std::string filename );

};

#endif
