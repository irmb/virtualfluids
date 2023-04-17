#ifndef  EnstrophyAnalyzer_H
#define  EnstrophyAnalyzer_H

#include <vector>
#include <string>



#include "PointerDefinitions.h"
#include "DataTypes.h"
#include "VirtualFluids_GPU_export.h"

class Parameter;

class VIRTUALFLUIDS_GPU_EXPORT EnstrophyAnalyzer
{
private:

	SPtr<Parameter> para;

    uint analyzeIter;

    std::vector<real> enstrophyTimeSeries;

public:

    EnstrophyAnalyzer( SPtr<Parameter> para, uint analyzeIter );

    bool run( uint iter );

	void writeToFile( std::string filename );

};

#endif
