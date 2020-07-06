#ifndef  EnstrophyAnalyzer_H
#define  EnstrophyAnalyzer_H

#include <vector>
#include <string>

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"

class Parameter;

class VF_PUBLIC EnstrophyAnalyzer
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
