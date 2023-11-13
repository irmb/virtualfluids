#ifndef  EnstrophyAnalyzer_H
#define  EnstrophyAnalyzer_H

#include <vector>
#include <string>



#include <basics/PointerDefinitions.h>
#include <basics/DataTypes.h>


class Parameter;

class EnstrophyAnalyzer
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
