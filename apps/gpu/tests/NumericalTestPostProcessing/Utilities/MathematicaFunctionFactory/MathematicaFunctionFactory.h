#ifndef MATHEMATICA_FUNCTION_FACTORY_H
#define MATHEMATICA_FUNCTION_FACTORY_H

#include <memory>
#include <string>
#include <vector>

class DataPoint;
class MathematicaFile;
class MathematicaListPlot;
class MathematicaListOfLists;
class MathematicaPointList;

class MathematicaFunctionFactory
{
public:
	virtual std::shared_ptr<MathematicaPointList> makeMathematicaPointList(std::shared_ptr<MathematicaFile> file, std::string listName, std::vector<std::shared_ptr<DataPoint> > plotData) = 0;
	virtual std::shared_ptr<MathematicaListPlot> makeMathematicaListPlot(std::shared_ptr<MathematicaFile> file, std::vector<std::shared_ptr<MathematicaPointList>> pointList, std::string plotType, std::string xAxes, std::string yAxes) = 0;
	virtual std::shared_ptr<MathematicaListOfLists> makeMathematicaListOfLists(std::shared_ptr<MathematicaFile> file, std::string listName, std::vector<std::vector<double>> listOfLists) = 0;
};
#endif