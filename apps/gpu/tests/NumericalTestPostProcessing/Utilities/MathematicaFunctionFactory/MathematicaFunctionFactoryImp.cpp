#include "MathematicaFunctionFactoryImp.h"

#include "Utilities/MathematicaFile/MathematicaFile.h"

#include "Utilities/MathematicaFunction/LinePlot/MathematicaListPlotImp.h"
#include "Utilities/MathematicaFunction/PointList/MathematicaPointListImp.h"
#include "Utilities/MathematicaFunction/ListOfLists/MathematicaListOfListsImp.h"

std::shared_ptr<MathematicaFunctionFactory> MathematicaFunctionFactoryImp::getNewInstance()
{
	return std::shared_ptr<MathematicaFunctionFactory>(new MathematicaFunctionFactoryImp());
}

std::shared_ptr<MathematicaPointList> MathematicaFunctionFactoryImp::makeMathematicaPointList(std::shared_ptr<MathematicaFile> file, std::string listName, std::vector<std::shared_ptr<DataPoint>> plotData)
{
	std::shared_ptr<MathematicaPointList> mathPointList = MathematicaPointListImp::getNewInstance(listName, plotData);
	file->addMathematicaFunction(mathPointList);
	return mathPointList;
}

std::shared_ptr<MathematicaListPlot> MathematicaFunctionFactoryImp::makeMathematicaListPlot(std::shared_ptr<MathematicaFile> file, std::vector<std::shared_ptr<MathematicaPointList>> pointList, std::string plotType, std::string xAxes, std::string yAxes)
{
	std::shared_ptr<MathematicaListPlot> listLinePlot = MathematicaListPlotImp::getNewInstance(pointList, plotType, xAxes, yAxes);
	file->addMathematicaFunction(listLinePlot);
	return listLinePlot;
}

std::shared_ptr<MathematicaListOfLists> MathematicaFunctionFactoryImp::makeMathematicaListOfLists(std::shared_ptr<MathematicaFile> file, std::string listName, std::vector<std::vector<double>> listOfLists)
{
	std::shared_ptr<MathematicaListOfLists> listOfList = MathematicaListOfListsImp::getNewInstance(listName, listOfLists);
	file->addMathematicaFunction(listOfList);
	return listOfList;
}

MathematicaFunctionFactoryImp::MathematicaFunctionFactoryImp()
{
}
