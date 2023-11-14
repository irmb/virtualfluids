#ifndef MATHEMATICA_FUNCTION_FACTORY_IMP_H
#define MATHEMATICA_FUNCTION_FACTORY_IMP_H

#include "MathematicaFunctionFactory.h"

class MathematicaFunctionFactoryImp : public MathematicaFunctionFactory
{
public:
    static std::shared_ptr<MathematicaFunctionFactory> getNewInstance();

    std::shared_ptr<MathematicaPointList> makeMathematicaPointList(std::shared_ptr<MathematicaFile> file, std::string listName, std::vector<std::shared_ptr<DataPoint> > plotData);
    std::shared_ptr<MathematicaListPlot> makeMathematicaListPlot(std::shared_ptr<MathematicaFile> file, std::vector<std::shared_ptr<MathematicaPointList>> pointList, std::string plotType, std::string xAxes, std::string yAxes);
    std::shared_ptr<MathematicaListOfLists> makeMathematicaListOfLists(std::shared_ptr<MathematicaFile> file, std::string listName, std::vector<std::vector<double>> listOfLists);

private:
    MathematicaFunctionFactoryImp();

};
#endif