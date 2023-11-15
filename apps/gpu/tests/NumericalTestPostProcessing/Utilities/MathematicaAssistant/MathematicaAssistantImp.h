#ifndef MATHEMATICA_ASSISTANT_IMP_H
#define MATHEMATICA_ASSISTANT_IMP_H

#include "MathematicaAssistant.h"

class MathematicaFunctionFactory;

class MathematicaAssistantImp : public MathematicaAssistant
{
public:
    virtual void makeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile) = 0;

protected:
    MathematicaAssistantImp();
    MathematicaAssistantImp(std::shared_ptr<MathematicaFunctionFactory> functionFactory);

    std::vector<std::string> finalizeListNames(std::vector<std::string> basicListNames, std::string testName, std::string dataToCalc);
    std::vector<std::string> finalizeListNames(std::vector<std::string> basicListNames, std::string testName, std::string dataToCalc, std::string normalizeData);
    std::vector<std::string> finalizeListNames(std::vector<std::string> basicListNames, std::string testName, std::string dataToCalc, std::string normalizeData, std::string timeStep);
    std::vector<std::string> finalizeListNames(std::vector<std::string> basicListNames, std::string testName, std::string dataToCalc, std::string normalizeData, std::vector<std::string> timeSteps);
    std::vector<std::string> finalizeListNames(std::vector<std::string> basicListNames, std::string testName, std::string dataToCalc, std::string normalizeData, std::vector<std::string> timeSteps, std::vector<std::string> basicKernels);
    std::string finalizeListName(std::string basicListName, std::string testName, std::string dataToCalc);

    void addListLogLogPlotToMathematicaFile(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::vector<std::string> listNames, std::vector<std::vector<double> > xAxesData, std::vector<std::vector<double> > yAxesData, std::string labelXAxes, std::string labelYAxes);
    void addListOfListsToMathematicaFile(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::string listName, std::vector<std::vector<double> > listOfLists);

    void addSecondOrderOfAccuracyRef(std::vector<std::vector<double> > &xAxesData, std::vector<std::vector<double> > &yAxesData, std::vector<std::string> &listNames);
    void addFourthOrderOfAccuracyRef(std::vector<std::vector<double> > &xAxesData, std::vector<std::vector<double> > &yAxesData, std::vector<std::string> &listNames);

    std::shared_ptr<MathematicaFunctionFactory> functionFactory;
};
#endif 