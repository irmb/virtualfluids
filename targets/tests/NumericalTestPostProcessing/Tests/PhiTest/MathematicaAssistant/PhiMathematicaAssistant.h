#ifndef PHI_MATHEMATICA_ASSISTANT_H
#define PHI_MATHEMATICA_ASSISTANT_H

#include "Utilities/MathematicaAssistant/MathematicaAssistantImp.h"

class MathematicaFunctionFactory;
class PhiLogFileData;

struct SortedDataPhi {
	std::vector<std::vector<std::shared_ptr<PhiLogFileData> > > testLogFileData;
	std::vector<std::vector<std::string> > basicListNames;
};

class PhiMathematicaAssistant : public MathematicaAssistantImp
{
public:
	static std::shared_ptr<PhiMathematicaAssistant> getNewInstance(std::shared_ptr<MathematicaFunctionFactory> functionFactory);

	void makeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile);

	
private:
	PhiMathematicaAssistant();
	PhiMathematicaAssistant(std::shared_ptr<MathematicaFunctionFactory> functionFactory);

	bool checkTestParameter(std::shared_ptr<PhiLogFileData> logFileData1, std::shared_ptr<PhiLogFileData> logFileData2);
	std::shared_ptr<SortedDataPhi> sortLogFileData(std::shared_ptr<LogFileDataGroup> logFileData);
	void addListLogLogPlotToMathematicaFile(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::vector<std::string> listNames, std::vector<std::vector<double> > xAxesData, std::vector<std::vector<double> > yAxesData, std::string labelXAxes, std::string labelYAxes);
	void addListOfListsToMathematicaFile(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::string listNames, std::vector<std::vector<double> > listOfLists);
	std::vector<std::string> finalizeListNames(std::vector<std::string> basicListNames, std::string dataToCalc, std::string testName);
	void makePhiDiffMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataPhi> sortedData);
	void makeOrderOfAccuracyMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataPhi> sortedData);

	std::shared_ptr<MathematicaFunctionFactory> functionFactory;
};
#endif 