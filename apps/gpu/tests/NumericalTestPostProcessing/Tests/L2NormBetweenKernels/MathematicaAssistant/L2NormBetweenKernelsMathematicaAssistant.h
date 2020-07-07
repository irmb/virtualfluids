#ifndef L2NORM_BK_MATHEMATICA_ASSISTANT_H
#define L2NORM_BK_MATHEMATICA_ASSISTANT_H

#include "Utilities/MathematicaAssistant/MathematicaAssistantImp.h"

class MathematicaFunctionFactory;
class L2NormBetweenKernelsLogFileData;

struct SortedDataL2NormBetweenKernels {
	std::vector<std::vector<std::shared_ptr<L2NormBetweenKernelsLogFileData> > > testLogFileData;
	std::vector<std::vector<std::string> > basicListNames;
};

class L2NormBetweenKernelsMathematicaAssistant : public MathematicaAssistantImp
{
public:
	static std::shared_ptr<L2NormBetweenKernelsMathematicaAssistant> getNewInstance(std::shared_ptr<MathematicaFunctionFactory> functionFactory);

	void makeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile);

	
private:
	L2NormBetweenKernelsMathematicaAssistant();
	L2NormBetweenKernelsMathematicaAssistant(std::shared_ptr<MathematicaFunctionFactory> functionFactory);

	bool checkTestParameter(std::shared_ptr<L2NormBetweenKernelsLogFileData> logFileData1, std::shared_ptr<L2NormBetweenKernelsLogFileData> logFileData2);
	std::shared_ptr<SortedDataL2NormBetweenKernels> sortLogFileData(std::shared_ptr<LogFileDataGroup> logFileData);

	void makeL2NormMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataL2NormBetweenKernels> sortedData);
	void makeL2NormBetweenKernelsMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataL2NormBetweenKernels> sortedData);

	
};
#endif 