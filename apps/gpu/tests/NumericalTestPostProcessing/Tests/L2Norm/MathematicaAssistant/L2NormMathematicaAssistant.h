#ifndef L2NORM_MATHEMATICA_ASSISTANT_H
#define L2NORM_MATHEMATICA_ASSISTANT_H

#include "Utilities/MathematicaAssistant/MathematicaAssistantImp.h"

class MathematicaFunctionFactory;
class L2NormLogFileData;

struct SortedDataL2Norm {
    std::vector<std::vector<std::shared_ptr<L2NormLogFileData> > > testLogFileData;
    std::vector<std::vector<std::string> > basicListNames;
};

class L2NormMathematicaAssistant : public MathematicaAssistantImp
{
public:
    static std::shared_ptr<L2NormMathematicaAssistant> getNewInstance(std::shared_ptr<MathematicaFunctionFactory> functionFactory);

    void makeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile);

    
private:
    L2NormMathematicaAssistant();
    L2NormMathematicaAssistant(std::shared_ptr<MathematicaFunctionFactory> functionFactory);

    bool checkTestParameter(std::shared_ptr<L2NormLogFileData> logFileData1, std::shared_ptr<L2NormLogFileData> logFileData2);
    std::shared_ptr<SortedDataL2Norm> sortLogFileData(std::shared_ptr<LogFileDataGroup> logFileData);

    void makeL2NormDiffMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataL2Norm> sortedData);
    void makeL2NormAllTimeStepsMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataL2Norm> sortedData);
    void makeL2NormBasicTimeStepMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataL2Norm> sortedData);
    void makeL2NormDivergentTimeStepMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataL2Norm> sortedData);

    
};
#endif 