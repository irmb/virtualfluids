#ifndef NY_MATHEMATICA_ASSISTANT_H
#define NY_MATHEMATICA_ASSISTANT_H

#include "Utilities/MathematicaAssistant/MathematicaAssistantImp.h"

class MathematicaFunctionFactory;
class NyLogFileData;

struct SortedDataNy {
    std::vector<std::vector<std::shared_ptr<NyLogFileData> > > testLogFileData;
    std::vector<std::vector<std::string> > basicListNames;
};

class NyMathematicaAssistant : public MathematicaAssistantImp
{
public:
    static std::shared_ptr<NyMathematicaAssistant> getNewInstance(std::shared_ptr<MathematicaFunctionFactory> functionFactory);

    void makeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile);

    
private:
    NyMathematicaAssistant();
    NyMathematicaAssistant(std::shared_ptr<MathematicaFunctionFactory> functionFactory);

    bool checkTestParameter(std::shared_ptr<NyLogFileData> logFileData1, std::shared_ptr<NyLogFileData> logFileData2);
    std::shared_ptr<SortedDataNy> sortLogFileData(std::shared_ptr<LogFileDataGroup> logFileData);

    void makeNyDiffMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataNy> sortedData);
    void makeOrderOfAccuracyMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataNy> sortedData);

};
#endif 