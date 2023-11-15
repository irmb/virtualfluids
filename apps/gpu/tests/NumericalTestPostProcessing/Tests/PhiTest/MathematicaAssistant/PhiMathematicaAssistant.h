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
    
    
    void makePhiDiffMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataPhi> sortedData);
    void makeOrderOfAccuracyMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataPhi> sortedData);

};
#endif 