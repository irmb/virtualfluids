#ifndef TIME_MATHEMATICA_ASSISTANT_H
#define TIME_MATHEMATICA_ASSISTANT_H

#include "Utilities/MathematicaAssistant/MathematicaAssistantImp.h"

class MathematicaFunctionFactory;
class MathematicaFile;

class TimeMathematicaAssistant : public MathematicaAssistantImp
{
public:
    static std::shared_ptr<TimeMathematicaAssistant> getNewInstance(std::shared_ptr<MathematicaFunctionFactory> functionFactory);

    void makeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile);

private:
    TimeMathematicaAssistant();
    TimeMathematicaAssistant(std::shared_ptr<MathematicaFunctionFactory> functionFactory);

    void makeSimulationTimeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile);
    void makeTestTimeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile);
    void makeAnalyticalWriteTimeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile);

    void makeTimeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile, std::vector<std::vector<double> > times, std::string timeName);

};
#endif 