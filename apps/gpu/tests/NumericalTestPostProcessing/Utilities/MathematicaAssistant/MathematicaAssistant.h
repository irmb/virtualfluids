#ifndef MATHEMATICA_ASSISTANT_H
#define MATHEMATICA_ASSISTANT_H

#include <memory>
#include <vector>
#include <string>

enum Assistant{Phi, Ny, L2Norm, L2NormBetweenKernels, Time };

class LogFileDataGroup;
class MathematicaFile;

class MathematicaAssistant
{
public:
	virtual void makeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile) = 0;
};
#endif 