#ifndef MATHEMATICA_ASSISTANT_H
#define MATHEMATICA_ASSISTANT_H

#include <memory>
#include <vector>

class LogFileDataGroup;
class MathematicaFile;

class MathematicaAssistant
{
public:
	virtual void makeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile) = 0;
};
#endif 