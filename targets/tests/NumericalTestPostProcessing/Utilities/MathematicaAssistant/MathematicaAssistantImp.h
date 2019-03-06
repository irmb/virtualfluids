#ifndef MATHEMATICA_ASSISTANT_IMP_H
#define MATHEMATICA_ASSISTANT_IMP_H

#include "MathematicaAssistant.h"

class MathematicaAssistantImp : public MathematicaAssistant
{
public:
	virtual void makeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile) = 0;

protected:
	MathematicaAssistantImp();
};
#endif 