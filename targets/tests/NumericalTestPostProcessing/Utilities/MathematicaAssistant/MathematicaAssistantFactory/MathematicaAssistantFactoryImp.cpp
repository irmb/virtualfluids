#include "MathematicaAssistantFactoryImp.h"

#include "Tests/PhiTest/MathematicaAssistant/PhiMathematicaAssistant.h"
#include "Tests/NyTest/MathematicaAssistant/NyMathematicaAssistant.h"
#include "Tests/L2Norm/MathematicaAssistant/L2NormMathematicaAssistant.h"
#include "Tests/L2NormBetweenKernels/MathematicaAssistant/L2NormBetweenKernelsMathematicaAssistant.h"

#include "Utilities/MathematicaAssistant/TimeAssistant/TimeMathematicaAssistant.h"

std::shared_ptr<MathematicaAssistantFactory> MathematicaAssistantFactoryImp::getNewInstance()
{
	return std::shared_ptr<MathematicaAssistantFactory>(new MathematicaAssistantFactoryImp());
}

std::vector<std::shared_ptr<MathematicaAssistant>> MathematicaAssistantFactoryImp::makeMathematicaAssistants(std::vector<Assistant> types, std::shared_ptr<MathematicaFunctionFactory> functionFactory)
{
	std::vector<std::shared_ptr<MathematicaAssistant>> myAssistants;

	for(int i = 0; i < types.size(); i++){
		switch (types.at(i))
		{
		case Phi:
			myAssistants.push_back(PhiMathematicaAssistant::getNewInstance(functionFactory));
			break;
		case Ny:
			myAssistants.push_back(NyMathematicaAssistant::getNewInstance(functionFactory));
			break;
		case L2Norm:
			myAssistants.push_back(L2NormMathematicaAssistant::getNewInstance(functionFactory));
			break;
		case L2NormBetweenKernels:
			myAssistants.push_back(L2NormBetweenKernelsMathematicaAssistant::getNewInstance(functionFactory));
			break;
		case Time:
			myAssistants.push_back(TimeMathematicaAssistant::getNewInstance(functionFactory));
			break;
		default:
			break;
		}

	}


	return myAssistants;
}

MathematicaAssistantFactoryImp::MathematicaAssistantFactoryImp()
{
}



