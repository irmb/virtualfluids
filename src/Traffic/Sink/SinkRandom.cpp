#include "SinkRandom.h"

#include "Utilities/invalidInput_error.h"

SinkRandom::SinkRandom(uint sinkIndex, real sinkBlockedPossibility)
{
	data.sinkIndex = sinkIndex;

	try {
		if (sinkBlockedPossibility >= 0 && sinkBlockedPossibility <= 1) {
			data.sinkBlockedPossibility = sinkBlockedPossibility;
		}
		else {
			throw invalidInput_error("possibility of the sink being blocked should be between 0 and 1");
		}
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cin.get();
		exit(EXIT_FAILURE);
	};
}


bool SinkRandom::carCanEnter()
{	
	return  !(distFloat(engine) < data.sinkBlockedPossibility);
}


real SinkRandom::getPossibilityBeingBlocked() const
{
	return data.sinkBlockedPossibility;
}


uint SinkRandom::getIndex() const
{
	return data.sinkIndex;
}

