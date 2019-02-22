#include "SinkSimple.h"



SinkSimple::SinkSimple(unsigned int sinkIndex, float sinkBlockedPossibility)
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
	catch (const exception& e) { 
		cerr << e.what() << endl;
		cin.get();
		exit(EXIT_FAILURE);
	};
}


bool SinkSimple::carCanEnter()
{	
	return  !(distFloat(engine) < data.sinkBlockedPossibility);
}


float SinkSimple::getPossibilityBeingBlocked() const
{
	return data.sinkBlockedPossibility;
}


unsigned int SinkSimple::getIndex() const
{
	return data.sinkIndex;
}

