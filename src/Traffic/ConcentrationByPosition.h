#pragma once
#include <VirtualFluidsDefinitions.h>

#include "ConcentrationOutwriter.h"

class VF_PUBLIC ConcentrationByPosition:
	public ConcentrationOutwriter
{
public:
	ConcentrationByPosition(unsigned int roadlength);
	~ConcentrationByPosition();

	virtual void writeToArray(const std::vector<int> & currentCarDistribution);

private:
	void dispConcentration();

private:
	std::vector<double> concentration;
};

