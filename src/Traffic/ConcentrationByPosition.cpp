#include "ConcentrationByPosition.h"



ConcentrationByPosition::ConcentrationByPosition(unsigned int roadlength)
{
	concentration.resize(roadlength);
}


ConcentrationByPosition::~ConcentrationByPosition()
{
}

void ConcentrationByPosition::writeToArray(const std::vector<int>& currentCarDistribution)
{
	for (unsigned int i = 0; i < currentCarDistribution.size(); i++) {
		if (currentCarDistribution[i] >= 0) 
			concentration[i] = 1.0;
		else
			concentration[i] = 0.0;
	}

	//dispConcentration();
}

void ConcentrationByPosition::dispConcentration()
{
	for (auto cell : concentration) {
		std::cout << std::setw(4) << cell;
	}
	std::cout << std::endl;
}

