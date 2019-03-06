#include "ConcentrationOutwriter.h"

void ConcentrationOutwriter::dispConcentration()
{
	for (auto cell : concentration) {
		if(cell>0)
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 12);
		else
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 8);
		std::cout << std::setw(4) << cell;
	}
	std::cout << std::endl;
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 7);
}


void ConcentrationOutwriter::resetConcVector()
{
	std::fill(concentration.begin(), concentration.end(), 0);
}


const std::vector<float> & ConcentrationOutwriter::getConcentrations()
{
	return concentration;
}