#include "ConcentrationOutwriter.h"

#include <iostream>
#include <iomanip>	//formatting output streams

#include "Utilities/ConsoleColor.h"

void ConcentrationOutwriter::resetConcentrations()
{
	if (useLBMConcArray)
		for (real* p = concArrayStart; p < concArrayStart + roadLength; ++p)
			*p = 0.0;
	else
		std::fill(concentration.begin(), concentration.end(), 0.0f);
}


void ConcentrationOutwriter::putConcIntoArrayOrVector(uint index, real conc)
{
	if (useLBMConcArray) {
		real *pos = concArrayStart + index;
		*pos = conc;
	}
	else
		concentration[index] = conc;
}


void ConcentrationOutwriter::addConcToArrayOrVector(uint index, real conc)
{
	if (useLBMConcArray) {
		real *pos = concArrayStart + index;
		//if ((*pos + conc) > 1.0) *pos = 1.0;
		//else
		*pos += conc;
	}
	else
		//	if (concentration[index] + conc > 1.0) concentration[index] = 1.0;
		//	else	
		concentration[index] += conc;
}


void ConcentrationOutwriter::dispCurrentConcentrations()
{
	if (useLBMConcArray)
		for (real* p = concArrayStart; p < concArrayStart + roadLength; ++p)
			dispSingleConcentration(*p);
	else
		for (auto cell : concentration)
			dispSingleConcentration(cell);


	std::cout << std::endl;

	ConsoleColor::setDefaultWhite();
}


void ConcentrationOutwriter::dispSingleConcentration(real conc)
{
	if (conc > 0)
		ConsoleColor::setBrightRed();
	else
		ConsoleColor::setDarkGrey();
	std::cout << std::setw(4) << conc;
}
