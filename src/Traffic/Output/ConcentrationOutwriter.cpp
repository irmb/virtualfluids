#include "ConcentrationOutwriter.h"

#include <iostream>
#include <iomanip>	//formatting output streams
#include <windows.h> //for colourful console output

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


void ConcentrationOutwriter::dispConcentrations()
{
	if (useLBMConcArray)
		for (real* p = concArrayStart; p < concArrayStart + roadLength; ++p)
			dispSingleConcentration(*p);
	else
		for (auto cell : concentration)
			dispSingleConcentration(cell);


	std::cout << std::endl;
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 7);
}



void ConcentrationOutwriter::dispSingleConcentration(real conc)
{
	if (conc > 0)
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 12);
	else
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 8);
	std::cout << std::setw(4) << conc;
}
