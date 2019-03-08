#include "ConcentrationOutwriter.h"


void ConcentrationOutwriter::resetConcentrations()
{
	if (useLBMConcArray)
		for (float* p = concArrayStart; p < concArrayStart + concArraySize; ++p)
			*p = 0.0;
	else
		std::fill(concentration.begin(), concentration.end(), 0.0);
}


void ConcentrationOutwriter::putConcIntoArrayOrVector(uint index, float conc)
{
	if (useLBMConcArray) {
		float *pos = concArrayStart + index;
		*pos = conc;
	}
	else
		concentration[index] = conc;
}


void ConcentrationOutwriter::addConcToArrayOrVector(uint index, float conc)
{
	if (useLBMConcArray) {
		float *pos = concArrayStart + index;
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
		for (float* p = concArrayStart; p < concArrayStart + concArraySize; ++p)
			dispSingleConcentration(*p);
	else
		for (auto cell : concentration)
			dispSingleConcentration(cell);


	std::cout << std::endl;
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 7);
}



void ConcentrationOutwriter::dispSingleConcentration(float conc)
{
	if (conc > 0)
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 12);
	else
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 8);
	std::cout << std::setw(4) << conc;
}
