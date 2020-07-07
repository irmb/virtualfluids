#ifndef L2NORM_CALCULATOR_H
#define L2NORM_CALCULATOR_H

#include <vector>
#include <string>


class L2NormCalculator
{
public:
	virtual double calc(std::vector<double> basicData, std::vector<double> divergentData, std::vector<unsigned int> level, double lx, double lz, double l0) = 0;
	virtual std::string getErrorMessage() = 0;

};
#endif