#include "L2NormCalculatorImp.h"

#include "Utilities/Calculator/AlmostEquals.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

L2NormCalculatorImp::L2NormCalculatorImp(std::string errorMessage) : errorMessage(errorMessage)
{
}

bool L2NormCalculatorImp::equalDouble(double num1, double num2)
{
    const FloatingPoint<double> lhs(num1), rhs(num2);

    if (lhs.AlmostEquals(rhs))
        return true;
    return false;
}

double L2NormCalculatorImp::calcCounter(std::vector<double> basicData, std::vector<double> divergentData, std::vector<unsigned int> level, double lx, double lz)
{
    double counter = 0.0;
    for (int i = 0; i < basicData.size(); i++) {
        double area = (1 / pow(2.0, level.at(i))) * (1 / pow(2.0, level.at(i)));
        counter += ((divergentData.at(i) - basicData.at(i))*(divergentData.at(i) - basicData.at(i))) * area;
    }
    return counter;
}

std::string L2NormCalculatorImp::getErrorMessage()
{
    return errorMessage;
}

L2NormCalculatorImp::L2NormCalculatorImp()
{
    
}