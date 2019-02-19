#include "LogFileDataAssistantStrategyImp.h"

#include "Utilities/AlmostEquals.h"

bool LogFileDataAssistantStrategyImp::equalDouble(double num1, double num2)
{
	const FloatingPoint<double> lhs(num1), rhs(num2);

	if (lhs.AlmostEquals(rhs))
		return true;
	return false;
}
