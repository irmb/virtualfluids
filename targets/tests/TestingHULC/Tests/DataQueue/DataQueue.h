#ifndef DATASTORAGEFORTEST_H
#define DATASTORAGEFORTEST_H

#include <string>

struct DataQueue
{
	double a, b;
	int la, lb;
	bool expected;
	int numberOfTests;
	std::string testName;
	std::string valueName;
	double orderOfAccuracy, minOrderOfAccuracy;

	struct DataQueue()
	{
		a = 0.0;
		b = 0.0;
		la = 0;
		lb = 0;
		testName = "noTestName";
		valueName = "noValueName";
		expected = false;
		numberOfTests = 0;
		orderOfAccuracy = 0.0;
		minOrderOfAccuracy = 0.0;
	}
};
#endif