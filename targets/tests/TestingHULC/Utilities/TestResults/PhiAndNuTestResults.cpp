#include "PhiAndNuTestResults.h"

std::shared_ptr<PhiAndNuTestResults> PhiAndNuTestResults::getNewInstance(std::string aTestName)
{
	return std::shared_ptr<PhiAndNuTestResults>(new PhiAndNuTestResults(aTestName));
}

void PhiAndNuTestResults::evaluate()
{
	if (nuDiff.size() > 1 && phiDiff.size() > 1) {
		const int sizeNuDiff = nuDiff.size();
		double ordOfAccNuDiff = log(nuDiff.end()[-1] / nuDiff.back()) / log(lx.back() / lx.end()[-1]);
		double ordOfAccPhiDiff = log(phiDiff.end()[-1] / phiDiff.back()) / log(lx.back()/ lx.end()[-1]);
		orderOfAccuracyNuDiff.push_back(ordOfAccNuDiff);
		orderOfAccuracyPhiDiff.push_back(ordOfAccPhiDiff);
		makeLastTestOutput();
	}
}

void PhiAndNuTestResults::add(double phiDiff, double nuDiff, double lx)
{
	this->nuDiff.push_back(nuDiff);
	this->phiDiff.push_back(phiDiff);
	this->lx.push_back(lx);
}

PhiAndNuTestResults::PhiAndNuTestResults(std::string aTestName) :testName(aTestName)
{
	phiDiff.resize(0);
	nuDiff.resize(0);
	lx.resize(0);
	orderOfAccuracyPhiDiff.resize(0);
	orderOfAccuracyNuDiff.resize(0);
}

void PhiAndNuTestResults::makeLastTestOutput()
{
	for (int i = 0; i < orderOfAccuracyNuDiff.size(); i++) {
		std::cout << orderOfAccuracyNuDiff.at(i) << std::endl;
	}
}
