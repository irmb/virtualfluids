#include "PhiAndNuTest.h"

#include "Utilities/TestCout/TestCoutImp.h"
#include "PhiAndNuTest.h"

std::shared_ptr<PhiAndNuTest> PhiAndNuTest::getNewInstance(std::string aTestName, double minOrderOfAccuracy, std::shared_ptr<TestCout> testOut)
{
	return std::shared_ptr<PhiAndNuTest>(new PhiAndNuTest(aTestName, minOrderOfAccuracy, testOut));
}

void PhiAndNuTest::evaluate()
{
	orderOfAccuracyNuDiff = calcOrderOfAccuracy(nuDiff);
	orderOfAccuracyPhiDiff = calcOrderOfAccuracy(phiDiff);

	nuDiffTestPassed = checkTestPassed(orderOfAccuracyNuDiff);
	phiDiffTestPassed = checkTestPassed(orderOfAccuracyPhiDiff);
	
	if (orderOfAccuracyNuDiff.size() > 0 && orderOfAccuracyPhiDiff.size() > 0)
		makeLastTestOutput();
}

void PhiAndNuTest::makeFinalOutput()
{
	for (int i = 1; i < lx.size(); i++) {
		testOut->makeTestOutput(nuDiffTestPassed.at(i - 1), testName, lx.at(i - 1), lx.at(i), "NuDiff", "NuDiff", "OrderOfAccuracy", nuDiff.at(i - 1), nuDiff.at(i), orderOfAccuracyNuDiff.at(i - 1));
		testOut->makeTestOutput(phiDiffTestPassed.at(i - 1), testName, lx.at(i - 1), lx.at(i), "PhiDiff", "PhiDiff", "OrderOfAccuracy", phiDiff.at(i - 1), phiDiff.at(i), orderOfAccuracyPhiDiff.at(i - 1));
	}
}

void PhiAndNuTest::add(double phiDiff, double nuDiff, double lx)
{
	this->nuDiff.push_back(nuDiff);
	this->phiDiff.push_back(phiDiff);
	this->lx.push_back(lx);
}

std::string PhiAndNuTest::getOutput()
{
	std::ostringstream oss;
	oss << "#################################################" << std::endl;
	oss << "#" << std::setfill(' ') << std::right << std::setw(24 + testName.length() / 2) << testName << std::setw(24 - testName.length() / 2) << "#" << std::endl;
	oss << "#################################################" << std::endl;

	oss << "L" << "\t" << std::setfill(' ') << std::left << std::setw(15) << "NuDiff" << "Order of Accuracy" << std::endl;
	oss << lx.at(0) << "\t" << nuDiff.at(0) << std::endl;
	for (int i = 0; i < nuDiff.size() - 1; i++) {
		oss << std::setfill(' ') << std::setw(23) << " " << orderOfAccuracyNuDiff.at(i) << std::endl;
		oss << lx.at(i + 1) << "\t" << nuDiff.at(i + 1) << std::endl;
	}
	oss << std::endl;

	oss << "L" << "\t" << std::setfill(' ') << std::left << std::setw(15) << "PhiDiff" << "Order of Accuracy" << std::endl;
	oss << lx.at(0) << "\t" << phiDiff.at(0) << std::endl;
	for (int i = 0; i < phiDiff.size() - 1; i++) {
		oss << std::setfill(' ') << std::setw(23) << " " << orderOfAccuracyPhiDiff.at(i) << std::endl;
		oss << lx.at(i + 1) << "\t" << phiDiff.at(i + 1) << std::endl;
	}
	oss << std::endl;

	return oss.str();
}

PhiAndNuTest::PhiAndNuTest(std::string aTestName, double minOrderOfAccuracy, std::shared_ptr<TestCout> testOut) : testName(aTestName), minOrderOfAccuracy(minOrderOfAccuracy), testOut(testOut)
{
	phiDiff.resize(0);
	nuDiff.resize(0);
	lx.resize(0);
	orderOfAccuracyPhiDiff.resize(0);
	orderOfAccuracyNuDiff.resize(0);
	nuDiffTestPassed.resize(0);
	phiDiffTestPassed.resize(0);
}

void PhiAndNuTest::makeLastTestOutput()
{
	testOut->makeTestOutput(nuDiffTestPassed.back(), testName, lx.at(lx.size() - 1), lx.back(), "NuDiff", "NuDiff", "OrderOfAccuracy", nuDiff.at(nuDiff.size() - 1), nuDiff.back(), orderOfAccuracyNuDiff.at(orderOfAccuracyNuDiff.size() - 1));
	testOut->makeTestOutput(phiDiffTestPassed.back(), testName, lx.at(lx.size() - 1), lx.back(), "PhiDiff", "PhiDiff", "OrderOfAccuracy", phiDiff.at(phiDiff.size() - 1), phiDiff.back(), orderOfAccuracyPhiDiff.at(orderOfAccuracyPhiDiff.size() - 1));
}

std::vector<double> PhiAndNuTest::calcOrderOfAccuracy(std::vector<double> data)
{
	std::vector<double> result;
	for (int i = 1; i < lx.size(); i++) {
		double ooa = log(data.at(i - 1) / data.at(i)) / log(lx.at(i) / lx.at(i - 1));
		result.push_back(ooa);
	}
	return result;
}

std::vector<bool> PhiAndNuTest::checkTestPassed(std::vector<double> orderOfAccuracy)
{
	std::vector<bool> result;
	for (int i = 0; i < orderOfAccuracy.size(); i++) 
		result.push_back(orderOfAccuracy.at(i) > minOrderOfAccuracy);
	return result;
}
