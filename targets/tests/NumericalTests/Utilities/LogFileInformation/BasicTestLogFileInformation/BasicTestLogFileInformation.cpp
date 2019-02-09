#include "BasicTestLogFileInformation.h"

std::shared_ptr<BasicTestLogFileInformation> BasicTestLogFileInformation::getInstance()
{
	static std::shared_ptr<BasicTestLogFileInformation> uniqueInstance;
	if (!uniqueInstance)
		uniqueInstance = std::shared_ptr<BasicTestLogFileInformation>(new BasicTestLogFileInformation());
	return uniqueInstance;
}

std::string BasicTestLogFileInformation::getOutput()
{
	if (!outputBuild)
		buildOutput();
	return oss.str();
}

void BasicTestLogFileInformation::addTest(std::string testName, bool testRun)
{
	bool isRegistered = false;
	for (int i = 0; i < this->testName.size(); i++) {
		if (this->testName.at(i) == testName)
			isRegistered = true;
	}
	if (!isRegistered) {
		this->testName.push_back(testName);
		this->testRun.push_back(testRun);
	}
}

BasicTestLogFileInformation::BasicTestLogFileInformation()
{
	testName.resize(0);
	testRun.resize(0);
	outputBuild = false;
}

void BasicTestLogFileInformation::buildOutput()
{
	makeCenterHead("Basic Test Information");

	for (int i = 0; i < testName.size(); i++)
		oss << testName.at(i) << "=" << std::boolalpha << testRun.at(i) << std::endl;
	oss << std::endl;

	outputBuild = true;
}
