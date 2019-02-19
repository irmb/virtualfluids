#include "BasicTestLogFileInformation.h"

std::shared_ptr<BasicTestLogFileInformation> BasicTestLogFileInformation::getNewInstance()
{
	return std::shared_ptr<BasicTestLogFileInformation>(new BasicTestLogFileInformation());;
}

std::string BasicTestLogFileInformation::getOutput()
{
	if (!outputBuild) {
		buildOutput();
		outputBuild = true;
	}
	return oss.str();
}

void BasicTestLogFileInformation::addTest(std::string testName, bool testRun)
{
	this->testName.push_back(testName);
	this->testRun.push_back(testRun);
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
