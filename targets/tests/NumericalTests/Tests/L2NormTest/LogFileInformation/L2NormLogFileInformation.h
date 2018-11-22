#ifndef L2NORM_LOGFILE_INFORMATION_H
#define L2NORM_LOGFILE_INFORMATION_H

#include "Utilities\LogFileInformation\TestLogFileInformation\TestLogFileInformation.h"

#include <memory>
#include <vector>

class L2NormTest;

class L2NormInformation : public TestLogFileInformation
{
public:
	static std::shared_ptr< L2NormInformation> getNewInstance(std::vector< std::shared_ptr< L2NormTest>> tests, unsigned int basicTimeStep, unsigned int divergentTimeStep);

	std::string getOutput();

private:
	L2NormInformation() {};
	L2NormInformation(std::vector< std::shared_ptr< L2NormTest>> tests, unsigned int basicTimeStep, unsigned int divergentTimeStep);

	std::vector< std::shared_ptr< L2NormTest>> tests;

	unsigned int basicTimeStep, divergentTimeStep;
};
#endif