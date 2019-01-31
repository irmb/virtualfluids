#ifndef LOGFILE_TIME_INFORMATION_H
#define LOGFILE_TIME_INFORMATION_H

#include "../LogFileInformationImp.h"

#include <memory>
#include <vector>

class TestSimulation;

class LogFileTimeInformation : public LogFileInformationImp
{
public:
	static std::shared_ptr<LogFileTimeInformation> getNewInstance(std::vector<std::shared_ptr<TestSimulation> > testSimulation, bool fileWriting);
	std::string getOutput();

private:
	LogFileTimeInformation();
	LogFileTimeInformation(std::vector<std::shared_ptr<TestSimulation> > testSimulation, bool fileWriting);

	std::vector<std::shared_ptr< TestSimulation> > testSimulation;
	bool fileWriting;
};
#endif