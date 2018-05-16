#ifndef LOGFILE_INFORMATION_OUTPUT_H
#define LOGFILE_INFORMATION_OUTPUT_H

#include "../LogFileInformationImp.h"

#include <memory>
#include <vector>

class LogFileInformationOutput : public LogFileInformationImp
{
public:
	static std::shared_ptr<LogFileInformation> getNewInstance(std::vector<int> devices);
	std::string getOutput();

private:
	void calcDateAndTime();
	LogFileInformationOutput() {};
	LogFileInformationOutput(std::vector<int> devices);

	std::vector<int> devices;
	time_t now;
	struct tm nowLocal;
};
#endif 