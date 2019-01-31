#ifndef LOGFILE_HEAD_H
#define LOGFILE_HEAD_H

#include "../LogFileInformationImp.h"

#include <memory>
#include <vector>

class LogFileHead : public LogFileInformationImp
{
public:
	static std::shared_ptr<LogFileHead> getNewInstance(std::vector<int> devices);
	std::string getOutput();

private:
	void calcDateAndTime();
	LogFileHead() {};
	LogFileHead(std::vector<int> devices);

	std::vector<int> devices;
	time_t now;
	struct tm nowLocal;
};
#endif 