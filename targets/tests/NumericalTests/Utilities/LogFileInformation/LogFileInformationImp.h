#ifndef LOGFILE_INFORMATION_IMP_H
#define LOGFILE_INFORMATION_IMP_H

#include "LogFileInformation.h"

#include <sstream>

class LogFileInformationImp : public LogFileInformation
{
public:

protected:
	void makeHashLine();
	void makeCenterHead(std::string head);

	std::ostringstream oss;
private:

};
#endif