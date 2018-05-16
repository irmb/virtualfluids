#include "LogFileInformationImp.h"

#include <iomanip>

void LogFileInformationImp::makeHashLine()
{
	oss << "#################################################" << std::endl;
}

void LogFileInformationImp::makeCenterHead(std::string head)
{
	makeHashLine();
	oss << "#" << std::setfill(' ') << std::right << std::setw(24 + head.length() / 2) << head << std::setw(24 - head.length() / 2) << "#" << std::endl;
	makeHashLine();
}
