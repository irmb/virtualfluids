#include "LogFileInformationOutput.h"

#include <iomanip>
#include <ctime>

#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>

std::shared_ptr<LogFileInformation> LogFileInformationOutput::getNewInstance(std::vector<int> devices)
{
	return std::shared_ptr<LogFileInformation>(new LogFileInformationOutput(devices));
}

std::string LogFileInformationOutput::getOutput()
{
	calcDateAndTime();

	makeCenterHead("LogFile Information");
	oss << "Date: " << std::setw(2) << std::setfill('0') << nowLocal.tm_mday << "." << std::setw(2) << nowLocal.tm_mon + 1 << "." << nowLocal.tm_year + 1900 << std::endl;
	oss << "Time: " << std::setw(2) << std::setfill('0') << nowLocal.tm_hour << ":" << std::setw(2) << nowLocal.tm_min << ":" << std::setw(2) << nowLocal.tm_sec << std::endl;
	oss << std::endl;

	for (int i = 0; i < devices.size(); i++) {
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, devices.at(i));
		oss << "GPU Device " << devices.at(i) << ": " << prop.name << std::endl;
	}
	oss << std::endl;

	return oss.str();
}

void LogFileInformationOutput::calcDateAndTime()
{
	now = time(NULL);
	nowLocal = *localtime(&now);
}

LogFileInformationOutput::LogFileInformationOutput(std::vector<int> devices) : devices(devices)
{

}
