#include "LogFileWriter.h"

#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>

#include <iomanip>
#include <ctime>

#include "Utilities\EvaluationParameter\EvaluationParameter.h"
#include "Utilities\TestInformation\TestInformation.h"
#include "Tests\DataQueue\DataQueue.h"

LogFileWriter::LogFileWriter(std::vector<std::shared_ptr<EvaluationParameter>> evaPara, std::shared_ptr<TestInformation> testinfo)
{
	std::ostringstream oss;
	oss << evaPara.at(0)->getLogFilePath() << "\\logFile_" << calcDateAndTime() << ".txt";
	this->logFilePath = oss.str();

	logFile.open(logFilePath, std::ios::out);

	makeCenterHead("LogFile Information");

	logFile << "Date: "<< std::setw(2) << std::setfill('0') << nowLocal.tm_mday << "." << std::setw(2) << nowLocal.tm_mon + 1 << "." << nowLocal.tm_year + 1900 << std::endl;
	logFile << "Time: " << std::setw(2) << std::setfill('0') << nowLocal.tm_hour << ":" << std::setw(2) << nowLocal.tm_min << ":" << std::setw(2) << nowLocal.tm_sec << std::endl;
	logFile << std::endl;


	int numberOfCudaDevices;
	cudaGetDeviceCount(&numberOfCudaDevices);
	for (int i = 0; i < numberOfCudaDevices; i++) {
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, i);
		logFile <<"GPU Device " << i + 1 << ": " << prop.name << std::endl;
	}
	logFile << std::endl;

	logFile << testinfo->getInformation();

	makeCenterHead("Test Time Information");
	logFile << "FileWriting: " << std::boolalpha << evaPara.at(0)->getWriteFiles() << std::endl;
	logFile << std::endl;
	logFile << "TestName \t \t \t" << " L\t\t" << "Time for Test" << std::endl;
	logFile << std::endl;
	for (int i = 0; i < evaPara.size(); i++)
		logFile << std::left << std::setfill(' ') << std::setw(17) << evaPara.at(i)->getTestName() << "\t" << std::right << std::setw(3) << evaPara.at(i)->getLx() << "\t\t" << std::setw(9) << evaPara.at(i)->getTestTime() << " sec" << std::endl;
	logFile << std::endl;
}

void LogFileWriter::makeDataQueueOutput(DataQueue* data, int arraySize)
{
	if (data[0].expected) {
		std::ostringstream oss;
		oss << data->testName << " " << data->valueName << " Test";
		makeCenterHead(oss.str());

		logFile << "L" << "\t" << std::setfill(' ') << std::left << std::setw(15) << data->valueName << "Order of Accuracy" << std::endl << std::endl;

		logFile << data[0].la << "\t" << data[0].a << std::endl;
		for (int i = 0; i < arraySize; i++) {
			if (data[i].expected) {
				logFile << std::setfill(' ') << std::setw(23) << " " << data[i].orderOfAccuracy << std::endl;
				logFile << data[i].lb << "\t" << data[i].b << std::endl;
			}
		}


		logFile << std::endl;
	}
}

std::string LogFileWriter::calcDateAndTime()
{
	std::ostringstream oss;
	now = time(NULL);
	nowLocal = *localtime(&now);
	oss << std::setfill('0')  << nowLocal.tm_year + 1900 << std::setw(2) << nowLocal.tm_mon + 1 << std::setw(2) << nowLocal.tm_mday << "_" << std::setw(2) << nowLocal.tm_hour << std::setw(2) << nowLocal.tm_min << std::setw(2) << nowLocal.tm_sec;
	return oss.str();
}

void LogFileWriter::makeHastTags()
{
	logFile << "#################################################" << std::endl;
}

void LogFileWriter::makeCenterHead(std::string output)
{
	makeHastTags();
	logFile << "#" << std::setfill(' ') << std::right << std::setw(24 + output.length()/2) << output << std::setw(24 - output.length() / 2) << "#" << std::endl;
	makeHastTags();
}
