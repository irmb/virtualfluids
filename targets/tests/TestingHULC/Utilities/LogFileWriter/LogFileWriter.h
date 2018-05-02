#ifndef LOGFILEWRITER_H
#define LOGFILEWRITER_H

#include <string>
#include <fstream>
#include <vector>
#include <memory>

class EvaluationParameter;
class DataQueue;
class TestInformation;

class LogFileWriter
{
public:
	LogFileWriter(std::vector<std::shared_ptr<EvaluationParameter>> evaPara, std::shared_ptr<TestInformation> testinfo);
	void makeDataQueueOutput(DataQueue* data, int arraySize);

private:
	std::string calcDateAndTime();
	void makeHastTags();
	void makeCenterHead(std::string output);
	
	std::fstream logFile;
	std::string logFilePath;

	time_t now;
	struct tm nowLocal;
};
#endif 
