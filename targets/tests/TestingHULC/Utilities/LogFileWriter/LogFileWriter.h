#ifndef LOG_FILE_WRITER_H
#define LOG_FILE_WRITER_H

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
	static std::shared_ptr<LogFileWriter> getNewInstance(std::string filePath);
	void makeOutput(std::string output);

private:
	LogFileWriter(std::string filePath);
	std::string calcDateAndTime();

	std::fstream logFile;
	std::string logFilePath;
	time_t now;
	struct tm nowLocal;
};
#endif 
