#include "LogFileReader.h"

#include "Utilities\LogFileData\LogFileDataImp.h"

#include "utilities/input/Input.h"
#include "utilities/StringUtil/StringUtil.h"

#include <filesystem>
#include <fstream>

std::shared_ptr<LogFileReader> LogFileReader::getInstance()
{
	static std::shared_ptr<LogFileReader> uniqueInstance;
	if (!uniqueInstance)
		uniqueInstance = std::shared_ptr<LogFileReader>(new LogFileReader());
	return uniqueInstance;
}

std::shared_ptr<LogFileData> LogFileReader::readLogFileToLogFileData(std::string filePath)
{
	std::shared_ptr<LogFileDataImp> logFileData = LogFileDataImp::getNewInstance();

	std::ifstream stream;
	stream.open(filePath.c_str(), std::ios::in);
	if (stream.fail())
		throw "can not open config file!\n";

	std::unique_ptr<input::Input> input = input::Input::makeInput(stream, "config");

	logFileData->setBasisTimeStepLength(StringUtil::toInt(input->getValue("BasisTimeStepLength")));
	logFileData->setKernel(StringUtil::toString(input->getValue("Kernel")));
	logFileData->setNumberOfTimeSteps(StringUtil::toInt(input->getValue("NumberOfTimeSteps")));
	logFileData->setViscosity(StringUtil::toDouble(input->getValue("Viscosity")));

	std::string test = StringUtil::toString(input->getValue("GPU_Device_1"));

	return logFileData;
}

std::vector<std::shared_ptr<LogFileData>> LogFileReader::readLogFilesInDirectoryToLogFileData(std::string directory)
{
	std::vector< std::shared_ptr< LogFileData> > logFileData;

	std::vector< std::string> filePaths = getAllFilesInDir(directory, ".txt");
	for (int i = 0; i < filePaths.size(); i++)
		logFileData.push_back(readLogFileToLogFileData(filePaths.at(i)));

	return logFileData;
}

LogFileReader::LogFileReader()
{
}

std::vector<std::string> LogFileReader::getAllFilesInDir(const std::string &dirPath, const std::string &fileExtension)
{
	std::vector<std::string> listOfFiles;
	std::experimental::filesystem::path myPath = dirPath;
	if (std::experimental::filesystem::exists(myPath) && std::experimental::filesystem::is_directory(myPath))
	{
		for (auto& item : std::experimental::filesystem::recursive_directory_iterator(myPath))
		{
			if (std::experimental::filesystem::is_regular_file(item.path()) && item.path().extension() == fileExtension)
				listOfFiles.push_back(item.path().string());
		}
	}
	return listOfFiles;
}