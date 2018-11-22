#ifndef CONFIG_FILE_READER_H
#define CONFIG_FILE_READER_H

#include "ConfigData.h"

#include <memory>
#include <string>

class ConfigFileReader
{
public:
	static std::shared_ptr< ConfigFileReader> getNewInstance(const std::string aFilePath);
	std::shared_ptr< ConfigDataStruct> getConfigData();
	
private:
	ConfigFileReader() {};
	ConfigFileReader(const std::string aFilePath);
	void readConfigFile(const std::string aFilePath);
	void checkConfigFileData();
	std::shared_ptr< ConfigDataStruct> configData;
};
#endif