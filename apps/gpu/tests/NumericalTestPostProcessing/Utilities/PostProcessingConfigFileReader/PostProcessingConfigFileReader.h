#ifndef POSTPROCESSING_CONFIG_FILE_READER_H
#define POSTPROCESSING_CONFIG_FILE_READER_H

#include <memory>
#include <string>

class PostProcessingConfigData;

class PostProcessingConfigFileReader
{
public:
	virtual std::shared_ptr<PostProcessingConfigData> readConfigFile(std::string filePath) = 0;
};
#endif