#ifndef POSTPROCESSING_CONFIG_FILE_READER_IMP_H
#define POSTPROCESSING_CONFIG_FILE_READER_IMP_H

#include "PostProcessingConfigFileReader.h"

class PostProcessingConfigFileReaderImp : public PostProcessingConfigFileReader
{
public:
	static std::shared_ptr<PostProcessingConfigFileReader> getNewInstance();

	std::shared_ptr<PostProcessingConfigData> readConfigFile(std::string filePath) override;

private:
	PostProcessingConfigFileReaderImp() = default;

};
#endif