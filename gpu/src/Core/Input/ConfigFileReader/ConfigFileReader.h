#ifndef CONFIGFILEREADER_H
#define CONFIGFILEREADER_H

#include "../Input.h"

#include <memory>

class ConfigData;

class  ConfigFileReader
{
public:
    VF_PUBLIC static std::shared_ptr<ConfigFileReader> getNewInstance();
    VF_PUBLIC virtual ~ConfigFileReader(void);

    VF_PUBLIC std::shared_ptr<ConfigData> readConfigFile(const std::string &filePath) const;

private:
	ConfigFileReader();
};
#endif
