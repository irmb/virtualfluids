#ifndef CONFIGFILEREADER_H
#define CONFIGFILEREADER_H

#include "../Input.h"

#include <memory>

class ConfigData;

class  ConfigFileReader
{
public:
    BASICS_EXPORT static std::shared_ptr<ConfigFileReader> getNewInstance();
    BASICS_EXPORT virtual ~ConfigFileReader();

    BASICS_EXPORT std::shared_ptr<ConfigData> readConfigFile(const std::string &filePath) const;

private:
	ConfigFileReader();
};
#endif
