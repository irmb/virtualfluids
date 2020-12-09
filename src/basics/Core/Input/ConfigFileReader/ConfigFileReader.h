#ifndef CONFIGFILEREADER_H
#define CONFIGFILEREADER_H


#include <memory>
#include <string>

#include "basics_export.h"

class ConfigData;

class ConfigFileReader
{
public:
    BASICS_EXPORT static std::shared_ptr<ConfigFileReader> getNewInstance();

    BASICS_EXPORT std::shared_ptr<ConfigData> readConfigFile(const char* filePath) const;

private:
    ConfigFileReader() = default;
};
#endif
