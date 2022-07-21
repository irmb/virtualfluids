#ifndef CONFIG_FILE_READER_H
#define CONFIG_FILE_READER_H

#include "Utilities/Structs/ConfigDataStruct.h"

#include <memory>
#include <string>
#include <vector>

namespace vf::basics 
{
class ConfigurationFile;
}

namespace vf::gpu::tests
{
    std::shared_ptr<ConfigDataStruct> readConfigFile(const std::string aFilePath);
}
#endif