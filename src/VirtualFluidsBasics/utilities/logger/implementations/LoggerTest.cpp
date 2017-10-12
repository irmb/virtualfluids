#include "gmock/gmock.h"
#include <fstream>
#include <memory>

#include "../Logger.h"


TEST(LoggerTest, logStringWithoutSettingLevels_WillPutTheLogMesssageIntoTheStream)
{
    std::ostringstream stream;
    logging::Logger::setStream(stream);

    *logging::out << "Hello World\n";

    EXPECT_THAT(stream.str(), "Hello World\n");
}

TEST(LoggerTest, logStringWithHighDebugLevel_logOnlyHighLevelMessages)
{
    std::ostringstream stream;
    logging::Logger::setStream(stream);

    logging::Logger::setDebugLevel(logging::Logger::HIGH);
    *logging::out << logging::Logger::LOW << "Low Debug Message\n" << logging::Logger::HIGH << "HIGH Debug Message\n";

    EXPECT_THAT(stream.str(), "HIGH Debug Message\n");
}