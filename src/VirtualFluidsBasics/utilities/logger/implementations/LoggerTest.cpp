#include "gmock/gmock.h"
#include <fstream>
#include <memory>

#include "../Logger.h"


TEST(LoggerTest, logStringWithoutSettingLevels_WillPutTheLogMesssageIntoTheStream)
{
    std::ostringstream stream;
    logging::Logger::setStream(&stream);

    *logging::out << logging::Logger::LOW << "Hello World\n";

    EXPECT_THAT(stream.str(), "[LOW]\tHello World\n");
}

TEST(LoggerTest, logStringWithHighDebugLevel_logOnlyHighLevelMessages)
{
    std::ostringstream stream;
    logging::Logger::setStream(&stream);

    logging::Logger::setDebugLevel(logging::Logger::HIGH);
    *logging::out << logging::Logger::LOW << "Low Debug Message\n" << logging::Logger::HIGH << "HIGH Debug Message\n";

    EXPECT_THAT(stream.str(), "[HIGH]\tHIGH Debug Message\n");
}

TEST(LoggerTest, addTwoStreams_shouldWriteToBoth)
{
    logging::Logger::resetStreams();

    std::ostringstream stream1, stream2;
    logging::out->addStream(&stream1);
    logging::out->addStream(&stream2);
    logging::Logger::setDebugLevel(logging::Logger::ERROR);

    *logging::out << logging::Logger::LOW <<"Hello World\n";

    EXPECT_THAT(stream1.str(), "[LOW]\tHello World\n");
    EXPECT_THAT(stream2.str(), "[LOW]\tHello World\n");
}

TEST(LoggerTest, splittetOutputShouldHaveDebugInformationOnce)
{
    std::ostringstream stream;
    logging::Logger::setStream(&stream);

    *logging::out << logging::Logger::LOW << "Hello" << " World\n";

    EXPECT_THAT(stream.str(), "[LOW]\tHello World\n");
}
