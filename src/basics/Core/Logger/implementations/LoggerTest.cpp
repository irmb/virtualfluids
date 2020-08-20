#include "gmock/gmock.h"
#include <fstream>
#include <memory>

#include "../Logger.h"


TEST(DISABLED_LoggerTest, logStringWithoutSettingLevels_WillPutTheLogMesssageIntoTheStream)
{
    std::ostringstream stream;
    logging::Logger::setStream(&stream);

    *logging::out << logging::Logger::INFO_LOW << "Hello World\n";

    EXPECT_THAT(stream.str(), "[LOW] Hello World\n");
}

TEST(DISABLED_LoggerTest, logStringWithHighDebugLevel_logOnlyHighLevelMessages)
{
    std::ostringstream stream;
    logging::Logger::setStream(&stream);

    logging::Logger::setDebugLevel(logging::Logger::INFO_HIGH);
    *logging::out << logging::Logger::INFO_LOW << "Low Debug Message\n" << logging::Logger::INFO_HIGH << "HIGH Debug Message\n";

    EXPECT_THAT(stream.str(), "[HIGH] HIGH Debug Message\n");
}

TEST(DISABLED_LoggerTest, addTwoStreams_shouldWriteToBoth)
{
    logging::Logger::resetStreams();

    std::ostringstream stream1, stream2;
    logging::out->addStream(&stream1);
    logging::out->addStream(&stream2);
    logging::Logger::setDebugLevel(logging::Logger::INFO_LOW);

    *logging::out << logging::Logger::INFO_LOW <<"Hello World\n";

    EXPECT_THAT(stream1.str(), "[LOW] Hello World\n");
    EXPECT_THAT(stream2.str(), "[LOW] Hello World\n");
}

TEST(DISABLED_LoggerTest, splittetOutputShouldHaveDebugInformationOnce)
{
    std::ostringstream stream;
    logging::Logger::setStream(&stream);

    *logging::out << logging::Logger::INFO_LOW << "Hello" << " World\n";

    EXPECT_THAT(stream.str(), "[LOW] Hello World\n");
}

TEST(DISABLED_LoggerTest, enableTimeStampInOutput)
{
    std::ostringstream stream;
    logging::Logger::setStream(&stream);
    logging::Logger::timeStamp(logging::Logger::TimeStamp::ENABLE);

    *logging::out << logging::Logger::INFO_LOW << "Hello" << " World\n";
    
    EXPECT_THAT(stream.str(), testing::StrNe("[LOW] Hello World\n"));
}

