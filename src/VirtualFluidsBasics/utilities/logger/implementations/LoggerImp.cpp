#include "LoggerImp.h"

#include "mpi.h"
#include <sstream>


logging::LoggerImp::LoggerImp(std::ostream* stream) : logging::Logger(stream)
{
    levelString[Level::WARNING] = "[WARNING]";
    levelString[Level::ERROR] = "[ERROR]";
    levelString[Level::LOW] = "[LOW]";
    levelString[Level::INTERMEDIATE] = "[INTERMEDIATE]";
    levelString[Level::HIGH] = "[HIGH]";
}

logging::LoggerImp::~LoggerImp()
{

}

logging::Logger& logging::LoggerImp::operator<<(const Level &level)
{
    localLogLevel = level;
    return *this;
}


logging::Logger& logging::LoggerImp::operator<<(const std::string &message)
{
    return this->log(message);
}

logging::Logger& logging::LoggerImp::operator<<(const int &message)
{
    return this->log(std::to_string(message));
}

logging::Logger& logging::LoggerImp::operator<<(const float &message)
{
    return this->log(std::to_string(message));
}

logging::Logger& logging::LoggerImp::operator<<(const double &message)
{
    return this->log(std::to_string(message));
}

logging::Logger& logging::LoggerImp::log(const std::string &message)
{
    if (isLocalLogLevel_greateEqual_GlobalLevel())
    {
        std::string modifiedMessage = message;
        addDebugInformation(modifiedMessage);
        for(auto stream : streams)
            *stream << modifiedMessage;
    }
    return *this;
}

bool logging::LoggerImp::isLocalLogLevel_greateEqual_GlobalLevel()
{
    return localLogLevel >= globalLogLevel;
}

void logging::LoggerImp::addDebugInformation(std::string& message)
{
    std::stringstream os;
    os << levelString[this->localLogLevel] << "\t" << message;
    message = os.str();
}

std::string logging::LoggerImp::getRankString()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return printRankNumber ? "[" + std::to_string(rank) + "] " : "";
}

