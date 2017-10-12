#include "Logger.h"
#include "mpi.h"
#include <memory>
#include <iostream>

#include "implementations/LoggerImp.h"



namespace logging {

    std::unique_ptr<Logger> out = std::unique_ptr<Logger>(new LoggerImp(std::cout));;

    logging::Logger::Level logging::Logger::globalLogLevel = logging::Logger::INTERMEDIATE;
    logging::Logger::Level logging::Logger::localLogLevel = logging::Logger::INTERMEDIATE;
    bool logging::Logger::printRankNumber = false;
   

    logging::Logger::Logger(std::ostream &stream) : stream(stream)
    {
        
    }

    logging::Logger::~Logger()
    {

    }

    void logging::Logger::setStream(std::ostream &stream)
    {
        out = std::unique_ptr<Logger>(new LoggerImp(stream));
    }

    void logging::Logger::setDebugLevel(const Level &level)
    {
        globalLogLevel = level;
    }

    void logging::Logger::enablePrintedRankNumbers(bool print)
    {
        printRankNumber = print;
    }

}
