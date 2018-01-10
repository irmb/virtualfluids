#include "Logger.h"
#include "mpi.h"
#include <memory>
#include <iostream>

#include "implementations/LoggerImp.h"



namespace logging {

    std::shared_ptr<Logger> out = std::make_shared<LoggerImp>(std::cout);

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
        out = std::make_shared<LoggerImp>(stream);
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
