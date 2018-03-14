#include "Logger.h"
#include "mpi.h"
#include <memory>
#include <iostream>

#include "implementations/LoggerImp.h"



namespace logging {

    std::shared_ptr<Logger> out = nullptr;

    logging::Logger::Level logging::Logger::globalLogLevel = logging::Logger::LOW;
    logging::Logger::Level logging::Logger::localLogLevel = logging::Logger::LOW;
    bool logging::Logger::printRankNumber = false;
    bool logging::Logger::timeStampEnabled = false;
   

    logging::Logger::Logger(std::ostream* stream)
    {
        streams.push_back(stream);
    }

    logging::Logger::~Logger()
    {

    }

    void Logger::addStreamToList(std::ostream* stream)
    {
        streams.push_back(stream);
    }

    void Logger::resetStreamList()
    {
        streams.clear();
    }


//-----------static methods----------------//
    void logging::Logger::resetStreams()
    {
        if (!out)
            out = std::make_shared<LoggerImp>(&std::cout);

        out->resetStreamList();
    }

    void logging::Logger::setStream(std::ostream* stream)
    {
        out = std::make_shared<LoggerImp>(stream);
    }

    void logging::Logger::addStream(std::ostream* stream)
    {
        if (!out)
            out = std::make_shared<LoggerImp>(stream);
        else
            out->addStreamToList(stream);
    }

    void logging::Logger::timeStamp(TimeStamp timeStamp)
    {
        switch(timeStamp)
        {
        case ENABLE:
            timeStampEnabled = true;
            break;
        case DISABLE:
            timeStampEnabled = false;
            break;
        }
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
