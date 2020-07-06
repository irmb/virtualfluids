#ifndef Logger_H
#define Logger_H

#include "basics_export.h"

#include <string>
#include <memory>
#include <ostream>
#include <vector>

namespace logging 
{
    class BASICS_EXPORT Logger
    {
    protected:
        Logger(std::ostream* stream);

    public:
        virtual ~Logger();

		enum Level
		{
			INFO_LOW = 3,
			INFO_INTERMEDIATE = 2,
			INFO_HIGH = 1,
			WARNING = 0,
			LOGGER_ERROR = -1
		};

        enum TimeStamp
        {
            ENABLE,
            DISABLE
		};

        static void setStream(std::ostream* stream);
        static void addStream(std::ostream* stream);
        static void resetStreams();

        static void timeStamp(TimeStamp timeStamp);

        static void setDebugLevel(const Level &level = Level::LOGGER_ERROR);
        static void enablePrintedRankNumbers(bool printRankNumbers);

        virtual Logger& operator<<(const Level &level) = 0;
        virtual Logger& operator<<(const std::string &log) = 0;
        virtual Logger& operator<<(const int &log) = 0;
        virtual Logger& operator<<(const unsigned int &log) = 0;
        virtual Logger& operator<<(const unsigned long& log) = 0;
        virtual Logger& operator<<(const float &log) = 0;
        virtual Logger& operator<<(const double &log) = 0;

    protected:
        void addStreamToList(std::ostream* stream);
        void resetStreamList();

        std::vector<std::ostream*> streams;

        static Level globalLogLevel;
        static Level localLogLevel;
        static bool printRankNumber;
        static bool timeStampEnabled;

    };
    extern VF_SHARED_LIB_IMPORT std::shared_ptr<Logger> out;
}



#endif
