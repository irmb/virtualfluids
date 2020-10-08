#ifndef LoggerImp_H
#define LoggerImp_H


#include "basics_export.h"

#include <string>

#include "../Logger.h"
#include <map>

namespace logging
{

    class BASICS_EXPORT LoggerImp : public Logger
    {
    public:
        LoggerImp(std::ostream* stream);
        ~LoggerImp() override;

        Logger& operator<<(const Level &level) override;
        Logger& operator<<(const std::string &message) override;
        Logger& operator<<(const int &message) override;
        Logger& operator<<(const unsigned int &message) override;
        Logger& operator<<(const unsigned long& log) override;
        Logger& operator<<(const float &message) override;
        Logger& operator<<(const double &message) override;


    private:
        std::string getRankString();
        static bool shouldBeLogged();

        static std::string getTimeStamp();
        void addDebugInformation(std::string& message);
        logging::Logger& log(const std::string &message);

    private:
        std::map<Logger::Level, std::string> levelString;
        bool newLoggingLine = true;
    };

}


#endif
