#ifndef LoggerImp_H
#define LoggerImp_H


#include <VirtualFluidsDefinitions.h>


#include <string>
#include <memory>
#include <ostream>

#include "../Logger.h"

namespace logging
{

    class __declspec(dllexport) LoggerImp : public Logger
    {
    public:
        LoggerImp(std::ostream* stream);
        virtual ~LoggerImp();

        Logger& operator<<(const Level &level) override;
        Logger& operator<<(const std::string &message) override;
        Logger& operator<<(const int &message) override;
        Logger& operator<<(const float &message) override;
        Logger& operator<<(const double &message) override;


    private:
        std::string getRankString();
        static bool isLocalLogLevel_greateEqual_GlobalLevel();

        logging::Logger& log(const std::string &message);


    };

}


#endif
