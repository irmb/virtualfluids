#ifndef LoggerImp_H
#define LoggerImp_H


#include <VirtualFluidsDefinitions.h>


#include <string>
#include <memory>
#include <ostream>

#include "../Logger.h"

namespace logging
{

    class VF_PUBLIC LoggerImp : public Logger
    {
    public:
        LoggerImp(std::ostream &stream);
        virtual ~LoggerImp();

        Logger& operator<<(const Level &level);
        Logger& operator<<(const std::string &message);
        Logger& operator<<(const int &message);
        Logger& operator<<(const float &message);
        Logger& operator<<(const double &message);


    private:
        std::string getRankString();
        bool isLocalLogLevelHighEnough();

        logging::Logger& log(const std::string &message);


    };

}


#endif
