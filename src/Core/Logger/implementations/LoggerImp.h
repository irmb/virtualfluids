//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __         
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |        
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |        
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |        
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____    
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|   
//      \    \  |    |   ________________________________________________________________    
//       \    \ |    |  |  ______________________________________________________________|   
//        \    \|    |  |  |         __          __     __     __     ______      _______    
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)   
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______    
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \   
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/   
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file LoggerImp.h
//! \ingroup Logger
//! \author Soeren Peters
//=======================================================================================
#ifndef LoggerImp_H
#define LoggerImp_H


#include <VirtualFluidsDefinitions.h>

#include <string>

#include "../Logger.h"
#include <map>

namespace logging
{

    class VF_PUBLIC LoggerImp : public Logger
    {
    public:
        LoggerImp(std::ostream* stream);
        virtual ~LoggerImp();

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
