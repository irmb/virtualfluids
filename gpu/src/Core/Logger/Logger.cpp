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
//! \file Logger.cpp
//! \ingroup Logger
//! \author Soeren Peters
//=======================================================================================
#include "Logger.h"
#include <memory>
#include <iostream>

#include "implementations/LoggerImp.h"



namespace logging {

    std::shared_ptr<Logger> out = nullptr;

    logging::Logger::Level logging::Logger::globalLogLevel = logging::Logger::INFO_LOW;
    logging::Logger::Level logging::Logger::localLogLevel = logging::Logger::INFO_LOW;
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
