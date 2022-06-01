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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
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
//! \author Soeren Peters
//=======================================================================================
#ifndef VF_LOGGER_H
#define VF_LOGGER_H

// VirtualFluids is using the spdlog logger https://github.com/gabime/spdlog
#include <spdlog/spdlog.h>
// To initialize spdlog initalizeLogger() must be called.
// spdlog supports 5 log level, which can be changed at runtime e.g.:
// spdlog::set_level(spdlog::level::debug)
// The default log level is set to trace. Supported levels: trace < debug < info < warning < critical
// 
// The logging is realized in 3 different log sinks:
// 1. colorded console output
// 2. a daily log file
// 3. a log file from the last run of VirtualFluids
// The default file path is relativ to executed command logs/
// File path can be changed via changeLogPath()

#define VF_LOG_TRACE(...) spdlog::trace(__VA_ARGS__)
#define VF_LOG_DEBUG(...) spdlog::debug(__VA_ARGS__)
#define VF_LOG_INFO(...) spdlog::info(__VA_ARGS__)
#define VF_LOG_WARNING(...) spdlog::warn(__VA_ARGS__)
#define VF_LOG_CRITICAL(...) spdlog::critical(__VA_ARGS__)


namespace vf::logging
{
    class Logger
    {
    public:
        // initalizing the above named logger
        static void initalizeLogger();

        // changing the path of the log files
        static void changeLogPath(const std::string& path);

    private:
        static void updateDefaultLogger();

        static std::string logPath;
    };
}

#endif
