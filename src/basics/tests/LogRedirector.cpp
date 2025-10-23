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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \{
//=======================================================================================

#include "LogRedirector.h"

#include <spdlog/sinks/ostream_sink.h>

#include <logger/Logger.h>

namespace testing::vf
{
    
LogRedirector::LogRedirector()
{
    this->redirectDefaultLoggerToString();
}

LogRedirector::~LogRedirector()
{
    logger->sinks()[0] = std::move(oldSink);
    logger->set_level(oldLevel);
}

void LogRedirector::redirectDefaultLoggerToString()
{
    redirectLoggerToString(spdlog::default_logger(), spdlog::default_logger()->level());
}

void LogRedirector::redirectLoggerToString(std::shared_ptr<spdlog::logger> logger, spdlog::level::level_enum newLevel)
{
    this->logger = logger;
    oldLevel = logger->level();
    std::vector<spdlog::sink_ptr>& sinks { logger->sinks() };
    assert(sinks.size() == 1);

    oldSink = std::move(sinks[0]);

    sinks[0] = std::make_shared<spdlog::sinks::ostream_sink_st>(oss);
    logger->set_pattern("[%l] %v");
    logger->set_level(newLevel);
}

bool LogRedirector::logContainsWarning()
{
    return getLoggerOutput().find("warning") != std::string::npos;
}

std::string LogRedirector::getLoggerOutput()
{
    return oss.str();
}

} // namespace testing::vf

//! \}
