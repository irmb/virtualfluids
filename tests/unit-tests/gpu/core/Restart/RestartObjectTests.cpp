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
//! \addtogroup gpu_Restart Restart
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#include "tests/LogRedirector.h"
#include <gmock/gmock.h>

#include <gpu/core/Restart/RestartObject.h>

using namespace vf::gpu;

template <typename Type>
void saveAndLoad()
{
    std::shared_ptr<RestartObject> write_object = std::make_shared<Type>();

    write_object->fs = std::vector<std::vector<real>> {
                { 1,2,3 },
                { 4,5,6 }
            };

    const std::string name {"test_do_check_point"};
    write_object->serialize_internal(name);

    std::shared_ptr<RestartObject> read_object = std::make_shared<Type>();
    read_object->deserialize_internal(name);

    EXPECT_THAT(write_object->fs, ::testing::ContainerEq(read_object->fs));
    write_object->delete_restart_file(name);
}

TEST(RestartObjectTests, saveAndLoad_ascii)
{
    saveAndLoad<ASCIIRestartObject>();
}

TEST(RestartObjectTests, saveAndLoad_binary)
{
    saveAndLoad<BinaryRestartObject>();
}

template <typename Type>
void failToDeleteRestartFile()
{
    std::shared_ptr<RestartObject> write_object = std::make_shared<Type>();
    testing::vf::LogRedirector logRedirector;
    write_object->delete_restart_file("does_not_exist");
    EXPECT_TRUE(logRedirector.logContainsWarning()) << "recorded log: " << logRedirector.getLoggerOutput();
}

TEST(RestartObjectTests, noAsciiRestartFile_tryDeleteRestartFile_displaysWarning)
{
    failToDeleteRestartFile<ASCIIRestartObject>();
}

TEST(RestartObjectTests, noBinaryRestartFile_tryDeleteRestartFile_displaysWarning)
{
    failToDeleteRestartFile<BinaryRestartObject>();
}

//! \}
