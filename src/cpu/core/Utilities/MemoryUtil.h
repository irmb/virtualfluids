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
//! \addtogroup cpu_Utilities Utilities
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef _MEMORYUTIL_H_
#define _MEMORYUTIL_H_

#if defined(_WIN32) || defined(_WIN64)
#define MEMORYUTIL_WINDOWS
#include "windows.h"
#include "psapi.h"
#pragma comment(lib, "psapi.lib")
#elif defined __APPLE__
#define MEMORYUTIL_APPLE
#include "sys/sysctl.h"
#include "sys/types.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#elif (defined(__amd64) || defined(__amd64__) || defined(__unix__) || defined(__CYGWIN__)) && !defined(__AIX__)
#define MEMORYUTIL_LINUX
#include "cstdio"
#include "cstdlib"
#include "cstring"
#include "sys/sysinfo.h"
#include "sys/types.h"
#else
#error "MemoryUtil::UnknownMachine"
#endif

#if defined(__CYGWIN__)
#define MEMORYUTIL_CYGWIN
#endif

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "Grid3D.h"
#include "lbm/constants/D3Q27.h"

//////////////////////////////////////////////////////////////////////////
// MemoryUtil
//////////////////////////////////////////////////////////////////////////
namespace Utilities
{
//////////////////////////////////////////////////////////////////////////
static long long getTotalPhysMem()
{
#if defined(MEMORYUTIL_WINDOWS)
    MEMORYSTATUSEX memInfo;
    memInfo.dwLength = sizeof(MEMORYSTATUSEX);
    GlobalMemoryStatusEx(&memInfo);
    DWORDLONG totalPhysMem = memInfo.ullTotalPhys;
#elif defined(MEMORYUTIL_LINUX)
    struct sysinfo memInfo;
    sysinfo(&memInfo);
    long long totalPhysMem = memInfo.totalram;
    // Multiply in next statement to avoid int overflow on right hand side...
    totalPhysMem *= memInfo.mem_unit;
#elif defined(MEMORYUTIL_APPLE)
    int mib[] = { CTL_HW, HW_MEMSIZE };
    int64_t totalPhysMem;
    size_t length = sizeof(totalPhysMem);

    if (sysctl(mib, 2, &totalPhysMem, &length, nullptr, 0) == -1)
        return 0;
#else
#error "MemoryUtil::getTotalPhysMem - UnknownMachine"
#endif

    return (long long)totalPhysMem;
}
//////////////////////////////////////////////////////////////////////////
static long long getPhysMemUsed()
{
#if defined(MEMORYUTIL_WINDOWS)
    MEMORYSTATUSEX memInfo;
    memInfo.dwLength = sizeof(MEMORYSTATUSEX);
    GlobalMemoryStatusEx(&memInfo);
    DWORDLONG physMemUsed = memInfo.ullTotalPhys - memInfo.ullAvailPhys;
#elif defined(MEMORYUTIL_LINUX)
    struct sysinfo memInfo;
    sysinfo(&memInfo);
    long long physMemUsed = memInfo.totalram - memInfo.freeram;
    // Multiply in next statement to avoid int overflow on right hand side...
    physMemUsed *= memInfo.mem_unit;
#elif defined(MEMORYUTIL_APPLE)
    long long physMemUsed = 0;
#else
#error "MemoryUtil::getPhysMemUsed - UnknownMachine"
#endif

    return (long long)physMemUsed;
}
//////////////////////////////////////////////////////////////////////////
#if defined(MEMORYUTIL_LINUX) || defined(MEMORYUTIL_APPLE) || defined(MEMORYUTIL_CYGWIN)
static int parseLine(char *line)
{
    int i = strlen(line);
    while (*line < '0' || *line > '9')
        line++;
    line[i - 3] = '\0';
    i           = atoi(line);
    return i;
}

static int getValue()
{ // Note: this value is in KB!
    FILE *file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL) {
        if (strncmp(line, "VmRSS:", 6) == 0) {
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}
#endif
//////////////////////////////////////////////////////////////////////////
static long long getPhysMemUsedByMe()
{
#if defined(MEMORYUTIL_WINDOWS) && !defined(__CYGWIN__)
    PROCESS_MEMORY_COUNTERS pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
    SIZE_T physMemUsedByMe = pmc.WorkingSetSize;
#elif defined(MEMORYUTIL_LINUX) || defined(MEMORYUTIL_APPLE) || defined(MEMORYUTIL_CYGWIN)
    long long physMemUsedByMe = (long long)getValue() * (long long)1024;
#else
#error "MemoryUtil::getPhysMemUsedByMe - UnknownMachine"
#endif

    return (long long)physMemUsedByMe;
}
//////////////////////////////////////////////////////////////////////////

static std::string toString(SPtr<Grid3D> grid, int numberOfProcesses)
{
    unsigned long long numberOfBlocks = (unsigned long long)grid->getNumberOfBlocks();
    int ghostLayer = grid->getGhostLayerWidth()*2+1;
    UbTupleInt3 blockNx = grid->getBlockNX();

    unsigned long long numberOfNodesPerBlock = (unsigned long long)(val<1>(blockNx)) *
                                               (unsigned long long)(val<2>(blockNx)) *
                                               (unsigned long long)(val<3>(blockNx));
    unsigned long long numberOfNodes = numberOfBlocks * numberOfNodesPerBlock;
    unsigned long long numberOfNodesPerBlockWithGhostLayer = numberOfBlocks * (val<1>(blockNx) + ghostLayer) *
                                                             (val<2>(blockNx) + ghostLayer) *
                                                             (val<3>(blockNx) + ghostLayer);
    real needMemAll = real(numberOfNodesPerBlockWithGhostLayer*(27*sizeof(real)+sizeof(int)+sizeof(real)*4));
    real needMem = needMemAll / real(numberOfProcesses);
    
    std::ostringstream out;
    out << "Grid information:" << std::endl;
    out << "###################################################" << std::endl;
    out << "# Number of blocks = " << numberOfBlocks << std::endl;
    out << "# Number of nodes  = " << numberOfNodes << std::endl;
    int minInitLevel = grid->getCoarsestInitializedLevel();
    int maxInitLevel = grid->getFinestInitializedLevel();
    for (int level = minInitLevel; level<=maxInitLevel; level++)
    {
        int nobl = grid->getNumberOfBlocks(level);
        out << "# Number of blocks for level " << level << " = " << nobl << std::endl;
        out << "# Number of nodes for level " << level << " = " << nobl * numberOfNodesPerBlock << std::endl;
    }
    out << "# Necessary memory  = " << needMemAll << " bytes" << std::endl;
    out << "# Necessary memory per process = " << needMem << " bytes" << std::endl;
    out << "# Available memory per process = " << (real)getTotalPhysMem() << " bytes" << std::endl;
    out << "###################################################" << std::endl;

    return out.str();
}

} // namespace Utilities

#endif

//! \}
