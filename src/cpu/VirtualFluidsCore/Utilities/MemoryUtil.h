#ifndef _MEMORYUTIL_H_
#define _MEMORYUTIL_H_

#if defined(_WIN32) || defined(_WIN64)
   #define MEMORYUTIL_WINDOWS
   #include "windows.h"
   #include "psapi.h"
   #pragma comment(lib, "psapi.lib")
#elif defined __APPLE__
#define MEMORYUTIL_APPLE
   #include "sys/types.h"
   #include "sys/sysctl.h"
   #include "stdlib.h"
   #include "stdio.h"
   #include "string.h"
#elif (defined(__amd64) || defined(__amd64__) || defined(__unix__) || defined(__CYGWIN__)) && !defined(__AIX__) 
   #define MEMORYUTIL_LINUX
   #include "sys/types.h"
   #include "sys/sysinfo.h"
   #include "stdlib.h"
   #include "stdio.h"
   #include "string.h"
#else
   #error "MemoryUtil::UnknownMachine"
#endif
//////////////////////////////////////////////////////////////////////////
//MemoryUtil
//////////////////////////////////////////////////////////////////////////
namespace Utilities
{
//////////////////////////////////////////////////////////////////////////
   static long long getTotalPhysMem()
   {
      #if defined MEMORYUTIL_WINDOWS
         MEMORYSTATUSEX memInfo;
         memInfo.dwLength = sizeof(MEMORYSTATUSEX);
         GlobalMemoryStatusEx(&memInfo);
         DWORDLONG totalPhysMem = memInfo.ullTotalPhys;            
      #elif defined(MEMORYUTIL_LINUX)
         struct sysinfo memInfo;
         sysinfo (&memInfo);
         long long totalPhysMem = memInfo.totalram;
         //Multiply in next statement to avoid int overflow on right hand side...
         totalPhysMem *= memInfo.mem_unit;
      #elif defined(MEMORYUTIL_APPLE)
        int mib [] = { CTL_HW, HW_MEMSIZE };
        int64_t totalPhysMem;
        size_t length = sizeof(totalPhysMem);

        if(sysctl(mib, 2, &totalPhysMem, &length, nullptr, 0) == -1)
            return 0;
      #else
      #error "MemoryUtil::getTotalPhysMem - UnknownMachine"
      #endif

      return (long long)totalPhysMem;
   }
//////////////////////////////////////////////////////////////////////////
   static long long getPhysMemUsed()
   {
      #if defined MEMORYUTIL_WINDOWS
         MEMORYSTATUSEX memInfo;
         memInfo.dwLength = sizeof(MEMORYSTATUSEX);
         GlobalMemoryStatusEx(&memInfo);
         DWORDLONG physMemUsed = memInfo.ullTotalPhys - memInfo.ullAvailPhys;          
      #elif defined(MEMORYUTIL_LINUX)
         struct sysinfo memInfo;
         sysinfo (&memInfo);
         long long physMemUsed = memInfo.totalram - memInfo.freeram;
         //Multiply in next statement to avoid int overflow on right hand side...
         physMemUsed *= memInfo.mem_unit;
         #elif defined(MEMORYUTIL_APPLE)
         long long physMemUsed = 0;
      #else
      #error "MemoryUtil::getPhysMemUsed - UnknownMachine"
      #endif

      return (long long)physMemUsed;
   }
//////////////////////////////////////////////////////////////////////////
#if defined(MEMORYUTIL_LINUX) || defined(MEMORYUTIL_APPLE)
   static int parseLine(char* line){
      int i = strlen(line);
      while (*line < '0' || *line > '9') line++;
      line[i-3] = '\0';
      i = atoi(line);
      return i;
   }

   static int getValue(){ //Note: this value is in KB!
      FILE* file = fopen("/proc/self/status", "r");
      int result = -1;
      char line[128];


      while (fgets(line, 128, file) != NULL){
         if (strncmp(line, "VmRSS:", 6) == 0){
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
      #if defined MEMORYUTIL_WINDOWS
         PROCESS_MEMORY_COUNTERS pmc;
         GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
         SIZE_T physMemUsedByMe = pmc.WorkingSetSize;          
      #elif defined(MEMORYUTIL_LINUX) || defined(MEMORYUTIL_APPLE)
         long long physMemUsedByMe = (long long)getValue() * (long long)1024;
      #else
         #error "MemoryUtil::getPhysMemUsedByMe - UnknownMachine"
      #endif

      return (long long)physMemUsedByMe;
   }
//////////////////////////////////////////////////////////////////////////

}

#endif

