#include "buildInfo.h"

#include "basics_export.h"


namespace buildInfo
{
    BASICS_EXPORT const char *gitCommitHash() { return "@git_commit_hash@";  }
    BASICS_EXPORT const char *gitBranch()     { return "@git_branch@";       }
    BASICS_EXPORT const char *buildType()     { return "@CMAKE_BUILD_TYPE@"; }
    BASICS_EXPORT const char *compilerFlags() { return "@COMPILER_FLAGS@";   }
    BASICS_EXPORT const char *buildMachine()  { return "@BUILD_computerName@";}
    BASICS_EXPORT const char *projectDir()    { return "@CMAKE_SOURCE_DIR@"; }
    BASICS_EXPORT const char *binaryDir()     { return "@CMAKE_BINARY_DIR@"; }
}
