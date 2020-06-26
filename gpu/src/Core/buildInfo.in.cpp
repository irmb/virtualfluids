#include "VirtualFluidsDefinitions.h"

namespace buildInfo
{
    VF_PUBLIC const char *gitCommitHash() { return "@git_commit_hash@";  }
    VF_PUBLIC const char *gitBranch()     { return "@git_branch@";       }
    VF_PUBLIC const char *buildType()     { return "@CMAKE_BUILD_TYPE@"; }
    VF_PUBLIC const char *compilerFlags() { return "@COMPILER_FLAGS@";   }
    VF_PUBLIC const char *buildMachine()  { return "@BUILD_computerName@";}
    VF_PUBLIC const char *projectDir()    { return "@CMAKE_SOURCE_DIR@"; }
    VF_PUBLIC const char *binaryDir()     { return "@CMAKE_BINARY_DIR@"; }
}
