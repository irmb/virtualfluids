#ifndef buildInfo_H
#define buildInfo_H


namespace buildInfo
{
    const char * gitCommitHash();
    const char * gitBranch();
    const char * buildType();
    const char * compilerFlags();
    const char * buildMachine();
    const char * projectDir();
    const char * binaryDir();
}

#define GIT_COMMIT_HASH   buildinfo::gitCommitHash()
#define GIT_BRANCH        buildinfo::gitBranch()
#define BUILD_MACHINE     buildinfo::buildMachine()
#define PROJECT_DIR       buildinfo::projectDir()
#define BINARY_DIR        buildinfo::binaryDir()
#define COMPILER_FLAGS    buildinfo::compilerFlags()
#define BUILD_TYPE        buildinfo::buildType()



#endif
