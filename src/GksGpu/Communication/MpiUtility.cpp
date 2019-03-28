#include "MpiUtility.h"

#include <exception>

#define ENV_LOCAL_RANK	     "OMPI_COMM_WORLD_LOCAL_RANK"
#define ENV_COMM_WORLD_SIZE  "OMPI_COMM_WORLD_SIZE"

int MpiUtility::getMpiRankBeforeInit()
{
    char * localRankStr = NULL;

    // We extract the local rank initialization using an environment variable
    if ((localRankStr = getenv(ENV_LOCAL_RANK)) != NULL)
    {
        return atoi(localRankStr);
    }
    else
    {
        return 0;
    }
}

int MpiUtility::getMpiWorldSizeBeforeInit()
{
    char * mpiWorldSizeStr = NULL;

    // We extract the local rank initialization using an environment variable
    if ((mpiWorldSizeStr = getenv(ENV_COMM_WORLD_SIZE)) != NULL)
    {
        return atoi(mpiWorldSizeStr);
    }
    else
    {
        return 1;
    }
}
