#ifndef MULTIPLE_GRID_BUILDER_H
#define MULTIPLE_GRID_BUILDER_H

#include "GridGenerator/global.h"

#include <vector>
#include <string>
#include <memory>

template<typename Grid>
class MultipleGridBuilder
{
private:
    VF_PUBLIC MultipleGridBuilder();

public:
    VF_PUBLIC static SPtr<MultipleGridBuilder> makeShared();

    VF_PUBLIC void addGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, bool periodicityX, bool periodicityY, bool periodicityZ);
    VF_PUBLIC uint getNumberOfLevels() const;

    
};

#endif

