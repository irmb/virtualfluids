#include "MultipleGridBuilder.h"

#include "../GridMocks.h"

template<typename Grid>
MultipleGridBuilder<Grid>::MultipleGridBuilder()
{
}

template<typename Grid>
SPtr<MultipleGridBuilder<Grid> > MultipleGridBuilder<Grid>::makeShared()
{
    return SPtr<MultipleGridBuilder>(new MultipleGridBuilder<Grid>());
}

template <typename Grid>
void MultipleGridBuilder<Grid>::addGrid(real startX, real startY, real startZ, real endX, real endY, real endZ,
    bool periodicityX, bool periodicityY, bool periodicityZ)
{

}

template <typename Grid>
uint MultipleGridBuilder<Grid>::getNumberOfLevels() const
{
    return 0;
}

template class MultipleGridBuilder<GridDummy>;
