#include "GridInterface.cuh"

GridInterface::GridInterface(const Grid& finerGrid)
{
    startCFCx = finerGrid.startX - finerGrid.delta * 0.5;
    startCFCy = finerGrid.startY - finerGrid.delta * 0.5;
    startCFCz = finerGrid.startZ - finerGrid.delta * 0.5;

    endCFCx = finerGrid.endX - finerGrid.delta * 1.5;
    endCFCy = finerGrid.endY - finerGrid.delta * 1.5;
    endCFCz = finerGrid.endZ - finerGrid.delta * 1.5;
}
