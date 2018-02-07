#ifndef GRID_CPU_STRATEGY_H
#define GRID_CPU_STRATEGY_H

#include "GridGenerator/global.h"

#include "../GridStrategy.h"

#include "core/PointerDefinitions.h"


class Grid;
struct Geometry;

class VF_PUBLIC GridCpuStrategy : public GridStrategy
{
public:
	virtual ~GridCpuStrategy() {};

    void allocateGridMemory(SPtr<Grid> grid) override;

    void initalNodes(SPtr<Grid> grid) override;
    void mesh(SPtr<Grid> grid, Geometry &geom) override;

    void createGridInterface(SPtr<Grid> grid, SPtr<Grid> finerGrid) override;

    void freeMemory(SPtr<Grid> grid) override;


    void deleteSolidNodes(SPtr<Grid> grid) override;

    virtual void copyDataFromGPU() {};

protected:
    static void findInvalidNodes(SPtr<Grid> grid);
    static void findForNeighborsNewIndices(SPtr<Grid> grid);
    static void findForGridInterfaceNewIndices(SPtr<Grid> grid);

    static void findGridInterface(SPtr<Grid> grid, SPtr<Grid> finerGrid);

};

#endif
