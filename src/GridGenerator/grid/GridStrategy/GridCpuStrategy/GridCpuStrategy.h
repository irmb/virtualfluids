#ifndef GRID_CPU_STRATEGY_H
#define GRID_CPU_STRATEGY_H

#include "GridGenerator/global.h"

#include "../GridStrategy.h"

#include "core/PointerDefinitions.h"


struct Grid;
struct Geometry;

class VF_PUBLIC GridCpuStrategy : public GridStrategy
{
public:
	virtual ~GridCpuStrategy() {};

    void allocateGridMemory(SPtr<Grid> grid) override;

    void initalNodes(SPtr<Grid> grid) override;
    void mesh(SPtr<Grid> grid, Geometry &geom) override;

    void freeMemory(SPtr<Grid> grid) override;


    void deleteSolidNodes(SPtr<Grid> grid) override;

    virtual void copyDataFromGPU() {};

protected:
    void findInvalidNodes(SPtr<Grid> grid);
    void findNeighborIndices(SPtr<Grid> grid);

};

#endif
