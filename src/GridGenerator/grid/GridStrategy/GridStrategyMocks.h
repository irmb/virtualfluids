#ifndef GRID_STRATEGYMOCKS_H
#define GRID_STRATEGYMOCKS_H

#include <VirtualFluidsDefinitions.h>
#include "core/PointerDefinitions.h"

#include "GridStrategy.h"

struct Geometry;
struct Grid;

class VF_PUBLIC GridStrategyDummy : public GridStrategy
{
public:
    virtual ~GridStrategyDummy() {};

    virtual void allocateGridMemory(SPtr<Grid> grid) override {};

    virtual void initalNodes(SPtr<Grid> grid) override {};
    virtual void mesh(SPtr<Grid> grid, Geometry &geom) override {};

    virtual void createGridInterface(SPtr<Grid> grid, SPtr<Grid> finerGrid) override {};

    virtual void deleteSolidNodes(SPtr<Grid> grid) override {};

    virtual void freeMemory(SPtr<Grid> grid) override {};

};

#endif
