#ifndef FIELD_H
#define FIELD_H

#include "GridGenerator/global.h"

struct Vertex;
class GridStrategy;

class VF_PUBLIC Field : public enableSharedFromThis<Field>
{
public:
    HOST Field(SPtr<GridStrategy> gridStrategy, uint size);
    HOSTDEVICE Field();
    HOSTDEVICE ~Field();
    HOST void allocateMemory();
    HOST void freeMemory();

    HOSTDEVICE uint getSize() const;
    HOSTDEVICE char getFieldEntry(uint index) const;

    HOSTDEVICE bool is(uint index, char type) const;
    HOSTDEVICE bool isCoarseToFineNode(uint index) const;
    HOSTDEVICE bool isFineToCoarseNode(uint index) const;
	HOSTDEVICE bool isFluid(uint index) const;
	HOSTDEVICE bool isInvalidSolid(uint index) const;
	HOSTDEVICE bool isQ(uint index) const;
    HOSTDEVICE bool isBoundaryConditionNode(uint index) const;
    HOSTDEVICE bool isInvalidCoarseUnderFine(uint index) const;
    HOSTDEVICE bool isStopperOutOfGrid(uint index) const;
    HOSTDEVICE bool isStopperCoarseUnderFine(uint index) const;
	HOSTDEVICE bool isStopperSolid(uint index) const;
	HOSTDEVICE bool isStopper(uint index) const;
    HOSTDEVICE bool isInvalidOutOfGrid(uint index) const;

    HOSTDEVICE void setFieldEntry(uint index, char val);
	HOSTDEVICE void setFieldEntryToFluid(uint index);
	HOSTDEVICE void setFieldEntryToInvalidSolid(uint index);
    HOSTDEVICE void setFieldEntryToStopperOutOfGrid(uint index);
    HOSTDEVICE void setFieldEntryToStopperCoarseUnderFine(uint index);
    HOSTDEVICE void setFieldEntryToInvalidCoarseUnderFine(uint index);
    HOSTDEVICE void setFieldEntryToInvalidOutOfGrid(uint index);

private:
    SPtr<GridStrategy> gridStrategy;

    char *field;
    uint size;

    friend class GridGpuStrategy;
    friend class GridCpuStrategy;
};

#endif
