#ifndef FIELD_H
#define FIELD_H

#include "GridGenerator/global.h"

struct Vertex;
class GridStrategy;

class VF_PUBLIC Field : public enableSharedFromThis<Field>
{
public:
    HOST Field(SPtr<GridStrategy> gridStrategy, uint size);
    Field();
    ~Field();
    HOST void allocateMemory();
    HOST void freeMemory();

    HOSTDEVICE uint getSize() const;
    HOSTDEVICE char getFieldEntry(uint index) const;

    HOSTDEVICE bool is(uint index, char type) const;
    HOSTDEVICE bool isCoarseToFineNode(uint index) const;
    HOSTDEVICE bool isFineToCoarseNode(uint index) const;
	HOSTDEVICE bool isFluid(uint index) const;
	HOSTDEVICE bool isSolid(uint index) const;
	HOSTDEVICE bool isQ(uint index) const;
    HOSTDEVICE bool isRb(uint index) const;
    HOSTDEVICE bool isInvalid(uint index) const;
    HOSTDEVICE bool isStopperEndOfGrid(uint index) const;
    HOSTDEVICE bool isStopperOverlapGrid(uint index) const;
    HOSTDEVICE bool isOutOfGrid(uint index) const;

    HOSTDEVICE void setFieldEntry(uint index, char val);
	HOSTDEVICE void setFieldEntryToFluid(uint index);
	HOSTDEVICE void setFieldEntryToSolid(uint index);
    HOSTDEVICE void setFieldEntryToStopperEndOfGrid(uint index);
    HOSTDEVICE void setFieldEntryToStopperOverlapGrid(uint index);
    HOSTDEVICE void setFieldEntryToInvalid(uint index);
    HOSTDEVICE void setFieldEntryToOutOfGrid(uint index);

private:
    SPtr<GridStrategy> gridStrategy;

    char *field;
    uint size;

    friend class GridGpuStrategy;
    friend class GridCpuStrategy;
};

#endif
