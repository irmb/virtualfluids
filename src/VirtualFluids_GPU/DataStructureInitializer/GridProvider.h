#ifndef GridReader_H
#define GridReader_H

#include <string>
#include <vector>
#include <memory>

#include "VirtualFluids_GPU_EXPORT.h"

class Parameter;
class GridBuilder;
class CudaMemoryManager;

class VirtualFluids_GPU_EXPORT GridProvider
{
public:
    static std::shared_ptr<GridProvider> makeGridGenerator(std::shared_ptr<GridBuilder> builder, std::shared_ptr<Parameter> para);
    static std::shared_ptr<GridProvider> makeGridReader(bool readBinaryFiles, std::shared_ptr<Parameter> para);

	virtual void allocArrays_CoordNeighborGeo() = 0;
	virtual void allocArrays_BoundaryValues() = 0;
	virtual void allocArrays_BoundaryQs() = 0;
	virtual void setDimensions() = 0;
	virtual void setBoundingBox() = 0;
	virtual void initPeriodicNeigh(std::vector<std::vector<std::vector<unsigned int> > > periodV, std::vector<std::vector<unsigned int> > periodIndex, std::string way) = 0;

    virtual void allocAndCopyForcing();
    virtual void freeMemoryOnHost();
    virtual void cudaCopyDataToHost(int level);

	virtual ~GridProvider() {};

protected:
	void setNumberOfNodes(const int numberOfNodes, const int level) const;
	void setInitalNodeValues(const int numberOfNodes, const int level) const;

	void setPressSizePerLevel(int level, int sizePerLevel) const;
	void setVelocitySizePerLevel(int level, int sizePerLevel) const;
	void setOutflowSizePerLevel(int level, int sizePerLevel) const;

    std::shared_ptr<Parameter> para;
    std::shared_ptr<CudaMemoryManager> cudaMemoryManager;
};

#endif
