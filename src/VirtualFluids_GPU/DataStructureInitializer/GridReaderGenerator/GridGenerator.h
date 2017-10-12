#ifndef GridReaderGenerator_H
#define GridReaderGenerator_H

#include "../GridProvider.h"
#include "VirtualFluids_GPU_EXPORT.h"
#include <vector>
#include <string>
#include <memory>

#include "LBM/LB.h"

class Parameter;
class GridBuilder;

class GridGenerator
	: public GridProvider
{
private:
	std::vector<std::string> channelDirections;
	std::vector<std::string> channelBoundaryConditions;

	std::shared_ptr<GridBuilder> builder;

public:
	VirtualFluids_GPU_EXPORT GridGenerator(std::shared_ptr<GridBuilder> builder, std::shared_ptr<Parameter> para);
	VirtualFluids_GPU_EXPORT virtual ~GridGenerator();

	void setUnstructuredGridBuilder(std::shared_ptr<GridBuilder> builder);

	virtual void allocArrays_CoordNeighborGeo();
	virtual void allocArrays_BoundaryValues();
	virtual void allocArrays_BoundaryQs();

	virtual void setDimensions();
	virtual void setBoundingBox();

	virtual void initPeriodicNeigh(std::vector<std::vector<std::vector<unsigned int> > > periodV, std::vector<std::vector<unsigned int> > periodIndex, std::string way);
	
private:
	void setPressureValues(int channelSide) const;
	void setPressRhoBC(int sizePerLevel, int level, int channelSide) const;

	void setVelocityValues(int channelSide) const;
	void setVelocity(int level, int sizePerLevel, int channelSide) const;

	void setOutflowValues(int channelSide) const;
	void setOutflow(int level, int sizePerLevel, int channelSide) const;

	void setPressQs(int channelSide) const;
	void setVelocityQs(int channelSide) const;
	void setOutflowQs(int channelSide) const;
	void setNoSlipQs(int channelSide) const;
	void setGeoQs() const;
	void modifyQElement(int channelSide, unsigned int level) const;

	void initalQStruct(QforBoundaryConditions& Q,int channelSide, unsigned int level) const;
	void printQSize(std::string bc,int channelSide, unsigned int level) const;
	void setSizeNoSlip(int channelSide, unsigned int level) const;
	void setSizeGeoQs(unsigned int level) const;
	void setQ27Size(QforBoundaryConditions &Q, doubflo* QQ, unsigned int sizeQ) const;
	bool hasQs(int channelSide, unsigned int level) const;
};

#endif
