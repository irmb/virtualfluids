#ifndef GridReaderFiles_H
#define GridReaderFiles_H

#include <VirtualFluidsDefinitions.h>

#include "../GridProvider.h"

#include <vector>
#include <string>
#include <memory>

#include "LBM/LB.h"

class Parameter;
class BoundaryValues;
class BoundaryQs;
class CoordNeighborGeoV;

class GridReader
	: public GridProvider
{
private:
	bool binaer;
	std::vector<std::string> channelDirections;
	std::vector<std::string> channelBoundaryConditions;
	std::shared_ptr<CoordNeighborGeoV> neighX, neighY, neighZ, neighWSB;
	std::vector<std::shared_ptr<BoundaryValues> > BC_Values;

public:
	VF_PUBLIC GridReader(bool binaer, std::shared_ptr<Parameter> para);
    VF_PUBLIC ~GridReader();
	void allocArrays_CoordNeighborGeo()override;
	void allocArrays_BoundaryValues()override;
    void allocArrays_OffsetScale() override;

	void initalValuesDomainDecompostion(int level);

	void setChannelBoundaryCondition();

	void allocArrays_BoundaryQs()override;
	bool getBinaer();
	void setDimensions();
	void setBoundingBox();
	void initPeriodicNeigh(std::vector<std::vector<std::vector<unsigned int> > > periodV, std::vector<std::vector<unsigned int> > periodIndex, std::string way);

private:
	void makeReader(std::shared_ptr<Parameter> para);
	void makeReader(std::vector<std::shared_ptr<BoundaryQs> > &BC_Qs, std::shared_ptr<Parameter> para);

	void setPressureValues(int channelSide) const;
	void setPressRhoBC(int sizePerLevel, int level, int channelSide) const;

	void setVelocityValues(int channelSide) const;
	void setVelocity(int level, int sizePerLevel, int channelSide) const;

	void setOutflowValues(int channelSide) const;
	void setOutflow(int level, int sizePerLevel, int channelSide) const;


	void setPressQs(std::shared_ptr<BoundaryQs> boundaryQ) const;
	void setVelocityQs(std::shared_ptr<BoundaryQs> boundaryQ) const;
	void setOutflowQs(std::shared_ptr<BoundaryQs> boundaryQ) const;
	void setNoSlipQs(std::shared_ptr<BoundaryQs> boundaryQ) const;
	void setGeoQs(std::shared_ptr<BoundaryQs> boundaryQ) const;
	void modifyQElement(std::shared_ptr<BoundaryQs> boundaryQ, unsigned int level) const;

	void initalQStruct(QforBoundaryConditions& Q, std::shared_ptr<BoundaryQs> boundaryQ, unsigned int level) const;
	void printQSize(std::string bc, std::shared_ptr<BoundaryQs> boundaryQ, unsigned int level) const;
	void setSizeNoSlip(std::shared_ptr<BoundaryQs> boundaryQ, unsigned int level) const;
	void setSizeGeoQs(std::shared_ptr<BoundaryQs> boundaryQ, unsigned int level) const;
	void setQ27Size(QforBoundaryConditions &Q, real* QQ, unsigned int sizeQ) const;
	bool hasQs(std::shared_ptr<BoundaryQs> boundaryQ, unsigned int level) const;
public:
    void initalGridInformations() override;
};

#endif
