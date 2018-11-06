#ifndef GridReaderGenerator_H
#define GridReaderGenerator_H

#include "../GridProvider.h"
#include <VirtualFluidsDefinitions.h>
#include <vector>
#include <string>
#include <memory>

#include "Core/DataTypes.h"

class Parameter;
class GridBuilder;
struct QforBC;
typedef QforBC QforBoundaryConditions;

class GridGenerator
	: public GridProvider
{
private:
	std::vector<std::string> channelDirections;
	std::vector<std::string> channelBoundaryConditions;

	std::shared_ptr<GridBuilder> builder;

public:
    VF_PUBLIC static std::shared_ptr<GridProvider> make(std::shared_ptr<GridBuilder> builder, std::shared_ptr<Parameter> para);


    VF_PUBLIC GridGenerator(std::shared_ptr<GridBuilder> builder, std::shared_ptr<Parameter> para);
	VF_PUBLIC virtual ~GridGenerator();

	void allocArrays_CoordNeighborGeo() override;
	void allocArrays_BoundaryValues() override;
	void allocArrays_BoundaryQs() override;
    void allocArrays_OffsetScale() override;

	virtual void setDimensions();
	virtual void setBoundingBox();

	virtual void initPeriodicNeigh(std::vector<std::vector<std::vector<unsigned int> > > periodV, std::vector<std::vector<unsigned int> > periodIndex, std::string way);
	
    void initalGridInformations() override;


private:
    std::string verifyNeighborIndices(int level) const;
    std::string verifyNeighborIndex(int level, int index, int &invalidNodes, int &stopperNodes, int &wrongNeighbors) const;
    std::string checkNeighbor(int level, real x, real y, real z, int index, int& numberOfWrongNeihgbors, int neighborIndex, real neighborX, real neighborY, real neighborZ, std::string direction) const;

};

#endif
