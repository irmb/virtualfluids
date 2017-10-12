#ifndef BoundaryValues_H
#define BoundaryValues_H

#include <string>
#include <vector>
#include <fstream>
#include <memory>

#include "LBM/LB.h"

class Parameter;

class BoundaryValues
{
private:
	std::string boundaryCondition;
	bool procNeighbor;

	std::vector< std::vector<std::vector<doubflo> > >values;
	std::vector< std::vector<unsigned int> >indices;
	std::vector<unsigned int> levelSizes;

	std::ifstream file;
	unsigned int maxLevel;
	
public:
	BoundaryValues(std::string path);
	BoundaryValues(std::string path, std::shared_ptr<Parameter> para, std::string str);
	BoundaryValues(int neighbor, std::shared_ptr<Parameter> para, std::string sor, std::string dir);
	~BoundaryValues();

	unsigned int getLevel();
	unsigned int getSize(unsigned int level);
	std::string getBoundaryCondition();
	void setBoundarys(std::vector<std::vector<std::vector<doubflo> > > &qs) const;
	void initIndex(/*unsigned*/ int *ptr, unsigned int level);

	void setProcNeighbor(bool pN);
	bool getProcNeighbor();

	void setPressValues(doubflo *RhoBC, int* kN, int level) const;
	void setVelocityValues(doubflo *vx, doubflo *vy, doubflo *vz, int level) const;
	void setOutflowValues(doubflo *RhoBC, int* kN, int level) const;

private:
	void init();
	int getNumberOfColumns();
	void initalVectors(unsigned int maxColumn);
	void readData(unsigned int level, int index, unsigned int maxColumn);
	void resizeVectorsPerLevel(unsigned int level, unsigned int maxColumn);
	void skipLine();
	void readLevelSize(unsigned int level);
	void initalVectorsWithSingleZero();
	void readNumberOfLevels();
	void readBC();
	void resizeVectors();


};

#endif
