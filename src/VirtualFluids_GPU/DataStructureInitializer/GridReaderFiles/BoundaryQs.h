#ifndef BoundaryQs_H
#define BoundaryQs_H

#include <vector>
#include <string>
#include <fstream>
#include <memory>

#include "LBM/LB.h"

class Parameter;

class BoundaryQs
{
public:
	BoundaryQs();
	BoundaryQs(std::string q, bool binaer);
	BoundaryQs(std::string q, std::shared_ptr<Parameter> para, std::string str, bool binaer);
	~BoundaryQs(void);

public:
	unsigned int getSize(unsigned int level);
	unsigned int getLevel();


private:	
	void checkFileStatus(std::string path);
	void init();
	void resizeVectors();
	void resizeVectorsPerLevel(unsigned int level, std::vector<uint32_t> &vec1D_code);

	void init_Binary();

public:
	void setIndex(int *indices, unsigned int level) const;
	void setValues(doubflo** q27, unsigned int level) const;

private:
	std::vector< std::vector<std::vector<doubflo> > >values;
	std::vector< std::vector<unsigned int> >indices;

	std::ifstream file;
	unsigned int maxLevel;
	std::vector<unsigned int> levelSizes;

};

#endif