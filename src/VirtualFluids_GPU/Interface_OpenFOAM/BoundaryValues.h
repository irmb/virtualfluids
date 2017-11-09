#pragma once
#include "Standard.h"
#include "BoundaryQs.h"

class BoundaryValues
	: public BoundaryQs
{
private:
	string bo;
	bool procNeighbor;

public:
	BoundaryValues(string ad);
	BoundaryValues(string ad, Parameter* para, string str);
	BoundaryValues(int neighbor, Parameter* para, string sor);
	BoundaryValues(int neighbor, Parameter* para, string sor, string dir);
	~BoundaryValues(void);
	void init();
	vector<unsigned int> getIndex(unsigned int level);
	vector<doubflo> getVec(unsigned int level, unsigned int column);
	string getWay();

	void setProcNeighbor(bool pN);
	bool getProcNeighbor();
};

