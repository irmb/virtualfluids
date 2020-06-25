#pragma once
#include "Standard.h"
#include "CoordNeighborGeoV.h"
#include "Parameter/Parameter.h"

using namespace std;

class BoundaryQs
	:public CoordNeighborGeoV
{
protected:
	vector<doubflo> vec1D_column;
	vector< vector<doubflo> >vec2D_lvl;
	vector< vector<vector<doubflo> > >vec3D_data;


	vector<unsigned int> vec1D_index;
	vector< vector<unsigned int> >vec2D_index;

	string boun;
	int column;

public:
	BoundaryQs();
	BoundaryQs(string q, bool binaer);
	BoundaryQs(string q, Parameter* para, string str, bool binaer);
	~BoundaryQs(void);
	void init();
	void init_Binary();

	void init_shortQs();

	void initArray(doubflo* ptr, unsigned int level, unsigned int column);
	void initProcNeighbor(int* ptr, unsigned int level, unsigned int column);
	void initIndex(/*unsigned*/ int *ptr, unsigned int level);

	vector<vector<vector<doubflo> > > setgetBoundarys(vector<vector<vector<doubflo> > > qs);
	vector<vector<vector<unsigned int> > > setgetBoundarys(vector<vector<vector<unsigned int> > > qs);
	vector<vector<unsigned int> > setgetIndex(vector<vector<unsigned int> > index);
	int getcolumn();

	void sout();
};

