#pragma once
#include "Standard.h"
//#include "VTKwriter.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "Temperature/FindTemperature.h"
#include "OffsetScale.h"
#include "BoundaryValues.h"
#include <cstdlib>


using namespace std;
class Interface
{
private:
	bool binaer;
	string *system;
	string *way;
	CoordNeighborGeoV *neighX, *neighY, *neighZ, *neighWSB;
public:
	Interface(bool binaer);
	~Interface(void);
	void allocArrays_CoordNeighborGeo(Parameter* para);
	void allocArrays_OffsetScale(Parameter* para);
	void allocArrays_BoundaryValues(Parameter* para);
	void allocArrays_BoundaryQs(Parameter* para);
	bool getBinaer();
	void setDimensions(Parameter* para);
	void setBoundingBox(Parameter* para);
	void sortSystem(BoundaryValues **BC_Values, Parameter* para, int t);
	void sortSystem(BoundaryQs **BC_Qs, Parameter* para, int t);

	void initPeriodicNeigh(vector<vector<vector<unsigned int> > > periodV, vector<vector<unsigned int> > periodIndex, string way);

	void rearrangeGeometry(Parameter* para, int lev);
};

