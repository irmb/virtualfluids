#ifndef D3Q27INTEGRATEVALUESHELPER_H
#define D3Q27INTEGRATEVALUESHELPER_H

#include "Grid3D.h"
#include "D3Q27System.h"
#include "Communicator.h"
#include "GbCuboid3D.h"

struct CalcNodes 
{
	Block3DPtr block;
	std::vector<UbTupleInt3> nodes;
};

class D3Q27IntegrateValuesHelper;
typedef boost::shared_ptr<D3Q27IntegrateValuesHelper> D3Q27IntegrateValuesHelperPtr;

class D3Q27IntegrateValuesHelper
{
public:
	D3Q27IntegrateValuesHelper(Grid3DPtr grid, CommunicatorPtr comm, 
		double minX1, double minX2, double minX3, 
		double  maxX1, double maxX2, double maxX3);
	virtual ~D3Q27IntegrateValuesHelper();

	void calculateMQ();
	void calculateAV();
	void clearData();

	LBMReal getRho() {return sRho;}
	LBMReal getVx1() {return sVx1;} 
	LBMReal getVx2() {return sVx2;}
	LBMReal getVx3() {return sVx3;}
   LBMReal getCellsVolume() { return sCellVolume; }
 //  LBMReal getVm() { return sVm; }
	//LBMReal getPress() {return sPress;}
	LBMReal getAvVx1(){return sAvVx1;}
	LBMReal getAvVx2(){return sAvVx2;}
	LBMReal getAvVx3(){return sAvVx3;}
	LBMReal getTSx1(){return sTSx1;}
	LBMReal getTSx2(){return sTSx2;}
	LBMReal getTSx3(){return sTSx3;}
	LBMReal getTSx1x3(){return sTSx1x3;}
	LBMReal getNumberOfFluidsNodes();
	LBMReal getNumberOfSolidNodes();
	GbCuboid3DPtr getBoundingBox();
   std::vector<CalcNodes> getCNodes();

protected:
private:
	void init();
	Grid3DPtr grid;
   LBMReal sVx1, sVx2, sVx3, sRho, sCellVolume;// sPress, sVm;
	LBMReal numberOfFluidsNodes, numberOfSolidNodes;
	LBMReal sAvVx1, sAvVx2, sAvVx3, sTSx1, sTSx2, sTSx3, sTSx1x3;
	std::vector<CalcNodes> cnodes;
	GbCuboid3DPtr boundingBox;
	CommunicatorPtr comm;
	enum Values{AvVx = 0, AvVy = 1, AvVz = 2, AvVxx = 3, AvVyy = 4, AvVzz = 5, AvVxy = 6, AvVyz = 7, AvVxz = 8};

};

#endif
