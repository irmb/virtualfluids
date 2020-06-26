#ifndef INTEGRATEVALUESHELPER_H
#define INTEGRATEVALUESHELPER_H

#include <PointerDefinitions.h>

#include "Grid3D.h"
#include "D3Q27System.h"
#include "Communicator.h"
#include "GbCuboid3D.h"
#include "CbArray2D.h"
#include "Block3D.h"

//struct CalcNodes 
//{
//	SPtr<Block3D> block;
//	std::vector<UbTupleInt3> nodes;
//};
//
//struct Nodes
//{
//   SPtr<Block3D> block;
//   UbTupleInt3 nodes;
//};


class IntegrateValuesHelper
{
public:
   struct CalcNodes
   {
      SPtr<Block3D> block;
      std::vector<UbTupleInt3> nodes;
   };

   struct Node
   {
      SPtr<Block3D> block;
      UbTupleInt3 node;
   };
public:
	IntegrateValuesHelper(SPtr<Grid3D> grid, SPtr<Communicator> comm, 
		double minX1, double minX2, double minX3, 
		double  maxX1, double maxX2, double maxX3);
   IntegrateValuesHelper(SPtr<Grid3D> grid, SPtr<Communicator> comm,
      double minX1, double minX2, double minX3,
      double  maxX1, double maxX2, double maxX3, int level);
	virtual ~IntegrateValuesHelper();

	void calculateMQ();
	void calculateAV();
	void clearData();

	double getRho() {return sRho;}
	double getVx1() {return sVx1;} 
	double getVx2() {return sVx2;}
	double getVx3() {return sVx3;}
   double getCellsVolume() { return sCellVolume; }
 //  LBMReal getVm() { return sVm; }
	//LBMReal getPress() {return sPress;}
	double getAvVx1(){return sAvVx1;}
	double getAvVx2(){return sAvVx2;}
	double getAvVx3(){return sAvVx3;}
	double getTSx1(){return sTSx1;}
	double getTSx2(){return sTSx2;}
	double getTSx3(){return sTSx3;}
	double getTSx1x3(){return sTSx1x3;}

	LBMReal getNumberOfFluidsNodes();
	LBMReal getNumberOfSolidNodes();
	GbCuboid3DPtr getBoundingBox();
   std::vector<CalcNodes> getCNodes();

protected:
private:
	void init(int level);

   bool root;
	SPtr<Grid3D> grid;
   double sVx1, sVx2, sVx3, sRho, sCellVolume;// sPress, sVm;
   double numberOfFluidsNodes, numberOfSolidNodes;
   double sAvVx1, sAvVx2, sAvVx3, sTSx1, sTSx2, sTSx3, sTSx1x3;
	std::vector<CalcNodes> cnodes;
	GbCuboid3DPtr boundingBox;
	SPtr<Communicator> comm;
   CbArray2D<Node> cnodes2DMatrix;
	enum Values{AvVx = 0, AvVy = 1, AvVz = 2, AvVxx = 3, AvVyy = 4, AvVzz = 5, AvVxy = 6, AvVyz = 7, AvVxz = 8};
};

#endif