/*
*  PathLinePostprocessor.h
*
*  Created on: 04.24.2012
*  Author: Fard
*/


#ifndef PATHLINEPOSTPROCESSOR_H
#define PATHLINEPOSTPROCESSOR_H

# define dirx 0
# define diry 1
# define dirz 2
# define dirxy 3 
# define dirxz 4 
# define diryz 5
# define dirxyz 6
#include <mpi.h>
#include "Postprocessor.h"
#include "Grid3D.h"
#include "Block3D.h"
#include "LBMUnitConverter.h"
#include "Communicator.h"
#include "D3Q27InterpolationProcessor.h"
#include "LBMKernelETD3Q27.h"
#include "D3Q27ETBCProcessor.h"
#include <boost/shared_ptr.hpp>
#include "WbWriter.h"

class D3Q27PathLinePostprocessor;
typedef boost::shared_ptr<D3Q27PathLinePostprocessor> D3Q27PathLinePostprocessorPtr;

class D3Q27PathLinePostprocessor : public Postprocessor
{
public:
   D3Q27PathLinePostprocessor(Grid3DPtr grid, const std::string& path, WbWriter* const writer,
      LBMUnitConverterPtr conv, UbSchedulerPtr s, CommunicatorPtr comm, 
      double x1, double x2, double x3,  LBMReal nue, D3Q27InterpolationProcessorPtr iProcessor);
   ~D3Q27PathLinePostprocessor();
   void update(double step);

protected:
   void collectPostprocessData();
   void initialMovement(Block3DPtr block, LBMKernelETD3Q27Ptr& kernel, DistributionArray3DPtr& distributions,  BCArray3D<D3Q27BoundaryCondition>& bcArray);
   void interpolMovement(Block3DPtr block, LBMKernelETD3Q27Ptr& kernel, DistributionArray3DPtr& distributions,  BCArray3D<D3Q27BoundaryCondition>& bcArray,int &isExtrapolation,double tau[]);
   void extrapolMovement(Block3DPtr block, LBMKernelETD3Q27Ptr& kernel, DistributionArray3DPtr& distributions,  BCArray3D<D3Q27BoundaryCondition>& bcArray,double tau[]);
   void clearData();
   void outICell(D3Q27ICell& iCell);
   void updateDistributedValues();
   void printPoint(double tau[]);
   void findPlane(int ix1,int ix2,int ix3,Grid3DPtr grid,Block3DPtr block,double &A,double &B,double &C,double &D,double &ii);
   void finddisPointToCentercell(double x1,double x2,double x3,double xp000, double yp000, double zp000,double &dis,double &dx,double &deltaT,Block3DPtr &block,int level);

   //void getWallPosition(Grid3DPtr grid,Block3DPtr &block,double xp000,double yp000,double zp000,double &Xw,double &Yw,double &Zw,double x1,double x2,double x3,
   //  double &xp000T,double &yp000T,double &zp000T,double A,double B,double C,double D,double normalDis,double di,double dj,double dk,double &cofWeightx
   //,double &cofWeighty,double &cofWeightz,double dismax,double dis,int level);
   void getWallPosition(double input[]);
   void findInitialCell(double s,double di,double dj,double dk,double dx,double ix1ph,double ix2ph,double ix3ph,double &xp000,double &yp000,double &zp000);
   void rearangeInt( int dis[7]);
   void rearangedDouble( double dis[7]);
   void addCorrection(double &vx1,double &vx2,double &vx3,double &tauxx,double &tauyy,double &tauzz,double &tauxy,double &tauxz,double &tauyz,
      double dx,D3Q27ICell iCell,double Xw, double Yw, double Zw,double  omega,double cofWeightx,double cofWeighty,double cofWeightz,double ii,double x1LB,double x2LB,double x3LB,double di,double dj,double dk);
   void collWall(double A,double B,double C,double D,double &x1,double &x2,double &x3,double x1old,double x2old,double x3old,double dx,double &vx1,double &vx2,double &vx3,double ii);
   void CalcVelParticle(double dx,double &vx1,double &vx2,double &vx3,double vx1old,double vx2old,double vx3old,double &vx1oldf,double &vx2oldf,double &vx3oldf);
   void CalcVelParticle2(double dx,double &vx1,double &vx2,double &vx3,double vx1old,double vx2old,double vx3old,double &vx1oldf,double &vx2oldf,double &vx3oldf);

   // void getWallPosition(double dixmax,double disx,double disy,double disz,double disxy,double disxz,double disyz,double disxyz,double &Xw,double &Yw,double &Zw,double x1,double x2,double x3,double A,double B,double C,double D,double normalDis,double di,double dj,double dk,
   //   double &cofWeightx,double &cofWeighty,double &cofWeightz);
   void getRankWithPosition(Grid3DPtr grid,double xp000, double yp000, double zp000,int &RankNeigber,int level);
   void getRankWithPosition(Grid3DPtr grid,double xp000, double yp000, double zp000,int &RankNeigber);
   void getRankWithPosition2(Grid3DPtr grid,double xp000, double yp000, double zp000,int &RankNeigber);
   void getAllTag(int allRank[7],int allTag[7]);
   void getAllInfoForCompare(double input[],double *Vect);
   void sendMaxtagToAllProccess(int allRank[],int allTag[],int numberpro,int &numberproccessused,int *mainproccess);
   void getAllRankWithAllPositions(double ix1ph,double ix2ph,double ix3ph,double xp000,double yp000,double zp000,int &Rankx,int &Ranky,int &Rankz,int &Rankxy,int &Rankxz,int &Rankyz,int &Rankxyz,int level);
   void getAllRankWithAllPositions2(double ix1ph,double ix2ph,double ix3ph,double xp000,double yp000,double zp000,int &Rankx,int &Ranky,int &Rankz,int &Rankxy,int &Rankxz,int &Rankyz,int &Rankxyz);
   void getVectfromiCell(D3Q27ICell iCell,double *Vect);
   void getiCellfromVect(D3Q27ICell &iCell,double *Vect);

   void getAfterCompare(double dis[],double input [7][25],double AllCell[7][8*27],double &xp000,double &yp000,double &zp000,
      double &Xw,double &Yw,double &Zw,double &cofWeightx,double &cofWeighty,double &cofWeightz,double &dx,double &omega,D3Q27ICell &iCell);
   void SendAndReceiveData(double input [7][25],double AllCell[7][8*27],int allRank[7],int allTag[7]);	
   void sendtoOtherIsOk(int isPointSuitable);
private:
   bool checkNodes(Block3DPtr block);
   bool checkNodes2( Block3DPtr block,double x11,double x22,double x33);
   void checkLevel(Block3DPtr& block, LBMKernelETD3Q27Ptr& kernel, DistributionArray3DPtr& distributions,  BCArray3D<D3Q27BoundaryCondition>& bcArray,bool &isPointSuitable2);

   std::vector<UbTupleDouble3> nodes;
   std::vector<std::string> datanames;
   std::vector<std::vector<double> > data;
   std::string path;
   WbWriter* writer;
   LBMUnitConverterPtr conv;
   CommunicatorPtr comm;
   D3Q27InterpolationProcessorPtr interpolation;
   double x1, x2, x3;
   double x1old, x2old, x3old;
   LBMReal vx1old,vx2old,vx3old;
   LBMReal vx1oldf,vx2oldf,vx3oldf;
   bool particleHasMass;
   int isExtrapolation;
   int istep;
   int stepcheck;
   LBMReal nue;
   D3Q27InterpolationProcessorPtr iProcessor;
   D3Q27InterpolationHelperPtr iHelper;
   int rank;
   int root;
   int maxtag;
   std::vector<std::vector<double> > Positions;
};
#endif

