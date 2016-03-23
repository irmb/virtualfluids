/*
*  PathLinePostprocessor.h
*
*  Created on: 3.6.2012
*  Author: Fard
*/


#ifndef PATHLINEPOSTPROCESSORMCPART_H
#define PATHLINEPOSTPROCESSORMCPART_H

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
//#include <boost/tr1/memory.hpp>
#include "Particles.h"
#include "basics/writer/WbWriterVtkXmlASCII.h"

class D3Q27PathLinePostprocessorMcpart;
typedef boost::shared_ptr<D3Q27PathLinePostprocessorMcpart> D3Q27PathLinePostprocessorMcpartPtr;

class D3Q27PathLinePostprocessorMcpart : public Postprocessor  
{
public:
   //typedef std::map<int,ParticlesPtr> ParticlesMap;
   std::map<int,int> neighbors;
   //std::map<ofstream,std::vector<int> > pathOfstream;
   D3Q27PathLinePostprocessorMcpart(Grid3DPtr grid, const std::string& path, WbWriter* const writer,
      LBMUnitConverterPtr conv, UbSchedulerPtr s, CommunicatorPtr comm, 
      std::vector<UbTupleDouble3 > Positions, LBMReal nue, D3Q27InterpolationProcessorPtr iProcessor);
   ~D3Q27PathLinePostprocessorMcpart();
   void update(double step);
   void getNeighborsRank();

protected:
   void collectPostprocessData();
   //void getParticles(std::vector<Particles> *particles,std::vector<UbTupleDouble3> x);
   void initializeParticle(double x, double y,double z,int i,int level);
   void getParticlePerProccess(std::vector<UbTupleDouble3 >);
   void getParticle(Particles &particle, double dx,double x, double y,double z,int level,int i);
   void initialMovement();
   void updateDistributedValues(std::vector<Particles> *particles,std::vector<UbTupleDouble3> &x);
   void receiveStatusOfPoint(bool &status,int rankRoot,double x1,double x2,double x3,double level);
   void sendStatusOfPoint(bool status,double x1,double x2,double x3,double level,int i);
   void receiveBool(bool &variable,int rankBlock);
   void receiveInt(int &rvalue,int root);
   void sendBool(bool _bool,int distinationRank,int i );
   void sendInt(int svalue,int distinationRank ,int tag);
   void getVectorFromParticles(std::vector<double> &particlesInfomation,ParticlesPtr particles);
   void getParticlesFromVector(std::vector<double> particlesInfomation,int numberOFVariable);
   void updateinfo(std::vector<Particles> *particlesVec,std::vector<UbTupleDouble3> &x);
   void allGatherDoubles(std::vector<double>& svalues, std::vector<double>& rvalues);
   void fillOutPositions(std::vector<double> particlesInfomation,std::vector<UbTupleDouble3> &x,double numberOFVariable);
   void updateParticles();
   void checkParticles();
   void finalMovement(ParticlesPtr particle,double tau[]);
   void interpolMovement(Block3DPtr block,DistributionArray3DPtr& distributions,ParticlesPtr particle,double tau[],int ix1, int ix2, int ix3);
   void extrapolMovement(Block3DPtr block, ParticlesPtr particle,double tau[],int ix1, int ix2, int ix3);
   void findPlane(int ix1,int ix2,int ix3,Grid3DPtr grid,Block3DPtr block,double &A,double &B,double &C,double &D,double &ii);
   void CalcVelParticle2(double dx,double &vx1,double &vx2,double &vx3,double vx1old,double vx2old,double vx3old,double &vx1oldf,double &vx2oldf,double &vx3oldf);
   void CalcVelParticle(double dx,double &vx1,double &vx2,double &vx3,double vx1old,double vx2old,double vx3old,double &vx1oldf,double &vx2oldf,double &vx3oldf);
   void findInitialCell(double s,double di,double dj,double dk,double dx,double ix1ph,double ix2ph,double ix3ph,double &xp000,double &yp000,double &zp000);
   void finddisPointToCentercell(double x1,double x2,double x3,double xp000, double yp000, double zp000,double &dis,int level);
   void collWall(double A,double B,double C,double D,double &x1,double &x2,double &x3,double x1old,double x2old,double x3old,double dx,double &vx1,double &vx2,double &vx3,double ii);
   void addCorrection(double &vx1,double &vx2,double &vx3,double &tauxx,double &tauyy,double &tauzz,double &tauxy,double &tauxz,double &tauyz,
      double dx,D3Q27ICell iCell,double Xw, double Yw, double Zw,double  omega,double cofWeightx,double cofWeighty,double cofWeightz,double ii,double x1LB,double x2LB,double x3LB,double di,double dj,double dk);
   void getAllDis(double dis[],double x,double y,double z,double ix1ph,double ix2ph,double ix3ph,double xp000,double yp000,double zp000,int level);
   void rearangedDouble( double dis[7]);
   void getWallPosition(double dx,double xp000,double yp000,double zp000,double &Xw,double &Yw,double &Zw,double x1,double x2,double x3,double A
      ,double B,double C,double D,double normalDis, double di,double dj,double dk, double &cofWeightx,double &cofWeighty,double &cofWeightz,int dir);
   void getIcell(D3Q27ICell &iCell,double xp000, double yp000, double zp000,int level);
   void getAfterCompare(double dis[],double dx,double ix1ph,double ix2ph,double ix3ph,double &xp000,double &yp000,double &zp000,double &Xw,
      double &Yw,double &Zw,double x1,double x2,double x3,double A,double B,double C,double D,double normalDis, double di,double dj,double dk, double &cofWeightx
      ,double &cofWeighty,double &cofWeightz,D3Q27ICell &iCell,int level);
   void printParticle(ParticlesPtr particle);
   void printParticle(int index);
   void initializeForPrinting(int number);
   void gatherData(ParticlesPtr particle);
private:
   bool checkNodes(Block3DPtr block,double _x1,double _x2,double _x3);
   CommunicatorPtr comm;
   //std::map<int ,std::vector<UbTupleDouble3> > nodes;
   //std::map<int ,std::vector<std::string> >datanames;
   //std::map<int ,std::vector<std::vector<double> > >data;

   std::vector<UbTupleDouble3> nodes;
   std::vector<std::string> datanames;
   std::vector<std::vector<double> > data;

   std::vector<boost::shared_ptr<std::ofstream> > files;
   std::string path;
   WbWriter* writer;
   LBMUnitConverterPtr conv;

   D3Q27InterpolationProcessorPtr interpolation;
   //std::vector<UbTupleDouble3> x;
   //std::vector<UbTupleDouble3> xold;
   //std::vector<UbTupleDouble3 >vxold;
   //std::vector<UbTupleDouble3 > vxoldf;
   bool particleHasMass;
   int isExtrapolation;
   int istep;
   int stepcheck;
   LBMReal nue;
   D3Q27InterpolationProcessorPtr iProcessor;
   D3Q27InterpolationHelperPtr iHelper;
   //Particles particles;
   int rank;
   int root;
   int maxtag;
   std::list<ParticlesPtr> particles;
   //std::vector<std::vector<UbTupleDouble3> > Positions;
   //Particles particles;

   //ParticlesMap particles;

};
#endif

