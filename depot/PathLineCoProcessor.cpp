#include "PathLineCoProcessor.h"
#include "LBMKernelETD3Q27.h"
#include "SimulationParameters.h"
#include "D3Q27ETBCProcessor.h"
#include <vector>
#include <string>
#include <boost/foreach.hpp>
#include "WbWriter.h"

using namespace std;


PathLineCoProcessor::PathLineCoProcessor(Grid3DPtr grid, const std::string& path,
                                                       WbWriter* const writer, LBMUnitConverterPtr conv,
                                                       UbSchedulerPtr s, CommunicatorPtr comm, 
                                                       double x1, double x2, double x3, LBMReal nue, D3Q27InterpolationProcessorPtr iProcessor)
                                                       : CoProcessor(grid, s),
                                                       path(path),
                                                       comm(comm),
                                                       writer(writer),
                                                       conv(conv),
                                                       x1(x1),
                                                       x2(x2),
                                                       x3(x3),
                                                       nue(nue),
                                                       iProcessor(iProcessor),
                                                       vx1old(0.0),
                                                       vx2old(0.0),
                                                       vx3old(0.0),
                                                       vx1oldf(0.1),
                                                       vx2oldf(0.1),
                                                       vx3oldf(0.1),
                                                       isExtrapolation(0),
                                                       istep(0),
                                                       maxtag(-1),
                                                       particleHasMass(true),
                                                       rank(comm->getProcessID())
{
   iHelper = D3Q27InterpolationHelperPtr(new D3Q27InterpolationHelper(iProcessor));
}
//////////////////////////////////////////////////////////////////////////
PathLineCoProcessor::~PathLineCoProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::process(double step)
{
   this->stepcheck=(int)step;
   double aa=scheduler->getMinEnd();
   if (step >= scheduler->getMaxBegin() && step <= scheduler->getMinEnd())
   {
      collectPostprocessData();
   }

   UBLOG(logDEBUG3, "PathLinePostprocessor::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::collectPostprocessData()
{
   //for (int n=0;n<Positions[][0].size();n++)
   //{
   //}

   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   int maxInitLevel = this->grid->getFinestInitializedLevel();
   bool goinside=true;

   LBMKernelETD3Q27Ptr kernel;
   DistributionArray3DPtr distributions;
   BCArray3D<D3Q27BoundaryCondition> bcArray;
   Block3DPtr block;
   Block3DPtr block2;
   int isPointSuitable=0;//0 means false
   ///////send tau for printing to proccess 0
   double tau[6];//={tauxx,tauyy,tauzz,tauxy,tauxz,tauyz};
   for(int level = minInitLevel; level<=maxInitLevel; level++)
   {	      
      UbTupleInt3 blockIndexes = grid->getBlockIndexes(x1, x2, x3,level);
      block = grid->getBlock(val<1>(blockIndexes), val<2>(blockIndexes), val<3>(blockIndexes), level);
      if(!block) continue; 
      if(block->isActive()) 
      {
         this->root=block->getRank();
         if(rank == block->getRank())
         { 
            if(!checkNodes(block))
            {		   
               isPointSuitable=0;//0 meas false
               if(comm->getNumberOfProcesses() > 1)   {  sendtoOtherIsOk(isPointSuitable);}
               continue;
            }  
            initialMovement(block, kernel, distributions,  bcArray);
            isPointSuitable=1;//0 meas false
            if(comm->getNumberOfProcesses() > 1)   {  sendtoOtherIsOk(isPointSuitable);}
            break;
         }	
         else
         {
            if(comm->getNumberOfProcesses() > 1)  
            { 
               MPI_Status status; 
               MPI_Recv(&isPointSuitable,1, MPI_INT,this->root,this->rank,MPI_COMM_WORLD,&status);
               if (isPointSuitable)
               {
                  break;
               }
            }
         } 
      }
   } 
   if (isPointSuitable)
   {
      if(comm->getNumberOfProcesses() > 1) updateDistributedValues();
      bool isPointSuitable2;
      checkLevel(block, kernel, distributions,  bcArray,isPointSuitable2);
      if (isPointSuitable2)
      {  
         if(block->isActive())
         { 
            if (rank == block->getRank())
            {
               interpolMovement(block, kernel, distributions,  bcArray,isExtrapolation,tau);
               if (isExtrapolation)
               {	
                  extrapolMovement(block, kernel, distributions,  bcArray,tau);
               }
            }
            else
            {
               if(comm->getNumberOfProcesses() > 1)  
               {  
                  //recive isExtrapolation
                  if (this->rank!=this->root)
                  {
                     MPI_Status status; 
                     MPI_Recv(&isExtrapolation,1, MPI_INT,this->root,this->rank,MPI_COMM_WORLD,&status);
                  }
               }
            }
         }   
         if(comm->getNumberOfProcesses() > 1) 
         {
            if (isExtrapolation)
            {	
               if (this->rank!=this->root)
               {
                  MPI_Status status;   
                  MPI_Recv(&maxtag,1, MPI_INT,this->root,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
               }	   
            }
            if (isExtrapolation)
            {	
               if (this->rank!=this->root)
               {
                  MPI_Status status; 
                  double input [7][25];
                  double AllCell[7][8*27];
                  int sizeinput=25;
                  for (int i=0;i<maxtag;i++)
                  {  
                     double vectforTrans[25]; 
                     double *Cell=(double*)malloc(8*27*sizeof(double));
                     MPI_Recv(vectforTrans,25, MPI_DOUBLE_PRECISION,this->root,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
                     getAllInfoForCompare(vectforTrans,Cell);

                     for (int j=0;j<sizeinput;j++)  {  input[i][j]=vectforTrans[j]; }
                     for (int k=0;k<8*27;k++)       {  AllCell[i][k]=Cell[k];   	  }
                     free (Cell);
                  }
                  for (int i=0;i<maxtag;i++)
                  { 
                     double vectforTrans[25];
                     double *Cell=(double*)malloc(8*27*sizeof(double));
                     for (int j=0;j<sizeinput;j++)   {	  vectforTrans[j]=input[i][j]; 	  }
                     for (int k=0;k<8*27;k++)       {    Cell[k]=AllCell[i][k];   	  }

                     MPI_Send(vectforTrans,25, MPI_DOUBLE_PRECISION,this->root,i,MPI_COMM_WORLD);
                     MPI_Send(Cell,8*27, MPI_DOUBLE_PRECISION,this->root,i,MPI_COMM_WORLD);
                     free (Cell);
                  }
               }
            }
            updateDistributedValues();
         }

         if(comm->getNumberOfProcesses() > 1)  MPI_Bcast(tau,6, MPI_DOUBLE, this->root, MPI_COMM_WORLD);	   
         if (rank==0)
         {
            MPI_Status status;   
            printPoint(tau);
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::printPoint(double tau[])
{
   double tauxx=tau[0];
   double tauyy=tau[1];
   double tauzz=tau[2];
   double tauxy=tau[3];
   double tauxz=tau[4];
   double tauyz=tau[5];
   datanames.resize(0);
   datanames.push_back("TimeStep");
   datanames.push_back("Vx");
   datanames.push_back("Vy");
   datanames.push_back("Vz");
   datanames.push_back("tauxx");
   datanames.push_back("tauyy");
   datanames.push_back("tauzz");
   datanames.push_back("tauxy");
   datanames.push_back("tauxz");
   datanames.push_back("tauyz");
   /*datanames.push_back("xoff");
   datanames.push_back("yoff");
   datanames.push_back("zoff");
   */
   data.resize(datanames.size());

   int index = 0;
   data[index++].push_back(istep++);
   data[index++].push_back(vx1old);
   data[index++].push_back(vx2old);
   data[index++].push_back(vx3old);
   data[index++].push_back(tauxx);
   data[index++].push_back(tauyy);
   data[index++].push_back(tauzz);
   data[index++].push_back(tauxy);
   data[index++].push_back(tauxz);
   data[index++].push_back(tauyz);
   //data[index++].push_back(xoff);
   //data[index++].push_back(yoff);
   //data[index++].push_back(zoff);

   nodes.push_back( makeUbTuple(  double(x1), double(x2), double(x3)) );

   //if(scheduler->isDue((double)istep) )
   if (istep%1 ==0)
   {
      string test = writer->writeNodesWithNodeDataDouble(path,nodes,datanames,data);
   }
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::clearData()
{
   nodes.clear();
   datanames.clear();
   data.clear();
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::initialMovement(Block3DPtr block, LBMKernelETD3Q27Ptr& kernel, DistributionArray3DPtr& distributions,  BCArray3D<D3Q27BoundaryCondition>& bcArray)
{
   //UBLOG(logINFO, "addPostprocessData1");
   //root = rank;

   double dx = grid->getDeltaX(block);
   kernel = boost::dynamic_pointer_cast<LBMKernelETD3Q27>(block->getKernel());
   distributions = kernel->getDataSet()->getFdistributions();
   bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();

   x1old=x1;
   x2old=x2;
   x3old=x3;

   x1=x1old+vx1old*conv->getFactorTimeLbToW(grid->getDeltaX(block));
   x2=x2old+vx2old*conv->getFactorTimeLbToW(grid->getDeltaX(block));
   x3=x3old+vx3old*conv->getFactorTimeLbToW(grid->getDeltaX(block));

}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::checkLevel(Block3DPtr& block, LBMKernelETD3Q27Ptr& kernel, DistributionArray3DPtr& distributions,  BCArray3D<D3Q27BoundaryCondition>& bcArray,bool &isPointSuitable2)
{
   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   int maxInitLevel = this->grid->getFinestInitializedLevel(); 
   bool goinside=true;
   int isPointSuitable=0;//0 means false
   for(int level = minInitLevel; level<=maxInitLevel; level++)
   {	      
      UbTupleInt3 blockIndexes = grid->getBlockIndexes(x1, x2, x3,level);
      block = grid->getBlock(val<1>(blockIndexes), val<2>(blockIndexes), val<3>(blockIndexes), level);
      if(!block) continue; 
      if(block->isActive()) 
      {
         this->root=block->getRank();
         if(rank == block->getRank())
         { 
            if(!checkNodes(block))
            {		   
               isPointSuitable=0;//0 meas false
               if(comm->getNumberOfProcesses() > 1)   {  sendtoOtherIsOk(isPointSuitable);}
               continue;
            }  
            kernel = boost::dynamic_pointer_cast<LBMKernelETD3Q27>(block->getKernel());
            distributions = kernel->getDataSet()->getFdistributions();
            bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();
            isPointSuitable=1;//0 meas false
            if(comm->getNumberOfProcesses() > 1)   {  sendtoOtherIsOk(isPointSuitable);}
            break;
         }	
         else
         {
            if(comm->getNumberOfProcesses() > 1)  
            { 
               MPI_Status status; 
               MPI_Recv(&isPointSuitable,1, MPI_INT,this->root,this->rank,MPI_COMM_WORLD,&status);
               if (isPointSuitable)
               {
                  break;
               }
            }
         } 
      }
   }
   if (isPointSuitable){isPointSuitable2=true;	} 
   else isPointSuitable2=false;
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::interpolMovement(Block3DPtr block, LBMKernelETD3Q27Ptr& kernel, DistributionArray3DPtr& distributions,  BCArray3D<D3Q27BoundaryCondition>& bcArray, int &isExtrapolation,double tau[])
{
   LBMReal tauxx,	   tauyy,	   tauzz,	   tauxy,	   tauxz,	   tauyz;
   LBMReal vx1,vx2,vx3,rho;
   UbTupleInt3 nodeIndexes = grid->getNodeIndexes(block, x1, x2, x3);
   LBMReal dx = grid->getDeltaX(block); 
   D3Q27ICell iCell;
   LBMReal f[D3Q27System::ENDF+1];
   int ix1 = val<1>(nodeIndexes);
   int ix2 = val<2>(nodeIndexes);
   int ix3 = val<3>(nodeIndexes);
   double xoff=0.0; 
   double yoff=0.0;
   double zoff=0.0;
   //typedef void (*CalcMacrosFct)(const LBMReal* const& /*feq[27]*/,LBMReal& /*(d)rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);
   /*
   CalcMacrosFct calcMacros = NULL;

   if(block->getKernel()->getCompressible())
   {
   calcMacros = &D3Q27System::calcCompMacroscopicValues;
   }
   else
   {
   calcMacros = &D3Q27System::calcIncompMacroscopicValues;
   }*/

   double Ucor=0,Vcor=0,Wcor=0;
   double tauxxC=00,tauyyC=0,tauzzC=0,tauxyC=0,tauxzC=0,tauyzC=0;
   double X_Wcor,Y_Wcor,Z_Wcor;
   double X_Pcor,Y_Pcor,Z_Pcor;
   double xp000,yp000,zp000;
   double Xw;
   double Yw;
   double Zw;

   UbTupleDouble3 orgNodeRW =  grid->getNodeCoordinates(block,  ix1, ix2, ix3);
   double ix1p= ( val<1>(orgNodeRW));
   double ix2p= ( val<2>(orgNodeRW));
   double ix3p= ( val<3>(orgNodeRW));
   if(!iProcessor->iCellHasSolid(bcArray, ix1, ix2, ix3))
   {
      iProcessor->readICell(distributions, iCell, ix1, ix2, ix3);
      double x1LB, x2LB, x3LB;
      UbTupleDouble3 orgNodeRW =  grid->getNodeCoordinates(block,  ix1, ix2, ix3);
      double offsetX1RW = abs(x1 - val<1>(orgNodeRW));
      double offsetX2RW = abs(x2 - val<2>(orgNodeRW));
      double offsetX3RW = abs(x3 - val<3>(orgNodeRW));
      x1LB = offsetX1RW / dx;
      x2LB = offsetX2RW / dx;
      x3LB = offsetX3RW / dx;
      x1LB -= 0.5-xoff;
      x2LB -= 0.5-yoff;
      x3LB -= 0.5-zoff;
      LBMReal omega = LBMSystem::calcCollisionFactor(nue, block->getLevel());
      /*iProcessor->interpolate8to1(iCell, f, x1LB, x2LB, x3LB, omega);
      calcMacros(f,rho,vx1,vx2,vx3);*/
      //iProcessor->interpolate8to1WithVelocity(iCell, f, x1LB, x2LB, x3LB, omega,vx1,vx2,vx3);
      iHelper->interpolate8to1WithVelocityWithShearStress(iCell, x1LB, x2LB, x3LB, omega,vx1,vx2,vx3,tauxx,tauyy,tauzz,tauxy,tauxz,tauyz);
      vx1 = vx1 *  conv->getFactorVelocityLbToW();
      vx2 = vx2 *  conv->getFactorVelocityLbToW();
      vx3 = vx3 *  conv->getFactorVelocityLbToW();
      tau[0]=tauxx/grid->getDeltaX(block);
      tau[1]=tauyy/grid->getDeltaX(block);
      tau[2]=tauzz/grid->getDeltaX(block);
      tau[3]=tauxy/grid->getDeltaX(block);
      tau[4]=tauxz/grid->getDeltaX(block);
      tau[5]=tauyz/grid->getDeltaX(block);  

      //if (particleHasMass) {   CalcVelParticle(dx,vx1,vx2,vx3,vx1old,vx2old,vx3old,vx1oldf,vx2oldf,vx3oldf);}
      if (particleHasMass) {   CalcVelParticle2(dx,vx1,vx2,vx3,vx1old,vx2old,vx3old,vx1oldf,vx2oldf,vx3oldf);}
      //heuns method
      x1 = (x1old) + 1.0/2.0*(vx1 + vx1old)*conv->getFactorTimeLbToW(dx);
      x2 = (x2old) + 1.0/2.0*(vx2 + vx2old)*conv->getFactorTimeLbToW(dx);
      x3 = (x3old) + 1.0/2.0*(vx3 + vx3old)*conv->getFactorTimeLbToW(dx);

      vx1old = vx1;
      vx2old = vx2;
      vx3old = vx3;

      //send isExtrapolation false
      isExtrapolation=0;//1 means true
      if(comm->getNumberOfProcesses() > 1)  
      { 
         int size=comm->getNumberOfProcesses();
         for (int i=0;i<size;i++)
         {
            if (i!=rank)
            {
               MPI_Send(&isExtrapolation,1, MPI_INT,i,i,MPI_COMM_WORLD); 
            }	
         }
      }
   }
   else
   {
      //send isExtrapolation true
      isExtrapolation=1;//1 means true
      if(comm->getNumberOfProcesses() > 1)  
      { 
         int size=comm->getNumberOfProcesses();
         for (int i=0;i<size;i++)
         {
            if (i!=rank)
            {
               MPI_Send(&isExtrapolation,1, MPI_INT,i,i,MPI_COMM_WORLD); 
            }	
         }
      }

   }
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::extrapolMovement(Block3DPtr block, LBMKernelETD3Q27Ptr& kernel, DistributionArray3DPtr& distributions,  BCArray3D<D3Q27BoundaryCondition>& bcArray,double tau[])
{
   LBMReal tauxx,	   tauyy,	   tauzz,	   tauxy,	   tauxz,	   tauyz;
   LBMReal vx1,vx2,vx3,rho;
   UbTupleInt3 nodeIndexes = grid->getNodeIndexes(block, x1, x2, x3);
   LBMReal dx = grid->getDeltaX(block);
   int    level = block->getLevel();
   double deltaT= 1.0/(double)(1<<level);	 
   LBMReal omega = LBMSystem::calcCollisionFactor(nue,level);
   D3Q27ICell iCell;
   LBMReal f[D3Q27System::ENDF+1];
   int ix1 = val<1>(nodeIndexes);
   int ix2 = val<2>(nodeIndexes);
   int ix3 = val<3>(nodeIndexes);
   double xoff=0.0; 
   double yoff=0.0;
   double zoff=0.0;
   double Ucor=0,Vcor=0,Wcor=0;
   double tauxxC=00,tauyyC=0,tauzzC=0,tauxyC=0,tauxzC=0,tauyzC=0;
   double X_Wcor,Y_Wcor,Z_Wcor;
   double X_Pcor,Y_Pcor,Z_Pcor;
   double xp000,yp000,zp000;
   double Xw=999.0;
   double Yw=999.0;
   double Zw=999.0;

   UbTupleDouble3 orgNodeRW =  grid->getNodeCoordinates(block,  ix1, ix2, ix3);
   double ix1p= ( val<1>(orgNodeRW));
   double ix2p= ( val<2>(orgNodeRW));
   double ix3p= ( val<3>(orgNodeRW));
   {
      UbTupleDouble3 orgNodeRW =  grid->getNodeCoordinates(block,  ix1, ix2, ix3);
      double ix1ph = ( val<1>(orgNodeRW));
      double ix2ph = ( val<2>(orgNodeRW));
      double ix3ph = ( val<3>(orgNodeRW)); 
      double cofWeightx=999.0,cofWeighty=999.0,cofWeightz=999.0;
      double A,B,C,D,ii=0.0;
      findPlane(ix1,ix2,ix3,grid,block,A,B,C,D,ii);
      double s = A*x1 + B*x2 + C*x3 + D;//The sign of s = Ax + By + Cz + D determines which side the point (x,y,z) lies with respect to the plane. If s > 0 then the point lies on the same side as the normal (A,B,C). If s < 0 then it lies on the opposite side, if s = 0 then the point (x,y,z) lies on the plane.
      if (s>0){s=1;} else if (s<0){s=-1;}else {s=0;}
      double normalDis=((A*x1 + B*x2 + C*x3 + D)/sqrt(A*A+B*B+C*C));///distance point to plane xp-Xw=distance
      double di=A/sqrt(A*A+B*B+C*C);      double dj=B/sqrt(A*A+B*B+C*C);       double dk=C/sqrt(A*A+B*B+C*C);
      double XXw=x1-di*normalDis;    double YYw=x2-dj*normalDis;    double ZZw=x3-dk*normalDis;
      findInitialCell(s,di,dj,dk,dx,ix1ph,ix2ph,ix3ph,xp000,yp000,zp000);//find initial cell with xp000,yp000,zp000
      /////////////////////////////////////////////////////////////////////////
      int level=block->getLevel();
      int Rankx,Ranky,Rankz,Rankxy,Rankxz,Rankyz,Rankxyz;
      double disx=9990,disy=9990,disz=9990,disxy=9990,disxz=9990,disyz=9990,disxyz=9990; 
      getAllRankWithAllPositions(ix1ph,ix2ph,ix3ph,xp000,yp000,zp000,Rankx,Ranky,Rankz,Rankxy,Rankxz,Rankyz,Rankxyz,level);
      //getAllRankWithAllPositions2(ix1ph,ix2ph,ix3ph,xp000,yp000,zp000,Rankx,Ranky,Rankz,Rankxy,Rankxz,Rankyz,Rankxyz);

      int allRank[7]={ Rankx,Ranky,Rankz,Rankxy,Rankxz,Rankyz,Rankxyz};
      int allTag[7];
      getAllTag(allRank,allTag);
      int numberproccessused=0;
      int numberpro=comm->getNumberOfProcesses();
      int *mainproccess=new int[numberpro] ; 
      sendMaxtagToAllProccess(allRank,allTag,numberpro,numberproccessused,mainproccess);

      //{0=xp000, 1=ix2ph,  2=ix3ph,3=x1,4=x2,5=x3,6=Xw,7=Yw,8=Zw,9=A,10=B,11=C,12=D,
      //13=normalDis,14=di,15=dj,16=dk,17=cofWeightx,18=cofWeighty,19=cofWeightz,,20=dx,21=omega,22=level,23=disx,24=dir  };
      const int sizeinput=25;
      double input [7][sizeinput]=
      {{xp000, ix2ph,  ix3ph,x1,x2,x3,Xw,Yw,Zw,A,B,C,D,normalDis,di,dj,dk,cofWeightx,cofWeighty,cofWeightz,dx,omega,level,disx  ,dirx  }
      ,{ix1ph, yp000,  ix3ph,x1,x2,x3,Xw,Yw,Zw,A,B,C,D,normalDis,di,dj,dk,cofWeightx,cofWeighty,cofWeightz,dx,omega,level,disy  ,diry  }
      ,{ix1ph, ix2ph,  zp000,x1,x2,x3,Xw,Yw,Zw,A,B,C,D,normalDis,di,dj,dk,cofWeightx,cofWeighty,cofWeightz,dx,omega,level,disz  ,dirz  }
      ,{xp000, yp000,  ix3ph,x1,x2,x3,Xw,Yw,Zw,A,B,C,D,normalDis,di,dj,dk,cofWeightx,cofWeighty,cofWeightz,dx,omega,level,disxy ,dirxy }
      ,{xp000, ix2ph,  zp000,x1,x2,x3,Xw,Yw,Zw,A,B,C,D,normalDis,di,dj,dk,cofWeightx,cofWeighty,cofWeightz,dx,omega,level,disxz ,dirxz }
      ,{ix1ph, yp000,  zp000,x1,x2,x3,Xw,Yw,Zw,A,B,C,D,normalDis,di,dj,dk,cofWeightx,cofWeighty,cofWeightz,dx,omega,level,disyz ,diryz }
      ,{xp000, yp000,  zp000,x1,x2,x3,Xw,Yw,Zw,A,B,C,D,normalDis,di,dj,dk,cofWeightx,cofWeighty,cofWeightz,dx,omega,level,disxyz,dirxyz}};

      double AllCell[7][8*27];
      int index;
      SendAndReceiveData(input,AllCell,allRank,allTag);

      double dis[7]={input[0][23],input[1][23],input[2][23],input[3][23],input[4][23],input[5][23],input[6][23]};                                                                  
      rearangedDouble(dis);
      UbTupleInt3 blockIndexes;

      getAfterCompare(dis,input,AllCell,xp000, yp000,zp000,Xw,Yw,Zw,cofWeightx,cofWeighty,cofWeightz,dx,omega,iCell);

      double x1LB, x2LB, x3LB;
      double offsetX1RW000 = (x1 - xp000);
      double offsetX2RW000 = (x2 - yp000);
      double offsetX3RW000 = (x3 - zp000);

      x1LB = offsetX1RW000 / dx;
      x2LB = offsetX2RW000 / dx;
      x3LB = offsetX3RW000 / dx;
      x1LB -= 0.5;
      x2LB -= 0.5;
      x3LB -= 0.5;
      // outICell(iCell);
      //iProcessor->interpolate8to1WithVelocity(iCell, f, x1LB, x2LB, x3LB, omega,vx1,vx2,vx3);
      iHelper->interpolate8to1WithVelocityWithShearStress(iCell, x1LB, x2LB, x3LB, omega,vx1,vx2,vx3,tauxx,tauyy,tauzz,tauxy,tauxz,tauyz);	    
      addCorrection(vx1,vx2,vx3,tauxx,tauyy,tauzz,tauxy,tauxz,tauyz,dx,iCell,Xw, Yw, Zw, omega,cofWeightx,cofWeighty,cofWeightz,ii, x1LB, x2LB, x3LB,di,dj,dk);
      tau[0]=tauxx;
      tau[1]=tauyy;
      tau[2]=tauzz;
      tau[3]=tauxy;
      tau[4]=tauxz;
      tau[5]=tauyz;
      // if (particleHasMass) {   CalcVelParticle(dx,vx1,vx2,vx3,vx1old,vx2old,vx3old,vx1oldf,vx2oldf,vx3oldf);}
      if (particleHasMass) {   CalcVelParticle2(dx,vx1,vx2,vx3,vx1old,vx2old,vx3old,vx1oldf,vx2oldf,vx3oldf);}
      //heuns method
      x1 = (x1old) + 1.0/2.0*(vx1 + vx1old)*conv->getFactorTimeLbToW(dx);
      x2 = (x2old) + 1.0/2.0*(vx2 + vx2old)*conv->getFactorTimeLbToW(dx);
      x3 = (x3old) + 1.0/2.0*(vx3 + vx3old)*conv->getFactorTimeLbToW(dx);
      collWall(A,B,C,D,x1,x2,x3,x1old,x2old,x3old,dx,vx1,vx2,vx3,ii);
      vx1old = vx1;
      vx2old = vx2;
      vx3old = vx3;
   }
}
//////////////////////////////////////////////////////////////////////////
bool PathLineCoProcessor::checkNodes( Block3DPtr block)
{
   bool result = true;
   LBMKernelETD3Q27Ptr kernel = boost::dynamic_pointer_cast<LBMKernelETD3Q27>(block->getKernel());
   BCArray3D<D3Q27BoundaryCondition> bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();
   D3Q27BoundaryConditionPtr bcPtr;

   double x1_ch = x1;// + vx1old*conv->getfactorTimeShouldMultiplebyDx()*grid->getDeltaX(block);
   double x2_ch = x2;// + vx2old*conv->getfactorTimeShouldMultiplebyDx()*grid->getDeltaX(block);
   double x3_ch = x3;// + vx3old*conv->getfactorTimeShouldMultiplebyDx()*grid->getDeltaX(block);

   int nx1 = static_cast<int>(bcArray.getNX1());
   int nx2 = static_cast<int>(bcArray.getNX2());
   int nx3 = static_cast<int>(bcArray.getNX3());

   UbTupleInt3 nodeIndexes = grid->getNodeIndexes(block, x1_ch, x2_ch, x3_ch);

   int ix1 = val<1>(nodeIndexes);
   int ix2 = val<2>(nodeIndexes);
   int ix3 = val<3>(nodeIndexes);
   bool outlet=false;

   for (int xx3 = ix3; xx3 <= ix3 + 1 && ix3 + 1 < nx3; xx3++)
      for(int xx2 = ix2; xx2 <= ix2 + 1 && ix2 + 1 < nx2; xx2++)
         for(int xx1 = ix1; xx1 <= ix1 + 1 && ix1 + 1 < nx1; xx1++)
         {
            bcPtr=bcArray.getBC(xx1,xx2,xx3);
            for(int fdir=D3Q27System::STARTF; fdir<=D3Q27System::ENDF; fdir++)
            {
               if ( bcPtr != NULL)
               {
                  if (bcPtr->hasDensityBoundaryFlag(fdir))
                  {	
                     if(outlet)break; 
                     outlet=true;
                     break;
                  }
               }
            }
            result = result && !bcArray.isUndefined(xx1, xx2, xx3) &&  !outlet;
         }

         return result;

}
//////////////////////////////////////////////////////////////////////////
bool PathLineCoProcessor::checkNodes2( Block3DPtr block,double x11,double x22,double x33)
{
   bool result = true;
   LBMKernelETD3Q27Ptr kernel = boost::dynamic_pointer_cast<LBMKernelETD3Q27>(block->getKernel());
   BCArray3D<D3Q27BoundaryCondition> bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();
   D3Q27BoundaryConditionPtr bcPtr;

   double x1_ch = x11;// + vx1old*conv->getfactorTimeShouldMultiplebyDx()*grid->getDeltaX(block);
   double x2_ch = x22;// + vx2old*conv->getfactorTimeShouldMultiplebyDx()*grid->getDeltaX(block);
   double x3_ch = x33;// + vx3old*conv->getfactorTimeShouldMultiplebyDx()*grid->getDeltaX(block);

   int nx1 = static_cast<int>(bcArray.getNX1());
   int nx2 = static_cast<int>(bcArray.getNX2());
   int nx3 = static_cast<int>(bcArray.getNX3());

   UbTupleInt3 nodeIndexes = grid->getNodeIndexes(block, x1_ch, x2_ch, x3_ch);

   int ix1 = val<1>(nodeIndexes);
   int ix2 = val<2>(nodeIndexes);
   int ix3 = val<3>(nodeIndexes);
   bool outlet=false;

   for (int xx3 = ix3; xx3 <= ix3 + 1 && ix3 + 1 < nx3; xx3++)
      for(int xx2 = ix2; xx2 <= ix2 + 1 && ix2 + 1 < nx2; xx2++)
         for(int xx1 = ix1; xx1 <= ix1 + 1 && ix1 + 1 < nx1; xx1++)
         {
            bcPtr=bcArray.getBC(xx1,xx2,xx3);
            for(int fdir=D3Q27System::STARTF; fdir<=D3Q27System::ENDF; fdir++)
            {
               if ( bcPtr != NULL)
               {
                  if (bcPtr->hasDensityBoundaryFlag(fdir))
                  {	
                     if(outlet)break; 
                     outlet=true;
                     break;
                  }
               }
            }
            result = result && !bcArray.isUndefined(xx1, xx2, xx3) &&  !outlet;
         }

         return result;

}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::outICell(D3Q27ICell& iCell)
{
   std::ofstream ostr;
   string fname = path+"_iCell";
   ostr.open(fname.c_str(), std::ios_base::out);
   if(!ostr)
   { 
      ostr.clear();
      string path = UbSystem::getPathFromString(fname);
      if(path.size()>0){ UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out);}
      if(!ostr) throw UbException(UB_EXARGS,"couldn't open file "+fname);
   }
   ostr << "step = "<< stepcheck << endl;
   for (int i = 0; i < 27; i++)
   {
      ostr << "iCell.BNE[" << i <<"] = "<< iCell.BNE[i] << "\t"<<"iCell.BNW[" << i <<"] = "<< iCell.BNW[i] << "\t"<<"iCell.BSE[" << i <<"] = "<< iCell.BSE[i]<<"\t"<<"iCell.BSW[" << i <<"] = "<< iCell.BSW[i]<<"\t"<<"iCell.TNE[" << i <<"] = "<< iCell.TNE[i] << "\t"<<"iCell.TNW[" << i <<"] = "<< iCell.TNW[i] << "\t"<<"iCell.TSE[" << i <<"] = "<< iCell.TSE[i]<<"\t"<<"iCell.TSW[" << i <<"] = "<< iCell.TSW[i]<<"\t"<< endl;

   }
   ostr.close();
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::updateDistributedValues()
{
   int size;	
   MPI_Status status; 
   MPI_Comm_size(MPI_COMM_WORLD,&size);	
   const int row=12;
   double values[row];//
   if (this->root == this->rank)
   {
      values[0]=x1;
      values[1]=x2;
      values[2]=x3;
      values[3]=x1old;
      values[4]=x2old;
      values[5]=x3old;
      values[6]=vx1old;
      values[7]=vx2old;
      values[8]=vx3old;
      values[9]=vx1oldf;
      values[10]=vx2oldf;
      values[11]=vx3oldf;
      for (int i=0;i<size;i++)
      {
         if (i!=this->root)
         {
            MPI_Send(&values,row, MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD);
         }	
      }
   }
   else{ 
      MPI_Recv(&values,row, MPI_DOUBLE_PRECISION,this->root,this->rank,MPI_COMM_WORLD,&status);
      x1     = values[0];
      x2     = values[1];
      x3     = values[2];
      x1old  = values[3];
      x2old  = values[4];
      x3old  = values[5];
      vx1old = values[6];
      vx2old = values[7];
      vx3old = values[8];
      vx1oldf =values[9];
      vx2oldf=values[10];
      vx3oldf=values[11];
   }
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::findPlane(int ix1,int ix2,int ix3,Grid3DPtr grid,Block3DPtr block,double &A,double &B,double &C,double &D,double &ii)
{
   double x1plane=0.0,y1plane=0.0,z1plane=0.0;
   double x2plane=0.0,y2plane=0.0,z2plane=0.0;
   double x3plane=0.0,y3plane=0.0,z3plane=0.0;
   D3Q27BoundaryConditionPtr bcPtr;
   double dx = grid->getDeltaX(block);
   LBMKernelETD3Q27Ptr kernel = boost::dynamic_pointer_cast<LBMKernelETD3Q27>(block->getKernel());
   DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();
   BCArray3D<D3Q27BoundaryCondition> bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();

   for (int k = ix3; k <= ix3 + 1; k++){
      for(int j = ix2; j <= ix2 + 1; j++){
         for(int i = ix1; i <= ix1 + 1; i++)
         {
            UbTupleDouble3 pointplane1 =  grid->getNodeCoordinates(block,  i,	  j,	 k);

            double   iph=val<1>(pointplane1);
            double   jph=val<2>(pointplane1);
            double   kph=val<3>(pointplane1);

            if(!bcArray.isSolid(i, j, k))
            {
               bcPtr=bcArray.getBC(i,j,k);
               if(bcPtr)
               {	 
                  for(int fdir=D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
                  {
                     if( ii<=2)
                     {
                        LBMReal q = bcPtr->getQ(fdir);
                        if (q!=999.00000)
                        {

                           if     ( fdir==D3Q27System::E ) 
                           {
                              if (i+q<=ix1+1)
                              {
                                 if      (ii==0)	    {	x1plane=iph+q*dx;	  y1plane=jph;	 z1plane=kph;    	ii++;	} 
                                 else if (ii==1) 	{	x2plane=iph+q*dx;     y2plane=jph;	 z2plane=kph;	    if  (x1plane!=x2plane||y1plane!=y2plane||z1plane!=z2plane) 	ii++;	}
                                 else if(ii==2) 	    {   x3plane=iph+q*dx;     y3plane=jph;	 z3plane=kph;       if ((x3plane!=x1plane||y3plane!=y1plane||z3plane!=z1plane)&&(x2plane!=x3plane||y2plane!=y3plane||z2plane!=z3plane)) ii++;} 
                              }
                           }
                           if     ( fdir==D3Q27System::W ) 
                           {
                              if (i-q>=ix1)
                              {
                                 if      (ii==0)	    {	x1plane=iph-q*dx;	  y1plane=jph;	 z1plane=kph;    	ii++;	} 
                                 else if (ii==1) 	    {	x2plane=iph-q*dx;    y2plane=jph;	 z2plane=kph;	    if  (x1plane!=x2plane||y1plane!=y2plane||z1plane!=z2plane) 	ii++;	}
                                 else if(ii==2) 	    {   x3plane=iph-q*dx;    y3plane=jph;	 z3plane=kph;       if ((x3plane!=x1plane||y3plane!=y1plane||z3plane!=z1plane)&&(x2plane!=x3plane||y2plane!=y3plane||z2plane!=z3plane)) ii++;} 
                              }
                           }
                           if     ( fdir==D3Q27System::N ) 
                           {
                              if(j+q<=ix2+1)
                              {
                                 if      (ii==0)	    {	x1plane=iph;	y1plane=jph+q*dx;	 z1plane=kph; 	    ii++;	} 
                                 else if (ii==1) 	    {	x2plane=iph;    y2plane=jph+q*dx;	 z2plane=kph;	    if  (x1plane!=x2plane||y1plane!=y2plane||z1plane!=z2plane) 	ii++;	}
                                 else if (ii==2) 	    {   x3plane=iph;    y3plane=jph+q*dx;	 z3plane=kph;       if ((x3plane!=x1plane||y3plane!=y1plane||z3plane!=z1plane)&&(x2plane!=x3plane||y2plane!=y3plane||z2plane!=z3plane)) ii++;} 
                              }
                           }
                           if     ( fdir==D3Q27System::S ) 
                           {
                              if (j-q>=ix2)
                              {
                                 if      (ii==0)	    {	x1plane=iph;	y1plane=jph-q*dx;	 z1plane=kph; 	ii++;	} 
                                 else if (ii==1) 	    {	x2plane=iph;    y2plane=jph-q*dx;	 z2plane=kph;	if  (x1plane!=x2plane||y1plane!=y2plane||z1plane!=z2plane) 	ii++;	}
                                 else if (ii==2) 	    {   x3plane=iph;    y3plane=jph-q*dx;	 z3plane=kph;   if ((x3plane!=x1plane||y3plane!=y1plane||z3plane!=z1plane)&&(x2plane!=x3plane||y2plane!=y3plane||z2plane!=z3plane)) ii++;} 
                              }
                           }

                           if     ( fdir==D3Q27System::T ) 
                           {
                              if(k+q<=ix3+1)
                              {
                                 if      (ii==0)	    {	x1plane=iph;	y1plane=jph;	 z1plane=kph+q*dx; 	    ii++;	} 
                                 else if (ii==1) 	    {	x2plane=iph;    y2plane=jph;	 z2plane=kph+q*dx;	    if  (x1plane!=x2plane||y1plane!=y2plane||z1plane!=z2plane) 	ii++;	}
                                 else if (ii==2) 	    {   x3plane=iph;    y3plane=jph;	 z3plane=kph+q*dx;      if ((x3plane!=x1plane||y3plane!=y1plane||z3plane!=z1plane)&&(x2plane!=x3plane||y2plane!=y3plane||z2plane!=z3plane)) ii++;} 
                              }
                           }
                           if     ( fdir==D3Q27System::B ) 
                           {
                              if (k-q>=ix3)
                              {
                                 if      (ii==0)	    {	x1plane=iph;	y1plane=jph;	 z1plane=kph-q*dx; 	ii++;	} 
                                 else if (ii==1) 	    {	x2plane=iph;    y2plane=jph;	 z2plane=kph-q*dx;	  if  (x1plane!=x2plane||y1plane!=y2plane||z1plane!=z2plane) 	ii++;	}
                                 else if (ii==2) 	    {   x3plane=iph;    y3plane=jph;	 z3plane=kph-q*dx;     if ((x3plane!=x1plane||y3plane!=y1plane||z3plane!=z1plane)&&(x2plane!=x3plane||y2plane!=y3plane||z2plane!=z3plane)) ii++;} 
                              }
                           }

                        }
                     }
                  }

               }

            }
         }
      }
   }

   A =   y1plane* (z2plane - z3plane) + y2plane*(z3plane - z1plane) + y3plane* (z1plane - z2plane);   
   B =   z1plane* (x2plane - x3plane) + z2plane*(x3plane - x1plane) + z3plane* (x1plane - x2plane) ;      
   C =   x1plane* (y2plane - y3plane) + x2plane*(y3plane - y1plane) + x3plane* (y1plane - y2plane) ;       
   D =-( x1plane*(y2plane*z3plane - y3plane*z2plane)+x2plane*(y3plane*z1plane - y1plane*z3plane) + x3plane* (y1plane* z2plane - y2plane* z1plane));
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::finddisPointToCentercell(double x1,double x2,double x3,double xp000, double yp000, double zp000,double &dis,double &dx,
                                                          double &deltaT,Block3DPtr &block,int level)
{	
   UbTupleInt3 blockIndexes = grid->getBlockIndexes(xp000, yp000,  zp000,level);
   block = grid->getBlock(val<1>(blockIndexes), val<2>(blockIndexes), val<3>(blockIndexes), level);
   LBMKernelETD3Q27Ptr kernel = boost::dynamic_pointer_cast<LBMKernelETD3Q27>(block->getKernel());
   if (!kernel)	  cout<<__FILE__<<" "<<__LINE__<<endl;
   else
   {
      DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();
      BCArray3D<D3Q27BoundaryCondition> bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();
      if(block && block->isActive())
      {	     		     
         deltaT= 1.0/(double)(1<<level);
         dx = grid->getDeltaX(block);
         UbTupleInt3 nodeIndexes = grid->getNodeIndexes(block,  xp000, yp000,  zp000);
         int ix1Se = val<1>(nodeIndexes);
         int ix2Se = val<2>(nodeIndexes);
         int ix3Se = val<3>(nodeIndexes);  
         if     (!iProcessor->iCellHasSolid(bcArray, ix1Se, ix2Se,  ix3Se  ))        
         {				 
            UbTupleDouble3 pointCell =  grid->getNodeCoordinates(block,  ix1Se,	  ix2Se,	 ix3Se);	 double   iph=val<1>(pointCell)+dx/2.0;	 double   jph=val<2>(pointCell)+dx/2.0;	 double   kph=val<3>(pointCell)+dx/2.0;
            dis=sqrt((iph-x1)*(iph-x1)+(jph-x2)*(jph-x2)+(kph-x3)*(kph-x3));
         }
      }	 	 	
   }
}
//////////////////////////////////////////////////////////////////////////
//void D3Q27PathLinePostprocessor::finddisPointToCentercell(double x1,double x2,double x3,double xp000, double yp000, double zp000,double &dis,double &dx,
// double &deltaT,Block3DPtr &block)
//{

// int minInitLevel = this->grid->getCoarsestInitializedLevel();
// int maxInitLevel = this->grid->getFinestInitializedLevel();
// bool goinside=true;
// Block3DPtr block2;
// for(int level = minInitLevel; level<=maxInitLevel; level++)
// {	   
//	 if (goinside )
//	 {
//		 UbTupleInt3 blockIndexes = grid->getBlockIndexes(xp000, yp000,  zp000,level);
//		 block = grid->getBlock(val<1>(blockIndexes), val<2>(blockIndexes), val<3>(blockIndexes), level);
//		 LBMKernelETD3Q27Ptr kernel = boost::dynamic_pointer_cast<LBMKernelETD3Q27>(block->getKernel());
//		 if (kernel)
//		 {
//			 DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();
//			 BCArray3D<D3Q27BoundaryCondition> bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();
//			 if(block && block->isActive())
//			 {	     		     
//				 if(!checkNodes(block))continue;
//				 goinside=false;
//				 deltaT= 1.0/(double)(1<<level);
//				 dx = grid->getDeltaX(block);
//				 UbTupleInt3 nodeIndexes = grid->getNodeIndexes(block,  xp000, yp000,  zp000);
//				 int ix1Se = val<1>(nodeIndexes);
//				 int ix2Se = val<2>(nodeIndexes);
//				 int ix3Se = val<3>(nodeIndexes);  
//				 if     (!iProcessor->iCellHasSolid(bcArray, ix1Se, ix2Se,  ix3Se  ))        
//				 {				 
//					 UbTupleDouble3 pointCell =  grid->getNodeCoordinates(block,  ix1Se,	  ix2Se,	 ix3Se);	 double   iph=val<1>(pointCell)+dx/2.0;	 double   jph=val<2>(pointCell)+dx/2.0;	 double   kph=val<3>(pointCell)+dx/2.0;
//					 dis=sqrt((iph-x1)*(iph-x1)+(jph-x2)*(jph-x2)+(kph-x3)*(kph-x3));
//					 block2=block;
//				 }
//			 }	 	 
//		 }
//	 }
// }

//}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::getWallPosition(double input[])
{ 

   //{0=xp000, 1=ix2ph,  2=ix3ph,3=x1,4=x2,5=x3,6=Xw,7=Yw,8=Zw,9=A,10=B,11=C,12=D,13=normalDis,14=di,15=dj,16=dk,
   //17=cofWeightx,18=cofWeighty,19=cofWeightz,,20=dx,21=omega,22=nue,23=disx,24=dir  };
   double dx=input[20];double xp000=input[0];double yp000=input[1];double zp000=input[2];double Xw=input[6];double Yw=input[7];double Zw=input[8];
   double x1=input[3] ;double x2=input[4] ;double x3=input[5]; double A=input[9];double B=input[10];double C=input[11];double D=input[12];
   double normalDis=input[13]; double di=input[14];double dj=input[15];double dk=input[16]; double cofWeightx=input[17];double cofWeighty=input[18];
   double cofWeightz=input[19];int dir=(int)input[24];

   //char dirx,diry,dirz,dirxy,dirxz,diryz,dirxyz;
   if      (dir==dirx  )  {  Xw=(B*x2+C*x3+D)/(-A); Yw=x2;                 Zw=x3;	         	      cofWeightx=1;  cofWeighty=0;cofWeightz=0;}     
   else if (dir==diry  )  {  Xw=x1;                 Yw=(A*x1+C*x3+D)/(-B); Zw=x3;		              cofWeightx=0;  cofWeighty=1;cofWeightz=0;}
   else if (dir==dirz  )  {  Xw=x1;                 Yw=x3;	             Zw=(B*x2+A*x1+D)/(-C);   cofWeightx=0;  cofWeighty=0;cofWeightz=1;    }
   else if (dir==dirxy )  {  Xw=x1-di*normalDis;    Yw=x2-dj*normalDis;    Zw=x3;        		      cofWeightx=1;  cofWeighty=1;cofWeightz=0;}
   else if (dir==dirxz )  {  Xw=x1-di*normalDis;    Yw=x2;                 Zw=x3-dk*normalDis;	  cofWeightx=1;  cofWeighty=0;cofWeightz=1;    }
   else if (dir==diryz )  {  Xw=x1;                 Yw=x2-dj*normalDis;    Zw=x3-dk*normalDis;	  cofWeightx=0;  cofWeighty=1;cofWeightz=1;    }
   else if (dir==dirxyz)  {  Xw=x1-di*normalDis;    Yw=x2-dj*normalDis;    Zw=x3-dk*normalDis;	  cofWeightx=1;  cofWeighty=1;cofWeightz=1;    }  

   Xw=(Xw-xp000)/dx-0.5;       
   Yw=(Yw-yp000)/dx-0.5;       
   Zw=(Zw-zp000)/dx-0.5; 

   input[6]=Xw;
   input[7]=Yw;
   input[8]=Zw;
   input[17]=cofWeightx;
   input[18]=cofWeighty;
   input[19]=cofWeightz;

}
////////////////////////////////////////////////////////////////////////////
// void D3Q27PathLinePostprocessor::getWallPosition(double dixmax,double disx,double disy,double disz,double disxy,double disxz,double disyz,double disxyz
//  ,double &Xw,double &Yw,double &Zw,double x1,double x2,double x3,double A,double B,double C,double D,double normalDis,double di,double dj,double dk,
//  double &cofWeightx,double &cofWeighty,double &cofWeightz)
// {
//  if      (dixmax==disx)  {  Xw=(B*x2+C*x3+D)/(-A); Yw=x2;                 Zw=x3;	         	  cofWeightx=1;  cofWeighty=0;cofWeightz=0;}     
//  else if (dixmax==disy)  {  Xw=x1;                 Yw=(A*x1+C*x3+D)/(-B); Zw=x3;		          cofWeightx=0;  cofWeighty=1;cofWeightz=0;}
//  else if (dixmax==disz)  {  Xw=x1;                 Yw=x3;	               Zw=(B*x2+A*x1+D)/(-C); cofWeightx=0;  cofWeighty=0;cofWeightz=1;}
//  else if (dixmax==disxy) {  Xw=x1-di*normalDis;    Yw=x2-dj*normalDis;    Zw=x3;        		  cofWeightx=1;  cofWeighty=1;cofWeightz=0;}
//  else if (dixmax==disxz) {  Xw=x1-di*normalDis;    Yw=x2;                 Zw=x3-dk*normalDis;	  cofWeightx=1;  cofWeighty=0;cofWeightz=1;}
//  else if (dixmax==disyz) {  Xw=x1;                 Yw=x2-dj*normalDis;    Zw=x3-dk*normalDis;	  cofWeightx=0;  cofWeighty=1;cofWeightz=1;}
//  else if (dixmax==disxyz){  Xw=x1-di*normalDis;    Yw=x2-dj*normalDis;    Zw=x3-dk*normalDis;	  cofWeightx=1;  cofWeighty=1;cofWeightz=1;}
// }
////////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::findInitialCell(double s,double di,double dj,double dk,double dx,double ix1ph,double ix2ph,double ix3ph,double &xp000,double &yp000,double &zp000)
{
   if (s*di>0)         {            xp000=ix1ph+dx;    }
   else if(s*di<0)     {            xp000=ix1ph-dx;    }
   else                {            xp000=ix1ph;      }
   if (s*dj>0)         {            yp000=ix2ph+dx;    }
   else if(s*dj<0)     {            yp000=ix2ph-dx;    }
   else                {            yp000=ix2ph;      }
   if (s*dk>0)         {            zp000=ix3ph+dx;    }
   else if(s*dk<0)     {            zp000=ix3ph-dx;    }
   else                {            zp000=ix3ph;      }
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::rearangedDouble( double dis[7])
{
   double item;
   for (int i=6;i>0;i--){   
      for(int j=0;j<i;j++){ 
         if (dis[j]>dis[j+1])
         {	  
            item=dis[j];  
            dis[j]=dis[j+1]; 
            dis[j+1]=item;					   
         }
      }
   } 
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::rearangeInt( int dis[7])
{
   double item;
   for (int i=6;i>0;i--){   
      for(int j=0;j<i;j++){ 
         if (dis[j]>dis[j+1])
         {	  
            item=dis[j];  
            dis[j]=dis[j+1]; 
            dis[j+1]=(int)item;					   
         }
      }
   } 
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::addCorrection(double &vx1,double &vx2,double &vx3,double &tauxx,double &tauyy,double &tauzz,double &tauxy,double &tauxz,double &tauyz,
                                               double dx,D3Q27ICell iCell,double Xw, double Yw, double Zw,double  omega,double cofWeightx,double cofWeighty,double cofWeightz,double ii,double x1LB,double x2LB,double x3LB,double di,double dj,double dk)
{
   LBMReal f[D3Q27System::ENDF+1];
   double  Udif,Vdif,Wdif;	  
   iHelper->interpolate8to1WithVelocity(iCell, Xw, Yw, Zw, omega,Udif,Vdif,Wdif);
   double Ucx,Ucy,Ucz;
   double Vcx,Vcy,Vcz;
   double Wcx,Wcy,Wcz;

   double dxUcx,dyUcy,dzUcz;
   double dxVcx,dyVcy,dzVcz;
   double dxWcx,dyWcy,dzWcz;
   //double di,dj,dk;

   if (cofWeightx==0){di=0;}
   if (cofWeighty==0){dj=0;}
   if (cofWeightz==0){dk=0;}

   double weightx=(di*di)/(di*di+dj*dj+dk*dk);
   double weighty=(dj*dj)/(di*di+dj*dj+dk*dk);
   double weightz=(dk*dk)/(di*di+dj*dj+dk*dk);
   // /* if (Xw<0.5||Xw>-0.5){weightx=0;}
   //  if (Yw<0.5||Yw>-0.5){weighty=0;}
   //  if (Zw<0.5||Zw>-0.5){weightz=0;}
   //*/
   if (weightx==0){Ucx=0;Vcx=0;Wcx=0;dxUcx=0;dxVcx=0;dxWcx=0;}else{
      Ucx=-Udif*weightx/(Xw*Xw*Xw-Xw/4.0)*(x1LB*x1LB*x1LB-x1LB/4.0);Vcx=-Vdif*weightx/(Xw*Xw*Xw-Xw/4.0)*(x1LB*x1LB*x1LB-x1LB/4.0);Wcx=-Wdif*weightx/(Xw*Xw*Xw-Xw/4.0)*(x1LB*x1LB*x1LB-x1LB/4.0);
      dxUcx=-Udif*weightx/(Xw*Xw*Xw-Xw/4.0)*(3.0*x1LB*x1LB-1.0/4.0);dxVcx=-Vdif*weightx/(Xw*Xw*Xw-Xw/4.0)*(3.0*x1LB*x1LB-1.0/4.0);dxWcx=-Wdif*weightx/(Xw*Xw*Xw-Xw/4.0)*(3.0*x1LB*x1LB-1.0/4.0);
   }
   if (weighty==0){Ucy=0;Vcy=0;Wcy=0;dyUcy=0;dyVcy=0;dyWcy=0;}else{
      Ucy=-Udif*weighty/(Yw*Yw*Yw-Yw/4.0)*(x2LB*x2LB*x2LB-x2LB/4.0);Vcy=-Vdif*weighty/(Yw*Yw*Yw-Yw/4.0)*(x2LB*x2LB*x2LB-x2LB/4.0);Wcy=-Wdif*weighty/(Yw*Yw*Yw-Yw/4.0)*(x2LB*x2LB*x2LB-x2LB/4.0);
      dyUcy=-Udif*weighty/(Yw*Yw*Yw-Yw/4.0)*(3.0*x2LB*x2LB-1.0/4.0);dyVcy=-Vdif*weighty/(Yw*Yw*Yw-Yw/4.0)*(3.0*x2LB*x2LB-1.0/4.0);dyWcy=-Wdif*weighty/(Yw*Yw*Yw-Yw/4.0)*(3.0*x2LB*x2LB-1.0/4.0);
   }
   if (weightz==0){Ucz=0;Vcz=0;Wcz=0;dzUcz=0;dzVcz=0;dzWcz=0;}else{
      Ucz=-Udif*weightz/(Zw*Zw*Zw-Zw/4.0)*(x3LB*x3LB*x3LB-x3LB/4.0);Vcz=-Vdif*weightz/(Zw*Zw*Zw-Zw/4.0)*(x3LB*x3LB*x3LB-x3LB/4.0);Wcz=-Wdif*weightz/(Zw*Zw*Zw-Zw/4.0)*(x3LB*x3LB*x3LB-x3LB/4.0);
      dzUcz=-Udif*weightz/(Zw*Zw*Zw-Zw/4.0)*(3.0*x3LB*x3LB-1.0/4.0);dzVcz=-Vdif*weightz/(Zw*Zw*Zw-Zw/4.0)*(3.0*x3LB*x3LB-1.0/4.0);dzWcz=-Wdif*weightz/(Zw*Zw*Zw-Zw/4.0)*(3.0*x3LB*x3LB-1.0/4.0);
   }

   double Ucor=Ucx+Ucy+Ucz;
   double Vcor=Vcx+Vcy+Vcz;
   double Wcor=Wcx+Wcy+Wcz;

   double tauxxC=dxUcx;
   double tauyyC=dyVcy;
   double tauzzC=dzWcz;
   double tauxyC=0.5*(dyUcy+dxVcx);
   double tauxzC=0.5*(dzUcz+dxWcx);
   double tauyzC=0.5*(dzVcz+dyWcy);

   if (ii!=3) 
   {
      //std::cout<<"there are not 3point for making plane"<<endl<<"Ucor="<<Ucor<<endl<<"vx1="<<vx1;
      Ucor  =0.0;
      Vcor  =0.0;
      Wcor  =0.0;
      tauxxC=0.0;
      tauyyC=0.0;
      tauzzC=0.0;
      tauxyC=0.0;
      tauxzC=0.0;
      tauyzC=0.0;
   }

   vx1+=Ucor;
   vx2+=Vcor;
   vx3+=Wcor;
   tauxx+=tauxxC;
   tauyy+=tauyyC;
   tauzz+=tauzzC;
   tauxy+=tauxyC;
   tauxz+=tauxzC;
   tauyz+=tauyzC;

   vx1 = vx1 *  conv->getFactorVelocityLbToW();
   vx2 = vx2 *  conv->getFactorVelocityLbToW();
   vx3 = vx3 *  conv->getFactorVelocityLbToW();
   tauxx=tauxx/dx;
   tauyy=tauyy/dx;
   tauzz=tauzz/dx;
   tauxy=tauxy/dx;
   tauxz=tauxz/dx;
   tauyz=tauyz/dx;
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::collWall(double A,double B,double C,double D,double &x1,double &x2,double &x3,double x1old,double x2old,double x3old,double dx,double &vx1,double &vx2,double &vx3,double ii)
{
   if (particleHasMass)
   {
      double s = A*x1 + B*x2 + C*x3 + D;//The sign of s = Ax + By + Cz + D determines which side the point (x,y,z) lies with respect to the plane. If s > 0 then the point lies on the same side as the normal (A,B,C). If s < 0 then it lies on the opposite side, if s = 0 then the point (x,y,z) lies on the plane.
      if (s>0){s=1;} else if (s<0){s=-1;}else {s=0;}
      double normalDis=((A*x1 + B*x2 + C*x3 + D)/sqrt(A*A+B*B+C*C));///distance point to plane xp-Xw=distance
      double di=A/sqrt(A*A+B*B+C*C);      double dj=B/sqrt(A*A+B*B+C*C);       double dk=C/sqrt(A*A+B*B+C*C);
      A*=s;B*=s;C*=s;
      if (abs(normalDis)<0.05*dx)
      {
         if (ii==3)
         {
            ///reflection without e factor
            /*  double vn=(A*vx1 + B*vx2 + C*vx3)/(A*A+B*B+C*C);
            double  vx1L=vx1-2*vn*A;
            double  vx2L=vx2-2*vn*B;
            double  vx3L=vx3-2*vn*C;*/
            double vnx=A*(A*vx1 + B*vx2 + C*vx3)/(A*A+B*B+C*C);
            double vny=B*(A*vx1 + B*vx2 + C*vx3)/(A*A+B*B+C*C);
            double vnz=C*(A*vx1 + B*vx2 + C*vx3)/(A*A+B*B+C*C);
            ////do collision  
            double CollitionFactor_n=0.01;
            double CollitionFactor_bt=0.01;
            vnx=-1*CollitionFactor_n*vnx;
            vny=-1*CollitionFactor_n*vny;
            vnz=-1*CollitionFactor_n*vnz;
            double vbtx=(B*B*vx1 + C*C*vx1 - A*B*vx2 - A*C*vx3)/(A*A+B*B+C*C);
            double vbty=(-(A*B*vx1) + A*A*vx2 + C*C*vx2 - B*C*vx3)/(A*A+B*B+C*C);
            double vbtz=(-(A*C*vx1) - B*C*vx2 + A*A*vx3 + B*B*vx3)/(A*A+B*B+C*C);
            vbtx=vbtx*CollitionFactor_bt;
            vbty=vbty*CollitionFactor_bt;
            vbtz=vbtz*CollitionFactor_bt;
            vx1=vnx+vbtx;
            vx2=vny+vbty;
            vx3=vnz+vbtz;

            x1 = (x1old);
            x2 = (x2old);
            x3 = (x3old);
         } 
         else   std::cout<<"there is no plane to reflect";
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::CalcVelParticle(double dx,double &vx1,double &vx2,double &vx3,double vx1old,double vx2old,
                                                 double vx3old,double &vx1oldf,double &vx2oldf,double &vx3oldf)
{
   double signDragx,signDragy,signDragz;
   double dragCoff=1.0;
   if     (vx1old>=vx1oldf)    {  signDragx=-1.0;   } 
   else if(vx1old<vx1oldf)     {  signDragx= 1.0;    }
   if     (vx2old>=vx2oldf)    {  signDragy=-1.0;   } 
   else if(vx2old< vx2oldf)    {  signDragy= 1.0;    }
   if     (vx3old>=vx3oldf)    {  signDragz=-1.0;   } 
   else if(vx3old< vx3oldf)    {  signDragz= 1.0;    }
   double vx1p=vx1old+signDragx*(vx1old-vx1oldf)*(vx1old-vx1oldf)*conv->getFactorTimeLbToW(dx);
   double vx2p=vx2old+signDragy*(vx2old-vx2oldf)*(vx2old-vx2oldf)*conv->getFactorTimeLbToW(dx);
   double vx3p=vx3old+signDragz*(vx3old-vx3oldf)*(vx3old-vx3oldf)*conv->getFactorTimeLbToW(dx);   

   double signDragNextx,signDragNexty,signDragNextz;
   if     (vx1p>=vx1)    {  signDragNextx=-1.0;   } 
   else if(vx1p< vx1)    {  signDragNextx= 1.0;    }
   if     (vx2p>=vx2)    {  signDragNexty=-1.0;   } 
   else if(vx2p< vx2)    {  signDragNexty= 1.0;    }
   if     (vx3p>=vx3)    {  signDragNextz=-1.0;   } 
   else if(vx3p< vx3)    {  signDragNextz= 1.0;    }
   ////////////////////velocity of particle////////////////////////////////////////////////////// 
   double velx1oldf=vx1;
   double velx2oldf=vx2;
   double velx3oldf=vx3;

   vx1=vx1old+1.0/2.0*( signDragNextx*(vx1p-vx1)*(vx1p-vx1)+signDragx*(vx1old-vx1oldf)*(vx1old-vx1oldf) )*conv->getFactorTimeLbToW(dx);
   vx2=vx2old+1.0/2.0*( signDragNexty*(vx2p-vx2)*(vx2p-vx2)+signDragy*(vx2old-vx2oldf)*(vx2old-vx2oldf) )*conv->getFactorTimeLbToW(dx);
   vx3=vx3old+1.0/2.0*( signDragNextz*(vx3p-vx3)*(vx3p-vx3)+signDragz*(vx3old-vx3oldf)*(vx3old-vx3oldf) )*conv->getFactorTimeLbToW(dx);

   vx1oldf=velx1oldf;
   vx2oldf=velx2oldf;
   vx3oldf=velx3oldf;
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::CalcVelParticle2(double dx,double &vx1,double &vx2,double &vx3,double vx1old,double vx2old,
                                                  double vx3old,double &vx1oldf,double &vx2oldf,double &vx3oldf)
{
   double signDragx,signDragy,signDragz;

   LBMReal muRE = 1.002*1e-3;//m2/s
   double mass=2.3*1e-17;// kg
   double Diameter=600*1e-9;//m
   double deltatime=conv->getFactorTimeLbToW(dx)*0.001;
   double Coff=3*PI*Diameter*muRE*deltatime/mass;
   double exCoff=exp(-Coff);


   //////////////////////velocity of particle////////////////////////////////////////////////////// 
   double velx1oldf=vx1;
   double velx2oldf=vx2;
   double velx3oldf=vx3;

   //air

   vx1=vx1old*exCoff+(velx1oldf-(velx1oldf)*exCoff);
   vx2=vx2old*exCoff+(velx2oldf-(velx2oldf)*exCoff);
   vx3=vx3old*exCoff+(velx3oldf-(velx3oldf)*exCoff);


   vx1oldf=velx1oldf;
   vx2oldf=velx2oldf;
   vx3oldf=velx3oldf;
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::getRankWithPosition(Grid3DPtr grid,double xp000, double yp000, double zp000,int &RankNeigber,int level)
{
   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   int maxInitLevel = this->grid->getFinestInitializedLevel();
   bool goinside=true;

   LBMKernelETD3Q27Ptr kernel;
   DistributionArray3DPtr distributions;
   BCArray3D<D3Q27BoundaryCondition> bcArray;
   Block3DPtr block;

   //for(int level = minInitLevel; level<=maxInitLevel; level++)
   {	   
      UbTupleInt3 blockIndexes = grid->getBlockIndexes(xp000,yp000, zp000,level);
      block = grid->getBlock(val<1>(blockIndexes), val<2>(blockIndexes), val<3>(blockIndexes), level);
      //LBMKernelETD3Q27Ptr kernel = boost::dynamic_pointer_cast<LBMKernelETD3Q27>(block->getKernel());
      // if(block )
      {
         if(block->isActive())
         {	     		     
            // if(!checkNodes2(block,xp000,yp000,zp000))continue;
            RankNeigber= block->getRank();   
            //break;
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::getRankWithPosition2(Grid3DPtr grid,double xp000, double yp000, double zp000,int &RankNeigber)
{
   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   int maxInitLevel = this->grid->getFinestInitializedLevel();
   bool goinside=true;

   LBMKernelETD3Q27Ptr kernel;
   DistributionArray3DPtr distributions;
   BCArray3D<D3Q27BoundaryCondition> bcArray;
   Block3DPtr block;

   for(int level = minInitLevel; level<=maxInitLevel; level++)
   {	   
      UbTupleInt3 blockIndexes = grid->getBlockIndexes(xp000,yp000, zp000,level);
      block = grid->getBlock(val<1>(blockIndexes), val<2>(blockIndexes), val<3>(blockIndexes), level);
      LBMKernelETD3Q27Ptr kernel = boost::dynamic_pointer_cast<LBMKernelETD3Q27>(block->getKernel());
      // if(kernel )
      {
         if(block->isActive())
         {	     		     
            // if(!checkNodes2(block,xp000,yp000,zp000))continue;
            RankNeigber= block->getRank();   
            break;
         }
      }
   }
}


////////////////////////////////////////////////////////////////////////// 
void PathLineCoProcessor::getAllTag(int allRank[7],int allTag[7])
{

   for (int i=0;i<7;i++)
   {
      int tagme=0;
      for (int j=0;j<=i;j++)
      {
         if(allRank[j]==allRank[i])
         {
            tagme++;
         }
         allTag[i]=tagme;
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::getAllInfoForCompare(double input[],double *Vect)
{
   //{0=xp000, 1=ix2ph,  2=ix3ph,3=x1,4=x2,5=x3,6=Xw,7=Yw,8=Zw,9=A,10=B,11=C,12=D,13=normalDis,14=di,15=dj,16=dk,
   //17=cofWeightx,18=cofWeighty,19=cofWeightz,,20=dx,21=omega,22=level,23=disx,24=dir  };

   double deltaT;
   Block3DPtr block;
   finddisPointToCentercell(input[3],input[4],input[5],input[0],input[1],input[2],input[23],input[20],deltaT,block,(int)input[22]);//get double &dis,double &dx,double &deltaT
   getWallPosition(input);	 //get wall position,Xw,Yw,Zw and cofWeightx,cofWeighty,cofWeightz	   

   double xp000=input[0];
   double yp000=input[1];
   double zp000=input[2]; 
   // double nue=input[22];
   LBMReal omega = LBMSystem::calcCollisionFactor(nue,block->getLevel());
   input[21]=omega;
   UbTupleInt3 nodeIndexes000 = grid->getNodeIndexes(block,  xp000, yp000,  zp000);
   double xp000id = val<1>(nodeIndexes000);
   double yp000id = val<2>(nodeIndexes000);
   double zp000id = val<3>(nodeIndexes000);
   LBMKernelETD3Q27Ptr kernel = boost::dynamic_pointer_cast<LBMKernelETD3Q27>(block->getKernel());
   if (kernel)
   {
      DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();
      D3Q27ICell iCell;
      iProcessor->readICell(distributions, iCell, (int)xp000id, (int)yp000id, (int)zp000id);
      getVectfromiCell(iCell,Vect);
   } 
   else
   {
      cout<<__FILE__<<" "<<__LINE__<<endl;
      for (int k=0;k<8*27;k++)       {  Vect[k]=999.0;   	  }
   }


}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::sendMaxtagToAllProccess(int allRank[],int allTag[],int numberpro,int &numberproccessused,int *mainproccess)
{

   int allRank_[7];
   for (int i=0;i<7;i++)
   {
      allRank_[i]=allRank[i];
   }
   rearangeInt(allRank_);

   for (int i=0;i<7;i++)
   {
      int maxtag_=0;
      if (i==6)
      {
         if (allRank_[i]!=this->root)
         {			   
            for (int k=0;k<=i;k++)
            {
               if (allRank_[i]==allRank_[k])
               {
                  maxtag_++;
               }
            }
            MPI_Send(&maxtag_,1, MPI_INT,allRank_[i],i,MPI_COMM_WORLD);
         }
         numberproccessused++;
         mainproccess[numberproccessused-1]=allRank_[i];
      } 
      else 
      {		
         if(allRank_[i]!=allRank_[i+1])			   
         {				   
            if (allRank_[i]!=this->root)
            { 
               for (int k=0;k<=i;k++)
               {
                  if (allRank_[i]==allRank_[k])
                  {
                     maxtag_++;
                  }
               }
               MPI_Send(&maxtag_,1, MPI_INT,allRank_[i],i,MPI_COMM_WORLD);  
            }
            numberproccessused++;
            mainproccess[numberproccessused-1]=allRank_[i];
         }   
      }
   }
   for (int i=0;i<numberpro;i++)
   {
      bool goinside=true;
      for (int j=0;j<numberproccessused;j++)
      {
         if (i==mainproccess[j]){goinside=false;break;}
      }
      if(goinside)
      {
         int maxtag_=0;
         if (i!=this->root)
         {
            MPI_Send(&maxtag_,1, MPI_INT,i,i,MPI_COMM_WORLD);
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::getAllRankWithAllPositions(double ix1ph,double ix2ph,double ix3ph,double xp000,double yp000,double zp000,
                                                            int &Rankx,int &Ranky,int &Rankz,int &Rankxy,int &Rankxz,int &Rankyz,int &Rankxyz,int level)
{

   getRankWithPosition(grid,xp000, ix2ph,  ix3ph,Rankx  ,level);
   getRankWithPosition(grid,ix1ph, yp000,  ix3ph,Ranky  ,level);
   getRankWithPosition(grid,ix1ph, ix2ph,  zp000,Rankz  ,level);
   getRankWithPosition(grid,xp000, yp000,  ix3ph,Rankxy ,level);
   getRankWithPosition(grid,xp000, ix2ph,  zp000,Rankxz ,level);
   getRankWithPosition(grid,ix1ph, yp000,  zp000,Rankyz ,level);
   getRankWithPosition(grid,xp000, yp000,  zp000,Rankxyz,level);  	

}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::getAllRankWithAllPositions2(double ix1ph,double ix2ph,double ix3ph,double xp000,double yp000,double zp000,
                                                             int &Rankx,int &Ranky,int &Rankz,int &Rankxy,int &Rankxz,int &Rankyz,int &Rankxyz)
{

   getRankWithPosition2(grid,xp000, ix2ph,  ix3ph,Rankx  );
   getRankWithPosition2(grid,ix1ph, yp000,  ix3ph,Ranky  );
   getRankWithPosition2(grid,ix1ph, ix2ph,  zp000,Rankz  );
   getRankWithPosition2(grid,xp000, yp000,  ix3ph,Rankxy );
   getRankWithPosition2(grid,xp000, ix2ph,  zp000,Rankxz );
   getRankWithPosition2(grid,ix1ph, yp000,  zp000,Rankyz );
   getRankWithPosition2(grid,xp000, yp000,  zp000,Rankxyz);  	

}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::getVectfromiCell(D3Q27ICell iCell,double *Vect)
{

   for (int i=0*27;i<1*27;i++)         {   Vect[i]=iCell.TSW[i-0*27];   }
   for (int i=1*27;i<2*27;i++)		   {   Vect[i]=iCell.TNW[i-1*27];   }
   for (int i=2*27;i<3*27;i++)		   {   Vect[i]=iCell.TNE[i-2*27];   }
   for (int i=3*27;i<4*27;i++)         {   Vect[i]=iCell.TSE[i-3*27];   }
   for (int i=4*27;i<5*27;i++)		   {   Vect[i]=iCell.BSW[i-4*27];   }
   for (int i=5*27;i<6*27;i++)		   {   Vect[i]=iCell.BNW[i-5*27];   }
   for (int i=6*27;i<7*27;i++)		   {   Vect[i]=iCell.BNE[i-6*27];   }
   for (int i=7*27;i<8*27;i++)		   {   Vect[i]=iCell.BSE[i-7*27];   }
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::getiCellfromVect(D3Q27ICell &iCell,double *Vect)
{

   for (int i=0*27;i<1*27;i++)         {   iCell.TSW[i-0*27]=Vect[i];   }
   for (int i=1*27;i<2*27;i++)		   {   iCell.TNW[i-1*27]=Vect[i];   }
   for (int i=2*27;i<3*27;i++)		   {   iCell.TNE[i-2*27]=Vect[i];   }
   for (int i=3*27;i<4*27;i++)         {   iCell.TSE[i-3*27]=Vect[i];   }
   for (int i=4*27;i<5*27;i++)		   {   iCell.BSW[i-4*27]=Vect[i];   }
   for (int i=5*27;i<6*27;i++)		   {   iCell.BNW[i-5*27]=Vect[i];   }
   for (int i=6*27;i<7*27;i++)		   {   iCell.BNE[i-6*27]=Vect[i];   }
   for (int i=7*27;i<8*27;i++)		   {   iCell.BSE[i-7*27]=Vect[i];   }
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::getAfterCompare(double dis[],double input [7][25],double AllCell[7][8*27],double &xp000,double &yp000,double &zp000,
                                                 double &Xw,double &Yw,double &Zw,double &cofWeightx,double &cofWeighty,double &cofWeightz,double &dx,double &omega,D3Q27ICell &iCell)
{

   double *Cell=(double*)malloc(8*27*sizeof(double));   	   
   bool goinside=true;

   for (int i=0;i<7;i++)
   {	
      double dist=dis[i];
      if(goinside )
      {  //17=cofWeightx,18=cofWeighty,19=cofWeightz,,20=dx,21=omega
         if      (dist==input[0][23]){  Xw=input[0][6];	Yw=input[0][7]; Zw=input[0][8]; xp000=input[0][0];	yp000=input[0][1]; zp000=input[0][2];  cofWeightx=input[0][17];	cofWeighty=input[0][18]; cofWeightz=input[0][19]; dx=input[0][20];	omega=input[0][21]; for (int i=0;i<8*27;i++){Cell[i]=AllCell[0][i];} } //get suitable block from suitable cell
         else if (dist==input[1][23]){  Xw=input[1][6];	Yw=input[1][7]; Zw=input[1][8]; xp000=input[1][0];	yp000=input[1][1]; zp000=input[1][2];  cofWeightx=input[1][17];	cofWeighty=input[1][18]; cofWeightz=input[1][19]; dx=input[1][20];	omega=input[1][21]; for (int i=0;i<8*27;i++){Cell[i]=AllCell[1][i];} }
         else if (dist==input[2][23]){  Xw=input[2][6];	Yw=input[2][7]; Zw=input[2][8]; xp000=input[2][0];	yp000=input[2][1]; zp000=input[2][2];  cofWeightx=input[2][17];	cofWeighty=input[2][18]; cofWeightz=input[2][19]; dx=input[2][20];	omega=input[2][21]; for (int i=0;i<8*27;i++){Cell[i]=AllCell[2][i];} }
         else if (dist==input[3][23]){  Xw=input[3][6];   Yw=input[3][7]; Zw=input[3][8]; xp000=input[3][0];  yp000=input[3][1]; zp000=input[3][2];  cofWeightx=input[3][17]; cofWeighty=input[3][18]; cofWeightz=input[3][19]; dx=input[3][20];  omega=input[3][21]; for (int i=0;i<8*27;i++){Cell[i]=AllCell[3][i];} }
         else if (dist==input[4][23]){  Xw=input[4][6];	Yw=input[4][7]; Zw=input[4][8]; xp000=input[4][0];	yp000=input[4][1]; zp000=input[4][2];  cofWeightx=input[4][17];	cofWeighty=input[4][18]; cofWeightz=input[4][19]; dx=input[4][20];	omega=input[4][21]; for (int i=0;i<8*27;i++){Cell[i]=AllCell[4][i];} }
         else if (dist==input[5][23]){  Xw=input[5][6];	Yw=input[5][7]; Zw=input[5][8]; xp000=input[5][0];	yp000=input[5][1]; zp000=input[5][2];  cofWeightx=input[5][17];	cofWeighty=input[5][18]; cofWeightz=input[5][19]; dx=input[5][20];	omega=input[5][21]; for (int i=0;i<8*27;i++){Cell[i]=AllCell[5][i];} }
         else if (dist==input[6][23]){  Xw=input[6][6];	Yw=input[6][7]; Zw=input[6][8]; xp000=input[6][0];	yp000=input[6][1]; zp000=input[6][2];  cofWeightx=input[6][17];	cofWeighty=input[6][18]; cofWeightz=input[6][19]; dx=input[6][20];	omega=input[6][21]; for (int i=0;i<8*27;i++){Cell[i]=AllCell[6][i];} }
         if (Xw<2.5&&Xw>-2.5&&Yw<2.50&&Yw>-2.50&&Zw<2.50&&Zw>-2.50){goinside=false;}      
      }
   }

   getiCellfromVect(iCell,Cell);
   free (Cell);
}

//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::SendAndReceiveData(double input[7][25],double AllCell[7][8*27],int allRank[7],int allTag[7])
{

   int sizeinput =25;
   double *Cell=(double*)malloc(8*27*sizeof(double));	


   for (int i=0;i<7;i++)
   {
      if (allRank[i]==this->root)
      {
         double vectforTrans[25];
         for (int j=0;j<sizeinput;j++)
         {
            vectforTrans[j]=input[i][j];
         }
         getAllInfoForCompare(vectforTrans,Cell);
         int index=(int)vectforTrans[24];
         for (int k=0;k<8*27;k++)
         {
            AllCell[index][k]=Cell[k];
         }
         for (int j=0;j<sizeinput;j++)
         {
            input[index][j]=vectforTrans[j];
         }
      }
      else
      {
         double vectforTrans[25];
         for (int j=0;j<sizeinput;j++)
         {
            vectforTrans[j]=input[i][j];
         }
         MPI_Send(vectforTrans,sizeinput, MPI_DOUBLE_PRECISION,allRank[i],allTag[i],MPI_COMM_WORLD);
      }
   }
   for (int i=0;i<7;i++)
   {
      if (allRank[i]!=this->root)
      {
         MPI_Status status; 
         double receiveVector[25];
         double *receiceCell=(double*)malloc(8*27*sizeof(double));
         MPI_Recv(receiveVector,25, MPI_DOUBLE_PRECISION,allRank[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
         MPI_Recv(receiceCell,8*27, MPI_DOUBLE_PRECISION,allRank[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);				  
         int index=(int)receiveVector[24];
         for (int j=0;j<sizeinput;j++)
         {
            input[index][j]=receiveVector[j];
         }
         for (int k=0;k<8*27;k++)
         {
            AllCell[index][k]=receiceCell[k];
         }
         free (receiceCell);
      }
   }
   free (Cell);
}
//////////////////////////////////////////////////////////////////////////
void PathLineCoProcessor::sendtoOtherIsOk(int isPointSuitable)
{
   if(comm->getNumberOfProcesses() > 1)  
   { 
      int size=comm->getNumberOfProcesses();
      for (int i=0;i<size;i++)
      {
         if (i!=rank)
         {
            MPI_Send(&isPointSuitable,1, MPI_INT,i,i,MPI_COMM_WORLD); 
         }	
      }
   }
}
