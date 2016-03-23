#include "D3Q27PathLinePostprocessorMcpart.h"
#include "LBMKernelETD3Q27.h"
#include "SimulationParameters.h"
#include "D3Q27ETBCProcessor.h"
#include <vector>
#include <string>
#include <boost/foreach.hpp>
#include "basics/writer/WbWriterVtkXmlASCII.h"
//#include "D3Q27ETFCoffVectorConnector.h"
using namespace std;


D3Q27PathLinePostprocessorMcpart::D3Q27PathLinePostprocessorMcpart(Grid3DPtr grid, const std::string& path, WbWriter* const writer,
                                                                   LBMUnitConverterPtr conv, UbSchedulerPtr s, CommunicatorPtr comm, 
                                                                   std::vector<UbTupleDouble3 > _Positions, LBMReal nue, D3Q27InterpolationProcessorPtr iProcessor)
                                                                   : Postprocessor(grid, s),
                                                                   path(path),
                                                                   comm(comm),
                                                                   writer(writer),
                                                                   conv(conv),
                                                                   nue(nue),
                                                                   iProcessor(iProcessor),
                                                                   istep(0),
                                                                   particleHasMass(true),
                                                                   rank(comm->getProcessID())
{
   iHelper = D3Q27InterpolationHelperPtr(new D3Q27InterpolationHelper(iProcessor));
   getParticlePerProccess(_Positions);
   getNeighborsRank();
   initializeForPrinting((int)_Positions.size());
}
//////////////////////////////////////////////////////////////////////////
D3Q27PathLinePostprocessorMcpart::~D3Q27PathLinePostprocessorMcpart()
{
}
//////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::update(double step)
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
void D3Q27PathLinePostprocessorMcpart::collectPostprocessData()
{

   int ran=rank;
   updateParticles();
   checkParticles();
   initialMovement();
   updateParticles();
   checkParticles();
   BOOST_FOREACH(ParticlesPtr particle, particles)
   { 
      double tau[6];//={tauxx,tauyy,tauzz,tauxy,tauxz,tauyz};
      finalMovement(particle,tau);	
      //printParticle(particle);
      gatherData(particle);
   }		

   if(scheduler->isDue((double)istep) )
   {
      string partName = writer->writeNodesWithNodeDataDouble(path+ UbSystem::toString(rank),nodes,datanames,data);
      size_t found=partName.find_last_of("//");
      string piece = partName.substr(found+1);

      vector<string> cellDataNames;

      vector<string> pieces = comm->gather(piece);
      if (comm->getProcessID() == comm->getRoot())
      {
         string pname = WbWriterVtkXmlASCII::getInstance()->writeParallelFile(path,pieces,datanames,cellDataNames);

         //vector<string> filenames;
         //filenames.push_back(pname);
         //if (istep == Postprocessor::scheduler->getMinBegin())
         //{
         //   WbWriterVtkXmlASCII::getInstance()->writeCollection(path+"_collection",filenames,0,false);
         //} 
         //else
         //{
         //   WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(path+"_collection",filenames,0,false);
         //}
         UBLOG(logINFO,"D3Q27MacroscopicQuantitiesPostprocessor step: " << istep);
      }
   }

   istep++;
}
////////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::initializeForPrinting(int _number)
{
   for (int i=0;i<_number;i++)
   {
      string fileName=path + UbSystem::toString(i)+ ".txt";
      files.push_back(boost::shared_ptr<std::ofstream>(new std::ofstream(fileName.c_str())));	
   }
}
//////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::getParticlePerProccess(std::vector<UbTupleDouble3 >particlesVec)
{
   for (int i=0;i<particlesVec.size();i++)
   {
      int minInitLevel = this->grid->getCoarsestInitializedLevel();
      int maxInitLevel = this->grid->getFinestInitializedLevel();
      LBMKernelETD3Q27Ptr kernel;
      DistributionArray3DPtr distributions;
      BCArray3D<D3Q27BoundaryCondition> bcArray;
      Block3DPtr block;

      for(int level = minInitLevel; level<=maxInitLevel; level++)
      {	      
         UbTupleInt3 blockIndexes = grid->getBlockIndexes(val<1>(particlesVec[i]), val<2>(particlesVec[i]), val<3>(particlesVec[i]),level);
         block = grid->getBlock(val<1>(blockIndexes), val<2>(blockIndexes), val<3>(blockIndexes), level);
         if(!block) continue; 
         if(block->isActive()) 
         {
            if (block->hasInterpolationFlagCF() || block->hasInterpolationFlagFC())
            {
               if(rank == block->getRank())
               { 					
                  if(!checkNodes(block,val<1>(particlesVec[i]), val<2>(particlesVec[i]), val<3>(particlesVec[i])))
                  {	
                     sendStatusOfPoint(false,val<1>(particlesVec[i]), val<2>(particlesVec[i]), val<3>(particlesVec[i]),level,i);
                     continue;
                  }
                  initializeParticle ( val<1>(particlesVec[i]), val<2>(particlesVec[i]), val<3>(particlesVec[i]),i ,level);
                  //getParticle(particle,grid->getDeltaX(block),particle.x,particle.y,particle.z,level,i);
                  sendStatusOfPoint(true,val<1>(particlesVec[i]), val<2>(particlesVec[i]), val<3>(particlesVec[i]),level,i);
                  break;
               }
               else							
               {								
                  bool isPointSuitable=true;
                  receiveStatusOfPoint(isPointSuitable,block->getRank(),val<1>(particlesVec[i]), val<2>(particlesVec[i]), val<3>(particlesVec[i]),level)	;	
                  if (isPointSuitable){break;}							
               }		
            }
            else
            {
               if(rank == block->getRank())
               { 
                  if(!checkNodes(block,val<1>(particlesVec[i]), val<2>(particlesVec[i]), val<3>(particlesVec[i])))
                  {		   
                     continue;
                  }  
                  initializeParticle ( val<1>(particlesVec[i]), val<2>(particlesVec[i]), val<3>(particlesVec[i]),i,level );
                  Particles particle;
                  //getParticle(particle,grid->getDeltaX(block),particle.x,particle.y,particle.z,level,i);
                  //particlesVec.push_back(particle);
                  //particles.erase(t.first);
                  break;	
               }

            }
         }		
      }
   } 			
}
////////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::initializeParticle(double x, double y,double z,int i,int level)
{
   ParticlesPtr particle = ParticlesPtr(new Particles() );

   particle->xold=x;	                    particle->yold=y;  	                particle->zold=z;
   particle->x   =x;	                    particle->y   =y;	                particle->z   =z;			 	  
   particle->vxold =0.1;		            particle->vyold =0.1;	            particle->vzold=0.1;	 
   particle->vxoldf=0.0;		            particle->vyoldf=0.0;	            particle->vzoldf=0.0;
   particle->rankOfParticle=rank;          particle->ID=i;                     particle->level=level;
   particles.push_back(particle);
}
//////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::getNeighborsRank()
{
   int gridRank = grid->getRank();
   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   int maxInitLevel = this->grid->getFinestInitializedLevel();

   for(int level = minInitLevel; level<=maxInitLevel;level++)
   {
      vector<Block3DPtr> blockVector;
      grid->getBlocks(level, gridRank, blockVector);
      BOOST_FOREACH(Block3DPtr block, blockVector)
      {
         int blockRank = block->getRank();
         if (gridRank == blockRank && block->isActive())
         {				
            //std::vector<Block3DPtr> neighbors; 
            int ix1 = block->getX1();
            int ix2 = block->getX2();
            int ix3 = block->getX3();
            int level = block->getLevel();
            //grid->getAllNeighbors(ix1, ix2, ix3, level, level, neighbors);
            int dirs=26;//number of directions
            for( int dir = 0; dir < dirs; dir++)
            {
               Block3DPtr neighBlock = grid->getNeighborBlock(dir, ix1, ix2, ix3, level);
               if(neighBlock)
               {
                  int neighBlockRank = neighBlock->getRank();			
                  if(blockRank != neighBlockRank && neighBlock->isActive())
                  {
                     if (neighbors.size()==0)
                     {
                        neighbors.insert((std::make_pair(neighBlockRank,0)));
                     } 
                     else
                     {
                        for ( std::map< int, int >::const_iterator iter = neighbors.begin();iter != neighbors.end(); ++iter )
                        {
                           if (iter->first!=neighBlockRank)
                           {
                              neighbors.insert((std::make_pair(neighBlockRank,0)));
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

//////////////////////////////////////////////////////////////////////////
//void D3Q27PathLinePostprocessorMcpart::getParticle(Particles &particle, double dx,double x, double y,double z,int level,int i)
//{
//	particle.particleNumber=i;
//	particle.xold=val<1>(xold[i]);	        particle.yold=val<2>(xold[i]);  	particle.zold=val<3>(xold[i]);
//	particle.x   =x;	                    particle.y   =y;	                particle.z   =z;			 
//	particle.level=level;
//	particle.dx=dx;		 
//	particle.vxold=val<1>(vxold[i]);		particle.vyold=val<2>(vxold[i]);	particle.vzold=val<3>(vxold[i]);	 
//	particle.vxoldf=val<1>(vxoldf[i]);		particle.vyoldf=val<2>(vxoldf[i]);	particle.vzoldf=val<3>(vxoldf[i]);
//
//}
//////////////////////////////////////////////////////////////////////////
bool D3Q27PathLinePostprocessorMcpart::checkNodes( Block3DPtr block,double _x1,double _x2,double _x3)
{
   bool result = true;
   LBMKernelETD3Q27Ptr kernel = boost::dynamic_pointer_cast<LBMKernelETD3Q27>(block->getKernel());
   BCArray3D<D3Q27BoundaryCondition> bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();
   D3Q27BoundaryConditionPtr bcPtr;

   double x1_ch = _x1;// + vx1old*conv->getfactorTimeShouldMultiplebyDx()*grid->getDeltaX(block);
   double x2_ch = _x2;// + vx2old*conv->getfactorTimeShouldMultiplebyDx()*grid->getDeltaX(block);
   double x3_ch = _x3;// + vx3old*conv->getfactorTimeShouldMultiplebyDx()*grid->getDeltaX(block);

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
//	//////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::initialMovement()
{
   BOOST_FOREACH(ParticlesPtr particle, particles)
   {
      double dx=grid->getDeltaX(particle->level);
      particle->xold=particle->x;
      particle->yold=particle->y;
      particle->zold=particle->z;
      particle->x   =particle->xold+particle->vxold*dx*conv->getFactorTimeLbToW(dx);
      particle->y   =particle->yold+particle->vyold*dx*conv->getFactorTimeLbToW(dx);
      particle->z   =particle->zold+particle->vzold*dx*conv->getFactorTimeLbToW(dx);	
   }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void D3Q27PathLinePostprocessorMcpart::updateDistributedValues(std::vector<Particles> *particlesVec,std::vector<UbTupleDouble3> &x)
//{
//	
//	int numberOFVariable=18;
//	std::vector<double> sendParticlesInfomation;	
//	std::vector<double> receiveParticlesInfomation;	
//	getVectorFromParticles(sendParticlesInfomation,particlesVec,numberOFVariable);
//	allGatherDoubles(sendParticlesInfomation, receiveParticlesInfomation);
//	fillOutPositions(receiveParticlesInfomation,x);
//	//getParticlesFromVector(receiveParticlesInfomation,particlesVec,numberOFVariable);
//	//std::vector<Particles> *particlesVec;
//}
//////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::sendBool(bool _bool,int distinationRank ,int i)
{
   int variable;
   if (_bool==false)
   {
      variable=0;
   }
   else variable=1;
   if(comm->getNumberOfProcesses() > 1)  
   { 		
      MPI_Send(&variable,1, MPI_INT,distinationRank,i,MPI_COMM_WORLD); 
      //std::cout<<"\n"<<variable<<"send from rank "<<rank<<" to rank "<< distinationRank;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::sendInt(int svalue,int distinationRank ,int tag)
{
   if(comm->getNumberOfProcesses() > 1)  
   { 		
      MPI_Send(&svalue,1, MPI_INT,distinationRank,tag,MPI_COMM_WORLD); 
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::receiveBool(bool &variable,int rankBlock)
{
   if(comm->getNumberOfProcesses() > 1)  
   { 	
      int _bool;
      MPI_Status status; 
      MPI_Recv(&_bool,1, MPI_INT,rankBlock,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      //std::cout<<"\n"<<_bool<<"received in rank "<<rank<<" from rank "<< rankBlock;
      if (_bool==false)
      {
         variable=false;
      }
      else variable=true;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::receiveInt(int &rvalue,int root)
{
   if(comm->getNumberOfProcesses() > 1)  
   { 	
      MPI_Status status; 
      MPI_Recv(&rvalue,1, MPI_INT,root,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
   }
}
////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::sendStatusOfPoint(bool status,double x1,double x2,double x3,double level,int i)
{
   for(int lev = (int)level-1; lev<=(int)level+1; lev++)
   {	      	
      UbTupleInt3 blockIndexesDistination = grid->getBlockIndexes(x1,x2,x3,lev);
      Block3DPtr blockDistination= grid->getBlock(val<1>(blockIndexesDistination), val<2>(blockIndexesDistination), val<3>(blockIndexesDistination), lev);
      if(!blockDistination) continue; 
      if(blockDistination->isActive()) 	
      {
         if(rank!=blockDistination->getRank())
         {
            sendBool(status,blockDistination->getRank(),i );  
         }

      }
   }
}
//////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::receiveStatusOfPoint(bool &status,int rankRoot,double x1,double x2,double x3,double level)
{
   for(int lev = (int)level-1; lev<=(int)level+1; lev++)
   {	      
      UbTupleInt3 blockIndexes = grid->getBlockIndexes(x1,x2,x3,lev);
      Block3DPtr block= grid->getBlock(val<1>(blockIndexes), val<2>(blockIndexes), val<3>(blockIndexes), lev);
      if (block)
      {
         if(rankRoot!=block->getRank())
         {
            LBMKernelETD3Q27Ptr kernel = boost::dynamic_pointer_cast<LBMKernelETD3Q27>(block->getKernel());
            if(kernel) 
            {
               if(checkNodes(block,x1,x2,x3))
               {	
                  receiveBool(status,rankRoot);	
               }
            }

         }

      }
   }
}
//////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::getVectorFromParticles(std::vector<double> &particlesInfomation,ParticlesPtr particle)
{
   particlesInfomation.push_back(particle->x     ); 
   particlesInfomation.push_back(particle->y     );	
   particlesInfomation.push_back(particle->z     );
   particlesInfomation.push_back(particle->xold  ); 
   particlesInfomation.push_back(particle->yold  );	
   particlesInfomation.push_back(particle->zold  );
   particlesInfomation.push_back(particle->vxold ); 
   particlesInfomation.push_back(particle->vyold );	
   particlesInfomation.push_back(particle->vzold );
   particlesInfomation.push_back(particle->vxoldf); 
   particlesInfomation.push_back(particle->vyoldf);	
   particlesInfomation.push_back(particle->vzoldf);
   particlesInfomation.push_back(particle->rankOfParticle); 
   particlesInfomation.push_back(particle->ID);
   particlesInfomation.push_back(particle->level);

}
////////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::getParticlesFromVector(std::vector<double> particlesInfomation,int numberOFVariable)
{
   int numParticles= (int)particlesInfomation.size()/numberOFVariable;
   for (int i=0;i<numParticles;i++)
   {
      ParticlesPtr particle = ParticlesPtr(new Particles() );

      particle->x                 =particlesInfomation[i*numberOFVariable+0]; 
      particle->y     			=particlesInfomation[i*numberOFVariable+1];
      particle->z     			=particlesInfomation[i*numberOFVariable+2];
      particle->xold   		    =particlesInfomation[i*numberOFVariable+3];
      particle->yold  			=particlesInfomation[i*numberOFVariable+4];
      particle->zold  			=particlesInfomation[i*numberOFVariable+5];
      particle->vxold  		    =particlesInfomation[i*numberOFVariable+6];
      particle->vyold 			=particlesInfomation[i*numberOFVariable+7];
      particle->vzold 			=particlesInfomation[i*numberOFVariable+8];
      particle->vxoldf 		    =particlesInfomation[i*numberOFVariable+9];
      particle->vyoldf			=particlesInfomation[i*numberOFVariable+10];
      particle->vzoldf			=particlesInfomation[i*numberOFVariable+11];
      particle->rankOfParticle    =(int)particlesInfomation[i*numberOFVariable+12];
      particle->ID	            =(int)particlesInfomation[i*numberOFVariable+13];
      particle->level	            =(int)particlesInfomation[i*numberOFVariable+14];
      particles.push_back(particle);
   }
}
////////////////////////////////////////////////////////////////////////////
//void D3Q27PathLinePostprocessorMcpart::updateinfo(std::vector<Particles> *particlesVec,std::vector<UbTupleDouble3> &x)
//{
//	std::vector<Particles>&particles=*particlesVec;
//	for (int i=0;i<particlesVec->size();i++)
//	{
//		Particles &particle=particles[i];
//		int index=particle.particleNumber;
//		val<1>(x[index])=particle.xold;	        
//		val<2>(x[index])=particle.yold; 
//		val<3>(x[index])=particle.zold;
//	}
//}
////////////////////////////////////////////////////////////////////////////
//void D3Q27PathLinePostprocessorMcpart::allGatherDoubles(std::vector<double>& svalues, std::vector<double>& rvalues)
//{
//	int scount;
//	vector<int> displs, rcounts;
//
//	scount = (int)(svalues.size());
//	rcounts.resize(comm->getNumberOfProcesses());
//	MPI_Allgather(&scount, 1, MPI_INT, &rcounts[0], 1, MPI_INT, MPI_COMM_WORLD);
//	displs.resize(comm->getNumberOfProcesses());
//
//	displs[0] = 0; 
//	for (int i=1; i<comm->getNumberOfProcesses(); ++i) { 
//		displs[i] = displs[i-1]+rcounts[i-1]; 
//	}
//
//	rvalues.resize(displs[comm->getNumberOfProcesses()-1]+rcounts[comm->getNumberOfProcesses()-1]); 
//
//	if(rvalues.size() == 0)
//	{
//		rvalues.resize(1);
//		rvalues[0] = -999;
//	}
//	if(scount == 0)
//	{
//		svalues.resize(1);
//		svalues[0] = -999;
//	}
//
//	MPI_Allgatherv(&svalues[0], scount, MPI_DOUBLE, &rvalues[0], &rcounts[0], &displs[0], MPI_DOUBLE, MPI_COMM_WORLD); 
//}
////////////////////////////////////////////////////////////////////////////
//void D3Q27PathLinePostprocessorMcpart::fillOutPositions(std::vector<double> particlesInfomation,std::vector<UbTupleDouble3> &x,double numberOFVariable)
//{
//	if (particlesInfomation.size()!=x.size()*numberOFVariable)
//	{
//		std::cout<<"number of particle is wrong";
//	}
//	for (int i=0;i<x.size();i++)
//	{
//		Particles particle;
//		particle.x              =particlesInfomation[i*numberOFVariable+0]; 
//		particle.y     			=particlesInfomation[i*numberOFVariable+1];
//		particle.z     			=particlesInfomation[i*numberOFVariable+2];
//		particle.xold   		=particlesInfomation[i*numberOFVariable+3];
//		particle.yold  			=particlesInfomation[i*numberOFVariable+4];
//		particle.zold  			=particlesInfomation[i*numberOFVariable+5];
//		particle.vx     		=particlesInfomation[i*numberOFVariable+6];
//		particle.vy    			=particlesInfomation[i*numberOFVariable+7];
//		particle.vz    			=particlesInfomation[i*numberOFVariable+8];
//		particle.vxold  		=particlesInfomation[i*numberOFVariable+9];
//		particle.vyold 			=particlesInfomation[i*numberOFVariable+10];
//		particle.vzold 			=particlesInfomation[i*numberOFVariable+11];
//		particle.vxoldf 		=particlesInfomation[i*numberOFVariable+12];
//		particle.vyoldf			=particlesInfomation[i*numberOFVariable+13];
//		particle.vzoldf			=particlesInfomation[i*numberOFVariable+14];
//		particle.dx     		=particlesInfomation[i*numberOFVariable+15];
//		particle.level 			=particlesInfomation[i*numberOFVariable+16];
//		particle.particleNumber	=particlesInfomation[i*numberOFVariable+17];
//
//		int index=particle.particleNumber;
//		val<1>(x[index])=particle.x;	        val<2>(x[index])=particle.y;        	val<3>(x[index])=particle.z;
//		val<1>(xold[index])=particle.xold;	    val<2>(xold[index])=particle.yold;  	val<3>(xold[index])=particle.zold;
//		val<1>(vxold[index])=particle.vxold;	val<2>(vxold[index])=particle.vyold;  	val<3>(vxold[index])=particle.vzold;
//		val<1>(vxoldf[index])=particle.vxoldf;	val<2>(vxoldf[index])=particle.vyoldf;  	val<3>(vxoldf[index])=particle.vzoldf;
//	}
//
//}
//
//
//

//////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::updateParticles()
{
   //find particles which want to go to another process
   BOOST_FOREACH(ParticlesPtr particle, particles)
   {
      int minInitLevel = this->grid->getCoarsestInitializedLevel();
      int maxInitLevel = this->grid->getFinestInitializedLevel();
      LBMKernelETD3Q27Ptr kernel;
      DistributionArray3DPtr distributions;
      BCArray3D<D3Q27BoundaryCondition> bcArray;
      Block3DPtr block;

      for(int level = minInitLevel; level<=maxInitLevel; level++)
      {	      
         UbTupleInt3 blockIndexes = grid->getBlockIndexes(particle->x,particle->y,particle->z,level);
         block = grid->getBlock(val<1>(blockIndexes), val<2>(blockIndexes), val<3>(blockIndexes), level);
         if(!block) continue; 
         if(block->isActive()) 
         {
            if(rank != block->getRank())
            { 					
               std::map<int,int>::iterator p;
               p=neighbors.find(block->getRank());
               p->second++;
               particle->rankOfParticle=block->getRank();
               particle->level=level;
               break;//break should be here
            }
            //break;//break should be here
         }			
      }
   }
   //send number of particle moving to another process 
   for ( std::map< int, int >::const_iterator iter = neighbors.begin();iter != neighbors.end(); ++iter )
   {					
      int svalue=iter->second;
      if(comm->getNumberOfProcesses() > 1)  
      { 		
         MPI_Send(&svalue,1, MPI_INT,iter->first,iter->first,MPI_COMM_WORLD); 
         //std::cout<<"\n"<<svalue<<" send from rank "<<rank<<" to rank "<< iter->first;
      }
   }	
   std::map<int,int> receiveNeighbors;
   //receive number of particle coming form another process
   for ( std::map< int, int >::const_iterator iter = neighbors.begin();iter != neighbors.end(); ++iter )
   {
      int rvalue;	
      if(comm->getNumberOfProcesses() > 1)  
      { 	
         MPI_Status status; 
         MPI_Recv(&rvalue,1, MPI_INT,iter->first,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
         //std::cout<<"\n"<<"rank "<< rank<<" receive "<<rvalue<<" from rank "<<iter->first;
         receiveNeighbors.insert((std::make_pair(iter->first,rvalue)));
      }
   }
   //int numOfvarClass=18;

   ///send data
   for ( std::map< int, int >::const_iterator iter = neighbors.begin();iter != neighbors.end(); ++iter )
   {					
      if (iter->second!=0)
      {
         std::vector<double> sendParticlesInfomation;
         int lengthSending=numOfvarClass*iter->second;

         BOOST_FOREACH(ParticlesPtr particle, particles)
         {
            if (particle->rankOfParticle==iter->first)
            {
               getVectorFromParticles(sendParticlesInfomation,particle);
            }
         }
         MPI_Send(&sendParticlesInfomation[0],lengthSending, MPI_DOUBLE,iter->first,iter->first,MPI_COMM_WORLD); 
         std::cout<<"\n send from rank "<<rank<<" to rank "<< iter->first;
      }
   }	
   //delete particle
   std::list<ParticlesPtr> newparticles;
   BOOST_FOREACH(ParticlesPtr particle, particles)
   {
      if (particle->rankOfParticle==rank)
      {
         newparticles.push_back(particle);
      }
      else
      {
         //	int idNumber =particle->ID;
         //	printParticle(idNumber);
         //delete data
      }
   }
   particles=newparticles;
   ///receive data
   for ( std::map< int, int >::const_iterator iter = receiveNeighbors.begin();iter != receiveNeighbors.end(); ++iter )
   {					
      if (iter->second!=0)
      {
         int lengthReceiving=numOfvarClass*iter->second;	
         std::vector<double> rvalue(lengthReceiving);	
         MPI_Status status; 
         MPI_Recv(&rvalue[0],lengthReceiving, MPI_DOUBLE,iter->first,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
         std::cout<<"\n"<<"rank "<< rank<<" receive  from rank "<<iter->first;
         getParticlesFromVector(rvalue,numOfvarClass);
      }
   }
   //reset neighbors
   std::map<int,int> Neighbors2;
   Neighbors2=neighbors;
   neighbors.clear();
   for ( std::map< int, int >::const_iterator iter = Neighbors2.begin();iter != Neighbors2.end(); ++iter )
   {
      neighbors.insert((std::make_pair(iter->first,0)));	
   }	



}
//////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::checkParticles()
{
   std::list<ParticlesPtr> newparticles;
   BOOST_FOREACH(ParticlesPtr particle, particles)
   { 
      int minInitLevel = this->grid->getCoarsestInitializedLevel();
      int maxInitLevel = this->grid->getFinestInitializedLevel();
      Block3DPtr block;

      for(int level = minInitLevel; level<=maxInitLevel; level++)
      {	      
         UbTupleInt3 blockIndexes = grid->getBlockIndexes(particle->x,particle->y,particle->z,level);
         block = grid->getBlock(val<1>(blockIndexes), val<2>(blockIndexes), val<3>(blockIndexes), level);
         if(!block) continue; 
         if(block->isActive()) 
         {
            LBMKernelETD3Q27Ptr kernel = boost::dynamic_pointer_cast<LBMKernelETD3Q27>(block->getKernel());
            if(!kernel) continue; 
            if(rank == block->getRank())
            { 					
               if(!checkNodes(block,particle->x,particle->y,particle->z))
               {	
                  std::cout<<"particle number "<<particle->ID <<"is gone";
                  continue;
               }
               newparticles.push_back(particle);
               break;
            }
         }
      }
   }
   particles=newparticles;
}
//////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::finalMovement(ParticlesPtr particle,double tau[])
{
   Block3DPtr block;
   UbTupleInt3 blockIndexes = grid->getBlockIndexes(particle->x, particle->y, particle->z,particle->level);
   block = grid->getBlock(val<1>(blockIndexes), val<2>(blockIndexes), val<3>(blockIndexes), particle->level);
   LBMKernelETD3Q27Ptr kernel = boost::dynamic_pointer_cast<LBMKernelETD3Q27>(block->getKernel());
   DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();
   BCArray3D<D3Q27BoundaryCondition> bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();

   UbTupleInt3 nodeIndexes = grid->getNodeIndexes(block, particle->x, particle->y, particle->z);
   int ix1 = val<1>(nodeIndexes);
   int ix2 = val<2>(nodeIndexes);
   int ix3 = val<3>(nodeIndexes);

   if(!iProcessor->iCellHasSolid(bcArray, ix1, ix2, ix3))
   {
      interpolMovement(block,distributions,particle,tau,ix1, ix2, ix3);
   }
   else
   {
      extrapolMovement( block, particle, tau, ix1,  ix2,  ix3);
   }

}
//////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::interpolMovement(Block3DPtr block,DistributionArray3DPtr& distributions,ParticlesPtr particle,double tau[],int ix1, int ix2, int ix3)
{
   D3Q27ICell iCell;
   //LBMReal f[D3Q27System::ENDF+1];
   LBMReal dx =grid->getDeltaX(particle->level);
   LBMReal vx1,vx2,vx3;
   iProcessor->readICell(distributions, iCell, ix1, ix2, ix3);
   double x1LB, x2LB, x3LB;
   UbTupleDouble3 orgNodeRW =  grid->getNodeCoordinates(block,  ix1, ix2, ix3);
   double offsetX1RW = abs(particle->x - val<1>(orgNodeRW));
   double offsetX2RW = abs(particle->y - val<2>(orgNodeRW));
   double offsetX3RW = abs(particle->z - val<3>(orgNodeRW));
   x1LB = offsetX1RW / dx;
   x2LB = offsetX2RW / dx;
   x3LB = offsetX3RW / dx;

   x1LB -= 0.5;
   x2LB -= 0.5;
   x3LB -= 0.5;
   LBMReal omega = LBMSystem::calcCollisionFactor(nue, block->getLevel());
   LBMReal tauxx,	   tauyy,	   tauzz,	   tauxy,	   tauxz,	   tauyz;
   iHelper->interpolate8to1WithVelocityWithShearStress(iCell, x1LB, x2LB, x3LB, omega,vx1,vx2,vx3,tauxx,tauyy,tauzz,tauxy,tauxz,tauyz);
   vx1 = vx1 *  conv->getFactorVelocityLbToW();
   vx2 = vx2 *  conv->getFactorVelocityLbToW();
   vx3 = vx3 *  conv->getFactorVelocityLbToW();
   tau[0]=tauxx/dx;
   tau[1]=tauyy/dx;
   tau[2]=tauzz/dx;
   tau[3]=tauxy/dx;
   tau[4]=tauxz/dx;
   tau[5]=tauyz/dx;  

   if (particleHasMass) {   CalcVelParticle(dx,vx1,vx2,vx3,particle->vxold,particle->vyold,particle->vzold,particle->vxoldf,particle->vyoldf,particle->vzoldf);}
   //if (particleHasMass) {   CalcVelParticle2(dx,vx1,vx2,vx3,particle->vxold,particle->vyold,particle->vzold,particle->vxoldf,particle->vyoldf,particle->vzoldf);}
   //heuns method
   particle->x= (particle->xold) + 1.0/2.0*(vx1 + particle->vxold)*conv->getFactorTimeLbToW(dx);
   particle->y= (particle->yold) + 1.0/2.0*(vx2 + particle->vyold)*conv->getFactorTimeLbToW(dx);
   particle->z= (particle->zold) + 1.0/2.0*(vx3 + particle->vzold)*conv->getFactorTimeLbToW(dx);

   particle->vxold = vx1;
   particle->vyold = vx2;
   particle->vzold = vx3;
}
//////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::extrapolMovement(Block3DPtr block,ParticlesPtr particle,double tau[],int ix1, int ix2, int ix3)
{

   LBMReal vx1,vx2,vx3;
   LBMReal dx =grid->getDeltaX(particle->level);
   int level =particle->level;
   LBMReal omega = LBMSystem::calcCollisionFactor(nue,level);
   D3Q27ICell iCell;
   //LBMReal f[D3Q27System::ENDF+1];
   double xp000,yp000,zp000;
   double Xw=999.0,Yw=999.0,Zw=999.0;
   double cofWeightx=999.0,cofWeighty=999.0,cofWeightz=999.0;
   double A,B,C,D,ii=0.0;
   double disx=9990,disy=9990,disz=9990,disxy=9990,disxz=9990,disyz=9990,disxyz=9990;

   UbTupleDouble3 orgNodeRW =  grid->getNodeCoordinates(block,  ix1, ix2, ix3);
   double ix1ph = ( val<1>(orgNodeRW));
   double ix2ph = ( val<2>(orgNodeRW));
   double ix3ph = ( val<3>(orgNodeRW)); 
   findPlane(ix1,ix2,ix3,grid,block,A,B,C,D,ii);
   double s = A*particle->x + B*particle->y + C*particle->z + D;//The sign of s = Ax + By + Cz + D determines which side the point (x,y,z) lies with respect to the plane. If s > 0 then the point lies on the same side as the normal (A,B,C). If s < 0 then it lies on the opposite side, if s = 0 then the point (x,y,z) lies on the plane.
   if (s>0){s=1;} else if (s<0){s=-1;}else {s=0;}
   double normalDis=((A*particle->x + B*particle->y + C*particle->z + D)/sqrt(A*A+B*B+C*C));///distance point to plane xp-Xw=distance
   double di=A/sqrt(A*A+B*B+C*C);      double dj=B/sqrt(A*A+B*B+C*C);       double dk=C/sqrt(A*A+B*B+C*C);
   double XXw=particle->x-di*normalDis;    double YYw=particle->y-dj*normalDis;    double ZZw=particle->z-dk*normalDis;
   findInitialCell(s,di,dj,dk,dx,ix1ph,ix2ph,ix3ph,xp000,yp000,zp000);//find initial cell with xp000,yp000,zp000

   double dis[7]={disx,disy,disz,disxy,disxz,disyz,disxyz};   
   getAllDis(dis,particle->x,particle->y,particle->z,ix1ph,ix2ph,ix3ph,xp000,yp000,zp000,level);
   rearangedDouble(dis);
   getAfterCompare(dis,dx,ix1ph,ix2ph,ix3ph,xp000,yp000,zp000,Xw,Yw,Zw,particle->x,particle->y,particle->z,A,B,C, D, normalDis,di,dj,dk,cofWeightx,cofWeighty,cofWeightz,iCell,level);

   double x1LB, x2LB, x3LB;
   double offsetX1RW000 = (particle->x - xp000);
   double offsetX2RW000 = (particle->y - yp000);
   double offsetX3RW000 = (particle->z - zp000);

   x1LB = offsetX1RW000 / dx;
   x2LB = offsetX2RW000 / dx;
   x3LB = offsetX3RW000 / dx;
   x1LB -= 0.5;
   x2LB -= 0.5;
   x3LB -= 0.5;
   // outICell(iCell);
   //iProcessor->interpolate8to1WithVelocity(iCell, f, x1LB, x2LB, x3LB, omega,vx1,vx2,vx3);
   LBMReal tauxx,	   tauyy,	   tauzz,	   tauxy,	   tauxz,	   tauyz;
   iHelper->interpolate8to1WithVelocityWithShearStress(iCell, x1LB, x2LB, x3LB, omega,vx1,vx2,vx3,tauxx,tauyy,tauzz,tauxy,tauxz,tauyz);	
   addCorrection(vx1,vx2,vx3,tauxx,tauyy,tauzz,tauxy,tauxz,tauyz,dx,iCell,Xw, Yw, Zw, omega,cofWeightx,cofWeighty,cofWeightz,ii, x1LB, x2LB, x3LB,di,dj,dk);
   tau[0]=tauxx;
   tau[1]=tauyy;
   tau[2]=tauzz;
   tau[3]=tauxy;
   tau[4]=tauxz;
   tau[5]=tauyz;
   // if (particleHasMass) {   CalcVelParticle(dx,vx1,vx2,vx3,vx1old,vx2old,vx3old,vx1oldf,vx2oldf,vx3oldf);}
   if (particleHasMass) {   CalcVelParticle2(dx,vx1,vx2,vx3,particle->vxold,particle->vyold,particle->vzold,particle->vxoldf,particle->vyoldf,particle->vzoldf);}
   //heuns method
   particle->x= (particle->xold) + 1.0/2.0*(vx1 + particle->vxold)*conv->getFactorTimeLbToW(dx);
   particle->y= (particle->yold) + 1.0/2.0*(vx2 + particle->vyold)*conv->getFactorTimeLbToW(dx);
   particle->z= (particle->zold) + 1.0/2.0*(vx3 + particle->vzold)*conv->getFactorTimeLbToW(dx);
   collWall(A,B,C,D,particle->x,particle->y,particle->z,particle->xold,particle->yold,particle->zold,dx,vx1,vx2,vx3,ii);
   particle->vxold = vx1;
   particle->vyold = vx2;
   particle->vzold = vx3;

}
//////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::CalcVelParticle(double dx,double &vx1,double &vx2,double &vx3,double vx1old,double vx2old,
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
void D3Q27PathLinePostprocessorMcpart::CalcVelParticle2(double dx,double &vx1,double &vx2,double &vx3,double vx1old,double vx2old,
                                                        double vx3old,double &vx1oldf,double &vx2oldf,double &vx3oldf)
{
   LBMReal muRE = 1.002*1e-3;//m2/s
   double mass=2.3*1e-17;// kg
   double Diameter=600*1e-9;//m
   double deltatime=conv->getFactorTimeLbToW(dx)*0.001;
   double Coff=3*PI*Diameter*muRE*deltatime/mass;
   double exCoff=exp(-Coff);
   //velocity of particle/// 
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
void D3Q27PathLinePostprocessorMcpart::findPlane(int ix1,int ix2,int ix3,Grid3DPtr grid,Block3DPtr block,double &A,double &B,double &C,double &D,double &ii)
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
void D3Q27PathLinePostprocessorMcpart::findInitialCell(double s,double di,double dj,double dk,double dx,double ix1ph,double ix2ph,double ix3ph,double &xp000,double &yp000,double &zp000)
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
void D3Q27PathLinePostprocessorMcpart::finddisPointToCentercell(double x1,double x2,double x3,double xp000, double yp000, double zp000,double &dis,int level)
{	
   UbTupleInt3 blockIndexes = grid->getBlockIndexes(xp000, yp000,  zp000,level);
   Block3DPtr block = grid->getBlock(val<1>(blockIndexes), val<2>(blockIndexes), val<3>(blockIndexes), level);
   if(block&& block->isActive())
   {
      LBMKernelETD3Q27Ptr kernel = boost::dynamic_pointer_cast<LBMKernelETD3Q27>(block->getKernel());
      if (kernel)
      {
         {
            DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();
            BCArray3D<D3Q27BoundaryCondition> bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();
            if(block && block->isActive())
            {	
               double dx = grid->getDeltaX(block);
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

   }
}
//////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::collWall(double A,double B,double C,double D,double &x1,double &x2,double &x3,double x1old,double x2old,double x3old,double dx,double &vx1,double &vx2,double &vx3,double ii)
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
void D3Q27PathLinePostprocessorMcpart::addCorrection(double &vx1,double &vx2,double &vx3,double &tauxx,double &tauyy,double &tauzz,double &tauxy,double &tauxz,double &tauyz,
                                                     double dx,D3Q27ICell iCell,double Xw, double Yw, double Zw,double  omega,double cofWeightx,double cofWeighty,double cofWeightz,double ii,double x1LB,double x2LB,double x3LB,double di,double dj,double dk)
{
   //LBMReal f[D3Q27System::ENDF+1];
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
void D3Q27PathLinePostprocessorMcpart::getAllDis(double dis[],double x,double y,double z,double ix1ph,double ix2ph,double ix3ph,double xp000,double yp000,double zp000,int level)
{
   finddisPointToCentercell(x,y,z,xp000, ix2ph,  ix3ph,dis[0],level);//disx  
   finddisPointToCentercell(x,y,z,ix1ph, yp000,  ix3ph,dis[1],level);//disy  
   finddisPointToCentercell(x,y,z,ix1ph, ix2ph,  zp000,dis[2],level);//disz  
   finddisPointToCentercell(x,y,z,xp000, yp000,  ix3ph,dis[3],level);//disxy 
   finddisPointToCentercell(x,y,z,xp000, ix2ph,  zp000,dis[4],level);//disxz 
   finddisPointToCentercell(x,y,z,ix1ph, yp000,  zp000,dis[5],level);//disyz 
   finddisPointToCentercell(x,y,z,xp000, yp000,  zp000,dis[6],level);//disxyz  	   
}
//////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::getAfterCompare(double dis[],double dx,double ix1ph,double ix2ph,double ix3ph,double &xp000,double &yp000,double &zp000,double &Xw,
                                                       double &Yw,double &Zw,double x1,double x2,double x3,double A,double B,double C,double D,double normalDis, double di,double dj,double dk, double &cofWeightx
                                                       ,double &cofWeighty,double &cofWeightz,D3Q27ICell &iCell,int level)
{
   for (int i=0;i<7;i++)
   {	
      double dist=dis[i];
      if      (dist==dirx  ){ Xw=(B*x2+C*x3+D)/(-A); Yw=x2;                 Zw=x3;	         	cofWeightx=1; cofWeighty=0; cofWeightz=0;  xp000=xp000;	yp000=ix2ph; zp000=ix3ph; getIcell(iCell,xp000, ix2ph, ix3ph,level); } //get suitable block from suitable cell
      else if (dist==diry  ){ Xw=x1;                 Yw=(A*x1+C*x3+D)/(-B); Zw=x3;		        cofWeightx=0; cofWeighty=1; cofWeightz=0;  xp000=ix1ph;	yp000=yp000; zp000=ix3ph; getIcell(iCell,ix1ph, yp000, ix3ph,level); } 
      else if (dist==dirz  ){ Xw=x1;                 Yw=x3;	              Zw=(B*x2+A*x1+D)/(-C);cofWeightx=0; cofWeighty=0; cofWeightz=1;  xp000=ix1ph;	yp000=ix2ph; zp000=zp000; getIcell(iCell,ix1ph, ix2ph, zp000,level); } 
      else if (dist==dirxy ){ Xw=x1-di*normalDis;    Yw=x2-dj*normalDis;    Zw=x3;        		cofWeightx=1; cofWeighty=1; cofWeightz=0;  xp000=xp000; yp000=yp000; zp000=ix3ph; getIcell(iCell,xp000, yp000, ix3ph,level); } 
      else if (dist==dirxz ){ Xw=x1-di*normalDis;    Yw=x2;                 Zw=x3-dk*normalDis;	cofWeightx=1; cofWeighty=0; cofWeightz=1;  xp000=xp000;	yp000=ix2ph; zp000=zp000; getIcell(iCell,xp000, ix2ph, zp000,level); } 
      else if (dist==diryz ){ Xw=x1;                 Yw=x2-dj*normalDis;    Zw=x3-dk*normalDis;	cofWeightx=0; cofWeighty=1; cofWeightz=1;  xp000=ix1ph;	yp000=yp000; zp000=zp000; getIcell(iCell,ix1ph, yp000, zp000,level); } 
      else if (dist==dirxyz){ Xw=x1-di*normalDis;    Yw=x2-dj*normalDis;    Zw=x3-dk*normalDis;	cofWeightx=1; cofWeighty=1; cofWeightz=1;  xp000=xp000;	yp000=yp000; zp000=zp000; getIcell(iCell,xp000, yp000, zp000,level); } 

      Xw=(Xw-xp000)/dx-0.5;       
      Yw=(Yw-yp000)/dx-0.5;       
      Zw=(Zw-zp000)/dx-0.5;
      if (Xw<2.5&&Xw>-2.5&&Yw<2.50&&Yw>-2.50&&Zw<2.50&&Zw>-2.50){break;}      
   }
}

//////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::rearangedDouble( double dis[7])
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
void D3Q27PathLinePostprocessorMcpart::getIcell(D3Q27ICell &iCell,double xp000, double yp000, double zp000,int level)
{ 
   UbTupleInt3 blockIndexes = grid->getBlockIndexes(xp000, yp000, zp000,level);
   Block3DPtr block= grid->getBlock(val<1>(blockIndexes), val<2>(blockIndexes), val<3>(blockIndexes), level);
   LBMKernelETD3Q27Ptr kernel = boost::dynamic_pointer_cast<LBMKernelETD3Q27>(block->getKernel());
   DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();
   iProcessor->readICell(distributions, iCell, (int)xp000, (int)yp000, (int)zp000);
}
////////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::printParticle(ParticlesPtr particle)

{
   /*BOOST_FOREACH(ParticlesPtr particle, particles)
   {*/
   int index=particle->ID;
   (*files[index]).precision (std::numeric_limits<double>::digits10 + 1);
   (*files[index])<<index<<" "<<istep<<" "<<particle->x<<"  "<< particle->y<<"  "<<particle->z<<" "<<particle->vxold<<"  "<<particle->vyold<<"  "<<particle->vzold<< std::endl;
   //}
}

///////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::printParticle(int index)

{

}
//////////////////////////////////////////////////////////////////////////
void D3Q27PathLinePostprocessorMcpart::gatherData(ParticlesPtr particle)
{
   //int index=particle->ID;
   //std::vector<std::vector<double> > tempData;
   //tempData.resize(7);
   //if (!data.count(index ))	
   //{	
   //	data.insert((std::make_pair(index,tempData)));	
   //}
   //tempData=data.at(index);
   //int t=0;
   //tempData[t++].push_back(istep);
   //tempData[t++].push_back(particle->x);
   //tempData[t++].push_back(particle->y);
   //tempData[t++].push_back(particle->z);
   //tempData[t++].push_back(particle->vxold);
   //tempData[t++].push_back(particle->vyold);
   //tempData[t++].push_back(particle->vzold);

   //  data.at(index)=tempData;

   datanames.resize(0);
   datanames.push_back("ID");
   datanames.push_back("TimeStep");
   datanames.push_back("Vx");
   datanames.push_back("Vy");
   datanames.push_back("Vz");
   //datanames.push_back("tauxx");
   //datanames.push_back("tauyy");
   //datanames.push_back("tauzz");
   //datanames.push_back("tauxy");
   //datanames.push_back("tauxz");
   //datanames.push_back("tauyz");
   //datanames.push_back("xoff");
   //datanames.push_back("yoff");
   //datanames.push_back("zoff");

   data.resize(datanames.size());

   int index = 0;
   data[index++].push_back(particle->ID);
   data[index++].push_back(istep);
   data[index++].push_back(particle->vxold);
   data[index++].push_back(particle->vyold);
   data[index++].push_back(particle->vzold);
   //data[index++].push_back(tauxx);
   //data[index++].push_back(tauyy);
   //data[index++].push_back(tauzz);
   //data[index++].push_back(tauxy);
   //data[index++].push_back(tauxz);
   //data[index++].push_back(tauyz);
   //data[index++].push_back(xoff);
   //data[index++].push_back(yoff);
   //data[index++].push_back(zoff);

   nodes.push_back( makeUbTuple(  double(particle->x), double(particle->y), double(particle->z)) );

}
