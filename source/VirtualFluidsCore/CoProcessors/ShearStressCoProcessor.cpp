#include "ShearStressCoProcessor.h"
#include "BCProcessor.h"
#include "WbWriterVtkXmlASCII.h"


ShearStressCoProcessor::ShearStressCoProcessor(Grid3DPtr grid, const std::string& path, 
                                                             WbWriter* const writer,
                                                             UbSchedulerPtr s,UbSchedulerPtr rs)
                                                             : CoProcessor(grid, s),
                                                             Resetscheduler(rs),
                                                             path(path),
                                                             writer(writer)
{
   CommunicatorPtr comm = Communicator::getInstance();
   normals.push_back(0);
   normals.push_back(0);
   normals.push_back(1);
   gridRank  = grid->getRank();
   minInitLevel = this->grid->getCoarsestInitializedLevel();
   maxInitLevel = this->grid->getFinestInitializedLevel();

   blockVector.resize(maxInitLevel+1);
   for(int level = minInitLevel; level<=maxInitLevel;level++)
   {
      grid->getBlocks(level, gridRank, true, blockVector[level]);
      BOOST_FOREACH(Block3DPtr block, blockVector[level])
      {
         UbTupleInt3 nx = grid->getBlockNX();
         ShearStressValuesArray3DPtr shearStressValues = ShearStressValuesArray3DPtr(new ShearStressValuesArray3D(14, val<1>(nx)+1, val<2>(nx)+1, val<3>(nx)+1, 0.0));
         block->getKernel()->getDataSet()->setShearStressValues(shearStressValues);
      }
   }
}
//////////////////////////////////////////////////////////////////////////
ShearStressCoProcessor::~ShearStressCoProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
void ShearStressCoProcessor::process( double step )
{
   if (step==0)
   {
      initDistance();
   }
   calculateShearStress(step);
   if(scheduler->isDue(step) )
      collectData(step);
   UBLOG(logDEBUG3, "D3Q27ShearStressCoProcessor::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void ShearStressCoProcessor::collectData(double step)
{
   using namespace std;

   int istep = int(step);
   addData();

   //string partName = writer->writeNodesWithNodeData(path+ UbSystem::toString(gridRank)+ "_" + UbSystem::toString(istep),nodes,datanames,data);
   //size_t found=partName.find_last_of("//");
   //string piece = partName.substr(found+1);

   //vector<string> cellDataNames;

   //CommunicatorPtr comm = Communicator::getInstance();
   //vector<string> pieces = comm->gatherStrings(piece);
   //if (comm->getProcessID() == comm->getRoot())
   //{
   //   string pname = WbWriterVtkXmlASCII::getInstance()->writeParallelFile(path+"_"+UbSystem::toString(istep),pieces,datanames,cellDataNames);

   //   vector<string> filenames;
   //   filenames.push_back(pname);
   //   if (step == CoProcessor::scheduler->getMinBegin())
   //   {
   //      WbWriterVtkXmlASCII::getInstance()->writeCollection(path+"__Shear_collection",filenames,istep,false);
   //   } 
   //   else
   //   {
   //      WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(path+"__Shear_collection",filenames,istep,false);
   //   }
   //   UBLOG(logINFO,"D3Q27ShearStressCoProcessor step: " << istep);
   //}

   string pfilePath, partPath, subfolder, cfilePath;
   subfolder = "shs"+UbSystem::toString(istep);
   pfilePath =  path+"/shs/"+subfolder;
   cfilePath =  path+"/shs/shs_collection";
   partPath = pfilePath+"/shs"+UbSystem::toString(gridRank)+ "_" + UbSystem::toString(istep);

   string partName = writer->writeNodesWithNodeData(partPath,nodes,datanames,data);
   size_t found=partName.find_last_of("/");
   string piece = partName.substr(found+1);
   piece = subfolder + "/" + piece;

   vector<string> cellDataNames;
   CommunicatorPtr comm = Communicator::getInstance();
   vector<string> pieces = comm->gather(piece);
   if (comm->getProcessID() == comm->getRoot())
   {
      string pname = WbWriterVtkXmlASCII::getInstance()->writeParallelFile(pfilePath,pieces,datanames,cellDataNames);
      found=pname.find_last_of("/");
      piece = pname.substr(found+1);

      vector<string> filenames;
      filenames.push_back(piece);
      if (step == CoProcessor::scheduler->getMinBegin())
      {
         WbWriterVtkXmlASCII::getInstance()->writeCollection(cfilePath,filenames,istep,false);
      } 
      else
      {
         WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(cfilePath,filenames,istep,false);
      }
      UBLOG(logINFO,"D3Q27ShearStressCoProcessor step: " << istep);
   }

   clearData();
}
//////////////////////////////////////////////////////////////////////////
void ShearStressCoProcessor::clearData()
{
   nodes.clear();
   datanames.clear();
   data.clear();
}
//////////////////////////////////////////////////////////////////////////
void ShearStressCoProcessor::calculateShearStress(double timeStep)
{
   using namespace D3Q27System;

   LBMReal f[27];
   LBMReal vx, vy, vz, sxx, syy, szz, sxy, syz, sxz;

   BOOST_FOREACH(D3Q27InteractorPtr interactor, interactors)
   {
      typedef std::map<Block3DPtr, std::set< std::vector<int> > > TransNodeIndicesMap;
      BOOST_FOREACH(TransNodeIndicesMap::value_type t, interactor->getTransNodeIndicesMap())
      {
         Block3DPtr block = t.first;
         std::set< std::vector<int> >& transNodeIndicesSet = t.second;

         LBMKernelPtr kernel = block->getKernel();
         BCArray3D& bcArray = kernel->getBCProcessor()->getBCArray();          
         DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions(); 
         ShearStressValuesArray3DPtr ssv = kernel->getDataSet()->getShearStressValues();

         int ghostLayer = kernel->getGhostLayerWidth();
         LBMReal collFactor = kernel->getCollisionFactor();

         int minX1 = ghostLayer;
         int maxX1 = (int)bcArray.getNX1() - 1 - ghostLayer;
         int minX2 = ghostLayer;
         int maxX2 = (int)bcArray.getNX2() - 1 - ghostLayer;
         int minX3 = ghostLayer;
         int maxX3 = (int)bcArray.getNX3() - 1 - ghostLayer;

         BOOST_FOREACH(std::vector<int> node, transNodeIndicesSet)
         {
            int ix1 = node[0];
            int ix2 = node[1];
            int ix3 = node[2];

            //without ghost nodes
            if (ix1 < minX1 || ix1 > maxX1 || ix2 < minX2 || ix2 > maxX2 ||ix3 < minX3 || ix3 > maxX3 ) continue;

            if(bcArray.isFluid(ix1,ix2,ix3)) 
            {
               double q=(*ssv)(normalq,ix1,ix2,ix3) ;
               double numPoint=(*ssv)(numberOfPoint,ix1,ix2,ix3) ;
               if (q==0||numPoint!=3)continue;
               // if (q==0)continue;
               //////////////////////////////////////////////////////////////////////////
               //read distribution
               ////////////////////////////////////////////////////////////////////////////
               distributions->getDistribution(f, ix1, ix2, ix3);
               //////////////////////////////////////////////////////////////////////////
               //compute velocity
               //////////////////////////////////////////////////////////////////////////
               vx = ((((f[TNE]-f[BSW]) + (f[TSE]-f[BNW])) + ((f[BSE]-f[TNW]) + (f[BNE]-f[TSW]))) +
                  (((f[BE]-f[TW]) + (f[TE]-f[BW])) + ((f[SE]-f[NW]) + (f[NE]-f[SW]))) +
                  (f[E]-f[W])); 

               vy = ((((f[TNE]-f[BSW]) + (f[BNW]-f[TSE])) + ((f[TNW]-f[BSE]) + (f[BNE]-f[TSW]))) +
                  (((f[BN]-f[TS]) + (f[TN]-f[BS])) + ((f[NW]-f[SE]) + (f[NE]-f[SW]))) +
                  (f[N]-f[S])); 

               vz = ((((f[TNE]-f[BSW]) + (f[TSE]-f[BNW])) + ((f[TNW]-f[BSE]) + (f[TSW]-f[BNE]))) +
                  (((f[TS]-f[BN]) + (f[TN]-f[BS])) + ((f[TW]-f[BE]) + (f[TE]-f[BW]))) +
                  (f[T]-f[B]));


               sxy = 3.0 * collFactor/(collFactor - 1.0)* ( ((f[TNE] + f[BSW])-(f[TSE]+f[BNW]))+(-(f[BSE]+f[TNW])+ (f[TSW]+f[BNE]))
                  +(((f[NE] + f[SW]) - (f[SE] + f[NW])))  -vx*vy);

               sxz = 3.0 * collFactor/(collFactor-1.0)*(((f[TNE] + f[BSW])+(f[TSE]+f[BNW]))+(-(f[BSE]+f[TNW])- (f[TSW]+f[BNE]))
                  +((f[TE] + f[BW])-(f[BE]+ f[TW])) -vx*vz);

               syz = 3.0 * collFactor/(collFactor-1.0)*(((f[TNE] + f[BSW])-(f[TSE]+f[BNW]))+((f[BSE]+f[TNW])- (f[TSW]+f[BNE]))
                  +(-(f[BN] + f[TS]) + (f[TN] + f[BS])) -vy*vz);                 

               LBMReal dxxMyy =3.0/2.0 * collFactor/(collFactor-1.0)* (((f[TE] + f[BW])+(f[BE]+ f[TW]))
                  -((f[BN] + f[TS]) + (f[TN] + f[BS]))+((f[E] + f[W])-(f[N] + f[S]))-vx*vx+vy*vy);

               LBMReal dxxMzz =3.0/2.0 * collFactor/(collFactor-1.0)*((((f[NE] + f[SW]) + (f[SE] + f[NW]))
                  -((f[BN] + f[TS]) + (f[TN] + f[BS])))+((f[E] + f[W])-(f[T] + f[B])) -vx*vx +vz*vz);

               // LBMReal dyyMzz =3.0/2.0 *collFactor/(collFactor-1.0)*((((f[NE] + f[SW]) + (f[SE] + f[NW]))-((f[TE] + f[BW])+(f[BE]+ f[TW])))
               //    +((f[N] + f[S])-(f[T] + f[B])) -vy*vy +vz*vz);

               sxx=(dxxMyy+dxxMzz)/3.0; // weil dxxPyyPzz=0

               syy=(dxxMzz-2*dxxMyy)/3.0;

               szz=(dxxMyy-2*dxxMzz)/3.0;

               //////////////////////////////////////////////////////////////////////////
               //compute average values
               //////////////////////////////////////////////////////////////////////////
               (*ssv)(AvVx,ix1,ix2,ix3) = ((*ssv)(AvVx,ix1,ix2,ix3)*timeStep + vx)/(timeStep+1.0);
               (*ssv)(AvVy,ix1,ix2,ix3) = ((*ssv)(AvVy,ix1,ix2,ix3)*timeStep + vy)/(timeStep+1.0);
               (*ssv)(AvVz,ix1,ix2,ix3) = ((*ssv)(AvVz,ix1,ix2,ix3)*timeStep + vz)/(timeStep+1.0);

               (*ssv)(AvSxx,ix1,ix2,ix3) = ((*ssv)(AvSxx,ix1,ix2,ix3)*timeStep + sxx)/(timeStep+1.0);
               (*ssv)(AvSyy,ix1,ix2,ix3) = ((*ssv)(AvSyy,ix1,ix2,ix3)*timeStep + syy)/(timeStep+1.0);
               (*ssv)(AvSzz,ix1,ix2,ix3) = ((*ssv)(AvSzz,ix1,ix2,ix3)*timeStep + szz)/(timeStep+1.0);
               (*ssv)(AvSxy,ix1,ix2,ix3) = ((*ssv)(AvSxy,ix1,ix2,ix3)*timeStep + sxy)/(timeStep+1.0);
               (*ssv)(AvSyz,ix1,ix2,ix3) = ((*ssv)(AvSyz,ix1,ix2,ix3)*timeStep + syz)/(timeStep+1.0);
               (*ssv)(AvSxz,ix1,ix2,ix3) = ((*ssv)(AvSxz,ix1,ix2,ix3)*timeStep + sxz)/(timeStep+1.0);
            }
         }

      }
   }
}
//////////////////////////////////////////////////////////////////////////
void ShearStressCoProcessor::addData()
{
   //Diese Daten werden geschrieben:
   datanames.resize(0);
   datanames.push_back("y^plus");
   datanames.push_back("u_tau");
   //datanames.push_back("yPlusFD");

   data.resize(datanames.size());

   BOOST_FOREACH(D3Q27InteractorPtr interactor, interactors)
   {
      typedef std::map<Block3DPtr, std::set< std::vector<int> > > TransNodeIndicesMap;
      BOOST_FOREACH(TransNodeIndicesMap::value_type t, interactor->getTransNodeIndicesMap())
      {
         Block3DPtr block = t.first;
         std::set< std::vector<int> >& transNodeIndicesSet = t.second;

         UbTupleDouble3 org          = grid->getBlockWorldCoordinates(block);
         UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
         UbTupleDouble3 nodeOffset   = grid->getNodeOffset(block);
         double         dx           = grid->getDeltaX(block);

         LBMKernelPtr kernel = block->getKernel();
         BCArray3D& bcArray = kernel->getBCProcessor()->getBCArray();          
         DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions(); 
         ShearStressValuesArray3DPtr ssv = kernel->getDataSet()->getShearStressValues();

         int ghostLayer = kernel->getGhostLayerWidth();
         LBMReal collFactor = kernel->getCollisionFactor();

         int minX1 = ghostLayer;
         int maxX1 = (int)bcArray.getNX1() - 1 - ghostLayer;
         int minX2 = ghostLayer;
         int maxX2 = (int)bcArray.getNX2() - 1 - ghostLayer;
         int minX3 = ghostLayer;
         int maxX3 = (int)bcArray.getNX3() - 1 - ghostLayer;

         int level=block->getLevel();
         if(level==1)
         {
            int le=0;
         }
         BOOST_FOREACH(std::vector<int> node, transNodeIndicesSet)
         {
            int ix1 = node[0];
            int ix2 = node[1];
            int ix3 = node[2];

            //without ghost nodes
            if (ix1 < minX1 || ix1 > maxX1 || ix2 < minX2 || ix2 > maxX2 ||ix3 < minX3 || ix3 > maxX3 ) continue;

            if(bcArray.isFluid(ix1,ix2,ix3)) 
            {
               double q=(*ssv)(normalq,ix1,ix2,ix3) ;
               double numPoint=(*ssv)(numberOfPoint,ix1,ix2,ix3) ;
               if (q==0||numPoint!=3)continue;
               // if (q==0)continue;

               int index = 0;
               nodes.push_back( makeUbTuple(float(val<1>(org) - val<1>(nodeOffset) + ix1*dx),
                  float(val<2>(org) - val<2>(nodeOffset) + ix2*dx),
                  float(val<3>(org) - val<3>(nodeOffset) + ix3*dx)) );

               //////get normal and distance//////
               double A,B,C;
               A=	(*ssv)(normalX1,ix1,ix2,ix3)	;
               B=	(*ssv)(normalX2,ix1,ix2,ix3)	;
               C=	(*ssv)(normalX3,ix1,ix2,ix3)	;

               ///////////
               //compute y plus
               //double vtxSonja, vtySonja, vtzSonja; //tangent velocity
               // double temp = (*av)(ix1,ix2,ix3,AvVx)*A+(*av)(ix1,ix2,ix3,AvVy)*B+(*av)(ix1,ix2,ix3,AvVz)*C;
               // vtxSonja = (*av)(ix1,ix2,ix3,AvVx)-normals[0]*temp;
               // vtySonja = (*av)(ix1,ix2,ix3,AvVy)-normals[1]*temp;
               // vtzSonja = (*av)(ix1,ix2,ix3,AvVz)-normals[2]*temp;

               double vtx=  (B*B*(*ssv)(AvVx,ix1,ix2,ix3) +  C*C*(*ssv)(AvVx,ix1,ix2,ix3) - A*B*(*ssv)(AvVy,ix1,ix2,ix3) - A*C*(*ssv)(AvVy,ix1,ix2,ix3))/(A*A+B*B+C*C);
               double vty=(-(A*B*(*ssv)(AvVx,ix1,ix2,ix3)) + A*A*(*ssv)(AvVy,ix1,ix2,ix3) + C*C*(*ssv)(AvVy,ix1,ix2,ix3) - B*C*(*ssv)(AvVz,ix1,ix2,ix3))/(A*A+B*B+C*C);
               double vtz=(-(A*C*(*ssv)(AvVx,ix1,ix2,ix3)) - B*C*(*ssv)(AvVy,ix1,ix2,ix3) + A*A*(*ssv)(AvVz,ix1,ix2,ix3) + B*B*(*ssv)(AvVz,ix1,ix2,ix3))/(A*A+B*B+C*C);

               double normVt = sqrt(vtx*vtx+vty*vty+vtz*vtz)+1e-100;
               double nvtx = vtx / normVt;
               double nvty = vty / normVt;
               double nvtz = vtz / normVt;

               double sx=0.5*((*ssv)(AvSxx,ix1,ix2,ix3)*nvtx+(*ssv)(AvSxy,ix1,ix2,ix3)*nvty+(*ssv)(AvSxz,ix1,ix2,ix3)*nvtz);
               double sy=0.5*((*ssv)(AvSxy,ix1,ix2,ix3)*nvtx+(*ssv)(AvSyy,ix1,ix2,ix3)*nvty+(*ssv)(AvSyz,ix1,ix2,ix3)*nvtz);
               double sz=0.5*((*ssv)(AvSxz,ix1,ix2,ix3)*nvtx+(*ssv)(AvSyz,ix1,ix2,ix3)*nvty+(*ssv)(AvSzz,ix1,ix2,ix3)*nvtz);
               double sabs=sqrt(sx*sx+sy*sy+sz*sz);

               double viscosity = (1.0/3.0)*(1.0/collFactor-0.5);
               double rho = 1.0;
               double utau=sqrt(viscosity/rho*sabs);

               // double q=(*av)(ix1,ix2,ix3,normalq) ;
               double yPlus = (utau*q)/viscosity;

               data[index++].push_back(yPlus);
               data[index++].push_back(utau);
            }
         }
      }
   }

}
//////////////////////////////////////////////////////////////////////////
void ShearStressCoProcessor::reset(double step)
{
   if(Resetscheduler->isDue(step) )
      resetData(step);

   UBLOG(logDEBUG3, "resetCoProcessor::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void ShearStressCoProcessor::resetData(double step)
{
   for(int level = minInitLevel; level<=maxInitLevel;level++)
   {
      BOOST_FOREACH(Block3DPtr block, blockVector[level])
      {
         if (block)
         {
            UbTupleDouble3 org          = grid->getBlockWorldCoordinates(block);
            UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
            UbTupleDouble3 nodeOffset   = grid->getNodeOffset(block);
            double         dx           = grid->getDeltaX(block);

            LBMKernelPtr kernel = block->getKernel();
            BCArray3D& bcArray = kernel->getBCProcessor()->getBCArray();          
            DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions(); 
            ShearStressValuesArray3DPtr ssv = kernel->getDataSet()->getShearStressValues();

            int minX1 = 0;
            int minX2 = 0;
            int minX3 = 0;

            int maxX1 = int(distributions->getNX1());
            int maxX2 = int(distributions->getNX2());
            int maxX3 = int(distributions->getNX3());

            for(int ix3=minX3; ix3<maxX3-1; ix3++)
            {
               for(int ix2=minX2; ix2<maxX2-1; ix2++)
               {
                  for(int ix1=minX1; ix1<maxX1-1; ix1++)
                  {
                     if(!bcArray.isUndefined(ix1,ix2,ix3) && !bcArray.isSolid(ix1,ix2,ix3))
                     {
                        //////////////////////////////////////////////////////////////////////////
                        //compute average values
                        //////////////////////////////////////////////////////////////////////////
                        (*ssv)(AvVx,ix1,ix2,ix3) = 0.0;
                        (*ssv)(AvVy,ix1,ix2,ix3) = 0.0;
                        (*ssv)(AvVz,ix1,ix2,ix3) = 0.0;

                        (*ssv)(AvSxx,ix1,ix2,ix3) = 0.0;
                        (*ssv)(AvSyy,ix1,ix2,ix3) = 0.0;
                        (*ssv)(AvSzz,ix1,ix2,ix3) = 0.0;
                        (*ssv)(AvSxy,ix1,ix2,ix3) = 0.0;
                        (*ssv)(AvSyz,ix1,ix2,ix3) = 0.0;
                        (*ssv)(AvSxz,ix1,ix2,ix3) = 0.0;
                        //////////////////////////////////////////////////////////////////////////
                     }
                  }
               }
            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void ShearStressCoProcessor::addInteractor( D3Q27InteractorPtr interactor )
{
   interactors.push_back(interactor);
}
//////////////////////////////////////////////////////////////////////////
void ShearStressCoProcessor::findPlane(int ix1,int ix2,int ix3,Grid3DPtr grid,Block3DPtr block,double &A,double &B,double &C,double &D,double &ii)
{
   double x1plane=0.0,y1plane=0.0,z1plane=0.0;
   double x2plane=0.0,y2plane=0.0,z2plane=0.0;
   double x3plane=0.0,y3plane=0.0,z3plane=0.0;
   BoundaryConditionsPtr bcPtr;
   double dx = grid->getDeltaX(block);
   LBMKernelPtr kernel = block->getKernel();
   DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();
   BCArray3D& bcArray = kernel->getBCProcessor()->getBCArray();
   bcPtr=bcArray.getBC(ix1,ix2,ix3);
   int x,y,z;


   if(InterpolationProcessor::iCellHasSolid(bcArray, ix1, ix2, ix3)) { x = ix1;y = ix2;z = ix3;}   
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1  , ix2-1, ix3  )) { x = ix1+0; y = ix2-1; z = ix3+0;}//S
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1  , ix2  , ix3-1)) { x = ix1+0; y = ix2+0; z = ix3-1;}//B		   
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1-1, ix2  , ix3  )) { x = ix1-1; y = ix2+0; z = ix3+0;}//w
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1  , ix2-1, ix3-1)) { x = ix1+0; y = ix2-1; z = ix3-1;}//BS
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1-1, ix2  , ix3-1)) { x = ix1-1; y = ix2+0; z = ix3-1;}//BW
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1-1, ix2-1, ix3  )) { x = ix1-1; y = ix2-1; z = ix3+0;}//SW

   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1-1, ix2-1, ix3-1)) { x = ix1-1; y = ix2-1; z = ix3-1;}//BSW	  
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1+1, ix2  , ix3  )) { x = ix1+1; y = ix2+0; z = ix3+0;}//E
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1  , ix2+1, ix3  )) { x = ix1+0; y = ix2+1; z = ix3+0;}//N
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1  , ix2  , ix3+1)) { x = ix1+0; y = ix2+0; z = ix3+1;}//T	
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1+1, ix2+1, ix3  )) { x = ix1+1; y = ix2+1; z = ix3+0;}//NE
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1+1, ix2  , ix3+1)) { x = ix1+1; y = ix2+0; z = ix3+1;}//TE
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1  , ix2+1, ix3+1)) { x = ix1+0; y = ix2+1; z = ix3+1;}//TN
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1+1, ix2+1, ix3+1)) { x = ix1+1; y = ix2+1; z = ix3+1;}//TNE

   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1+1, ix2-1, ix3  )) { x = ix1+1; y = ix2-1; z = ix3+0;}//SE
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1-1, ix2+1, ix3  )) { x = ix1-1; y = ix2+1; z = ix3+0;}//NW
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1+1, ix2  , ix3-1)) { x = ix1+1; y = ix2+0; z = ix3-1;}//BE
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1-1, ix2  , ix3+1)) { x = ix1-1; y = ix2+0; z = ix3+1;}//TW
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1+0, ix2+1, ix3-1)) { x = ix1+0; y = ix2+1; z = ix3-1;}//BN
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1+0, ix2-1, ix3+1)) { x = ix1+0; y = ix2-1; z = ix3+1;}//TS

   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1-1, ix2+1, ix3+1)) { x = ix1-1; y = ix2+1; z = ix3+1;}//TNW 
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1+1, ix2-1, ix3+1)) { x = ix1+1; y = ix2-1; z = ix3+1;}//TSE 
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1-1, ix2-1, ix3+1)) { x = ix1-1; y = ix2-1; z = ix3+1;}//TSW 
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1+1, ix2+1, ix3-1)) { x = ix1+1; y = ix2+1; z = ix3-1;}//BNE 
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1-1, ix2+1, ix3-1)) { x = ix1-1; y = ix2+1; z = ix3-1;}//BNW 
   else if(InterpolationProcessor::iCellHasSolid(bcArray, ix1+1, ix2-1, ix3-1)) { x = ix1+1; y = ix2-1; z = ix3-1;}//BSE 


   else   {{UB_THROW( UbException(UB_EXARGS,"there is no cell  ix1="+UbSystem::toString(ix1)+"ix2="+UbSystem::toString(ix2)+"ix3="+UbSystem::toString(ix3)+"GlobalID="+UbSystem::toString(block->getGlobalID())+"dx="+UbSystem::toString(dx)
      +"T="+UbSystem::toString(bcPtr->getQ(D3Q27System::T))+"B="+UbSystem::toString(bcPtr->getQ(D3Q27System::B))
      +"E="+UbSystem::toString(bcPtr->getQ(D3Q27System::E))+"W="+UbSystem::toString(bcPtr->getQ(D3Q27System::W))+"N="+UbSystem::toString(bcPtr->getQ(D3Q27System::N))	
      +"S="+UbSystem::toString(bcPtr->getQ(D3Q27System::S))+"NE="+UbSystem::toString(bcPtr->getQ(D3Q27System::NE))+"SW="+UbSystem::toString(bcPtr->getQ(D3Q27System::SW))
      +"SE="+UbSystem::toString(bcPtr->getQ(D3Q27System::SE))+"NW="+UbSystem::toString(bcPtr->getQ(D3Q27System::NW))    
      +"TE="+
      UbSystem::toString(bcPtr->getQ(D3Q27System::TE))+"BW="+UbSystem::toString(bcPtr->getQ(D3Q27System::BW))+"BE="+UbSystem::toString(bcPtr->getQ(D3Q27System::BE))+"TW="+
      UbSystem::toString(bcPtr->getQ(D3Q27System::TW))+"TN="+UbSystem::toString(bcPtr->getQ(D3Q27System::TN))+"BS="+UbSystem::toString(bcPtr->getQ(D3Q27System::BS))+"BN="+
      UbSystem::toString(bcPtr->getQ(D3Q27System::BN))+"TS="+UbSystem::toString(bcPtr->getQ(D3Q27System::TS))+"TNE="+UbSystem::toString(bcPtr->getQ(D3Q27System::TNE))+"TNW="+
      UbSystem::toString(bcPtr->getQ(D3Q27System::TNW))+"TSE="+UbSystem::toString(bcPtr->getQ(D3Q27System::TSE))+"TSW="+UbSystem::toString(bcPtr->getQ(D3Q27System::TSW))+"BNE="+
      UbSystem::toString(bcPtr->getQ(D3Q27System::BNE))+"BNW="+UbSystem::toString(bcPtr->getQ(D3Q27System::BNW))+"BSE="+UbSystem::toString(bcPtr->getQ(D3Q27System::BSE))+"BSW="+UbSystem::toString(bcPtr->getQ(D3Q27System::BSW)*dx)
      ) ) ;}} 


   if(InterpolationProcessor::iCellHasSolid(bcArray, x, y, z))
   {  
      for(int i = x; i <= x + 1; i++){
         for(int j = y; j <= y + 1; j++){
            for (int k = z; k <= z + 1; k++)
            {
               UbTupleDouble3 pointplane1 =  grid->getNodeCoordinates(block,  i,	  j,	 k);

               double   iph=val<1>(pointplane1);
               double   jph=val<2>(pointplane1);
               double   kph=val<3>(pointplane1);

               if(!bcArray.isSolid(i, j, k))
               {
                  BoundaryConditionsPtr bcPtrIn=bcArray.getBC(i,j,k);
                  if(bcPtrIn)
                  {	 
                     for(int fdir=D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
                     {
                        if( ii<=2)
                        {
                           LBMReal q = bcPtrIn->getQ(fdir);
                           if (q!=999.00000)
                           {	  
                              if     ( fdir==D3Q27System::E ) 
                              {
                                 //if(!bcArray.isSolid(i, j, k))continue;
                                 if (i+q<=x+1)
                                 {
                                    if      (ii==0)	    {	x1plane=iph+q*dx;	  y1plane=jph;	 z1plane=kph;    	ii++;	} 
                                    else if (ii==1) 	{	x2plane=iph+q*dx;     y2plane=jph;	 z2plane=kph;	    if  (x1plane!=x2plane||y1plane!=y2plane||z1plane!=z2plane) 	ii++;	}
                                    else if(ii==2) 	    {   x3plane=iph+q*dx;     y3plane=jph;	 z3plane=kph;       if ((x3plane!=x1plane||y3plane!=y1plane||z3plane!=z1plane)&&(x2plane!=x3plane||y2plane!=y3plane||z2plane!=z3plane)) ii++;} 
                                 }
                              }
                              if     ( fdir==D3Q27System::W ) 
                              {
                                 //if(!bcArray.isSolid(i, j, k))continue;
                                 if (i-q>=x)
                                 {
                                    if      (ii==0)	    {	x1plane=iph-q*dx;	  y1plane=jph;	 z1plane=kph;    	ii++;	} 
                                    else if (ii==1) 	    {	x2plane=iph-q*dx;    y2plane=jph;	 z2plane=kph;	    if  (x1plane!=x2plane||y1plane!=y2plane||z1plane!=z2plane) 	ii++;	}
                                    else if(ii==2) 	    {   x3plane=iph-q*dx;    y3plane=jph;	 z3plane=kph;       if ((x3plane!=x1plane||y3plane!=y1plane||z3plane!=z1plane)&&(x2plane!=x3plane||y2plane!=y3plane||z2plane!=z3plane)) ii++;} 
                                 }
                              }
                              if     ( fdir==D3Q27System::N ) 
                              {
                                 //if(!bcArray.isSolid(i, j, k))continue;
                                 if(j+q<=y+1)
                                 {
                                    if      (ii==0)	    {	x1plane=iph;	y1plane=jph+q*dx;	 z1plane=kph; 	    ii++;	} 
                                    else if (ii==1) 	    {	x2plane=iph;    y2plane=jph+q*dx;	 z2plane=kph;	    if  (x1plane!=x2plane||y1plane!=y2plane||z1plane!=z2plane) 	ii++;	}
                                    else if (ii==2) 	    {   x3plane=iph;    y3plane=jph+q*dx;	 z3plane=kph;       if ((x3plane!=x1plane||y3plane!=y1plane||z3plane!=z1plane)&&(x2plane!=x3plane||y2plane!=y3plane||z2plane!=z3plane)) ii++;} 
                                 }
                              }
                              if     ( fdir==D3Q27System::S ) 
                              {
                                 //if(!bcArray.isSolid(i, j, k))continue;
                                 if (j-q>=y)
                                 {
                                    if      (ii==0)	    {	x1plane=iph;	y1plane=jph-q*dx;	 z1plane=kph; 	ii++;	} 
                                    else if (ii==1) 	    {	x2plane=iph;    y2plane=jph-q*dx;	 z2plane=kph;	if  (x1plane!=x2plane||y1plane!=y2plane||z1plane!=z2plane) 	ii++;	}
                                    else if (ii==2) 	    {   x3plane=iph;    y3plane=jph-q*dx;	 z3plane=kph;   if ((x3plane!=x1plane||y3plane!=y1plane||z3plane!=z1plane)&&(x2plane!=x3plane||y2plane!=y3plane||z2plane!=z3plane)) ii++;} 
                                 }
                              }

                              if     ( fdir==D3Q27System::T ) 
                              {
                                 //if(!bcArray.isSolid(i, j, k))continue;
                                 if(k+q<=z+1)
                                 {
                                    if      (ii==0)	    {	x1plane=iph;	y1plane=jph;	 z1plane=kph+q*dx; 	    ii++;	} 
                                    else if (ii==1) 	    {	x2plane=iph;    y2plane=jph;	 z2plane=kph+q*dx;	    if  (x1plane!=x2plane||y1plane!=y2plane||z1plane!=z2plane) 	ii++;	}
                                    else if (ii==2) 	    {   x3plane=iph;    y3plane=jph;	 z3plane=kph+q*dx;      if ((x3plane!=x1plane||y3plane!=y1plane||z3plane!=z1plane)&&(x2plane!=x3plane||y2plane!=y3plane||z2plane!=z3plane)) ii++;} 
                                 }
                              }
                              if     ( fdir==D3Q27System::B ) 
                              {
                                 //if(!bcArray.isSolid(i, j, k))continue;
                                 if (k-q>=z)
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
   if(ii!=3){

      {{UB_THROW( UbException(UB_EXARGS,"ii is="+UbSystem::toString(ii)+"  ix1="+UbSystem::toString(ix1)+" ix2="+UbSystem::toString(ix2)+" ix3="+UbSystem::toString(ix3)+" Block3D::GlobalID="+UbSystem::toString(block->getGlobalID())+" dx="+UbSystem::toString(dx)
         +" T="+UbSystem::toString(bcPtr->getQ(D3Q27System::T))+" B="+UbSystem::toString(bcPtr->getQ(D3Q27System::B))
         +" E="+UbSystem::toString(bcPtr->getQ(D3Q27System::E))+" W="+UbSystem::toString(bcPtr->getQ(D3Q27System::W))+" N="+UbSystem::toString(bcPtr->getQ(D3Q27System::N))	
         +" S="+UbSystem::toString(bcPtr->getQ(D3Q27System::S))+" NE="+UbSystem::toString(bcPtr->getQ(D3Q27System::NE))+" SW="+UbSystem::toString(bcPtr->getQ(D3Q27System::SW))
         +" SE="+UbSystem::toString(bcPtr->getQ(D3Q27System::SE))+" NW="+UbSystem::toString(bcPtr->getQ(D3Q27System::NW))    
         +" TE="+
         UbSystem::toString(bcPtr->getQ(D3Q27System::TE))+" BW="+UbSystem::toString(bcPtr->getQ(D3Q27System::BW))+" BE="+UbSystem::toString(bcPtr->getQ(D3Q27System::BE))+" TW="+
         UbSystem::toString(bcPtr->getQ(D3Q27System::TW))+" TN="+UbSystem::toString(bcPtr->getQ(D3Q27System::TN))+" BS="+UbSystem::toString(bcPtr->getQ(D3Q27System::BS))+" BN="+
         UbSystem::toString(bcPtr->getQ(D3Q27System::BN))+" TS="+UbSystem::toString(bcPtr->getQ(D3Q27System::TS))+" TNE="+UbSystem::toString(bcPtr->getQ(D3Q27System::TNE))+" TNW="+
         UbSystem::toString(bcPtr->getQ(D3Q27System::TNW))+" TSE="+UbSystem::toString(bcPtr->getQ(D3Q27System::TSE))+" TSW="+UbSystem::toString(bcPtr->getQ(D3Q27System::TSW))+" BNE="+
         UbSystem::toString(bcPtr->getQ(D3Q27System::BNE))+" BNW="+UbSystem::toString(bcPtr->getQ(D3Q27System::BNW))+" BSE="+UbSystem::toString(bcPtr->getQ(D3Q27System::BSE))+" BSW="+UbSystem::toString(bcPtr->getQ(D3Q27System::BSW))
         ) ) ;}}

   }	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
bool ShearStressCoProcessor::checkUndefindedNodes( BCArray3D& bcArray,int ix1,int ix2,int ix3)
{
   for(int i = ix1; i <= ix1 + 1; i++){
      for(int j = ix2; j <= ix2 + 1; j++){
         for (int k = ix3; k <= ix3 + 1; k++)
         {	
            if(bcArray.isUndefined(i, j, k)) return true;
         }  
      }
   }
   return false;
}
//////////////////////////////////////////////////////////////////////////////////////
void ShearStressCoProcessor::initDistance()
{
   BOOST_FOREACH(D3Q27InteractorPtr interactor, interactors)
   {
      typedef std::map<Block3DPtr, std::set< std::vector<int> > > TransNodeIndicesMap;
      BOOST_FOREACH(TransNodeIndicesMap::value_type t, interactor->getTransNodeIndicesMap())
      {
         Block3DPtr block = t.first;
         std::set< std::vector<int> >& transNodeIndicesSet = t.second;

         UbTupleDouble3 org          = grid->getBlockWorldCoordinates(block);
         UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
         UbTupleDouble3 nodeOffset   = grid->getNodeOffset(block);
         double         dx           = grid->getDeltaX(block);

         LBMKernelPtr kernel = block->getKernel();
         BCArray3D& bcArray = kernel->getBCProcessor()->getBCArray();          
         DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions(); 
         ShearStressValuesArray3DPtr ssv = kernel->getDataSet()->getShearStressValues();

         int ghostLayer = kernel->getGhostLayerWidth();
         LBMReal collFactor = kernel->getCollisionFactor();

         int minX1 = ghostLayer;
         int maxX1 = (int)bcArray.getNX1() - 1 - ghostLayer;
         int minX2 = ghostLayer;
         int maxX2 = (int)bcArray.getNX2() - 1 - ghostLayer;
         int minX3 = ghostLayer;
         int maxX3 = (int)bcArray.getNX3() - 1 - ghostLayer;

         BOOST_FOREACH(std::vector<int> node, transNodeIndicesSet)
         {
            int ix1 = node[0];
            int ix2 = node[1];
            int ix3 = node[2];

            //without ghost nodes
            if (ix1 < minX1 || ix1 > maxX1 || ix2 < minX2 || ix2 > maxX2 ||ix3 < minX3 || ix3 > maxX3 ) continue;

            if(bcArray.isFluid(ix1,ix2,ix3) )
            {
               BoundaryConditionsPtr bc = bcArray.getBC(ix1,ix2,ix3);
               if((bc->hasDensityBoundary()||bc->hasVelocityBoundary()))continue;
               int numberOfCorner=0;

               if(bc->getQ(D3Q27System::T)  !=999.000){ numberOfCorner++;}
               if(bc->getQ(D3Q27System::B)  !=999.000){ numberOfCorner++;}
               if(bc->getQ(D3Q27System::E)  !=999.000){ numberOfCorner++;}
               if(bc->getQ(D3Q27System::W)  !=999.000){ numberOfCorner++;}
               if(bc->getQ(D3Q27System::N)  !=999.000){ numberOfCorner++;}
               if(bc->getQ(D3Q27System::S)  !=999.000){ numberOfCorner++;}
               // if(bc->hasVelocityBoundary()||bc->hasDensityBoundary())continue;
               if(numberOfCorner>1)continue;
               if(checkUndefindedNodes( bcArray, ix1,ix2,ix3))continue;


               //////get normal and distance//////
               double A,B,C,D,ii=0.0;
               findPlane(ix1,ix2,ix3,grid,block,A,B,C,D,ii);
               UbTupleDouble3 pointplane1 =  grid->getNodeCoordinates(block, ix1,ix2,ix3);
               double   ix1ph=val<1>(pointplane1);
               double   ix2ph=val<2>(pointplane1);
               double   ix3ph=val<3>(pointplane1);
               double normalDis;
               if(ii!=3)
               {
                  UB_THROW( UbException(UB_EXARGS,"not enough points to create plane"+UbSystem::toString(ii)));
               }
               else
               {    
                  double s = A*ix1ph + B*ix2ph + C*ix3ph + D;//The sign of s = Ax + By + Cz + D determines which side the point (x,y,z) lies with respect to the plane. If s > 0 then the point lies on the same side as the normal (A,B,C). If s < 0 then it lies on the opposite side, if s = 0 then the point (x,y,z) lies on the plane.
                  if (s>0){s=1;} else if (s<0){s=-1;}else {s=0;}

                  normalDis=((A*ix1ph + B*ix2ph + C*ix3ph + D)/sqrt(A*A+B*B+C*C));///distance point to plane xp-Xw=distance
                  normalDis*=s;

                  (*ssv)(normalX1,ix1,ix2,ix3) = A;
                  (*ssv)(normalX2,ix1,ix2,ix3) = B;
                  (*ssv)(normalX3,ix1,ix2,ix3) = C;
                  (*ssv)(normalq,ix1,ix2,ix3) = normalDis;
                  (*ssv)(numberOfPoint,ix1,ix2,ix3) = ii;

               }
            }
         }
      }
   }
}
