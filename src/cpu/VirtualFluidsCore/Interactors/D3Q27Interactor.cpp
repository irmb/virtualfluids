#include "D3Q27Interactor.h"
#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbLogger.h>

#include <basics/writer/WbWriterVtkXmlBinary.h>

#include <geometry3d/GbCuboid3D.h>
#include <geometry3d/GbLine3D.h>
#include "Block3D.h"
#include "Grid3D.h"
#include "BCArray3D.h"
#include "BoundaryConditions.h"
#include "LBMKernel.h"
#include "BCProcessor.h"
#include "BCAdapter.h"

//#include <3rdParty/MarchingCubes/MarchingCubes.h>




using namespace std;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
D3Q27Interactor::D3Q27Interactor() : Interactor3D()
{
   this->reinitWithStoredQsFlag = false;
   this->initRayVectors();
}
//////////////////////////////////////////////////////////////////////////
D3Q27Interactor::D3Q27Interactor(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, int type)
: Interactor3D(geoObject3D, grid, type), relevantForForces(false)
{
   this->reinitWithStoredQsFlag = false;
   this->initRayVectors();
}
//////////////////////////////////////////////////////////////////////////
D3Q27Interactor::D3Q27Interactor(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, SPtr<BCAdapter> bcAdapter,  int type)
   :   Interactor3D(geoObject3D, grid, type), relevantForForces(false)
{
   this->reinitWithStoredQsFlag = false;
   this->addBCAdapter(bcAdapter); 
   this->initRayVectors();
}
//////////////////////////////////////////////////////////////////////////
D3Q27Interactor::D3Q27Interactor(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, SPtr<BCAdapter> bcAdapter,  int type, Interactor3D::Accuracy a)
   :   Interactor3D(geoObject3D, grid, type, a), relevantForForces(false)
{
   this->reinitWithStoredQsFlag = false;
   this->addBCAdapter(bcAdapter); 
   this->initRayVectors();
}
//////////////////////////////////////////////////////////////////////////
D3Q27Interactor::~D3Q27Interactor()
{

}
//////////////////////////////////////////////////////////////////////////
void D3Q27Interactor::initRayVectors()
{
   int fdir; 
   double c1oS2 = UbMath::one_over_sqrt2;
   double c1oS3 = UbMath::one_over_sqrt3;
   fdir = D3Q27System::E;  rayX1[fdir] =  1.0;   rayX2[fdir] =  0.0;   rayX3[fdir] =  0.0;
   fdir = D3Q27System::W;  rayX1[fdir] = -1.0;   rayX2[fdir] =  0.0;   rayX3[fdir] =  0.0;
   fdir = D3Q27System::N;  rayX1[fdir] =  0.0;   rayX2[fdir] =  1.0;   rayX3[fdir] =  0.0;
   fdir = D3Q27System::S;  rayX1[fdir] =  0.0;   rayX2[fdir] = -1.0;   rayX3[fdir] =  0.0;
   fdir = D3Q27System::T;  rayX1[fdir] =  0.0;   rayX2[fdir] =  0.0;   rayX3[fdir] =  1.0;
   fdir = D3Q27System::B;  rayX1[fdir] =  0.0;   rayX2[fdir] =  0.0;   rayX3[fdir] = -1.0;
   fdir = D3Q27System::NE; rayX1[fdir] =  c1oS2; rayX2[fdir] =  c1oS2; rayX3[fdir] =  0.0;
   fdir = D3Q27System::SW; rayX1[fdir] = -c1oS2; rayX2[fdir] = -c1oS2; rayX3[fdir] =  0.0;
   fdir = D3Q27System::SE; rayX1[fdir] =  c1oS2; rayX2[fdir] = -c1oS2; rayX3[fdir] =  0.0;
   fdir = D3Q27System::NW; rayX1[fdir] = -c1oS2; rayX2[fdir] =  c1oS2; rayX3[fdir] =  0.0;
   fdir = D3Q27System::TE; rayX1[fdir] =  c1oS2; rayX2[fdir] = 0.0;    rayX3[fdir] =  c1oS2;
   fdir = D3Q27System::BW; rayX1[fdir] = -c1oS2; rayX2[fdir] = 0.0;    rayX3[fdir] = -c1oS2;
   fdir = D3Q27System::BE; rayX1[fdir] =  c1oS2; rayX2[fdir] = 0.0;    rayX3[fdir] = -c1oS2;
   fdir = D3Q27System::TW; rayX1[fdir] = -c1oS2; rayX2[fdir] = 0.0;    rayX3[fdir] =  c1oS2;
   fdir = D3Q27System::TN; rayX1[fdir] =  0.0;   rayX2[fdir] = c1oS2;  rayX3[fdir] =  c1oS2;
   fdir = D3Q27System::BS; rayX1[fdir] =  0.0;   rayX2[fdir] =-c1oS2;  rayX3[fdir] = -c1oS2;
   fdir = D3Q27System::BN; rayX1[fdir] =  0.0;   rayX2[fdir] = c1oS2;  rayX3[fdir] = -c1oS2;
   fdir = D3Q27System::TS; rayX1[fdir] =  0.0;   rayX2[fdir] =-c1oS2;  rayX3[fdir] =  c1oS2;

   fdir = D3Q27System::TNW; rayX1[fdir] = -c1oS3; rayX2[fdir] = c1oS3; rayX3[fdir] =   c1oS3;
   fdir = D3Q27System::TNE; rayX1[fdir] =  c1oS3; rayX2[fdir] = c1oS3; rayX3[fdir] =   c1oS3;
   fdir = D3Q27System::TSW; rayX1[fdir] = -c1oS3; rayX2[fdir] =-c1oS3; rayX3[fdir] =   c1oS3;
   fdir = D3Q27System::TSE; rayX1[fdir] =  c1oS3; rayX2[fdir] =-c1oS3; rayX3[fdir] =   c1oS3;
   fdir = D3Q27System::BNW; rayX1[fdir] = -c1oS3; rayX2[fdir] = c1oS3; rayX3[fdir] =  -c1oS3;
   fdir = D3Q27System::BNE; rayX1[fdir] =  c1oS3; rayX2[fdir] = c1oS3; rayX3[fdir] =  -c1oS3;
   fdir = D3Q27System::BSW; rayX1[fdir] = -c1oS3; rayX2[fdir] =-c1oS3; rayX3[fdir] =  -c1oS3;
   fdir = D3Q27System::BSE; rayX1[fdir] =  c1oS3; rayX2[fdir] =-c1oS3; rayX3[fdir] =  -c1oS3;
}
//////////////////////////////////////////////////////////////////////////
void D3Q27Interactor::initInteractor(const double& timeStep)
{
   UBLOG(logDEBUG5,"D3Q27Interactor::initInteractor - "<<" for timestep = "<<timeStep);
   
   //////////////////////////////////////////////////////////////////////////
   //init bcs
   int nofAdapter = (int)bcAdapters.size();
   if(nofAdapter==0) UBLOG(logWARNING,"WARNING - D3Q27Interactor::initInteractor Warning - no nodeAdapter available");
   bool needTimeDependence = false;
   for(int pos=0; pos<nofAdapter; ++pos)
   {
      bcAdapters[pos]->init(this,timeStep);
      if(bcAdapters[pos]->isTimeDependent()) needTimeDependence = true;
   }
   if(needTimeDependence) this->setTimeDependent();
   else                   this->unsetTimeDependent();

   Interactor3D::initInteractor(timeStep);

  ////calcForces arbeitet nicht korrekt, wenn Geo mit Bloecken 
  // //unterschiedlicher Leveltiefe diskretisiert -> exception
  // //abfrage steht hier, weil es theoretisch sein kann, dass bei parallelen rechnungen
  // //genau der block mit dem anderen level auf einem anderen prozess liegt...
  // //Update: es kann u.U. passieren, dass Bl�cke in der Liste nicht aktiv sin
  // //(falls diese z.B. duch andere Interactoren solid gesetzt wurden)
  // //diese werden nicht ber�cksichtigt (auch nicht beid er kraftauswertung sp�ter)
  // if(this->isRelevantForForces() )
  // {
  //    int level = -1;     
  //    for( vector<Block3D*>::const_iterator pos = transBlockSet->begin(); pos!=transBlockSet->end(); ++pos)
  //       if( (*pos)->isActive() )
  //       {
  //          level = (*pos)->getLevel();
  //          break;
  //       }
  //    
  //    bool check = false;
  //    for( vector<Block3D*>::const_iterator pos = transBlockSet->begin(); pos!=transBlockSet->end(); ++pos)
  //       if( (*pos)->isActive() && (*pos)->getLevel()!=level)
  //       {
  //          throw UbException(UB_EXARGS,"interactor is relevant for forces,"
  //                                  +(string)" but has transblocks with different levels"
  //                                  +(string)" -> not supportet by D3Q27Interactor::getForces()"
  //                                  +(string)" -> increase refineWidth");
  //       }
  // }
}
//////////////////////////////////////////////////////////////////////////
void D3Q27Interactor::updateInteractor(const double& timestep)
{
   //UB_THROW( UbException(UB_EXARGS,"toDo") );
   //if(this->isFluid()) return;

   UBLOG(logDEBUG5,"D3Q27Interactor::updateInteractor - for timestep = "<<timestep);

   //////////////////////////////////////////////////////////////////////////
   //update bcs
   int nofAdapter = (int)bcAdapters.size();
   if(nofAdapter==0) UBLOG(logERROR,"WARNING - D3Q27Interactor::updateInteractor Warning - no nodeAdapter available for ");

   bool needTimeDependence = false;

   for(int pos=0; pos<nofAdapter; ++pos)
   {
      bcAdapters[pos]->update(this,timestep);
      if(bcAdapters[pos]->isTimeDependent()) needTimeDependence = true;
   }
   if(needTimeDependence) this->setTimeDependent();
   else                   this->unsetTimeDependent();
   
   //UBLOG(logINFO, "transNodeIndicesMap = "<<transNodeIndicesMap.size());
   
   for(BcNodeIndicesMap::value_type t : bcNodeIndicesMap)
   {
      SPtr<Block3D> block = t.first;
      std::set< std::vector<int> >& transNodeIndicesSet = t.second;

      //UBLOG(logINFO, "transNodeIndicesSet = "<<transNodeIndicesSet.size());

      if(block->isNotActive() || !block) continue;

      SPtr<ILBMKernel> kernel = block->getKernel();
      SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();

      set< std::vector<int> >::iterator setPos;

      for(setPos=transNodeIndicesSet.begin(); setPos!=transNodeIndicesSet.end();  ++setPos)
      {
         int    x1      = (*setPos)[0];
         int    x2      = (*setPos)[1];
         int    x3      = (*setPos)[2];
         Vector3D coords = grid.lock()->getNodeCoordinates(block, x1, x2, x3);
         double worldX1 = coords[0];
         double worldX2 = coords[1];
         double worldX3 = coords[2];

         SPtr<BoundaryConditions> bc = bcArray->getBC(x1,x2,x3);
         if(bc) //kann sein, dass die BC durch das solid setzen eines andern interactors geloescht wurde
         {
            for(size_t i=0; i<bcAdapters.size(); i++)
               bcAdapters[i]->adaptBC(*this,bc,worldX1,worldX2,worldX3,timestep);
         }
         //else 
         //{
         //   UBLOG(logERROR,"D3Q27Interactor.updateInteractor (Z."<<__LINE__<<"): da ist kein BC dran ... kann aber korrekt sein s.h. code");
         //}
      }
   }
}
//////////////////////////////////////////////////////////////////////////
// Berechnung findet im realen Koordinatensystem statt !!!
// nicht im normierten !
//x1,x2,x3 sind die Koordinaten unten links vom "System"
//extendedBoundingGeoOfGeoObject MUSS bereits um delta_x_level in jede richtung vergroesert worden sein fuer SOLID
bool D3Q27Interactor::setDifferencesToGbObject3D(const SPtr<Block3D> block/*,const double& orgX1,const double& orgX2,const double& orgX3,const double& blockLengthX1,const double& blockLengthX2,const double& blockLengthX3, const double& timestep*/)
{
   if(!block) return false;

   if(block->isNotActive()) return false;//continue;

   bcNodeIndicesMap[block] = set< std::vector<int> >();
   set< std::vector<int> >& transNodeIndices = bcNodeIndicesMap[block];
   solidNodeIndicesMap[block] = set< UbTupleInt3 >();
   set< UbTupleInt3 >& solidNodeIndices = solidNodeIndicesMap[block];


   double timestep = 0;
   bool oneEntryGotBC = false; //ob ueberhaupt ein eintrag ein BC zugewiesen wurde
   bool gotQs         = false; //true, wenn "difference" gesetzt wurde
   SPtr<BoundaryConditions> bc;

   SPtr<ILBMKernel> kernel = block->getKernel();
   SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();

   double internX1,internX2,internX3;


   //width of ghost layer 
   int gl = kernel->getGhostLayerWidth();

   int startIX1 = 0;
   int startIX2 = 0;
   int startIX3 = 0; 
   int stopIX1  = (int)bcArray->getNX1();
   int stopIX2  = (int)bcArray->getNX2();
   int stopIX3  = (int)bcArray->getNX3(); 

   double         dx       = grid.lock()->getDeltaX(block);
   UbTupleDouble3 orgDelta = grid.lock()->getNodeOffset(block);

   //anderes boundingRect als in init, da hier das boundrect um ein dx vergroessert werden muss
   GbCuboid3D extendedBoundingGeoOfGeoObject(  geoObject3D->getX1Minimum()-1.02*dx 
                                             , geoObject3D->getX2Minimum()-1.02*dx 
                                             , geoObject3D->getX3Minimum()-1.02*dx 
                                             , geoObject3D->getX1Maximum()+1.02*dx 
                                             , geoObject3D->getX2Maximum()+1.02*dx 
                                             , geoObject3D->getX3Maximum()+1.02*dx );

   double deltaX1=dx, deltaX2 =dx, deltaX3=dx;

   if(   geoObject3D->hasRaytracing() 
        || (this->isInverseSolid() && geoObject3D->raytracingSupportsPointsInside() ) )
   {
      //wenn deltaX1==deltaX2==deltaX3 (muss fuer LB!!)
      if(!UbMath::zero( deltaX1-deltaX2 + deltaX1-deltaX3 + deltaX2-deltaX3 ) )
         throw UbException(UB_EXARGS,"fuer den bei LB nicht vorkommenden Fall deltaX1!=deltaX2!=deltaX3  nicht implementiert ");

      vector<double> distNeigh(D3Q27System::FENDDIR+1, UbMath::sqrt2*deltaX1);
      distNeigh[D3Q27System::E] = distNeigh[D3Q27System::W] = distNeigh[D3Q27System::N] = deltaX1;
      distNeigh[D3Q27System::S] = distNeigh[D3Q27System::T] = distNeigh[D3Q27System::B] = deltaX1;
      distNeigh[D3Q27System::NE]  = distNeigh[D3Q27System::NW]  = distNeigh[D3Q27System::SW]  = distNeigh[D3Q27System::SE]  = UbMath::sqrt2*deltaX1;
      distNeigh[D3Q27System::TE]  = distNeigh[D3Q27System::TN]  = distNeigh[D3Q27System::TW]  = distNeigh[D3Q27System::TS]  = UbMath::sqrt2*deltaX1;
      distNeigh[D3Q27System::BE]  = distNeigh[D3Q27System::BN]  = distNeigh[D3Q27System::BW]  = distNeigh[D3Q27System::BS]  = UbMath::sqrt2*deltaX1;
      distNeigh[D3Q27System::TNE] = distNeigh[D3Q27System::TNW] = distNeigh[D3Q27System::TSE] = distNeigh[D3Q27System::TSW] = UbMath::sqrt3*deltaX1;
      distNeigh[D3Q27System::BNE] = distNeigh[D3Q27System::BNW] = distNeigh[D3Q27System::BSE] = distNeigh[D3Q27System::BSW] = UbMath::sqrt3*deltaX1;
      double q;
      bool pointOnBoundary = false;

//#ifdef _OPENMP
//      #pragma omp parallel for private(internX1,internX2,internX3,gotQs,bc,q )
//#endif
      for(int ix3=startIX3; ix3<stopIX3; ix3++)
      {
         for(int ix2=startIX2; ix2<stopIX2; ix2++)
         {
            for(int ix1=startIX1; ix1<stopIX1; ix1++)
            {
               //TODO weiter untersuchen, ob das nicht ein Fehler ist
               if(bcArray->isUndefined(ix1, ix2, ix3)) continue;

               Vector3D coords = grid.lock()->getNodeCoordinates(block, ix1, ix2, ix3);
               internX1 = coords[0];
               internX2 = coords[1];
               internX3 = coords[2];

               //point in object test ist ueberfluessig, weil die start und stop indices bereits zuvor 
               //ermittelt werden -> es werden nur point-IN-cube indizes betrachtet
               if(extendedBoundingGeoOfGeoObject.isPointInGbObject3D(internX1,internX2,internX3))
               {
                  if(this->isSolid() )
                  {
                     if(this->geoObject3D->isPointInGbObject3D(internX1, internX2, internX3))
                     {
//#ifdef _OPENMP
//                        #pragma omp critical (SOLIDNODE_SET_CHANGE)
//#endif
                        {
                           solidNodeIndices.insert(UbTupleInt3(ix1, ix2, ix3));
                           bcArray->setSolid(ix1,ix2,ix3); 
                        }
                        continue;
                     }
                  }
                  else if( this->isInverseSolid()  )
                  {
                     //bei inverse solid sind alle Knoten AUSSERHALB und auf der boundary SOLID
                     if( !this->geoObject3D->isPointInGbObject3D(internX1, internX2, internX3, pointOnBoundary) 
                        || pointOnBoundary == true )
                     {
//#ifdef _OPENMP
//                        #pragma omp critical (SOLID_SET_CHANGE)
//#endif
                        {
                           solidNodeIndices.insert(UbTupleInt3(ix1, ix2, ix3));
                           bcArray->setSolid(ix1,ix2,ix3); 
                        }
                        continue;
                     }
                  }

                  //evtl wurde node von anderen interactoren solid gesetzt (muss hie rein->sonst evtl bei
                  //ueberschneidender geo -> solidNodeIndicesMap unvollstaendig)
                  if(bcArray->isSolid(ix1,ix2,ix3)) 
                     continue;

                  gotQs = false;

                  //TODO: prüfen was passiert wenn ein Geoobjekt zwischen zwei knoten rausguckt
                  //  * /
                  //<
                  //  * \  //
                  //sollen dann keine qs gesetzt werden
                  for(int fdir=D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
                  {
                     q = geoObject3D->getIntersectionRaytraceFactor(internX1,internX2,internX3,rayX1[fdir],rayX2[fdir],rayX3[fdir]);
                     q /= distNeigh[fdir];

                     //assert(UbMath::lessEqual(q, 1.0));

                     if( UbMath::inClosedInterval(q, 1.0, 1.0) ) q = 1.0;
                     if( UbMath::greater(q, 0.0) && UbMath::lessEqual(q, 1.0) )
                     {
                        //#pragma omp critical (BC_CHANGE)
                        {
                           bc = bcArray->getBC(ix1,ix2,ix3);
                           if(!bc)
                           {
                              //bc = bvd->createD3Q27BoundaryCondition(); //= new D3Q27BoundaryCondition();
                              bc = SPtr<BoundaryConditions>(new BoundaryConditions);
                              bcArray->setBC(ix1,ix2,ix3,bc);
                           }
                     //TODO: man muss ueberlegen, wie kann man, dass die Geschwindigkeit auf 0.0 gesetzt werden, vermeiden
                     //das folgt zu unguenstigen "design rules"
                      //SG 31.08.2010 das Problem - bewegter Interactor angrenzend an stehenden noslip interactor
                      // hier sollte die Geschwindigkeit auf 0.0 gesetzt werden
                           if(bc->hasNoSlipBoundary())
                           {
                              bc->setBoundaryVelocityX1(0.0);
                              bc->setBoundaryVelocityX2(0.0);
                              bc->setBoundaryVelocityX3(0.0);
                           }
                       //SG 31.08.2010
                        
                           for(int index=(int)bcAdapters.size()-1; index>=0; --index)
                              bcAdapters[index]->adaptBCForDirection(*this,bc,internX1,internX2,internX3,q,fdir, timestep);
                        }

                        gotQs=true;
                     }
                  }

                  if(gotQs)
                  {
//#ifdef _OPENMP
//                     #pragma omp critical (TRANSNODE_SET_CHANGE)
//#endif
                     {
                        oneEntryGotBC = true;

                        std::vector<int> p(3);
                        p[0]=ix1; p[1]=ix2; p[2]=ix3;
                        transNodeIndices.insert(p);
                     
                        for(int index=(int)bcAdapters.size()-1; index>=0; --index)
                           bcAdapters[index]->adaptBC(*this,bc,internX1,internX2,internX3, timestep);
                     }
                  }
               }
               else if( this->isInverseSolid()  )
               {
                  //bei inverse solid sind alle Knoten AUSSERHALB und auf der boundary SOLID
                  if( !this->geoObject3D->isPointInGbObject3D(internX1, internX2, internX3, pointOnBoundary) 
                     || pointOnBoundary == true )
                 {
//#ifdef _OPENMP 
//   #pragma omp critical (SOLID_SET_CHANGE)
//#endif
                     {
                        solidNodeIndices.insert(UbTupleInt3(ix1, ix2, ix3));
                        bcArray->setSolid(ix1,ix2,ix3); 
                     }
                     continue;
                  }
               }
            }
         }
      }
   }
   else  //clipping -> langsamer (wird derzeit auch fuer alle inverseSolid objekte verwendet deren raytracing nicht fuer nodes INNERHALB der geo funzt)
   {
      bool pointOnBoundary = false;
      //#pragma omp parallel for private(internX1,internX2,internX3,gotQs,bc,pointOnBoundary )
      for(int ix1=startIX1; ix1<stopIX1; ix1++)
      {
         for(int ix2=startIX2; ix2<stopIX2; ix2++)
         {
            for(int ix3=startIX3; ix3<stopIX3; ix3++)
            {
               if(bcArray->isSolid(ix1,ix2,ix3) || bcArray->isUndefined(ix1, ix2, ix3)) continue;

               Vector3D coords = grid.lock()->getNodeCoordinates(block, ix1, ix2, ix3);
               internX1 = coords[0];
               internX2 = coords[1];
               internX3 = coords[2];

               if(extendedBoundingGeoOfGeoObject.isPointInGbObject3D(internX1,internX2,internX3))
               {
                  if( this->isSolid() && this->geoObject3D->isPointInGbObject3D(internX1, internX2, internX3) )
                  {
                     //#pragma omp critical (SOLID_SET_CHANGE)
                     {
                        solidNodeIndices.insert(UbTupleInt3(ix1, ix2, ix3));
                        bcArray->setSolid(ix1,ix2,ix3); 
                     }
                     continue;
                  }
                  else if( this->isInverseSolid()  )
                  {
                     //bei inverse solid sind alle Knoten AUSSERHALB und auf der boundary SOLID
                     if( !this->geoObject3D->isPointInGbObject3D(internX1, internX2, internX3, pointOnBoundary) 
                         || pointOnBoundary == true )
                     {
                        //#pragma omp critical (SOLID_SET_CHANGE)
                        {
                           solidNodeIndices.insert(UbTupleInt3(ix1, ix2, ix3));
                           bcArray->setSolid(ix1,ix2,ix3); 
                        }
                        continue;
                     }
                  }

                  gotQs = false;

                  GbPoint3D pointA(internX1,internX2,internX3);
                  for(int fdir=D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
                  {
                     double x1B = internX1+D3Q27System::DX1[fdir]*deltaX1;
                     double x2B = internX2+D3Q27System::DX2[fdir]*deltaX2;
                     double x3B = internX3+D3Q27System::DX3[fdir]*deltaX3;

                     GbPoint3D pointB(x1B,x2B,x3B);
                     GbLine3D* clippedLine = this->geoObject3D->createClippedLine3D(pointA, pointB);

                     if(clippedLine)
                     {
                        double q=0.0;
                        if( !this->isInverseSolid() )  //A liegt auf jeden Fall aussen
                        {
                           double distanceAB = pointA.getDistance(&pointB); //pointA to B
                           double distanceAP = UbMath::min(pointA.getDistance(clippedLine->getPoint1()),
                                                           pointA.getDistance(clippedLine->getPoint2()) );
                           q = distanceAP/distanceAB;
                        }
                        else
                        {
                           bool pointIsOnBoundary = false;
                           if(   !clippedLine->getPoint1()->equals(&pointB)
                              && !clippedLine->getPoint2()->equals(&pointB) )
                           {
                              //A liegt auf jeden Fall drinnen, clipped line darf B nicht enthalten
                              double distanceAB = pointA.getDistance(&pointB); //pointA to B
                              double distanceAP = clippedLine->getLength();
                              q = distanceAP/distanceAB;
                           }
                           else if(  this->geoObject3D->isPointInGbObject3D( pointB.getX1Coordinate()
                                                                            ,pointB.getX2Coordinate()
                                                                            ,pointB.getX3Coordinate()
                                                                            ,pointIsOnBoundary )
                                   && pointIsOnBoundary )
                           {
                              //A liegt auf jeden Fall drinnen, B liegt genau auf ObjektBoundary => q=1.0
                              q=1.0;
                           }
                           else
                           {
                              q = 0.0;
                           }
                        }

                        if(UbMath::inClosedInterval(q, 1.0, 1.0)) q = 1.0;
                        if(UbMath::lessEqual(q, 1.0) && UbMath::greater(q, 0.0))
                        {
                           //#pragma omp critical (BC_CHANGE)
                           {
                              bc = bcArray->getBC(ix1,ix2,ix3);
                              if(!bc)
                              {
                                 //bc = bvd->createD3Q27BoundaryCondition(); //= new D3Q27BoundaryCondition();
                                 bc = SPtr<BoundaryConditions>(new BoundaryConditions);
                                 bcArray->setBC(ix1,ix2,ix3,bc);
                              }
                              for(int index=(int)bcAdapters.size()-1; index>=0; --index)
                                 bcAdapters[index]->adaptBCForDirection(*this,bc,internX1,internX2,internX3,q,fdir,timestep);
                           }
                                          
                           gotQs=true;
                        }

                        clippedLine->deletePoint1();
                        clippedLine->deletePoint2();
                        delete clippedLine;
                     }
                  }

                  if(gotQs)
                  {
                     //#pragma omp critical (TRANSNODE_SET_CHANGE)
                     {
                        oneEntryGotBC = true;

                        std::vector<int> p(3);
                        p[0]=ix1; p[1]=ix2; p[2]=ix3;
                        transNodeIndices.insert(p);

                        for(int index=(int)bcAdapters.size()-1; index>=0; --index)
                           bcAdapters[index]->adaptBC(*this,bc,internX1,internX2,internX3,timestep);
                     }
                  }
               }
            }
         }
      }
   }

   return oneEntryGotBC;
}
//////////////////////////////////////////////////////////////////////////
void D3Q27Interactor::addQsLineSet(std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt2 >& lines)
{
      for(SPtr<Block3D> block : bcBlocks)
      {
         if(!block) continue;

         double         dx       = grid.lock()->getDeltaX(block);
         UbTupleDouble3 orgDelta = grid.lock()->getNodeOffset(block);

         SPtr<ILBMKernel> kernel = block->getKernel();
         SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();

         map<SPtr<Block3D>, set< std::vector<int> > >::iterator pos = bcNodeIndicesMap.find(block);
         if(pos==bcNodeIndicesMap.end()) 
         {
            UB_THROW( UbException(UB_EXARGS,"block nicht in indizes map!!!") );
         }
         set< std::vector<int> >& transNodeIndicesSet = pos->second;
         set< std::vector<int> >::iterator setPos;

         std::size_t node1Index, node2Index;

         UbTupleDouble3 blockOrg = grid.lock()->getBlockWorldCoordinates(block);

         for(setPos=transNodeIndicesSet.begin(); setPos!=transNodeIndicesSet.end();  ++setPos)
         {
            int ix1 = (*setPos)[0];
            int ix2 = (*setPos)[1];
            int ix3 = (*setPos)[2];

            if(bcArray->isFluid(ix1,ix2,ix3)) //es kann sein, dass der node von einem anderen interactor z.B. als solid gemarkt wurde!!!
            {
               if( !bcArray->hasBC(ix1,ix2,ix3) ) continue;
               SPtr<BoundaryConditions> bc = bcArray->getBC(ix1,ix2,ix3);

               double x1a = val<1>(blockOrg) - val<1>(orgDelta) + ix1 * dx;
               double x2a = val<2>(blockOrg) - val<2>(orgDelta) + ix2 * dx;
               double x3a = val<3>(blockOrg) - val<3>(orgDelta) + ix3 * dx;
               nodes.push_back( makeUbTuple( (float)x1a, (float)x2a, (float)x3a ) );
               node1Index = nodes.size()-1;

               for(int dir = D3Q27System::FSTARTDIR; dir<=D3Q27System::FENDDIR; dir++)
               {
                  if (bc->hasBoundaryConditionFlag(D3Q27System::INVDIR[dir]))
                  {
                     double x1b, x2b, x3b, q = bc->getQ(dir);
                     switch(dir)
                     {
                     case D3Q27System::E : x1b = x1a+q*dx; x2b = x2a         ; x3b = x3a         ; break;
                     case D3Q27System::N : x1b = x1a         ; x2b = x2a+q*dx; x3b = x3a         ; break;
                     case D3Q27System::W : x1b = x1a-q*dx; x2b = x2a         ; x3b = x3a         ; break;
                     case D3Q27System::S : x1b = x1a         ; x2b = x2a-q*dx; x3b = x3a         ; break;
                     case D3Q27System::NE: x1b = x1a+q*dx; x2b = x2a+q*dx; x3b = x3a         ; break;
                     case D3Q27System::NW: x1b = x1a-q*dx; x2b = x2a+q*dx; x3b = x3a         ; break;
                     case D3Q27System::SW: x1b = x1a-q*dx; x2b = x2a-q*dx; x3b = x3a         ; break;
                     case D3Q27System::SE: x1b = x1a+q*dx; x2b = x2a-q*dx; x3b = x3a         ; break;
                     case D3Q27System::T : x1b = x1a         ; x2b = x2a         ; x3b = x3a+q*dx; break;
                     case D3Q27System::TE: x1b = x1a+q*dx; x2b = x2a         ; x3b = x3a+q*dx; break;
                     case D3Q27System::TN: x1b = x1a         ; x2b = x2a+q*dx; x3b = x3a+q*dx; break;
                     case D3Q27System::TW: x1b = x1a-q*dx; x2b = x2a         ; x3b = x3a+q*dx; break;
                     case D3Q27System::TS: x1b = x1a         ; x2b = x2a-q*dx; x3b = x3a+q*dx; break;
                     case D3Q27System::B : x1b = x1a         ; x2b = x2a         ; x3b = x3a-q*dx; break;
                     case D3Q27System::BE: x1b = x1a+q*dx; x2b = x2a         ; x3b = x3a-q*dx; break;
                     case D3Q27System::BN: x1b = x1a         ; x2b = x2a+q*dx; x3b = x3a-q*dx; break;
                     case D3Q27System::BW: x1b = x1a-q*dx; x2b = x2a         ; x3b = x3a-q*dx; break;
                     case D3Q27System::BS: x1b = x1a         ; x2b = x2a-q*dx; x3b = x3a-q*dx; break;
                     case D3Q27System::TNE : x1b = x1a+q*dx  ; x2b = x2a+q*dx; x3b = x3a+q*dx; break;
                     case D3Q27System::BSW : x1b = x1a-q*dx  ; x2b = x2a-q*dx; x3b = x3a-q*dx; break;
                     case D3Q27System::BNE : x1b = x1a+q*dx  ; x2b = x2a+q*dx; x3b = x3a-q*dx; break;
                     case D3Q27System::TSW : x1b = x1a-q*dx  ; x2b = x2a-q*dx; x3b = x3a+q*dx; break;
                     case D3Q27System::TSE : x1b = x1a+q*dx  ; x2b = x2a-q*dx; x3b = x3a+q*dx; break;
                     case D3Q27System::BNW : x1b = x1a-q*dx  ; x2b = x2a+q*dx; x3b = x3a-q*dx; break;
                     case D3Q27System::BSE : x1b = x1a+q*dx  ; x2b = x2a-q*dx; x3b = x3a-q*dx; break;
                     case D3Q27System::TNW : x1b = x1a-q*dx  ; x2b = x2a+q*dx; x3b = x3a+q*dx; break;
                     default: throw UbException(UB_EXARGS,"unknown direction");
                     }

                     nodes.push_back( makeUbTuple( (float)x1b, (float)x2b, (float)x3b ) );
                     node2Index = nodes.size()-1;

                     lines.push_back( makeUbTuple( (int)node1Index, (int)node2Index) );
                  }
               }
            }
         }
      }
}
////////////////////////////////////////////////////////////////////////////
vector< pair<GbPoint3D,GbPoint3D> >  D3Q27Interactor::getQsLineSet()
{
   vector< pair<GbPoint3D,GbPoint3D> >  QsLineSet;
   pair<GbPoint3D,GbPoint3D> pointpair;

   UbTupleInt3 blocknx = grid.lock()->getBlockNX();

   int blocknx1 = val<1>(blocknx);
   int blocknx2 = val<2>(blocknx);
   int blocknx3 = val<3>(blocknx);
   //   vector<double> deltaT = grid->getD3Q27Calculator()->getDeltaT();

   for(SPtr<Block3D> block : bcBlocks)
   {
       SPtr<ILBMKernel> kernel = block->getKernel();
      SPtr<BCArray3D> bcMatrix = kernel->getBCProcessor()->getBCArray();
      UbTupleDouble3 nodeOffset   = grid.lock()->getNodeOffset(block);

      //double collFactor = ((LbD3Q27Calculator*)grid->getCalculator())->getCollisionsFactors()[block->getLevel()];
      //checken ob obere reihe doppelt im system vorhanden oder nicht
      bool include_N_Face  = false; //x1=[0..blocknx1[ && x3=[0..blocknx3[
      bool include_E_Face  = false; //x2=[0..blocknx2[ && x3=[0..blocknx3[
      bool include_T_Face  = false; //x1=[0..blocknx1[ && x2=[0..blocknx2[
      bool include_NE_Edge = false; //(x1/x2/x3)=(blocknx1/blocknx2/[0..blocknx3[)
      bool include_TN_Edge = false; //(x1/x2/x3)=([0..blocknx1[/blocknx2/blocknx1)
      bool include_TE_Edge = false; //(x1/x2/x3)=(blocknx1/[0..blocknx2[/blocknx2)
      if(block)
      {
         if( !block->getConnector(D3Q27System::N ) ) include_N_Face  = true;
         if( !block->getConnector(D3Q27System::E ) ) include_E_Face  = true;
         if( !block->getConnector(D3Q27System::T ) ) include_T_Face  = true;
         if( !block->getConnector(D3Q27System::NE) && include_N_Face && include_E_Face ) include_NE_Edge = true;
         if( !block->getConnector(D3Q27System::TN) && include_T_Face && include_N_Face ) include_TN_Edge = true;
         if( !block->getConnector(D3Q27System::TE) && include_T_Face && include_E_Face ) include_TE_Edge = true;
      }

      //      double dT = deltaT[block->getLevel()];

      map<SPtr<Block3D>, set< std::vector<int> > >::iterator pos = bcNodeIndicesMap.find(block);
      if(pos==bcNodeIndicesMap.end()) throw UbException(UB_EXARGS,"block nicht in indizes map!!!"+block->toString());
      set< std::vector<int> >& transNodeIndicesSet = pos->second;
      set< std::vector<int> >::iterator setPos;

      double x1,x2,x3,dx;
      grid.lock()->calcStartCoordinatesAndDelta(block,x1,x2,x3,dx);

      //cout<<"getQs: "<<transBlockSet->size()<<" "<<transNodeIndicesVec.size()<<endl;

      for(setPos=transNodeIndicesSet.begin(); setPos!=transNodeIndicesSet.end();  ++setPos)
      {
         int ix1 = (*setPos)[0];
         int ix2 = (*setPos)[1];
         int ix3 = (*setPos)[2];

         if(   ( ix1<blocknx1 && ix2<blocknx2 && ix3<blocknx3 ) //std fall
            || ( include_E_Face  && ix1==blocknx1 && ix2<blocknx2  && ix3<blocknx3  )
            || ( include_N_Face  && ix2==blocknx2 && ix1<blocknx1  && ix3<blocknx3  )
            || ( include_T_Face  && ix3==blocknx3 && ix1<blocknx1  && ix2<blocknx2  )
            || ( include_NE_Edge && ix1==blocknx1 && ix2==blocknx2 )
            || ( include_TN_Edge && ix2==blocknx2 && ix3==blocknx3 )
            || ( include_TE_Edge && ix1==blocknx1 && ix3==blocknx3 ) ) //ansonsten doppelt im kraftwert
         {
            if(bcMatrix->isFluid(ix1,ix2,ix3)) //es kann sein, dass der node von einem anderen interactor z.B. als solid gemarkt wurde!!!
            {
               if( !bcMatrix->hasBC(ix1,ix2,ix3) ) continue;
               SPtr<BoundaryConditions> bc = bcMatrix->getBC(ix1,ix2,ix3);
               double x1a = x1-val<1>(nodeOffset)+dx * ix1;
               double x2a = x2-val<2>(nodeOffset)+dx * ix2;
               double x3a = x3-val<3>(nodeOffset)+dx * ix3;
               pointpair.first.setX1(x1a);
               pointpair.first.setX2(x2a);
               pointpair.first.setX3(x3a);
               for(int dir = D3Q27System::FSTARTDIR; dir<=D3Q27System::FENDDIR; dir++)
               {
                  if (bc->hasBoundaryConditionFlag(D3Q27System::INVDIR[dir]))
                  {
                     double x1b, x2b, x3b, q = bc->getQ(dir);
                     switch(dir)
                     {
                     case D3Q27System::E : x1b = x1a+q*dx; x2b = x2a         ; x3b = x3a         ; break;
                     case D3Q27System::N : x1b = x1a         ; x2b = x2a+q*dx; x3b = x3a         ; break;
                     case D3Q27System::W : x1b = x1a-q*dx; x2b = x2a         ; x3b = x3a         ; break;
                     case D3Q27System::S : x1b = x1a         ; x2b = x2a-q*dx; x3b = x3a         ; break;
                     case D3Q27System::NE: x1b = x1a+q*dx; x2b = x2a+q*dx; x3b = x3a         ; break;
                     case D3Q27System::NW: x1b = x1a-q*dx; x2b = x2a+q*dx; x3b = x3a         ; break;
                     case D3Q27System::SW: x1b = x1a-q*dx; x2b = x2a-q*dx; x3b = x3a         ; break;
                     case D3Q27System::SE: x1b = x1a+q*dx; x2b = x2a-q*dx; x3b = x3a         ; break;
                     case D3Q27System::T : x1b = x1a         ; x2b = x2a         ; x3b = x3a+q*dx; break;
                     case D3Q27System::TE: x1b = x1a+q*dx; x2b = x2a         ; x3b = x3a+q*dx; break;
                     case D3Q27System::TN: x1b = x1a         ; x2b = x2a+q*dx; x3b = x3a+q*dx; break;
                     case D3Q27System::TW: x1b = x1a-q*dx; x2b = x2a         ; x3b = x3a+q*dx; break;
                     case D3Q27System::TS: x1b = x1a         ; x2b = x2a-q*dx; x3b = x3a+q*dx; break;
                     case D3Q27System::B : x1b = x1a         ; x2b = x2a         ; x3b = x3a-q*dx; break;
                     case D3Q27System::BE: x1b = x1a+q*dx; x2b = x2a         ; x3b = x3a-q*dx; break;
                     case D3Q27System::BN: x1b = x1a         ; x2b = x2a+q*dx; x3b = x3a-q*dx; break;
                     case D3Q27System::BW: x1b = x1a-q*dx; x2b = x2a         ; x3b = x3a-q*dx; break;
                     case D3Q27System::BS: x1b = x1a         ; x2b = x2a-q*dx; x3b = x3a-q*dx; break;
                     case D3Q27System::TNE : x1b = x1a+q*dx  ; x2b = x2a+q*dx; x3b = x3a+q*dx; break;
                     case D3Q27System::BSW : x1b = x1a-q*dx  ; x2b = x2a-q*dx; x3b = x3a-q*dx; break;
                     case D3Q27System::BNE : x1b = x1a+q*dx  ; x2b = x2a+q*dx; x3b = x3a-q*dx; break;
                     case D3Q27System::TSW : x1b = x1a-q*dx  ; x2b = x2a-q*dx; x3b = x3a+q*dx; break;
                     case D3Q27System::TSE : x1b = x1a+q*dx  ; x2b = x2a-q*dx; x3b = x3a+q*dx; break;
                     case D3Q27System::BNW : x1b = x1a-q*dx  ; x2b = x2a+q*dx; x3b = x3a-q*dx; break;
                     case D3Q27System::BSE : x1b = x1a+q*dx  ; x2b = x2a-q*dx; x3b = x3a-q*dx; break;
                     case D3Q27System::TNW : x1b = x1a-q*dx  ; x2b = x2a+q*dx; x3b = x3a+q*dx; break;
                     default: throw UbException(UB_EXARGS,"unknown direction");
                     }
                     pointpair.second.setX1(x1b);
                     pointpair.second.setX2(x2b);
                     pointpair.second.setX3(x3b);
                     QsLineSet.push_back(pointpair);
                  }
               }

            }
         }
      }
   }
   //cout<<"getQs: "<<QsLineSet.size()<<endl;
   return QsLineSet;
}

void D3Q27Interactor::writeValidationAVSFile(string filename)
{
   UBLOG(logINFO,"D3Q27Interactor::writeValidationAVSFile("<<filename<<") - start ");
   ofstream out(filename.c_str(),ios::out);
   if(!out) throw UbException(UB_EXARGS,"couldn't open file "+filename);

   int numpoints, numlines;
   vector< pair<GbPoint3D,GbPoint3D> > qsLineSet = this->getQsLineSet();
   numlines  = (unsigned)qsLineSet.size();
   numpoints = numlines*2;

   out<<"# UCD-File created by D3Q27Interactor\n";
   out<<numpoints<<" "<<numlines<<" 0 0 0 "<<endl;
   int nr=1;
   for (int i=0; i<numlines; i++)
   {
      out<<nr++<<" "<<qsLineSet[i].first.getX1Coordinate() <<" " <<qsLineSet[i].first.getX2Coordinate() <<" " <<qsLineSet[i].first.getX3Coordinate() <<" \n";
      out<<nr++<<" "<<qsLineSet[i].second.getX1Coordinate()<<" " <<qsLineSet[i].second.getX2Coordinate()<<" " <<qsLineSet[i].second.getX3Coordinate() <<" \n";
   }
   nr = 1;
   for (int i=0; i<numlines; i++)
   {
      int el = nr+1;
      out<<i+1<<" "<<2<<" line "<<nr<<" "<<el<<" "<<endl;
      nr=el+1;
   }
   UBLOG(logINFO,"D3Q27Interactor::writeValidationAVSFile("<<filename<<") - end");
}
