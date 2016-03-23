#include "D3Q27ETBCProcessor.h"
#include "D3Q27EsoTwist3DSplittedVector.h"

D3Q27ETBCProcessor::D3Q27ETBCProcessor()
{
   
}
//////////////////////////////////////////////////////////////////////////
D3Q27ETBCProcessor::D3Q27ETBCProcessor(LBMKernel3DPtr kernel)
{
   distributions = boost::dynamic_pointer_cast<EsoTwist3D>(kernel->getDataSet()->getFdistributions());
   ghostLayerWitdh = kernel->getGhostLayerWidth();
   collFactor = kernel->getCollisionFactor();
   bcArray.resize( distributions->getNX1(), distributions->getNX2(), distributions->getNX3(), BCArray3D<D3Q27BoundaryCondition>::FLUID);
   this->compressible = kernel->getCompressible();

   minX1 = 0;
   minX2 = 0;
   minX3 = 0;
   maxX1 = (int)bcArray.getNX1();
   maxX2 = (int)bcArray.getNX2();
   maxX3 = (int)bcArray.getNX3();
}
//////////////////////////////////////////////////////////////////////////
D3Q27ETBCProcessor::~D3Q27ETBCProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
//void D3Q27ETBCProcessor::init(LBMKernel3DPtr kernel)
//{
//   distributions = boost::dynamic_pointer_cast<EsoTwist3D>(kernel->getDataSet()->getFdistributions());
//   ghostLayerWitdh = kernel->getGhostLayerWidth();
//   collFactor = kernel->getCollisionFactor();
//   bcArray.resize( distributions->getNX1(), distributions->getNX2(), distributions->getNX3(), BCArray3D<D3Q27BoundaryCondition>::FLUID);
//   this->compressible = kernel->getCompressible();
//   CalcFeqForDirFct calcFeqsForDirFct = NULL;
//   CalcMacrosFct    calcMacrosFct     = NULL;
//   CalcFeqFct       calcFeqFct        = NULL;
//
//   if(this->compressible)
//   {
//      rho0 = REAL_CAST(0.0);
//      calcFeqsForDirFct = &D3Q27System::getCompFeqForDirection; 
//      calcMacrosFct     = &D3Q27System::calcCompMacroscopicValues;
//      calcFeqFct        = &D3Q27System::calcCompFeq; 
//   }
//   else
//   {
//      rho0 = REAL_CAST(1.0);
//      calcFeqsForDirFct = &D3Q27System::getIncompFeqForDirection; 
//      calcMacrosFct     = &D3Q27System::calcIncompMacroscopicValues;
//      calcFeqFct        = &D3Q27System::calcIncompFeq; 
//   }
//
//   minX1 = 0;
//   minX2 = 0;
//   minX3 = 0;
//   maxX1 = (int)bcArray.getNX1();
//   maxX2 = (int)bcArray.getNX2();
//   maxX3 = (int)bcArray.getNX3();
//}
//////////////////////////////////////////////////////////////////////////
BCProcessorPtr D3Q27ETBCProcessor::clone(LBMKernel3DPtr kernel)
{
   BCProcessorPtr bcProcessor(new D3Q27ETBCProcessor(kernel));
   //boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(bcProcessor)->init(kernel);
   BoundaryConditionPtr velocity = this->getBC(BoundaryCondition::Velocity);
   BoundaryConditionPtr density  = this->getBC(BoundaryCondition::Density);
   BoundaryConditionPtr noSlip   = this->getBC(BoundaryCondition::NoSlip);
   BoundaryConditionPtr slip     = this->getBC(BoundaryCondition::Slip);
   if (velocity)bcProcessor->addBC(velocity->clone());
   if (density) bcProcessor->addBC(density->clone());
   if (noSlip)  bcProcessor->addBC(noSlip->clone());
   if (slip)    bcProcessor->addBC(slip->clone());
   return bcProcessor;
}
//////////////////////////////////////////////////////////////////////////
BCArray3D<D3Q27BoundaryCondition>& D3Q27ETBCProcessor::getBCArray()
{ 
   return this->bcArray; 
}
//////////////////////////////////////////////////////////////////////////
void D3Q27ETBCProcessor::applyBC()
{
   typedef void      (*CalcMacrosFct)(const LBMReal* const& /*f[27]*/, LBMReal& /*rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);
   typedef LBMReal (*CalcFeqForDirFct)(const int& /*direction*/,const LBMReal& /*(d)rho*/,const LBMReal& /*vx1*/,const LBMReal& /*vx2*/,const LBMReal& /*vx3*/);
   typedef  void     (*CalcFeqFct)(LBMReal* const& /*feq/*[27]*/,const LBMReal& /*rho*/,const LBMReal& /*vx1*/,const LBMReal& /*vx2*/,const LBMReal& /*vx3*/);

   CalcFeqForDirFct calcFeqsForDirFct = NULL;
   CalcMacrosFct    calcMacrosFct     = NULL;
   CalcFeqFct       calcFeqFct        = NULL;

   if(this->compressible)
   {
      rho0 = REAL_CAST(0.0);
      calcFeqsForDirFct = &D3Q27System::getCompFeqForDirection; 
      calcMacrosFct     = &D3Q27System::calcCompMacroscopicValues;
      calcFeqFct        = &D3Q27System::calcCompFeq; 
   }
   else
   {
      rho0 = REAL_CAST(1.0);
      calcFeqsForDirFct = &D3Q27System::getIncompFeqForDirection; 
      calcMacrosFct     = &D3Q27System::calcIncompMacroscopicValues;
      calcFeqFct        = &D3Q27System::calcIncompFeq; 
   }



   //////////////////////////////////////////////////////////////////////////
   //EsoTwist3DPtr distributionsTemp(new D3Q27EsoTwist3DSplittedVector(maxX1, maxX2, maxX2, -999.0));
   //printf("create new distributionsTemp\n");
   //////////////////////////////////////////////////////////////////////////

   for(int x3 = minX3; x3 < maxX3; x3++)      
   {
      for(int x2 = minX2; x2 < maxX2; x2++)      
      {
         for(int x1 = minX1; x1 < maxX1; x1++)  
         {
            if(!bcArray.isSolid(x1,x2,x3) && !bcArray.isUndefined(x1,x2,x3)) 
            {
               if( (bcPtr=bcArray.getBC(x1,x2,x3)) != NULL)
               {
                  //velocity boundary condition
                  if(bcPtr->hasVelocityBoundary())
                  {
                     //distributions->getDistributionInv(f,x1,x2,x3);
                     //int nx1 = x1;
                     //int nx2 = x2;
                     //int nx3 = x3;

                     //rho = rho0;
                     //if(compressible) rho += D3Q27System::getDensity(f); 

                     //for(int fdir=D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
                     //{
                     //   if(bcPtr->hasVelocityBoundaryFlag(fdir))
                     //   {
                     //      int invDir = D3Q27System::INVDIR[fdir];
                     //      LBMReal ftemp = distributions->getDistributionInvForDirection(x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
                     //      ftemp += rho * bcPtr->getBoundaryVelocity(fdir);
                     //      distributions->setDistributionForDirection(ftemp, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
                     //   }
                     //}
                     distributions->getDistributionInv(ftemp,x1,x2,x3);
                     calcMacrosFct(ftemp,rho,vx1,vx2,vx3);
                     calcFeqFct(feq,rho,vx1,vx2,vx3);

                     for(int fdir=D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
                     {
                        if(bcPtr->hasVelocityBoundaryFlag(fdir))
                        {
                           const int invDir = D3Q27System::INVDIR[fdir];
                           LBMReal q = 1.0;//bcPtr->getQ(invDir);// m+m q=0 stabiler
                           LBMReal velocity = bcPtr->getBoundaryVelocity(invDir);
                           //LBMReal fReturn = ((1.0-q)/(1.0+q))*((ftemp[invDir]-feq[invDir]*collFactor)/(1.0-collFactor))+((q*(ftemp[invDir]+ftemp[fdir])-velocity)/(1.0+q));
                           LBMReal fReturn = ((1.0-q)/(1.0+q))*((ftemp[invDir]-feq[invDir])/(1.0-collFactor)+feq[invDir])+((q*(ftemp[invDir]+ftemp[fdir])-velocity)/(1.0+q));
                           distributions->setDistributionForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
                        }
                     }
                  }
///////////////////velocity bc for non reflecting pressure bc
                  //if (bcPtr->hasVelocityBoundary())
                  //{
                  //   distributions->getDistributionInv(ftemp, x1, x2, x3);
                  //   calcMacrosFct(ftemp, rho, vx1, vx2, vx3);
                  //   calcFeqFct(feq, rho, vx1, vx2, vx3);

                  //   for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++)
                  //   {
                  //      if (bcPtr->hasVelocityBoundaryFlag(fdir))
                  //      {
                  //         const int invDir = D3Q27System::INVDIR[fdir];
                  //         LBMReal q = 1.0;//bcPtr->getQ(invDir);// m+m q=0 stabiler
                  //         LBMReal velocity = bcPtr->getBoundaryVelocity(invDir);
                  //         //LBMReal fReturn = ((1.0-q)/(1.0+q))*((ftemp[invDir]-feq[invDir]*collFactor)/(1.0-collFactor))+((q*(ftemp[invDir]+ftemp[fdir])-velocity)/(1.0+q));
                  //         LBMReal fReturn = ((1.0 - q) / (1.0 + q))*((ftemp[invDir] - feq[invDir]) / (1.0 - collFactor) + feq[invDir]) + ((q*(ftemp[invDir] + ftemp[fdir]) - velocity) / (1.0 + q))-rho*D3Q27System::WEIGTH[invDir];
                  //         distributions->setDistributionForDirection(fReturn, x1 + D3Q27System::DX1[invDir], x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], fdir);
                  //      }
                  //   }
                  //}
                  //density boundary condition
                  else if(bcPtr->hasDensityBoundary())
                  {
                     distributions->getDistributionInv(f,x1,x2,x3);
                     int nx1 = x1;
                     int nx2 = x2;
                     int nx3 = x3;
                     int direction = -1;

                     //flag points in direction of fluid
                     if     ( bcPtr->hasDensityBoundaryFlag(D3Q27System::E) ) {nx1 -= 1; direction=D3Q27System::E; }
                     else if( bcPtr->hasDensityBoundaryFlag(D3Q27System::W) ) {nx1 += 1; direction=D3Q27System::W; }
                     else if( bcPtr->hasDensityBoundaryFlag(D3Q27System::N) ) {nx2 -= 1; direction=D3Q27System::N; }
                     else if( bcPtr->hasDensityBoundaryFlag(D3Q27System::S) ) {nx2 += 1; direction=D3Q27System::S; }
                     else if( bcPtr->hasDensityBoundaryFlag(D3Q27System::T) ) {nx3 -= 1; direction=D3Q27System::T; }
                     else if( bcPtr->hasDensityBoundaryFlag(D3Q27System::B) ) {nx3 += 1; direction=D3Q27System::B; }
                     else UB_THROW( UbException(UB_EXARGS,"Danger...no orthogonal BC-Flag on density boundary...") );
                        
                     #ifdef _DEBUG
                     if( nx1<0 || nx1>maxX1 ) UB_THROW( UbException(UB_EXARGS,"nx1<0 || nx1>=lengthX1") );
                     if( nx2<0 || nx2>maxX2 ) UB_THROW( UbException(UB_EXARGS,"nx2<0 || nx2>=lengthX2") );
                     if( nx3<0 || nx3>maxX3 ) UB_THROW( UbException(UB_EXARGS,"nx3<0 || nx3>=lengthX3") );
                     #endif

                     //if(bcPtr->getDensitySecondaryOption(direction)!=0) 
                     //   UB_THROW( UbException(UB_EXARGS,"invalid bc->getDensitySecondaryOption()="+UbSystem::toString(bcPtr->getDensitySecondaryOption(direction))) );

                     calcMacrosFct(f,rho,vx1,vx2,vx3);
                     LBMReal rhoBC = bcPtr->getBoundaryDensity();
                     for(int fdir=D3Q27System::STARTF; fdir<=D3Q27System::ENDF; fdir++)
                      {
                         if(bcPtr->hasDensityBoundaryFlag(fdir))
                         {
                            short densitySO = bcPtr->getDensitySecondaryOption(fdir);
                            if ( densitySO == 0)
                            {
                               // Martins NEQ ADDON
                               ////original: 15.2.2013:
                               LBMReal ftemp = calcFeqsForDirFct(fdir,rho,vx1,vx2,vx3);             
                               ftemp = calcFeqsForDirFct(fdir,rhoBC,vx1,vx2,vx3)+f[fdir]-ftemp;
                               distributions->setDistributionForDirection(ftemp, nx1, nx2, nx3, fdir);
                            }
                            else if(densitySO == 1)
                            {
                               //Ehsan: 15.2.2013:
                               LBMReal ftemp = calcFeqsForDirFct(fdir,rhoBC,vx1,vx2,vx3);
                               distributions->setDistributionForDirection(ftemp, nx1, nx2, nx3, fdir);
                            }
                         }
                     }
                  }
///////////////non reflecting pressure baundary condition             
               //else if(bcPtr->hasDensityBoundary())
               //   {
//                     distributions->swap();
//                     distributions->getDistribution(f, x1, x2, x3);
//                     int nx1 = x1;
//                     int nx2 = x2;
//                     int nx3 = x3;
//                     int direction = -1;
//
//                     //flag points in direction of fluid
//                     if (bcPtr->hasDensityBoundaryFlag(D3Q27System::E))      { nx1 += 1; direction = D3Q27System::E; }
//                     else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::W)) { nx1 -= 1; direction = D3Q27System::W; }
//                     else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::N)) { nx2 += 1; direction = D3Q27System::N; }
//                     else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::S)) { nx2 -= 1; direction = D3Q27System::S; }
//                     else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::T)) { nx3 += 1; direction = D3Q27System::T; }
//                     else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::B)) { nx3 -= 1; direction = D3Q27System::B; }
//                     else UB_THROW(UbException(UB_EXARGS, "Danger...no orthogonal BC-Flag on density boundary..."));
//
//#ifdef _DEBUG
//                     if (nx1<0 || nx1>maxX1) UB_THROW(UbException(UB_EXARGS, "nx1<0 || nx1>=lengthX1"));
//                     if (nx2<0 || nx2>maxX2) UB_THROW(UbException(UB_EXARGS, "nx2<0 || nx2>=lengthX2"));
//                     if (nx3<0 || nx3>maxX3) UB_THROW(UbException(UB_EXARGS, "nx3<0 || nx3>=lengthX3"));
//#endif
//                     distributions->getDistribution(ftemp, nx1, nx2, nx3);
//                     calcMacrosFct(f, rho, vx1, vx2, vx3);
//
//                     distributions->swap();
//
//                     double dim = 0.01;
//
//                     f[D3Q27System::W] = ftemp[D3Q27System::W] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[D3Q27System::W] - rho*dim*D3Q27System::WEIGTH[D3Q27System::W];
//                     f[D3Q27System::NW] = ftemp[D3Q27System::NW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[D3Q27System::NW] - rho*dim*D3Q27System::WEIGTH[D3Q27System::NW];
//                     f[D3Q27System::SW] = ftemp[D3Q27System::SW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[D3Q27System::SW] - rho*dim*D3Q27System::WEIGTH[D3Q27System::SW];
//                     f[D3Q27System::TW] = ftemp[D3Q27System::TW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[D3Q27System::TW] - rho*dim*D3Q27System::WEIGTH[D3Q27System::TW];
//                     f[D3Q27System::BW] = ftemp[D3Q27System::BW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[D3Q27System::BW] - rho*dim*D3Q27System::WEIGTH[D3Q27System::BW];
//                     f[D3Q27System::TNW] = ftemp[D3Q27System::TNW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[D3Q27System::TNW] - rho*dim*D3Q27System::WEIGTH[D3Q27System::TNW];
//                     f[D3Q27System::TSW] = ftemp[D3Q27System::TSW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[D3Q27System::TSW] - rho*dim*D3Q27System::WEIGTH[D3Q27System::TSW];
//                     f[D3Q27System::BNW] = ftemp[D3Q27System::BNW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[D3Q27System::BNW] - rho*dim*D3Q27System::WEIGTH[D3Q27System::BNW];
//                     f[D3Q27System::BSW] = ftemp[D3Q27System::BSW] * (one_over_sqrt3 - vx1) + (1.0 - one_over_sqrt3 + vx1)*f[D3Q27System::BSW] - rho*dim*D3Q27System::WEIGTH[D3Q27System::BSW];
//                                                                                                                                             
//                     distributions->setDistributionInvForDirection(f[D3Q27System::W], x1, x2, x3, D3Q27System::W);
//                     distributions->setDistributionInvForDirection(f[D3Q27System::NW], x1, x2, x3, D3Q27System::NW);
//                     distributions->setDistributionInvForDirection(f[D3Q27System::SW], x1, x2, x3, D3Q27System::SW);
//                     distributions->setDistributionInvForDirection(f[D3Q27System::TW], x1, x2, x3, D3Q27System::TW);
//                     distributions->setDistributionInvForDirection(f[D3Q27System::BW], x1, x2, x3, D3Q27System::BW);
//                     distributions->setDistributionInvForDirection(f[D3Q27System::TNW], x1, x2, x3, D3Q27System::TNW);
//                     distributions->setDistributionInvForDirection(f[D3Q27System::TSW], x1, x2, x3, D3Q27System::TSW);
//                     distributions->setDistributionInvForDirection(f[D3Q27System::BNW], x1, x2, x3, D3Q27System::BNW);
//                     distributions->setDistributionInvForDirection(f[D3Q27System::BSW], x1, x2, x3, D3Q27System::BSW);
                   
                  //}
               else if(bcPtr->hasNoSlipBoundary())
               {
                  distributions->getDistributionInv(ftemp,x1,x2,x3);
                  calcMacrosFct(ftemp,rho,vx1,vx2,vx3);
                  calcFeqFct(feq,rho,vx1,vx2,vx3);

                  for(int fdir=D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
                  {
                     if(bcPtr->hasNoSlipBoundaryFlag(fdir))
                     {
                        switch(bcPtr->getNoSlipSecondaryOption(fdir))
                        {
                        case 0:
                           //simple bounce back - do nothing
                           break;
                        case 1:
                           {
                              //quadratic bounce back
                              const int invDir = D3Q27System::INVDIR[fdir];
                              LBMReal q = bcPtr->getQ(invDir);
                              //LBMReal fReturn = ((1.0-q)/(1.0+q))*((ftemp[invDir]-feq[invDir]*collFactor)/(1.0-collFactor))+((q/(1.0+q))*(ftemp[invDir]+ftemp[fdir]));
                              LBMReal fReturn = ((1.0-q)/(1.0+q))*((ftemp[invDir]-feq[invDir])/(1.0-collFactor)+feq[invDir])+((q/(1.0+q))*(ftemp[invDir]+ftemp[fdir]));
                              distributions->setDistributionForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
                           }
                           break;
                           //////////////////////////////////////////////////////////////////////////
                        ////HighViscosityNoSlipBoundaryCondition
                        //case 1:
                        //{
                        //   //quadratic bounce back
                        //   const int invDir = D3Q27System::INVDIR[fdir];
                        //   LBMReal q = bcPtr->getQ(invDir);
                        //   LBMReal fReturn = ((1.0-q)*ftemp[invDir]+q*((ftemp[invDir]+ftemp[fdir])*(1-collFactor)+collFactor*(feq[invDir]+feq[fdir])))/(1.0+q);
                        //   distributions->setDistributionForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
                        //}
                        //break;
                        //////////////////////////////////////////////////////////////////////////
                        //case 2:
                        //   {
                        //      //quadratic bounce back 2nd choice
                        //      const int invDir = D3Q27System::INVDIR[fdir];
                        //      LBMReal q = bcPtr->getQ(invDir);
                        //      LBMReal fReturn = ((1.0-q)/(1.0+q))*0.5*(ftemp[invDir]-ftemp[fdir]+(ftemp[invDir]+ftemp[fdir]-collFactor*(feq[fdir]+feq[invDir]))/(1.0-collFactor)); 
                        //      distributions->setDistributionForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
                        //   }
                        //   break;
                        case 2:
                           {
                              //quadratic bounce back with for thin walls
                              const int invDir = D3Q27System::INVDIR[fdir];
                              LBMReal q = bcPtr->getQ(invDir);
                              //LBMReal fReturn = ((1.0-q)/(1.0+q))*((ftemp[invDir]-feq[invDir])/(1.0-collFactor)+feq[invDir])+((q/(1.0+q))*(ftemp[invDir]+ftemp[fdir]));
                              LBMReal fReturn = ((1.0-q)/(1.0+q))*0.5*(ftemp[invDir]-ftemp[fdir]+(ftemp[invDir]+ftemp[fdir]-collFactor*(feq[fdir]+feq[invDir]))/(1.0-collFactor)); 
                              //printf("set distributionsTemp - start\n");
                              distributionsTemp->setDistributionForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
                              //printf("set distributionsTemp - end\n");
                           }
                           break;
                        }
                     }
                  }
               }
               else if(bcPtr->hasSlipBoundary()) 
               {
                  distributions->getDistributionInv(ftemp,x1,x2,x3);
                  calcMacrosFct(ftemp,rho,vx1,vx2,vx3);
                  calcFeqFct(feq,rho,vx1,vx2,vx3);

                  UbTupleFloat3 normale = bcPtr->getNormalVector();
                  LBMReal amp = vx1*val<1>(normale)+vx2*val<2>(normale)+vx3*val<3>(normale);

                  vx1 = vx1 - amp * val<1>(normale); //normale zeigt von struktur weg!
                  vx2 = vx2 - amp * val<2>(normale); //normale zeigt von struktur weg!
                  vx3 = vx3 - amp * val<3>(normale); //normale zeigt von struktur weg!

                  for(int fdir=D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
                  {
                     if(bcPtr->hasSlipBoundaryFlag(fdir))
                     {


                        switch( bcPtr->getSlipSecondaryOption(fdir) )
                        {
                        case 0:
                           //doesn't work
                           //simple bounce back - do nothing
                           break;
                        case 1:
                           {
                              //quadratic bounce back
                              const int invDir = D3Q27System::INVDIR[fdir];
                              LBMReal q = 1.0;//bcPtr->getQ(invDir);// m+m q=0 stabiler
                              //vx3=0;
                              LBMReal velocity = 0.0;
                              switch(invDir)
                              {
                              case D3Q27System::E : velocity = ( UbMath::c4o9*(+vx1) ); break;      //(2/cs^2)(=6)*rho_0(=1 bei imkompr)*wi*u*ei mit cs=1/sqrt(3)
                              case D3Q27System::W : velocity = ( UbMath::c4o9*(-vx1) ); break;      //z.B. aus paper manfred MRT LB models in three dimensions (2002)   
                              case D3Q27System::N : velocity = ( UbMath::c4o9*(+vx2) ); break;   
                              case D3Q27System::S : velocity = ( UbMath::c4o9*(-vx2) ); break;
                              case D3Q27System::T : velocity = ( UbMath::c4o9*(+vx3) ); break;
                              case D3Q27System::B : velocity = ( UbMath::c4o9*(-vx3) ); break;
                              case D3Q27System::NE: velocity = ( UbMath::c1o9*(+vx1+vx2             ) ); break;
                              case D3Q27System::SW: velocity = ( UbMath::c1o9*(-vx1-vx2             ) ); break;
                              case D3Q27System::SE: velocity = ( UbMath::c1o9*(+vx1-vx2             ) ); break;
                              case D3Q27System::NW: velocity = ( UbMath::c1o9*(-vx1+vx2             ) ); break;
                              case D3Q27System::TE: velocity = ( UbMath::c1o9*(+vx1             +vx3) ); break;
                              case D3Q27System::BW: velocity = ( UbMath::c1o9*(-vx1             -vx3) ); break;
                              case D3Q27System::BE: velocity = ( UbMath::c1o9*(+vx1             -vx3) ); break;
                              case D3Q27System::TW: velocity = ( UbMath::c1o9*(-vx1             +vx3) ); break;
                              case D3Q27System::TN: velocity = ( UbMath::c1o9*(             +vx2+vx3) ); break;
                              case D3Q27System::BS: velocity = ( UbMath::c1o9*(             -vx2-vx3) ); break;
                              case D3Q27System::BN: velocity = ( UbMath::c1o9*(             +vx2-vx3) ); break;
                              case D3Q27System::TS: velocity = ( UbMath::c1o9*(             -vx2+vx3) ); break;
                              case D3Q27System::TNE: velocity = ( UbMath::c1o36*(+vx1+vx2+vx3) ); break;
                              case D3Q27System::BSW: velocity = ( UbMath::c1o36*(-vx1-vx2-vx3) ); break;
                              case D3Q27System::BNE: velocity = ( UbMath::c1o36*(+vx1+vx2-vx3) ); break;
                              case D3Q27System::TSW: velocity = ( UbMath::c1o36*(-vx1-vx2+vx3) ); break;
                              case D3Q27System::TSE: velocity = ( UbMath::c1o36*(+vx1-vx2+vx3) ); break;
                              case D3Q27System::BNW: velocity = ( UbMath::c1o36*(-vx1+vx2-vx3) ); break;
                              case D3Q27System::BSE: velocity = ( UbMath::c1o36*(+vx1-vx2-vx3) ); break;
                              case D3Q27System::TNW: velocity = ( UbMath::c1o36*(-vx1+vx2+vx3) ); break; 
                              default: throw UbException(UB_EXARGS,"unknown error");
                              }
                              //LBMReal fReturn = ((1.0-q)/(1.0+q))*((ftemp[invDir]-feq[invDir]*collFactor)/(1.0-collFactor))+((q*(ftemp[invDir]+ftemp[fdir])-velocity)/(1.0+q));
                              LBMReal fReturn = ((1.0-q)/(1.0+q))*((ftemp[invDir]-feq[invDir])/(1.0-collFactor)+feq[invDir])+((q*(ftemp[invDir]+ftemp[fdir])-velocity)/(1.0+q));
                              distributions->setDistributionForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
                           }
                           break;
                        default:  
                           UB_THROW( UbException(UB_EXARGS,"unknown secondary option for SlipBC") );
                        }
                     }
                  }
               }
               }
            }
         }   
      }
   }


//   for(int x3 = minX3; x3 < maxX3; x3++)      
//   {
//      for(int x2 = minX2; x2 < maxX2; x2++)      
//      {
//         for(int x1 = minX1; x1 < maxX1; x1++)  
//         {
//            if(!bcArray.isSolid(x1,x2,x3) && !bcArray.isUndefined(x1,x2,x3)) 
//            {
//               if( (bcPtr=bcArray.getBC(x1,x2,x3)) != NULL)
//               {
//                  if(bcPtr->hasNoSlipBoundary())
//                  {
//                     //distributions->getDistributionInv(ftemp,x1,x2,x3);
//                     //calcMacrosFct(ftemp,rho,vx1,vx2,vx3);
//                     //calcFeqFct(feq,rho,vx1,vx2,vx3);
//
//                     for(int fdir=D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
//                     {
//                        if(bcPtr->hasNoSlipBoundaryFlag(fdir))
//                        {
//                           switch(bcPtr->getNoSlipSecondaryOption(fdir))
//                           {
//                           case 2:
//                              {
//                                 //quadratic bounce back with for thin walls
//                                const int invDir = D3Q27System::INVDIR[fdir];
///*                                 LBMReal q = bcPtr->getQ(invDir);
//                                 LBMReal fReturn = ((1.0-q)/(1.0+q))*0.5*(ftemp[invDir]-ftemp[fdir]+(ftemp[invDir]+ftemp[fdir]-collFactor*(feq[fdir]+feq[invDir]))/(1.0-collFactor));*/ 
//                                //printf("get distributionsTemp - start\n");
//                                //printf("x1=%d,x2=%d,x3=%d,D3Q27System::DX1[invDir]=%d,D3Q27System::DX2[invDir]=%d,D3Q27System::DX3[invDir]%d\n",x1,x2,x3,D3Q27System::DX1[invDir],D3Q27System::DX2[invDir],D3Q27System::DX3[invDir]);
//                                ftemp[fdir] = distributionsTemp->getDistributionInvForDirection(x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
//                                //printf("get distributionsTemp - end\n"); 
//                                distributions->setDistributionForDirection(ftemp[fdir], x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
//                              }
//                              break;
//                           }
//                        }
//                     }
//                  }
//               }
//            }
//         }
//      }
//   }

}
//////////////////////////////////////////////////////////////////////////
//void D3Q27ETBCProcessor::init()
//{
//   BoundaryConditionPtr velocity = getBC(BoundaryCondition::Velocity);
//   BoundaryConditionPtr density  = getBC(BoundaryCondition::Density);
//   BoundaryConditionPtr noSlip   = getBC(BoundaryCondition::NoSlip);
//   BoundaryConditionPtr slip     = getBC(BoundaryCondition::Slip);
//
//   if (velocity) velocity->setCompressible(compressible);
//   if (density) density->setCompressible(compressible);
//   if (noSlip) noSlip->setCompressible(compressible);
//   if (slip) slip->setCompressible(compressible);
//
//   int vCount = 0;
//   int dCount = 0;
//   int nsCount = 0;
//   int sCount = 0;
//
//   for (int x3 = minX3; x3 < maxX3; x3++)
//   {
//      for (int x2 = minX2; x2 < maxX2; x2++)
//      {
//         for (int x1 = minX1; x1 < maxX1; x1++)
//         {
//            if (!bcArray.isSolid(x1, x2, x3) && !bcArray.isUndefined(x1, x2, x3))
//            {
//               if ((bcPtr = bcArray.getBC(x1, x2, x3)) != NULL)
//               {
//                  //velocity boundary condition
//                  if (bcPtr->hasVelocityBoundary() && velocity)
//                  {
//                     velocity->addNode(x1, x2, x3);
//                     velocity->addBcPointer(bcPtr);
//                     vCount++;
//                  }
//                  //density boundary condition
//                  else if (bcPtr->hasDensityBoundary() && density)
//                  {
//                     density->addNode(x1, x2, x3);
//                     density->addBcPointer(bcPtr);
//                     dCount++;
//                  }
//                  //no slip boundary condition
//                  else if (bcPtr->hasNoSlipBoundary() && noSlip)
//                  {
//                     noSlip->addNode(x1, x2, x3);
//                     noSlip->addBcPointer(bcPtr);
//                     nsCount++;
//                  }
//                  //slip boundary condition
//                  else if (bcPtr->hasSlipBoundary() && slip)
//                  {
//                     slip->addNode(x1, x2, x3);
//                     slip->addBcPointer(bcPtr);
//                     sCount++;
//                  }
//               }
//            }
//         }
//      }
//   }
//   if (vCount > 0)
//   {
//      velocity->addNnode(vCount);
//      velocity->addDistributions(distributions);
//      velocity->setCollFactor(collFactor);
//   }
//   if (dCount > 0)
//   {
//      density->addNnode(dCount);
//      density->addDistributions(distributions);
//      density->setCollFactor(collFactor);
//   }
//   if (nsCount > 0)
//   {
//      noSlip->addNnode(nsCount);
//      noSlip->addDistributions(distributions);
//      noSlip->setCollFactor(collFactor);
//   }
//   if (sCount > 0)
//   {
//      slip->addNnode(sCount);
//      slip->addDistributions(distributions);
//      slip->setCollFactor(collFactor);
//   }
//}
//////////////////////////////////////////////////////////////////////////
void D3Q27ETBCProcessor::addBC(BoundaryConditionPtr bc)
{
   BoundaryCondition::Type type = bc->getType();
   if (bc->isPreCollision())
   {
      preBC.push_back(bc);
   }
   else
   {
      postBC.push_back(bc);
   }
}
//////////////////////////////////////////////////////////////////////////
BoundaryConditionPtr D3Q27ETBCProcessor::getBC(BoundaryCondition::Type type)
{
   BOOST_FOREACH(BoundaryConditionPtr bc, preBC)
   {
      if (bc->getType() == type)
      {
         return bc;
      }
   }

   BOOST_FOREACH(BoundaryConditionPtr bc, postBC)
   {
      if (bc->getType() == type)
      {
         return bc;
      }
   }

   return BoundaryConditionPtr();
}
//////////////////////////////////////////////////////////////////////////
void D3Q27ETBCProcessor::applyPreCollisionBC()
{
   BOOST_FOREACH(BoundaryConditionPtr bc, preBC)
   {
      bc->apply();
   }
}
//////////////////////////////////////////////////////////////////////////
void D3Q27ETBCProcessor::applyPostCollisionBC()
{
   BOOST_FOREACH(BoundaryConditionPtr bc, postBC)
   {
      bc->apply();
   }
}
