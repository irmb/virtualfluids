#include "LBMKernelETD3Q27Cascaded.h"
#include "D3Q27System.h"
#include "D3Q27NoSlipBCAdapter.h"
#include "D3Q27DensityBCAdapter.h"
#include "D3Q27VelocityBCAdapter.h"
#include "SimulationParameters.h"
#include "D3Q27InterpolationProcessor.h"
#include "D3Q27EsoTwist3DSplittedVector.h"


//BOOST_CLASS_EXPORT(LBMKernelETD3Q27Cascaded)

//////////////////////////////////////////////////////////////////////////
LBMKernelETD3Q27Cascaded::LBMKernelETD3Q27Cascaded()
{

}
//////////////////////////////////////////////////////////////////////////
LBMKernelETD3Q27Cascaded::LBMKernelETD3Q27Cascaded(int nx1, int nx2, int nx3) 
   : LBMKernelETD3Q27(nx1, nx2, nx3)
{
   this->compressible = true;
}
//////////////////////////////////////////////////////////////////////////
LBMKernelETD3Q27Cascaded::~LBMKernelETD3Q27Cascaded(void)
{

}
//////////////////////////////////////////////////////////////////////////
void LBMKernelETD3Q27Cascaded::init()
{
   //DistributionArray3DPtr d(new D3Q27EsoTwist3DSplittedVector(nx1+ghostLayerWitdh*2, nx2+ghostLayerWitdh*2, nx3+ghostLayerWitdh*2, -999.0));
   DistributionArray3DPtr d(new D3Q27EsoTwist3DSplittedVector(nx1+2, nx2+2, nx3+2, -999.0));
   dataSet->setFdistributions(d);
}
//////////////////////////////////////////////////////////////////////////
LBMKernel3DPtr LBMKernelETD3Q27Cascaded::clone()
{
   LBMKernel3DPtr kernel(new LBMKernelETD3Q27Cascaded(nx1, nx2, nx3));
   boost::dynamic_pointer_cast<LBMKernelETD3Q27Cascaded>(kernel)->init();
   kernel->setCollisionFactor(this->collFactor);
   kernel->setBCProcessor(bcProcessor->clone(kernel));
   kernel->setWithForcing(withForcing);
   kernel->setForcingX1(muForcingX1);
   kernel->setForcingX2(muForcingX2);
   kernel->setForcingX3(muForcingX3);
   kernel->setIndex(ix1, ix2, ix3);
   kernel->setDeltaT(deltaT);
   return kernel;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernelETD3Q27Cascaded::calculate()
{
   timer.resetAndStart();
   collideAll();
   timer.stop();
}
//////////////////////////////////////////////////////////////////////////
void LBMKernelETD3Q27Cascaded::collideAll()
{
   using namespace D3Q27System;

   //initializing of forcing stuff 
   if (withForcing)
   {
      muForcingX1.DefineVar("x1",&muX1); muForcingX1.DefineVar("x2",&muX2); muForcingX1.DefineVar("x3",&muX3);
      muForcingX2.DefineVar("x1",&muX1); muForcingX2.DefineVar("x2",&muX2); muForcingX2.DefineVar("x3",&muX3);
      muForcingX3.DefineVar("x1",&muX1); muForcingX3.DefineVar("x2",&muX2); muForcingX3.DefineVar("x3",&muX3);

      muDeltaT = deltaT;

      muForcingX1.DefineVar("dx",&muDeltaT);
      muForcingX2.DefineVar("dx",&muDeltaT);
      muForcingX3.DefineVar("dx",&muDeltaT);

      muNue = (1.0/3.0)*(1.0/collFactor - 1.0/2.0);

      muForcingX1.DefineVar("nue",&muNue);
      muForcingX2.DefineVar("nue",&muNue);
      muForcingX3.DefineVar("nue",&muNue);

      LBMReal forcingX1 = 0;
      LBMReal forcingX2 = 0;
      LBMReal forcingX3 = 0;
   }
   /////////////////////////////////////

   s9 = - collFactor;
   c1o27=1.0/27.0;
   c2o3=2.0/3.0;
   w2=-1.0; //MXXpMYYpMZZ bulk viscosity
   w7=-1.0;//s9; //ORDER 4 Isotropic
   w9=-1.0;
   w10=-1.0;//s9;//-1.0; // ORDER 6 Isotropic
   w1=s9;
   // wenn es mal an den Ecken nicht gut aussieht -2.0-s9 probieren
   w3=-1.0;//-2.0-s9;//-1.0;//MXXYpMYZZ
   w4=-1.0;//-2.0-s9;//-1.0;//MXXYmMYZZ
   w5=-1.0;//-2.0-s9;//-1.0;//MYXZ
   w6=-1.0; //MXXYYpm2p
   w8=-1.0; //M_zXXYZ 

   localDistributions = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getLocalDistributions();
   nonLocalDistributions = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getNonLocalDistributions();
   zeroDistributions = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getZeroDistributions();

   BCArray3D<D3Q27BoundaryCondition>& bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(this->getBCProcessor())->getBCArray();

   const int bcArrayMaxX1 = (int)bcArray.getNX1();
   const int bcArrayMaxX2 = (int)bcArray.getNX2();
   const int bcArrayMaxX3 = (int)bcArray.getNX3();
   
   int minX1 = ghostLayerWidth;
   int minX2 = ghostLayerWidth;
   int minX3 = ghostLayerWidth;
   int maxX1 = bcArrayMaxX1-ghostLayerWidth;
   int maxX2 = bcArrayMaxX2-ghostLayerWidth;
   int maxX3 = bcArrayMaxX3-ghostLayerWidth;

   for(int x3 = minX3; x3 < maxX3; x3++)
   {
      for(int x2 = minX2; x2 < maxX2; x2++)
      {
         for(int x1 = minX1; x1 < maxX1; x1++)
         {
            if(!bcArray.isSolid(x1,x2,x3) && !bcArray.isUndefined(x1,x2,x3))
            {
               //////////////////////////////////////////////////////////////////////////
               //read distribution
               ////////////////////////////////////////////////////////////////////////////
               f[ZERO] = (*this->zeroDistributions)(x1,x2,x3);

               f[E] = (*this->localDistributions)(D3Q27System::ET_E, x1,x2,x3);
               f[N] = (*this->localDistributions)(D3Q27System::ET_N,x1,x2,x3);  
               f[T] = (*this->localDistributions)(D3Q27System::ET_T,x1,x2,x3);
               f[NE] = (*this->localDistributions)(D3Q27System::ET_NE,x1,x2,x3);
               f[NW] = (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,x3);
               f[TE] = (*this->localDistributions)(D3Q27System::ET_TE,x1,x2,x3);
               f[TW] = (*this->localDistributions)(D3Q27System::ET_TW, x1+1,x2,x3);
               f[TN] = (*this->localDistributions)(D3Q27System::ET_TN,x1,x2,x3);
               f[TS] = (*this->localDistributions)(D3Q27System::ET_TS,x1,x2+1,x3);
               f[TNE] = (*this->localDistributions)(D3Q27System::ET_TNE,x1,x2,x3);
               f[TNW] = (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,x3);
               f[TSE] = (*this->localDistributions)(D3Q27System::ET_TSE,x1,x2+1,x3);
               f[TSW] = (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3);

               f[W ] = (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,x3  );
               f[S ] = (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,x2+1,x3  );
               f[B ] = (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,x2,x3+1  );
               f[SW] = (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3 );
               f[SE] = (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,x2+1,x3 );
               f[BW] = (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,x3+1 );
               f[BE] = (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,x2,x3+1 );
               f[BS] = (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,x2+1,x3+1 );
               f[BN] = (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,x2,x3+1 );
               f[BSW] = (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1);
               f[BSE] = (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,x2+1,x3+1);
               f[BNW] = (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,x3+1);
               f[BNE] = (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,x2,x3+1);
               //////////////////////////////////////////////////////////////////////////


               //////////////////////////////////////////////////////////////////////////
               //compute collision
               //////////////////////////////////////////////////////////////////////////
               rho = f[ZERO] + f[E] + f[W] + f[N] + f[S] + f[T] + f[B] 
               + f[NE] + f[SW] + f[SE] + f[NW] + f[TE] + f[BW] + f[BE]
               + f[TW] + f[TN] + f[BS] + f[BN] + f[TS] + f[TNE] + f[TSW]
               + f[TSE] + f[TNW] + f[BNE] + f[BSW] + f[BSE] + f[BNW];

               vx = f[E] - f[W] + f[NE] - f[SW] + f[SE] - f[NW] + f[TE] - f[BW]
               + f[BE] - f[TW] + f[TNE] - f[TSW] + f[TSE] - f[TNW] + f[BNE] - f[BSW]
               + f[BSE] - f[BNW]; 

               vy = f[N] - f[S] + f[NE] - f[SW] - f[SE] + f[NW] + f[TN] - f[BS] + f[BN]
               - f[TS] + f[TNE] - f[TSW] - f[TSE] + f[TNW] + f[BNE] - f[BSW] - f[BSE] 
               + f[BNW]; 

               vz = f[T] - f[B] + f[TE] - f[BW] - f[BE] + f[TW] + f[TN] - f[BS] - f[BN] 
               + f[TS] + f[TNE] + f[TSW] + f[TSE] + f[TNW] - f[BNE] - f[BSW] - f[BSE] 
               - f[BNW];

               vx/=rho;
               vy/=rho;
               vz/=rho;
                  
               durchrho=1.0/rho;

               //forcing 
               /////////////////////////////////////////////////////////////////////////////////////////////
               if (withForcing)
               {
                  muX1 = static_cast<double>(x1+ix1*bcArrayMaxX1);
                  muX2 = static_cast<double>(x2+ix2*bcArrayMaxX2);
                  muX3 = static_cast<double>(x3+ix3*bcArrayMaxX3);

                  forcingX1 = muForcingX1.Eval();
                  forcingX2 = muForcingX2.Eval();
                  forcingX3 = muForcingX3.Eval();

                  //mu100 += forcingX1; 
                  vx += forcingX1*0.5; // X
                  //mu010 += forcingX2; 
                  vy += forcingX2*0.5; // Y
                  //mu001 += forcingX3; 
                  vz += forcingX3*0.5; // Z
               }
               /////////////////////////////////////////////////////////////////////////////////////////////


               M_zXX=0.0;     M_zYY=0.0;     M_zZZ=0.0;     M_zXY=0.0;    M_zXZ=0.0;  M_zYZ=0.0;
               M_zXXY=0.0;    M_zXYY=0.0;    M_zXXZ=0.0;    M_zXZZ=0.0;   M_zYYZ=0.0;  M_zYZZ=0.0;  M_zXYZ=0.0;
               M_zXXYY=0.0;   M_zXXZZ=0.0;   M_zYYZZ=0.0;   M_zXXYZ=0.0;  M_zXYYZ=0.0;  M_zXYZZ=0.0;
               M_zXXYYZ=0.0;  M_zXXYZZ=0.0;  M_zXYYZZ=0.0;  M_zXXYYZZ=0.0;

               vx_sq=vx*vx;
               vy_sq=vy*vy;
               vz_sq=vz*vz;
               vx_vy=vx*vy;
               vx_vz=vx*vz;
               vy_vz=vy*vz;
               vx_vy_vz=vx_vy*vz;

               mu200=(f[E]+f[W]+f[NE]+f[SW]+f[SE]+f[NW]+f[TE]+f[BW]+f[BE]+f[TW]+f[TNE]+f[BNE]+f[TSE]+f[BSE]+f[TNW]+f[BNW]+f[TSW]+f[BSW])*durchrho;
               mu020=(f[N]+f[S]+f[NE]+f[SW]+f[SE]+f[NW]+f[TN]+f[BS]+f[BN]+f[TS]+f[TNE]+f[BNE]+f[TSE]+f[BSE]+f[TNW]+f[BNW]+f[TSW]+f[BSW])*durchrho;
               mu002=(f[T]+f[B]+f[TE]+f[BW]+f[BE]+f[TW]+f[TN]+f[BS]+f[BN]+f[TS]+f[TNE]+f[BNE]+f[TSE]+f[BSE]+f[TNW]+f[BNW]+f[TSW]+f[BSW])*durchrho;
               mu110=(f[NE]-f[SE]+f[SW]-f[NW]+f[TNE]-f[TSE]+f[BNE]-f[BSE]+f[TSW]-f[TNW]+f[BSW]-f[BNW])*durchrho;
               mu101=(f[TE]+f[BW]-f[BE]-f[TW]+f[TNE]-f[BNE]+f[TSE]-f[BSE]-f[TNW]+f[BNW]-f[TSW]+f[BSW])*durchrho;
               mu011=(f[TN]+f[BS]-f[BN]-f[TS]+f[TNE]-f[BNE]-f[TSE]+f[BSE]+f[TNW]-f[BNW]-f[TSW]+f[BSW])*durchrho;
               mu210=(f[NE]-f[SW]-f[SE]+f[NW]+f[TNE]+f[BNE]-f[TSE]-f[BSE]+f[TNW]+f[BNW]-f[TSW]-f[BSW])*durchrho;
               mu120=(f[NE]-f[SW]+f[SE]-f[NW]+f[TNE]-f[TNW]+f[BNE]-f[BNW]+f[TSE]-f[TSW]+f[BSE]-f[BSW])*durchrho;
               mu102=(f[TE]-f[BW]+f[BE]-f[TW]+f[TNE]-f[TNW]+f[BNE]-f[BNW]+f[TSE]-f[TSW]+f[BSE]-f[BSW])*durchrho;
               mu111=(f[TNE]-f[BNE]-f[TSE]+f[BSE]-f[TNW]+f[BNW]+f[TSW]-f[BSW])*durchrho;
               mu201=(f[TE]-f[BW]-f[BE]+f[TW]+f[TNE]-f[BNE]+f[TSE]-f[BSE]+f[TNW]-f[BNW]+f[TSW]-f[BSW])*durchrho;
               mu021=(f[TN]-f[BS]-f[BN]+f[TS]+f[TNE]-f[BNE]+f[TSE]-f[BSE]+f[TNW]-f[BNW]+f[TSW]-f[BSW])*durchrho;
               mu012=(f[TN]-f[BS]+f[BN]-f[TS]+f[TNE]+f[BNE]-f[TSE]-f[BSE]+f[TNW]+f[BNW]-f[TSW]-f[BSW])*durchrho;
               mu220=(f[NE]+f[SW]+f[SE]+f[NW]+f[TNE]+f[BNE]+f[TSE]+f[BSE]+f[TNW]+f[BNW]+f[TSW]+f[BSW])*durchrho;
               mu121=(f[TNE]-f[BNE]+f[TSE]-f[BSE]-f[TNW]+f[BNW]-f[TSW]+f[BSW])*durchrho;
               mu202=(f[TE]+f[BW]+f[BE]+f[TW]+f[TNE]+f[BNE]+f[TSE]+f[BSE]+f[TNW]+f[BNW]+f[TSW]+f[BSW])*durchrho;
               mu211=(f[TNE]-f[BNE]-f[TSE]+f[BSE]+f[TNW]-f[BNW]-f[TSW]+f[BSW])*durchrho;
               mu112=(f[TNE]+f[BNE]-f[TSE]-f[BSE]-f[TNW]-f[BNW]+f[TSW]+f[BSW])*durchrho;
               mu022=(f[TN]+f[BS]+f[BN]+f[TS]+f[TNE]+f[BNE]+f[TSE]+f[BSE]+f[TNW]+f[BNW]+f[TSW]+f[BSW])*durchrho;
               mu221=(f[TNE]-f[BNE]+f[TSE]-f[BSE]+f[TNW]-f[BNW]+f[TSW]-f[BSW])*durchrho;
               mu122=(f[TNE]-f[TNW]+f[BNE]-f[BNW]+f[TSE]-f[TSW]+f[BSE]-f[BSW])*durchrho;
               mu212=(f[TNE]+f[BNE]-f[TSE]-f[BSE]+f[TNW]+f[BNW]-f[TSW]-f[BSW])*durchrho;
               mu222=(f[TNE]+f[BNE]+f[TSE]+f[BSE]+f[TNW]+f[BNW]+f[TSW]+f[BSW])*durchrho;
               mu000=1.0;
               mu100=vx;
               mu010=vy; 
               mu001=vz;

               M_zYZ= -vy_vz  +mu011 ;
               M_zZZ= -vz_sq  +mu002 ; 
               M_zYY= -vy_sq  +mu020 ;
               M_zXX= -vx_sq  +mu200 ;
               M_zXY= -vx_vy  +mu110 ;          
               M_zXZ =-vx_vz  +mu101 ;

               M_zYYZ=mu021-(mu020*vz + 2.0*vy*M_zYZ);          
               M_zYZZ=mu012-(vy*mu002 + 2.0*vz*M_zYZ);
               M_zXZZ=mu102-(vx*mu002 + 2.0*vz*M_zXZ);
               M_zXXZ=mu201-(mu200*vz + 2.0*vx*M_zXZ );
               M_zXXY=mu210-(mu200*vy + 2.0*vx*M_zXY );
               M_zXYY=mu120-(vx*mu020 + 2.0*vy*M_zXY);
               M_zXYZ=mu111-(vx_vy_vz + vx*M_zYZ + vy*M_zXZ + vz*M_zXY);

               M_zXYYZ=mu121-(vz*mu120+ 2.0*vx_vy*M_zYZ+vx*M_zYYZ+vy_sq*M_zXZ+2.0*vy*M_zXYZ); 
               M_zXYZZ=mu112-(vy*mu102+ 2.0*vx_vz*M_zYZ+vx*M_zYZZ+vz_sq*M_zXY+2.0*vz*M_zXYZ);
               M_zXXYZ=mu211-(vz*mu210+ 2.0*vx_vy*M_zXZ+vy*M_zXXZ+vx_sq*M_zYZ+2.0*vx*M_zXYZ); 


               M_zXXYY=mu220-( vx_sq*mu020+ 4.0*vx_vy*M_zXY + 2.0*vx*M_zXYY + vy_sq*M_zXX + 2.0*vy*M_zXXY); 
               M_zYYZZ=mu022-( vy_sq*mu002+ 4.0*vy_vz*M_zYZ + 2.0*vy*M_zYZZ + vz_sq*M_zYY + 2.0*vz*M_zYYZ);
               M_zXXZZ=mu202-( vx_sq*mu002+ 4.0*vx_vz*M_zXZ + 2.0*vx*M_zXZZ + vz_sq*M_zXX + 2.0*vz*M_zXXZ); 

               M_zXXYYZ=mu221-(vz*mu220+
                  2.0*vx_sq*vy*M_zYZ  + vx_sq*M_zYYZ + 2.0*vx*vy_sq*M_zXZ  +4.0*vx_vy*M_zXYZ+ 
                  2.0*vx*M_zXYYZ  +    vy_sq*M_zXXZ  + 2.0*vy*M_zXXYZ);
               M_zXYYZZ=mu122-(vx*mu022 
                  + 2.0*vy_sq*vz*M_zXZ + vy_sq*M_zXZZ + 
                  2.0*vy*vz_sq*M_zXY + 4.0*vy_vz*M_zXYZ + 2.0*vy*M_zXYZZ +vz_sq*M_zXYY + 2.0*vz*M_zXYYZ);
               M_zXXYZZ=mu212-(vy*mu202
                  + 2.0*vx_sq*vz*M_zYZ +vx_sq*M_zYZZ  + 2.0*vx*vz_sq*M_zXY + 4.0*vx_vz*M_zXYZ 
                  + 2.0*vx*M_zXYZZ  + vz_sq*M_zXXY +  2.0*vz*M_zXXYZ);

               M_zXXYYZZ=mu222-(vx_sq*mu022 
                  + 4.0* vx*vy_sq*vz*M_zXZ   + 
                  2.0* vx*vy_sq*M_zXZZ + 4.0* vx_vy*vz_sq*M_zXY+ 
                  8.0* vx_vy_vz*M_zXYZ + 4.0* vx_vy*M_zXYZZ + 
                  2.0* vx*vz_sq*M_zXYY + 4.0* vx_vz*M_zXYYZ + 2.0*vx*M_zXYYZZ + 
                  vy_sq*vz_sq*M_zXX + 2.0* vy_sq*vz*M_zXXZ  + 
                  vy_sq*M_zXXZZ     + 2.0* vy*vz_sq*M_zXXY  + 4.0*vy_vz*M_zXXYZ + 
                  2.0* vy*M_zXXYZZ + vz_sq*M_zXXYY +2.0* vz*M_zXXYYZ);

               ////////lin kombi bilden:
               MXXpMYYpMZZ = M_zXX + M_zYY + M_zZZ;
               MXXmMYY = M_zXX - M_zYY;
               MXXmMZZ = M_zXX - M_zZZ;

               MXXYpMYZZ  =  M_zXXY+M_zYZZ;
               MXXYmMYZZ  =  M_zXXY-M_zYZZ;
               MXXZpMYYZ  =  M_zXXZ+M_zYYZ;
               MXXZmMYYZ  =  M_zXXZ-M_zYYZ;
               MXYYpMXZZ  =  M_zXYY+M_zXZZ;
               MXYYmMXZZ  =  M_zXYY-M_zXZZ;

               MXXYYppp  = M_zXXYY + M_zXXZZ + M_zYYZZ;
               MXXYYpm2p = M_zXXYY - 2.0*M_zXXZZ + M_zYYZZ;
               MXXYYppm2 = M_zXXYY + M_zXXZZ - 2.0*M_zYYZZ;


               //relaxation: 
               MXXpMYYpMZZ -= w2*(1.0-MXXpMYYpMZZ);
               MXXmMYY +=     w1*(    MXXmMYY); 
               MXXmMZZ +=     w1*(    MXXmMZZ); 

               MXXYpMYZZ  +=  w3*(MXXYpMYZZ); 
               MXXYmMYZZ  +=  w4*(MXXYmMYZZ); 
               MXXZpMYYZ  +=  w3*(MXXZpMYYZ); 
               MXXZmMYYZ  +=  w4*(MXXZmMYYZ); 
               MXYYpMXZZ  +=  w3*(MXYYpMXZZ); 
               MXYYmMXZZ  +=  w4*(MXYYmMXZZ); 

               MXXYYppp  -= w7*(c1o3 -MXXYYppp );
               MXXYYpm2p += w6*(MXXYYpm2p);
               MXXYYppm2 += w6*(MXXYYppm2);


               M_zXXYYZZ-=w10*(c1o27-M_zXXYYZZ);
               M_zXZ+=w1*(  M_zXZ    );
               M_zYZ+=w1*(  M_zYZ    );
               M_zXY+=w1*(  M_zXY    );

               M_zXYZ+=w5*(  M_zXYZ   );

               M_zXYYZ+=w8*(  M_zXYYZ  );
               M_zXYZZ+=w8*(  M_zXYZZ  );
               M_zXXYZ+=w8*(  M_zXXYZ  );

               M_zXXYYZ+=w9*( M_zXXYYZ );
               M_zXXYZZ+=w9*( M_zXXYZZ );
               M_zXYYZZ+=w9*( M_zXYYZZ );

               ////////von Lin Kombis zurueck:
               M_zXX = c1o3 *MXXmMYY + c1o3 *MXXmMZZ +  c1o3* MXXpMYYpMZZ; 
               M_zYY = -c2o3* MXXmMYY+ c1o3* MXXmMZZ +  c1o3* MXXpMYYpMZZ; 
               M_zZZ = c1o3   * MXXmMYY -c2o3 *MXXmMZZ +  c1o3* MXXpMYYpMZZ; 

               M_zXXY = (MXXYmMYZZ + MXXYpMYZZ)*0.5; 
               M_zYZZ = 0.5*(-MXXYmMYZZ + MXXYpMYZZ); 
               M_zXYY =(MXYYmMXZZ + MXYYpMXZZ)*0.5; 
               M_zXZZ = 0.5*(-MXYYmMXZZ + MXYYpMXZZ); 
               M_zXXZ = (MXXZmMYYZ + MXXZpMYYZ)*0.5; 
               M_zYYZ = 0.5*(-MXXZmMYYZ + MXXZpMYYZ);

               M_zXXYY =  c1o3* MXXYYpm2p + c1o3* MXXYYppm2 +  0.33333333*MXXYYppp; 
               M_zXXZZ = -c1o3* MXXYYpm2p + c1o3* MXXYYppp; 
               M_zYYZZ = -c1o3* MXXYYppm2 + c1o3* MXXYYppp;

               //forcing 
               /////////////////////////////////////////////////////////////////////////////////////////////
               if (withForcing)
               {
                  mu100 += forcingX1*0.5; 
                  vx += forcingX1*0.5; // X
                  mu010 += forcingX2*0.5; 
                  vy += forcingX2*0.5; // Y
                  mu001 += forcingX3*0.5; 
                  vz += forcingX3*0.5; // Z
                  vx_sq=vx*vx;
                  vy_sq=vy*vy;
                  vz_sq=vz*vz;
                  vx_vy=vx*vy;
                  vx_vz=vx*vz;
                  vy_vz=vy*vz;
                  vx_vy_vz=vx_vy*vz;
               }
               /////////////////////////////////////////////////////////////////////////////////////////////
                  
               mu011=  vy_vz + M_zYZ;
               mu002= vz_sq + M_zZZ; 
               mu020= vy_sq + M_zYY ;
               mu200=   vx_sq + M_zXX;
               mu110= vx_vy+ M_zXY;          
               mu101 =vx_vz + M_zXZ;

               mu021=mu020*vz + 2.0*vy*M_zYZ  + M_zYYZ;          
               mu012=vy*mu002 + 2.0*vz*M_zYZ + M_zYZZ;
               mu102=vx*mu002 + 2.0*vz*M_zXZ + M_zXZZ;
               mu112=vy*mu102 + 2.0*vx_vz*M_zYZ+vx*M_zYZZ+vz_sq*M_zXY+2.0*vz*M_zXYZ+M_zXYZZ;

               mu201=mu200*vz + 2.0*vx*M_zXZ  + M_zXXZ;
               mu210=mu200*vy + 2.0*vx*M_zXY  + M_zXXY;
               mu211=vz*mu210 + 2.0*vx_vy*M_zXZ+vy*M_zXXZ+vx_sq*M_zYZ+2.0*vx*M_zXYZ+M_zXXYZ; 

               mu120=vx*mu020 + 2.0*vy*M_zXY + M_zXYY;
               mu121=vz*mu120 + 2.0*vx_vy*M_zYZ+vx*M_zYYZ+vy_sq*M_zXZ+2.0*vy*M_zXYZ+ M_zXYYZ; 

               mu111=vx_vy_vz + vx*M_zYZ + vy*M_zXZ + vz*M_zXY+M_zXYZ;

               mu220=vx_sq*mu020 + 4.0*vx_vy*M_zXY + 2.0*vx*M_zXYY  + vy_sq*M_zXX + 2.0*vy*M_zXXY +  M_zXXYY; 

               mu221=vz*mu220 + 2.0*vx_sq*vy*M_zYZ  + vx_sq*M_zYYZ + 2.0*vx*vy_sq*M_zXZ  +4.0*vx_vy*M_zXYZ+ 
                  2.0*vx*M_zXYYZ  +    vy_sq*M_zXXZ  + 2.0*vy*M_zXXYZ  + M_zXXYYZ;

               mu022= vy_sq*mu002 + 4.0*vy_vz*M_zYZ  +2.0*vy*M_zYZZ+ vz_sq*M_zYY+ 2.0*vz*M_zYYZ +  M_zYYZZ;

               mu122=vx*mu022 + 2.0*vy_sq*vz*M_zXZ + vy_sq*M_zXZZ + 
                  2.0*vy*vz_sq*M_zXY + 4.0*vy_vz*M_zXYZ + 2.0*vy*M_zXYZZ +vz_sq*M_zXYY + 2.0*vz*M_zXYYZ + M_zXYYZZ;

               mu202= vx_sq*mu002 + 4.0*vx_vz*M_zXZ+ 2.0*vx*M_zXZZ + vz_sq*M_zXX+ 2.0*vz*M_zXXZ + M_zXXZZ; 

               mu212=vy*mu202 + 2.0*vx_sq*vz*M_zYZ +vx_sq*M_zYZZ  + 2.0*vx*vz_sq*M_zXY + 4.0*vx_vz*M_zXYZ 
                  + 2.0*vx*M_zXYZZ  + vz_sq*M_zXXY +  2.0*vz*M_zXXYZ + M_zXXYZZ;

               mu222=vx_sq*mu022 
                  + 4.0* vx*vy_sq*vz*M_zXZ   + 
                  2.0* vx*vy_sq*M_zXZZ + 4.0* vx_vy*vz_sq*M_zXY+ 
                  8.0* vx_vy_vz*M_zXYZ + 4.0* vx_vy*M_zXYZZ + 
                  2.0* vx*vz_sq*M_zXYY + 4.0* vx_vz*M_zXYYZ + 2.0*vx*M_zXYYZZ + 
                  vy_sq*vz_sq*M_zXX + 2.0* vy_sq*vz*M_zXXZ  + 
                  vy_sq*M_zXXZZ     + 2.0* vy*vz_sq*M_zXXY  + 4.0*vy_vz*M_zXXYZ + 
                  2.0* vy*M_zXXYZZ + vz_sq*M_zXXYY +2.0* vz*M_zXXYYZ + M_zXXYYZZ;

               f[E] =0.5* (mu200 -  mu220 + mu222 - mu202 - mu120 + mu122 - mu102 +mu100) *(rho);
               f[W] =0.5* (mu200 - mu220 + mu222 - mu202 + mu120 - mu122 + mu102 -mu100)*  (rho);
               f[N] =0.5* (-mu210 - mu220 + mu222 + mu212 + mu020 - mu022 - mu012 +mu010)* (rho);
               f[S] =0.5* (mu210 -  mu220 + mu222 - mu212 + mu020 - mu022 + mu012 -mu010)* (rho);
               f[T] =0.5* (mu221 +  mu222 - mu201 - mu202 - mu021 - mu022 + mu002 +mu001)* (rho);
               f[B] =0.5* (-mu221 + mu222 + mu201  - mu202 + mu021 - mu022 + mu002-mu001)* (rho);
               f[NE] =0.25*( mu210  + mu220- mu222 - mu212 + mu110+ mu120- mu122 -mu112)*  (rho); 
               f[SW] =0.25*(-mu210 + mu220- mu222 + mu212 + mu110- mu120+ mu122 -mu112)*  (rho); 
               f[SE] =0.25*(-mu210 + mu220- mu222 + mu212 - mu110+ mu120- mu122 +mu112)*  (rho); 
               f[NW] =0.25*( mu210  + mu220- mu222 - mu212 - mu110- mu120+ mu122 + mu112)* (rho); 
               f[TE]=0.25*(-mu221 - mu222 + mu201 + mu202 - mu121 - mu122 + mu101 + mu102)*(rho); 
               f[BW]=0.25*( mu221  -mu222 - mu201 + mu202 - mu121 + mu122 + mu101 - mu102)*(rho);

               f[BE]=0.25*(mu221 - mu222 - mu201 + mu202 + mu121 - mu122 - mu101 +mu102)*  (rho);
               f[TW]=0.25*(-mu221 - mu222 + mu201 + mu202 + mu121 + mu122 - mu101 -mu102)* (rho); 
               f[TN]=0.25*(-mu221 - mu222 - mu211 - mu212 + mu021 + mu022 + mu011+mu012)*  (rho);
               f[BS]=0.25*( mu221 - mu222 - mu211 + mu212 - mu021 + mu022 + mu011 - mu012)*(rho);
               f[BN]=0.25*( mu221 - mu222 + mu211 - mu212 - mu021 + mu022 - mu011 + mu012)*(rho);
               f[TS]=0.25*(-mu221 - mu222 + mu211 + mu212 + mu021 + mu022 - mu011 -mu012)* (rho); 
               f[ZERO]=     (-mu200 + mu220 - mu222 + mu202 - mu020 + mu022 - mu002 + mu000)*(rho);
               f[TNE]=0.125*( mu221 + mu222 + mu211 + mu212 + mu121 + mu122 + mu111 + mu112)*  (rho); 
               f[BNE]=0.125*(-mu221 + mu222 -mu211 + mu212 -mu121 + mu122 -mu111 + mu112)*     (rho);
               f[TSE]=0.125*( mu221 + mu222 - mu211 - mu212 + mu121 + mu122 - mu111 - mu112)*  (rho); 
               f[BSE]=0.125*(-mu221 + mu222 +mu211 - mu212 -mu121 + mu122 +mu111 - mu112)*     (rho); 
               f[TNW]=0.125*( mu221 + mu222 + mu211 + mu212 - mu121 - mu122 - mu111 - mu112)*  (rho); 
               f[BNW]=0.125*(-mu221 + mu222 -mu211 + mu212 +mu121 - mu122 +mu111 - mu112)*     (rho); 
               f[TSW]=0.125*( mu221 + mu222 - mu211 - mu212 - mu121 - mu122 + mu111 + mu112)*  (rho); 
               f[BSW]=0.125*(-mu221 + mu222+mu211 - mu212+mu121 - mu122-mu111 + mu112)*        (rho);
               //////////////////////////////////////////////////////////////////////////
               //proof correctness
               //////////////////////////////////////////////////////////////////////////
               LBMReal rho_post = f[ZERO] + f[E] + f[W] + f[N] + f[S] + f[T] + f[B] 
               + f[NE] + f[SW] + f[SE] + f[NW] + f[TE] + f[BW] + f[BE]
               + f[TW] + f[TN] + f[BS] + f[BN] + f[TS] + f[TNE] + f[TSW]
               + f[TSE] + f[TNW] + f[BNE] + f[BSW] + f[BSE] + f[BNW];
               //LBMReal dif = fabs(rho - rho_post);
               LBMReal dif = rho - rho_post;
#ifdef SINGLEPRECISION
               if(dif > 10.0E-7 || dif < -10.0E-7)
#else
               if(dif > 10.0E-15 || dif < -10.0E-15)
#endif
               {
                  UB_THROW(UbException(UB_EXARGS,"rho is not correct for node "+UbSystem::toString(x1)+","+UbSystem::toString(x2)+","+UbSystem::toString(x3)));
                  //UBLOG(logERROR,"LBMKernelETD3Q27Cascaded::collideAll(): rho is not correct for node "+UbSystem::toString(x1)+","+UbSystem::toString(x2)+","+UbSystem::toString(x3));
                  //exit(EXIT_FAILURE);
               }
               //////////////////////////////////////////////////////////////////////////
               //write distribution
               //////////////////////////////////////////////////////////////////////////
               (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3) = f[D3Q27System::INV_E];
               (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3) = f[D3Q27System::INV_N];
               (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3) = f[D3Q27System::INV_T];
               (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3) = f[D3Q27System::INV_NE];
               (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3) = f[D3Q27System::INV_NW];
               (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3) = f[D3Q27System::INV_TE];
               (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3) = f[D3Q27System::INV_TW];
               (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3) = f[D3Q27System::INV_TN];
               (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3) = f[D3Q27System::INV_TS];
               (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3) = f[D3Q27System::INV_TNE];
               (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3) = f[D3Q27System::INV_TNW];
               (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3) = f[D3Q27System::INV_TSE];
               (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3) = f[D3Q27System::INV_TSW];

               (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = f[D3Q27System::INV_W ];
               (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = f[D3Q27System::INV_S ];
               (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = f[D3Q27System::INV_B ];
               (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = f[D3Q27System::INV_SW];
               (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = f[D3Q27System::INV_SE];
               (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = f[D3Q27System::INV_BW];
               (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = f[D3Q27System::INV_BE];
               (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = f[D3Q27System::INV_BS];
               (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = f[D3Q27System::INV_BN];
               (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = f[D3Q27System::INV_BSW];
               (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = f[D3Q27System::INV_BSE];
               (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = f[D3Q27System::INV_BNW];
               (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1) = f[D3Q27System::INV_BNE];

               (*this->zeroDistributions)(x1,x2,x3) = f[D3Q27System::ZERO];
               //////////////////////////////////////////////////////////////////////////
               
            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
double LBMKernelETD3Q27Cascaded::getCallculationTime()
{
   //return timer.getDuration();
   return timer.getTotalTime();
}
