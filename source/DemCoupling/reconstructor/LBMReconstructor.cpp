#include "LBMReconstructor.h"

#include "ILBMKernel.h"
#include "D3Q27System.h"
#include "DataSet3D.h"
#include "BCProcessor.h"
#include "BCArray3D.h"

#include "PhysicsEngineGeometryAdapter.h"

using namespace D3Q27System;

LBMReconstructor::LBMReconstructor(bool compressible)
{
   if (compressible)
   {
      calcMacrosFct = &D3Q27System::calcCompMacroscopicValues;
   }
   else
   {
      calcMacrosFct = &D3Q27System::calcIncompMacroscopicValues;
   }
}

void LBMReconstructor::reconstructNode(const int& x1, const int& x2, const int& x3,
                                               const Vector3D& worldCoordinates, std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry, std::shared_ptr<ILBMKernel> kernel) const
{
   LBMReal pdf[D3Q27System::ENDF + 1];

   LBMReal rho, vx1, vx2, vx3;
   calcMacrosFct(pdf, rho, vx1, vx2, vx3);

   LBMReal rho_dif = 1; 

   while (rho_dif > 1e-7)
   {
      for (int fDir = D3Q27System::FSTARTDIR; fDir <= D3Q27System::FENDDIR; fDir++)
      {

         UbTupleInt3 neighbor(x1 + D3Q27System::DX1[fDir], x2 + D3Q27System::DX2[fDir], x3 + D3Q27System::DX3[fDir]);

         if (!kernel->getBCProcessor()->getBCArray()->isFluid(val<1>(neighbor), val<2>(neighbor), val<3>(neighbor)))
         {
            LBMReal pdfNeighbor[D3Q27System::ENDF + 1];
            SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
            const int invDir = D3Q27System::INVDIR[fDir];
            distributions->getDistributionForDirection(pdfNeighbor[invDir], val<1>(neighbor), val<2>(neighbor), val<3>(neighbor));
            distributions->setDistributionInvForDirection(pdf[invDir], x1, x2, x3, invDir);
         }


      }
   }



   LBMReal collFactor = kernel->getCollisionFactor();
   collide(pdf, collFactor);




}


void LBMReconstructor::collide(LBMReal* f, LBMReal collFactor)
{

   LBMReal drho, vx1, vx2, vx3;
   LBMReal feq[D3Q27System::ENDF+1];


   drho = ((f[TNE]+f[BSW])+(f[TSE]+f[BNW]))+((f[BSE]+f[TNW])+(f[TSW]+f[BNE]))
      +(((f[NE]+f[SW])+(f[SE]+f[NW]))+((f[TE]+f[BW])+(f[BE]+f[TW]))
         +((f[BN]+f[TS])+(f[TN]+f[BS])))+((f[E]+f[W])+(f[N]+f[S])
            +(f[T]+f[B]))+f[ZERO];

   vx1 = ((((f[TNE]-f[BSW])+(f[TSE]-f[BNW]))+((f[BSE]-f[TNW])+(f[BNE]-f[TSW])))+
      (((f[BE]-f[TW])+(f[TE]-f[BW]))+((f[SE]-f[NW])+(f[NE]-f[SW])))+
      (f[E]-f[W]));

   vx2 = ((((f[TNE]-f[BSW])+(f[BNW]-f[TSE]))+((f[TNW]-f[BSE])+(f[BNE]-f[TSW])))+
      (((f[BN]-f[TS])+(f[TN]-f[BS]))+((f[NW]-f[SE])+(f[NE]-f[SW])))+
      (f[N]-f[S]));

   vx3 = ((((f[TNE]-f[BSW])+(f[TSE]-f[BNW]))+((f[TNW]-f[BSE])+(f[TSW]-f[BNE])))+
      (((f[TS]-f[BN])+(f[TN]-f[BS]))+((f[TW]-f[BE])+(f[TE]-f[BW])))+
      (f[T]-f[B]));

   LBMReal cu_sq = 1.5*(vx1*vx1+vx2*vx2+vx3*vx3);

   feq[ZERO] = c8o27*(drho-cu_sq);
   feq[E] = c2o27*(drho+3.0*(vx1)+c9o2*(vx1)*(vx1)-cu_sq);
   feq[W] = c2o27*(drho+3.0*(-vx1)+c9o2*(-vx1)*(-vx1)-cu_sq);
   feq[N] = c2o27*(drho+3.0*(vx2)+c9o2*(vx2)*(vx2)-cu_sq);
   feq[S] = c2o27*(drho+3.0*(-vx2)+c9o2*(-vx2)*(-vx2)-cu_sq);
   feq[T] = c2o27*(drho+3.0*(vx3)+c9o2*(vx3)*(vx3)-cu_sq);
   feq[B] = c2o27*(drho+3.0*(-vx3)+c9o2*(-vx3)*(-vx3)-cu_sq);
   feq[NE] = c1o54*(drho+3.0*(vx1+vx2)+c9o2*(vx1+vx2)*(vx1+vx2)-cu_sq);
   feq[SW] = c1o54*(drho+3.0*(-vx1-vx2)+c9o2*(-vx1-vx2)*(-vx1-vx2)-cu_sq);
   feq[SE] = c1o54*(drho+3.0*(vx1-vx2)+c9o2*(vx1-vx2)*(vx1-vx2)-cu_sq);
   feq[NW] = c1o54*(drho+3.0*(-vx1+vx2)+c9o2*(-vx1+vx2)*(-vx1+vx2)-cu_sq);
   feq[TE] = c1o54*(drho+3.0*(vx1+vx3)+c9o2*(vx1+vx3)*(vx1+vx3)-cu_sq);
   feq[BW] = c1o54*(drho+3.0*(-vx1-vx3)+c9o2*(-vx1-vx3)*(-vx1-vx3)-cu_sq);
   feq[BE] = c1o54*(drho+3.0*(vx1-vx3)+c9o2*(vx1-vx3)*(vx1-vx3)-cu_sq);
   feq[TW] = c1o54*(drho+3.0*(-vx1+vx3)+c9o2*(-vx1+vx3)*(-vx1+vx3)-cu_sq);
   feq[TN] = c1o54*(drho+3.0*(vx2+vx3)+c9o2*(vx2+vx3)*(vx2+vx3)-cu_sq);
   feq[BS] = c1o54*(drho+3.0*(-vx2-vx3)+c9o2*(-vx2-vx3)*(-vx2-vx3)-cu_sq);
   feq[BN] = c1o54*(drho+3.0*(vx2-vx3)+c9o2*(vx2-vx3)*(vx2-vx3)-cu_sq);
   feq[TS] = c1o54*(drho+3.0*(-vx2+vx3)+c9o2*(-vx2+vx3)*(-vx2+vx3)-cu_sq);
   feq[TNE] = c1o216*(drho+3.0*(vx1+vx2+vx3)+c9o2*(vx1+vx2+vx3)*(vx1+vx2+vx3)-cu_sq);
   feq[BSW] = c1o216*(drho+3.0*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
   feq[BNE] = c1o216*(drho+3.0*(vx1+vx2-vx3)+c9o2*(vx1+vx2-vx3)*(vx1+vx2-vx3)-cu_sq);
   feq[TSW] = c1o216*(drho+3.0*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
   feq[TSE] = c1o216*(drho+3.0*(vx1-vx2+vx3)+c9o2*(vx1-vx2+vx3)*(vx1-vx2+vx3)-cu_sq);
   feq[BNW] = c1o216*(drho+3.0*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
   feq[BSE] = c1o216*(drho+3.0*(vx1-vx2-vx3)+c9o2*(vx1-vx2-vx3)*(vx1-vx2-vx3)-cu_sq);
   feq[TNW] = c1o216*(drho+3.0*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);

   //Relaxation
   f[ZERO] += (feq[ZERO]-f[ZERO])*collFactor;
   f[E] += (feq[E]-f[E])*collFactor;
   f[W] += (feq[W]-f[W])*collFactor;
   f[N] += (feq[N]-f[N])*collFactor;
   f[S] += (feq[S]-f[S])*collFactor;
   f[T] += (feq[T]-f[T])*collFactor;
   f[B] += (feq[B]-f[B])*collFactor;
   f[NE] += (feq[NE]-f[NE])*collFactor;
   f[SW] += (feq[SW]-f[SW])*collFactor;
   f[SE] += (feq[SE]-f[SE])*collFactor;
   f[NW] += (feq[NW]-f[NW])*collFactor;
   f[TE] += (feq[TE]-f[TE])*collFactor;
   f[BW] += (feq[BW]-f[BW])*collFactor;
   f[BE] += (feq[BE]-f[BE])*collFactor;
   f[TW] += (feq[TW]-f[TW])*collFactor;
   f[TN] += (feq[TN]-f[TN])*collFactor;
   f[BS] += (feq[BS]-f[BS])*collFactor;
   f[BN] += (feq[BN]-f[BN])*collFactor;
   f[TS] += (feq[TS]-f[TS])*collFactor;

   f[TNE] += (feq[TNE]-f[TNE])*collFactor;
   f[BSW] += (feq[BSW]-f[BSW])*collFactor;
   f[BNE] += (feq[BNE]-f[BNE])*collFactor;
   f[TSW] += (feq[TSW]-f[TSW])*collFactor;
   f[TSE] += (feq[TSE]-f[TSE])*collFactor;
   f[BNW] += (feq[BNW]-f[BNW])*collFactor;
   f[BSE] += (feq[BSE]-f[BSE])*collFactor;
   f[TNW] += (feq[TNW]-f[TNW])*collFactor;
}