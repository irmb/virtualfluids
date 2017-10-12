#include "SlipBCAdapter.h"
#include "LBM/D3Q27System.h"
#include "Interactors//D3Q27Interactor.h"
#include "numerics/geometry3d/GbCuboid3D.h"
#include <boost/pointer_cast.hpp>

//*==========================================================*/
//ObObject* D3Q27SlipBCAdapterCreator::createObObject()
//{
//   return new D3Q27SlipBCAdapter; 
//}
//*==========================================================*/
//ObObjectCreator* D3Q27SlipBCAdapter::getCreator()
//{
//   return D3Q27SlipBCAdapterCreator::getInstance();
//}
//*==========================================================*/
void SlipBCAdapter::adaptBC(const D3Q27Interactor& interactor, BoundaryConditionsPtr bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& time)
{
   //////////////////////////////////////////////////////////////////////////
   //>>> nur workaround! -> Hendrick nach normalen berechnung aus qs fragen
   
   GbCuboid3DPtr geo = boost::dynamic_pointer_cast<GbCuboid3D>(interactor.getGbObject3D());
   if(!geo) throw UbException(UB_EXARGS,"derzeit nur fuer Cubes valide");

   if     ( bc->hasSlipBoundaryFlag(D3Q27System::E) ) bc->setNormalVector( 1.0, 0.0, 0.0);  
   else if( bc->hasSlipBoundaryFlag(D3Q27System::W) ) bc->setNormalVector(-1.0, 0.0, 0.0);  
   else if( bc->hasSlipBoundaryFlag(D3Q27System::N) ) bc->setNormalVector( 0.0, 1.0, 0.0);  
   else if( bc->hasSlipBoundaryFlag(D3Q27System::S) ) bc->setNormalVector( 0.0,-1.0, 0.0);  
   else if( bc->hasSlipBoundaryFlag(D3Q27System::T) ) bc->setNormalVector( 0.0, 0.0, 1.0);  
   else if( bc->hasSlipBoundaryFlag(D3Q27System::B) ) bc->setNormalVector( 0.0, 0.0,-1.0);  

   bc->setBcAlgorithmType(algorithmType);
}

