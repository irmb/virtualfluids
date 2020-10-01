//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __         
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |        
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |        
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |        
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____    
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|   
//      \    \  |    |   ________________________________________________________________    
//       \    \ |    |  |  ______________________________________________________________|   
//        \    \|    |  |  |         __          __     __     __     ______      _______    
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)   
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______    
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/   
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file VelocityBCAdapter.cpp
//! \ingroup BoundarConditions
//! \author Sören Freudiger
//=======================================================================================

#include "VelocityBCAdapter.h"
#include "basics/utilities/UbLogger.h"
#include "basics/utilities/UbMath.h"
#include "basics/utilities/UbTuple.h"

using namespace std;


VelocityBCAdapter::VelocityBCAdapter(const bool& vx1, const bool& vx2, const bool& vx3, const BCFunction& velVxBC)
{
   if(vx1) this->vx1BCs.push_back(velVxBC);
   if(vx2) this->vx2BCs.push_back(velVxBC);
   if(vx3) this->vx3BCs.push_back(velVxBC);
   this->init();
}
/*==========================================================*/
VelocityBCAdapter::VelocityBCAdapter(const bool& vx1, const bool& vx2, const bool& vx3, const mu::Parser& function, const double& startTime, const double& endTime )
{
   if(vx1) this->vx1BCs.push_back(BCFunction(function,startTime,endTime));
   if(vx2) this->vx2BCs.push_back(BCFunction(function,startTime,endTime));
   if(vx3) this->vx3BCs.push_back(BCFunction(function,startTime,endTime));
   this->init();
}
/*==========================================================*/
VelocityBCAdapter::VelocityBCAdapter(const bool& vx1, const bool& vx2, const bool& vx3, const mu::Parser& function1, const mu::Parser& function2, const mu::Parser& function3, const double& startTime, const double& endTime )
{
   if(vx1) this->vx1BCs.push_back(BCFunction(function1,startTime,endTime));
   if(vx2) this->vx2BCs.push_back(BCFunction(function2,startTime,endTime));
   if(vx3) this->vx3BCs.push_back(BCFunction(function3,startTime,endTime));
   this->init();
}
/*==========================================================*/
VelocityBCAdapter::VelocityBCAdapter(const bool& vx1, const bool& vx2, const bool& vx3, const string& functionstring, const double& startTime, const double& endTime )
{
   if(vx1) this->vx1BCs.push_back(BCFunction(functionstring,startTime,endTime));
   if(vx2) this->vx2BCs.push_back(BCFunction(functionstring,startTime,endTime));
   if(vx3) this->vx3BCs.push_back(BCFunction(functionstring,startTime,endTime));
   this->init();
}
/*==========================================================*/
VelocityBCAdapter::VelocityBCAdapter(const BCFunction& velBC, bool x1Dir, bool x2Dir, bool x3Dir)
{
   if(x1Dir) this->vx1BCs.push_back(velBC);
   if(x2Dir) this->vx2BCs.push_back(velBC);
   if(x3Dir) this->vx3BCs.push_back(velBC);
   this->init();
}
/*==========================================================*/
VelocityBCAdapter::VelocityBCAdapter(const BCFunction& velVx1BC, const BCFunction& velVx2BC, const BCFunction& velVx3BC)
{
   if( velVx1BC.getEndTime()!=-Ub::inf ) this->vx1BCs.push_back(velVx1BC);
   if( velVx2BC.getEndTime()!=-Ub::inf ) this->vx2BCs.push_back(velVx2BC);
   if( velVx3BC.getEndTime()!=-Ub::inf ) this->vx3BCs.push_back(velVx3BC);
   this->init();
}
/*==========================================================*/
VelocityBCAdapter::VelocityBCAdapter(const vector< BCFunction >& velVx1BCs, const vector< BCFunction >& velVx2BCs, const vector< BCFunction >& velVx3BCs)
{
   this->vx1BCs = velVx1BCs;
   this->vx2BCs = velVx2BCs;
   this->vx3BCs = velVx3BCs;
   this->init();
}
/*==========================================================*/
VelocityBCAdapter::VelocityBCAdapter(const double& vx1, const double& vx1StartTime, const double& vx1EndTime,
                                               const double& vx2, const double& vx2StartTime, const double& vx2EndTime,
                                               const double& vx3, const double& vx3StartTime, const double& vx3EndTime )
{
   this->vx1BCs.push_back(BCFunction(vx1,vx1StartTime,vx1EndTime));
   this->vx2BCs.push_back(BCFunction(vx2,vx2StartTime,vx2EndTime));
   this->vx3BCs.push_back(BCFunction(vx3,vx3StartTime,vx3EndTime));
   this->init();
}
/*==========================================================*/
VelocityBCAdapter::VelocityBCAdapter(const string& vx1Function, const double& vx1StartTime, const double& vx1EndTime,
                                               const string& vx2Function, const double& vx2StartTime, const double& vx2EndTime,
                                               const string& vx3Function, const double& vx3StartTime, const double& vx3EndTime ) 
{
   if(vx1Function.size()) this->vx1BCs.push_back(BCFunction(vx1Function,vx1StartTime,vx1EndTime));
   if(vx2Function.size()) this->vx2BCs.push_back(BCFunction(vx2Function,vx2StartTime,vx2EndTime));
   if(vx3Function.size()) this->vx3BCs.push_back(BCFunction(vx3Function,vx3StartTime,vx3EndTime));
   this->init();
}
/*==========================================================*/
void VelocityBCAdapter::setNewVelocities(const double& vx1, const double& vx1StartTime, const double& vx1EndTime,
                                              const double& vx2, const double& vx2StartTime, const double& vx2EndTime,
                                              const double& vx3, const double& vx3StartTime, const double& vx3EndTime )
{
   this->clear();
   this->vx1BCs.push_back(BCFunction(vx1,vx1StartTime,vx1EndTime));
   this->vx2BCs.push_back(BCFunction(vx2,vx2StartTime,vx2EndTime));
   this->vx3BCs.push_back(BCFunction(vx3,vx3StartTime,vx3EndTime));
   this->init();
}
/*==========================================================*/
void VelocityBCAdapter::init()
{
   this->unsetTimeDependent();
   
   this->timeStep = 0.0;

   this->x1 = 0.0;
   this->x2 = 0.0;
   this->x3 = 0.0;

   this->tmpVx1Function = NULL;
   this->tmpVx2Function = NULL;
   this->tmpVx3Function = NULL;

   try //initilialization and validation of functions
   {
      this->init(vx1BCs);
      this->init(vx2BCs);
      this->init(vx3BCs);
   }
   catch(mu::Parser::exception_type& e){ stringstream error; error<<"mu::parser exception occurs, message("<<e.GetMsg()<<"), formula("<<e.GetExpr()+"), token("+e.GetToken()<<")"
                                          <<", pos("<<e.GetPos()<<"), error code("<<e.GetCode(); throw UbException(error.str()); }
   catch(...)                          { throw UbException(UB_EXARGS,"unknown exception" ); }
}
/*==========================================================*/
void VelocityBCAdapter::init(std::vector<BCFunction>& vxBCs)
{
   for(size_t pos=0; pos<vxBCs.size(); ++pos)
   {
      if( !(    UbMath::equal( BCFunction::INFCONST, vxBCs[pos].getEndTime() )
             && UbMath::greaterEqual( this->timeStep,  vxBCs[pos].getStartTime()  ) ) )
      {
         this->setTimeDependent();
      }

      vxBCs[pos].getFunction().DefineVar("t" , &this->timeStep);
      vxBCs[pos].getFunction().DefineVar("x1", &this->x1      );
      vxBCs[pos].getFunction().DefineVar("x2", &this->x2      );
      vxBCs[pos].getFunction().DefineVar("x3", &this->x3      );

      vxBCs[pos].getFunction().Eval(); //<-- validation
   }
}
/*==========================================================*/
void VelocityBCAdapter::init(const D3Q27Interactor* const& interactor, const double& time)
{
   this->timeStep       = time;
   this->tmpVx1Function = this->tmpVx2Function = this->tmpVx3Function = NULL;

   //aktuelle velocityfunction bestimmen
   double maxEndtime = -Ub::inf;
   
   for(size_t pos=0; pos<vx1BCs.size(); ++pos)
   {
      if( UbMath::equal(vx1BCs[pos].getEndTime(),BCFunction::INFTIMEDEPENDENT) ) maxEndtime=Ub::inf;
      maxEndtime = UbMath::max(maxEndtime,vx1BCs[pos].getStartTime(),vx1BCs[pos].getEndTime()); //startTime abfragen, da  INFCONST=-10
      
      if( UbMath::greaterEqual(this->timeStep,vx1BCs[pos].getStartTime()) ) 
      {
          if(   UbMath::lessEqual( this->timeStep     , vx1BCs[pos].getEndTime()     )
             || UbMath::equal(     vx1BCs[pos].getEndTime(), (double)BCFunction::INFCONST        )
             || UbMath::equal(     vx1BCs[pos].getEndTime(), (double)BCFunction::INFTIMEDEPENDENT)  )
         {
            tmpVx1Function = &vx1BCs[pos].getFunction();
            break;
         }
      }
   }
   for(size_t pos=0; pos<vx2BCs.size(); ++pos)
   {
      if( UbMath::equal(vx2BCs[pos].getEndTime(),BCFunction::INFTIMEDEPENDENT)) maxEndtime=Ub::inf;
      maxEndtime = UbMath::max(maxEndtime,vx2BCs[pos].getStartTime(),vx2BCs[pos].getEndTime()); //startTime abfragen, da  INFCONST=-10

      if( UbMath::greaterEqual(this->timeStep,vx2BCs[pos].getStartTime()) ) 
      {
         if(   UbMath::lessEqual( this->timeStep     , vx2BCs[pos].getEndTime()      )
            || UbMath::equal(     vx2BCs[pos].getEndTime(), (double)BCFunction::INFCONST         )
            || UbMath::equal(     vx2BCs[pos].getEndTime(), (double)BCFunction::INFTIMEDEPENDENT )  )
         {
            tmpVx2Function = &vx2BCs[pos].getFunction();
            break;
         }
      }
   }
   for(size_t pos=0; pos<vx3BCs.size(); ++pos)
   {
      if( UbMath::equal(vx3BCs[pos].getEndTime(),BCFunction::INFTIMEDEPENDENT)) maxEndtime=Ub::inf;
      maxEndtime = UbMath::max(maxEndtime,vx3BCs[pos].getStartTime(),vx3BCs[pos].getEndTime()); //startTime abfragen, da  INFCONST=-10

      if( UbMath::greaterEqual(this->timeStep,vx3BCs[pos].getStartTime()) ) 
      {
         if(   UbMath::lessEqual( this->timeStep     , vx3BCs[pos].getEndTime()      )
            || UbMath::equal(     vx3BCs[pos].getEndTime(), (double)BCFunction::INFCONST         )
            || UbMath::equal(     vx3BCs[pos].getEndTime(), (double)BCFunction::INFTIMEDEPENDENT )  )
         {
            tmpVx3Function = &vx3BCs[pos].getFunction();
            break;
         }
      }
   }

   if( UbMath::greaterEqual(time,maxEndtime) ) 
   {
      if( !this->isTimePeriodic() ) this->unsetTimeDependent();
      else //bei peridoic die interavalle neu setzen:
      {
         if( UbMath::equal(maxEndtime,BCFunction::INFCONST) )
            for(size_t pos=0; pos<vx1BCs.size(); ++pos)
            {
               vx1BCs[pos].setStartTime( vx1BCs[pos].getStartTime() + timeStep );
               vx1BCs[pos].setEndTime( vx1BCs[pos].getEndTime() + timeStep );
            }
            if( UbMath::equal(maxEndtime,BCFunction::INFCONST) )
            for(size_t pos=0; pos<vx2BCs.size(); ++pos)
            {
               vx2BCs[pos].setStartTime( vx2BCs[pos].getStartTime() + timeStep );
               vx2BCs[pos].setEndTime( vx2BCs[pos].getEndTime() + timeStep );
            }
         if( UbMath::equal(maxEndtime,BCFunction::INFCONST) )
            for(size_t pos=0; pos<vx3BCs.size(); ++pos)
            {
               vx3BCs[pos].setStartTime( vx3BCs[pos].getStartTime() + timeStep );
               vx3BCs[pos].setEndTime( vx3BCs[pos].getEndTime() + timeStep );
            }
        this->init(interactor,time);
      }
   }

   UBLOG(logDEBUG4,"D3Q27VelocityBCAdapter::init(time="<<time<<") "
                   <<", vx1= \""<<(tmpVx1Function ? tmpVx1Function->GetExpr() : "-")<<"\""
                   <<", vx2= \""<<(tmpVx2Function ? tmpVx2Function->GetExpr() : "-")<<"\""
                   <<", vx3= \""<<(tmpVx3Function ? tmpVx3Function->GetExpr() : "-")<<"\""
                   <<", timedependent="<<boolalpha<<this->isTimeDependent()   );
}
/*==========================================================*/
void VelocityBCAdapter::update( const D3Q27Interactor* const& interactor, const double& time ) 
{
   this->init(interactor,time);
}
/*==========================================================*/
void VelocityBCAdapter::adaptBCForDirection( const D3Q27Interactor& interactor, SPtr<BoundaryConditions> bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& q, const int& fdirection, const double& time )
{
   bc->setVelocityBoundaryFlag(D3Q27System::INVDIR[fdirection],secondaryBcOption);
   bc->setQ((float)q,fdirection);
}
/*==========================================================*/
void VelocityBCAdapter::adaptBC( const D3Q27Interactor& interactor, SPtr<BoundaryConditions> bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& time ) 
{
   this->setNodeVelocity(interactor,bc,worldX1,worldX2,worldX3,time);
   bc->setBcAlgorithmType(algorithmType);
}
/*==========================================================*/
void VelocityBCAdapter::setNodeVelocity( const D3Q27Interactor& interactor, SPtr<BoundaryConditions> bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& timestep) 
{
   //Geschwindigkeiten setzen
   try
   {
      //PunktKoordinaten bestimmen
      this->x1 = worldX1;
      this->x2 = worldX2;
      this->x3 = worldX3;
      this->timeStep = timestep;

      if(tmpVx1Function) bc->setBoundaryVelocityX1((LBMReal)tmpVx1Function->Eval());  
      if(tmpVx2Function) bc->setBoundaryVelocityX2((LBMReal)tmpVx2Function->Eval());
      if(tmpVx3Function) bc->setBoundaryVelocityX3((LBMReal)tmpVx3Function->Eval());
   }
   catch(mu::Parser::exception_type& e){ stringstream error; error<<"mu::parser exception occurs, message("<<e.GetMsg()<<"), formula("<<e.GetExpr()+"), token("+e.GetToken()<<")"
                                         <<", pos("<<e.GetPos()<<"), error code("<<e.GetCode(); throw UbException(error.str()); }
   catch(...)                          { throw UbException(UB_EXARGS,"unknown exception" ); }
}
/*==========================================================*/
UbTupleDouble3 VelocityBCAdapter::getVelocity(const double& x1, const double& x2, const double& x3, const double& timeStep) const
{
	double vx1 = 0.0;
	double vx2 = 0.0;
	double vx3 = 0.0;
   this->x1 = x1;
   this->x2 = x2;
   this->x3 = x3;
   this->timeStep = timeStep;
	
	if(tmpVx1Function) vx1 = tmpVx1Function->Eval();  
   if(tmpVx2Function) vx2 = tmpVx2Function->Eval();
   if(tmpVx3Function) vx3 = tmpVx3Function->Eval();
    
   return UbTupleDouble3(vx1,vx2,vx3);
}
/*==========================================================*/
string VelocityBCAdapter::toString()
{
   stringstream info;
   info<<"D3Q27VelocityBCAdapter:\n";
   info<<" #vx1-functions = "<<(int)vx1BCs.size()<<endl;
   info<<" #vx2-functions = "<<(int)vx2BCs.size()<<endl;
   info<<" #vx3-functions = "<<(int)vx3BCs.size()<<endl;
   info<<" protected variables: x1, x2, x3, t"<<endl;
   
   const vector<BCFunction>* bcvecs[3] = { &vx1BCs, &vx2BCs, &vx3BCs };
   for(int i=0; i<3; i++)
   {
      for(size_t pos=0; pos<bcvecs[i]->size(); ++pos)
      {
         info<<"\n   vx"<<(i+1)<<"-function nr."<<pos<<":"<<endl;
         info<<(*bcvecs[i])[pos]<<endl;
      }
   }
   return info.str();
}


