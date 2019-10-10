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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \   
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
//! \file BCFunction.h
//! \ingroup BoundarConditions
//! \author Sören Freudiger
//=======================================================================================

#ifndef D3Q27BCFUNCTION_H
#define D3Q27BCFUNCTION_H

#include <basics/utilities/UbInfinity.h>

#include <MuParser/include/muParser.h>

//! A class implements function parcer for boundary conditions  
class BCFunction
{
public:
   static const double INFTIMEDEPENDENT;
   static const double INFCONST;

public:
   BCFunction() 
      : starttime(-Ub::inf ), endtime(-Ub::inf ) 
   {

   }
   BCFunction( const mu::Parser& function, const double& starttime, const double& endtime )
      : function(function), starttime(starttime), endtime(endtime)
   {

   }
   BCFunction( const std::string& functionstring, const double& starttime, const double& endtime )
      : starttime(starttime), endtime(endtime)
   {
      this->setFunction(functionstring); 
   }
   BCFunction( const double& velocity, const double& starttime, const double& endtime )
      : starttime(starttime), endtime(endtime)
   {
      this->setFunction(velocity); 
   }

   void setFunction(const mu::Parser& function) { this->function = function; }
   void setFunction(const std::string& functionstring) { this->function.SetExpr(functionstring); }
   void setFunction(const double& constVelocity) { std::stringstream dummy; dummy<<constVelocity; function.SetExpr(dummy.str());  }
   void setStartTime(const double& starttime) {this->starttime = starttime; }
   void setEndTime(const double& starttime) {this->endtime = endtime; }

   mu::Parser&        getFunction()        { return function;  }
   const mu::Parser&  getFunction()  const { return function;  }
   const double&      getStartTime() const { return starttime; }
   const double&      getEndTime()   const { return endtime;   }

   std::string toString() const
   {
      std::stringstream info;
      if     (starttime==INFTIMEDEPENDENT) info<<"start=inf. timedep., ";
      else if(starttime==INFCONST        ) info<<"start=inf. const., ";
      else                                 info<<"start="<<starttime<<", ";
      if     (endtime==INFTIMEDEPENDENT) info<<"end=inf. timedep."<<std::endl;
      else if(endtime==INFCONST        ) info<<"end=inf. const."<<std::endl;
      else                               info<<"end="<<endtime<<std::endl;
      info<<"expr="<<function.GetExpr()<<std::endl;
      info<<"with constants: ";
      mu::valmap_type cmap = function.GetConst();
      for(mu::valmap_type::const_iterator item = cmap.begin(); item!=cmap.end(); ++item)
         info<<item->first<<"="<<item->second<<", ";
      return info.str();
   }
   /*==========================================================*/
   friend inline std::ostream& operator << (std::ostream& os, const BCFunction& bc) 
   {
      os<<bc.toString();
      return os;
   }

protected:
   mu::Parser function;
   double starttime;
   double endtime;

private:

};

#endif //D3Q27BCFUNCTION_H
