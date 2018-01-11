//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef D3Q27BCFUNCTION_H
#define D3Q27BCFUNCTION_H

#include <basics/utilities/UbInfinity.h>

#include <MuParser/include/muParser.h>


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
