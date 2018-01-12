//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBSCHEDULER_H
#define UBSCHEDULER_H

#include <iostream>
#include <string>
#include <limits>
#include <cassert> 
#include <sstream>
#include <iomanip>
#include <algorithm>

#include <basics/utilities/UbSystem.h>
#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbInfinity.h>
#include <basics/utilities/UbComparators.h>
#include <basics/utilities/UbFileOutput.h>
#include <basics/utilities/UbFileInput.h>

/*=========================================================================*/
/*  UbScheduler                                                            */
/*                                                                         */
/**
namespace for global system-functions
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@author <A HREF="mailto:hegewald@cab.bau.tu-bs.de">J. Hegewald</A>
@version 1.0 - 06.09.06
@version 1.1 - 09.09.06
@version 1.2 - 03.07.08 - nun auch isDue(t) mehrmals fuer dasselbe t moeglich
                          isDue(t) auch fuer t < lastUsedT
                          bug entfernt, der bei Schedule (5,0,500) auch 505 als Due zurückgibt!
*/ 

/*
usage: ...
*/

// this class is not thread save
//

class UbScheduler
{
public:
   class UbSchedule
   {
      friend class UbScheduler;
   public:
      UbSchedule() :  step(Ub::inf), begin(Ub::inf), end(Ub::inf) { }
      UbSchedule(const double& step, const double& begin=0.0, const double& end=Ub::inf) 
         : step(step), begin(begin), end(end) 
      {  
      }
      double getStep()  const { return this->step;  }
      double getBegin() const { return this->begin; }
      double getEnd()   const { return this->end;   }
      
      /*==========================================================*/
      std::string toString() { std::stringstream text; text<<*this; return text.str(); }
      /*==========================================================*/
      friend inline std::ostream& operator << (std::ostream& os, const UbSchedule& schedule) 
      {
         os<<"Schedule[start,end,step]=["<<schedule.begin<<", "<<schedule.end<<", "<<schedule.step<<"]";
         return os;
      }

      //------------- implements CAB serialization ----- start
      virtual void write(UbFileOutput* out)
      {
         out->writeDouble( begin );
         out->writeDouble( end );
         out->writeDouble( step );
      }
      virtual void read(UbFileInput* in)
      {
         begin = in->readDouble();
         end   = in->readDouble();
         step  = in->readDouble();
      }
  

   private:
      double step, begin, end;
   };

public:
   UbScheduler() 
   {
      this->initVals();
   }
   /*==========================================================*/                         
   UbScheduler(const double& step,const double& begin=0, const double& end=Ub::inf ) 
   {
      this->initVals();
      this->addSchedule(step,begin,end);
   }
   /*==========================================================*/
   UbScheduler(const UbSchedule& schedule) 
   {
      this->initVals();
      this->addSchedule(schedule);
   }
   /*==========================================================*/
   virtual ~UbScheduler() {}
   /*==========================================================*/
   inline void addSchedule(const UbSchedule& schedule)
   {
      this->addSchedule(schedule.step, schedule.begin, schedule.end);
   }
   /*==========================================================*/
   bool addSchedule(const double& step, const double& begin, double end)
   {
      if( UbMath::zero(step) || begin>end )
      { 
         std::cerr<<"UbScheduler::addSchedule - invalid Schedule:\n\t"<<UbSchedule(step, begin, end)<<std::endl;
         return false; 
      }
      
      if( UbMath::less( end, (double)Ub::inf )  )
      {
         //es kann vorkommen, dass man mit dem intervall nicht genau auf den letzten wert kommt
         //(z.B. step=2; start=0; end=9; -> ende wird angepasst)
         //also wenn end-begin>Ub::inf ist, dann geht es halt nicht.. ein cast in long double half hier nichts
         double multiplier=0.0;
         double fractpart =  modf( (end-begin)/step, &multiplier);
         if( !UbMath::zero(fractpart) )
         {
            //tmp-speicherung (fuer cerr)
            fractpart = end;
            //neues ende
            end = begin+multiplier*step;
            
            std::cerr<<"Warning: UbScheduler::addSchedule - "
                      <<"end of schedule was adapted to intervall \n\t"
                      <<"from "<< UbSchedule(step, begin, fractpart) <<" to "<< UbSchedule(step, begin, end) <<std::endl;
         }
      }

      //nu aber:
      schedules.push_back(UbSchedule(step, begin, end));

      if( end>maxT ) maxT = end;

      double potentialDueTime;
      if(   calcNextDueTimeForSchedule(schedules.back(), lastUsedT, potentialDueTime)
         && potentialDueTime < nextDueTime   )
      {
         nextDueTime = potentialDueTime;
      }

      return true;
   }
   /*==========================================================*/
   //returns true if scheduler contains schedules
   bool   hasSchedules() const { return !schedules.empty(); }
   /*==========================================================*/
   //time bei dem das letzte mal isDue(time) true war
   double getLastDueTime() const { return lastDueTime; }
   /*==========================================================*/
   //time bei dem das naechste mal isDue(time) true ergibt
   double getNextDueTime() const { return nextDueTime; }
   /*==========================================================*/
   //maxDueTime (maxTime der Schedules!
   double getMaxDueTime()  const { return this->maxT; }
   /*==========================================================*/
   bool isDue(const double& t)
   {
      lastUsedT = t;
      if( UbMath::greaterEqual(t,nextDueTime) ) 
      {
         //groesser maxT is nicht
         if( UbMath::greater(t,maxT) )  return false;
         
         //temp var
         double actDueTime = nextDueTime;

         //um Suche nach nextDueTime bei "Zukunfts-t" zu optimieren, setzt man die "start"-suchzeit auf "t-1":
         nextDueTime = t-1; //t-1 deshlab, damit falls z.B. while Schleife nicht durchlaufen wird
                            //die folgende if Abfrage nicht faelschlicher Weise true ist!
         while( UbMath::greaterEqual(t,nextDueTime) && !UbMath::equal(nextDueTime, maxT) )
         {
            double tmpNextDueTime = maxT, potentialDueTime=-1.0;
            for(std::size_t i=0; i<schedules.size(); i++)
            {
               if(   calcNextDueTimeForSchedule(schedules[i], nextDueTime, potentialDueTime)
                  && potentialDueTime < tmpNextDueTime                 )
               {
                  assert( nextDueTime < potentialDueTime );
                  tmpNextDueTime = potentialDueTime;
               }
            }
            actDueTime  = nextDueTime;
            nextDueTime = tmpNextDueTime;
         } 

         //wenn t = der aktuuellen oder gar schon der nächstmöglichen ist (hierbei wurde
         //zuvor actDueTime und nextDueTime ggf. angepasst)
         //Bsp.: nextDuTime war 5, aber für t=400 gilt andere schedule -> Bsp actDue=350 und nextDue 405
         if(    UbMath::equal(t,actDueTime)    
             || UbMath::equal(t,nextDueTime) ) 
         {
            lastDueTime = t;
            return true;
         }
      }
      else if( UbMath::lessEqual(t, lastDueTime) ) 
      {
         if(UbMath::equal(t, lastDueTime) ) return true; //braucht man, wenn man für dasselbe t isDue(t) aufruft
         else  
         {
            //Fall: Zeit liegt faktisch in der Vergangenheit -> neu initialsisieren
            double tmpNextDueTime = maxT, potentialDueTime=-1.0;
            for(size_t i=0; i<schedules.size(); i++)
            {
               if(   calcNextDueTimeForSchedule(schedules[i], t-1, potentialDueTime)
                  && potentialDueTime < tmpNextDueTime                 )
               {
                  tmpNextDueTime = potentialDueTime;
               }
            }
            nextDueTime = tmpNextDueTime;

            return UbMath::equal(t, nextDueTime);
         }
      }

      return false;
   }
   /*==========================================================*/
   inline double getMinBegin( ) const
   {
      if( schedules.empty() ) return Ub::inf;
      return std::min_element(schedules.begin(), schedules.end(),UbComparators::membercomp(&UbSchedule::getBegin) )->getBegin();
   }
   /*==========================================================*/
   inline double getMaxBegin( ) const
   {
      if( schedules.empty() ) return Ub::inf;
      return std::max_element(schedules.begin(), schedules.end(),UbComparators::membercomp(&UbSchedule::getBegin) )->getBegin();
   }
   /*==========================================================*/
   inline double getMinEnd( ) const
   {
      if( schedules.empty() ) return Ub::inf;
      return std::min_element(schedules.begin(), schedules.end(),UbComparators::membercomp(&UbSchedule::getEnd) )->getEnd();
   }
   /*==========================================================*/
   inline double getMaxEnd( ) const
   {
      if( schedules.empty() ) return Ub::inf;
      return std::max_element(schedules.begin(), schedules.end(),UbComparators::membercomp(&UbSchedule::getEnd) )->getEnd();
   }
   /*==========================================================*/
   inline double getMinStep( ) const
   {
      if( schedules.empty() ) return Ub::inf;
      return std::min_element(schedules.begin(), schedules.end(),UbComparators::membercomp(&UbSchedule::getStep) )->getStep();
   }
   /*==========================================================*/
   inline double getMaxStep( ) const
   {
      if( schedules.empty() ) return Ub::inf;
      return std::max_element(schedules.begin(), schedules.end(),UbComparators::membercomp(&UbSchedule::getStep) )->getStep();
   }
   /*==========================================================*/
   inline std::string toString() const
   {
      std::stringstream text;
      text<<*this;
      return text.str();
   }
   /*==========================================================*/
   friend inline std::ostream& operator << (std::ostream& os, const UbScheduler& scheduler) 
   {
      os<<"UbScheduler\n";
      os<<"Schedule |       start       |        end        |     intervall     "<<std::endl;
      for(std::size_t i=0; i<scheduler.schedules.size(); i++)
         os<<std::setw(9)<<i<<"|"
           <<std::setw(19)<<scheduler.schedules[i].getBegin()<<"|"
           <<std::setw(19)<<scheduler.schedules[i].getEnd()  <<"|"
           <<std::setw(19)<<scheduler.schedules[i].getStep() <<std::endl;
      return os;
   }

   //------------- implements CAB serialization ----- start
   virtual void write(UbFileOutput* out)
   {
      out->writeSize_t( schedules.size() );
      
      for(std::size_t i=0; i<schedules.size(); i++)
         schedules[i].write(out);
   }
   virtual void read(UbFileInput* in)
   {
      this->initVals();

      std::size_t nofSchedules = in->readSize_t();
      for(std::size_t i=0; i<nofSchedules; i++)
      {
         UbSchedule schedule;
         schedule.read(in);
         this->addSchedule(schedule);
      }
   }

protected:
   /*==========================================================*/
   void initVals()
   {
      lastUsedT   = -Ub::inf; 
      lastDueTime = -Ub::inf;
      nextDueTime =  Ub::inf;
      maxT        = -Ub::inf;
   }
   /*==========================================================*/
   // calculates next due time for a schedule 
   // with  nextDueTime > searchStart
   bool calcNextDueTimeForSchedule(const UbSchedule& schedule, const double& searchStart, double& nextDueTime )
   {
      if     ( UbMath::greater(searchStart, schedule.end  ) ) return false;
      else if( UbMath::less(   searchStart, schedule.begin) ) nextDueTime = schedule.begin;
      else                            
      {
         nextDueTime = schedule.begin + ((int)((searchStart-schedule.begin)/schedule.step)+1)*schedule.step;
         if(   UbMath::less(   nextDueTime, searchStart )
            || UbMath::greater(nextDueTime, schedule.end) ) 
         {
            return false;
         }
      }
      return true;
   }

protected:
   double lastUsedT;
   double lastDueTime;
   double nextDueTime;
   double maxT;
   
   std::vector<UbSchedule> schedules;
};

typedef UbScheduler::UbSchedule UbSchedule;
// inline std::ostream& operator<<( std::ostream& os, const UbScheduler& scheduler )
// {
//    os<<"UbScheduler\n";
//    os<<"Schedule |       start       |        end        |     intervall     "<<std::endl;
//    for(std::size_t i=0; i<scheduler.schedules.size(); i++)
//       os<<std::setw(9)<<i<<"|"
//         <<std::setw(19)<<scheduler.schedules[i].getBegin()<<"|"
//         <<std::setw(19)<<scheduler.schedules[i].getEnd()  <<"|"
//         <<std::setw(19)<<scheduler.schedules[i].getStep() <<std::endl;
//    return os;
// }

#endif //UBSCHEDULER_H



//int main(int argc, char** argv)            
//{   
//	UbScheduler writeSchedule;
////	writeSchedule.addSchedule(0,2000,100);
////	writeSchedule.addSchedule(3005,4500,300);
////	writeSchedule.addSchedule(0,10,1);
////	writeSchedule.addSchedule(0,100001,100);
//	writeSchedule.addSchedule(0,2,1);
//	writeSchedule.addSchedule(0,100001,200);
//
//	for(int t = 0; t < 1001; t++)
//	{
//		if(writeSchedule.isDue(t))
//		{
//			cout<<"due@ "<<t<<endl;
//		}
//	}
//	return 0;
//}

