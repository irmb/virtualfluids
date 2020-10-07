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
//! \file UbTiming.h
//! \ingroup utilities
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef UBTIMING_H
#define UBTIMING_H

#include <string>
#include <limits>
#include <iostream>
#include <sstream>
#include <vector>
#include <ctime>

#ifdef VF_MPI
   #include <mpi.h>
   #include <basics/parallel/PbMpi.h>
#endif //VF_MPI

class UbTiming
{
public:
	UbTiming()
   {
      this->duration		= 0.0;
      this->deltaT		= 0.0;
      this->startTime	= 0;
      this->name        = "noname";
   }
   /*==========================================================*/
   UbTiming(const std::string& name)
   {
      this->duration		= 0.0;
      this->deltaT		= 0.0;
      this->startTime	= 0;
      this->name        = name;
   }
   /*==========================================================*/
   virtual ~UbTiming() = default;  
   /*==========================================================*/
   virtual void initTiming()
   {
      this->duration = 0.0;	
   }
   /*==========================================================*/
   virtual void startTiming()
   {
      #if defined(VF_MPI) && !defined(CAB_RUBY)
         this->startTime = PbMpi::Wtime();
      #else
         this->startTime = (double)clock();	
      #endif //VF_MPI 
   }
   /*==========================================================*/
   virtual void initAndStartTiming()
   {
      this->initTiming();
      this->startTiming();
   }
   /*==========================================================*/
   virtual void endTiming()
   {
      this->stopTiming();
   }
   /*==========================================================*/
   virtual void stopTiming()
   {
      #if defined(VF_MPI) && !defined(CAB_RUBY)
            this->deltaT   = PbMpi::Wtime()-this->startTime;
      #else
         this->deltaT   = ((double)clock()-this->startTime)/(double)CLOCKS_PER_SEC;
      #endif //VF_MPI 

      this->duration += this->deltaT;
   }
   /*==========================================================*/
   virtual double getDuration() const
   {
      return this->duration;
   }
   /*==========================================================*/
   virtual void setName(const std::string& name)
   {
      this->name = name;
   }
   /*==========================================================*/
   virtual std::string getName() const
   { 
      return this->name; 
   }
   /*==========================================================*/
   void start()
   {
      this->duration = 0.0;

      #if defined(VF_MPI) && !defined(CAB_RUBY)
         this->startTime = PbMpi::Wtime();
      #else
         this->startTime = (double)clock();
      #endif //VF_MPI 
   }
   /*==========================================================*/
   void pause()
   {
      #if defined(VF_MPI) && !defined(CAB_RUBY)
         this->duration += PbMpi::Wtime()-this->startTime;
      #else
         this->duration +=((double)clock()-this->startTime)/(double)CLOCKS_PER_SEC;
      #endif //VF_MPI 
   }
   /*==========================================================*/
   void unpause()
   {
      #if defined(VF_MPI) && !defined(CAB_RUBY)
         this->startTime   = PbMpi::Wtime();
      #else
         this->startTime = (double)clock();
      #endif //VF_MPI 
   }
   /*==========================================================*/
   void stop()
   {
      #if defined(VF_MPI) && !defined(CAB_RUBY)
         this->duration += PbMpi::Wtime()-this->startTime;
      #else
         this->duration +=((double)clock()-this->startTime)/(double)CLOCKS_PER_SEC;
      #endif //VF_MPI 
   }
   /*==========================================================*/
   double getTicks() const            
   { 
      #if defined(VF_MPI) && !defined(CAB_RUBY)
         return PbMpi::Wtick();
      #else
         return double(1.0)/double(CLOCKS_PER_SEC);
      #endif  //VF_MPI 
   }

protected:
   std::string name;

   double startTime;
   double duration;
	double deltaT;
};

#include <basics/utilities/UbSystem.h> //for definitons of system/OS type

#ifdef UBSYSTEM_APPLE   //Apple hack
   #include <mach/mach_time.h>  
   #include <time.h>  
   #include <stdio.h> 
   inline void mach_absolute_difference(const uint64_t& end, const uint64_t& start, struct timespec *tp) 
   {  
         uint64_t difference = end - start;  
         static mach_timebase_info_data_t info = {0,0};  
   
         if (info.denom == 0)  
                 mach_timebase_info(&info);  
   
         uint64_t elapsednano = difference * (info.numer / info.denom);  
   
         tp->tv_sec = elapsednano * 1e-9;  
         tp->tv_nsec = elapsednano - (tp->tv_sec * 1e9);  
   } 
#elif defined(UBSYSTEM_LINUX) || defined(UBSYSTEM_AIX)
   #include <ctime>
   #include <unistd.h> // for sysconf
   #include <pthread.h>
#endif

/*=========================================================================*/
//! \brief Time Measuring                                              
//! \details                                                                         
//! example:
//! \code
//! t=0  start 
//! t=1 
//! t=2  stop  -> return 2; getLapTime=2; getTotalTime 2; getLapTimes:  2
//! t=3 
//! t=4 
//! t=5  stop  -> return 3; getLapTime=3; getTotalTime 5; getLapTimes:  2,3
//! t=6  stop  -> return 1; getLapTime=1; getTotalTime 6; getLapTimes:  2,3,1
//! t=7  
//! t=8  start ->no consideration of time 7 and 8 
//! t=9  
//! t=10 stop  -> return 2; getLapTime=2; getTotalTime 8; getLapTimes:  2,3,1,2
//! t=11 resetAndStart -> Timer is reset and restarted
//! t=12
//! t=13 
//! t=14 stop  -> return 3; getLapTime=3; getTotalTime 3; getLapTimes:  3
//! \endcode

class UbTimer
{
public:
   UbTimer(const bool& storeLapTimes = false) 
      :  name("unamed"),  storeLapTimes(storeLapTimes)
        
   {

   }
   /*==========================================================*/
   UbTimer(const std::string& name, const bool& storeLapTimes = false) 
      :  name(name), isMeasuring(false), storeLapTimes(storeLapTimes)
       , startTime(0.0), totalTime(0.0), lapTime(0.0)
   {

   }
   /*==========================================================*/
   virtual ~UbTimer() = default;  
   /*==========================================================*/
   double              getLapTime() const               { return this->lapTime;  }
   std::vector<double> getLapTimes() const              { return this->lapTimes; }
   void                setName(const std::string& name) { this->name = name;     }
   std::string         getName() const                  { return this->name;     }
   bool                isRunning() const                { return isMeasuring;    }
   bool                isStoringLapTimes() const        { return storeLapTimes;  }
   /*==========================================================*/
   void setStoreLapTimes(const bool& storeLapTimes) { this->storeLapTimes = storeLapTimes; }
   /*==========================================================*/
   void start()
   {
      this->isMeasuring = true;

      #if defined(VF_MPI) && !defined(CAB_RUBY)
          this->startTime = PbMpi::Wtime();
      #elif defined(UBSYSTEM_APPLE)
    	 this->startTime = mach_absolute_time();  
      #elif defined(UBSYSTEM_LINUX) || defined(UBSYSTEM_AIX)
         timespec tp;
         clock_gettime(CLOCK_REALTIME,&tp);
         this->startTime = (double)(tp.tv_sec)*1.0e9 + (double)(tp.tv_nsec);
      #else
         this->startTime = (double)clock();
      #endif //VF_MPI
   }
   /*==========================================================*/
   void resetAndStart() { this->reset(); this->start(); }
   /*==========================================================*/
   //stop: - stops the calculation and returns the time elapsed since last start/stop
   //      - timing continues
   double stop()
   {
      //if start() was never activated before:
      if(!isMeasuring) return 0.0; 
      
      #if defined(VF_MPI) && !defined(CAB_RUBY)
         double actTime = PbMpi::Wtime();
         this->lapTime  = actTime-this->startTime;
      #elif defined(UBSYSTEM_APPLE)
    	 double actTime = mach_absolute_time();  
         timespec tp;  
         mach_absolute_difference(actTime, this->startTime, &tp);
         this->lapTime  =  tp.tv_sec + tp.tv_nsec*1e-9;
	  #elif defined(UBSYSTEM_LINUX) || defined(UBSYSTEM_AIX)
         timespec tp;
         clock_gettime(CLOCK_REALTIME,&tp);
         double actTime = (double)(tp.tv_sec)*1.0e9 + (double)(tp.tv_nsec);
         this->lapTime  = (actTime-this->startTime)*1.0e-9;
      #else
         double actTime = (double)clock();
         this->lapTime  = (actTime-this->startTime)/(double)CLOCKS_PER_SEC;
      #endif //VF_MPI 
      
      this->startTime  = actTime;
      this->totalTime += this->lapTime;
      if(storeLapTimes) lapTimes.push_back(this->lapTime);

      return lapTime;
   }
   /*==========================================================*/
   void reset()
   {
      this->isMeasuring = false;
      
      this->startTime   = 0.0;
      this->totalTime   = 0.0;
      this->lapTime     = 0.0;

      lapTimes.resize(0);
   }
   /*==========================================================*/
   double getCurrentLapTime() const
   {
     //if start() was never activated before:
      if(!isMeasuring) return 0.0; 
      
      #if defined(VF_MPI) && !defined(CAB_RUBY)
         return PbMpi::Wtime() - this->startTime;
      #elif defined(UBSYSTEM_APPLE)
         timespec tp;  
         mach_absolute_difference(mach_absolute_time(), this->startTime, &tp);
         return tp.tv_sec + tp.tv_nsec*1e-9;
      #elif defined(UBSYSTEM_LINUX) || defined(UBSYSTEM_AIX)
         timespec tp;
         clock_gettime(CLOCK_REALTIME,&tp);
         return ((double)(tp.tv_sec)*1.0e9 + (double)(tp.tv_nsec) - this->startTime)*1.0e-9;
      #else
         return ( (double)clock() - this->startTime ) / (double)CLOCKS_PER_SEC;
      #endif //VF_MPI 
      
   }
   /*==========================================================*/
   double getTotalTime() const
   {
      return this->totalTime;
   }
   /*==========================================================*/
   std::string toString()
   {
      std::stringstream text;
      text<<*this;
      return text.str();
   }

   //ueberladene Operatoren
   /*==========================================================*/
   friend inline std::ostream& operator << (std::ostream& os, const UbTimer& timer) 
   {
       os<<"UbTimer[totalTime="<<timer.totalTime<<"sec, lapTimes(";
       for(std::size_t i=0; i<timer.lapTimes.size(); i++) os<<timer.lapTimes[i]<<",";
       os<<")]";
       return os;
   }


protected:
   std::string name;
   bool        isMeasuring{false};
   bool        storeLapTimes;

   double      startTime{0.0};
   double      totalTime{0.0};
   double      lapTime{0.0};
   
   std::vector<double> lapTimes;
};


/*=========================================================================*/
//! \brief Time Measuring                                              
//! 
//! \details UbProressTimer measures the time from its instantiation to destruction and spend the elapsed time on "os" in [s]
//! example:
//! \code
//!  {
//!     UbProgressTimer timer;
//!     UbSystem::sleepS(10);
//!  } //--> 10s
//! \endcode

class UbProgressTimer : public UbTimer
{
private:
	UbProgressTimer(const UbProgressTimer& rhs);
public:
  explicit UbProgressTimer( std::ostream & os = std::cout )
     : UbTimer(),os(os) 
  {
  	  this->start();
  }
  /*==========================================================*/
  ~UbProgressTimer() override
  {
  //  A) Throwing an exception from a destructor is a Bad Thing.
  //  B) The progress_timer destructor does output which may throw.
  //  C) A progress_timer is usually not critical to the application.
  //  Therefore, wrap the I/O in a try block, catch and ignore all exceptions.
    try
    {
      // use istream instead of ios_base to workaround GNU problem (Greg Chicares)
      std::istream::fmtflags old_flags = os.setf( std::istream::fixed,
                                                  std::istream::floatfield );
      std::streamsize old_prec = os.precision( 2 );
      os << stop() << " s" << std::endl;
      os.flags( old_flags );
      os.precision( old_prec );
    }
    catch (...) {} // eat any exceptions
  } 

private:
  std::ostream & os;
};


#endif //UBTIMING_H
