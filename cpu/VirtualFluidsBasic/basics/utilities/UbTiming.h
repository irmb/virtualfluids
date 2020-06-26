//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBTIMING_H
#define UBTIMING_H

#include <string>
#include <limits>
#include <iostream>
#include <sstream>
#include <vector>
#include <ctime>

#ifdef CAB_RCF
   #include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif //CAB_RCF

#ifdef VF_MPI
   #include <mpi.h>
   #include <basics/parallel/PbMpi.h>
#endif //VF_MPI

/*=========================================================================*/
//  UbTiming - Time Measuring                                              
//                                                                         
//
//
//This Class provides the base for ...
//<BR><BR>
//@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
//@author <A HREF="mailto:geller@cab.bau.tu-bs.de">S. Geller</A>
//@version 1.1 - 14.02.06
// 

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
   virtual ~UbTiming() {}  
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

/*=========================================================================*/
//  UbTimer - Time Measuring                                              
//                                                                         
//
//
//This Class provides the base for ...
//<BR><BR>
//@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
//@version 1.0 - 16.08.2007
// 
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

//example:
//t=0  start 
//t=1 
//t=2  stop  -> return 2; getLapTime=2; getTotalTime 2; getLapTimes:  2
//t=3 
//t=4 
//t=5  stop  -> return 3; getLapTime=3; getTotalTime 5; getLapTimes:  2,3
//t=6  stop  -> return 1; getLapTime=1; getTotalTime 6; getLapTimes:  2,3,1
//t=7  
//t=8  start ->no consideration of time 7 and 8 
//t=9  
//t=10 stop  -> return 2; getLapTime=2; getTotalTime 8; getLapTimes:  2,3,1,2
//t=11 resetAndStart timer wird zurueckgestellt und neu gestaret
//t=12
//t=13 
//t=14 stop  -> return 3; getLapTime=3; getTotalTime 3; getLapTimes:  3

class UbTimer
{
public:
   UbTimer(const bool& storeLapTimes = false) 
      :  name("unamed"), isMeasuring(false), storeLapTimes(storeLapTimes)
       , startTime(0.0), totalTime(0.0), lapTime(0.0)
   {

   }
   /*==========================================================*/
   UbTimer(const std::string& name, const bool& storeLapTimes = false) 
      :  name(name), isMeasuring(false), storeLapTimes(storeLapTimes)
       , startTime(0.0), totalTime(0.0), lapTime(0.0)
   {

   }
   /*==========================================================*/
   virtual ~UbTimer() {}  
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


#ifdef CAB_RCF
   template<class Archive>
   void SF_SERIALIZE(Archive & ar)
   {
      ar & name;
      ar & isMeasuring;
      ar & startTime;
      ar & totalTime;
      ar & lapTime;
      ar & lapTimes;
      ar & storeLapTimes;
   }
#endif //CAB_RCF

protected:
   std::string name;
   bool        isMeasuring;
   bool        storeLapTimes;

   double      startTime;
   double      totalTime;
   double      lapTime;
   
   std::vector<double> lapTimes;
};


/*=========================================================================*/
//  UbProgressTimer - Time Measuring                                              
//                                                                         
//
//
//UbProressTimer misst die Zeit von seiner Instantiierung bis zur Zerst�rung
//und gib die verstrichene Zeit auf "os" in [s] aus
//example:
// {
//    UbProgressTimer timer;
//    UbSystem::sleepS(10);
// } //--> 10s
//<BR><BR>
//@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
//@version 1.0 - 10.03.2008
// 
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
  ~UbProgressTimer()
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