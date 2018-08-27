#ifndef UBNUPSTIMER_H
#define UBNUPSTIMER_H

#include <basics/utilities/UbTiming.h>
#include <sstream>
#include <vector>


/*=========================================================================*/
/*  UbNupsTimer                                                             */
/*                                                                         */
/**
This Class provides the base for ...
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 01.11.04
*/ 
class UbNupsTimer : public UbTiming
{
public:
   UbNupsTimer() : UbTiming()
   {
      mTempNodes = 0.0;
      mNofNodes.resize(0);
      mDurations.resize(0);
   }
   /*==========================================================*/
   UbNupsTimer(std::string name) : UbTiming(name)
   {
      mNofNodes.resize(0);
      mDurations.resize(0);
      mTempNodes = 0.0;
   }
   /*==========================================================*/
   void initTiming()
   {
      UbTiming::initTiming();
      mNofNodes.resize(0);
      mDurations.resize(0);
      mTempNodes   = 0.0;
   }
   /*==========================================================*/
   void startTiming(double nofNodes)
   {
      mTempNodes=nofNodes;
      UbTiming::startTiming();
   }
   /*==========================================================*/
   void endTiming()
   {
      UbTiming::endTiming();
      //save #node and time informations
      mNofNodes.push_back(mTempNodes);
      mDurations.push_back(UbTiming::getDuration());
      //reset internal timecounter
      UbTiming::initTiming();
   }
   /*==========================================================*/
   double getAverageNups()
   { 
      double averageNups = 0.0;
      for(int i=0;i<(int)mNofNodes.size();i++)
         averageNups+=mNofNodes.at(i)/mDurations.at(i);
      
      return averageNups/(double)mNofNodes.size(); 
   }
   /*==========================================================*/
   double getSumOfDuration()
   {
      double duration = 0.0;
      for(int i=0;i<(int)mDurations.size();i++) duration+=mDurations.at(i);
      return duration;
   }
   /*==========================================================*/
   std::string getNupsString()
   {
      std::stringstream ss;
      ss<<"saved nups informations"<<std::endl;
      for(int i=0;i<(int)mNofNodes.size();i++)
         ss<<mNofNodes.at(i)<<"nodes/"<<mDurations.at(i)<<"sec="<<mNofNodes.at(i)/mDurations.at(i)<<"nups\n";
      return ss.str();
   }

protected:

private:
   std::vector<double> mNofNodes;
   std::vector<double> mDurations;
   
   double mTempNodes;
};

#endif
