//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef BLOCK2DCONNECTOR_H
#define BLOCK2DCONNECTOR_H

#include <vector>
#include <string>

#include <basics/utilities/UbTuple.h>

#include <boost/shared_ptr.hpp>

class Block3DConnector;
typedef boost::shared_ptr<Block3DConnector> Block3DConnectorPtr;

class Block3DConnector
{
   //FETOL
#ifdef VF_FETOL
public:
   enum TransmitterType{NONE=0, MPI=1, BOND=2};
#endif
public:
   Block3DConnector() 
      : sendDir(-1)
      , invStep(false)
#ifdef VF_FETOL
      , ttype(NONE)
#endif
   {}

   Block3DConnector(const int& sendDir) 
      : sendDir(sendDir)
      , invStep(false)
#ifdef VF_FETOL
      , ttype(NONE)
#endif
   {}

   virtual ~Block3DConnector() {}

   virtual void init()=0;

   //for synchronize the send- and receive-bufferlength
   virtual void sendTransmitterDataSize()=0;  
   virtual void receiveTransmitterDataSize()=0;

   //send operations (should be called in given order!!!)
   virtual void prepareForSend()=0;
   virtual void fillSendVectors()=0;
   virtual void sendVectors()=0;
   
   //receive operations (should be called in given order!!!)
   virtual void prepareForReceive()=0;
   virtual void receiveVectors()=0;
   virtual void distributeReceiveVectors()=0;

   //info section
   virtual bool isLocalConnector()  = 0;
   virtual bool isRemoteConnector() = 0;
   virtual bool isInterpolationConnectorCF() = 0;
   virtual bool isInterpolationConnectorFC() = 0;

   //grid refinement
   virtual void setInvStep(bool step) {invStep = step;}
   virtual int getSendDir() const { return sendDir; } 

   //virtual double getSendRecieveTime() = 0;
   
   virtual void prepareForSendX1() = 0;
   virtual void prepareForSendX2() = 0;
   virtual void prepareForSendX3() = 0;

   virtual void sendVectorsX1() = 0;
   virtual void sendVectorsX2() = 0;
   virtual void sendVectorsX3() = 0;

   virtual void prepareForReceiveX1() = 0;
   virtual void prepareForReceiveX2() = 0;
   virtual void prepareForReceiveX3() = 0;

   virtual void receiveVectorsX1() = 0;
   virtual void receiveVectorsX2() = 0;
   virtual void receiveVectorsX3() = 0;

   //FETOL
#ifdef VF_FETOL
   void setTransmitterType(TransmitterType ttype) { this->ttype=ttype;}
   TransmitterType getTransmitterType() {return ttype;}
#endif

protected:
   int  sendDir;
   bool invStep;
   //FETOL 
#ifdef VF_FETOL
   TransmitterType ttype;
#endif
};

#endif //BLOCK2DCONNECTOR_H
