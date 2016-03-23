#ifndef CREATETRANSMITTERSHELPER_H 
#define CREATETRANSMITTERSHELPER_H

#include "Block3D.h"
#include "Communicator.h"

#include <basics/transmitter/TbTransmitter.h>
#include <basics/transmitter/TbTransmitterMpiPool.h>
#include <basics/container/CbVector.h>

//! \brief The class helps to create Transmitters. 
//! \details It is created two types of Transmitters: MPI and BOND
//! \author K. Kucher
class CreateTransmittersHelper
{
public:
   //! Switch between same level and interpolation Connectors. NONE - for same level Connectors; SW, NW, NE, SE - source/target fine blocks in grid interface   
   enum IBlock {NONE, SW, NW, NE, SE};
   //! Switch between MPI and BOND Transmitters
   enum TransmitterType {MPI, BOND, MPI2BOND};
public:
   typedef CbVector <LBMReal> DataType;
   typedef boost::shared_ptr< TbTransmitter< DataType > > TransmitterPtr;
public:
   CreateTransmittersHelper();
   void createTransmitters(const Block3DPtr sblock, const Block3DPtr tblock, int dir, IBlock ib,
                                 TransmitterPtr& sender, TransmitterPtr& receiver, CommunicatorPtr comm, TransmitterType tType);
protected:
private:
   unsigned int generatePoolKey(int srcRank, int srcLevel, int tgtRank, int tgtLevel);
   unsigned int  generateVectorKey(int x1, int x2, int x3,/*int id,*/ int dir, IBlock ib);
   int generateMPITag(int srcLevel, int tgtLevel);
   static unsigned int vKey;
};

#endif
