#ifndef CREATETRANSMITTERSHELPER_H 
#define CREATETRANSMITTERSHELPER_H

#include "Block3D.h"
#include "Communicator.h"

#include "LBMSystem.h"

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
   typedef SPtr< TbTransmitter< DataType > > TransmitterPtr;
public:
   CreateTransmittersHelper();
   void createTransmitters(const SPtr<Block3D> sblock, const SPtr<Block3D> tblock, int dir, IBlock ib,
                                 TransmitterPtr& sender, TransmitterPtr& receiver, SPtr<Communicator> comm, TransmitterType tType);
protected:
private:
   std::string generatePoolKey(int srcRank, int srcLevel, int tgtRank, int tgtLevel);
   std::string  generateVectorKey(int x1, int x2, int x3,/*int id,*/ int dir, IBlock ib);
   int generateMPITag(int srcLevel, int tgtLevel);
   static unsigned int vKey;
};

#endif
