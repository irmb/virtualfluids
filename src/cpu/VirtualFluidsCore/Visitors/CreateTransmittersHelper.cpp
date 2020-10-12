#include "CreateTransmittersHelper.h"
#include <D3Q27System.h>
#include <Communicator.h>
#include <string>

#ifdef VF_FETOL
   #include <FETOLTransmitterBondPool.h>
#endif
#include <MathUtil.hpp>

unsigned CreateTransmittersHelper::vKey = 0;

using namespace std;

//////////////////////////////////////////////////////////////////////////
CreateTransmittersHelper::CreateTransmittersHelper()
= default;
//////////////////////////////////////////////////////////////////////////
void CreateTransmittersHelper::createTransmitters(SPtr<Block3D> sblock, SPtr<Block3D> tblock, int dir, IBlock ib,
                                                        TransmitterPtr& sender, TransmitterPtr& receiver, SPtr<Communicator> comm, TransmitterType tType)
{
   //SourceBlock
   int srcLevel = sblock->getLevel();
//   int srcID    = sblock->getGlobalID();
 
   //TargetBlock 
   int tgtLevel = tblock->getLevel();
//   int tgtID    = tblock->getGlobalID();

   int invDir = D3Q27System::INVDIR[dir];

   if( srcLevel != tgtLevel ) invDir = dir;

   int srcRank = 0;
   int tgtRank = 0;

   if (tType == MPI)
   {
      srcRank  = sblock->getRank();
      tgtRank  = tblock->getRank();
   } 
#ifdef VF_FETOL
   else if (tType == MPI2BOND)
   {
      srcRank  = sblock->getLocalRank();
      tgtRank  = tblock->getLocalRank();
   }
#endif

   if (tType == MPI 
#ifdef VF_FETOL
      || tType == MPI2BOND
#endif
      )
   {
      string sendPoolKey = generatePoolKey(srcRank, srcLevel, tgtRank, tgtLevel);
      string receivePoolKey = generatePoolKey(tgtRank, tgtLevel, srcRank, srcLevel);

      TbCbVectorMpiPool <LBMReal>::MpiPoolPtr sendPool = TbCbVectorMpiPool <LBMReal>::getTbCbVectorMpiPool(sendPoolKey   );
      TbCbVectorMpiPool <LBMReal>::MpiPoolPtr recvPool = TbCbVectorMpiPool <LBMReal>::getTbCbVectorMpiPool(receivePoolKey);

      MPI_Comm mpi_comm = *((MPI_Comm*) comm->getNativeCommunicator());

      if( !sendPool ) sendPool = TbCbVectorMpiPool <LBMReal>::createTbCbVectorMpiPool(sendPoolKey   ,tgtRank, generateMPITag(srcLevel, tgtLevel), mpi_comm);
      if( !recvPool ) recvPool = TbCbVectorMpiPool <LBMReal>::createTbCbVectorMpiPool(receivePoolKey,tgtRank, generateMPITag(tgtLevel, srcLevel), mpi_comm);

      TbCbVectorMpiPool <LBMReal>::CbVectorKey keyOfSendCbVectorKey = generateVectorKey(sblock->getX1(), sblock->getX2(), sblock->getX3()/*tgtID*/, dir, ib);
      TbCbVectorMpiPool <LBMReal>::CbVectorKey keyOfRecvCbVectorKey = generateVectorKey(tblock->getX1(), tblock->getX2(), tblock->getX3()/*srcID*/, invDir, ib);

      ////////////////////////////////////////////////////////
      //DEBUG
      //int myid = comm->getProcessID();
      //FILE * file;
      ////char * name = "d:/temp/sendPoolKey.csv";
      //std::string name = "d:/temp/VectorKey" + UbSystem::toString(myid) + ".csv";
      //file = fopen(name.c_str(), "a");
      //fprintf(file, "%d;%d%;%d;%d;%d;%u;%d;%d%;%d;%d;%d;%u\n", sblock->getX1(), sblock->getX2(), sblock->getX3()/*tgtID*/, dir, ib, keyOfSendCbVectorKey, tblock->getX1(), tblock->getX2(), tblock->getX3()/*srcID*/, invDir, ib, keyOfRecvCbVectorKey);
      //fclose(file);
      ////////////////////////////////////////////////////////

      //create sender-/receiver
      sender   = TransmitterPtr( new TbCbVectorSenderMpiPool< LBMReal >(keyOfSendCbVectorKey,sendPool.get()) );
      receiver = TransmitterPtr( new TbCbVectorReceiverMpiPool< LBMReal >(keyOfRecvCbVectorKey,recvPool.get()) );
   }
#ifdef VF_FETOL
   if (tType == BOND)
   {
      int srcBondRank  = sblock->getRank();
      int tgtBondRank  = tblock->getRank();

      int sendBondPoolKey    = generatePoolKey(srcBondRank,srcLevel,tgtBondRank,tgtLevel);
      int receiveBondPoolKey = generatePoolKey(tgtBondRank,tgtLevel,srcBondRank,srcLevel);

      TbCbVectorBondPool <LBMReal>::BondPoolPtr sendPool = TbCbVectorBondPool <LBMReal>::getTbCbVectorBondPool(sendBondPoolKey   );
      TbCbVectorBondPool <LBMReal>::BondPoolPtr recvPool = TbCbVectorBondPool <LBMReal>::getTbCbVectorBondPool(receiveBondPoolKey);

      if( !sendPool ) sendPool = TbCbVectorBondPool <LBMReal>::createTbCbVectorBondPool(sendBondPoolKey   ,tgtBondRank, generateMPITag(srcLevel, tgtLevel));
      if( !recvPool ) recvPool = TbCbVectorBondPool <LBMReal>::createTbCbVectorBondPool(receiveBondPoolKey,tgtBondRank, generateMPITag(tgtLevel, srcLevel));

      TbCbVectorBondPool <LBMReal>::CbVectorKey keyOfSendCbVectorKey = generateVectorKey(tgtID, dir, ib);     
      TbCbVectorBondPool <LBMReal>::CbVectorKey keyOfRecvCbVectorKey = generateVectorKey(srcID, invDir, ib);  

      //create sender-/receiver 
      sender   = TransmitterPtr( new TbCbVectorSenderBondPool< LBMReal >(keyOfSendCbVectorKey,sendPool.get()) );
      receiver = TransmitterPtr( new TbCbVectorReceiverBondPool< LBMReal >(keyOfRecvCbVectorKey,recvPool.get()) );
   }
#endif
}
//////////////////////////////////////////////////////////////////////////
int CreateTransmittersHelper::generateMPITag(int srcLevel, int tgtLevel)
{
   //The MPI standard guarantees that integers 0-32767 can be used as tags
   if (srcLevel == tgtLevel)
   {
      return srcLevel;
   }
   else
   {
      srcLevel++;
      tgtLevel++;
      std::string str = UbSystem::toString<int>(srcLevel) + UbSystem::toString<int>(tgtLevel);
      int r = UbSystem::stringTo<int>(str);
      return r;
   }
}
//////////////////////////////////////////////////////////////////////////
string CreateTransmittersHelper::generatePoolKey(int srcRank, int srcLevel, int tgtRank, int tgtLevel)
{
   std::string str;
   str = UbSystem::toString<int>(srcLevel);
   str += "#";
   str += UbSystem::toString<int>(tgtLevel);
   str += "#";
   str += UbSystem::toString<int>(srcRank);
   str += "#";
   str += UbSystem::toString<int>(tgtRank);

   return str;
}
//////////////////////////////////////////////////////////////////////////
string CreateTransmittersHelper::generateVectorKey(int x1, int x2, int x3, int dir, IBlock ib)
{
   std::string str;

   str += UbSystem::toString<int>(x1);
   str += "#";
   str += UbSystem::toString<int>(x2);
   str += "#";
   str += UbSystem::toString<int>(x3);
   str += "#";
   str += UbSystem::toString<int>(dir);
   str += "#";
   str += UbSystem::toString<int>(ib);

   return str;
}


