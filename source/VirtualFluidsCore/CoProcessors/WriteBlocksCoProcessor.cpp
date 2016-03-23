#include "WriteBlocksCoProcessor.h"
#include "basics/writer/WbWriterVtkXmlASCII.h"
#include <boost/foreach.hpp>
#include "D3Q27System.h"

WriteBlocksCoProcessor::WriteBlocksCoProcessor(Grid3DPtr grid, UbSchedulerPtr s, 
                                         const std::string& path, WbWriter* const writer, 
                                         CommunicatorPtr comm) :
                                         CoProcessor(grid, s),
                                         path(path),
                                         writer(writer),
                                         comm(comm)
{

}
//////////////////////////////////////////////////////////////////////////
WriteBlocksCoProcessor::~WriteBlocksCoProcessor() 
{
}
//////////////////////////////////////////////////////////////////////////
void WriteBlocksCoProcessor::process(double step)
{
   if(scheduler->isDue(step) )
      collectData(step);
}
//////////////////////////////////////////////////////////////////////////
void WriteBlocksCoProcessor::collectData(double step)
{
   if (comm->getProcessID() == comm->getRoot())
   {
      int istep = int(step);
      std::vector<std::string> filenames;
      std::vector< UbTupleFloat3 > nodes;
      std::vector< UbTupleInt8 > cells;
      std::vector<std::string> celldatanames;

      celldatanames.push_back("isActive");
      celldatanames.push_back("rank");
      celldatanames.push_back("interface");
      celldatanames.push_back("ID");
      celldatanames.push_back("part");
      celldatanames.push_back("level");
      //celldatanames.push_back("connectorCF");
      //celldatanames.push_back("connectorFC");
#if defined VF_FETOL
      celldatanames.push_back("bundle");
#endif

      std::vector< std::vector< double > > celldata(celldatanames.size());

      int nr=0;
      int minInitLevel = this->grid->getCoarsestInitializedLevel();
      int maxInitLevel = this->grid->getFinestInitializedLevel();

      for(int level = minInitLevel; level<=maxInitLevel;level++)
      {
         std::vector<Block3DPtr> blockVector;
         grid->getBlocks(level, blockVector);
         BOOST_FOREACH(Block3DPtr block, blockVector)
         {
            UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);
            UbTupleDouble3 blockLength = grid->getBlockLengths(block);

            nodes.push_back(makeUbTuple((float)(val<1>(org)), (float)(val<2>(org)), (float)(val<3>(org))));
            nodes.push_back(makeUbTuple((float)(val<1>(org)+val<1>(blockLength)), (float)(val<2>(org)), (float)(val<3>(org))));
            nodes.push_back(makeUbTuple((float)(val<1>(org)+val<1>(blockLength)), (float)(val<2>(org)+val<2>(blockLength)), (float)(val<3>(org))));
            nodes.push_back(makeUbTuple((float)(val<1>(org)), (float)(val<2>(org)+val<2>(blockLength)), (float)(val<3>(org))));
            nodes.push_back(makeUbTuple((float)(val<1>(org)), (float)(val<2>(org)), (float)(val<3>(org)+val<3>(blockLength))));
            nodes.push_back(makeUbTuple((float)(val<1>(org)+val<1>(blockLength)), (float)(val<2>(org)), (float)(val<3>(org)+val<3>(blockLength))));
            nodes.push_back(makeUbTuple((float)(val<1>(org)+val<1>(blockLength)), (float)(val<2>(org)+val<2>(blockLength)), (float)(val<3>(org)+val<3>(blockLength))));
            nodes.push_back(makeUbTuple((float)(val<1>(org)), (float)(val<2>(org)+val<2>(blockLength)), (float)(val<3>(org)+val<3>(blockLength))));
            cells.push_back(makeUbTuple(nr, nr+1, nr+2, nr+3, nr+4, nr+5, nr+6, nr+7));
            nr += 8;

            //data
            celldata[0].push_back((double)block->isActive());
            celldata[1].push_back((double)block->getRank());
            celldata[2].push_back((double)block->hasInterpolationFlag());
            celldata[3].push_back((double)block->getGlobalID());
            celldata[4].push_back((double)block->getPart());
            celldata[5].push_back((double)block->getLevel());

            //bool flag = false;
            //std::vector<Block3DConnectorPtr> connectors;

            //block->pushBackLocalInterpolationConnectorsCF(connectors);
            //for (std::size_t i = 0; i<connectors.size(); i++)
            //   if (connectors[i])
            //   {
            //      if (connectors[i]->getSendDir() == D3Q27System::BS)
            //      {

            //         flag = true;
            //      }
            //   }

            //if (flag)
            //{
            //   celldata[6].push_back(1);
            //   UBLOG(logINFO, "CF: "+block->toString());
            //}
            //else
            //{
            //   celldata[6].push_back(0);
            //}

            //flag = false;
            //connectors.resize(0);
            //block->pushBackLocalInterpolationConnectorsFC(connectors);
            //for (std::size_t i = 0; i<connectors.size(); i++)
            //   if (connectors[i])
            //   {
            //      if (connectors[i]->getSendDir() == D3Q27System::BS)
            //      {

            //         flag = true;
            //      }
            //   }

            //if (flag)
            //{
            //   celldata[7].push_back(1);
            //   UBLOG(logINFO, "FC: "+block->toString());
            //}
            //else
            //{
            //   celldata[7].push_back(0);
            //}

#ifdef VF_FETOL            
            celldata[6].push_back( (double)block->getBundle());
#endif
         }
      }

      filenames.push_back(writer->writeOctsWithCellData(path+"/blocks/blocks_" + UbSystem::toString(grid->getRank()) + "_" + UbSystem::toString(istep),nodes,cells,celldatanames,celldata));

      if (istep == CoProcessor::scheduler->getMinBegin())
      {
         WbWriterVtkXmlASCII::getInstance()->writeCollection(path+"/blocks/blocks_collection",filenames,istep,false);
      } 
      else
      {
         WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(path + "/blocks/blocks_collection", filenames, istep, false);
      }

      UBLOG(logINFO,"BlocksPostprocessor step: " << istep);
   }
}
