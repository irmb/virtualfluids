/*
*  WritePeBlocksCoProcessor.h
*
*  Created on: 07.09.2018
*  Author: K. Kutscher
*/

#ifndef WritePeBlocksCoProcessor_H_
#define WritePeBlocksCoProcessor_H_

#include <PointerDefinitions.h>
#include <string>

#include "CoProcessor.h"

#include <pe/basic.h>

class Communicator;
class Grid3D;
class UbScheduler;
class WbWriter;

class WritePeBlocksCoProcessor : public CoProcessor
{
public:
   WritePeBlocksCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string& path, WbWriter* const writer, SPtr<Communicator> comm, SPtr<walberla::blockforest::BlockForest> forest);
   virtual ~WritePeBlocksCoProcessor();

   void process(double step) override;

protected:
   void collectData(double step);

   std::string path;
   WbWriter* writer;
   SPtr<Communicator>  comm;
   SPtr<walberla::blockforest::BlockForest> forest;
};



#endif 

