/*
*  BlocksPostprocessor.h
*
*  Created on: 24.09.2012
*  Author: K. Kucher
*/

#ifndef BlocksPostprocessor_H_
#define BlocksPostprocessor_H_

#include "CoProcessor.h"
#include "Communicator.h"
#include "WbWriter.h"

#include <boost/shared_ptr.hpp>
class WriteBlocksCoProcessor;
typedef boost::shared_ptr<WriteBlocksCoProcessor> WriteBlocksCoProcessorPtr;

class WriteBlocksCoProcessor: public CoProcessor {
public:
   WriteBlocksCoProcessor(Grid3DPtr grid, UbSchedulerPtr s, const std::string& path, WbWriter* const writer, CommunicatorPtr comm);
   virtual ~WriteBlocksCoProcessor();
   void process(double step);
protected:
   void collectData(double step);
   std::string path;
   WbWriter* writer;
   CommunicatorPtr comm;
};


#endif 
