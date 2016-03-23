/*
*  BlocksPostprocessor.h
*
*  Created on: 24.09.2012
*  Author: K. Kucher
*/

#ifndef BlocksPostprocessor_H_
#define BlocksPostprocessor_H_

#include "Postprocessor.h"
#include "Communicator.h"
#include "WbWriter.h"

#include <boost/shared_ptr.hpp>
class BlocksPostprocessor;
typedef boost::shared_ptr<BlocksPostprocessor> BlocksPostprocessorPtr;

class BlocksPostprocessor: public Postprocessor {
public:
   BlocksPostprocessor(Grid3DPtr grid, UbSchedulerPtr s, const std::string& path, WbWriter* const writer, CommunicatorPtr comm);
   virtual ~BlocksPostprocessor();
   void update(double step);
protected:
   void collectPostprocessData(double step);
   std::string path;
   WbWriter* writer;
   CommunicatorPtr comm;
};


#endif 
