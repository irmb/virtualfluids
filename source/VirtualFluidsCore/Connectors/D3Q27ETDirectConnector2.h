#ifndef D3Q27ETDirectConnector2_H
#define D3Q27ETDirectConnector2_H

#include <boost/weak_ptr.hpp>

#include "LocalBlock3DConnector.h"
#include "Block3D.h"

class D3Q27ETDirectConnector2 : public LocalBlock3DConnector
{
public:
   D3Q27ETDirectConnector2(Block3DPtr from, Block3DPtr to, const int& sendDir) 
      :  LocalBlock3DConnector(from, to, sendDir)
   {

   }
   void sendVectors();
};

#endif //D3Q27ETDirectConnector2_H

