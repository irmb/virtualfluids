/**
* @file D3Q27ETFullDirectConnector.h
* @brief Connector send and receive full distribution in shared memory
*         
* @author Kostyantyn Kucher
* @date 08.06.2011
*/
#ifndef D3Q27ETFULLDIRECTCONNECTOR_H
#define D3Q27ETFULLDIRECTCONNECTOR_H

#include <boost/weak_ptr.hpp>

#include "LocalBlock3DConnector.h"
#include "Block3D.h"
#include "D3Q27System.h"

class D3Q27ETFullDirectConnector : public LocalBlock3DConnector
{
public:
   D3Q27ETFullDirectConnector(Block3DPtr from, Block3DPtr to, int sendDir);
   virtual ~D3Q27ETFullDirectConnector();
   void init();                       
   void sendVectors();

 protected:
   EsoTwist3DPtr  fFrom;
   EsoTwist3DPtr  fTo;
   int maxX1;
   int maxX2;
   int maxX3;
   LBMReal f[D3Q27System::ENDF+1];

   void fillData(EsoTwist3DPtr  fFrom, const int& x1, const int& x2, const int& x3);
   void distributeData(EsoTwist3DPtr  fTo, const int& x1, const int& x2, const int& x3);
};


#endif 

