//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef TOTRANSMITTERLOCAL_H
#define TOTRANSMITTERLOCAL_H

#include <basics/utilities/UbException.h>
#include <basics/transmitter/TbTransmitter.h>

#include <boost/shared_ptr.hpp>

/*================================================================================*/
/*   TbLocalTransmitter, TbVectorSenderLocal, TbVectorReceiverLocal               */
/*                                                                                */
/**
This Class provides the base for exception handling.
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 08.11.07
*/ 

/*
usage: ...
*/

//////////////////////////////////////////////////////////////////////////
// LocalTransmitter lokalen Datenaustausch
// data = send- und zugleich receive-buffer
template<typename T>
class TbLocalTransmitter : public TbTransmitter<T>
{
public:
   typedef boost::shared_ptr< TbLocalTransmitter<T> > TbLocalTransmitterPtr;

   typedef T value_type;

public:
   TbLocalTransmitter() : TbTransmitter<T>() 
   {

   }
   
   bool isLocalTransmitter()  const { return true;                         }
   bool isRemoteTransmitter() const { return !this->isLocalTransmitter();  }

   //send buffer wird autom resized
   void sendDataSize()    { }
   //reiceive braucht nichts machen, da send==receive buffer ;-)
   void receiveDataSize() { } 

   void        sendData()    { }
   value_type& receiveData() { return this->data; }

   std::string toString()  const { return "TbLocalTransmitter"+(std::string)typeid(T).name(); }
};

//////////////////////////////////////////////////////////////////////////
// TbVectorSender/ReceiverLocal lokalen Datenaustausch ueber ZWEI vektoren
template<typename T>
class TbVectorReceiverLocal : public TbTransmitter<T>
{
public:
   typedef T value_type;

public:
   TbVectorReceiverLocal() : TbTransmitter<value_type>() 
   {

   }
   //virtual ~TbVectorReceiverLocal() { std::cout<<typeid(*this).name()<<" tot"<<std::endl;   }

   bool isLocalTransmitter()  const { return true;                         }
   bool isRemoteTransmitter() const { return !this->isLocalTransmitter();  }

   //send buffer wird autom resized
   void sendDataSize()    { UB_THROW( UbException(UB_EXARGS,"empfaengt nur") ); }
   //reiceive braucht nichts machen, das macht der sender :-)
   void receiveDataSize() { } 

   void         sendData()    { UB_THROW( UbException(UB_EXARGS,"empfaengt nur") ); }
   value_type&  receiveData() { return this->data; }

   std::string toString() const { return "TbVectorReceiverLocal<"+(std::string)typeid(T).name()+">"; }
};

template<typename T>
class TbVectorSenderLocal : public TbTransmitter<T>
{
public:
   typedef T value_type;

public:
   TbVectorSenderLocal(boost::shared_ptr< TbVectorReceiverLocal< value_type > > receiver) 
      : TbTransmitter< value_type >(), receiver(receiver) 
   {

   }
   //virtual ~TbVectorSenderLocal() { std::cout<<typeid(*this).name()<<" tot"<<std::endl;   }

   bool isLocalTransmitter()  const { return true;                         }
   bool isRemoteTransmitter() const { return !this->isLocalTransmitter();  }

   //send buffer wird autom resized
   void sendDataSize()  
   { 
      assert(receiver!=NULL); 
      receiver->getData().resize( this->data.size() ); 
   }
   //reiceive braucht nichts machen, da send==receive buffer ;-)
   void receiveDataSize()  { UB_THROW( UbException(UB_EXARGS,"sendet nur") ); } 
   
   void sendData()    
   { 
      assert( this->data.size() == receiver->getData().size() );
      receiver->getData() = this->data;
//       for(int i=(int)this->data.size()-1; i>=0; --i)
//          receiver->getData()[i]= this->data[i];
   }
   value_type& receiveData() { UB_THROW( UbException(UB_EXARGS,"sendet nur") ); }

   std::string toString() const { return "TbVectorSenderLocal<"+(std::string)typeid(T).name()+">"; }

protected:
   boost::shared_ptr< TbVectorReceiverLocal< value_type > > receiver; 
};
                                        
#endif //TOTRANSMITTERLOCAL_H 
