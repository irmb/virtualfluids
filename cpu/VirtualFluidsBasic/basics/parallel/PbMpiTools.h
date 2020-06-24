//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef PBMPITOOLS_H
#define PBMPITOOLS_H

#include <vector>
#include <sstream>

//#undef SEEK_SET
//#undef SEEK_CUR
//#undef SEEK_END
#include <mpi.h>
#include <basics/utilities/UbException.h>

#ifdef USE_MPI_CXX_SYNTAX

namespace PbMpiTools
{
   /*======================================================================*/  
   // send a single value "value" of MPI_Datatype
   template <class T>
   inline void sendSingleValue(const T& value, MPI_Datatype datatype, int dest, int tag, MPI::Intracomm comm);
   
   /*======================================================================*/  
   // receives a single value "value" of MPI_Datatype
   template <class T>
   inline void receiveSingleValue(T& value, MPI_Datatype datatype, int source, int tag, MPI::Intracomm comm);
   
   /*======================================================================*/  
   // receives and returns a single value of MPI_Datatype
   // expample: int value = PbMpiTools::receiveSingleValue<int>(MPI::INT,0,10,comm);
   template <class T>
   inline T receiveSingleValue(MPI_Datatype datatype, int source, int tag, MPI::Intracomm comm);
   
   /*======================================================================*/  
   // sends bool value (doesn't work with template, why ever... stupid MPI)
   inline void sendBoolValue(const bool& value,int dest, int tag, MPI::Intracomm comm);

   /*======================================================================*/  
   // receives bool value (doesn't work with template, why ever... stupid MPI)
   inline bool receiveBoolValue(int source, int tag, MPI::Intracomm comm);

   /*======================================================================*/  
   // sends bool value (doesn't work with template, why ever... stupid MPI)
   inline void sendStringValue(const std::string& value,int dest, int tag, MPI::Intracomm comm);

   /*======================================================================*/  
   // receives bool value (doesn't work with template, why ever... stupid MPI)
   inline std::string receiveStringValue(int source, int tag, MPI::Intracomm comm);

   /*======================================================================*/  
   // send a vector of MPI_Datatype
   template <class T>
	inline void sendVector(const std::vector<T>& v, MPI_Datatype datatype, int dest, int tag, MPI::Intracomm comm);
	
   /*======================================================================*/  
   // receive a std::vector of MPI_Datatype
   template <class T>
   inline void receiveVector(std::vector<T>& v, MPI_Datatype datatype, int source, int tag, MPI::Intracomm comm);

   /*======================================================================*/  
   // receive a vector of MPI_Datatype and adds this vector to existing vector
   // ans returns number of received elements
   template <class T>
   inline int receiveVectorAndAddToVector(std::vector<T>& v, MPI_Datatype datatype, int source, int tag, MPI::Intracomm comm);

   /*======================================================================*/  
   // send a std::vector of strings
   inline void sendStringVector(const std::vector<std::string>& v, int dest, int tag, MPI::Intracomm comm);

   /*======================================================================*/  
   // send a vector of strings
   inline void receiveStringVector(std::vector<std::string>& v, int dest, int tag, MPI::Intracomm comm);
};

/*======================================================================*/  
// send a single value of MPI_Datatype
template <class T>
void PbMpiTools::sendSingleValue(const T& value, MPI_Datatype datatype, int dest, int tag, MPI::Intracomm comm)
{
   try{ comm.Send(&value, 1, datatype, dest, tag); }
   catch(MPI::Exception& e){ std::stringstream ss; ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
                             UB_THROW( UbException(UB_EXARGS,"catched with info at send size\n"+ss.str()) ); }
   catch(...)              { UB_THROW( UbException(UB_EXARGS,"unknown exception at send size") ); }
}
/*======================================================================*/  
template <class T>
void PbMpiTools::receiveSingleValue(T& value, MPI_Datatype datatype, int source, int tag, MPI::Intracomm comm) 
{
   try { comm.Recv(&value, 1, datatype, source, tag); }
   catch(MPI::Exception& e){ std::stringstream ss;ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
                             UB_THROW( UbException(UB_EXARGS,"MPI:Exception catched with info at receive size\n"+ss.str()) );}
   catch(...)              { UB_THROW( UbException(UB_EXARGS,"unknown exception at receive size") ); }
}
/*======================================================================*/  
template <class T>
T PbMpiTools::receiveSingleValue(MPI_Datatype datatype, int source, int tag, MPI::Intracomm comm) 
{
   T value;
   try { comm.Recv(&value, 1, datatype, source, tag); }
   catch(MPI::Exception& e){ std::stringstream ss;ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
                             UB_THROW( UbException(UB_EXARGS,"MPI:Exception catched with info at receive size\n"+ss.str()) );}
   catch(...)              { UB_THROW( UbException(UB_EXARGS,"unknown exception at receive size") ); }

   return value;
}
/*======================================================================*/  
// send a bool value (bool doesn't work with template, why ever)
void PbMpiTools::sendBoolValue(const bool& value,int dest, int tag, MPI::Intracomm comm)
{
   short dummy;
   if(value) dummy=1;                  
   else      dummy=0;

   try{ comm.Send(&dummy, 1, MPI::SHORT, dest, tag); }
   catch(MPI::Exception& e){ std::stringstream ss; ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
                             UB_THROW( UbException(UB_EXARGS,"MPI:Exception catched with info at send size\n"+ss.str()) ); }
   catch(...)              { UB_THROW( UbException(UB_EXARGS,"unknown exception at send size") ); }
}
/*======================================================================*/  
bool PbMpiTools::receiveBoolValue(int source, int tag, MPI::Intracomm comm) 
{
   short dummy;
   try { comm.Recv(&dummy, 1, MPI::SHORT, source, tag); }
   catch(MPI::Exception& e){ std::stringstream ss;ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
                             UB_THROW( UbException(UB_EXARGS,"MPI:Exception catched with info at receive size\n"+ss.str()) );}
   catch(...)              { UB_THROW( UbException(UB_EXARGS,"unknown exception at receive size") ); }

   return (dummy==1);
}
/*======================================================================*/  
// sends bool value (doesn't work with template, why ever... stupid MPI)
void PbMpiTools::sendStringValue(const std::string& value,int dest, int tag, MPI::Intracomm comm)
{
   std::vector<char> vec;
   for(std::size_t i=0; i<value.size(); i++)
      vec.push_back(value[i]);
   PbMpiTools::sendVector(vec,MPI::CHAR,dest,tag,comm);
}

/*======================================================================*/  
// receives bool value (doesn't work with template, why ever... stupid MPI)
std::string PbMpiTools::receiveStringValue(int source, int tag, MPI::Intracomm comm)
{
   std::vector<char> vec;
   PbMpiTools::receiveVector(vec,MPI::CHAR,source,tag,comm);
   std::string str;
   for(std::size_t i=0; i<vec.size(); i++)
      str+=vec[i];

   return str;
}
/*======================================================================*/  
// send a vector of MPI_Datatype
template <class T>
void PbMpiTools::sendVector(const std::vector<T>& v, MPI_Datatype datatype, int dest, int tag, MPI::Intracomm comm)
{
   // send size
   int size = (int)v.size();
   
   try{ comm.Send(&size, 1, MPI::INT, dest, tag); }
   catch(MPI::Exception& e){ std::stringstream ss; ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
                             UB_THROW( UbException(UB_EXARGS,"MPI:Exception catched with info at send size\n"+ss.str()) ); }
   catch(...)              { UB_THROW( UbException(UB_EXARGS,"unknown exception at send size") ); }
   
   if(size>0)
	{
      try{ comm.Send(&v[0], size, datatype, dest, tag); }
      catch(MPI::Exception& e){ std::stringstream ss; ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
                                UB_THROW( UbException(UB_EXARGS,"MPI:Exception catched with info at send vector<T>\n"+ss.str()) );}
      catch(...)              { UB_THROW( UbException(UB_EXARGS,"unknown exception at send vector<T>") ); }
   }
}
/*======================================================================*/  
// receive a vector of MPI_Datatype
template <class T>
void PbMpiTools::receiveVector(std::vector<T>& v, MPI_Datatype datatype, int source, int tag, MPI::Intracomm comm) 
{
   int size;

   try { comm.Recv(&size, 1, MPI::INT, source, tag); }
   catch(MPI::Exception& e){ std::stringstream ss;ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
                             UB_THROW( UbException(UB_EXARGS,"MPI:Exception catched with info at receive size\n"+ss.str()) );}
   catch(...)              { UB_THROW( UbException(UB_EXARGS,"unknown exception at receive size") ); }

   v.resize(size);

   if( size>0 )
   {
      try{ comm.Recv(&v[0], size, datatype, source, tag); }
      catch(MPI::Exception& e){ std::stringstream ss; ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
                                UB_THROW( UbException(UB_EXARGS,"MPI:Exception catched with info at receive vector\n"+ss.str()) ); }
      catch(...)              { UB_THROW( UbException(UB_EXARGS,"unknown exception at receive vector") ); }
   }
}
/*======================================================================*/  
// receive a vector of MPI_Datatype and adds this vector to existing vector
// return value is size of received elements
template <class T>
int PbMpiTools::receiveVectorAndAddToVector(std::vector<T>& v, MPI_Datatype datatype, int source, int tag, MPI::Intracomm comm) 
{
   int incommingSize;

   try { comm.Recv(&incommingSize, 1, MPI::INT, source, tag); }
   catch(MPI::Exception& e){ std::stringstream ss;ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
                             UB_THROW( UbException(UB_EXARGS,"MPI:Exception catched with info at receive size\n"+ss.str()) );}
   catch(...)              { UB_THROW( UbException(UB_EXARGS,"unknown exception at receive size") ); }

   int oldSize = (int)v.size();
   v.resize(oldSize+incommingSize);

   if( incommingSize>0 )
   {
      try{ comm.Recv(&v[oldSize], incommingSize, datatype, source, tag); }
      catch(MPI::Exception& e){ std::stringstream ss; ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
                                UB_THROW( UbException(UB_EXARGS,"MPI:Exception catched with info at receive vector\n"+ss.str()) ); }
      catch(...)              { UB_THROW( UbException(UB_EXARGS,"unknown exception at receive vector") ); }
   }

   return incommingSize;
}
/*======================================================================*/  
// send a vector of strings
void PbMpiTools::sendStringVector(const std::vector<std::string>& v, int dest, int tag, MPI::Intracomm comm)
{
   // send size
   int stringVectorSize = (int)v.size();

   try{ comm.Send(&stringVectorSize, 1, MPI::INT, dest, tag); }
   catch(MPI::Exception& e){ std::stringstream ss; ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
                             UB_THROW( UbException(UB_EXARGS,"MPI:Exception catched with info at send size\n"+ss.str()) ); }
   catch(...)              { UB_THROW( UbException(UB_EXARGS,"unknown exception at send size") ); }

   if(stringVectorSize>0)
   {
      std::vector<int> singleStringSizes(stringVectorSize+1);
      int nofChars = 0;
      for(int i=0; i<stringVectorSize; i++)
         nofChars += singleStringSizes[i] = (int)v[i].length();
      singleStringSizes[stringVectorSize] = nofChars;

      try{ comm.Send(&singleStringSizes[0], stringVectorSize+1, MPI::INT, dest, tag); }
      catch(MPI::Exception& e){ std::stringstream ss; ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
                                UB_THROW( UbException(UB_EXARGS,"MPI:Exception catched with info at send vector<T>\n"+ss.str()) );}
      catch(...)              { UB_THROW( UbException(UB_EXARGS,"unknown exception at send vector<T>") ); }

      std::vector<char> charVector(nofChars);
      int pos = 0;
      for(int i=0; i<stringVectorSize; i++)
         for(int j=0; j<singleStringSizes[i]; j++)
            charVector[pos++] = v[i][j];      

      try{ comm.Send(&charVector[0], nofChars, MPI::CHAR, dest, tag); }
      catch(MPI::Exception& e){ std::stringstream ss; ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
                                UB_THROW( UbException(UB_EXARGS,"MPI:Exception catched with info at send vector<T>\n"+ss.str()) );}
      catch(...)              { UB_THROW( UbException(UB_EXARGS,"unknown exception at send vector<T>") ); }
   }
}
/*======================================================================*/  
// send a vector of strings
void PbMpiTools::receiveStringVector(std::vector<std::string>& v, int source, int tag, MPI::Intracomm comm)
{
   // send size
   int stringVectorSize;
   try { comm.Recv(&stringVectorSize, 1, MPI::INT, source, tag); }
   catch(MPI::Exception& e){ std::stringstream ss; ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
                             UB_THROW( UbException(UB_EXARGS,"MPI:Exception catched with info at send size\n"+ss.str()) ); }
   catch(...)              { UB_THROW( UbException(UB_EXARGS,"unknown exception at send size") ); }

   v.clear();
   v.resize(stringVectorSize);

   if(stringVectorSize>0)
   {
      std::vector<int> singleStringSizes(stringVectorSize+1);

      try{ comm.Recv(&singleStringSizes[0], stringVectorSize+1, MPI::INT, source, tag); }
      catch(MPI::Exception& e){ std::stringstream ss; ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
                                UB_THROW( UbException(UB_EXARGS,"MPI:Exception catched with info at send vector<T>\n"+ss.str()) );}
      catch(...)              { UB_THROW( UbException(UB_EXARGS,"unknown exception at send vector<T>") ); }

      int nofChars = singleStringSizes[stringVectorSize];
      std::vector<char> charVector(nofChars);

       try{ comm.Recv(&charVector[0], nofChars, MPI::CHAR, source, tag); }
       catch(MPI::Exception& e){ std::stringstream ss; ss<<"MPI::Exception error_string="<<e.Get_error_string()<<std::endl;
                                 UB_THROW( UbException(UB_EXARGS,"MPI:Exception catched with info at send vector<T>\n"+ss.str()) );}
       catch(...)              { UB_THROW( UbException(UB_EXARGS,"unknown exception at send vector<T>") ); }
      
      int pos=0;
      for(int i=0; i<stringVectorSize; i++)
         for(int j=0; j<singleStringSizes[i]; j++)
            v[i].push_back(charVector[pos++]);      
   }
}

#endif

#endif //PBMPITOOLS_H
