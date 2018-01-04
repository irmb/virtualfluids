//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBSYSTEM_H
#define UBSYSTEM_H

#if defined(_WIN32) || defined(_WIN64)
   #define UBSYSTEM_WINDOWS
   #include <process.h>
   #include <io.h>
   #include <direct.h>
   //#ifndef _WINSOCK2API_  //ansonsten gibt es mecker bei #include "Windows.h" und ::Sleep()
   //   #define _WINSOCK2API_
   //   #include<WinSock2.h> 
   //#endif
  #include <windows.h>
  //#include <Windows.h>
  //#include <tchar.h>
#elif defined(__APPLE__)
   #define UBSYSTEM_APPLE
   #include "dirent.h"
   #include "sys/stat.h"
   #include <sys/syscall.h>
   #include <sys/stat.h>
#elif (defined(__amd64) || defined(__amd64__) || defined(__unix__) || defined(__CYGWIN__)) && !defined(__AIX__) 
   #define UBSYSTEM_LINUX
   #include "dirent.h"
   #include "sys/stat.h"
   #include <sys/syscall.h>
   #include <sys/stat.h>
   #include <unistd.h>
   #include <string.h>
#elif defined(__AIX__)
   #define UBSYSTEM_AIX
   #include "dirent.h"
   #include <unistd.h>
   #include <sys/stat.h>
   #include <sys/types.h>
#else
   #error "UbSystem::UnknownMachine"
#endif



#if defined(min) || defined(max) //daruch kann man sich spaeter #undef min; #undef max erparen
#   error add NOMINMAX to preprocessor defines
#endif


#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <algorithm>
#include <typeinfo>
#include <cctype> //for toupper
#include <ctime>

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbLogger.h>

#if defined(CAB_BOOST)
#include <boost/thread.hpp>
#endif // CAB_BOOST

//DEFINE TO STRING
//e.g. #define FOO hallo
//     -> QUOTEME(FOO) == "hallo"
#define _QUOTEME(x) #x
#define QUOTEME(x) _QUOTEME(x)

//allg.:
//const int * C1        -> C1 is variable pointer to a constant integer
//int const * C2        -> C2 is variable pointer to a constant integer (same as above)
//int * const C3        -> C3 is constant pointer to a variable integer
//int const * const C4  -> C4 is constant pointer to a constant integer

//////////////////////////////////////////////////////////////////////////
//UbSystem
//////////////////////////////////////////////////////////////////////////
namespace UbSystem
{
   template<bool> struct ub_static_assert;     //deklaration (ub_xxx da static_assert in C++0x ein keyword werden wird)
   template<> struct ub_static_assert<true>{}; //deklaration + definition der spezialisierung fuer "true"
                                               //ub_static_assert<false> fuehrt zu compiler fehler, da dafuer
                                               //keine implementierung vorhanden!  //UB_STATIC_ASSERT(false)

   /*==========================================================*/
   inline void sleepMs(const unsigned int& msec)
   {
      #if defined UBSYSTEM_WINDOWS
         ::Sleep(  (msec==0) ? 1 : msec );  // +1 here causes a context switch if SleepMSec(0) is called
      #elif defined(UBSYSTEM_LINUX) || defined(UBSYSTEM_APPLE) || defined(UBSYSTEM_AIX)
         ::usleep(1000*msec);
      #else
         #error "UbSystem::sleepMSec - UnknownMachine"
      #endif
   }
   /*==========================================================*/
   inline void sleepS(const unsigned int& sec)
   {
      #if defined UBSYSTEM_WINDOWS
         ::Sleep( (sec==0) ? 1 : sec*1000 );  // +1 here causes a context switch if sleepS(0) is called
      #elif defined(UBSYSTEM_LINUX) || defined(UBSYSTEM_APPLE) || defined(UBSYSTEM_AIX)
         ::sleep(sec);
      #else
         #error "UbSystem::sleepS - UnknownMachine"
      #endif
   }
   /*==========================================================*/
   //checks if the bits of bitmask are set in value
   template<typename T>
   inline bool bitCheck(const T& value, const T& bitmask)
   {
      return  ( (value & bitmask) == bitmask);
   }
   /*==========================================================*/
   //checks if the bits of bitmask are set in value
   template<typename T>
   inline void setBit(T& value, const T& bitmask)
   {
      value |= bitmask;
   }
   /*==========================================================*/
   template<typename T>
   inline void unsetBit(T& value, const T& bitmask)
   {
      value &= ~bitmask;
   }
   /*==========================================================*/
   //returns bitmask as string e.g. 0001 0100 1101
   template<typename T>
   inline std::string getBitString(const T& value)
   {
      std::stringstream text;
      for(int i=sizeof(value)*8-1/*8 bits per byte*/; i>=0; i--)
      {
         text<<(char) ( ((value>>i) & 1) + '0');
         if(i%4 == 0 && i>0) text<<' ';
      }
      return text.str();
   }
   /*==========================================================*/
   //converts string to type T
   // usage: int x = stringTo<int>("123");
   template<typename T>
   inline T stringTo(const std::string& s)
   {
     std::istringstream iss(s);
     T x;
     iss >> x;
     if(!iss)
        UB_THROW( UbException(UB_EXARGS," cannot convert \""+s+"\" to type <"+static_cast<std::string>(typeid(x).name())+">") );

     return x;
   }
   /*==========================================================*/
   // usage: string s = toString(x);
   template<typename T>
   inline std::string toString(const T& x, int precision=15)
   {
     std::ostringstream oss;
     oss<<std::setprecision(precision);
     oss<<x;
     return oss.str();
   }
   /*==========================================================*/
   //e.g. str="iHcsnW" -> "IHCSNW"
   inline std::string toUpperString(const std::string& str)
   {
      std::string tmp(str);
      std::transform(tmp.begin(),tmp.end(),tmp.begin(), static_cast<int (*)(int)>(std::toupper));

      return tmp;
   }
   /*==========================================================*/
   //e.g. str="iHcsnW" -> "ihcsnw"
   inline std::string toLowerString(const std::string& str)
   {
      std::string tmp(str);
      std::transform(tmp.begin(),tmp.end(),tmp.begin(), static_cast<int (*)(int)>(std::tolower));

      return tmp;
   }
   /*==========================================================*/
   // usage: std::string s = replaceInString(str,"\\","/");
   //        std::string s = replaceInString(str,"ich","du");
   static std::string replaceInString(std::string original, const std::string& replace, const std::string& replaceWith )
   {
      size_t pos=0;
      while( (pos=original.find(replace,pos))!=std::string::npos )
      {
         original.replace(pos,replace.size(),replaceWith);
         pos+=replaceWith.size();
      }
      return original;
   }
   /*==========================================================*/
   //returns content of an enviroment variable
   inline std::string getEnv(const std::string& var)
   {
      char* str = getenv( var.c_str());
      if(  str == NULL  ) 
      {
         return std::string("");
      }
      
      return static_cast<std::string>( str );
   }
   /*==========================================================*/
   inline bool isDirectory(const std::string& dir, const unsigned& attemptions = 3)
   {
      if( dir.empty() ) 
         UB_THROW( UbException(UB_EXARGS,"dir is empty") );
      
      std::string path = UbSystem::replaceInString(dir,"\\","/");

      #if defined UBSYSTEM_WINDOWS
         #ifndef _UNICODE 
            if( _access(path.c_str(), 0  ) == -1 ) return false;
         #else
            if( _waccess(path.c_str(), 0 ) == -1 ) return false;
         #endif
      #elif defined(UBSYSTEM_LINUX) || defined(UBSYSTEM_APPLE) || defined(UBSYSTEM_AIX)
         struct stat stFileInfo;
         if( stat(path.c_str(),&stFileInfo) != 0) 
         {
            return false;
         } 
      #endif
      
      return true;
   }
   /*==========================================================*/
   // usage:  makeDirectory("c:/temp");
   //         makeDirectory("c:/temp/");
   // return: true  -> successful
   //         false -> failed
   #if defined(CAB_BOOST) 
      static boost::mutex mtx_makeDirectory;
   #endif
   inline bool makeDirectory(const std::string& dir, const unsigned& attemptions = 3)
   {
      UBLOG(logDEBUG5,"UbSystem::makeDirectory - start, dir="<<dir<<" #attemptions="<<attemptions);

      if( dir.empty() ) UB_THROW( UbException(UB_EXARGS,"dir is empty") );
      std::string path = UbSystem::replaceInString(dir,"\\","/");

      bool dirCreated = true;
      #if defined UBSYSTEM_WINDOWS
         if(path[path.size()-1] != '/') path+="/";
         size_t  pos = 0;
         while( ( pos=path.find("/",pos+1) ) != std::string::npos )
         {
            std::string tmpdir = path.substr(0,pos);
            #if defined(CAB_BOOST) 
            boost::mutex::scoped_lock lock(mtx_makeDirectory);
            #endif
            if( 
                #ifndef _UNICODE 
                 _access(tmpdir.c_str(), 0 ) == -1 && _mkdir(tmpdir.c_str() ) == -1
                #else
                 _waccess(tmpdir.c_str(), 0) == -1 && _wmkdir(tmpdir.c_str()) == -1
                #endif
               )
               {
                  UBLOG(logDEBUG5,"UbSystem::makeDirectory-  dir=\""<<tmpdir<<"\" doesn't exit or makedir failed");
                  dirCreated = false;
                  break;
               }
         }
      #elif defined(UBSYSTEM_LINUX) || defined(UBSYSTEM_APPLE) || defined(UBSYSTEM_AIX)
         std::string command = "mkdir -p \""+path+"\"";
         {
            #if defined(CAB_BOOST) 
               boost::mutex::scoped_lock lock(mtx_makeDirectory);
            #endif
            if(system(command.c_str())!=0)
            {
               UBLOG(logDEBUG5,"UbSystem::makeDirectory-  dir=\""<<path<<"\" doesn't exit or makedir failed");
               dirCreated = false;
            }
         }
      #else
         #error "UbSystem::makeDirectory - UnknownMachine"
      #endif

      if(!dirCreated && attemptions > 1)
      {
         UBLOG(logDEBUG5,"UbSystem::makeDirectory - internal call of UbSystem::makeDirectory");
         UbSystem::sleepMs(500);
         dirCreated = UbSystem::makeDirectory(path, attemptions-1);
      }
      
      UBLOG(logDEBUG5,"UbSystem::makeDirectory - end (success="<<dirCreated<<", attemptions = "<<attemptions<<")");
      return dirCreated;
   }
   /*==========================================================*/
#if defined(CAB_BOOST) 
   static boost::mutex mtx_removeDirectory;
#endif
   inline int removeDirectory(const std::string& dir)
   {
      #if defined(CAB_BOOST) 
         boost::mutex::scoped_lock lock(mtx_removeDirectory);
      #endif
      std::string command = "rmdir \""+dir+"\"";
      return std::system(command.c_str());
   }
   /*==========================================================*/
   // usage  : getPathFromString("c:/temp/foo.txt");
   //returns: "c:/temp"
   // usage  : getPathFromString("c:\\temp\\foo.txt");
   //returns: "c:/temp"
   // usage  : getPathFromString("foo.txt");
   // returns: ""
   inline std::string getPathFromString(const std::string& fileStringWithPath)
   {
      std::string tmp = UbSystem::replaceInString(fileStringWithPath,"\\","/");
      std::size_t last = tmp.rfind("/");
      if(last!=std::string::npos) tmp.resize(last);
      else                        tmp = "";
      return tmp;
   }
   /*==========================================================*/
   // usage  : getFilenameFromString("c:/temp/foo.txt");
   // returns: "foo.txt"
   // usage  : getFilenameFromString("c:/temp/foo.txt",false);
   // returns: "foo"
   // usage  : getFilenameFromString("c:/temp/");
   // returns: ""
   inline std::string getFilenameFromString(const std::string& fileStringWithPath, bool withExtension = true)
   {
      std::string tmp = UbSystem::replaceInString(fileStringWithPath,"\\","/");
      
      //remove path
      std::size_t last = tmp.rfind("/");
      if(last!=std::string::npos && (last+1)<tmp.size()) tmp.erase(0,last+1);
      
      //remove extension
      if(!withExtension)
      {
         last = tmp.rfind(".");
         if(last!=std::string::npos) tmp.erase(last);
      }

      return tmp;
   }
   /*==========================================================*/
   inline int getProcessID()
   {
      #if defined UBSYSTEM_WINDOWS
         return _getpid();
      #elif defined(UBSYSTEM_LINUX) || defined(UBSYSTEM_APPLE) || defined(UBSYSTEM_AIX)
         return getpid();
      #else
         #error "int UbSystem::getProcessID() - UnknownMachine"
      #endif
   }
   /*==========================================================*/
   inline unsigned long getCurrentThreadID()
   {
      #if defined UBSYSTEM_WINDOWS
         return (unsigned long)GetCurrentThreadId();
      #elif defined(UBSYSTEM_LINUX) || defined(UBSYSTEM_APPLE)
         return (unsigned long)syscall(SYS_gettid);
      #elif defined(UBSYSTEM_AIX)
         return (unsigned long) getpid(); //WORKAROUND for IBM (for get thread id is another function necessary) 
      #else
         #error "unsigned long UbSystem::getCurrentThreadID() - UnknownMachine"
      #endif
   }
   /*==========================================================*/
   inline bool isBigEndian()
   {
      short word = 0x4321;
      if((*(char*)& word) != 0x21 ) return true;
      else                           return false;
   }
   /*==========================================================*/
   inline bool isLittleEndian()
   {
      return !isBigEndian();
   }
   /*==========================================================*/
   inline std::string getTimeStamp()
   {
      time_t t = time(NULL);
      tm* localTime = localtime(&t); 	
      
      std::stringstream tmp;
      tmp.fill('0');
      
      tmp << localTime->tm_year+1900 
          << "." << std::setw(2) <<localTime->tm_mon+1
          << "." << std::setw(2) << localTime->tm_mday 
          << "@" << std::setw(2) << localTime->tm_hour  
          << "." << std::setw(2) << localTime->tm_min   
          << "." << std::setw(2) << localTime->tm_sec  ;

      return tmp.str();
   }
   /*==========================================================*/
   //swap Byte Order
   //usage: int test = 8;
   //       swapByteOrder((unsigned char*)& test, sizeof(int))
   //#define ByteSwap5(x) ByteSwap((unsigned char *) &x,sizeof(x))
   inline void swapByteOrder(unsigned char* toSwap, int length)
   {
      register int i = 0;
      register int j = length-1;
      while(i<j)
      {
         std::swap(toSwap[i], toSwap[j]);
         i++, j--;
      }
   }
   //////////////////////////////////////////////////////////////////////////
   //get host name
   inline std::string getMachineName()
   {
      char Name[150];
      int i = 0;

#ifdef UBSYSTEM_WINDOWS
      TCHAR infoBuf[150];
      DWORD bufCharCount = 150;
      memset(Name, 0, 150);
      if (GetComputerName(infoBuf, &bufCharCount))
      {
         for (i = 0; i<150; i++)
         {
            Name[i] = infoBuf[i];
         }
      }
      else
      {
         strcpy(Name, "Unknown_Host_Name");
      }
#else
      memset(Name, 0, 150);
      gethostname(Name, 150);
#endif
      return std::string(Name);
   }

   //////////////////////////////////////////////////////////////////////////
   // generic IfThenElse - start
   //////////////////////////////////////////////////////////////////////////
   // primary template: yield second or third argument depending on first argument
   template<bool C, typename Ta, typename Tb>
   class IfThenElse;

   // partial specialization: true yields second argument
   template<typename Ta, typename Tb>
   class IfThenElse<true, Ta, Tb> {
   public:
      typedef Ta ResultT;
   };

   // partial specialization: false yields third argument
   template<typename Ta, typename Tb>
   class IfThenElse<false, Ta, Tb> {
   public:
      typedef Tb ResultT;
   };
   //////////////////////////////////////////////////////////////////////////
   // generic IfThenElse - end
   //////////////////////////////////////////////////////////////////////////

   //////////////////////////////////////////////////////////////////////////
   //help struct for overloading methods in template classes for specific types
   //////////////////////////////////////////////////////////////////////////
   template< typename T>
   struct type2type
   {
      typedef T type;
   };


   //////////////////////////////////////////////////////////////////////////
   // pair selector
   //////////////////////////////////////////////////////////////////////////
   template <typename Pair>
   struct select1st
   {
      typedef Pair argument_type ;
      typedef typename Pair::first_type result_type ;

      const result_type&  operator()(const argument_type &p) const
      {
         return p.first ;
      }
   };

   template <typename Pair>
   struct select2nd
   {
      typedef Pair argument_type ;
      typedef typename Pair::second_type result_type ;

      const result_type& operator()(const argument_type &p) const
      {
         return p.second ;
      }
   };

};

#define UB_STATIC_ASSERT(expr) static_cast<void>(sizeof( UbSystem::ub_static_assert<expr> ));
//zum ueberpruefen von STATISCHEN ausdruecken waehrend der compile-zeit
//--> Ausdruecke muessen schon ZUR compilerzeit auswertbar sein !!!
//Anwendung z.B. zur Ueberpruefung von Funktionalitaeten, wie z.B. bei UbMath::getNegativeInfinity<double>();
//
//Grund fuer macro ist einfach, dass es besser anzuwenden ist in der praxis!
//ansonsten würde es so aussehen:
//     UbSystem::ub_static_assert< aaa == 1 > test();
//    da ist  UB_STATIC_ASSERT(aaa == 1); schoener
//
//um das zu vermeiden machtman hier diesen static_cast<void>(sizeof(...) )
//Code-Snippet:
// struct Test { const static bool m_const_bool = true; bool m_bool; };
// int main() {
//  UB_STATIC_ASSERT( Test::m_const_bool == true );
//  --> okay, assert bestanden
//  UB_STATIC_ASSERT( Test::m_const_bool == false); //:
//  --> assert nicht bestanden z.B. error C2027: use of undefined type 'UbSystem::ub_static_assert<__formal> with __formal = false --> funzt nicht. fehler im code
//  UB_STATIC_ASSERT( Test::m_bool == true );
//  --> nicht erlaubt, da m_bool nicht statisch und nicht const ist.
//}

#endif //UBSYSTEM_H
