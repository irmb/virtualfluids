//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBEXCEPTION_H
#define UBEXCEPTION_H

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>

#include "./UbTuple.h"

/*=========================================================================*/
/*  UbException                                                             */
/*                                                                         */
/**
This Class provides the base for exception handling.
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0  - 23.11.04
@version 1.5  - 14.03.08
@version 1.6  - 31.03.08 derivation from std::run_time_error
@version 1.6a - helper marco UB_EXARGS
*/ 

/*
usage: throw UbException("error message");
       throw UbException(__FILE__, __LINE__,"error message");
       throw UbException(__FILE__, __LINE__,UB_FUNCTION,"error message");
       throw UbException(UB_EXARGS,"error"); //same as above
*/

//Macro UB_FUNCTION: figures out the method/function name (platform dependant)
#if defined(__GNUC__) || (defined(__MWERKS__) && (__MWERKS__ >= 0x3000)) || (defined(__ICC) && (__ICC >= 600))
 # define UB_FUNCTION __PRETTY_FUNCTION__
#elif defined(__DMC__) && (__DMC__ >= 0x810)
 # define UB_FUNCTION __PRETTY_FUNCTION__
#elif defined(__FUNCSIG__)
 # define UB_FUNCTION __FUNCSIG__
#elif (defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 600)) || (defined(__IBMCPP__) && (__IBMCPP__ >= 500))
 # define UB_FUNCTION __FUNCTION__
#elif defined(__BORLANDC__) && (__BORLANDC__ >= 0x550)
 # define UB_FUNCTION __FUNC__
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901)
 # define UB_FUNCTION __func__
#else
 # define UB_FUNCTION "(unknown)"
#endif

//Helper Marco
#ifndef SWIG
#define UB_EXARGS __FILE__,__LINE__,UB_FUNCTION
#endif
class UbException : public std::runtime_error
{
public:
   typedef UbTuple< std::string, int, std::string, std::string > ExceptionData;
public:
   //////////////////////////////////////////////////////////////////////////
   //constructors
   UbException()
      : std::runtime_error("")
   { 
   }
   /*==========================================================*/
   UbException(const std::string& str)
      : std::runtime_error("")
   {
      this->addInfo(str);		
   }
   /*==========================================================*/
   UbException(const std::string& file, const int& line, const std::string& err_str)
      : std::runtime_error("")
   {
      this->addInfo(file,line,"unknown",err_str);		
   }
   /*==========================================================*/
   //UbException(const char* file, const int& line, const char* function, const std::string& err_str)
   UbException(const std::string& file, const int& line, const std::string& function, const std::string& err_str)
      : std::runtime_error("")
   {
      this->addInfo(file,line,function,err_str);		
   }
   //////////////////////////////////////////////////////////////////////////
   //destructor
   virtual ~UbException() throw() { }
   //////////////////////////////////////////////////////////////////////////
   //virtual public methods
   //returns  exception-string
   virtual const char* what() const throw()
   {
      exceptionString = this->toString();
      return exceptionString.c_str();  //ansonsten ist das Verhalten anschlieﬂend undefiniert!
   }
   /*==========================================================*/
   virtual void addInfo(const std::string& err_str)	 
   { 
      exceptionData.push_back( makeUbTuple( (std::string)"-", 0, (std::string)"unknown", err_str) ); 
   }
   /*==========================================================*/
   //add exception
   virtual void addInfo(const std::string& file, const int& line, const std::string& function, const std::string& err_str)	 
   { 
      exceptionData.push_back( makeUbTuple( file, line, function, err_str ) ); 
   }
   /*==========================================================*/
   //returns exception-string with all calles exceptions
   virtual const std::vector<std::string> getInfo() const
   { 
      std::vector<std::string> tmp;
      for(std::size_t i=0; i<exceptionData.size(); i++)
      {
         std::stringstream str;
         str << val<1>( exceptionData[i] ) << ", " 
             << val<2>( exceptionData[i] ) << ", " 
             << val<3>( exceptionData[i] ) << ", " 
             << val<4>( exceptionData[i] );
         tmp.push_back( str.str());
      }
      return tmp; 
   }
   /*==========================================================*/
   //returns exception-string with all calles exceptions and detailes informations
   virtual std::string toString() const
   { 
      std::stringstream str("UbExeption");
      
      for(std::size_t i=0; i<exceptionData.size(); i++)
         str<<(std::string)"caller[" << i << "]\n"
            <<"  - file:     "<< val<1>( exceptionData[i] )<<"\n"
            <<"  - line:     "<< val<2>( exceptionData[i] )<<"\n"
            <<"  - function: "<< val<3>( exceptionData[i] )<<"\n"
            <<"  - what:     "<< val<4>( exceptionData[i] )<< std::endl; 

      return str.str();
   }

protected:
   //////////////////////////////////////////////////////////////////////////
   //protected member
   std::vector< ExceptionData > exceptionData;
   mutable std::string exceptionString;
};

//overlading operator <<
inline std::ostream& operator<<(std::ostream& os, const UbException& e)
{
   return os<<e.toString();
}

#endif //UBEXCEPTION_H
