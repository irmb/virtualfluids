#ifndef LOGWRITER_H
#define LOGWRITER_H

#include <iostream>
#include <fstream>
//#include <string>

//#include "Utilities/StringUtil.hpp"


////////////////////////////////////////////////////////////////////////////////
class LogWriter
{
public:
   LogWriter()
   {  
      consoleOut = false;
   }
   LogWriter(std::string fname)
   {
      consoleOut = false;
      this->fname = fname;
   }
   void setName(std::string name)
   {
      this->fname = name;
   }
   void setConsoleOut(bool flag)
   {
      consoleOut = flag;
   }
   void clearLogFile()
   {
      ostr.open(fname.c_str(), std::ios_base::out);
      if (ostr.bad())
      {
         std::string exceptionText = "Error: Output file/directory not found! LogWriter::operator << \n";
         throw exceptionText;
      }
      ostr << "";
      ostr.close();
   }
   template <typename T>
   LogWriter&  operator << (const T& arg)
   {
      ostr.open(fname.c_str(), std::ios_base::app);
      if (ostr.bad())
      {
         //std::cout << "Error: Output file/directory not found! LogWriter::operator <<" << std::endl;
         //return *this;
         std::string exceptionText = "Error: Output file/directory not found! LogWriter::operator << \n";
         throw exceptionText;
      }
      ostr << arg;
      ostr.close();
      if(consoleOut) std::cout << arg << std::flush;
      return *this;
   }
protected:
private:
   std::string fname;
   std::ofstream ostr;
   bool consoleOut;
};

#endif	
