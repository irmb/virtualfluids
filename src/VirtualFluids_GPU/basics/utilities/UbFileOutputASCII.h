//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBFILEOUTPUTASCII_H
#define UBFILEOUTPUTASCII_H

#include <iomanip>
#include <fstream>
#include <iostream>

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbFileOutput.h>
#define NOMINMAX
/*=========================================================================*/
/*  UbFileOutputASCII                                                             */
/*                                                                         */
/**
...
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 23.11.04
*/ 

/*
usage: ...
*/

class UbFileOutputASCII : public UbFileOutput
{
public:
   UbFileOutputASCII() : UbFileOutput() {}
   UbFileOutputASCII(const std::string& filename, const bool& createPath=true,  const int& precision=15);             
   UbFileOutputASCII(const std::string& filename, CREATEOPTION opt, const bool& createPath=true, const int& precision=15);
   
   bool open(const std::string& filename, CREATEOPTION opt=OUTFILE);
   
   void writeBool(const bool& value, const int& width=0);
   void writeDouble(const double& value, const int& width=0);
	void writeFloat(const float& value, const int& width=0);
	void writeInteger(const int& value, const int& width=0);
   void writeLong(const long& value, const int& width=0);
   void writeLongLong(const long long& value, const int& width=0);
   void writeSize_t(const std::size_t& value, const int& width=0);
   void writeChar(const char& value, const int& width=0);
   void writeString(const std::string& value, const int& width=0);
   void writeStringOnly(const std::string& value);
   void writeLine(const std::string& value, const int& width=0);
   void writeLine();
  
   void setPrecision(const int& precision);
   int  getPrecision() { return (int)outfile.precision(); }

   void setCommentIndicator(char commentindicator) {this->commentindicator = commentindicator;} 
   
   void writeCommentLine(const std::string& line);
   void writeCommentLine(char indicator, const std::string& line);
   void writeCopyOfFile(const std::string& filename);

   FILETYPE getFileType() { return ASCII; }

   template< typename T >
   friend inline UbFileOutputASCII& operator<<(UbFileOutputASCII& file, const T& data) 
   {
      file.outfile<<data;
      return file;
   }
};

#endif


