//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBFILEOUTPUTBINARY_H
#define UBFILEOUTPUTBINARY_H

#include <iomanip>
#include <fstream>
#include <iostream>


#include <basics/utilities/UbException.h>
#include <basics/utilities/UbFileOutput.h>

/*=========================================================================*/
/*  UbFileOutputBinary                                                             */
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

class UbFileOutputBinary : public UbFileOutput
{
public:
   UbFileOutputBinary() : UbFileOutput() {}
   UbFileOutputBinary(const std::string& filename, const bool& createPath=true);
   UbFileOutputBinary(const std::string& filename, UbFileOutput::CREATEOPTION opt, const bool& createPath);
   
   bool open(const std::string& filename, UbFileOutput::CREATEOPTION opt=OUTFILE);

   void writeInteger(const int& value, const int& width=0);
   void writeDouble(const double& value, const int& width=0);
	void writeFloat(const float& value, const int& width=0);
	void writeBool(const bool& value, const int& width=0);
   void writeChar(const char& value, const int& width=0);
   void writeSize_t(const std::size_t& value, const int& width=0);
   void writeString(const std::string& value, const int& width=0);
   void writeStringOnly(const std::string& value);
   void writeLine(const std::string& value, const int& width=0);
   void writeLine();
   void writeCommentLine(const std::string& line);
   void writeCommentLine(char indicator, const std::string& line);
   void writeCopyOfFile(const std::string& filename);

   void setPrecision(const int& precision);
   int  getPrecision();

   FILETYPE getFileType() { return BINARY; }

   template< typename T >
   friend inline UbFileOutputBinary& operator<<(UbFileOutputBinary& file, const T& data) 
   {
      file.outfile.write((char*)&data,sizeof(T));
      return file;
   }
};

#endif


