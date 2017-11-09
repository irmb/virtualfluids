//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBFILEOUTPUT_H
#define UBFILEOUTPUT_H            

#include <iomanip>
#include <fstream>
#include <iostream>
#include <string>

#include <basics/utilities/UbException.h>
#define NOMINMAX
/*=========================================================================*/
/*  UbFileOutput                                                             */
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

class UbFileOutput
{
public:
   enum CREATEOPTION {OUTFILE=0, INANDOUTFILE=1, APPENDFILE=2};      
   enum FILETYPE {ASCII, BINARY};      

public:
   UbFileOutput() : filename(""), commentindicator('C') {  }
   UbFileOutput(const std::string& filename)  : filename(filename), commentindicator('C') { }             
   virtual ~UbFileOutput() { outfile.flush();outfile.close(); }

   virtual bool open(const std::string& filename, CREATEOPTION opt=OUTFILE) = 0;

   virtual bool operator!() { return !(outfile); }
   virtual bool isOpen()    { return !(!(outfile)); }
   
   virtual void flush() { outfile.flush(); }
   virtual void close() { outfile.close(); }
   
   std::fstream::pos_type tellp() { return outfile.tellp(); } 
   void seekp(std::fstream::off_type offset, std::ios::seekdir origin=std::ios::beg) { outfile.seekp(offset,origin); }  

   virtual void writeInteger(const int& value, const int& width=0)=0;
   virtual void writeLong(const long& value, const int& width=0)=0;
   virtual void writeLongLong(const long long& value, const int& width=0)=0;
   virtual void writeDouble(const double& value, const int& width=0)=0;
	virtual void writeFloat(const float& value, const int& width=0)=0;
   virtual void writeBool(const bool& value, const int& width=0)=0;
   virtual void writeSize_t(const std::size_t& value, const int& width=0)=0;
   virtual void writeChar(const char& value, const int& width=0)=0;
   virtual void writeString(const std::string& value, const int& width=0)=0;
   virtual void writeStringOnly(const std::string& value)=0;
	virtual void writeLine(const std::string& value, const int& width=0)=0;
	virtual void writeLine()=0;

	virtual void writeCommentLine(const std::string& line)=0;
	virtual void writeCommentLine(char indicator, const std::string& line)=0;
   virtual void writeCopyOfFile(const std::string& filename)=0;
	
   virtual void setCommentIndicator(char commentindicator) {this->commentindicator = commentindicator;} 
   
   virtual void setPrecision(const int& precision)=0;
   virtual int  getPrecision()=0;

   //returns "ASCII", "BINARY"
   virtual FILETYPE getFileType()=0;

   //returns file extension:
   //e.g. "./../test/ich.inp" -> "inp", "./../test/ich" -> ""
   virtual std::string getFileExtension()  
   {
	   std::size_t pos1 = filename.rfind("/");
      if(pos1==std::string::npos) pos1 = 0;
      std::size_t pos2 = filename.rfind(".");
      if(pos2!=std::string::npos && pos2>pos1)
         return filename.substr(pos2+1);

      return "";
   }				

   virtual std::string  getFileName() {return this->filename;}
protected:
   std::ofstream outfile;
   std::string   filename; 
   char          commentindicator; 
};

#endif //UBFILEOUTPUT_H
