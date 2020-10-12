#include <basics/utilities/UbFileOutputASCII.h>
#include <basics/utilities/UbSystem.h>
#include <basics/utilities/UbInfinity.h>
#include <basics/utilities/UbMath.h>
#include <cstring>

using namespace std;

UbFileOutputASCII::UbFileOutputASCII(const string& filename, const bool& createPath, const int&  /*precision*/)
   : UbFileOutput(filename)
{
	this->commentindicator = 'C'; 
	this->setPrecision(20);

   outfile.open(filename.c_str(),ios::out);
   
   if(!outfile && createPath) 
   {
      outfile.clear(); //flags ruecksetzen (ansonsten liefert utern if(!outfile) weiterhin true!!!
      string path = UbSystem::getPathFromString(filename);
      if(path.size()>0) 
      {
         UbSystem::makeDirectory(path);
         outfile.open(filename.c_str(),ios::out);
      }
   }

      if(!outfile) UB_THROW( UbException(UB_EXARGS,"couldn't open file:\n "+filename) );
}
/*==========================================================*/
UbFileOutputASCII::UbFileOutputASCII(const std::string& filename, CREATEOPTION opt, const bool& createPath, const int& precision)
   : UbFileOutput(filename)
{
	this->commentindicator = 'C'; 
   this->setPrecision(precision);
	
   if(!this->open(filename,opt) && createPath) 
   {
      outfile.clear(); //flags ruecksetzen (ansonsten liefert utern if(!outfile) weiterhin true!!!
      string path = UbSystem::getPathFromString(filename);
      if(path.size()>0) UbSystem::makeDirectory(path);

      this->open(filename,opt);
   }

   if(!outfile) UB_THROW( UbException(UB_EXARGS,"couldn't open file:\n "+filename) );
}
/*==========================================================*/
bool UbFileOutputASCII::open(const std::string& filename, CREATEOPTION opt)
{
   outfile.close();
   outfile.clear(); //setzt flags zurueck
   this->filename = filename;

   if     (opt==UbFileOutput::OUTFILE      ) outfile.open(this->filename.c_str(),ios::out); 
   else if(opt==UbFileOutput::INANDOUTFILE ) outfile.open(this->filename.c_str(),ios::out | ios::in);
   else if(opt==UbFileOutput::APPENDFILE   ) outfile.open(this->filename.c_str(),ios::app);
   else UB_THROW( UbException(UB_EXARGS,"undefined CREATEOPTION") );

   return outfile.is_open();
}
/*==========================================================*/
void UbFileOutputASCII::writeBool(const bool& value, const int& width)				
{
   outfile.width(width);
   if(value) outfile<<"true ";
   else      outfile<<"false ";
}
/*==========================================================*/
void UbFileOutputASCII::writeDouble(const double& value, const int& width)				
{
   outfile.width(width);
   //Problem: Ub::inf wird gerundet 
   //         -> beim Einlesen ist der Wert evtl zu gross und es kommt murks raus 
   //         -> max Laenge darstellen und gut ist
   if(UbMath::equal(value, (double)Ub::inf) )
   {
      ios_base::fmtflags flags = outfile.flags();
      outfile<<setprecision(std::numeric_limits<double>::digits10+2);
      outfile<<value<<" ";
      outfile.flags(flags);
      return;
   }
   outfile<<value<<" ";
}
/*==========================================================*/
void UbFileOutputASCII::writeFloat(const float& value, const int& width)				
{
   outfile.width(width);
   //Problem: Ub::inf wird gerundet 
   //         -> beim Einlesen ist der Wert evtl zu gross und es kommt murks raus 
   //         -> max Laenge darstellen und gut ist
   if(UbMath::equal(value, (float)Ub::inf) )
   {
      ios_base::fmtflags flags = outfile.flags();
      outfile<<setprecision(std::numeric_limits<float>::digits10+2);
      outfile<<value<<" ";
      outfile.flags(flags);
      return;
   }
   outfile<<value<<" ";
}
/*==========================================================*/
void UbFileOutputASCII::setPrecision(const int& precision)				
{
   outfile<<setprecision(precision);
}
/*==========================================================*/
void UbFileOutputASCII::writeInteger(const int& value, const int& width)				
{
   outfile.width(width);
   outfile<<value<<" ";
}
/*==========================================================*/
void UbFileOutputASCII::writeSize_t(const std::size_t& value, const int& width)
{
   outfile.width(width);
   outfile<<value<<" ";
}
/*==========================================================*/
void UbFileOutputASCII::writeChar(const char& value, const int& width)
{
   outfile.width(width);
   outfile<<(int)value<<" ";	
}
/*==========================================================*/
void UbFileOutputASCII::writeString(const string& value, const int& width)				
{
   outfile.width(width);
   outfile<<value.c_str()<<" ";	
}
/*==========================================================*/
void UbFileOutputASCII::writeStringOnly(const string& value)				
{
	outfile<<value.c_str();	
}

/*==========================================================*/
void UbFileOutputASCII::writeLine(const string& value, const int& width)				
{
   outfile.width(width);
   outfile<<value.c_str()<<endl;	
}
/*==========================================================*/
void UbFileOutputASCII::writeLine()				
{
	outfile<<endl;	
}
/*==========================================================*/
void UbFileOutputASCII::writeCommentLine(const string& line) 
{
   this->writeCommentLine(this->commentindicator, line); 
}
/*==========================================================*/
void UbFileOutputASCII::writeCommentLine(char indicator, const string& line) 
{
	this->outfile<<indicator<<line<<endl;
}
/*==========================================================*/
void UbFileOutputASCII::writeCopyOfFile(const string& filename)
{
   ifstream infile(filename.c_str());
   if(!infile) UB_THROW( UbException(UB_EXARGS,"couldn't open file:\n "+filename) );

   try
   {
      char c;
      while(infile.get(c)) 
      {
         outfile.put(c);  //=out<<c;
      }
      outfile.flush();
      infile.close();
   }
   catch(std::exception& e) { UB_THROW( UbException(UB_EXARGS,"catched std::exception, error: "+(std::string)e.what()) ); }
   catch(...)               { UB_THROW( UbException(UB_EXARGS,"unknown error") ); }
}
