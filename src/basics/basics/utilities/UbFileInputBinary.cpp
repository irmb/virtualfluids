#include <basics/utilities/UbFileInputBinary.h>
#include <cstring>

using namespace std;

/*==========================================================*/
UbFileInputBinary::UbFileInputBinary(string filename)
{
	this->filename = filename;
   infile.open(filename.c_str(), ios::in | ios::binary);
}
/*==========================================================*/
bool UbFileInputBinary::open(string filename)
{
   infile.close();
   infile.clear(); //setzt flags zurueck
   
   this->filename = filename;
   infile.open(this->filename.c_str(), ios::in | ios::binary);

   return infile.is_open();
}
/*==========================================================*/
int UbFileInputBinary::readInteger()				
{
   int dummy;
   infile.read((char*)&dummy,sizeof(int));
   return dummy; 
}
/*==========================================================*/
std::size_t UbFileInputBinary::readSize_t()				
{
   std::size_t dummy;
   infile.read((char*)&dummy,sizeof(std::size_t));
   return dummy;
}
/*==========================================================*/
double UbFileInputBinary::readDouble()	
{
   double dummy;
   infile.read((char*)&dummy,sizeof(double));
   return dummy; 
}
/*==========================================================*/
float UbFileInputBinary::readFloat()	
{
	float dummy;
	infile.read((char*)&dummy,sizeof(float));
	return dummy; 
}
/*==========================================================*/
char UbFileInputBinary::readChar()	
{
   char dummy;
   infile.read((char*)&dummy,sizeof(char));
   return dummy; 
}
/*==========================================================*/
string UbFileInputBinary::readString()	
{
   char c;
   infile.read(&c,sizeof(char));
   while(c==' ' || c=='\t') infile.read(&c,sizeof(char));  
   
   string dummy;
   dummy+=c;

   infile.read(&c,sizeof(char));
   while(c!='\0' && c!=' ' && c!='\t' && c!='\n')
   {
      dummy+=c;
      infile.read(&c,sizeof(char));
   }
   return dummy;
}
/*==========================================================*/
bool UbFileInputBinary::readBool()	
{
   bool dummy;
   infile.read((char*)&dummy,sizeof(bool));
   return dummy; 
}
/*==========================================================*/
void UbFileInputBinary::skipLine()				
{
   char c;
   do{
      infile.read(&c,sizeof(char));
   }while(c!='\n');
}
/*==========================================================*/
void UbFileInputBinary::readLine()				
{
   char c;
   infile.read(&c,sizeof(char));
   while(c!='\n') infile.read(&c,sizeof(char));
}
/*==========================================================*/
string UbFileInputBinary::readStringLine()				
{
   char c;
   string dummy;
   infile.read(&c,sizeof(char));
   while(c!='\n')
   {
      dummy+=c;
      infile.read(&c,sizeof(char));
   }
   return dummy;
}
/*==========================================================*/
string UbFileInputBinary::readLineTill(char  /*stop*/)				
{
   UB_THROW( UbException(UB_EXARGS,"method makes no sense for binary streams") );
}
/*==========================================================*/
string UbFileInputBinary::parseString()				
{
   UB_THROW( UbException(UB_EXARGS,"method makes no sense for binary streams") );
}
/*==========================================================*/
bool UbFileInputBinary::containsString(const string&  /*var*/)
{
   UB_THROW( UbException(UB_EXARGS,"method makes no sense for binary streams") );
}
/*==========================================================*/
void UbFileInputBinary::setPosAfterLineWithString(const string&  /*var*/)
{
   UB_THROW( UbException(UB_EXARGS,"method makes no sense for binary streams") );
}
/*==========================================================*/
int UbFileInputBinary::readIntegerAfterString(const string&  /*var*/)
{
   UB_THROW( UbException(UB_EXARGS,"method makes no sense for binary streams") );
}
/*==========================================================*/
double UbFileInputBinary::readDoubleAfterString(const string&  /*var*/)	
{
   UB_THROW( UbException(UB_EXARGS,"method makes no sense for binary streams") );
}
/*==========================================================*/
string UbFileInputBinary::readStringAfterString(const string&  /*var*/)	
{
   UB_THROW( UbException(UB_EXARGS,"method makes no sense for binary streams") );
}
/*==========================================================*/
bool UbFileInputBinary::readBoolAfterString(const string&  /*var*/)	
{
   UB_THROW( UbException(UB_EXARGS,"method makes no sense for binary streams") );
}
