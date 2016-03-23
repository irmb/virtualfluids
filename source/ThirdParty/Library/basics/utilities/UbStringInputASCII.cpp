#include <basics/utilities/UbStringInputASCII.h>
#include <cstring>

using namespace std;


UbStringInputASCII::UbStringInputASCII(string inputString) : UbFileInputASCII("")
{
	instream.str(inputString);
	
	
//	this->filename         = filename;
//   this->commentindicator = 'C'; 
//   
//   infile.open(filename.c_str());

}
/*==========================================================*/
int UbStringInputASCII::readInteger()				
{
	int dummy;
	instream>>dummy;
	return dummy;
}
/*==========================================================*/
std::size_t UbStringInputASCII::readSize_t()				
{
   std::size_t dummy;
   instream>>dummy;
   return dummy;
}
/*==========================================================*/
string UbStringInputASCII::getFileName()				
{
	return this->filename;
}

/*==========================================================*/
void UbStringInputASCII::skipLine()				
{
	string dummy;
	getline(instream, dummy);
}
/*==========================================================*/
void UbStringInputASCII::readLine()				
{
	string dummy;
	getline(instream, dummy);
}
/*==========================================================*/
string UbStringInputASCII::readStringLine()				
{
   string dummy;
   getline(instream, dummy);
   return dummy;
}
/*==========================================================*/
string UbStringInputASCII::readLineTill(char stop)				
{
	string dummy;
	getline(instream, dummy, stop);
	return dummy;
}
/*==========================================================*/
string UbStringInputASCII::parseString()				
{
	string dummy;
	getline(instream, dummy, ' ');
	return dummy;
}
/*==========================================================*/
double UbStringInputASCII::readDouble()	
{
   double dummy;
   instream>>dummy;
   return dummy;
}
/*==========================================================*/
float UbStringInputASCII::readFloat()	
{
   float dummy;
   instream>>dummy;
   return dummy;
}
/*==========================================================*/
char UbStringInputASCII::readChar()	
{
   char dummy;
   instream>>dummy;
   return dummy;
}
/*==========================================================*/
string UbStringInputASCII::readString()	
{
	string dummy;
	instream>>dummy;
	return dummy;
}
/*==========================================================*/
bool UbStringInputASCII::containsString(string var)
{
   instream.seekg(0L, ios::beg); //Positionszeiger der Datei auf den Anfang setzen
   char line[512];								
   do
   { 
      instream.getline(line,512);
      if(instream.eof()) return false;
   }while (strstr(line,var.c_str()) != line);		// Ende Schleife, wenn varname ganz in zeile vorkommt	
   
   return true;
}
/*==========================================================*/
void UbStringInputASCII::setPosAfterLineWithString(string var)
{
   instream.seekg(0L, ios::beg); //Positionszeiger der Datei auf den Anfang setzen
   char line[512];								
   do
   { 
      instream.getline(line,512);
      if(instream.eof()) UB_THROW( UbException(UB_EXARGS,var+" wasn't found in "+this->filename) );
   }while (strstr(line,var.c_str()) != line);		// Ende Schleife, wenn varname ganz in zeile vorkommt	
}
/*==========================================================*/
int UbStringInputASCII::readIntegerAfterString(string var)
// last change [10.3.2004] at [9:46] 
//suchts in einer Datei nach varname und gibt den dahinter stehenden int-Wert zurueck
//z.B. timesteps 9
{
   instream.seekg(0L, ios::beg); //Positionszeiger der Datei auf den Anfang setzen

   char line[512];								

   do
   { 
      instream.getline(line,512);
      if(instream.eof()) UB_THROW( UbException(UB_EXARGS,var+" wasn't found in "+this->filename) );
   }while (strstr(line,var.c_str()) != line);		// Ende Schleife, wenn varname ganz in zeile vorkommt	

   strcpy (line, (line+strlen(var.c_str())));	    // zeile um "varname" kuerzen 
   while ((line[0] == ' ') || (line[0] == '\t')) strcpy (line, (line+1));	// Whitespaces entfernen

   return(atoi(line));						// Umwandlung in int 					
}
/*==========================================================*/
// last change [10.3.2004] at [9:46] 
//sucht in einer Datei nach varname und gibt den dahinter stehenden int-Wert zurueck
//z.B. nue 9.5
double UbStringInputASCII::readDoubleAfterString(string var)	
{
   instream.seekg(0L, ios::beg); //Positionszeiger der Datei auf den Anfang setzen

   char line[512];								

   do
   { 
      instream.getline(line,512);
      if(instream.eof()) UB_THROW( UbException(UB_EXARGS,var+" wasn't found in "+this->filename) );
   }while (/*!strncmp(varname,line,sizeof(varname))==0*/strstr(line,var.c_str()) != line);		// Ende Schleife, wenn varname ganz in zeile vorkommt	


   strcpy (line, (line+strlen(var.c_str())));	    // zeile um "varname" kuerzen 
   while ((line[0] == ' ') || (line[0] == '\t')) strcpy (line, (line+1));	// Whitespaces entfernen

   return (atof(line));			// Umwandlung in double 					
}
/*==========================================================*/
//  [9.9.2002]
// liefert sring-Wert der hinter dem uebergebenen char feld in der datei instream steht
// zudem wird der wert in die uebergebene variable value uebertragen (falls man das ergebniss als char benoetig)
string UbStringInputASCII::readStringAfterString(string var)	
{
   instream.seekg(0L, ios::beg); //Positionszeiger der Datei auf den Anfang setzen

   char line[512];								
   //string line_copy[512];

   do{ 
      instream.getline(line,512);
      if(instream.eof()) UB_THROW( UbException(UB_EXARGS,var+" wasn't found in "+this->filename) );
   }while (strstr(line,var.c_str()) != line);		// Ende Schleife, wenn varname ganz in zeile vorkommt	

   strcpy (line, (line+strlen(var.c_str())));										// zeile um "varname" kuerzen 
   while ((line[0] == ' ') || (line[0] == '\t')) strcpy (line, (line+1));	// Whitespaces entfernen

   char *p;
   p=strtok(line," "); //schneidet alles "ab und inklusive space " nach namen ab
   p=strtok(line,"\t");//schneidet alles "ab und inklusive tab   " nach namen ab

   return (string)p;			// Umwandlung in string					
}
/*==========================================================*/
// last change [10.3.2004] at [9:46] 
//sucht in einer Datei nach varname und gibt den dahinter stehenden int-Wert zurueck
//z.B. nue 9.5
bool UbStringInputASCII::readBoolAfterString(string var)	
{
   if(this->readStringAfterString(var.c_str())      == "true" ) return true;
   else if(this->readStringAfterString(var.c_str()) == "false") return false;
   else UB_THROW( UbException(UB_EXARGS,var+" is not equal to 'true' or 'false' in "+this->filename) );
}
/*==========================================================*/
// last change [10.3.2004] at [9:46] 
//sucht in einer Datei nach varname und gibt den dahinter stehenden int-Wert zurueck
//z.B. nue 9.5
bool UbStringInputASCII::readBool()	
{
   string tmp = this->readString();
   if(     tmp == "true" ) return true;
   else if(tmp == "false") return false;
   else UB_THROW( UbException(UB_EXARGS,"expression is not equal to 'true' or 'false' in "+this->filename) );
}
