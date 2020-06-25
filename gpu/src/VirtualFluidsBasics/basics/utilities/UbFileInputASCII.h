//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBFILEINPUTASCII_H
#define UBFILEINPUTASCII_H

#include <fstream>
#include <iostream>
#include <string>

#include <VirtualFluidsDefinitions.h>

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbFileInput.h>

/*=========================================================================*/
/*  UbFileInputASCII                                                       */
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

class VF_PUBLIC UbFileInputASCII : public UbFileInput
{                               
public:
   UbFileInputASCII() : UbFileInput() { }
   UbFileInputASCII(std::string filename);
	
   bool open(std::string filename);

   std::string getFileName();				
	void	      skipLine();					   // Springt zur naechsten Zeile

   void        readLine();		 
   std::string readStringLine();				
	int		   readInteger();				   // Liest einen Int-Wert ein
    long long   readLongLong();                       // Liest einen long-Wert ein

   std::size_t readSize_t();
   double	   readDouble();				   // Liest einen double-Wert ein
	float 	   readFloat();				   // Liest einen float-Wert ein
	bool  	   readBool();				      // Liest einen bool-Wert ein
   char        readChar();                // Liest einen char-Wert ein
   std::string	readString();				   // Liest ein Wort ein
	std::string	readLineTill(char stop);	// Liest gesamte Zeile ein bis zu einem bestimmten Zeichen
	std::string	parseString();	


   bool        containsString(const std::string& var);
   void        setPosAfterLineWithString(const std::string& var);
   int		   readIntegerAfterString(const std::string& var);
   double	   readDoubleAfterString(const std::string& var);
   bool        readBoolAfterString(const std::string& var);
   std::string readStringAfterString(const std::string& var);

   FILETYPE getFileType() { return ASCII; }

   template< typename T >
   friend inline UbFileInputASCII& operator>>(UbFileInputASCII& file, T& data) 
   {
      file.infile>>data;
      return file;
   }
};

#endif //UBFILEINPUTASCII_H


