//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBSTRINGINPUTASCII_H
#define UBSTRINGINPUTASCII_H

#include <fstream>
#include <iostream>
#include <string>

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbFileInput.h>

#include <basics/utilities/UbFileInputASCII.h>

class UbStringInputASCII : public UbFileInputASCII
{                               
public:
	UbStringInputASCII(std::string inputString);
	
	std::string getFileName();				
	void	      skipLine();					   // Springt zur naechsten Zeile

   void        readLine();		 
   std::string readStringLine();				
   std::size_t readSize_t();				
   int		   readInteger();				   // Liest einen int-Wert ein
   long		   readLong();				      // Liest einen long-Wert ein
   double	   readDouble();				   // Liest einen double-Wert ein
	float 	   readFloat();				   // Liest einen float-Wert ein
	bool  	   readBool();				      // Liest einen bool-Wert ein
   char        readChar();                // Liest einen char-Wert ein
   std::string	readString();				   // Liest ein Wort ein
	std::string	readLineTill(char stop);	// Liest gesamte Zeile ein bis zu einem bestimmten Zeichen
	std::string	parseString();	

   bool        containsString(std::string var);
   void        setPosAfterLineWithString(std::string var);
   int		   readIntegerAfterString(std::string var);
   double	   readDoubleAfterString(std::string var);
   bool        readBoolAfterString(std::string var);
   std::string readStringAfterString(std::string var);

   FILETYPE getFileType() { return ASCII; }

private:
	std::istringstream instream;
};


#endif


