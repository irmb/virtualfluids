//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __         
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |        
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |        
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |        
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____    
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|   
//      \    \  |    |   ________________________________________________________________    
//       \    \ |    |  |  ______________________________________________________________|   
//        \    \|    |  |  |         __          __     __     __     ______      _______    
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)   
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______    
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/   
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file LogWriter.hpp
//! \ingroup Output
//! \author Martin Schoenherr
//=======================================================================================
#ifndef LOGWRITER_H
#define LOGWRITER_H

#include <iostream>
#include <fstream>

////////////////////////////////////////////////////////////////////////////////
//! \class LogWriter manages the generation and output of log files
class LogWriter
{
public:
	//! Class default constructor
   LogWriter()
   {  
      consoleOut = false;
   }
   //! Class 2nd constructor
   //! \param string with path an filename
   LogWriter(std::string fname)
   {
      consoleOut = false;
      this->fname = fname;
   }
   //! \brief sets the path and filename for the log file
   //! \param string with path an filename
   void setName(std::string fname)
   {
      this->fname = fname;
   }
   //! \brief sets the status, if the logging information should be displayed on the console
   //! \param flag is an boolean
   void setConsoleOut(bool flag)
   {
      consoleOut = flag;
   }
   //! \brief clear the log file
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
   //! \brief operator overloading for operator "<<"
   template <typename T>
   LogWriter&  operator << (const T& arg)
   {
      ostr.open(fname.c_str(), std::ios_base::app);
      if (ostr.bad())
      {
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
   //! \property fname holds a string with the path and filename for the log file
   std::string fname;
   //! \property ostr is the output stream
   std::ofstream ostr;
   //! \property consoleOut holds the information if the output should be displayed on the console
   bool consoleOut;
};

#endif	
