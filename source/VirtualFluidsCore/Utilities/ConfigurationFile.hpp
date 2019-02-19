#ifndef Configuration_h__
#define Configuration_h__

#include <map>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>

//! \brief  Simple configuration file
//! \details The Configuration class presented here can read and keep values of any configuration file written in a format like this:
//!#
//!# Simulation parameters
//!#
//!
//!nbDimensions    = 2
//!temperature     = 25.001
//!epsilon         = 1.013e-14
//!writeLogFile    = false      # NOTE: Set to "true" in debug mode only.
//!                             #       Logging slows down the program.
//!errorMessage    = the simulation failed
//!origin          = 0.0 0.0 0.0 # x, y, z of origin
//!
//!Example how to use it:
//!
//!ConfigurationFile   config;
//!config.load(configname);
//!
//!int            nbDimensions = config.getValue<int>("nbDimensions");
//!float          temperature  = config.getValue<float>("temperature");
//!double         epsilon      = config.getValue<double>("epsilon");
//!bool           writeLogFile = config.getValue<bool>("writeLogFile");
//!string         errorMessage = config.getValue<string>("errorMessage");
//!vector<double> origin       = config.getVector<double>("origin");
//!            
//! \author  Konstantin Kutscher

class ConfigurationFile
{
public:
   //! clear all values
   void clear();

   //! load a configuration file
   bool load(const std::string& File);

   //! check if value associated with given key exists
   bool contains(const std::string& key) const;

   //! get vector with key
   template<class T>
   std::vector<T> getVector(const std::string& key) const;

   //! get value with key
   template<class T>
   T getValue(const std::string& key) const;

private:
   //! the container
   std::map<std::string, std::string> data;

   //! get string with key
   std::string  getString(const std::string& key) const;

   //! remove leading and trailing tabs and spaces
   static std::string trim(const std::string& str);

   //! convert string to data type T
   template<class T>
   T fromString(const std::string& str) const;

   void split(std::vector<std::string>& lst, const std::string& input, const std::string& separators, bool remove_empty = true) const;
};


// ----------------------------------
// method implementations
// ----------------------------------

void ConfigurationFile::clear()
{
   data.clear();
}
//////////////////////////////////////////////////////////////////////////
bool ConfigurationFile::load(const std::string& file)
{
   std::ifstream inFile(file.c_str());

   if (!inFile.good())
   {
      UB_THROW(UbException(UB_EXARGS, "Cannot read configuration file "+file+"!"));
   }

   while (inFile.good() && ! inFile.eof())
   {
      std::string line;
      getline(inFile, line);

      // filter out comments
      if (!line.empty())
      {
         size_t pos = line.find('#');

         if (pos != std::string::npos)
         {
            line = line.substr(0, pos);
         }
      }

      // split line into key and value
      if (!line.empty())
      {
         size_t pos = line.find('=');

         if (pos != std::string::npos)
         {
            std::string key = trim(line.substr(0, pos));
            std::string value = trim(line.substr(pos + 1));

            if (!key.empty() && !value.empty())
            {
               data[key] = value;
            }
         }
      }
   }

   return true;
}
//////////////////////////////////////////////////////////////////////////
bool ConfigurationFile::contains(const std::string& key) const
{
   return data.find(key) != data.end();
}
//////////////////////////////////////////////////////////////////////////
std::string ConfigurationFile::getString(const std::string& key) const
{
   std::map<std::string, std::string>::const_iterator iter = data.find(key);

   if (iter != data.end())
   {
      std::string value = iter->second;
      return value;
   }
   else
   {
      UB_THROW(UbException(UB_EXARGS, "The parameter \"" + key + "\" is missing!"));
   }
}
//////////////////////////////////////////////////////////////////////////
std::string ConfigurationFile::trim(const std::string& str)
{
   size_t first = str.find_first_not_of(" \t\n\r");

   if (first != std::string::npos)
   {
      size_t last = str.find_last_not_of(" \t\n\r");

      return str.substr(first, last - first + 1);
   }
   else
   {
      return "";
   }
}
//////////////////////////////////////////////////////////////////////////
template<class T>
std::vector<T> ConfigurationFile::getVector(const std::string& key) const
{
   std::string str = getString(key);
   std::vector<T> v;
   std::vector<std::string> strings;
   split(strings, str, "\t\n\r;, ");
   for (std::vector<std::string>::iterator it = strings.begin(); it != strings.end(); ++it)
   {
      if (*it != "")
      {
         v.push_back(fromString<T>(*it));
      }
   }
   return v;
}
//////////////////////////////////////////////////////////////////////////
void ConfigurationFile::split(std::vector<std::string>& lst, const std::string& input, const std::string& separators, bool remove_empty) const
{
   std::ostringstream word;
   for (size_t n = 0; n < input.size(); ++n)
   {
      if (std::string::npos == separators.find(input[n]))
         word << input[n];
      else
      {
         if (!word.str().empty() || !remove_empty)
            lst.push_back(word.str());
         word.str("");
      }
   }
   if (!word.str().empty() || !remove_empty)
      lst.push_back(word.str());
}
//////////////////////////////////////////////////////////////////////////
template<class T>
T ConfigurationFile::fromString(const std::string& str) const
{
   //boolean hack
   if (str == "true")
      return true;
   else if (str == "false")
      return false;
   //////////////
   std::istringstream stream(str);
   T t;
   stream >> t;
   return t;
}
//////////////////////////////////////////////////////////////////////////
template<class T>
T ConfigurationFile::getValue(const std::string& key) const
{
   std::string str = getString(key);
   bool bFlag = false;
   if ((std::string)typeid(T).name() == (std::string)typeid(bool).name()) 
   {
      bFlag = true;
   }
      
   std::istringstream iss(str);
   T x;
   iss >> x;
   if (!iss && !bFlag)
      UB_THROW(UbException(UB_EXARGS, " cannot convert \"" + str + "\" to type <" + static_cast<std::string>(typeid(x).name()) + ">"));

   if (bFlag)
   {
      bool value = (str == "true");
      x = value;
   }

   return x;
}
#endif // Configuration_h__
