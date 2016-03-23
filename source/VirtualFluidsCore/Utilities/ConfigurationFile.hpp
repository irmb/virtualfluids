#ifndef Configuration_h__
#define Configuration_h__

#include <map>
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
//!Configuration   config;
//!config.load(configname);
//!
//!int            nbDimensions = config.getInt("nbDimensions");
//!float          temperature  = config.getFloat("temperature");
//!double         epsilon      = config.getDouble("epsilon");
//!bool           writeLogFile = config.getBool("writeLogFile");
//!string         errorMessage = config.getString("pathname");
//!vector<double> origin       = config.getVector<double>("origin");
//!            
//! \author  Konstantin Kutscher

class ConfigurationFile
{
public:
   // clear all values
   void clear();

   // load a configuration file
   bool load(const std::string& File);

   // check if value associated with given key exists
   bool contains(const std::string& key) const;

   // get value associated with given key
   
   int    getInt(const std::string& key) const;
   long   getLong(const std::string& key) const;
   float  getFloat(const std::string& key) const;
   double getDouble(const std::string& key) const;
   bool   getBool(const std::string& key) const;
   std::string  getString(const std::string& key) const;
   template<class T>
   std::vector<T> getVector(const std::string& Key) const;

   template<class T>
   T get(const std::string& Key) const;

private:
   // the container
   std::map<std::string, std::string> data;

   // remove leading and trailing tabs and spaces
   static std::string trim(const std::string& str);

   template<class T>
   T fromString(const std::string& str) const;
};


// ----------------------------------
// method implementations
// ----------------------------------

void ConfigurationFile::clear()
{
   data.clear();
}

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

bool ConfigurationFile::contains(const std::string& key) const
{
   return data.find(key) != data.end();
}

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

int ConfigurationFile::getInt(const std::string& key) const
{
   std::string str = getString(key);
   int value = std::atoi(str.c_str());
   return value;
}

long ConfigurationFile::getLong(const std::string& key) const
{
   std::string str = getString(key);
   long value = std::atol(str.c_str());
   return value;
}

float ConfigurationFile::getFloat(const std::string& key) const
{
   std::string str = getString(key);
   float value = (float)std::atof(str.c_str());
   return value;
}

double ConfigurationFile::getDouble(const std::string& key) const
{
   std::string str = getString(key);
   double value = std::atof(str.c_str());
   return value;
}

bool ConfigurationFile::getBool(const std::string& key) const
{
   std::string str = getString(key);
   bool value = (str == "true");
   return value;
}

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

template<class T>
std::vector<T> ConfigurationFile::getVector(const std::string& key) const
{
   std::string str = getString(key);
   std::vector<T> v;
   std::vector<std::string> strings;
   boost::algorithm::split(strings, str, boost::algorithm::is_any_of("\t\n\r;, "));
   BOOST_FOREACH(std::string s, strings)
   {
      if (s != "")
      {
         v.push_back(fromString<T>(s));
      }
   }
   return v;
}

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

#endif // Configuration_h__
