#ifndef STRINGUTIL_H
#define STRINGUTIL_H
#include <errno.h>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <string>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

template <class T>
bool from_string(T &t, const std::string &s, std::ios_base &(*f)(std::ios_base&))
{
   std::istringstream iss(s);
   return !(iss>>f>>t).fail();
}
class StringUtil 
{
public:
   StringUtil() {}
   ~StringUtil() {}
   // Find the given string in the source string and replace it with the
   // "replace" string, everywhere instances of that string exist.
   static void findandreplace( std::string& source, const std::string& find, const std::string& replace )
   {
      size_t j;
      for (;(j = source.find( find )) 
         != std::string::npos;)
      {
         source.replace( j, 
            find.length(), replace );
      }
   }
   // The following function returns a string with all-uppercase characters.
   static std::string makeUpper( const std::string& instring)
   {
      std::string temp=instring;
      transform( temp.begin(), temp.end(), temp.begin(), ::toupper );
      return temp;
   }
   // The following function returns a string with all-lowercase characters.
   static std::string makeLower( const std::string& instring)
   {
      std::string temp;
      transform( temp.begin(), temp.end(), temp.begin(), ::tolower );
      return temp;
   }
   static bool contains( const std::string& source, const char *find )
   {
      return ( 0!=strstr(source.c_str(),find) );
   }
   static std::string pad( const std::string& instring, char padchar, int length )
   {
      std::string outstring = instring;
      for ( int i=(int)outstring.length(); i<length; ++i )
         outstring += padchar;
      return outstring;
   }
   // Trim the given characters from the beginning and end of a string.
   // the default is to trim whitespace. If the string is empty or contains
   // only the trim characters, an empty string is returned.
   static std::string trim( const std::string &instring, const std::string &trimstring=std::string(" \t\n"))
   {
      if (trimstring.size()==0) 
         return instring;
      std::string temp="";
      std::string::size_type begpos=instring.find_first_not_of (trimstring);
      if (begpos==std::string::npos)
      {
         return temp;
      }
      else
      {
         std::string::size_type endpos=instring.find_last_not_of (trimstring);
         temp=instring.substr(begpos, endpos-begpos+1);
      }
      return temp;
   }
   // Convert the string to an int. Note that a string exception is thrown if
   // it is invalid.
   static int toInt(const std::string & myInString)
   {
      int i=0;
      std::string inString = trim(myInString);
      if( !from_string<int>(i, inString, std::dec) )
      {
         std::string exceptionText = "StringUtils::toInt() - Not an integer: " + inString;
         throw exceptionText;
      }
      // Time to run some more checks.
      for (unsigned int j=0; j < inString.length(); j++)
      {
         if ( !isNumeric(inString[j]) )
         {
            if (j==0 && inString[j] =='-')
            {
               continue;
            }
            else
            {
               std::string exceptionText = "StringUtils::toInt() - Not an integer: " + inString;
               throw exceptionText;
            }
         }
      }
      return (i);
   }
   // Convert the string to a float. Note: A string exception is thrown if
   // it is invalid.
   static float toFloat(const std::string & myInString)
   {
      float f=0;
      std::string inString = trim(myInString);
      if( !from_string<float>(f, inString, std::dec) )
      {
         std::string exceptionText = "StringUtils::toFloat() - Not a float: " + inString;
         throw exceptionText;
      }
      // Now it runs some more checks.
      int dec_count=0;
      int e_count=0;
      for (unsigned int j=0; j < inString.length(); j++)
      {
         if ( !isNumeric(inString[j]) )
         {
            if ((j==0 || inString[j-1] == 'e' || inString[j-1] == 'E') && inString[j] =='-')
            {
               continue;
            }
            else if (inString[j]=='.')
            {
               dec_count++;
               if (dec_count > 1)
               {
                  std::string exceptionText = "StringUtils::toFloat() - Not a float: " + inString;
                  throw exceptionText;
               }
               continue;
            }
            else if (inString[j] == 'e' || inString[j] == 'E')
            {
               e_count++;
               if (e_count > 1)
               {
                  std::string exceptionText = "StringUtils::toFloat() - Not a float: " + inString;
                  throw exceptionText;
               }
               continue;
            }
            else
            {
               std::string exceptionText = "StringUtils::toFloat() - Not a float: " + inString;
               throw exceptionText;
            }
         }
      }
      return (f);
   }
   // Convert the string to a double. Note: A string exception is thrown if
   // it is invalid.
   static double toDouble(const std::string & myInString)
   {
	   double d=0;
	   std::string inString = trim(myInString);
	   if( !from_string<double>(d, inString, std::dec) )
	   {
		   std::string exceptionText = "StringUtils::toDouble() - Not a double: " + inString;
		   throw exceptionText;
	   }
	   // Now it runs some more checks.
      int dec_count=0;
      int e_count=0;
      for (unsigned int j=0; j < inString.length(); j++)
      {
         if ( !isNumeric(inString[j]) )
         {
            if ((j==0 || inString[j-1] == 'e' || inString[j-1] == 'E') && inString[j] =='-')
            {
               continue;
            }
            else if (inString[j]=='.')
            {
               dec_count++;
               if (dec_count > 1)
               {
                  std::string exceptionText = "StringUtils::toDouble() - Not a double: " + inString;
                  throw exceptionText;
               }
               continue;
            }
            else if (inString[j] == 'e' || inString[j] == 'E')
            {
               e_count++;
               if (e_count > 1)
               {
                  std::string exceptionText = "StringUtils::toDouble() - Not a double: " + inString;
                  throw exceptionText;
               }
               continue;
            }
            else
            {
               std::string exceptionText = "StringUtils::toDouble() - Not a double: " + inString;
               throw exceptionText;
            }
         }
      }
	   return (d);
   }
   // Convert the string to a boolean. Note: A string exception is thrown if
   // it is invalid.
   static bool toBool(const std::string & myInString)
   {
      bool b=0;
      std::string inString = trim(myInString);
      if( !from_string<bool>(b, inString, std::boolalpha) )
      {
         std::string exceptionText = "StringUtils::toBool() - Not a bool: " + inString;
         throw exceptionText;
      }
      return (b);
   }
   // Returns true if the character is numeric.
   static bool isNumeric(char c)
   {
      return ('0' <= c && c <= '9');
   }
   // Replace environment variables in the string with their values.
   // Note: environment variables must be of the form ${ENVVAR}.
   //static std::string substituteEnvVar( const std::string &myInString )
   //{
   //   std::string outString="";
   //   char variable[512];
   //   const char *s = myInString.c_str();
   //   while(*s!=0)
   //   {
   //      if (*s=='$' && *(s+1)=='{')
   //      {
   //         // When you�ve found beginning of variable, find the end.
   //         //strcpy(variable,s+2);
   //         strcpy_s(variable, s+2);
   //         char *end = strchr (variable,'}');
   //         if (end)
   //         {
   //            *end='\0';
   //            //char *cp = (char *)getenv(variable);
   //            char *cp;
   //            size_t len;
   //            _dupenv_s(&cp, &len, variable );
   //            if (cp)
   //               //outString += (char *) getenv(variable);
   //               outString += (char *) cp;
   //            s = strchr(s,'}');
   //         }
   //         else
   //         {
   //            outString += *s;
   //         }
   //      }
   //      else
   //      {
   //         outString += *s;
   //      }
   //      s++;
   //   }
   //   return outString;
   //}
   template<class T>
   static std::vector<T> toVector(const std::string& s)
   {
      std::vector<T> v;
      std::vector<std::string> strings;
      boost::algorithm::split(strings, s, boost::is_any_of("\t\n "));
      BOOST_FOREACH(std::string s, strings){
         if (s != "")
         {
            v.push_back(fromString<T>(s));
         }
      }
      return v;
   }
   template<class T>
 	static std::string toString(const T& t)
	{
		std::ostringstream stream;
		stream << t;
		return stream.str();
	}
   template<class T>
   static T fromString(const std::string& s)
   {
      std::istringstream stream (s);
      T t;
      stream >> t;
      return t;
   }
};
#endif //STRINGUTIL_H

