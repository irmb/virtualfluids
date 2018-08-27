//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBCONVERTER_H 
#define UBCONVERTER_H 

#include <cstdlib> 
#include <ctime> 
#include <cassert> 
#include <string>

/*=========================================================================*/
/*  UBConverter                                                             */
/*                                                                         */
//
// encodes  vals to   e.g. base64
// dencodes vals from e.g. base64
// author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
// version 1.0 - 22.10.2007


class UbConverter
{
public:
   static std::string base64_encode(unsigned char const* , unsigned int len);
   static std::string base64_decode(std::string const& s);

   static inline bool is_base64(const unsigned char& c)
   {
      return (isalnum(c) || (c == '+') || (c == '/'));
   }

protected:
   UbConverter() {}
   ~UbConverter() {}

private:
   UbConverter(const UbConverter&);  // not implemented.
   void operator=(const UbConverter&);  //not implemented.

   static const std::string base64_chars;
};



#endif //UBCONVERTER_H
