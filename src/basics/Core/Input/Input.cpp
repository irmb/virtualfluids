#include "Input.h"

#ifdef BUILD_JSONCPP
#include "JsonInput/JsonInput.h"
#endif

#include "ConfigInput/ConfigInput.h"

namespace input
{

    std::unique_ptr<input::Input> Input::makeInput(std::istream &stream, const std::string & /*inputType*/)
    {
#ifdef BUILD_JSONCPP
        if(inputType == "json")
            return std::unique_ptr<Input>(new JsonInput(stream));
#endif
         
        // changed by St. Lenz: make_unique<...> is C++ 14 standard!
        //return std::make_unique<ConfigInput>(stream);
        return std::unique_ptr<ConfigInput>(new ConfigInput(stream));
    }

}

