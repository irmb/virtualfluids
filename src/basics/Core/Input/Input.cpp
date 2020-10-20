#include "Input.h"

#include <memory>

#ifdef BUILD_JSONCPP
#include "JsonInput/JsonInput.h"
#endif

#include "ConfigInput/ConfigInput.h"

namespace input
{

std::unique_ptr<input::Input> Input::makeInput(std::istream &stream, const std::string & /*inputType*/)
{
#ifdef BUILD_JSONCPP
    if (inputType == "json")
        return std::unique_ptr<Input>(new JsonInput(stream));
#endif

    return std::make_unique<ConfigInput>(stream);
}

} // namespace input
