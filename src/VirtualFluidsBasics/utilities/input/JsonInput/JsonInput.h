#ifdef BUILD_JSONCPP

#ifndef JsonInput_H
#define JsonInput_H


#include <string>
#include <json/json.h>


#include "VirtualFluidsBasics_EXPORT.h"


#include "../Input.h"

namespace input
{
	class JsonInput : public Input
	{
	public:
        VirtualFluidsBasics_EXPORT JsonInput(std::istream &stream);

        VirtualFluidsBasics_EXPORT virtual bool hasValue(const std::string &key) const override;
        VirtualFluidsBasics_EXPORT virtual std::string getValue(const std::string &key) override;

    private:
        Json::Value jsonValue;

    };
}

#endif

#endif
