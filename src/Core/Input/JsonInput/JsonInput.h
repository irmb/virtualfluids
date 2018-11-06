#ifdef BUILD_JSONCPP

#ifndef JsonInput_H
#define JsonInput_H


#include <string>
#include <json/json.h>


#include <VirtualFluidsDefinitions.h>


#include "../Input.h"

namespace input
{
	class JsonInput : public Input
	{
	public:
        VF_PUBLIC JsonInput(std::istream &stream);

        VF_PUBLIC virtual bool hasValue(const std::string &key) const override;
        VF_PUBLIC virtual std::string getValue(const std::string &key) override;

    private:
        Json::Value jsonValue;

    };
}

#endif

#endif
