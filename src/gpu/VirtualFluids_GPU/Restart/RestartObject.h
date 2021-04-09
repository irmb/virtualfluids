#ifndef RestartObject_H
#define RestartObject_H

#include <memory>
#include <string>
#include <vector>

#include <basics/Core/DataTypes.h>

class Parameter;

class RestartObject
{
public:
    virtual ~RestartObject() = default;

    void deserialize(const std::string &filename, std::shared_ptr<Parameter> para);
    void serialize(const std::string &filename, std::shared_ptr<Parameter> para);

    std::vector<std::vector<real>> fs;

    virtual void serialize_internal(const std::string &filename)   = 0;
    virtual void deserialize_internal(const std::string &filename) = 0;

private:
    void clear(std::shared_ptr<Parameter> para);
};


class ASCIIRestartObject : public RestartObject
{
private:
    virtual void serialize_internal(const std::string &filename);
    virtual void deserialize_internal(const std::string &filename);
};

class BinaryRestartObject : public RestartObject
{
private:
    virtual void serialize_internal(const std::string &filename);
    virtual void deserialize_internal(const std::string &filename);
};

#endif
