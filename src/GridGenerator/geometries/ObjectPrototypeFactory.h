//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef OBJECT_PROTOTYPE_FACTORY_H
#define OBJECT_PROTOTYPE_FACTORY_H

#include "PointerDefinitions.h"

#include "Sphere/Sphere.h"
#include "Cuboid/Cuboid.h"

enum class ObjectType
{
    Cuboid
};


class AbstractFactory
{
public:
    virtual Object* make(ObjectType type) = 0;
};

class VF_PUBLIC ObjectPrototypeFactory : public AbstractFactory
{
public:
    ~ObjectPrototypeFactory() = default;
    static SPtr<ObjectPrototypeFactory> make()
    {
        return SPtr<ObjectPrototypeFactory>(new ObjectPrototypeFactory());
    }

    Object* make(ObjectType type)
    {
        switch (type)
        {
        case ObjectType::Cuboid:
            return cuboid.clone();
        //case ObjectType::Sphere:
        //    return sphere.clone();
        }
    }

    void setCuboid(Cuboid cuboid)
    {
        this->cuboid = cuboid;
    }


private:
    ObjectPrototypeFactory() = default;
    ObjectPrototypeFactory(const ObjectPrototypeFactory& a) = default;


    Cuboid cuboid;
};




#endif
