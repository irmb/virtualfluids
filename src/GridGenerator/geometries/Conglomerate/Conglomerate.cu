#include "Conglomerate.h"
#include "utilities/logger/Logger.h"


Conglomerate::Conglomerate()
{
    objects = new Object*[MAX_NUMBER_OF_OBJECTS];
}

Conglomerate::~Conglomerate()
{

}

void Conglomerate::add(Object* object)
{
    if (numberOfObjects < MAX_NUMBER_OF_OBJECTS)
    {
        objects[numberOfObjects] = object;
        numberOfObjects++;
    } else
        *logging::out << logging::Logger::WARNING << "max numbers of " << MAX_NUMBER_OF_OBJECTS << "reached! Object was not added.\n";
}

Object* Conglomerate::clone() const
{
    return new Conglomerate();
}

double Conglomerate::getX1Centroid()
{
    return 0.0;
}

double Conglomerate::getX1Minimum()
{
    double minimum = objects[0]->getX1Minimum();
    for(int i = 1; i < numberOfObjects; i++)
    {
        const double compareMinimum = objects[i]->getX1Minimum();
        if (minimum > compareMinimum)
            minimum = compareMinimum;
    }
    return minimum;
}

double Conglomerate::getX1Maximum()
{
    return 0.0;
}

double Conglomerate::getX2Centroid()
{
    return 0.0;
}

double Conglomerate::getX2Minimum()
{
    return 0.0;
}	

double Conglomerate::getX2Maximum()
{
    return 0.0;
}

double Conglomerate::getX3Centroid()
{
    return 0.0;
}

double Conglomerate::getX3Minimum()
{
    return 0.0;
}	

double Conglomerate::getX3Maximum()
{
    return 0.0;
}


bool Conglomerate::isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset)
{

    return false;
}



void Conglomerate::scale(double delta)
{

}
