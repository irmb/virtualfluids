#include "Conglomerate.h"


Conglomerate::Conglomerate()
{
    addObjects = new Object*[MAX_NUMBER_OF_OBJECTS];
    subtractObjects = new Object*[MAX_NUMBER_OF_OBJECTS];
}

Conglomerate::~Conglomerate()
{
    for (uint i = 0; i < numberOfAddObjects; i++)
        delete addObjects[i];

    for (uint i = 0; i < numberOfSubtractObjects; i++)
        delete subtractObjects[i];

    delete[] addObjects;
    delete[] subtractObjects;
}

SPtr<Conglomerate> Conglomerate::makeShared()
{
    return SPtr<Conglomerate>(new Conglomerate());
}

void Conglomerate::add(Object* object)
{
    if (numberOfAddObjects < MAX_NUMBER_OF_OBJECTS)
    {
        addObjects[numberOfAddObjects] = object;
        numberOfAddObjects++;
    } else
        printf("[WARNING] max numbers of %d reached! Object was not added.\n", MAX_NUMBER_OF_OBJECTS);
}

void Conglomerate::subtract(Object* object)
{
    if (numberOfSubtractObjects < MAX_NUMBER_OF_OBJECTS)
    {
        subtractObjects[numberOfSubtractObjects] = object;
        numberOfSubtractObjects++;
    }
    else
        printf("[WARNING] max numbers of %d reached! Object was not added.\n", MAX_NUMBER_OF_OBJECTS);
}

Object* Conglomerate::clone() const
{
    auto conglomerate = new Conglomerate();
    for (uint i = 0; i < numberOfAddObjects; i++)
        conglomerate->add(addObjects[i]->clone());

    for (uint i = 0; i < numberOfSubtractObjects; i++)
        conglomerate->subtract(subtractObjects[i]->clone());

    return conglomerate;
}

double Conglomerate::getX1Centroid()
{
    return (getX1Minimum() + getX1Maximum()) * 0.5;
}

double Conglomerate::getX1Minimum()
{
    double minimum = addObjects[0]->getX1Minimum();
    for(uint i = 1; i < numberOfAddObjects; i++)
        minimum = getMinimum(minimum, addObjects[i]->getX1Minimum());
    return minimum;
}

double Conglomerate::getX1Maximum()
{
    double maximum = addObjects[0]->getX1Maximum();
    for (uint i = 1; i < numberOfAddObjects; i++)
        maximum = getMaximum(maximum, addObjects[i]->getX1Maximum());
    return maximum;
}

double Conglomerate::getX2Centroid()
{
    return (getX2Minimum() + getX2Maximum()) * 0.5;
}

double Conglomerate::getX2Minimum()
{
    double minimum = addObjects[0]->getX2Minimum();
    for (uint i = 1; i < numberOfAddObjects; i++)
        minimum = getMinimum(minimum, addObjects[i]->getX2Minimum());
    return minimum;
}	

double Conglomerate::getX2Maximum()
{
    double maximum = addObjects[0]->getX2Maximum();
    for (uint i = 1; i < numberOfAddObjects; i++)
        maximum = getMaximum(maximum, addObjects[i]->getX2Maximum());
    return maximum;
}

double Conglomerate::getX3Centroid()
{
    return (getX3Minimum() + getX3Maximum()) * 0.5;
}

double Conglomerate::getX3Minimum()
{
    double minimum = addObjects[0]->getX3Minimum();
    for (uint i = 1; i < numberOfAddObjects; i++)
        minimum = getMinimum(minimum, addObjects[i]->getX3Minimum());
    return minimum;
}	

double Conglomerate::getX3Maximum()
{
    double maximum = addObjects[0]->getX3Maximum();
    for (uint i = 1; i < numberOfAddObjects; i++)
        maximum = getMaximum(maximum, addObjects[i]->getX3Maximum());
    return maximum;
}

double Conglomerate::getMinimum(double val1, double val2)
{
    if (val1 > val2)
        return val2;
    return val1;
}

double Conglomerate::getMaximum(double val1, double val2)
{
    if (val1 < val2)
        return val2;
    return val1;
}

bool Conglomerate::isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset)
{
    for (uint i = 0; i < numberOfSubtractObjects; i++)
        if (subtractObjects[i]->isPointInObject(x1, x2, x3, minOffset, maxOffset))
            return false;

    for (uint i = 0; i < numberOfAddObjects; i++)
        if (addObjects[i]->isPointInObject(x1, x2, x3, minOffset, maxOffset))
            return true;

    return false;
}

void Conglomerate::scale(double delta)
{
    for (uint i = 0; i < numberOfAddObjects; i++)
        addObjects[i]->scale(delta);

    for (uint i = 0; i < numberOfSubtractObjects; i++)
        subtractObjects[i]->scale(delta);
}
