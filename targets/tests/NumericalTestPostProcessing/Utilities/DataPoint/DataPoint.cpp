#include "DataPoint.h"

DataPoint::DataPoint()
{
}

std::shared_ptr<DataPoint> DataPoint::getNewInstance(double x, double y)
{
	return std::shared_ptr<DataPoint>(new DataPoint(x, y));
}

double DataPoint::getX()
{
	return x;
}

double DataPoint::getY()
{
	return y;
}

DataPoint::DataPoint(double x, double y) : x(x), y(y)
{
}
