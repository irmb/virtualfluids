#include "RestartObject.h"
#include "Parameter/Parameter.h"
#include "basics/utilities/UbMath.h"

RestartObject::RestartObject()
{
}
//////////////////////////////////////////////////////////////////////////
RestartObject::~RestartObject()
{
}
//////////////////////////////////////////////////////////////////////////
void RestartObject::load(Parameter* para)
{
	for (int j = para->getCoarse(); j <= para->getFine(); j++)
	{
		std::vector<doubflo> vec;
		fs.push_back(vec);

		for (unsigned int i = 0; i < (para->getD3Qxx()*para->getParH(j)->size_Mat_SP); i++)
		{
			para->getParH(j)->d0SP.f[0][i] = fs[j][i];
		}
	}
}
//////////////////////////////////////////////////////////////////////////
void RestartObject::safe(Parameter* para)
{
	if (fs.size() > 0)
	{
		clear(para);
	}
	for (int j = para->getCoarse(); j <= para->getFine(); j++)
	{
		std::vector<doubflo> vec;
		fs.push_back(vec);

		for (unsigned int i = 0; i < (para->getD3Qxx()*para->getParH(j)->size_Mat_SP); i++)
		{
			if (UbMath::isNaN(para->getParH(j)->d0SP.f[0][i]))
			{
				fs[j].push_back((doubflo)0.0);
			}
			else
			{
				fs[j].push_back(para->getParH(j)->d0SP.f[0][i]);
			}
		}
	}
}
//////////////////////////////////////////////////////////////////////////
void RestartObject::clear(Parameter* para)
{
	for (int j = para->getCoarse(); j <= para->getFine(); j++)
	{
		fs[j].resize(0);
	}
	fs.resize(0);
}
//////////////////////////////////////////////////////////////////////////

