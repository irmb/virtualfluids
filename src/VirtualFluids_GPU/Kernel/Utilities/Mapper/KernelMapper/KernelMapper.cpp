#include "KernelMapper.h"

std::shared_ptr<KernelMapper> KernelMapper::getInstance()
{
	static std::shared_ptr<KernelMapper> uniqueInstance;
	if (!uniqueInstance)
		uniqueInstance = std::shared_ptr<KernelMapper>(new KernelMapper());
	return uniqueInstance;
}

std::string KernelMapper::getString(KernelType enumeration)
{
	return myEnumMapper.getString(enumeration);
}

KernelType KernelMapper::getEnum(std::string name)
{
	return myEnumMapper.getEnum(name);
}

KernelMapper::KernelMapper()
{
	myEnumMapper.addEnum(LB_BGKCompSP27,						"BGKCompSP27");
	myEnumMapper.addEnum(LB_BGKPlusCompSP27,					"BGKPlusCompSP27");
	myEnumMapper.addEnum(LB_CascadeCompSP27,					"CascadeCompSP27");
	myEnumMapper.addEnum(LB_CumulantCompSP27,					"CumulantCompSP27");
	myEnumMapper.addEnum(LB_CumulantAA2016CompSP27,				"CumulantAA2016CompSP27");
	myEnumMapper.addEnum(LB_CumulantAA2016CompBulkSP27,			"CumulantAA2016CompBulkSP27");
	myEnumMapper.addEnum(LB_CumulantAll4CompSP27,				"CumulantAll4CompSP27");
	myEnumMapper.addEnum(LB_CumulantF3CompSP27,					"CumulantF3CompSP27");
	myEnumMapper.addEnum(LB_CumulantF32018CompSP27,				"CumulantF32018CompSP27");
	myEnumMapper.addEnum(LB_CumulantOneCompSP27,				"CumulantOneCompSP27");
	myEnumMapper.addEnum(LB_CumulantOneCompBulkSP27,			"CumulantOneCompBulkSP27");
	myEnumMapper.addEnum(LB_CumulantOneCompSpongeSP27,			"CumulantOneCompSpongeSP27");
	myEnumMapper.addEnum(LB_MRTCompSP27,						"MRTCompSP27");
	myEnumMapper.addEnum(LB_BGKIncompSP27,						"BGKIncompSP27");
	myEnumMapper.addEnum(LB_BGKPlusIncompSP27,					"BGKPlusIncompSP27");
	myEnumMapper.addEnum(LB_CascadeIncompSP27,					"CascadeIncompSP27");
	myEnumMapper.addEnum(LB_Cumulant1hIncompSP27,				"Cumulant1hIncompSP27");
	myEnumMapper.addEnum(LB_CumulantIsoIncompSP27,				"CumulantIsoIncompSP27");
	myEnumMapper.addEnum(LB_CumulantOneIncompSP27,				"CumulantOneIncompSP27");
	myEnumMapper.addEnum(LB_MRTIncompSP27,						"MRTIncompSP27");
	myEnumMapper.addEnum(LB_PMCumulantOneCompSP27,				"PMCumulantOneCompSP27");
	myEnumMapper.addEnum(LB_WaleCumulantAA2016CompSP27,			"WaleCumulantAA2016CompSP27");
	myEnumMapper.addEnum(LB_WaleCumulantAA2016DebugCompSP27,	"WaleCumulantAA2016DebugCompSP27");
	myEnumMapper.addEnum(LB_WaleCumulantOneCompSP27,			"WaleCumulantOneCompSP27");
	myEnumMapper.addEnum(LB_WaleBySoniMalavCumulantOneCompSP27,	"WaleBySoniMalavCumulantOneCompSP27");
}