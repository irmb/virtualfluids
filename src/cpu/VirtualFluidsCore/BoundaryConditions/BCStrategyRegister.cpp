#include "BCStrategyRegister.h"
#include "BC.h"
#include "BCStrategy.h"

std::mutex BCStrategyRegister::instantiation_mutex = std::mutex();
std::shared_ptr<BCStrategyRegister> BCStrategyRegister::instance = std::shared_ptr<BCStrategyRegister>();

std::shared_ptr<BCStrategyRegister> BCStrategyRegister::getInstance()
{
    std::lock_guard<std::mutex> myLock(instantiation_mutex);
    if (!instance) {
        instance = std::shared_ptr<BCStrategyRegister>(new BCStrategyRegister);
    }
    return instance;
}

void BCStrategyRegister::setBCStrategy(char strategyKey, std::shared_ptr<BCStrategy> bcStrategy)
{
    bcMap.insert(std::make_pair(strategyKey, bcStrategy));
}

std::shared_ptr<BCStrategy> BCStrategyRegister::getBCStrategy(char strategyKey)
{
    return bcMap[strategyKey];
}