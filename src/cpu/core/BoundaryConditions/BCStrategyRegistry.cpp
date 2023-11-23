#include "BCStrategyRegistry.h"
#include "BC.h"
#include "BCStrategy.h"

std::mutex BCStrategyRegistry::instantiation_mutex = std::mutex();
std::shared_ptr<BCStrategyRegistry> BCStrategyRegistry::instance = std::shared_ptr<BCStrategyRegistry>();

std::shared_ptr<BCStrategyRegistry> BCStrategyRegistry::getInstance()
{
    std::lock_guard<std::mutex> myLock(instantiation_mutex);
    if (!instance) {
        instance = std::shared_ptr<BCStrategyRegistry>(new BCStrategyRegistry);
    }
    return instance;
}

void BCStrategyRegistry::setBCStrategy(char strategyKey, std::shared_ptr<BCStrategy> bcStrategy)
{
    bcMap.insert(std::make_pair(strategyKey, bcStrategy));
}

std::shared_ptr<BCStrategy> BCStrategyRegistry::getBCStrategy(char strategyKey)
{
    return bcMap[strategyKey];
}

void BCStrategyRegistry::clear()
{
    bcMap.clear();
}