#include <Incube.h>

Incubator::Incubator()
: mGenSize(50), mMutationRate(0.0015f) {}


Incubator::Incubator(const int& mutationRate, const int& genSize)
: mGenSize(genSize), mMutationRate(mutationRate){}

