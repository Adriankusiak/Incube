#ifndef INCUBE_H
#define INCUBE_H

#include <vector>

using std::vector

template <typename T>
class Imprintable{
public:
    virtual void imprint(const vector<T>& imprints) = 0;
    virtual const vector<T>& getFit() = 0;

    // Add fitness function substitution here
};


template <typename T>
class Seeder{
public:
    virtual vector<T>& getSeeds() const = 0;
};


template <typename T>
class Incubator{

public:

void seed(const Seeder<T>& geneSeeder){}
void generationStep(const Imprintable<T>& tester, const int stepCount){}
vector<T>& getGeneration(){}



private:

};


#endif
