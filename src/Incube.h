#ifndef INCUBE_H
#define INCUBE_H

#include <vector>
#include <random>

using std::vector;

template <typename T>
class Imprintable{
public:
    virtual void imprint(const vector<vector<T>>& imprints) = 0;
    virtual const vector<vector<T>>& getFit() = 0;

    // Add fitness function substitution here
};


template <typename T>
class Seeder{
public:
    virtual vector<vector<T>> getSeeds() const = 0;
};


template <typename T>
class Incubator{

public:

    enum CrossoverType{
         SINGLE_POINT,
         TWO_POINT,
    };

    enum MutationType{
        DESTRUCTIVE,
        GENERATIVE,
        ALLELE_SWAP
    };


    Incubator();
    Incubator(const int& mutationRate, const int& genSize, const CrossoverType& coType);

    void seed(const Seeder<T>& geneSeeder);
    void seed(const vector<vector<T>>& seed);
    void evolve(Imprintable<T>& tester, const int stepCount);
    auto getGeneration() const;
    void setGenSize(const int& num);
    void setMutationChance(const float& percentage);

    void setCrossoverType(CrossoverType type);
    void setMutationType(MutationType type);

private:

    void singlePoint(const vector<vector<T>>& breeders);
    void twoPoint(const vector<vector<T>>& breeders);
    void mutate();
    void alleleSwap(vector<T>& sequence);
    void alleleDestructive(vector<T>& sequence);
    void alleleGenerative(vector<T>& sequence);

    vector<vector<T>> mCurrentGen;
    CrossoverType mCrossoverType;
    MutationType mMutationType;
    float mMutationRate;
    int mGenSize;
    std::random_device mSeed;
    std::mt19937 mRNGen;
    std::uniform_real_distribution<float> mDist;

};

#endif
