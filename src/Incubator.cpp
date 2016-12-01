
#include <math.h>
#include <algorithm>
#include "Incube.h"


using std::floor;

template <typename T>
Incubator<T>::Incubator()
: mGenSize(50), mMutationRate(0.0015f), mSeed(), mRNGen(mSeed), mDist(0,100),
  mCrossoverType(Incubator::SINGLE_POINT) {}

template <typename T>
Incubator<T>::Incubator(const int& mutationRate, const int& genSize, const CrossoverType& coType)
: mGenSize(genSize), mMutationRate(mutationRate), mSeed(), mRNGen(mSeed), mDist(0,100), mCrossoverType(coType) {}

template <typename T>
void Incubator<T>::seed(const Seeder<T>& geneSeeder){
    mCurrentGen = geneSeeder.getSeeds();
    mGenSize  = mCurrentGen.size();
}

template <typename T>
void Incubator<T>::seed(const vector<vector<T>>& seed){
    mCurrentGen = seed;
    mGenSize  = mCurrentGen.size();
}

template <typename T>
void Incubator<T>::evolve(Imprintable<T>& tester, const int stepCount){
    tester.imprint(mCurrentGen);
    auto fit = tester.getFit();
    int newSize = mCurrentGen.size();
    if(fit.size() < newSize){
        for(auto f: fit){
            mCurrentGen.erase(std::find(mCurrentGen.begin(), mCurrentGen.end(), f));
        }
        fit.reserve(newSize);
        fit.insert(fit.end(), mCurrentGen.begin(), mCurrentGen.end());
    }
    switch(mCrossoverType){
        case SINGLE_POINT:
            singlePoint(fit);
            break;
        case TWO_POINT:
            twoPoint(fit);
            break;
    }
}

template <typename T>
void Incubator<T>::singlePoint(const vector<vector<T>>& breeders){

    vector<vector<T>> newGen;
    for(int i = 0; i < mGenSize+2/2; ++i){
        vector<T> specimen;
        vector<T> specimen2;
        int crossIndex = floor(breeders[0].size()*((mDist(mRNGen)/100)));
        specimen.insert(specimen.end(), breeders[i].begin(), breeders[i].begin()+crossIndex);
        specimen.insert(specimen.end(), breeders[i+1].begin()+crossIndex, breeders[i+1].end());
        specimen2.insert(specimen2.end(), breeders[i+1].begin(), breeders[i+1].begin()+crossIndex);
        specimen2.insert(specimen2.end(),  breeders[i].begin()+crossIndex, breeders[i].end());
        newGen.push_back(specimen);
        newGen.push_back(specimen2);
    }

    newGen.resize(mGenSize);
    mCurrentGen = newGen;
}

template <typename T>
void Incubator<T>::twoPoint(const vector<vector<T>>& breeders){

    vector<vector<T>> newGen;
    for(int i = 0; i < mGenSize+4/4; ++i){
        vector<T> specimen;
        vector<T> specimen2;
        vector<T> specimen3;
        vector<T> specimen4;

        int crossIndex = floor(breeders[0].size()*((mDist(mRNGen)/100)));
        int crossIndexTwo = floor(breeders[0].size()*((mDist(mRNGen)/100)));
        if(crossIndex > crossIndexTwo){
            int temp = crossIndex;
            crossIndex = crossIndexTwo;
            crossIndexTwo = temp;
        }
        specimen.insert(specimen.end(), breeders[i].begin(), breeders[i].begin()+crossIndex);
        specimen.insert(specimen.end(), breeders[i+1].begin()+crossIndex, breeders[i+1].end());
        specimen2.insert(specimen2.end(), breeders[i+1].begin(), breeders[i+1].begin()+crossIndex);
        specimen2.insert(specimen2.end(),  breeders[i].begin()+crossIndex, breeders[i].end());

        specimen3.insert(specimen3.end(), specimen.begin(), specimen.begin()+crossIndexTwo);
        specimen3.insert(specimen3.end(), specimen2.begin()+crossIndexTwo, specimen2.end());
        specimen4.insert(specimen4.end(), specimen2.begin(), specimen2.begin()+crossIndexTwo);
        specimen4.insert(specimen4.end(), specimen.begin()+crossIndexTwo, specimen.end());

        specimen.resize(breeders.size());
        specimen2.resize(breeders.size());
        specimen3.resize(breeders.size());
        specimen4.resize(breeders.size());

        newGen.push_back(specimen);
        newGen.push_back(specimen2);
        newGen.push_back(specimen3);
        newGen.push_back(specimen4);
    }

    newGen.resize(mGenSize);
    mCurrentGen = newGen;
}

template <typename T>
void Incubator<T>::setCrossoverType(CrossoverType coType){
    mCrossoverType = coType;
}

template <typename T>
auto Incubator<T>::getGeneration() const{
    return mCurrentGen;
}

template <typename T>
void Incubator<T>::setGenSize(const int& num){
    mGenSize = num;
}

template <typename T>
void Incubator<T>::setMutationChance(const float& percentage){
    mMutationRate = percentage;
}

