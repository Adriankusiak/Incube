
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

    mutate();
}


template <typename T>
void Incubator<T>::mutate() {
    for(auto sequence : mCurrentGen){
        if(mDist(mRNGen) < mMutationRate ? true : false){
            switch(mMutationType){
                case ALLELE_SWAP:
                    alleleSwap(sequence);
                    break;
                case DESTRUCTIVE:
                    alleleDestructive(sequence);
                    break;
                case GENERATIVE:
                    alleleGenerative(sequence);
                    break;
                default:
                    break;
            }
        }
    }
}

template <typename T>
void Incubator<T>::alleleDestructive(vector<T>& sequence){
    int randIndx = sequence.size()*(floor(mDist(mRNGen))/100);
    sequence.erase(sequence.begin()+randIndx);
}

template <typename T>
void Incubator<T>::alleleGenerative(vector<T>& sequence) {
    int randIndx = sequence.size()*(floor(mDist(mRNGen))/100);
    sequence.erase(sequence.begin()+randIndx);
}

template <typename T>
void Incubator<T>::alleleSwap(vector<T>& sequence) {
    int randIndx = sequence.size()*(floor(mDist(mRNGen))/100);
    int randIndx2 = sequence.size()*(floor(mDist(mRNGen))/100);

    while(randIndx==randIndx2) randIndx2 = sequence.size()*(floor(mDist(mRNGen))/100);
    T allele = sequence[randIndx];


    sequence[randIndx] = sequence[randIndx2];
    sequence[randIndx2] = allele;
}




template <typename T>
void Incubator<T>::singlePoint(const vector<vector<T>>& breeders){

    vector<vector<T>> newGen;
    for(int i = 0; i < mGenSize+2/2; ++i){
        vector<T> sequence;
        vector<T> sequence2;
        int crossIndex = floor(breeders[0].size()*((mDist(mRNGen)/100)));
        sequence.insert(sequence.end(), breeders[i].begin(), breeders[i].begin()+crossIndex);
        sequence.insert(sequence.end(), breeders[i+1].begin()+crossIndex, breeders[i+1].end());
        sequence2.insert(sequence2.end(), breeders[i+1].begin(), breeders[i+1].begin()+crossIndex);
        sequence2.insert(sequence2.end(),  breeders[i].begin()+crossIndex, breeders[i].end());
        newGen.push_back(sequence);
        newGen.push_back(sequence2);
    }

    newGen.resize(mGenSize);
    mCurrentGen = newGen;
}

template <typename T>
void Incubator<T>::twoPoint(const vector<vector<T>>& breeders){

    vector<vector<T>> newGen;
    for(int i = 0; i < mGenSize+4/4; ++i){
        vector<T> sequence;
        vector<T> sequence2;
        vector<T> sequence3;
        vector<T> sequence4;

        int crossIndex = floor(breeders[0].size()*((mDist(mRNGen)/100)));
        int crossIndexTwo = floor(breeders[0].size()*((mDist(mRNGen)/100)));
        if(crossIndex > crossIndexTwo){
            int temp = crossIndex;
            crossIndex = crossIndexTwo;
            crossIndexTwo = temp;
        }
        sequence.insert(sequence.end(), breeders[i].begin(), breeders[i].begin()+crossIndex);
        sequence.insert(sequence.end(), breeders[i+1].begin()+crossIndex, breeders[i+1].end());
        sequence2.insert(sequence2.end(), breeders[i+1].begin(), breeders[i+1].begin()+crossIndex);
        sequence2.insert(sequence2.end(),  breeders[i].begin()+crossIndex, breeders[i].end());

        sequence3.insert(sequence3.end(), sequence.begin(), sequence.begin()+crossIndexTwo);
        sequence3.insert(sequence3.end(), sequence2.begin()+crossIndexTwo, sequence2.end());
        sequence4.insert(sequence4.end(), sequence2.begin(), sequence2.begin()+crossIndexTwo);
        sequence4.insert(sequence4.end(), sequence.begin()+crossIndexTwo, sequence.end());

        newGen.push_back(sequence);
        newGen.push_back(sequence2);
        newGen.push_back(sequence3);
        newGen.push_back(sequence4);
    }

    newGen.resize(mGenSize);
    mCurrentGen = newGen;
}

template <typename T>
void Incubator<T>::setCrossoverType(CrossoverType coType){
    mCrossoverType = coType;
}

template <typename T>
void Incubator<T>::setMutationType(MutationType muType) {
    mMutationType = muType;
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

