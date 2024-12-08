//
//  multichain.cpp
//  gctb
//
//  Created by Jian Zeng on 24/11/2024.
//

#include "multichain.hpp"


void MultiChainParameter::getValues(){
    for (unsigned i=0; i<numChains; ++i) {
        perChainValue[i] = chainVec[i]->value;
    }
    value = perChainValue.mean();
}

void MultiChainParamSet::getValues(){
    for (unsigned i=0; i<numChains; ++i) {
        perChainValues.col(i) = chainVec[i]->values;
    }
    values = perChainValues.rowwise().mean();
}

void MultiChainParamVec::getValues(){
    for (unsigned i=0; i<numParams; ++i) {
        (*this)[i]->getValues();
    }
}

void MultiChainParamSetVec::getValues(){
    for (unsigned i=0; i<numParams; ++i) {
        (*this)[i]->getValues();
    }
}

void MultiChainSBayesR::NumBadSnps::output(){
    value = 0;
    for (unsigned i=0; i<numChains; ++i) {
        vector<unsigned> badSnpIdxVec = nBadSnpVec[i]->badSnpIdx;
        vector<string> badSnpNameVec = nBadSnpVec[i]->badSnpName;
        for (unsigned j=0; j<badSnpNameVec.size(); ++j) {
            if(badSnpSet.insert(badSnpNameVec[j]).second) {
                out << badSnpIdxVec[j] << "\t" << badSnpNameVec[j] << endl;
                ++value;
            }
        }
    }
    //value = badSnpSet.size();
}

void MultiChainSBayesR::NumHighPIPs::getValue(const VectorXf &PIP){
    value = 0;
    unsigned size = PIP.size();
    for (unsigned i=0; i<size; ++i) {
        //cout << i << " " << PIP[i] << " " << threshold << endl;
        if (PIP[i] > threshold) ++value;
    }
}

void MultiChainSBayesR::sampleUnknowns(const unsigned iter){

#pragma omp parallel for num_threads(numThreadLevel1)
    for (unsigned i=0; i<numChains; ++i) {
//        cout << "sampling chain " << i << " in " << numChains << " chains " << endl;
//        printf("Outer: Thread %d of %d\n", omp_get_thread_num(), omp_get_num_threads());

        // Restrict inner parallelism to numThreadLevel2 threads
        omp_set_num_threads(numThreadLevel2);

        chainVec[i]->sampleUnknowns(iter);
    }
    
    snpEffects.getValues();
    pip.getValues();
    deltaPi.getValues();
    hsq.getValues();
    numSnpMix.getValues();
    vgMix.getValues();
    
    nHighPips.getValue(pip.values);
    nBadSnps.output();
}

void MultiChainSBayesRC::NumBadSnps::output(){
    value = 0;
    for (unsigned i=0; i<numChains; ++i) {
        vector<unsigned> badSnpIdxVec = nBadSnpVec[i]->badSnpIdx;
        vector<string> badSnpNameVec = nBadSnpVec[i]->badSnpName;
        for (unsigned j=0; j<badSnpNameVec.size(); ++j) {
            if(badSnpSet.insert(badSnpNameVec[j]).second) {
                out << badSnpIdxVec[j] << "\t" << badSnpNameVec[j] << endl;
                ++value;
            }
        }
    }
    //value = badSnpSet.size();
}

void MultiChainSBayesRC::sampleUnknowns(const unsigned iter){
    
#pragma omp parallel for num_threads(numThreadLevel1)
    for (unsigned i=0; i<numChains; ++i) {
//        cout << "sampling chain " << i << " in " << numChains << " chains " << endl;
        omp_set_num_threads(numThreadLevel2);

        chainVec[i]->sampleUnknowns(iter);
    }
        
    snpEffects.getValues();
    pip.getValues();
    deltaPi.getValues();
    hsq.getValues();
    numSnpMix.getValues();
    vgMix.getValues();
    annoEffects.getValues();
    annoJointProb.getValues();
    annoTotalGenVar.getValues();
    annoPerSnpHsqEnrich.getValues();
    
    nHighPips.getValue(pip.values);
    nBadSnps.output();
}


void MultiModelSBayesR::sampleUnknowns(const unsigned iter){
        
#pragma omp parallel for num_threads(numThreadLevel1)
    for (unsigned i=0; i<numModels; ++i) {
        //cout << "sampling chain " << i << " in " << numChains << " chains " << endl;
        //printf("Outer: Thread %d of %d\n", omp_get_thread_num(), omp_get_num_threads());

        // Restrict inner parallelism to numThreadLevel2 threads
        omp_set_num_threads(numThreadLevel2);

        modelVec[i]->sampleUnknowns(iter);
    }    
}
