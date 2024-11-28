//
//  multichains.cpp
//  gctb
//
//  Created by Jian Zeng on 24/11/2024.
//

#include "multichains.hpp"


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
    for (unsigned i=0; i<numChains; ++i) {
        vector<unsigned> badSnpIdxVec = nBadSnpVec[i]->badSnpIdx;
        vector<string> badSnpNameVec = nBadSnpVec[i]->badSnpName;
        for (unsigned j=0; j<badSnpNameVec.size(); ++j) {
            if(badSnpSet.insert(badSnpNameVec[j]).second) {
                out << badSnpIdxVec[j] << "\t" << badSnpNameVec[j] << endl;
            }
        }
    }
    value = badSnpSet.size();
}

void MultiChainSBayesR::sampleUnknowns(){
    // Set the maximum number of nested parallelism levels
//    omp_set_max_active_levels(2);
    
//#pragma omp parallel num_threads(numChains)
    for (unsigned i=0; i<numChains; ++i) {
//        cout << "sampling chain " << i << " in " << numChains << " chains " << endl;
        chainVec[i]->sampleUnknowns();
    }
    
    snpEffects.getValues();
    pip.getValues();
    deltaPi.getValues();
    hsq.getValues();
    numSnpMix.getValues();
    vgMix.getValues();
    
    nBadSnps.output();
}

void MultiChainSBayesRC::NumBadSnps::output(){
    for (unsigned i=0; i<numChains; ++i) {
        vector<unsigned> badSnpIdxVec = nBadSnpVec[i]->badSnpIdx;
        vector<string> badSnpNameVec = nBadSnpVec[i]->badSnpName;
        for (unsigned j=0; j<badSnpNameVec.size(); ++j) {
            if(badSnpSet.insert(badSnpNameVec[j]).second) {
                out << badSnpIdxVec[j] << "\t" << badSnpNameVec[j] << endl;
            }
        }
    }
    value = badSnpSet.size();
}

void MultiChainSBayesRC::sampleUnknowns(){
    // Set the maximum number of nested parallelism levels
//    omp_set_max_active_levels(2);
    
//#pragma omp parallel num_threads(numChains)
    for (unsigned i=0; i<numChains; ++i) {
//        cout << "sampling chain " << i << " in " << numChains << " chains " << endl;
        chainVec[i]->sampleUnknowns();
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
    
    nBadSnps.output();
}
