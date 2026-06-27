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
                badSnpIdxSet.insert(badSnpIdxVec[j]);
                out << badSnpIdxVec[j] << "\t" << badSnpNameVec[j] << endl;
                ++value;
            }
        }
    }
    //value = badSnpSet.size();
    
    for (unsigned i=0; i<numChains; ++i) {
        nBadSnpVec[i]->badSnpIdx.resize(badSnpSet.size());
        nBadSnpVec[i]->badSnpName.resize(badSnpSet.size());
        set<string>::iterator it1;
        unsigned j=0;
        for (it1 = badSnpSet.begin(); it1 != badSnpSet.end(); ++it1) {
            nBadSnpVec[i]->badSnpName[j++] = *it1;
        }
        set<unsigned>::iterator it2;
        j=0;
        for (it2 = badSnpIdxSet.begin(); it2 != badSnpIdxSet.end(); ++it2) {
            nBadSnpVec[i]->badSnpIdx[j++] = *it2;
        }
    }
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
    if (mcHsqFixed) {
        mcHsqFixed->getValues();
        mcHsqRandom->getValues();
    }
    numSnpMix.getValues();
    vgMix.getValues();
    snpHsqPep.getValues();
    
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
                badSnpIdxSet.insert(badSnpIdxVec[j]);
                out << badSnpIdxVec[j] << "\t" << badSnpNameVec[j] << endl;
                ++value;
            }
        }
    }
    //value = badSnpSet.size();
    
    for (unsigned i=0; i<numChains; ++i) {
        nBadSnpVec[i]->badSnpIdx.resize(badSnpSet.size());
        nBadSnpVec[i]->badSnpName.resize(badSnpSet.size());
        set<string>::iterator it1;
        unsigned j=0;
        for (it1 = badSnpSet.begin(); it1 != badSnpSet.end(); ++it1) {
            nBadSnpVec[i]->badSnpName[j++] = *it1;
        }
        set<unsigned>::iterator it2;
        j=0;
        for (it2 = badSnpIdxSet.begin(); it2 != badSnpIdxSet.end(); ++it2) {
            nBadSnpVec[i]->badSnpIdx[j++] = *it2;
        }
    }
}

void MultiChainSBayesRC::SnpPIP::computeGelmanRubinStat(const unsigned iter){
    GelmanRubinStat.resize(0);
    if (iter % 10) return;
    ++cntSample;
    for (unsigned i=0; i<numChains; ++i) {
         perChainMean.col(i).array()    += (perChainValues.col(i) - perChainMean.col(i)).array()/cntSample;
         perChainSqrMean.col(i).array() += (perChainValues.col(i).array().square() - perChainSqrMean.col(i).array())/cntSample;
     }
    posteriorMean.array() += (values - posteriorMean).array()/cntSample;


    if (cntSample > 10) {
        VectorXf meanVarWithinChain(size);
        VectorXf varMeansBetweenChains(size);
        VectorXf posteriorVar(size);
        GelmanRubinStat.setZero(size);
        for (unsigned i=0; i<size; ++i) {
            meanVarWithinChain[i] = (perChainSqrMean.row(i).array() - perChainMean.row(i).array().square()).mean(); // W
            varMeansBetweenChains[i] = float(cntSample)*Gadget::calcVariance(perChainMean.row(i).transpose());                           // B
            posteriorVar[i] = (cntSample-1.0)*meanVarWithinChain[i]/float(cntSample) + varMeansBetweenChains[i]/float(cntSample);
            if (meanVarWithinChain[i]) {
                GelmanRubinStat[i] = sqrt(posteriorVar[i]/meanVarWithinChain[i]);
            } else {
                GelmanRubinStat[i] = 1.0;
            }
        }
    }
}

void MultiChainSBayesRC::SnpPIP::selectUnconvergedSnps(const vector<LDBlockInfo*> &keptLdBlockInfoVec){
    selectedSnpForTGS.resize(0);
    if (!GelmanRubinStat.size()) return;
    long nBlocks = keptLdBlockInfoVec.size();
    for(unsigned blk = 0; blk < nBlocks; blk++){
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[blk];
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd   = blockInfo->endSnpIdx;
        unsigned blockSize  = blockEnd - blockStart + 1;
        for(unsigned i = 0; i < blockSize; i++){
            unsigned snpIdx = i + blockStart;
            if (GelmanRubinStat[snpIdx] > 1.2 && posteriorMean[snpIdx] > 0.1) {
                vector <int> vec(3);
                vec[0] = blockInfo->chrom;
                vec[1] = blk;
                vec[2] = snpIdx;
                selectedSnpForTGS.push_back(vec);
            }
        }
    }
}

void MultiChainSBayesRC::AnnoPerSnpRsqEnrichment::compute(const VectorXf &snpCrossChainEffectProduct, const MatrixXf &annoMat, const vector<AnnoInfo*> &annoInfoVec){
    VectorXf rsqVec;
    rsqVec.setZero(size);
    
    for (unsigned k=0; k<size; ++k) {
        AnnoInfo *anno = annoInfoVec[k];
        for (unsigned j=0; j<anno->size; ++j) {
            SnpInfo *snp = anno->memberSnpVec[j];
            rsqVec[k] += annoMat(snp->index,k) * snpCrossChainEffectProduct[snp->index];
        }
        if (rsqVec[0] == 0.0f) {
            values[k] = 0.0f;
        } else if (anno->isBinary) {
            values[k] = rsqVec[k]/rsqVec[0] / anno->fraction;
        } else {
            MatrixXf XPX(2,2);
            XPX(0,0) = size;
            XPX(0,1) = XPX(1,0) = anno->sum;
            XPX(1,1) = anno->ssq;
            VectorXf XPy(2);
            XPy(0) = rsqVec[0];
            XPy(1) = rsqVec[k];
            VectorXf coef = XPX.householderQr().solve(XPy);
            values[k] = 1.0f + coef(1) / rsqVec[0];
        }
    }
}

void MultiChainSBayesRC::AnnoJointPerSnpRsqEnrichment::compute(const VectorXf &snpCrossChainEffectProduct, const MatrixXf &annoMat, const AnnoJointProb &annoJointProb, const VectorXf &gamma, const VectorXf &snpAnnoCntInv){
    unsigned numAnno = annoMat.cols();
    unsigned numSnps = annoMat.rows();
    unsigned numDist = annoJointProb.size();
    VectorXf rsqAnnoObs;
    VectorXf rsqAnnoNull;
    rsqAnnoObs.setZero(numAnno);
    rsqAnnoNull.setZero(numAnno);
    
    float betaSqNull = snpCrossChainEffectProduct.sum()/float(numSnps);
    VectorXf probAnnoObs(numAnno);
    for (unsigned j=0; j<numSnps; ++j) {
        probAnnoObs.setZero(numAnno);
        for (unsigned i=0; i<numAnno; ++i) {
            if (annoMat(j,i)) {
                for (unsigned k=1; k<numDist; ++k) {
                    probAnnoObs[i] += (*annoJointProb[k]).values[i] * gamma[k];
                }
            }
        }
        float probSum = probAnnoObs.sum();
        if (probSum == 0.0f) continue;
        probAnnoObs /= probSum;
        float betajSq = snpCrossChainEffectProduct[j];
        for (unsigned i=0; i<numAnno; ++i) {
            if (annoMat(j,i)) {
                rsqAnnoObs[i]  += probAnnoObs[i]  * betajSq;
                rsqAnnoNull[i] += snpAnnoCntInv[j] * betaSqNull;
            }
        }
    }
    
    values.setZero(numAnno);
    for (unsigned i=0; i<numAnno; ++i) {
        values[i] = rsqAnnoNull[i] == 0.0f ? 0.0f : rsqAnnoObs[i]/rsqAnnoNull[i];
    }
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
    if (rcMcHsqFixed) {
        rcMcHsqFixed->getValues();
        rcMcHsqRandom->getValues();
    }
    numSnpMix.getValues();
    vgMix.getValues();
    annoEffects.getValues();
    annoJointProb.getValues();
    annoTotalGenVar.getValues();
    annoJointPerSnpHsqEnrich.getValues();
    annoPerSnpHsqEnrich.getValues();
    if (estimateRsqEnrich) {
        VectorXf chainSum = snpEffects.perChainValues.rowwise().sum();
        VectorXf chainSumSq = snpEffects.perChainValues.array().square().rowwise().sum();
        VectorXf snpCrossChainEffectProduct = (chainSum.array().square() - chainSumSq.array()) / float(numChains * (numChains - 1));
        annoPerSnpRsqEnrich.compute(snpCrossChainEffectProduct, annoMat, annoInfoVec);
        annoJointPerSnpRsqEnrich.compute(snpCrossChainEffectProduct, annoMat, annoJointProb, gamma.values, snpAnnoCntInv);
    }
    snpHsqPep.getValues();

    nHighPips.getValue(pip.values);
    nBadSnps.output();
    
//    pip.computeGelmanRubinStat(iter);
//    pip.selectUnconvergedSnps(keptLdBlockInfoVec);
//    if (pip.selectedSnpForTGS.size()) {
//        for (unsigned i=0; i<numChains; ++i) {
//            chainVec[i]->sampleUnknownsTGS(pip.selectedSnpForTGS);
//        }
//        snpEffects.getValues();
//        pip.getValues();
//    }
}


void MultiModelSBayesR::sampleUnknowns(const unsigned iter){
    omp_set_max_active_levels(2);

#pragma omp parallel for num_threads(numThreadLevel1)
    for (unsigned i=0; i<numModels; ++i) {
        //cout << "sampling chain " << i << " in " << numChains << " chains " << endl;
        //printf("Outer: Thread %d of %d\n", omp_get_thread_num(), omp_get_num_threads());

        // Restrict inner parallelism to numThreadLevel2 threads
        omp_set_num_threads(numThreadLevel2);

        modelVec[i]->sampleUnknowns(iter);
    }    
}

void MultiChainSBayesRD::NumBadSnps::output(){
    value = 0;
    for (unsigned i=0; i<numChains; ++i) {
        vector<unsigned> badSnpIdxVec = nBadSnpVec[i]->badSnpIdx;
        vector<string> badSnpNameVec = nBadSnpVec[i]->badSnpName;
        for (unsigned j=0; j<badSnpNameVec.size(); ++j) {
            if(badSnpSet.insert(badSnpNameVec[j]).second) {
                badSnpIdxSet.insert(badSnpIdxVec[j]);
                out << badSnpIdxVec[j] << "\t" << badSnpNameVec[j] << endl;
                ++value;
            }
        }
    }
    //value = badSnpSet.size();
    
    for (unsigned i=0; i<numChains; ++i) {
        nBadSnpVec[i]->badSnpIdx.resize(badSnpSet.size());
        nBadSnpVec[i]->badSnpName.resize(badSnpSet.size());
        set<string>::iterator it1;
        unsigned j=0;
        for (it1 = badSnpSet.begin(); it1 != badSnpSet.end(); ++it1) {
            nBadSnpVec[i]->badSnpName[j++] = *it1;
        }
        set<unsigned>::iterator it2;
        j=0;
        for (it2 = badSnpIdxSet.begin(); it2 != badSnpIdxSet.end(); ++it2) {
            nBadSnpVec[i]->badSnpIdx[j++] = *it2;
        }
    }
}

void MultiChainSBayesRD::sampleUnknowns(const unsigned iter){
    
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
    annoJointPerSnpHsqEnrich.getValues();
    snpHsqPep.getValues();
    piAnno.getValues();
    piAnnos.getValues();
    annoPip.getValues();
    annoPips.getValues();

    nHighPips.getValue(pip.values);
    nBadSnps.output();
}

void MultiChainSBayesS::NumBadSnps::output(){
    value = 0;
    for (unsigned i=0; i<numChains; ++i) {
        vector<unsigned> badSnpIdxVec = nBadSnpVec[i]->badSnpIdx;
        vector<string> badSnpNameVec = nBadSnpVec[i]->badSnpName;
        for (unsigned j=0; j<badSnpNameVec.size(); ++j) {
            if(badSnpSet.insert(badSnpNameVec[j]).second) {
                badSnpIdxSet.insert(badSnpIdxVec[j]);
                out << badSnpIdxVec[j] << "\t" << badSnpNameVec[j] << endl;
                ++value;
            }
        }
    }
    
    for (unsigned i=0; i<numChains; ++i) {
        nBadSnpVec[i]->badSnpIdx.resize(badSnpSet.size());
        nBadSnpVec[i]->badSnpName.resize(badSnpSet.size());
        set<string>::iterator it1;
        unsigned j=0;
        for (it1 = badSnpSet.begin(); it1 != badSnpSet.end(); ++it1) {
            nBadSnpVec[i]->badSnpName[j++] = *it1;
        }
        set<unsigned>::iterator it2;
        j=0;
        for (it2 = badSnpIdxSet.begin(); it2 != badSnpIdxSet.end(); ++it2) {
            nBadSnpVec[i]->badSnpIdx[j++] = *it2;
        }
    }
}

void MultiChainSBayesS::NumHighPIPs::getValue(const VectorXf &PIP){
    value = 0;
    unsigned size = PIP.size();
    for (unsigned i=0; i<size; ++i) {
        if (PIP[i] > threshold) ++value;
    }
}

void MultiChainSBayesS::sampleUnknowns(const unsigned iter){

#pragma omp parallel for num_threads(numThreadLevel1)
    for (unsigned i=0; i<numChains; ++i) {
        // Restrict inner parallelism to numThreadLevel2 threads
        omp_set_num_threads(numThreadLevel2);

        chainVec[i]->sampleUnknowns(iter);
    }
    
    snpEffects.getValues();
    pip.getValues();
    hsq.getValues();
    S.getValues();
    pi.getValues();
    nnzSnp.getValues();
    sigmaSq.getValues();
    varg.getValues();
    vare.getValues();
    
    nHighPips.getValue(pip.values);
    nBadSnps.output();
}
