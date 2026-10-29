//
//  multichain.cpp
//  gctb
//
//  Created by Jian Zeng on 24/11/2024.
//

#include "multichain.hpp"
#include <exception>


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
                if (enabled) out << badSnpIdxVec[j] << "\t" << badSnpNameVec[j] << endl;
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

void MultiChainSBayesR::NumHighPIPs::getValue(const VectorXf &PIP, const unsigned iter,
                                               const set<unsigned> &skepticalSnpIdx){
    const bool isPosterior = iter > burnin;
    if (isPosterior && !posteriorPhase) {
        runningPip.setZero(PIP.size());
        cntSample = 0;
        value = 0;
        posteriorPhase = true;
    }

    const unsigned sampleInterval = isPosterior ? thin : previewThin;
    if (iter % sampleInterval) return;

    if (runningPip.size() != PIP.size()) {
        runningPip.setZero(PIP.size());
        cntSample = 0;
    }
    ++cntSample;
    runningPip.array() += (PIP - runningPip).array()/float(cntSample);

    value = 0;
    for (unsigned i=0; i<(unsigned)runningPip.size(); ++i) {
        if (runningPip[i] <= threshold) continue;
        if (!excludeSkepticalSnps || !skepticalSnpIdx.count(i)) ++value;
    }
}

void MultiChainSBayesR::sampleUnknowns(const unsigned iter){

    Stat::prepareParallelChainRng(numThreadLevel1);
    std::exception_ptr firstException;
#pragma omp parallel for num_threads(numThreadLevel1)
    for (unsigned i=0; i<numChains; ++i) {
//        cout << "sampling chain " << i << " in " << numChains << " chains " << endl;
//        printf("Outer: Thread %d of %d\n", omp_get_thread_num(), omp_get_num_threads());

        // Restrict inner parallelism to numThreadLevel2 threads
        omp_set_num_threads(numThreadLevel2);
        Stat::McmcChainRngScope chainRng(i);

        try {
            chainVec[i]->sampleUnknowns(iter);
        } catch (...) {
#pragma omp critical(gctb_multichain_r_exception)
            {
                if (!firstException) firstException = std::current_exception();
            }
        }
    }
    Stat::finishParallelChainRng();
    if (firstException) std::rethrow_exception(firstException);
    
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
    
    nBadSnps.output();
    nHighPips.getValue(pip.values, iter, nBadSnps.badSnpIdxSet);
}

void MultiChainSBayesRC::NumBadSnps::output(){
    value = 0;
    for (unsigned i=0; i<numChains; ++i) {
        vector<unsigned> badSnpIdxVec = nBadSnpVec[i]->badSnpIdx;
        vector<string> badSnpNameVec = nBadSnpVec[i]->badSnpName;
        for (unsigned j=0; j<badSnpNameVec.size(); ++j) {
            if(badSnpSet.insert(badSnpNameVec[j]).second) {
                badSnpIdxSet.insert(badSnpIdxVec[j]);
                if (enabled) out << badSnpIdxVec[j] << "\t" << badSnpNameVec[j] << endl;
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

void MultiChainSBayesRC::computePredictionAccuracyEnrichment(){
    const unsigned numBlocks = keptLdBlockInfoVec.size();
    const unsigned numSnps = annoMat.rows();
    const float pairDenom = float(numChains * (numChains - 1));

    vector<int> binaryAnnoSlot(annoInfoVec.size(), -1);
    vector<unsigned> binaryAnnoIdx;
    vector<unsigned> quantitativeAnnoIdx;
    for (unsigned k=0; k<annoInfoVec.size(); ++k) {
        if (annoInfoVec[k]->isBinary) {
            if (k) {
                binaryAnnoSlot[k] = binaryAnnoIdx.size();
                binaryAnnoIdx.push_back(k);
            }
        } else {
            quantitativeAnnoIdx.push_back(k);
        }
    }

    VectorXf blockRsq = VectorXf::Zero(numBlocks);
    MatrixXf blockAnnoStats = MatrixXf::Zero(numBlocks, annoInfoVec.size());

#pragma omp parallel for schedule(dynamic) num_threads(numThreadTotal)
    for (unsigned blk=0; blk<numBlocks; ++blk) {
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[blk];
        const QuantizedEigenBlock *Qq = quantizedRsqEnrich ? &quantizedEigenBlocks[blk] : nullptr;
        const MatrixXf *Q = quantizedRsqEnrich ? nullptr : &Qblocks[blk];
        const int rank = Qq ? Qq->k : Q->rows();
        const int numQCols = Qq ? Qq->m : Q->cols();
        const VectorXf *originalNormWeight = nullptr;
        if (contextOriginalNormWeights.size() == numBlocks
            && contextOriginalNormWeights[blk].size() == rank) {
            originalNormWeight = &contextOriginalNormWeights[blk];
        }

        MatrixXf what(rank, numChains);
        for (unsigned c=0; c<numChains; ++c)
            what.col(c) = chainVec[c]->whatBlocks[blk];

        const VectorXf whatSum = what.rowwise().sum();
        const ArrayXf crossProductByRow = whatSum.array().square()
            - what.array().square().rowwise().sum();
        blockRsq[blk] = (originalNormWeight
            ? (crossProductByRow * originalNormWeight->array()).sum()
            : crossProductByRow.sum()) / pairDenom;

        vector<MatrixXf> binaryWhat(binaryAnnoIdx.size());
        for (unsigned a=0; a<binaryAnnoIdx.size(); ++a)
            binaryWhat[a] = MatrixXf::Zero(rank, numChains);

        MatrixXf betaByQCol;
        if (quantitativeAnnoIdx.size())
            betaByQCol = MatrixXf::Zero(numQCols, numChains);

        vector<pair<unsigned, int> > mappedSnps;
        mappedSnps.reserve(blockInfo->snpInfoVec.size());
        for (unsigned t=0; t<blockInfo->snpInfoVec.size(); ++t) {
            SnpInfo *snp = blockInfo->snpInfoVec[t];
            if (!snp || snp->index < 0) continue;
            const unsigned snpIdx = snp->index;
            const int qcol = eigenQColIndex(blockInfo, snpIdx);
            if (qcol < 0 || qcol >= numQCols) continue;
            mappedSnps.push_back(make_pair(snpIdx, qcol));

            for (unsigned c=0; c<numChains; ++c) {
                const float beta = chainVec[c]->snpEffects.values[snpIdx];
                if (quantitativeAnnoIdx.size()) betaByQCol(qcol,c) = beta;
                if (beta == 0.0f) continue;

                for (unsigned a=0; a<snp->annoVec.size(); ++a) {
                    const int annoIdx = snp->annoVec[a]->idx;
                    if (annoIdx < 0 || annoIdx >= (int)binaryAnnoSlot.size()) continue;
                    const int slot = binaryAnnoSlot[annoIdx];
                    if (slot < 0) continue;
                    if (Qq) {
                        auto target = binaryWhat[slot].col(c);
                        Qq->addScaledQ(qcol, beta, target);
                    } else {
                        binaryWhat[slot].col(c).noalias() += Q->col(qcol) * beta;
                    }
                }
            }
        }

        for (unsigned a=0; a<binaryAnnoIdx.size(); ++a) {
            const VectorXf annoWhatSum = binaryWhat[a].rowwise().sum();
            const ArrayXf annoCrossProductByRow = annoWhatSum.array().square()
                - binaryWhat[a].array().square().rowwise().sum();
            blockAnnoStats(blk, binaryAnnoIdx[a]) = (originalNormWeight
                ? (annoCrossProductByRow * originalNormWeight->array()).sum()
                : annoCrossProductByRow.sum()) / pairDenom;
        }

        if (quantitativeAnnoIdx.size()) {
            MatrixXf rBetaByQCol(numQCols, numChains);
            if (Qq) {
                for (int qcol=0; qcol<numQCols; ++qcol) {
                    for (unsigned c=0; c<numChains; ++c)
                        rBetaByQCol(qcol,c) = Qq->dotQ(qcol, chainVec[c]->whatBlocks[blk].data());
                }
            } else {
                if (originalNormWeight) {
                    const MatrixXf originalScaleWhat =
                        (what.array().colwise() * originalNormWeight->array()).matrix();
                    rBetaByQCol.noalias() = Q->transpose() * originalScaleWhat;
                } else {
                    rBetaByQCol.noalias() = Q->transpose() * what;
                }
            }

            const VectorXf betaSum = betaByQCol.rowwise().sum();
            const VectorXf rBetaSum = rBetaByQCol.rowwise().sum();
            const VectorXf withinChain =
                (betaByQCol.array() * rBetaByQCol.array()).rowwise().sum().matrix();
            const VectorXf snpRsq =
                (betaSum.array() * rBetaSum.array() - withinChain.array()).matrix() / pairDenom;

            for (unsigned t=0; t<mappedSnps.size(); ++t) {
                const unsigned snpIdx = mappedSnps[t].first;
                const int qcol = mappedSnps[t].second;
                const float rsq = snpRsq[qcol];
                for (unsigned a=0; a<quantitativeAnnoIdx.size(); ++a) {
                    const unsigned annoIdx = quantitativeAnnoIdx[a];
                    blockAnnoStats(blk, annoIdx) += annoMat(snpIdx, annoIdx) * rsq;
                }
            }
        }
    }

    const float totalRsq = blockRsq.sum();
    const VectorXf annoStats = blockAnnoStats.colwise().sum().transpose();
    annoPerSnpRsqEnrich.values.setZero();

    if (numSnps && totalRsq != 0.0f) {
        for (unsigned k=0; k<annoInfoVec.size(); ++k) {
            AnnoInfo *anno = annoInfoVec[k];
            if (k == 0) {
                annoPerSnpRsqEnrich.values[k] = 1.0f;
            } else if (anno->isBinary) {
                annoPerSnpRsqEnrich.values[k] = anno->fraction == 0.0f
                    ? 0.0f
                    : annoStats[k] / totalRsq / anno->fraction;
            } else if (anno->ssq == 0.0f || anno->sd == 0.0f) {
                annoPerSnpRsqEnrich.values[k] = 1.0f;
            } else {
                Matrix2f XPX;
                XPX(0,0) = float(numSnps);
                XPX(0,1) = XPX(1,0) = anno->sum;
                XPX(1,1) = anno->ssq;
                const Vector2f XPy(totalRsq, annoStats[k]);
                const Vector2f coef = XPX.householderQr().solve(XPy);
                annoPerSnpRsqEnrich.values[k] =
                    1.0f + coef(1) * anno->sd / (totalRsq / float(numSnps));
            }
        }
    } else if (annoPerSnpRsqEnrich.size) {
        annoPerSnpRsqEnrich.values[0] = 1.0f;
    }

    for (unsigned c=0; c<numChains; ++c)
        annoPerSnpRsqEnrich.perChainValues.col(c) = annoPerSnpRsqEnrich.values;
}

void MultiChainSBayesRC::sampleUnknowns(const unsigned iter){

    Stat::prepareParallelChainRng(numThreadLevel1);
    std::exception_ptr firstException;
#pragma omp parallel for num_threads(numThreadLevel1)
    for (unsigned i=0; i<numChains; ++i) {
//        cout << "sampling chain " << i << " in " << numChains << " chains " << endl;
        omp_set_num_threads(numThreadLevel2);
        Stat::McmcChainRngScope chainRng(i);
        chainVec[i]->summaryRequired = summaryRequired;  // set by MCMC::run for this wrapper

        try {
            chainVec[i]->sampleUnknowns(iter);
        } catch (...) {
#pragma omp critical(gctb_multichain_rc_exception)
            {
                if (!firstException) firstException = std::current_exception();
            }
        }
    }
    Stat::finishParallelChainRng();
    if (firstException) std::rethrow_exception(firstException);
        
    snpEffects.getValues();
    pip.getValues();
    deltaPi.getValues();
    hsq.getValues();
    if (contextVarianceRatioE8) {
        contextVarianceRatioE8->getValues();
        contextVarianceE8->getValues();
    }
    if (rcMcHsqFixed) {
        rcMcHsqFixed->getValues();
        rcMcHsqRandom->getValues();
    }
    numSnpMix.getValues();
    vgMix.getValues();
    annoEffects.getValues();
    annoJointProb.getValues();
    annoTotalGenVar.getValues();
    annoPerSnpHsqEnrich.getValues();
    annoConditionalHsqEnrich.getValues();
    annoConditionalPolygenicityEnrich.getValues();
    if (estimateRsqEnrich && iter > rsqEnrichBurnin && !(iter % rsqEnrichFreq))
        computePredictionAccuracyEnrichment();
    snpHsqPep.getValues();

    nBadSnps.output();
    nHighPips.getValue(pip.values, iter, nBadSnps.badSnpIdxSet);
    
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

    Stat::prepareParallelChainRng(numThreadLevel1);
    std::exception_ptr firstException;
#pragma omp parallel for num_threads(numThreadLevel1)
    for (unsigned i=0; i<numModels; ++i) {
        try {
            // Restrict inner parallelism to numThreadLevel2 threads
            omp_set_num_threads(numThreadLevel2);
            Stat::McmcChainRngScope chainRng(i);
            modelVec[i]->sampleUnknowns(iter);
        } catch (...) {
#pragma omp critical(gctb_multimodel_exception)
            {
                if (!firstException) firstException = std::current_exception();
            }
        }
    }
    Stat::finishParallelChainRng();
    if (firstException) std::rethrow_exception(firstException);
}

void MultiChainSBayesRD::NumBadSnps::output(){
    value = 0;
    for (unsigned i=0; i<numChains; ++i) {
        vector<unsigned> badSnpIdxVec = nBadSnpVec[i]->badSnpIdx;
        vector<string> badSnpNameVec = nBadSnpVec[i]->badSnpName;
        for (unsigned j=0; j<badSnpNameVec.size(); ++j) {
            if(badSnpSet.insert(badSnpNameVec[j]).second) {
                badSnpIdxSet.insert(badSnpIdxVec[j]);
                if (enabled) out << badSnpIdxVec[j] << "\t" << badSnpNameVec[j] << endl;
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

    Stat::prepareParallelChainRng(numThreadLevel1);
    std::exception_ptr firstException;
#pragma omp parallel for num_threads(numThreadLevel1)
    for (unsigned i=0; i<numChains; ++i) {
//        cout << "sampling chain " << i << " in " << numChains << " chains " << endl;
        omp_set_num_threads(numThreadLevel2);
        Stat::McmcChainRngScope chainRng(i);
        chainVec[i]->summaryRequired = summaryRequired;  // set by MCMC::run for this wrapper

        try {
            chainVec[i]->sampleUnknowns(iter);
        } catch (...) {
#pragma omp critical(gctb_multichain_rd_exception)
            {
                if (!firstException) firstException = std::current_exception();
            }
        }
    }
    Stat::finishParallelChainRng();
    if (firstException) std::rethrow_exception(firstException);
        
    snpEffects.getValues();
    pip.getValues();
    deltaPi.getValues();
    hsq.getValues();
    numSnpMix.getValues();
    vgMix.getValues();
    annoEffects.getValues();
    annoMixtureWeights.getValues();
    annoDeltas.getValues();
    annoJointProb.getValues();
    annoTotalGenVar.getValues();
    annoPerSnpHsqEnrich.getValues();
    annoConditionalHsqEnrich.getValues();
    annoConditionalPolygenicityEnrich.getValues();
    snpHsqPep.getValues();
    piAnno.getValues();
    piAnnos.getValues();
    annoPip.getValues();
    annoPips.getValues();

    nBadSnps.output();
    nHighPips.getValue(pip.values, iter, nBadSnps.badSnpIdxSet);
}

void MultiChainSBayesS::NumBadSnps::output(){
    value = 0;
    for (unsigned i=0; i<numChains; ++i) {
        vector<unsigned> badSnpIdxVec = nBadSnpVec[i]->badSnpIdx;
        vector<string> badSnpNameVec = nBadSnpVec[i]->badSnpName;
        for (unsigned j=0; j<badSnpNameVec.size(); ++j) {
            if(badSnpSet.insert(badSnpNameVec[j]).second) {
                badSnpIdxSet.insert(badSnpIdxVec[j]);
                if (enabled) out << badSnpIdxVec[j] << "\t" << badSnpNameVec[j] << endl;
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

void MultiChainSBayesS::NumHighPIPs::getValue(const VectorXf &PIP, const unsigned iter,
                                               const set<unsigned> &skepticalSnpIdx){
    const bool isPosterior = iter > burnin;
    if (isPosterior && !posteriorPhase) {
        runningPip.setZero(PIP.size());
        cntSample = 0;
        value = 0;
        posteriorPhase = true;
    }

    const unsigned sampleInterval = isPosterior ? thin : previewThin;
    if (iter % sampleInterval) return;

    if (runningPip.size() != PIP.size()) {
        runningPip.setZero(PIP.size());
        cntSample = 0;
    }
    ++cntSample;
    runningPip.array() += (PIP - runningPip).array()/float(cntSample);

    value = 0;
    for (unsigned i=0; i<(unsigned)runningPip.size(); ++i) {
        if (runningPip[i] <= threshold) continue;
        if (!excludeSkepticalSnps || !skepticalSnpIdx.count(i)) ++value;
    }
}

void MultiChainSBayesS::sampleUnknowns(const unsigned iter){

    Stat::prepareParallelChainRng(numThreadLevel1);
#pragma omp parallel for num_threads(numThreadLevel1)
    for (unsigned i=0; i<numChains; ++i) {
        // Restrict inner parallelism to numThreadLevel2 threads
        omp_set_num_threads(numThreadLevel2);
        Stat::McmcChainRngScope chainRng(i);

        chainVec[i]->sampleUnknowns(iter);
    }
    Stat::finishParallelChainRng();
    
    snpEffects.getValues();
    pip.getValues();
    hsq.getValues();
    S.getValues();
    pi.getValues();
    nnzSnp.getValues();
    sigmaSq.getValues();
    varg.getValues();
    vare.getValues();
    
    nBadSnps.output();
    nHighPips.getValue(pip.values, iter, nBadSnps.badSnpIdxSet);
}
