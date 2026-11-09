//
//  multichain.cpp
//  gctb
//
//  Created by Jian Zeng on 24/11/2024.
//

#include "multichain.hpp"
#include <exception>
#include <algorithm>
#include <limits>
#include <memory>
#include <random>
#if defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#endif


void MultiChainParameter::getValues(){
    for (unsigned i=0; i<numChains; ++i) {
        perChainValue[i] = chainVec[i]->value;
    }
    value = perChainValue.mean();
}

void MultiChainParamSet::getValues(){
    // The per-SNP sets (millions of rows) in parallel over row ranges, with the same arithmetic per row.
    const long rows = perChainValues.rows(), rangeRows = 1 << 16;
    bool sizesMatch = rows > rangeRows;
    for (unsigned i=0; i<numChains && sizesMatch; ++i) sizesMatch = chainVec[i]->values.size() == rows;
    if (sizesMatch) {
        values.resize(rows);
        #pragma omp parallel for schedule(static)
        for (long first = 0; first < rows; first += rangeRows) {
            const long n = std::min(rangeRows, rows - first);
            for (unsigned i=0; i<numChains; ++i)
                perChainValues.col(i).segment(first, n) = chainVec[i]->values.segment(first, n);
            values.segment(first, n) = perChainValues.middleRows(first, n).rowwise().mean();
        }
        return;
    }
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

    // --batch-chains (changes the sampling, multichain.hpp): all chains' SNP effects in one sweep that reads each LD
    // block once; each chain then runs the rest of its iteration without its own sweep.
    if (batchChains && batchSupported < 0) {
        batchSupported = multichain::snpSweepSupported(multichain::sweepChains(chainVec));
        cout << "Batched multi-chain SNP sweep: " << (batchSupported ? "on" : "not supported for this run, off") << endl;
    }
    const bool batchedSweep = batchChains && batchSupported == 1;
    if (batchedSweep) {
        multichain::sampleSnpEffects(multichain::sweepChains(chainVec), numThreadTotal, batchSweepCache);
        for (ApproxBayesR *chain : chainVec) chain->snpSweepDone = true;
    }

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
    if (batchedSweep) for (ApproxBayesR *chain : chainVec) chain->snpSweepDone = false;
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

void MultiChainSBayesRC::sampleChains(const vector<ApproxBayesRC*> &chains, const unsigned iter){

    // --batch-chains (changes the sampling, multichain.hpp): all chains' SNP effects in one sweep that reads each LD
    // block once; each chain then runs the rest of its iteration without its own sweep.
    if (batchChains && batchSupported < 0) {
        batchSupported = multichain::snpSweepSupported(chains);
        cout << "Batched multi-chain SNP sweep: " << (batchSupported ? "on" : "not supported for this run, off") << endl;
    }
    const bool batchedSweep = batchChains && batchSupported == 1;
    if (batchedSweep) {
        // GCTB_PROFILE_RC: the LD-swap step and the sweep go into every chain's profile.
        const bool profile = chains[0]->profileTiming;
        const double t0 = profile ? gctbWallTime() : 0.0;
        if (chains[0]->algorithm != ApproxBayesR::gibbs) {  // TGS / TGS_thin: each chain's LD-swap step first
            std::exception_ptr swapException;
            Stat::prepareParallelChainRng(numThreadLevel1);
#pragma omp parallel for num_threads(numThreadLevel1)
            for (unsigned i=0; i<numChains; ++i) {
                omp_set_num_threads(numThreadLevel2);
                Stat::McmcChainRngScope chainRng(i);
                try {
                    chains[i]->sampleLdSwap(iter);
                } catch (...) {
#pragma omp critical(gctb_multichain_rc_swap_exception)
                    {
                        if (!swapException) swapException = std::current_exception();
                    }
                }
            }
            Stat::finishParallelChainRng();
            if (swapException) std::rethrow_exception(swapException);
        }
        const double t1 = profile ? gctbWallTime() : 0.0;
        multichain::sampleSnpEffects(multichain::sweepChains(chains), numThreadTotal, batchSweepCache);
        if (profile) {
            const double t2 = gctbWallTime();
            for (ApproxBayesRC *chain : chains) {
                chain->addSharedProfileTime(chain->profileTimeTgs, t1 - t0);
                chain->addSharedProfileTime(chain->profileTimeSnpEffects, t2 - t1);
            }
        }
        for (ApproxBayesRC *chain : chains) chain->snpSweepDone = true;
    }

    // Shared annotation passes (multichain.hpp): each chain's iteration runs in two stages around its annotation
    // sampler, with its random stream carried from stage 1 to stage 2, and the passes after each stage.
    const bool sharedAnno = multichain::sharedAnnoPasses(chains);
    vector<Stat::random_engine> chainEngine(sharedAnno ? numChains : 0);
    vector<Stat::normal_distribution> chainNormal(sharedAnno ? numChains : 0);
    std::exception_ptr firstException;
    auto runChains = [&](const int stage) {
        Stat::prepareParallelChainRng(numThreadLevel1);
#pragma omp parallel for num_threads(numThreadLevel1)
        for (unsigned i=0; i<numChains; ++i) {
//            cout << "sampling chain " << i << " in " << numChains << " chains " << endl;
            omp_set_num_threads(numThreadLevel2);
            Stat::McmcChainRngScope chainRng(i);
            if (stage == 2) {  // continue the chain's random stream from stage 1
                Stat::engine = chainEngine[i];
                Stat::snorm.distribution() = chainNormal[i];
            }
            chains[i]->summaryRequired = summaryRequired;  // set by MCMC::run for this wrapper
            chains[i]->annoPassStage = stage;

            try {
                chains[i]->sampleUnknowns(iter);
            } catch (...) {
#pragma omp critical(gctb_multichain_rc_exception)
                {
                    if (!firstException) firstException = std::current_exception();
                }
            }
            chains[i]->annoPassStage = 0;
            if (stage == 1) {
                chainEngine[i] = Stat::engine;
                chainNormal[i] = Stat::snorm.distribution();
            }
        }
        Stat::finishParallelChainRng();
    };
    if (!sharedAnno) {
        runChains(0);
    } else {
        // GCTB_PROFILE_RC: the shared passes go into every chain's profile, which is reported after them.
        const bool profile = chains[0]->profileTiming;
        double scores = 0.0, probabilities = 0.0, logPi = 0.0;
        runChains(1);
        if (!firstException) {
            const double t0 = profile ? gctbWallTime() : 0.0;
            multichain::annoFirstScores(chains);
            if (profile) scores = gctbWallTime() - t0;
            runChains(2);
        }
        if (!firstException) {
            const double t0 = profile ? gctbWallTime() : 0.0;
            multichain::annoProbabilities(chains);
            const double t1 = profile ? gctbWallTime() : 0.0;
            for (ApproxBayesRC *chain : chains) ApproxBayesRC::logPiFromP(chain->snpP, chain->snpLogPi);  // each in parallel
            if (profile) {
                probabilities = t1 - t0;
                logPi = gctbWallTime() - t1;
            }
        }
        if (profile && !firstException) {
            for (ApproxBayesRC *chain : chains) {
                chain->addSharedProfileTime(chain->profileTimeAnnoEffects, scores + probabilities);
                chain->addSharedProfileTime(chain->profileTimeLogPi, logPi);
                chain->reportProfile(iter);
            }
        }
    }
    if (batchedSweep) for (ApproxBayesRC *chain : chains) chain->snpSweepDone = false;
    if (firstException) std::rethrow_exception(firstException);
}

void MultiChainSBayesRC::sampleUnknowns(const unsigned iter){
    sampleChains(chainVec, iter);
        
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

    // --batch-chains (changes the sampling, multichain.hpp): all models' SNP effects in one sweep that reads each LD
    // block once; each model then runs the rest of its iteration without its own sweep.
    if (batchChains && batchSupported < 0) {
        batchSupported = multichain::snpSweepSupported(multichain::sweepChains(modelVec));
        cout << "Batched multi-model SNP sweep: " << (batchSupported ? "on" : "not supported for this run, off") << endl;
    }
    const bool batchedSweep = batchChains && batchSupported == 1;
    if (batchedSweep) {
        multichain::sampleSnpEffects(multichain::sweepChains(modelVec), numThreadTotal, batchSweepCache);
        for (ApproxBayesR *model : modelVec) model->snpSweepDone = true;
    }

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
    if (batchedSweep) for (ApproxBayesR *model : modelVec) model->snpSweepDone = false;
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
    // The chain loop of SBayesRC (with --batch-chains the batched SNP sweep); SBayesRD's chains have no annotation
    // pass stages (ApproxBayesRD::annoPassStagesSupported), so each chain runs its whole iteration.
    sampleChains(vector<ApproxBayesRC*>(chainVec.begin(), chainVec.end()), iter);
        
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


// ---------------------------------------------------------------------------------------------------------------------
// Work shared by the chains of a multi-chain run (namespace multichain, multichain.hpp).

bool multichain::sharedAnnoPasses(const vector<ApproxBayesRC*> &chains) {
    return chains.size() > 1 && chains[0]->annoPassStagesSupported() && chains[0]->estimatePi
        && !chains[0]->annoEffects.blockedGaussian;
}

void multichain::annoFirstScores(const vector<ApproxBayesRC*> &chains) {
    vector<const VectorXf*> latents;
    vector<VectorXd*> scores;
    for (ApproxBayesRC *chain : chains) {
        latents.push_back(&chain->annoFirstLatents);
        scores.push_back(&chain->annoFirstScores);
    }
    annotationScores(chains[0]->data.annoMat, latents, scores);
    for (ApproxBayesRC *chain : chains) chain->annoFirstLatents.resize(0);
}

void multichain::annoProbabilities(const vector<ApproxBayesRC*> &chains) {
    const unsigned numChains = chains.size();
    vector<MatrixXf> alpha(numChains);
    vector<const MatrixXf*> alphaPtrs;
    vector<MatrixXf*> predictorMeans, snpP;
    for (unsigned c=0; c<numChains; ++c) {
        ApproxBayesRC::AnnoEffects &effects = chains[c]->annoEffects;
        alpha[c].resize(effects.numAnno, effects.numComp);
        for (unsigned i=0; i<effects.numComp; ++i) alpha[c].col(i) = effects[i]->values;
        alphaPtrs.push_back(&alpha[c]);
        predictorMeans.push_back(&effects.predictorMeans);
        snpP.push_back(&chains[c]->snpP);
    }
    annotationProbabilities(chains[0]->annoEffects, chains[0]->data.annoMat, alphaPtrs, predictorMeans, snpP);
}

// ---------------------------------------------------------------------------------------------------------------------
// --batch-chains: the SNP effects of all chains in one sweep that reads each LD block's Q once. Unlike the annotation
// passes above, this changes the sampling: within a block the chains visit the same random order of chunks of 16
// consecutive Q columns, each chain visiting the chunk's SNPs in its own random order, with random numbers from streams
// keyed by chain, iteration and block. A chunk's terms q_j'w are computed for all chains in one product and, after each
// change, corrected with a column of the chunk's Gram matrix Q_chunk'Q_chunk (computed once per run), so each is q_j'w
// at the time of the update; the per-SNP sampling is that of SBayesR's and SBayesRC's sampleFromFC_eigen
// (ApproxBayesR::SnpEffects::sampleSnp).

namespace {
const unsigned sweepChunkSize = 16;

#if defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__)) && !defined(GCTB_BATCH_PORTABLE)
// R = Qc' W (Qc k x count and W k x nc, column-major with column stride k; R(j, h) at R[h * ldr + j]) by AVX2 + FMA,
// compiled for AVX2 by a target attribute on this function alone (GCTB's compiler flags stay as they are). Per group of
// four chains (the last repeated when fewer remain), the rows in tiles of 512 so that a tile of W stays in L1 cache while
// each pair of Qc's columns (the last repeated when count is odd) is multiplied with it, each pair's eight vector sums
// carrying over from tile to tile; then the rows after the last whole vector one by one. It sums in a different order
// from Eigen, so the rounding differs. count must be at most 16.
__attribute__((target("avx2,fma")))
void chunkProductAvx2(const float *Qc, const int k, const int count, const float *W, const int nc, float *R,
                      const int ldr) {
    const int tileRows = 512, vectorRows = k / 8 * 8, numPairs = (count + 1) / 2;
    __m256 acc[8][8];
    for (int h0 = 0; h0 < nc; h0 += 4) {
        const int nh = std::min(4, nc - h0);
        const float *w[4];
        for (int h = 0; h < 4; ++h) w[h] = W + (size_t)(h0 + std::min(h, nh - 1)) * k;
        for (int p = 0; p < numPairs; ++p) for (int i = 0; i < 8; ++i) acc[p][i] = _mm256_setzero_ps();
        for (int tile = 0; tile < vectorRows; tile += tileRows) {
            const int tileEnd = std::min(tile + tileRows, vectorRows);
            for (int p = 0; p < numPairs; ++p) {
                const float *q0 = Qc + (size_t)(2 * p) * k, *q1 = Qc + (size_t)std::min(2 * p + 1, count - 1) * k;
                __m256 a00 = acc[p][0], a01 = acc[p][1], a02 = acc[p][2], a03 = acc[p][3];
                __m256 a10 = acc[p][4], a11 = acc[p][5], a12 = acc[p][6], a13 = acc[p][7];
                for (int r = tile; r < tileEnd; r += 8) {
                    const __m256 x0 = _mm256_loadu_ps(q0 + r), x1 = _mm256_loadu_ps(q1 + r);
                    __m256 b = _mm256_loadu_ps(w[0] + r);
                    a00 = _mm256_fmadd_ps(x0, b, a00); a10 = _mm256_fmadd_ps(x1, b, a10);
                    b = _mm256_loadu_ps(w[1] + r);
                    a01 = _mm256_fmadd_ps(x0, b, a01); a11 = _mm256_fmadd_ps(x1, b, a11);
                    b = _mm256_loadu_ps(w[2] + r);
                    a02 = _mm256_fmadd_ps(x0, b, a02); a12 = _mm256_fmadd_ps(x1, b, a12);
                    b = _mm256_loadu_ps(w[3] + r);
                    a03 = _mm256_fmadd_ps(x0, b, a03); a13 = _mm256_fmadd_ps(x1, b, a13);
                }
                acc[p][0] = a00; acc[p][1] = a01; acc[p][2] = a02; acc[p][3] = a03;
                acc[p][4] = a10; acc[p][5] = a11; acc[p][6] = a12; acc[p][7] = a13;
            }
        }
        for (int j = 0; j < count; ++j)
            for (int h = 0; h < nh; ++h) {
                const __m256 v = acc[j / 2][(j % 2) * 4 + h];  // the eight lanes added pairwise, then the tail rows
                const __m128 v4 = _mm_add_ps(_mm256_castps256_ps128(v), _mm256_extractf128_ps(v, 1));
                const __m128 v2 = _mm_add_ps(v4, _mm_movehl_ps(v4, v4));
                float s = _mm_cvtss_f32(_mm_add_ss(v2, _mm_shuffle_ps(v2, v2, 1)));
                const float *qj = Qc + (size_t)j * k;
                for (int r = vectorRows; r < k; ++r) s += qj[r] * w[h][r];
                R[(size_t)(h0 + h) * ldr + j] = s;
            }
    }
}

// Whether the CPU has AVX2 and FMA; checked once, with a line in the log.
bool chunkKernelAvx2() {
    static const bool avx2 = [] {
        __builtin_cpu_init();
        const bool yes = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
        std::cout << "Batched multi-chain products: " << (yes ? "AVX2 kernel" : "Eigen") << std::endl;
        return yes;
    }();
    return avx2;
}
#endif

// R's first Qc.cols() rows = Qc' W: the AVX2 kernel where available, else Eigen's product (everywhere when built
// with -DGCTB_BATCH_PORTABLE).
void chunkProduct(const Map<const MatrixXf> &Qc, const MatrixXf &W, MatrixXf &R) {
#if defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__)) && !defined(GCTB_BATCH_PORTABLE)
    if (chunkKernelAvx2()) {
        chunkProductAvx2(Qc.data(), Qc.rows(), Qc.cols(), W.data(), W.cols(), R.data(), R.rows());
        return;
    }
#endif
    R.topRows(Qc.cols()).noalias() = Qc.transpose() * W;
}
}  // namespace

multichain::SweepChain multichain::sweepChain(ApproxBayesR &chain) {
    return {&chain, &chain.snpEffects, &chain.snpEffects.snpset, nullptr, &chain.snpStore, nullptr, nullptr};
}

multichain::SweepChain multichain::sweepChain(ApproxBayesRC &chain) {
    return {&chain, &chain.snpEffects, &chain.snpEffects.snpset, &chain.snpEffects.numSnpMix, nullptr, &chain.snpLogPi,
            &chain.snpEffects.fcMean};
}

bool multichain::snpSweepSupported(const vector<SweepChain> &chains) {
    if (chains.size() < 2) return false;
    for (const SweepChain &chain : chains) {
        if (!chain.model->lowRankModel) return false;
    }
    const Data &data = chains[0].model->data;
    return !data.coreContextLdmEigen && !data.integrateOutContext && !data.numFixedEffectSnps;
}

bool multichain::snpSweepSupported(const vector<ApproxBayesRC*> &chains) {
    for (const ApproxBayesRC *chain : chains) {
        if (chain->integrateOutContext) return false;
        if (chain->algorithm != ApproxBayesR::gibbs && chain->algorithm != ApproxBayesR::tgs
            && chain->algorithm != ApproxBayesR::tgs_thin) return false;
    }
    return snpSweepSupported(sweepChains(chains));
}

void multichain::sampleSnpEffects(const vector<SweepChain> &chains, const unsigned numThreads, SweepCache &cache) {
    const unsigned numChains = chains.size(), chunkSize = sweepChunkSize;
    const Data &data = chains[0].model->data;
    const vector<LDBlockInfo*> &blocks = data.keptLdBlockInfoVec;
    const unsigned nBlocks = blocks.size();
    const unsigned size = chains[0].effects->size;

    // Q's columns for a chunk: the dense block's, or with quantized LD (q16 / q8) the chunk's columns decoded into
    // buf (scale times quantized value, as QuantizedEigenBlock::decodeQColumn).
    const bool quantized = chains[0].model->quantizedLowRankModel;
    auto qRows = [&](const unsigned blk) { return quantized ? data.quantizedEigenBlocks[blk].k : (int)data.Qblocks[blk].rows(); };
    auto qCols = [&](const unsigned blk) { return quantized ? data.quantizedEigenBlocks[blk].m : (int)data.Qblocks[blk].cols(); };
    auto chunkColumns = [&](const unsigned blk, const unsigned first, const unsigned count, MatrixXf &buf) {
        const int k = qRows(blk);
        if (!quantized) return Map<const MatrixXf>(data.Qblocks[blk].data() + (size_t)first * k, k, count);
        const QuantizedEigenBlock &Qq = data.quantizedEigenBlocks[blk];
        const float *scale = Qq.qScale.data();
        buf.resize(k, count);
        for (unsigned j = 0; j < count; ++j) {
            float *out = buf.col(j).data();
            const size_t offset = (size_t)(first + j) * k;
            if (Qq.bits == 16) for (int r = 0; r < k; ++r) out[r] = scale[r] * (float)Qq.q16[offset + r];
            else for (int r = 0; r < k; ++r) out[r] = scale[r] * (float)Qq.q8[offset + r];
        }
        return Map<const MatrixXf>(buf.data(), k, count);
    };

    // Per block, once per run: the block's SNPs (positions) sorted by Q column, where each chunk of columns starts in
    // that list, and the SNPs without a column (eigenQColIndex < 0).
    typedef SweepCache::BlockLayout BlockLayout;
    vector<BlockLayout> &layout = cache.layout;
    if (layout.size() != nBlocks) {
        layout.assign(nBlocks, BlockLayout());
        #pragma omp parallel for schedule(dynamic) num_threads(numThreads)
        for (unsigned blk = 0; blk < nBlocks; ++blk) {
            const LDBlockInfo *blockInfo = blocks[blk];
            const unsigned blockStart = blockInfo->startSnpIdx;
            const unsigned blockSize = blockInfo->endSnpIdx - blockStart + 1;
            const int numCols = qCols(blk);
            BlockLayout &L = layout[blk];
            vector<std::pair<int, unsigned> > byCol;
            for (unsigned p = 0; p < blockSize; ++p) {
                const int qcol = eigenQColIndex(blockInfo, blockStart + p);
                if (qcol < 0) L.noColumn.push_back(p);
                else byCol.push_back(std::make_pair(qcol, p));
            }
            std::stable_sort(byCol.begin(), byCol.end(),
                             [](const std::pair<int, unsigned> &a, const std::pair<int, unsigned> &b) { return a.first < b.first; });
            for (const auto &x : byCol) {
                L.col.push_back(x.first);
                L.pos.push_back(x.second);
            }
            const unsigned numChunks = (numCols + chunkSize - 1) / chunkSize;
            L.chunkStart.resize(numChunks + 1);
            unsigned at = 0;
            for (unsigned c = 0; c <= numChunks; ++c) {
                while (at < L.col.size() && L.col[at] < (int)(c * chunkSize)) ++at;
                L.chunkStart[c] = at;
            }
            L.chunkStart[numChunks] = L.col.size();
        }
    }
    // Per chunk, once per run: G = Q_chunk' Q_chunk.
    vector<vector<float> > &chunkGram = cache.chunkGram;  // per block: chunk c's G at c * chunkSize^2, leading dimension chunkSize
    if (chunkGram.size() != nBlocks) {
        const double gramClock = omp_get_wtime();
        chunkGram.assign(nBlocks, vector<float>());
        #pragma omp parallel for schedule(dynamic) num_threads(numThreads)
        for (unsigned blk = 0; blk < nBlocks; ++blk) {
            const int numCols = qCols(blk);
            const unsigned numChunks = (numCols + chunkSize - 1) / chunkSize;
            MatrixXf buf;
            vector<float> &G = chunkGram[blk];
            G.assign((size_t)numChunks * chunkSize * chunkSize, 0.0f);
            for (unsigned c = 0; c < numChunks; ++c) {
                const unsigned first = c * chunkSize, count = std::min((int)chunkSize, numCols - (int)first);
                float *Gc = G.data() + (size_t)c * chunkSize * chunkSize;
                Map<MatrixXf, 0, OuterStride<> > Gmap(Gc, count, count, OuterStride<>(chunkSize));
                const Map<const MatrixXf> Qc = chunkColumns(blk, first, count, buf);
                Gmap.noalias() = Qc.transpose() * Qc;
            }
        }
        cout << "Batched sweep: chunk Gram matrices computed in " << omp_get_wtime() - gramClock << " s" << endl;
    }

    // Per chain: its draws and per-block tallies, and its sweep terms (prior variances, prior log pi, deltaPi); each chain
    // has its own number of components (the component selection's models).
    typedef ApproxBayesR::SnpEffects Effects;
    struct ChainState {
        std::unique_ptr<float[]> urnd, nrnd;  // drawn per block in the block loop
        vector<Effects::BlockTally> tallies;
    };
    vector<ChainState> state(numChains);
    vector<Effects::SweepTerms> terms;
    terms.reserve(numChains);
    vector<Stat::random_engine::result_type> chunkSeed(nBlocks);  // shared by the chains: keyed by iteration and block
    for (unsigned blk = 0; blk < nBlocks; ++blk) chunkSeed[blk] = Stat::parallelTaskSeed(blk, Stat::batchSweepChunkRngSalt);
    for (unsigned h = 0; h < numChains; ++h) {
        ApproxBayesR &chain = *chains[h].model;
        ChainState &s = state[h];
        s.urnd.reset(new float[size]);
        s.nrnd.reset(new float[size]);
        chain.whatBlocks.resize(nBlocks);
        for (unsigned blk = 0; blk < nBlocks; ++blk) chain.whatBlocks[blk].resize(chain.wcorrBlocks[blk].size());
        s.tallies = Effects::blockTallies(nBlocks, chain.gamma.values.size());
        if (chains[h].snpLogPi)
            terms.emplace_back(chain.gamma.values, chain.sigmaSq.value, chain.vargRandom.value, chain.hsqPercModel,
                               chain.deltaPi, *chains[h].snpLogPi, *chains[h].fcMean);
        else
            terms.emplace_back(chain.gamma.values, chain.sigmaSq.value, chain.vargRandom.value, chain.hsqPercModel,
                               chain.deltaPi, chain.Pis.values);
    }
    // The per-SNP outputs start at zero, as in sampleFromFC_eigen (the chains in parallel).
    #pragma omp parallel for num_threads(numThreads)
    for (unsigned h = 0; h < numChains; ++h) chains[h].effects->startSweep(chains[h].model->deltaPi);

    #pragma omp parallel for schedule(dynamic) num_threads(numThreads)
    for (unsigned blk = 0; blk < nBlocks; ++blk) {
        const LDBlockInfo *blockInfo = blocks[blk];
        const VectorXf &qDiag = blockInfo->qDiag;
        const unsigned blockStart = blockInfo->startSnpIdx;
        const int k = qRows(blk), numCols = qCols(blk);
        MatrixXf chunkBuf;  // quantized LD: the current chunk's columns, decoded
        const BlockLayout &L = layout[blk];

        MatrixXf W(k, numChains);  // the chains' residuals for this block, side by side
        for (unsigned h = 0; h < numChains; ++h) {
            W.col(h) = chains[h].model->wcorrBlocks[blk];
            chains[h].model->whatBlocks[blk].setZero();
        }
        vector<unsigned> chunks(L.chunkStart.size() - 1);
        for (unsigned c = 0; c < chunks.size(); ++c) chunks[c] = c;
        std::mt19937 chunkRng(chunkSeed[blk]);
        std::shuffle(chunks.begin(), chunks.end(), chunkRng);
        vector<std::mt19937> withinRng;  // each chain's orders within the chunks

        // Each chain's draws for this block, and the seed of its orders within the chunks, from streams keyed by
        // chain, iteration and block.
        for (unsigned h = 0; h < numChains; ++h) {
            Stat::McmcChainRngScope chainRng(h);
            withinRng.emplace_back(Stat::parallelTaskSeed(blk, Stat::batchSweepOrderRngSalt));
            Stat::ParallelTaskRngScope drawRng(blk, Stat::batchSweepDrawRngSalt);
            ChainState &s = state[h];
            for (unsigned i = blockStart; i <= blockInfo->endSnpIdx; ++i) {
                s.urnd[i] = Stat::ranf();
                s.nrnd[i] = Stat::snorm();
            }
        }
        // SNPs without a Q column get only the skip and bad-SNP updates, as in sampleFromFC_eigen.
        for (unsigned h = 0; h < numChains; ++h)
            for (unsigned p : L.noColumn)
                chains[h].effects->clearUnsampledSnp(blockStart + p, blockInfo->snpInfoVec[p], chains[h].model->deltaPi);

        MatrixXf R(chunkSize, numChains);  // q'w for the chunk's columns (rows) and the chains (columns)
        vector<unsigned> order;            // positions in L.pos, in a chain's order
        vector<Effects::SnpScratch> scratch;  // the per-SNP arrays, one set per chain
        scratch.reserve(numChains);
        for (unsigned h = 0; h < numChains; ++h) scratch.emplace_back(terms[h].invWtdSigmaSq.size());
        for (unsigned ci = 0; ci < chunks.size(); ++ci) {
            const unsigned c = chunks[ci];
            const unsigned begin = L.chunkStart[c], end = L.chunkStart[c + 1];
            if (begin == end) continue;
            const unsigned first = c * chunkSize, count = std::min((int)chunkSize, numCols - (int)first);
            const Map<const MatrixXf> Qc = chunkColumns(blk, first, count, chunkBuf);
            chunkProduct(Qc, W, R);
            for (unsigned h = 0; h < numChains; ++h) {
                ApproxBayesR &chain = *chains[h].model;
                Effects &effects = *chains[h].effects;
                ChainState &s = state[h];
                Ref<VectorXf> what = chain.whatBlocks[blk];
                const float invVareDn = data.nGWASblock[blk] / chain.vareBlk.values[blk];
                order.resize(end - begin);
                for (unsigned t = 0; t < order.size(); ++t) order[t] = begin + t;
                std::shuffle(order.begin(), order.end(), withinRng[h]);
                for (unsigned at : order) {
                    const unsigned blockIndex = L.pos[at], i = blockStart + blockIndex;
                    const int qcol = L.col[at];
                    const unsigned local = qcol - first;
                    if (effects.clearUnsampledSnp(i, blockInfo->snpInfoVec[blockIndex], chain.deltaPi)) continue;
                    const float oldSample = effects.values[i];
                    Ref<const VectorXf> Qi = Qc.col(local);
                    const unsigned delta = effects.sampleSnp(i, R(local, h), qDiag[qcol], invVareDn, terms[h],
                                                             s.urnd[i], s.nrnd[i], scratch[h], s.tallies[blk]);
                    if (delta) what += Qi * effects.values[i];
                    const float change = oldSample - effects.values[i];
                    if (change) {  // update this chain's residual and the chunk's remaining terms
                        W.col(h) += Qi * change;
                        const float *g = chunkGram[blk].data() + ((size_t)c * chunkSize + local) * chunkSize;
                        float *r = R.col(h).data();
                        for (unsigned row = 0; row < count; ++row) r[row] += g[row] * change;
                    }
                }
            }
        }
        for (unsigned h = 0; h < numChains; ++h) chains[h].model->wcorrBlocks[blk] = W.col(h);
    }

    // Tallies, as at the end of sampleFromFC_eigen (no fixed-effect SNPs here).
    for (unsigned h = 0; h < numChains; ++h) {
        if (chains[h].numSnpMix) chains[h].effects->finishSweep(state[h].tallies, *chains[h].numSnpMix, *chains[h].snpset);
        else chains[h].effects->finishSweep(state[h].tallies, *chains[h].snpStore, *chains[h].snpset);
    }
}
