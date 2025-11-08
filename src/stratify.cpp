//
//  stratify.cpp
//  gctb
//
//  Created by Jian Zeng on 25/06/2018.
//  Copyright © 2018 Jian Zeng. All rights reserved.
//

#include "stratify.hpp"

void StratApproxBayesS::VarEffectStratified::sampleFromFC(const VectorXf &snpEffSumSq, const VectorXf &numSnpEff) {
    for (unsigned i=0; i<size; ++i) {
        float dfTilde = df + numSnpEff[i];
        float scaleTilde = snpEffSumSq[i] + df*scales[i];
        values[i] = InvChiSq::sample(dfTilde, scaleTilde);
    }
}

void StratApproxBayesS::VarEffectStratified::sampleFromPrior() {
    for (unsigned i=0; i<size; ++i) {
        values[i] = InvChiSq::sample(df, scales[i]);
    }
}

void StratApproxBayesS::VarEffectEnrichment::compute(const VectorXf &sigmaSqStrat, const float sigmaSq) {
    for (unsigned i=0; i<size; ++i) {
        values[i] = sigmaSqStrat[i]/sigmaSq;
    }
}

void StratApproxBayesS::PiStratified::sampleFromFC(const vector<unsigned> &numSnps, const VectorXf &numSnpEff) {
    for (unsigned i=0; i<size; ++i) {
        float alphaTilde = numSnpEff[i] + alpha;
        float betaTilde  = numSnps[i] - numSnpEff[i] + beta;
        values[i] = Beta::sample(alphaTilde, betaTilde);
    }
}

void StratApproxBayesS::PiStratified::sampleFromPrior() {
    for (unsigned i=0; i<size; ++i) {
        values[i] = Beta::sample(alpha, beta);
    }
}

void StratApproxBayesS::PiEnrichment::compute(const VectorXf &nnzStrat, const float nnzTotal) {
    for (unsigned i=0; i<size; ++i) {
        float obs = nnzStrat[i]/nnzTotal;
        values[i] = obs/expectation[i];
    }
}

void StratApproxBayesS::HeritabilityStratified::compute(const VectorXf &sigmaSq, const VectorXf &sum2pqSplusOne, const float genVar, const float resVar) {
    values = sigmaSq.cwiseProduct(sum2pqSplusOne)/(genVar + resVar);
}

void StratApproxBayesS::HeritabilityStratified::compute(const vector<VectorXf> &snpEffectPerAnno, const vector<SpMat> &annowiseZPZsp, const vector<VectorXf> &annowiseZPZdiag, const float genVar, const float resVar) {
    float varp = genVar + resVar;

    long chunkSize = size/omp_get_max_threads();
#pragma omp parallel for schedule(dynamic, chunkSize)
    for (unsigned i=0; i<size; ++i) {
        values[i] = 2.0*snpEffectPerAnno[i].transpose()*annowiseZPZsp[i]*snpEffectPerAnno[i] + snpEffectPerAnno[i].cwiseProduct(snpEffectPerAnno[i]).dot(annowiseZPZdiag[i]);
        values[i] /= float(sampleSize);
        values[i] /= varp;
    }
}

void StratApproxBayesS::PerSnpHeritabilityEnrichment::compute(const VectorXf &hsqStrat, const VectorXf &expFrac, const float hsqTotal) {
    for (unsigned i=0; i<size; ++i) {
        values[i] = hsqStrat[i]/hsqTotal/expFrac[i];
    }
}

void StratApproxBayesS::PerNzHeritabilityEnrichment::compute(const VectorXf &hsqStrat, const VectorXf &nnzStrat, const float hsqTotal, const float nnzTotal) {
    float expectation = hsqTotal/float(nnzTotal);
    for (unsigned i=0; i<size; ++i) {
        float obs = nnzStrat[i] ? hsqStrat[i]/float(nnzStrat[i]) : 0;
        values[i] = obs/expectation;
    }
}

void StratApproxBayesS::SpStratified::sampleFromFC(const vector<VectorXf> &snpEffects, const VectorXf &numNonZeros, const VectorXf &sigmaSq, const VectorXf &hsq, const float genVar, const float resVar, const vector<AnnoInfo*> &annoInfoVec, VectorXf &scales, VectorXf &sum2pqSplusOneVec) {
    
    VectorXf varg = hsq.array()*(genVar + resVar);
    
    long chunkSize = size/omp_get_max_threads();
#pragma omp parallel for schedule(dynamic, chunkSize)
    for (unsigned i=0; i<size; ++i) {
        unsigned nnzi = numNonZeros[i];
        VectorXf snpEffectsAnnoi(nnzi);
        VectorXf snp2pqAnnoi(nnzi);
        VectorXf snp2pqLogAnnoi(nnzi);
        AnnoInfo *anno = annoInfoVec[i];
        unsigned idx = 0;
        for (unsigned j=0; j<anno->size; ++j) {
            if (snpEffects[i](j)) {
                snpEffectsAnnoi[idx] = snpEffects[i](j);
                snp2pqAnnoi[idx] = anno->snp2pq[j];
                snp2pqLogAnnoi[idx] = snp2pqLog[i](j);
                ++idx;
            }
        }

        if (nnzi < 3) {
            values[i] = Stat::snorm()*sqrtf(var);
            if (nnzi) {
                sum2pqSplusOneVec[i] = snp2pqAnnoi.sum();
                scales[i] = 0.5*varg[i]/sum2pqSplusOneVec[i];
            }
        } else {
            hmcSampler(i, snpEffectsAnnoi, snp2pqAnnoi, snp2pqLogAnnoi, sigmaSq[i], varg[i], scales[i], sum2pqSplusOneVec[i], values[i]);
        }
    }
}

void StratApproxBayesS::SpStratified::hmcSampler(const unsigned annoIdx, const VectorXf &snpEffects, const VectorXf &snp2pq, const VectorXf &snp2pqLog, const float sigmaSq, const float varg, float &scale, float &sum2pqSplusOne, float &value) {
    float snp2pqLogSum = snp2pqLog.sum();
    
    float curr = value;
    float curr_p = Stat::snorm();
    
    float cand = curr;
    // Make a half step for momentum at the beginning
    
    gradientU(curr, snpEffects, snp2pqLogSum, snp2pq, snp2pqLog, sigmaSq);
    
    float cand_p = curr_p - 0.5*stepSize[annoIdx] * gradientU(curr, snpEffects, snp2pqLogSum, snp2pq, snp2pqLog, sigmaSq);
    
    for (unsigned i=0; i<numSteps; ++i) {
        // Make a full step for the position
        cand += stepSize[annoIdx] * cand_p;
        if (i < numSteps-1) {
            // Make a full step for the momentum, except at end of trajectory
            cand_p -= stepSize[annoIdx] * gradientU(cand, snpEffects, snp2pqLogSum, snp2pq, snp2pqLog, sigmaSq);
        } else {
            // Make a half step for momentum at the end
            cand_p -= 0.5*stepSize[annoIdx] * gradientU(cand, snpEffects, snp2pqLogSum, snp2pq, snp2pqLog, sigmaSq);
        }
    }
    
    // Evaluate potential (negative log posterior) and kinetic energies at start and end of trajectory
    float curr_H = computeU(curr, snpEffects, snp2pqLogSum, snp2pq, sigmaSq) + 0.5*curr_p*curr_p;
    float cand_H = computeU(cand, snpEffects, snp2pqLogSum, snp2pq, sigmaSq) + 0.5*cand_p*cand_p;
    
    if (Stat::ranf() < exp(curr_H-cand_H)) {  // accept
        value = cand;
        sum2pqSplusOne = snp2pq.array().pow(value+1.0).sum();
        scale = 0.5*varg/sum2pqSplusOne;
        if (scale != scale) {
            cout << snp2pq << endl;
            cout << "pow sum: " << snp2pq.array().pow(value+1.0).sum() << endl;
            throw("Error: scale is nan!");
        }
        ars[annoIdx]->count(1, 0.5, 0.9);
    } else {
        ars[annoIdx]->count(0, 0.5, 0.9);
    }
    
    if (!(ars[annoIdx]->cnt % 10)) {
        if      (ars[annoIdx]->value < 0.6) stepSize[annoIdx] *= 0.8;
        else if (ars[annoIdx]->value > 0.8) stepSize[annoIdx] *= 1.2;
    }
    
    if (ars[annoIdx]->consecRej > 20) stepSize[annoIdx] *= 0.8;
}

float StratApproxBayesS::SpStratified::gradientU(const float S, const VectorXf &snpEffects, const float snp2pqLogSum, const VectorXf &snp2pq, const VectorXf &snp2pqLog, const float sigmaSq){
    // compute the first derivative of the negative log posterior
    return 0.5*snp2pqLogSum - 0.5/sigmaSq*(snpEffects.array().square()*snp2pqLog.array()/snp2pq.array().pow(S)).sum() + S/var;
}

float StratApproxBayesS::SpStratified::computeU(const float S, const VectorXf &snpEffects, const float snp2pqLogSum, const VectorXf &snp2pq, const float sigmaSq){
    // compute negative log posterior and scale
    return 0.5*S*snp2pqLogSum + 0.5/sigmaSq*(snpEffects.array().square()/snp2pq.array().pow(S)).sum() + 0.5*S*S/var;
}

void StratApproxBayesS::SpStratified::sampleFromPrior() {
    for (unsigned i=0; i<size; ++i) {
        values[i] = Stat::snorm()*sqrtf(var);
    }
}

void StratApproxBayesS::SpEnrichment::compute(const VectorXf &Sstrat, const float S) {
    for (unsigned i=0; i<size; ++i) {
        values[i] = Sstrat[i]/S;
    }
}

void StratApproxBayesS::SnpAnnoMembership::getValues(const VectorXf &snpAnnoVec) {
    long numSnps = snpAnnoVec.size();
    long pos = 0;
    values.setZero(size);
    for (unsigned i=0; i<numSnps; ++i) {
        values[snpAnnoVec[i] + pos] = 1;
        pos += numAnnoPerSnpVec[i];
        //if (i==0) cout << "i " << i << " snpAnnoVec[i] " << snpAnnoVec[i] << " numAnnoPerSnpVec[i] " << numAnnoPerSnpVec[i] << " pos " << pos << " size " << size << endl;
    }
}

void StratApproxBayesS::PerSnpPi::sampleFromFC(const MatrixXf &annoMat, const MatrixXf &APA, const VectorXf &snpEffects) {
    // MH random walk
    // sample beta candidates
    VectorXf betaCand(beta.size());
    for (unsigned i=0; i<beta.size(); ++i) {
        betaCand[i] = sample(beta[i], varProp);
    }
    
    // compute pi given beta candidates
    VectorXf curr = values;
    VectorXf cand = 1.0f/(1.0f + (-annoMat*betaCand).array().exp());  // expit function
    VectorXf logCurr = values.array().log();
    VectorXf logCand = cand.array().log();
    VectorXf logCurrComp = (1.0-values.array()).log();
    VectorXf logCandComp = (1.0-cand.array()).log();
    
    // compute log likelihood
    double logLikeCurr = 0.0;
    double logLikeCand = 0.0;
//    int nnz = 0, nz=0;
    for (unsigned i=0; i<snpEffects.size(); ++i) {
        if (snpEffects[i]) {
            logLikeCurr += logCurr[i];
            logLikeCand += logCand[i];
//            ++nnz;
        } else {
            logLikeCurr += logCurrComp[i];
            logLikeCand += logCandComp[i];
//            ++nz;
        }
    }
    
    // compute log posterior
    double logPostCurr = logLikeCurr - 0.5f*beta.squaredNorm()/sigmaSq;
    double logPostCand = logLikeCand - 0.5f*betaCand.squaredNorm()/sigmaSq;

    if (Stat::ranf() < exp(logPostCand-logPostCurr)) {  // accept
        values = cand;
        beta = betaCand;
        ar.count(1, 0.1, 0.5);
    } else {
        ar.count(0, 0.1, 0.5);
    }

//    cout << "betaCand " << betaCand.transpose() << endl;
//    cout << "betaCurr " << beta.transpose() << endl;
//    cout << "cand " << cand[0] << " " << cand.mean() << " logCand " << logCand[0] << " " << logCand.mean() << " logCandComp " << logCandComp[0] << " " << logCandComp.mean() << endl;
//    cout << "curr " << curr[0] << " " << curr.mean() << " logCurr " << logCurr[0] << " " << logCurr.mean() << " logCurrComp " << logCurrComp[0] << " " << logCurrComp.mean() << endl;
////    cout << "nnz " << nnz << " nz " << nz << endl;
//    cout << "logLikeCand " << logLikeCand << " logLikeCurr " << logLikeCurr << endl;
//    cout << "logPostCand " << logPostCand << " logPostCurr " << logPostCurr << " exp(logPostCand-logPostCurr) " << exp(logPostCand-logPostCurr) << " ar " << ar.value << endl;
        
    if (!(ar.cnt % 10)) {
        if      (ar.value < 0.2) varProp *= 0.8;
        else if (ar.value > 0.5) varProp *= 1.2;
    }

//    cout << " ar " << ar.value << " " << varProp << endl;
    
    
//
////    cout << "PerSnpPi begin" << endl;
//    VectorXf y = (values.array()/(1.0f-values.array())).log();
//    beta = APA.inverse()*(annoMat.transpose()*y);
//    VectorXf yhat = annoMat*beta;
//    values = yhat.array().exp()/(1.0f+yhat.array().exp());
////    cout << "y " << y.head(10).transpose() << endl;
////    cout << "mean y " << y.mean() << endl;
////    cout << "APA " << APA << endl;
////    cout << "annoMat " << annoMat.block(0, 0, 10, 10) << endl;
////    cout << "annoMat.transpose()*y " << (annoMat.transpose()*y).head(10).transpose() << endl;
////    cout << " beta " << beta.transpose() << endl;
////    cout << "PerSnpPi end" << endl;
}

void StratApproxBayesS::PerSnpPi::computeFromAnnoPi(const MatrixXf &annoMat, const VectorXf &piStrat) {
    for (unsigned i=0; i<size; ++i) {
        values[i] = annoMat.row(i).dot(piStrat)/annoMat.row(i).sum();
    }
}

void StratApproxBayesS::SnpEffects::sampleFromFCMixture(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag, const vector<ChromInfo *> &chromInfoVec, const vector<SnpInfo *> &incdSnpInfoVec, const VectorXf &snp2pq, const VectorXf &LDsamplVar, const unsigned numAnnos, const VectorXf &sigmaSq, const VectorXf &pi, const VectorXf &S, const float Sgw, const float varg, const float vare, const float ps, const float overdispersion) {

    long numChr = chromInfoVec.size();

    wtdSumSq = 0;
    numNonZeros = 0;
    sum2pqBetaSq = 0;
    wtdSumSqPerAnno.setZero(numAnnos);
    numNonZeroPerAnno.setZero(numAnnos);
    sum2pqBetaSqAnno.setZero(numAnnos);

    VectorXf logPi = pi.array().log();
    VectorXf logPiComp = (1.0-pi.array()).log();
    VectorXf invSigmaSq = sigmaSq.cwiseInverse();

    for (unsigned i=0; i<numAnnos; ++i) {
        valuesPerAnno[i].setZero(valuesPerAnno[i].size());
    }

    for (unsigned chr=0; chr<numChr; ++chr) {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned annoIdx, snpIdx;
        unsigned i, j;

        float oldSample;
        float rhs;
        float probDelta1;
        float varei;
        float sampleDiff;

        ArrayXf invLhs, uhat;
        ArrayXf snp2pqPowS;
        ArrayXf logDelta1Anno;
        ArrayXf logDelta0Anno;
        ArrayXf probDeltaAnno;

        SnpInfo *snp;

        for (i = chrStart; i <= chrEnd; ++i) {
            snp = incdSnpInfoVec[i];
            invLhs.resize(snp->numAnnos);
            uhat.resize(snp->numAnnos);
            snp2pqPowS.resize(snp->numAnnos);
            logDelta1Anno.resize(snp->numAnnos);
            logDelta0Anno.resize(snp->numAnnos);
            probDeltaAnno.resize(snp->numAnnos);

            oldSample = values[i];
            varei = LDsamplVar[i]*varg + vare + ps + overdispersion;

            rhs  = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;

            for (j=0; j<snp->numAnnos; ++j) {
                annoIdx = snp->annoVec[j]->idx;

                snp2pqPowS[j] = powf(snp2pq[i], S[annoIdx]);
                invLhs[j] = 1.0f/(ZPZdiag[i]/varei + invSigmaSq[annoIdx]/snp2pqPowS[j]);
                uhat[j] = invLhs[j]*rhs;

                logDelta1Anno[j]  = 0.5*(logf(invLhs[j]) - logf(snp2pqPowS[j]*sigmaSq[annoIdx]) + uhat[j]*rhs) + logPi[annoIdx];
                logDelta0Anno[j] = logPiComp[annoIdx];
            }

            for (j=0; j<snp->numAnnos; ++j) {
                probDeltaAnno[j] = (1.0+expf(logDelta0Anno[j]-logDelta1Anno[j]))/((logDelta1Anno-logDelta1Anno[j]).exp().sum() + (logDelta0Anno-logDelta1Anno[j]).exp().sum());
            }

            snpAnnoVec[i] = j = bernoulli.sample(probDeltaAnno);

            annoIdx = snp->annoVec[j]->idx;

            probDelta1 = 1.0f/(1.0f + expf(logDelta0Anno[j]-logDelta1Anno[j]));

            if (bernoulli.sample(probDelta1)) {
                values[i] = normal.sample(uhat[j], invLhs[j]);
                sampleDiff = oldSample - values[i];
                for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                    rcorr[it.index()] += it.value() * sampleDiff;
                }
                wtdSumSq += values[i]*values[i]/powf(snp2pq[i], Sgw);
                sum2pqBetaSq += snp2pq[i]*values[i]*values[i];
                ++numNonZeros;

                for (j=0; j<snp->numAnnos; ++j) {
                    annoIdx = snp->annoVec[j]->idx;
                    snpIdx = snp->annoIdx[j];
                    valuesPerAnno[annoIdx][snpIdx] = values[i];
                    wtdSumSqPerAnno[annoIdx] += values[i]*values[i]/snp2pqPowS[j];
                    sum2pqBetaSqAnno[annoIdx] += snp2pq[i]*values[i]*values[i];
                    ++numNonZeroPerAnno[annoIdx];
                }
            } else {
                if (oldSample) {
                    for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                        rcorr[it.index()] += it.value() * oldSample;
                    }
                }
                values[i] = 0.0;
            }
        }
    }
}

void StratApproxBayesS::SnpEffects::sampleFromFCMixture(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo *> &chromInfoVec, const vector<SnpInfo *> &incdSnpInfoVec, const VectorXf &snp2pq, const VectorXf &LDsamplVar, const unsigned numAnnos, const VectorXf &sigmaSq, const VectorXf &pi, const VectorXf &S, const float Sgw, const float varg, const float vare, const float ps, const float overdispersion) {
    
    long numChr = chromInfoVec.size();
    
    wtdSumSq = 0;
    numNonZeros = 0;
    sum2pqBetaSq = 0;
    wtdSumSqPerAnno.setZero(numAnnos);
    numNonZeroPerAnno.setZero(numAnnos);
    sum2pqBetaSqAnno.setZero(numAnnos);

    VectorXf logPi = pi.array().log();
    VectorXf logPiComp = (1.0-pi.array()).log();
    VectorXf invSigmaSq = sigmaSq.cwiseInverse();
    
    for (unsigned i=0; i<numAnnos; ++i) {
        valuesPerAnno[i].setZero(valuesPerAnno[i].size());
    }
    
    for (unsigned chr=0; chr<numChr; ++chr) {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned annoIdx, snpIdx;
        unsigned i, j;
        
        float oldSample;
        float rhs;
        float probDelta1;
        float varei;
        float sampleDiff;
        
        ArrayXf invLhs, uhat;
        ArrayXf snp2pqPowS;
        ArrayXf logDelta1Anno;
        ArrayXf logDelta0Anno;
        ArrayXf probDeltaAnno;
        
        SnpInfo *snp;
        
        for (i = chrStart; i <= chrEnd; ++i) {
            snp = incdSnpInfoVec[i];
            invLhs.resize(snp->numAnnos);
            uhat.resize(snp->numAnnos);
            snp2pqPowS.resize(snp->numAnnos);
            logDelta1Anno.resize(snp->numAnnos);
            logDelta0Anno.resize(snp->numAnnos);
            probDeltaAnno.resize(snp->numAnnos);
            
            oldSample = values[i];
            varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
            
            rhs  = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            
            for (j=0; j<snp->numAnnos; ++j) {
                annoIdx = snp->annoVec[j]->idx;
                
                snp2pqPowS[j] = powf(snp2pq[i], S[annoIdx]);
                invLhs[j] = 1.0f/(ZPZdiag[i]/varei + invSigmaSq[annoIdx]/snp2pqPowS[j]);
                uhat[j] = invLhs[j]*rhs;
                
                logDelta1Anno[j]  = 0.5*(logf(invLhs[j]) - logf(snp2pqPowS[j]*sigmaSq[annoIdx]) + uhat[j]*rhs) + logPi[annoIdx];
                logDelta0Anno[j] = logPiComp[annoIdx];
            }
            
            for (j=0; j<snp->numAnnos; ++j) {
                probDeltaAnno[j] = (1.0+expf(logDelta0Anno[j]-logDelta1Anno[j]))/((logDelta1Anno-logDelta1Anno[j]).exp().sum() + (logDelta0Anno-logDelta1Anno[j]).exp().sum());
            }
            
            snpAnnoVec[i] = j = bernoulli.sample(probDeltaAnno);
            
            annoIdx = snp->annoVec[j]->idx;
            
            probDelta1 = 1.0f/(1.0f + expf(logDelta0Anno[j]-logDelta1Anno[j]));
            
            if (bernoulli.sample(probDelta1)) {
                values[i] = normal.sample(uhat[j], invLhs[j]);
                sampleDiff = oldSample - values[i];
                rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*sampleDiff;
                wtdSumSq += values[i]*values[i]/powf(snp2pq[i], Sgw);
                sum2pqBetaSq += snp2pq[i]*values[i]*values[i];
                ++numNonZeros;
                
                for (j=0; j<snp->numAnnos; ++j) {
                    annoIdx = snp->annoVec[j]->idx;
                    snpIdx = snp->annoIdx[j];
                    valuesPerAnno[annoIdx][snpIdx] = values[i];
                    wtdSumSqPerAnno[annoIdx] += values[i]*values[i]/snp2pqPowS[j];
                    sum2pqBetaSqAnno[annoIdx] += snp2pq[i]*values[i]*values[i];
                    ++numNonZeroPerAnno[annoIdx];
                }
            } else {
                if (oldSample) {
                    rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*oldSample;
                }
                values[i] = 0.0;
            }
        }
    }
}

void StratApproxBayesS::SnpEffects::sampleFromFCLinear(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag, const vector<ChromInfo *> &chromInfoVec, const vector<SnpInfo *> &incdSnpInfoVec, const VectorXf &snp2pq, const VectorXf &LDsamplVar, const unsigned numAnnos, const VectorXf &sigmaSq, const VectorXf &pi, const VectorXf &S, const float Sgw, const float varg, const float vare, const float ps, const float overdispersion) {
    
    long numChr = chromInfoVec.size();
    
    wtdSumSq = 0;
    numNonZeros = 0;
    sum2pqBetaSq = 0;
    wtdSumSqPerAnno.setZero(numAnnos);
    numNonZeroPerAnno.setZero(numAnnos);
    sum2pqBetaSqAnno.setZero(numAnnos);
    
    VectorXf logPi = pi.array().log();
    VectorXf logPiComp = (1.0-pi.array()).log();
    
    for (unsigned chr=0; chr<numChr; ++chr) {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned annoIdx, snpIdx;
        
        float oldSample;
        float rhs;
        float logDelta1;
        float logDelta0;
        float probDelta1;
        float varei;
        float sigmaSqi;
        float sampleDiff;
        float invLhs, uhat;
        float rcorrAnno;
        
        ArrayXf snp2pqPowSAnno;
        ArrayXf sigmaSqAnno;
        ArrayXf betaAnno;
        
        SnpInfo *snp;
        
        for (unsigned i = chrStart; i <= chrEnd; ++i) {
            snp = incdSnpInfoVec[i];
            oldSample = values[i];
            
            varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
                        
            snp2pqPowSAnno.resize(snp->numAnnos);
            sigmaSqAnno.resize(snp->numAnnos);
            betaAnno.resize(snp->numAnnos);

            for (unsigned j=0; j<snp->numAnnos; ++j) {
                annoIdx = snp->annoVec[j]->idx;
                snpIdx = snp->annoIdx[j];
                snp2pqPowSAnno[j] = powf(snp2pq[i], S[annoIdx]);
                sigmaSqAnno[j] = sigmaSq[annoIdx];
                betaAnno[j] = valuesPerAnno[annoIdx][snpIdx];
            }
            
            sigmaSqi = snp2pqPowSAnno.cwiseProduct(sigmaSqAnno).sum();
            
            rhs  = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            invLhs = 1.0f/(ZPZdiag[i]/varei + 1.0f/sigmaSqi);
            uhat = invLhs*rhs;
            
            logDelta1 = 0.5*(logf(invLhs) - logf(sigmaSqi) + uhat*rhs) + logPi[i];
            logDelta0 = logPiComp[i];
            
            probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
            
            if (bernoulli.sample(probDelta1)) {
                // sample effect for each annotation
                rcorrAnno = rcorr[i];
                for (unsigned j=0; j<snp->numAnnos; ++j) {
                    annoIdx = snp->annoVec[j]->idx;
                    snpIdx = snp->annoIdx[j];
                    
                    rcorrAnno += ZPZdiag[i]*betaAnno[j];
                    rhs = rcorrAnno/varei;
                    invLhs = 1.0f/(ZPZdiag[i]/varei + 1.0f/(snp2pqPowSAnno[j]*sigmaSqAnno[j]));
                    uhat = invLhs*rhs;
                    betaAnno[j] = normal.sample(uhat, invLhs);
                    rcorrAnno -= ZPZdiag[i]*betaAnno[j];
                    
                    valuesPerAnno[annoIdx][snpIdx] = betaAnno[j];
                    wtdSumSqPerAnno[annoIdx] += betaAnno[j]*betaAnno[j]/snp2pqPowSAnno[j];
                    ++numNonZeroPerAnno[annoIdx];
                }
                
                // total effect size
                values[i] = betaAnno.sum();
                sampleDiff = oldSample - values[i];
                for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                    rcorr[it.index()] += it.value() * sampleDiff;
                }
                wtdSumSq += values[i]*values[i]/powf(snp2pq[i], Sgw);
                sum2pqBetaSq += snp2pq[i]*values[i]*values[i];
                for (unsigned j=0; j<snp->numAnnos; ++j) {
                    annoIdx = snp->annoVec[j]->idx;
                    sum2pqBetaSqAnno[annoIdx] += snp2pq[i]*values[i]*values[i];
                }
                ++numNonZeros;
            }
            else {
                if (oldSample) {
                    for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                        rcorr[it.index()] += it.value() * oldSample;
                    }
                }
                values[i] = 0.0;
                for (unsigned j=0; j<snp->numAnnos; ++j) {
                    annoIdx = snp->annoVec[j]->idx;
                    snpIdx = snp->annoIdx[j];
                    valuesPerAnno[annoIdx][snpIdx] = 0.0;
                }
            }
        }
    }
}

void StratApproxBayesS::SnpEffects::sampleFromFCLinear(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag,
                                                       const VectorXi &windStart, const VectorXi &windSize,
                                                       const vector<ChromInfo*> &chromInfoVec, const vector<SnpInfo*> &incdSnpInfoVec,
                                                       const VectorXf &snp2pq, const VectorXf &LDsamplVar, const unsigned numAnnos,
                                                       const VectorXf &sigmaSq, const VectorXf &pi, const VectorXf &S, const float Sgw,
                                                       const float varg, const float vare, const float ps, const float overdispersion) {
    
    long numChr = chromInfoVec.size();
    
    wtdSumSq = 0;
    numNonZeros = 0;
    sum2pqBetaSq = 0;
    wtdSumSqPerAnno.setZero(numAnnos);
    numNonZeroPerAnno.setZero(numAnnos);
    sum2pqBetaSqAnno.setZero(numAnnos);

    VectorXf logPi = pi.array().log();
    VectorXf logPiComp = (1.0-pi.array()).log();
    
    for (unsigned chr=0; chr<numChr; ++chr) {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned annoIdx, snpIdx;
        
        float oldSample;
        float rhs;
        float logDelta1;
        float logDelta0;
        float probDelta1;
        float varei;
        float sigmaSqi;
        float sampleDiff;
        float invLhs, uhat;
        float rcorrAnno;
        
        ArrayXf snp2pqPowSAnno;
        ArrayXf sigmaSqAnno;
        ArrayXf betaAnno;
        
        SnpInfo *snp;
        
        for (unsigned i = chrStart; i <= chrEnd; ++i) {
            snp = incdSnpInfoVec[i];
            oldSample = values[i];
            
            varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
                        
            snp2pqPowSAnno.resize(snp->numAnnos);
            sigmaSqAnno.resize(snp->numAnnos);
            betaAnno.resize(snp->numAnnos);

            for (unsigned j=0; j<snp->numAnnos; ++j) {
                annoIdx = snp->annoVec[j]->idx;
                snpIdx = snp->annoIdx[j];
                snp2pqPowSAnno[j] = powf(snp2pq[i], S[annoIdx]);
                sigmaSqAnno[j] = sigmaSq[annoIdx];
                betaAnno[j] = valuesPerAnno[annoIdx][snpIdx];
            }
            
            sigmaSqi = snp2pqPowSAnno.cwiseProduct(sigmaSqAnno).sum();
            
            rhs  = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            invLhs = 1.0f/(ZPZdiag[i]/varei + 1.0f/sigmaSqi);
            uhat = invLhs*rhs;
            
            logDelta1 = 0.5*(logf(invLhs) - logf(sigmaSqi) + uhat*rhs) + logPi[i];
            logDelta0 = logPiComp[i];
            
            probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
            
            if (bernoulli.sample(probDelta1)) {
                // sample effect for each annotation
                rcorrAnno = rcorr[i];
                for (unsigned j=0; j<snp->numAnnos; ++j) {
                    annoIdx = snp->annoVec[j]->idx;
                    snpIdx = snp->annoIdx[j];
                    
                    rcorrAnno += ZPZdiag[i]*betaAnno[j];
                    rhs = rcorrAnno/varei;
                    invLhs = 1.0f/(ZPZdiag[i]/varei + 1.0f/(snp2pqPowSAnno[j]*sigmaSqAnno[j]));
                    uhat = invLhs*rhs;
                    betaAnno[j] = normal.sample(uhat, invLhs);
                    rcorrAnno -= ZPZdiag[i]*betaAnno[j];
                    
                    valuesPerAnno[annoIdx][snpIdx] = betaAnno[j];
                    wtdSumSqPerAnno[annoIdx] += betaAnno[j]*betaAnno[j]/snp2pqPowSAnno[j];
                    ++numNonZeroPerAnno[annoIdx];
                }
                
                // total effect size
                values[i] = betaAnno.sum();
                sampleDiff = oldSample - values[i];
                rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*sampleDiff;
                wtdSumSq += values[i]*values[i]/powf(snp2pq[i], Sgw);
                sum2pqBetaSq += snp2pq[i]*values[i]*values[i];
                for (unsigned j=0; j<snp->numAnnos; ++j) {
                    annoIdx = snp->annoVec[j]->idx;
                    sum2pqBetaSqAnno[annoIdx] += snp2pq[i]*values[i]*values[i];
                }
                ++numNonZeros;
            }
            else {
                if (oldSample) {
                    rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*oldSample;
                }
                values[i] = 0.0;
                for (unsigned j=0; j<snp->numAnnos; ++j) {
                    annoIdx = snp->annoVec[j]->idx;
                    snpIdx = snp->annoIdx[j];
                    valuesPerAnno[annoIdx][snpIdx] = 0.0;
                }
            }
        }
    }
}

void StratApproxBayesS::sampleUnknowns() {
    unsigned cnt=0;
    
    if (model == mixture) {
        do {
            if (sparse) {
                snpEffects.sampleFromFCMixture(rcorr, data.ZPZsp, data.ZPZdiag, data.chromInfoVec, data.incdSnpInfoVec, data.snp2pq,
                                        data.LDsamplVar, data.numAnnos, sigmaSqStrat.values, piStrat.values, Sstrat.values, S.value,
                                        varg.value, vare.value, ps.value, overdispersion);
            } else {
                snpEffects.sampleFromFCMixture(rcorr, data.ZPZ, data.ZPZdiag, data.windStart, data.windSize, data.chromInfoVec, data.incdSnpInfoVec, data.snp2pq,
                                        data.LDsamplVar, data.numAnnos, sigmaSqStrat.values, piStrat.values, Sstrat.values, S.value,
                                        varg.value, vare.value, ps.value, overdispersion);
            }
            if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
        } while (snpEffects.numNonZeros == 0);
        
        snpAnnoMembership.getValues(snpEffects.snpAnnoVec);

    }
    else if (model == linear) {
        do {
            if (sparse) {
                snpEffects.sampleFromFCLinear(rcorr, data.ZPZsp, data.ZPZdiag, data.chromInfoVec, data.incdSnpInfoVec, data.snp2pq,
                                        data.LDsamplVar, data.numAnnos, sigmaSqStrat.values, perSnpPi.values, Sstrat.values, S.value,
                                        varg.value, vare.value, ps.value, overdispersion);
            } else {
                snpEffects.sampleFromFCLinear(rcorr, data.ZPZ, data.ZPZdiag, data.windStart, data.windSize, data.chromInfoVec, data.incdSnpInfoVec, data.snp2pq,
                                        data.LDsamplVar, data.numAnnos, sigmaSqStrat.values, perSnpPi.values, Sstrat.values, S.value,
                                        varg.value, vare.value, ps.value, overdispersion);
            }
            if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
        } while (snpEffects.numNonZeros == 0);

//        if (estimatePi) perSnpPi.sampleFromFC(data.annoMat, data.APA, snpEffects.values);

    }
    
    nnzSnp.getValue(snpEffects.numNonZeros);
    nnzStrat.getValues(snpEffects.numNonZeroPerAnno);
    propNnzStrat.compute(nnzStrat.values, nnzSnp.value);

    if (estimatePi) {
        pi.sampleFromFC(data.numIncdSnps, nnzSnp.value);
        piStrat.sampleFromFC(data.numSnpAnnoVec, nnzStrat.values);
        piEnrich.compute(nnzStrat.values, nnzSnp.value);
        if (model == linear) perSnpPi.computeFromAnnoPi(data.annoMat, piStrat.values);
    }
    
    sigmaSq.sampleFromFC(snpEffects.wtdSumSq, nnzSnp.value);
    sigmaSqStrat.sampleFromFC(snpEffects.wtdSumSqPerAnno, nnzStrat.values);
    sigmaSqEnrich.compute(sigmaSqStrat.values, sigmaSq.value);
    
//    if (model == mixture) {
//        sigmaSqG.compute(sigmaSq.value, snpEffects.sum2pqSplusOne);
//        varg.value = sigmaSqG.value;
//        //    covg.compute(data.ypy, snpEffects.values, data.ZPy, rcorr);
//        //    varg.compute(snpEffects.values, data.ZPy, rcorr, covg.value);
//        vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, covg.value);
//        hsq.compute(varg.value, vare.value);
//        hsqStrat.compute(sigmaSqStrat.values, snpEffects.sum2pqSplusOnePerAnno, varg.value, vare.value);
//        //    hsqStrat.compute(snpEffects.values, data.annowiseZPZsp, data.annowiseZPZdiag, data.annoInfoVec, varg.value, vare.value);
//    }
//    else if (model == linear) {
        varg.value = sigmaSqG.value = snpEffects.sum2pqBetaSq;
        vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, covg.value);
        hsq.compute(varg.value, vare.value);
        hsqStrat.values = snpEffects.sum2pqBetaSqAnno/(varg.value + vare.value);
//    }
    propHsqStrat.compute(hsqStrat.values, hsq.value);
    perSnpHsqEnrich.compute(hsqStrat.values, piEnrich.expectation, hsq.value);
    perNzHsqEnrich.compute(hsqStrat.values, nnzStrat.values, hsq.value, nnzSnp.value);
    
    S.sampleFromFC(snpEffects.wtdSumSq, nnzSnp.value, sigmaSq.value, snpEffects.values, data.snp2pq, snp2pqPowS, logSnp2pq, varg.value, sigmaSq.scale, snpEffects.sum2pqSplusOne);
    if (model == mixture)
        Sstrat.sampleFromFC(snpEffects.valuesPerAnno, nnzStrat.values, sigmaSqStrat.values, hsqStrat.values, varg.value, vare.value, data.annoInfoVec, sigmaSqStrat.scales, snpEffects.sum2pqSplusOnePerAnno);
    else if (model == linear)
        Sstrat.sampleFromFC(snpEffects.valuesPerAnno, nnzStrat.values, sigmaSqStrat.values, sigmaSqStrat.values.cwiseProduct(snpEffects.sum2pqSplusOnePerAnno), varg.value, vare.value, data.annoInfoVec, sigmaSqStrat.scales, snpEffects.sum2pqSplusOnePerAnno);
    Senrich.compute(Sstrat.values, S.value);
    
    if (modelPS) ps.compute(rcorr, data.ZPZdiag, data.LDsamplVar, varg.value, vare.value, data.chisq);
    
    if (sparse)
        rounding.computeRcorr(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    else
        rounding.computeRcorr(data.ZPy, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    
    scaleStrat.values = sigmaSqStrat.scales;
}

void StratApproxBayesS::sampleStartVal(){
    sigmaSqStrat.sampleFromPrior();
    if (estimatePi) piStrat.sampleFromPrior();
    Sstrat.sampleFromPrior();
    cout << "  Starting value for " << sigmaSq.label << ": " << sigmaSqStrat.values.transpose() << endl;
    if (estimatePi) cout << "  Starting value for " << pi.label << ": " << piStrat.values.transpose() << endl;
    cout << "  Starting value for " << S.label << ": " << Sstrat.values.transpose() << endl;
    cout << endl;
}




///// post hoc stratified analysis based on MCMC samples of SNP effects

void PostHocStratifyS::SnpEffects::getValues(const SparseVector<float> &snpEffects, const vector<SnpInfo*> &snpInfoVec, const vector<AnnoInfo*> &annoInfoVec, const VectorXf &snp2pq, const VectorXf &S, const float Sgw) {
    VectorXf beta(snpEffects);
    values = beta;
    long numAnnos = annoInfoVec.size();
    wtdSumSqPerAnno.setZero(numAnnos);
    numNonZeroPerAnno.setZero(numAnnos);
//    numNonZeros = snpEffects.nonZeros();   // THIS IS NOT CORRECT!!
    wtdSumSq = 0;
    for (unsigned i=0; i<numAnnos; ++i) {
        valuesPerAnno[i].setZero();
    }
    numNonZeros= 0;
    for (unsigned i=0; i<size; ++i) {
        if (!values[i]) continue;
        ++numNonZeros;
        wtdSumSq += values[i]*values[i]/powf(snp2pq[i], Sgw);
        SnpInfo *snp = snpInfoVec[i];
        for (unsigned j=0; j<snp->numAnnos; ++j) {
            unsigned annoIdx = snp->annoVec[j]->idx;
            unsigned snpIdx = snp->annoIdx[j];
            valuesPerAnno[annoIdx][snpIdx] = values[i];
            wtdSumSqPerAnno[annoIdx] += values[i]*values[i]/powf(snp2pq[i], S[annoIdx]);
            ++numNonZeroPerAnno[annoIdx];
        }
    }
//    long chunkSize = size/omp_get_max_threads();
//#pragma omp parallel for schedule(dynamic, chunkSize)
//    for (unsigned i=0; i<numAnnos; ++i) {
//        AnnoInfo *anno = annoInfoVec[i];
//        valuesPerAnno[i].setZero(anno->size);
//        for (unsigned j=0; j<anno->size; ++j) {
//            unsigned snpidx = anno->memberSnpVec[j]->index;
//            if (beta[snpidx]) {
//                valuesPerAnno[i][j] = beta[snpidx];
//                wtdSumSqPerAnno[i] += beta[snpidx]*beta[snpidx]/powf(snp2pq[snpidx], S[i]);
//                ++numNonZeroPerAnno[i];
//            }
//        }
//    }
}

void PostHocStratifyS::PiStratified::compute(const vector<unsigned int> &numSnps, const VectorXf &numSnpEff) {
    for (unsigned i=0; i<size; ++i) {
        values[i] = numSnpEff[i]/numSnps[i];
    }
}

void PostHocStratifyS::VarEffectStratified::compute(const VectorXf &snpEffSumSq, const VectorXf &numSnpEff) {
    for (unsigned i=0; i<size; ++i) {
        if (numSnpEff[i]) values[i] = snpEffSumSq[i]/numSnpEff[i];
    }
}

void PostHocStratifyS::sampleUnknowns() {
    snpEffects.getValues(snpEffectsMcmc.datMatSp.row(iter), data.incdSnpInfoVec, data.annoInfoVec, data.snp2pq, Sstrat.values, S.value);
    hsq.value = hsqMcmc.datMat.row(iter*thin)[0];
    varg.value = hsq.value;
    vare.value = 1.0 - hsq.value;

    nnzSnp.getValue(snpEffects.numNonZeros);
    nnzStrat.getValues(snpEffects.numNonZeroPerAnno);
    propNnzStrat.compute(nnzStrat.values, nnzSnp.value);

    pi.compute(data.numIncdSnps, nnzSnp.value);
    piStrat.compute(data.numSnpAnnoVec, nnzStrat.values);
    piEnrich.compute(nnzStrat.values, nnzSnp.value);
    
    sigmaSq.compute(snpEffects.wtdSumSq, nnzSnp.value);
    sigmaSqStrat.compute(snpEffects.wtdSumSqPerAnno, nnzStrat.values);
    sigmaSqEnrich.compute(sigmaSqStrat.values, sigmaSq.value);
    
//    hsq.compute(varg.value, vare.value);
//    hsqStrat.compute(sigmaSqStrat.values, snpEffects.sum2pqSplusOnePerAnno, varg.value, vare.value);
    hsqStrat.compute(snpEffects.valuesPerAnno, data.annowiseZPZsp, data.annowiseZPZdiag, varg.value, vare.value);
    propHsqStrat.compute(hsqStrat.values, hsq.value);
    perSnpHsqEnrich.compute(hsqStrat.values, piEnrich.expectation, hsq.value);
    perNzHsqEnrich.compute(hsqStrat.values, nnzStrat.values, hsq.value, nnzSnp.value);
    
    S.sampleFromFC(snpEffects.wtdSumSq, nnzSnp.value, sigmaSq.value, snpEffects.values, data.snp2pq, snp2pqPowS, logSnp2pq, varg.value, sigmaSq.scale, snpEffects.sum2pqSplusOne);
    Sstrat.sampleFromFC(snpEffects.valuesPerAnno, nnzStrat.values, sigmaSqStrat.values, hsqStrat.values, varg.value, vare.value, data.annoInfoVec, sigmaSqStrat.scales, snpEffects.sum2pqSplusOnePerAnno);
    Senrich.compute(Sstrat.values, S.value);
    
    ++iter;
}


void PostHocStratifySMix::DeltaS::getValues(const SparseVector<float> &deltaS, const vector<AnnoInfo*> &annoInfoVec){
    VectorXf delta(deltaS);
    values = delta;
    sum = values.sum();
    
    long numAnnos = annoInfoVec.size();
    sumPerAnno.setZero(numAnnos);
    
    long chunkSize = size/omp_get_max_threads();
#pragma omp parallel for schedule(dynamic, chunkSize)
    for (unsigned i=0; i<numAnnos; ++i) {
        AnnoInfo *anno = annoInfoVec[i];
        for (unsigned j=0; j<anno->size; ++j) {
            unsigned snpidx = anno->memberSnpVec[j]->index;
            sumPerAnno[i] += values[snpidx];
        }
    }
}


void PostHocStratifySMix::sampleUnknowns() {
    snpEffects.getValues(snpEffectsMcmc.datMatSp.row(iter), data.incdSnpInfoVec, data.annoInfoVec, data.snp2pq, Sstrat.values, S.value);
    deltaS.getValues(deltaSmcmc.datMatSp.row(iter), data.annoInfoVec);
    hsq.value = hsqMcmc.datMat.row(iter*thin)[0];
    varg.value = hsq.value;
    vare.value = 1.0 - hsq.value;
    
    nnzSnp.getValue(snpEffects.numNonZeros);
    nnzStrat.getValues(snpEffects.numNonZeroPerAnno);
    propNnzStrat.compute(nnzStrat.values, nnzSnp.value);
    
    pi.compute(data.numIncdSnps, nnzSnp.value);
    piStrat.compute(data.numSnpAnnoVec, nnzStrat.values);
    piEnrich.compute(nnzStrat.values, nnzSnp.value);
    
//    hsq.compute(varg.value, vare.value);
    hsqStrat.compute(snpEffects.valuesPerAnno, data.annowiseZPZsp, data.annowiseZPZdiag, varg.value, vare.value);
    propHsqStrat.compute(hsqStrat.values, hsq.value);
    perSnpHsqEnrich.compute(hsqStrat.values, piEnrich.expectation, hsq.value);
    perNzHsqEnrich.compute(hsqStrat.values, nnzStrat.values, hsq.value, nnzSnp.value);
    
    piS.compute(data.numIncdSnps, deltaS.sum);
    piSstrat.compute(data.numSnpAnnoVec, deltaS.sumPerAnno);
    piSenrich.compute(deltaS.sumPerAnno, deltaS.sum);
    
    ++iter;
}



