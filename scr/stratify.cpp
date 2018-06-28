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

void StratApproxBayesS::VarEffectEnrichment::compute(const VectorXf &sigmaSqStrat, const float sigmaSq) {
    for (unsigned i=0; i<size; ++i) {
        values[i] = sigmaSqStrat[i] > sigmaSq ? 1 : 0;
    }
}

void StratApproxBayesS::PiStratified::sampleFromFC(const vector<unsigned> &numSnps, const VectorXf &numSnpEff) {
    for (unsigned i=0; i<size; ++i) {
        float alphaTilde = numSnpEff[i] + alpha;
        float betaTilde  = numSnps[i] - numSnpEff[i] + beta;
        values[i] = Beta::sample(alphaTilde, betaTilde);
    }
}

void StratApproxBayesS::PiEnrichment::compute(const VectorXf &piStrat, const float pi) {
    for (unsigned i=0; i<size; ++i) {
        values[i] = piStrat[i] > pi ? 1 : 0;
    }
}

void StratApproxBayesS::HeritabilityStratified::compute(const float genVar, const float resVar, const VectorXf &snpEffects,
                                                        const vector<SparseVector<float> > &ZPZsp, const vector<AnnoInfo *> &annoInfoVec) {
    for (unsigned i=0; i<size; ++i) {
        AnnoInfo *anno = annoInfoVec[i];
        float varg = 0.0;
        unsigned jj, kk;
        for (unsigned j=0; j<anno->size; ++j) {
            jj = anno->memberSnpVec[j]->index;
            if (snpEffects[jj]) {
                for (unsigned k=j; k<anno->size; ++k) {
                    kk = anno->memberSnpVec[k]->index;
                    if (snpEffects[kk]) {
                        if (jj==kk) {
                            varg += snpEffects[jj]*snpEffects[kk]*ZPZsp[jj].coeff(kk);
                        } else {
                            varg += 2.0*snpEffects[jj]*snpEffects[kk]*ZPZsp[jj].coeff(kk);
                        }
                    }
                }
            }
        }
        varg /= float(sampleSize);
        values[i] = varg/(genVar + resVar);
        //        cout << "shsq " << i << " " << values[i] << endl;
    }
}

void StratApproxBayesS::HeritabilityStratified::compute(const VectorXf &snpEffects, const vector<SparseMatrix<float> > &annowiseZPZsp, const vector<VectorXf> &annowiseZPZdiag, const vector<AnnoInfo*> &annoInfoVec, const float genVar, const float resVar) {
    for (unsigned i=0; i<size; ++i) {
        AnnoInfo *anno = annoInfoVec[i];
        SnpInfo *snp;
        VectorXf beta(anno->size);
        for (unsigned j=0; j<anno->size; ++j) {
            snp = anno->memberSnpVec[j];
            beta[j] = snpEffects[snp->index];
        }
        float varg = beta.dot(annowiseZPZdiag[i].cwiseProduct(beta)) + 2.0*beta.dot(annowiseZPZsp[i]*beta);
        varg /= float(sampleSize);
        values[i] = varg/(genVar + resVar);
    }
}

void StratApproxBayesS::HeritabilityEnrichment::compute(const VectorXf &hsqStrat, const float hsqTotal) {
    float averageHsq = hsqTotal/float(numSnps);
    for (unsigned i=0; i<size; ++i) {
        values[i] = hsqStrat[i]/float(numSnpAnno[i]) > averageHsq ? 1 : 0;
    }
}

void StratApproxBayesS::SpStratified::sampleFromFC(const VectorXf &snpEffects, const VectorXf &numNonZeros, const VectorXf &snp2pq, const VectorXf &sigmaSq, const VectorXf &hsq, const float genVar, const float resVar, const vector<AnnoInfo*> &annoInfoVec, VectorXf &scales) {
    
    VectorXf varg = hsq.array()*(genVar + resVar);
    
    unsigned idx;
    for (unsigned i=0; i<size; ++i) {
        unsigned numNonZeroAnnoi = numNonZeros[i];
        if (numNonZeroAnnoi < 3) {   // do not estimate S if few non-zero effects in this annotation
            values[i] = 0.0;
            continue;
        }
        VectorXf snpEffectsAnnoi(numNonZeroAnnoi);
        VectorXf snp2pqAnnoi(numNonZeroAnnoi);
        VectorXf snp2pqLogAnnoi(numNonZeroAnnoi);
        AnnoInfo *anno = annoInfoVec[i];
        unsigned ii = 0;
        for (unsigned j=0; j<anno->size; ++j) {
            idx = anno->memberSnpVec[j]->index;
            if (snpEffects[idx]) {
                snpEffectsAnnoi[ii] = snpEffects[idx];
                snp2pqAnnoi[ii] = snp2pq[idx];
                snp2pqLogAnnoi[ii] = snp2pqLog[idx];
                ++ii;
            }
        }
        
        hmcSampler(i, snpEffectsAnnoi, snp2pqAnnoi, snp2pqLogAnnoi, sigmaSq[i], varg[i], scales[i], values[i]);
        
    }
}

void StratApproxBayesS::SpStratified::hmcSampler(const unsigned annoIdx, const VectorXf &snpEffects, const VectorXf &snp2pq, const VectorXf &snp2pqLog, const float sigmaSq, const float varg, float &scale, float &value) {
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
        scale = 0.5*varg/snp2pq.array().pow(value+1.0).sum();
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

void StratApproxBayesS::SpEnrichment::compute(const VectorXf &Sstrat, const float S) {
    for (unsigned i=0; i<size; ++i) {
        values[i] = Sstrat[i] < S ? 1 : 0;   // lower than S since mostly S is negative
    }
}

void StratApproxBayesS::SnpEffects::sampleFromFC(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag, const vector<ChromInfo *> &chromInfoVec, const vector<SnpInfo *> &incdSnpInfoVec, const VectorXf &snp2pq, const VectorXf &LDsamplVar, const unsigned numAnnos, const VectorXf &sigmaSq, const VectorXf &pi, const VectorXf &S, const float varg, const float vare, const float ps, const float overdispersion) {
    
    long numChr = chromInfoVec.size();
    
    wtdSumSq.setZero(numAnnos);
    numNonZeros.setZero(numAnnos);
    
    VectorXf logPi = pi.array().log();
    VectorXf logPiComp = (1.0-pi.array()).log();
    VectorXf invSigmaSq = sigmaSq.cwiseInverse();
    
    for (unsigned chr=0; chr<numChr; ++chr) {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        
        float oldSample;
        float rhs, invLhs, uhat;
        float logDelta0, logDelta1, probDelta1;
        float snp2pqPowS;
        float varei;
        float sampleDiff;
        
        unsigned annoIdx;
        
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            annoIdx = incdSnpInfoVec[i]->annoPtr[0]->idx;
            oldSample = values[i];
            
            varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
            snp2pqPowS = powf(snp2pq[i], S[annoIdx]);
            
            rhs  = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            invLhs = 1.0f/(ZPZdiag[i]/varei + invSigmaSq[annoIdx]/snp2pqPowS);
            uhat = invLhs*rhs;
            
            logDelta1 = 0.5*(logf(invLhs) - logf(snp2pqPowS*sigmaSq[annoIdx]) + uhat*rhs) + logPi[annoIdx];
            logDelta0 = logPiComp[annoIdx];
            probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
            
            if (bernoulli.sample(probDelta1)) {
                values[i] = normal.sample(uhat, invLhs);
                sampleDiff = oldSample - values[i];
                for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                    rcorr[it.index()] += it.value() * sampleDiff;
                }
                wtdSumSq[annoIdx] += values[i]*values[i]/snp2pqPowS;
                ++numNonZeros[annoIdx];
            } else {
                if(oldSample) {
                    for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                        rcorr[it.index()] += it.value() * oldSample;
                    }
                }
                values[i] = 0.0;
            }
        }
    }
}

void StratApproxBayesS::sampleUnknowns() {
    unsigned cnt=0;
    do {
        snpEffects.sampleFromFC(rcorr, data.ZPZsp, data.ZPZdiag, data.chromInfoVec, data.incdSnpInfoVec, data.snp2pq,
                                data.LDsamplVar, data.numAnnos, sigmaSqStrat.values, piStrat.values, Sstrat.values,
                                varg.value, vare.value, ps.value, overdispersion);
        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros.sum() == 0);
    
    if (estimatePi) {
        pi.sampleFromFC(data.numIncdSnps, snpEffects.numNonZeros.sum());
        piStrat.sampleFromFC(data.numSnpAnnoVec, snpEffects.numNonZeros);
        piEnrich.compute(piStrat.values, pi.value);
    }
    
    nnzSnp.getValue(snpEffects.numNonZeros.sum());
    nnzStrat.getValues(snpEffects.numNonZeros);
    
    sigmaSq.sampleFromFC(snpEffects.wtdSumSq.sum(), snpEffects.numNonZeros.sum());
    sigmaSqStrat.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
    sigmaSqEnrich.compute(sigmaSqStrat.values, sigmaSq.value);
    
    vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr);
    varg.compute(snpEffects.values, data.ZPy, rcorr);
    hsq.compute(varg.value, vare.value);
    hsqStrat.compute(snpEffects.values, annowiseZPZsp, annowiseZPZdiag, data.annoInfoVec, varg.value, vare.value);
    hsqEnrich.compute(hsqStrat.values, hsq.value);
    
    S.sampleFromFC(snpEffects.wtdSumSq.sum(), snpEffects.numNonZeros.sum(), sigmaSq.value, snpEffects.values, data.snp2pq, snp2pqPowS, logSnp2pq, varg.value, sigmaSq.scale, snpEffects.sum2pqOneMinusS);
    Sstrat.sampleFromFC(snpEffects.values, snpEffects.numNonZeros, data.snp2pq, sigmaSqStrat.values, hsqStrat.values,
                        varg.value, vare.value, data.annoInfoVec, sigmaSqStrat.scales);
    Senrich.compute(Sstrat.values, S.value);
    
    if (modelPS) ps.compute(rcorr, data.ZPZdiag, data.LDsamplVar, varg.value, vare.value);
    
    rounding.computeRcorr(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    
    scaleStrat.values = sigmaSqStrat.scales;
}

void StratApproxBayesS::makeAnnowiseSparseLDM(const vector<SparseVector<float> > &ZPZsp, const vector<AnnoInfo *> &annoInfoVec, const vector<SnpInfo*> &snpInfoVec) {
    if (myMPI::rank==0) {
        cout << " making annotation-wise sparse LD matrix ..." << endl;
    }
    long numAnnos = annoInfoVec.size();
    long numSnps = ZPZsp.size();
    annowiseZPZsp.resize(numAnnos);
    annowiseZPZdiag.resize(numAnnos);
    for (unsigned i=0; i<numAnnos; ++i) {
        AnnoInfo *anno = annoInfoVec[i];
        annowiseZPZsp[i].resize(anno->size, anno->size);
        //        VectorXi nnz(anno->size);
        //        for (unsigned j=0; j<anno->size; ++j) {
        //            nnz[j] = anno->memberSnpVec[j]->windSize;
        //        }
        //        annowiseZPZsp[i].reserve(nnz);   // wield - make it slower
        annowiseZPZdiag[i].resize(anno->size);
        SnpInfo *snpj, *snpk;
        unsigned r, c;
        float v;
        for (unsigned j=0; j<anno->size; ++j) {
            if (!(j%10000))
                cout << "  annotaion " << std::setw(2) << std::left << i+1 << " snp " << std::setw(6) << std::left << j << "\r";
            snpj = anno->memberSnpVec[j];
            r = j;
            c = 0;
            VectorXf dense;
            dense.setZero(numSnps);
            for (SparseVector<float>::InnerIterator it(ZPZsp[snpj->index]); it; ++it) {
                dense[it.index()] = it.value();
            }
            for (unsigned k=j+1; k<anno->size; ++k) {
                snpk = anno->memberSnpVec[k];
                if (snpj->chrom != snpk->chrom) continue;
                //                v = ZPZsp[snpj->index].coeff(snpk->index);
                v = dense[snpk->index];
                if (v) {
                    annowiseZPZsp[i].insert(c++, r) = v;
                }
            }
            annowiseZPZdiag[i][j] = dense[snpj->index];
            if (j==anno->size-1)
                cout << "  annotaion " << std::setw(2) << std::left << i+1 << " snp " << std::setw(6) << std::left << j+1 << " nonzeros " << annowiseZPZsp[i].nonZeros() << "\r";
        }
        if (myMPI::rank==0) cout << endl;
        annowiseZPZsp[i].makeCompressed();
    }
}

