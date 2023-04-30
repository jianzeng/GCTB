//
//  model.cpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "model.hpp"


void BayesC::FixedEffects::sampleFromFC(VectorXf &ycorr, const MatrixXf &X,
                                        const VectorXf &XPXdiag, const float vare){
    float rhs;
    for (unsigned i=0; i<size; ++i) {
        if (!XPXdiag[i]) continue;
        float oldSample = values[i];
        float rhs = X.col(i).dot(ycorr);
        rhs += XPXdiag[i]*oldSample;
        float invLhs = 1.0f/XPXdiag[i];
        float bhat = invLhs*rhs;
        values[i] = Normal::sample(bhat, invLhs*vare);
        ycorr += X.col(i) * (oldSample - values[i]);
    }
}

void BayesC::RandomEffects::sampleFromFC(VectorXf &ycorr, const MatrixXf &W, const VectorXf &WPWdiag, const VectorXf &Rsqrt, const bool weightedRes, const float sigmaSqRand, const float vare, VectorXf &rhat){
    rhat.setZero(ycorr.size());
    float invVare = 1.0f/vare;
    float invSigmaSqRand = 1.0f/sigmaSqRand;
    float rhs = 0.0;
    ssq = 0.0;
    for (unsigned i=0; i<size; ++i) {
        if (!WPWdiag[i]) continue;
        float oldSample = values[i];
        float rhs = W.col(i).dot(ycorr) + WPWdiag[i]*oldSample;
        rhs *= invVare;
        float invLhs = 1.0f/(WPWdiag[i]*invVare + invSigmaSqRand);
        float uhat = invLhs*rhs;
        values[i] = Normal::sample(uhat, invLhs);
        ssq += values[i]*values[i];
        if (weightedRes) rhat += W.col(i).cwiseProduct(Rsqrt) * values[i];
        else rhat  += W.col(i) * values[i];
        ycorr += W.col(i) * (oldSample - values[i]);
    }
}

void BayesC::VarRandomEffects::sampleFromFC(const float randEffSumSq, const unsigned int numRandEff){
    float dfTilde = df + numRandEff;
    float scaleTilde = randEffSumSq + df*scale;
    value = InvChiSq::sample(dfTilde, scaleTilde);    
}

void BayesC::SnpEffects::sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes,
const float sigmaSq, const float pi, const float vare, VectorXf &ghat){
    if (algorithm == gibbs) {
        gibbsSampler(ycorr, Z, ZPZdiag, Rsqrt, weightedRes, sigmaSq, pi, vare, ghat);
    } else if (algorithm == hmc) {
        hmcSampler(ycorr, Z, ZPZdiag, sigmaSq, pi, vare, ghat);
    }
}

void BayesC::SnpEffects::gibbsSampler(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes,
                                      const float sigmaSq, const float pi, const float vare, VectorXf &ghat){
    sumSq = 0.0;
    numNonZeros = 0;
    
    pip.setZero(size);
    ghat.setZero(ycorr.size());
    
    float oldSample;
    float rhs, invLhs, uhat;
    float logDelta0, logDelta1, probDelta1;
    float logPi = log(pi);
    float logPiComp = log(1.0-pi);
    float logSigmaSq = log(sigmaSq);
    float invVare = 1.0f/vare;
    float invSigmaSq = 1.0f/sigmaSq;
    
    for (unsigned i=0; i<size; ++i) {
        oldSample = values[i];
        rhs = Z.col(i).dot(ycorr);
        rhs += ZPZdiag[i]*oldSample;
        rhs *= invVare;
        invLhs = 1.0f/(ZPZdiag[i]*invVare + invSigmaSq);
        uhat = invLhs*rhs;
        logDelta1 = 0.5*(logf(invLhs) - logSigmaSq + uhat*rhs) + logPi;
        //logDelta1 = rhs*oldSample - 0.5*ZPZdiag[i]*oldSample*oldSample/vare + logPiComp;
        logDelta0 = logPiComp;
        probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
        pip[i] = probDelta1;
        
        //cout << i << " rhs " << rhs << " invLhs " << invLhs << " uhat " << uhat << endl;

        if (bernoulli.sample(probDelta1)) {
            values[i] = normal.sample(uhat, invLhs);
            ycorr += Z.col(i) * (oldSample - values[i]);
            if (weightedRes) ghat += Z.col(i).cwiseProduct(Rsqrt) * values[i];
            else ghat  += Z.col(i) * values[i];
            sumSq += values[i]*values[i];
            ++numNonZeros;
        } else {
            if (oldSample) ycorr += Z.col(i) * oldSample;
            values[i] = 0.0;
        }
    }
}

void BayesC::SnpEffects::hmcSampler(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const float sigmaSq, const float pi, const float vare, VectorXf &ghat){
    // Hamiltonian Monte Carlo
    // Only BayesC0 model available
    
    float stepSize = 0.1;
    unsigned numSteps = 10;
    
    ycorr += Z*values;
    
    static MatrixXf ZPZ;
    if (cnt==0) ZPZ = Z.transpose()*Z;
    VectorXf ypZ = ycorr.transpose()*Z;
    
    VectorXf curr = values;
    
    ArrayXf curr_p(size);
    for (unsigned i=0; i<size; ++i) {
        curr_p[i] = Stat::snorm();
    }
    
    VectorXf cand = curr;
    // Make a half step for momentum at the beginning
    ArrayXf cand_p = curr_p - 0.5*stepSize * gradientU(curr, ZPZ, ypZ, sigmaSq, vare);
    
    for (unsigned i=0; i<numSteps; ++i) {
        cand.array() += stepSize * cand_p;
        if (i < numSteps-1) {
            cand_p -= stepSize * gradientU(cand, ZPZ, ypZ, sigmaSq, vare);
        } else {
            cand_p -= 0.5*stepSize * gradientU(cand, ZPZ, ypZ, sigmaSq, vare);
        }
    }
    
    float curr_H = computeU(curr, ZPZ, ypZ, sigmaSq, vare) + 0.5*curr_p.matrix().squaredNorm();
    float cand_H = computeU(cand, ZPZ, ypZ, sigmaSq, vare) + 0.5*cand_p.matrix().squaredNorm();
    
    if (Stat::ranf() < exp(curr_H-cand_H)) {  // accept
        values = cand;
        ghat = Z*values;
        ++mhr;
    }
    
    if (!(++cnt % 100)) {
        float ar = mhr/float(cnt);
        if      (ar < 0.5) cout << "Warning: acceptance rate for SNP effects is too low "  << ar << endl;
        else if (ar > 0.9) cout << "Warning: acceptance rate for SNP effects is too high " << ar << endl;
    }
    
    numNonZeros = size;
    sumSq = values.squaredNorm();
    
    ycorr -= Z*values;
}

ArrayXf BayesC::SnpEffects::gradientU(const VectorXf &alpha, const MatrixXf &ZPZ, const VectorXf &ypZ, const float sigmaSq, const float vare){
    return 1.0/vare*(ZPZ*alpha - ypZ) + 1/sigmaSq*alpha;
}

float BayesC::SnpEffects::computeU(const VectorXf &alpha, const MatrixXf &ZPZ, const VectorXf &ypZ, const float sigmaSq, const float vare){
    return 0.5/vare*(alpha.transpose()*ZPZ*alpha + vare/sigmaSq*alpha.squaredNorm() - 2.0*ypZ.dot(alpha));
}

void BayesC::SnpEffects::sampleFromFC_omp(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag,
                                          const float sigmaSq, const float pi, const float vare, VectorXf &ghat){
    // speed-enhanced single site Gibbs sampling due to the use of parallel computing on SNPs with zero effect
    
    unsigned blockSize = 1; //omp_get_max_threads();
    //cout << blockSize << endl;
    
    sumSq = 0.0;
    numNonZeros = 0;
    
    ghat.setZero(ycorr.size());
    
    float oldSample;
    float logPi = log(pi);
    float logPiComp = log(1.0-pi);
    float logSigmaSq = log(sigmaSq);
    float invVare = 1.0f/vare;
    float invSigmaSq = 1.0f/sigmaSq;
    
    vector<int> deltaVec(blockSize);
    vector<float> invLhsVec(blockSize);
    vector<float> uhatVec(blockSize);
    
    unsigned blocki;
    unsigned i, j;
    bool breakFlag;
    
    for (i=0; i<size; ) {
        
        if (blockSize + i < size) {
            blocki = blockSize;
        } else {
            blocki = size - i;
            deltaVec.resize(blocki);
            invLhsVec.resize(blocki);
            uhatVec.resize(blocki);
        }
        
        #pragma omp parallel for
        for (j=0; j<blocki; ++j) {
            float rhsj = (Z.col(i+j).dot(ycorr) + ZPZdiag[i+j]*values[i+j])*invVare;
            invLhsVec[j] = 1.0f/(ZPZdiag[i+j]*invVare + invSigmaSq);
            uhatVec[j] = invLhsVec[j]*rhsj;
            float logDelta0minusDelta1j = logPiComp - (0.5f*(logf(invLhsVec[j]) - logSigmaSq + uhatVec[j]*rhsj) + logPi);
            deltaVec[j] = bernoulli.sample(1.0f/(1.0f + expf(logDelta0minusDelta1j)));
        }
        
        breakFlag = false;
        for (j=0; j<blocki; ++j) {
            if (values[i+j] || deltaVec[j]) {   // need to update ycorr for the first snp who is in the model at either last or this iteration
                i += j;
                breakFlag = true;
                break;
            }
        }
        
        if (breakFlag) {
            oldSample = values[i];
            if (deltaVec[j]) {
                values[i] = normal.sample(uhatVec[j], invLhsVec[j]);
                ycorr += Z.col(i) * (oldSample - values[i]);
                ghat  += Z.col(i) * values[i];
                sumSq += values[i]*values[i];
                ++numNonZeros;
            } else {
                if (oldSample) ycorr += Z.col(i) * oldSample;
                values[i] = 0.0;
            }
            ++i;
        }
        else {
            i += blocki;
        }
    }
}

void BayesC::VarEffects::sampleFromFC(const float snpEffSumSq, const unsigned numSnpEff){
    float dfTilde = df + numSnpEff;
    float scaleTilde = snpEffSumSq + df*scale;
    value = InvChiSq::sample(dfTilde, scaleTilde);
    //cout << "snpEffSumSq " << snpEffSumSq << " scale " << scale << " scaleTilde " << scaleTilde << " dfTilde " << dfTilde << " value " << value << endl;
}

void BayesC::VarEffects::sampleFromPrior(){
    value = InvChiSq::sample(df, scale);
}

void BayesC::VarEffects::computeScale(const float varg, const VectorXf &snp2pq, const float pi){
    if (noscale)
        scale = (df-2)/df * varg/(snp2pq.sum()*pi);
    else
        scale = (df-2)/df * varg/(snp2pq.size()*pi);
}

void BayesC::VarEffects::computeScale(const float varg, const float sum2pq){
        scale = (df-2)/df * varg/sum2pq;
}

void BayesC::VarEffects::compute(const float snpEffSumSq, const float numSnpEff){
    if (numSnpEff) value = snpEffSumSq/numSnpEff;
}

void BayesC::ScaleVar::sampleFromFC(const float sigmaSq, const float df, float &scaleVar){
    float shapeTilde = shape + 0.5*df;
    float scaleTilde = 1.0/(1.0/scale + 0.5*df/sigmaSq);
    value = Gamma::sample(shapeTilde, scaleTilde);
    scaleVar = value;
}

void BayesC::Pi::sampleFromFC(const unsigned numSnps, const unsigned numSnpEff){
    float alphaTilde = numSnpEff + alpha;
    float betaTilde  = numSnps - numSnpEff + beta;
    value = Beta::sample(alphaTilde, betaTilde);
}

void BayesC::Pi::sampleFromPrior(){
    value = Beta::sample(alpha, beta);
}

void BayesC::Pi::compute(const float numSnps, const float numSnpEff){
    value = numSnpEff/numSnps;
}

void BayesC::ResidualVar::sampleFromFC(VectorXf &ycorr){
    float sse = ycorr.squaredNorm();
    float dfTilde = df + nobs;
    float scaleTilde = sse + df*scale;
    value = InvChiSq::sample(dfTilde, scaleTilde);
}

void BayesC::GenotypicVar::compute(const VectorXf &ghat){
    //value = Gadget::calcVariance(ghat);
    float sum = ghat.sum();
    float ssq = ghat.squaredNorm();
    unsigned size = (unsigned)ghat.size();
    float mean = sum/size;
    value = ssq/size - mean*mean;
}

void BayesC::RandomVar::compute(const VectorXf &rhat){
    //value = Gadget::calcVariance(ghat);
    float sum = rhat.sum();
    float ssq = rhat.squaredNorm();
    unsigned size = (unsigned)rhat.size();
    float mean = sum/size;
    value = ssq/size - mean*mean;
}

void BayesC::Rounding::computeYcorr(const VectorXf &y, const MatrixXf &X, const MatrixXf &W, const MatrixXf &Z,
                                    const VectorXf &fixedEffects, const VectorXf &randomEffects, const VectorXf &snpEffects,
                                    VectorXf &ycorr){
    if (count++ % 100) return;
    VectorXf oldYcorr = ycorr;
    ycorr = y - X*fixedEffects;
    if (randomEffects.size()) ycorr -= W*randomEffects;
    for (unsigned i=0; i<snpEffects.size(); ++i) {
        if (snpEffects[i]) ycorr -= Z.col(i)*snpEffects[i];
    }
    float ss = (ycorr - oldYcorr).squaredNorm();
    value = sqrt(ss);
}

void BayesC::sampleUnknowns(){
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }
    unsigned cnt=0;
//    do {
        snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, data.Rsqrt, data.weightedRes, sigmaSq.value, pi.value, vare.value, ghat);
//        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
//    } while (snpEffects.numNonZeros == 0);
    snpPip.getValues(snpEffects.pip);
    sigmaSq.sampleFromFC(snpEffects.sumSq, snpEffects.numNonZeros);
    //scale.sampleFromFC(sigmaSq.value, sigmaSq.df, sigmaSq.scale);
    if (estimatePi) pi.sampleFromFC(snpEffects.size, snpEffects.numNonZeros);
    vare.sampleFromFC(ycorr);
    
    varg.compute(ghat);
    hsq.compute(varg.value, vare.value);
    
    rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);
    nnzSnp.getValue(snpEffects.numNonZeros);
}

void BayesC::sampleStartVal(){
    sigmaSq.sampleFromPrior();
    if (estimatePi) pi.sampleFromPrior();
        cout << "  Starting value for " << sigmaSq.label << ": " << sigmaSq.value << endl;
    if (estimatePi) cout << "  Starting value for " << pi.label << ": " << pi.value << endl;
        cout << endl;
}


void BayesB::SnpEffects::sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes,
                                      const VectorXf &sigmaSq, const float pi, const float vare, VectorXf &ghat){
    numNonZeros = 0;
    
    ghat.setZero(ycorr.size());
    
    float oldSample;
    float rhs, invLhs, uhat;
    float logDelta0, logDelta1, probDelta1;
    float logPi = log(pi);
    float logPiComp = log(1.0-pi);
    float invVare = 1.0f/vare;
    float beta;
    
    for (unsigned i=0; i<size; ++i) {
        oldSample = values[i];
        rhs = Z.col(i).dot(ycorr);
        rhs += ZPZdiag[i]*oldSample;
        rhs *= invVare;
        invLhs = 1.0f/(ZPZdiag[i]*invVare + 1.0f/sigmaSq[i]);
        uhat = invLhs*rhs;
        logDelta1 = 0.5*(logf(invLhs) - logf(sigmaSq[i]) + uhat*rhs) + logPi;
        logDelta0 = logPiComp;
        probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
        
        //cout << i << " rhs " << rhs << " invLhs " << invLhs << " uhat " << uhat << endl;
        
        if (bernoulli.sample(probDelta1)) {
            values[i] = normal.sample(uhat, invLhs);
            ycorr += Z.col(i) * (oldSample - values[i]);
            if (weightedRes) ghat += Z.col(i).cwiseProduct(Rsqrt) * values[i];
            else ghat  += Z.col(i) * values[i];
            betaSq[i] = values[i]*values[i];
            ++numNonZeros;
        } else {
            if (oldSample) ycorr += Z.col(i) * oldSample;
            beta = normal.sample(0, sigmaSq[i]);
            betaSq[i] = beta*beta;
            values[i] = 0.0;
        }
    }
}

void BayesB::VarEffects::sampleFromFC(const VectorXf &betaSq){
    float dfTilde = df + 1.0f;
    ArrayXf scaleTilde = betaSq.array() + df*scale;
    for (unsigned i=0; i<size; ++i) {
        values[i] = InvChiSq::sample(dfTilde, scaleTilde[i]);
    }
}

void BayesB::sampleUnknowns(){
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }
    unsigned cnt=0;
    do {
        snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, data.Rsqrt, data.weightedRes, sigmaSq.values, pi.value, vare.value, ghat);
        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros == 0);
    sigmaSq.sampleFromFC(snpEffects.betaSq);
    if (estimatePi) pi.sampleFromFC(snpEffects.size, snpEffects.numNonZeros);
    vare.sampleFromFC(ycorr);
    varg.compute(ghat);
    hsq.compute(varg.value, vare.value);
    rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);
    nnzSnp.getValue(snpEffects.numNonZeros);
}


void BayesN::SnpEffects::sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes,
                                      const float sigmaSq, const float pi, const float vare, VectorXf &ghat){
    sumSq = 0.0;
    numNonZeros = 0;
    numNonZeroWind = 0;
    
    ghat.setZero(ycorr.size());
    
    pip.setZero(size);
    windPip.setZero(numWindows);
    
    float oldSample;
    float rhs, invLhs, uhat;
    float logDelta0, logDelta1, probDelta1;
    float logPi = log(pi);
    float logPiComp = log(1.0-pi);
    float logSigmaSq = log(sigmaSq);
    float invVare = 1.0f/vare;
    float invSigmaSq = 1.0f/sigmaSq;
    float diffQuadSum;
    float logDelta0MinusLogDelta1;
    
    unsigned start, end;
    
    for (unsigned i=0; i<numWindows; ++i) {
        start = windStart[i];
        end = i+1 < numWindows ? windStart[i] + windSize[i] : size;
        
        // sample window delta
        diffQuadSum = 0.0;
        if (windDelta[i]) {
            for (unsigned j=start; j<end; ++j) {
                if (snpDelta[j]) {
                    rhs = Z.col(j).dot(ycorr);
                    diffQuadSum += 2.0f*beta[j]*rhs + beta[j]*beta[j]*ZPZdiag[j];
                }
            }
        } else {
            for (unsigned j=start; j<end; ++j) {
                if (snpDelta[j]) {
                    rhs = Z.col(j).dot(ycorr);
                    diffQuadSum += 2.0f*beta[j]*rhs - beta[j]*beta[j]*ZPZdiag[j];
                }
            }
        }
        
        diffQuadSum *= invVare;
        logDelta0MinusLogDelta1 = -0.5f*diffQuadSum + logPiComp - logPi;
        probDelta1 = 1.0f/(1.0f + expf(logDelta0MinusLogDelta1));
        windPip[i] = probDelta1;
        
        if (bernoulli.sample(probDelta1)) {
            if (!windDelta[i]) {
                for (unsigned j=start; j<end; ++j) {
                    if (snpDelta[j]) {
                        ycorr -= Z.col(j) * beta[j];
                    }
                }
            }
            windDelta[i] = 1.0;
            ++numNonZeroWind;
            
            for (unsigned j=start; j<end; ++j) {
                oldSample = beta[j]*snpDelta[j];
                rhs = Z.col(j).dot(ycorr);
                rhs += ZPZdiag[j]*oldSample;
                rhs *= invVare;
                invLhs = 1.0f/(ZPZdiag[j]*invVare + invSigmaSq);
                uhat = invLhs*rhs;
                logDelta1 = 0.5*(logf(invLhs) - logSigmaSq + uhat*rhs) + logLocalPi[i];
                logDelta0 = logLocalPiComp[i];
                probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
                pip[j] = probDelta1;
                if (bernoulli.sample(probDelta1)) {
                    values[j] = beta[j] = normal.sample(uhat, invLhs);
                    ycorr += Z.col(j) * (oldSample - values[j]);
                    if (weightedRes) ghat += Z.col(j).cwiseProduct(Rsqrt) * values[j];
                    else ghat  += Z.col(j) * values[j];
                    sumSq += values[j]*values[j];
                    snpDelta[j] = 1.0;
                    ++cumDelta[j];
                    ++numNonZeros;
                } else {
                    if (oldSample) ycorr += Z.col(j) * oldSample;
                    beta[j] = normal.sample(0.0, sigmaSq);
                    snpDelta[j] = 0.0;
                    values[j] = 0.0;
                }
                //sumSq += beta[j]*beta[j];
            }
        }
        else {
//            unsigned windSize = end-start;
//            float localSum = cumDelta.segment(start,windSize).sum();
            for (unsigned j=start; j<end; ++j) {
                beta[j] = normal.sample(0.0, sigmaSq);
                snpDelta[j] = bernoulli.sample(localPi[i]);
                pip[j] = localPi[i];
//                float seudopi = (localPi[i]/(windSize-1)+cumDelta[j])/(localPi[i]+localSum-cumDelta[j]);
//                snpDelta[j] = bernoulli.sample(seudopi);
                if (values[j]) ycorr += Z.col(j) * values[j];
                values[j] = 0.0;
                //sumSq += beta[j]*beta[j];
            }
            windDelta[i] = 0.0;
        }
    }
}

void BayesN::sampleUnknowns(){
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }
    unsigned cnt=0;
    do {
        snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, data.Rsqrt, data.weightedRes, sigmaSq.value, pi.value, vare.value, ghat);
        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros == 0);
    sigmaSq.sampleFromFC(snpEffects.sumSq, snpEffects.numNonZeros);
    //scale.sampleFromFC(sigmaSq.value, sigmaSq.df, sigmaSq.scale);
    if (estimatePi) pi.sampleFromFC(snpEffects.numWindows, snpEffects.numNonZeroWind);
    vare.sampleFromFC(ycorr);
    
    varg.compute(ghat);
    hsq.compute(varg.value, vare.value);
    
    rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);
    nnzSnp.getValue(snpEffects.numNonZeros);
    nnzWind.getValue(snpEffects.numNonZeroWind);
    windDelta.getValues(snpEffects.windDelta);
}

// ----------------------------------------------------------------------------------------
// Bayes R
// ----------------------------------------------------------------------------------------

void BayesR::ProbMixComps::sampleFromFC(const VectorXf &snpStore) {
	VectorXf dirx;
	dirx = snpStore + alphaVec;
    values = Dirichlet::sample(ndist, dirx);
    for (unsigned i=0; i<ndist; ++i) {
      (*this)[i]->value=values[i];  
    }
}

void BayesR::NumSnpMixComps::getValues(const VectorXf &snpStore) {
    values = snpStore;
    for (unsigned i=0; i<ndist; ++i) {
        (*this)[i]->value=values[i];
    }
}

void BayesR::VgMixComps::compute(const VectorXf &snpEffects, const MatrixXf &Z, const vector<vector<unsigned> > snpset, const float varg) {
    values.setZero(ndist);
    long nobs = Z.rows();
//    for (unsigned k=0; k<ndist; ++k) {
//        if (k!=zeroIdx && k!=minIdx) {
//            long numSnps = snpset[k].size();
//            unsigned idx;
//            VectorXf ghat;
//            ghat.setZero(nobs);
//            for (unsigned i=0; i<numSnps; ++i) {
//                idx = snpset[k][i];
//                ghat += snpEffects[idx]*Z.col(idx);
//            }
//            (*this)[k]->value = values[k] = Gadget::calcVariance(ghat)/varg;
//        }
//    }
//    float sum = values.sum();
//    (*this)[minIdx]->value = values[minIdx] = 1.0 - sum;

    for (unsigned k=0; k<ndist; ++k) {
        if (k!=zeroIdx) {
            long numSnps = snpset[k].size();
            unsigned idx;
            VectorXf ghat;
            ghat.setZero(nobs);
            for (unsigned i=0; i<numSnps; ++i) {
                idx = snpset[k][i];
                ghat += snpEffects[idx]*Z.col(idx);
            }
            values[k] = Gadget::calcVariance(ghat);
        }
    }
    float sum = values.sum();
    for (unsigned k=0; k<ndist; ++k) {
        (*this)[k]->value = values[k] = values[k]/sum;
    }

}

void BayesR::SnpEffects::sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes,
                                      const float sigmaSq, const VectorXf &pis, const VectorXf &gamma,
                                      const float vare, VectorXf &ghat, VectorXf &snpStore,
                                      const float varg, const bool originalModel){
    sumSq = 0.0;
    numNonZeros = 0;
    
    ghat.setZero(ycorr.size());
    float oldSample;
    float rhs;
    // -----------------------------------------
    // Initialise the parameters in MCMC sampler
    // -----------------------------------------
    // ----------------
    // Bayes R specific
    // ----------------
    int ndist, indistflag;
    double v1,  b_ls, ssculm, r;
    VectorXf gp, ll, ll2, pll, snpindist, var_b_ls;
    ndist = pis.size();
    snpStore.setZero(pis.size());
    pll.setZero(pis.size());
    // --------------------------------------------------------------------------------
    // Scale the variances in each of the normal distributions by the genetic variance
    // and initialise the class membership probabilities
    // --------------------------------------------------------------------------------
    if (originalModel)
        gp = gamma * 0.01 * varg;
    else
        gp = gamma * sigmaSq;
//    cout << varg << " " << gp.transpose() << endl;
    snpset.resize(ndist);
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
    }
    
    for (unsigned i=0; i<size; ++i) {
        // ------------------------------
        // Derived Bayes R implementation
        // ------------------------------
        // ----------------------------------------------------
        // Add back the content for the corrected rhs for SNP k
        // ----------------------------------------------------
        rhs = Z.col(i).dot(ycorr);
        oldSample = values[i];
        rhs += ZPZdiag[i] * oldSample;
        // ------------------------------------------------------
        // Calculate the beta least squares updates and variances
        // ------------------------------------------------------
        b_ls = rhs / ZPZdiag[i];
        var_b_ls = gp.array() + vare / ZPZdiag[i];
        // ------------------------------------------------------
        // Calculate the likelihoods for each distribution
        // ------------------------------------------------------
        // ll  = (-1.0 / 2.0) * var_b_ls.array().log()  - (b_ls * b_ls)  / (2 * var_b_ls.array());
        ll = (-1.0 / 2.0) * var_b_ls.array().log()  - (b_ls * b_ls)  / (2 * var_b_ls.array()) + pis.array().log();
        // --------------------------------------------------------------
        // Calculate probability that snp is in each of the distributions
        // in this iteration
        // --------------------------------------------------------------
        // pll = (ll.array().exp().cwiseProduct(pis.array())) / ((ll.array().exp()).cwiseProduct(pis.array())).sum();
        for (unsigned k=0; k<pis.size(); ++k) {
            pll[k] = 1.0 / (exp(ll.array() - ll[k])).sum();
        }
        // --------------------------------------------------------------
        // Sample the group based on the calculated probabilities
        // --------------------------------------------------------------
        ssculm = 0.0;
        r = Stat::ranf();
        indistflag = 1;
        for (int kk = 0; kk < ndist; kk++)
        {
            ssculm += pll(kk);
            if (r < ssculm)
            {
                indistflag = kk + 1;
                snpStore(kk) = snpStore(kk) + 1; 
                break;
            }
        }
        snpset[indistflag-1].push_back(i);
        // --------------------------------------------------------------
        // Sample the effect given the group and adjust the rhs
        // --------------------------------------------------------------
        if (indistflag != 1)
        {
            v1 = ZPZdiag[i] + vare / gp((indistflag - 1));
            values[i] = normal.sample(rhs / v1, vare / v1);
            ycorr += Z.col(i) * (oldSample - values[i]);
            if (weightedRes) ghat += Z.col(i).cwiseProduct(Rsqrt) * values[i];
            else ghat  += Z.col(i) * values[i];
            sumSq += (values[i] * values[i]) / gamma[indistflag - 1];
            ++numNonZeros;
        } else {
            if (oldSample) ycorr += Z.col(i) * oldSample;
            values[i] = 0.0;
        }
    }
}

void BayesR::VarEffects::computeScale(const float varg, const VectorXf &snp2pq, const VectorXf &gamma, const VectorXf &pi){
    if (noscale)
        scale = (df-2)/df * varg/(snp2pq.sum()*gamma.dot(pi));
    else
        scale = (df-2)/df * varg/(snp2pq.size()*gamma.dot(pi));
}

void BayesR::sampleUnknowns(){
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }
    unsigned cnt=0;
    do {
        snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, data.Rsqrt, data.weightedRes, sigmaSq.value, Pis.values, gamma.values, vare.value, ghat, snpStore, varg.value, originalModel);
        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros == 0);  
    sigmaSq.sampleFromFC(snpEffects.sumSq, snpEffects.numNonZeros);
    vare.sampleFromFC(ycorr);
    Pis.sampleFromFC(snpStore);
    numSnps.getValues(snpStore);
    varg.compute(ghat);
    hsq.compute(varg.value, vare.value);
    if (originalModel) Vgs.compute(snpEffects.values, data.Z, snpEffects.snpset, varg.value);
    rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);
    nnzSnp.getValue(snpEffects.numNonZeros);
}


void BayesS::AcceptanceRate::count(const bool state, const float lower, const float upper){
    accepted += state;
    value = accepted/float(++cnt);
    if (!state) ++consecRej;
    else consecRej = 0;
//    if (!(cnt % 100) && myMPI::rank==0) {
//        if      (value < lower) cout << "Warning: acceptance rate is too low  " << value << endl;
//        else if (value > upper) cout << "Warning: acceptance rate is too high " << value << endl;
//    }
}

void BayesS::Sp::sampleFromFC(const float snpEffWtdSumSq, const unsigned numNonZeros, float &sigmaSq, const VectorXf &snpEffects,
                              const VectorXf &snp2pq, ArrayXf &snp2pqPowS, const ArrayXf &logSnp2pq,
                              const float vg, float &scale, float &sum2pqSplusOne){
    if (algorithm == random_walk) {
        randomWalkMHsampler(snpEffWtdSumSq, numNonZeros, sigmaSq, snpEffects, snp2pq, snp2pqPowS, logSnp2pq, vg, scale, sum2pqSplusOne);
    } else if (algorithm == hmc) {
        hmcSampler(numNonZeros, sigmaSq, snpEffects, snp2pq, snp2pqPowS, logSnp2pq, vg, scale, sum2pqSplusOne);
    } else if (algorithm == reg) {
        regression(snpEffects, logSnp2pq, snp2pqPowS, sigmaSq);
    }
}

void BayesS::Sp::sampleFromPrior(){
    value = sample(mean, var);
}

void BayesS::Sp::randomWalkMHsampler(const float snpEffWtdSumSq, const unsigned numNonZeros, const float sigmaSq, const VectorXf &snpEffects,
                                     const VectorXf &snp2pq, ArrayXf &snp2pqPowS, const ArrayXf &logSnp2pq,
                                     const float vg, float &scale, float &sum2pqSplusOne){
    // Random walk Mentroplis-Hastings
    // note that the scale factor of sigmaSq will be simultaneously updated
    
    float curr = value;
    float cand = sample(value, varProp);
    
    float sumLog2pq = 0;
    float snpEffWtdSumSqCurr = snpEffWtdSumSq;
    float snpEffWtdSumSqCand = 0;
    float snp2pqCand = 0;
    float sum2pqCandPlusOne = 0;
    for (unsigned i=0; i<numSnps; ++i) {
        if (snpEffects[i]) {
            sumLog2pq += logf(snp2pq[i]);
            snp2pqCand = powf(snp2pq[i], cand);
            snpEffWtdSumSqCand += snpEffects[i]*snpEffects[i]/snp2pqCand;
            sum2pqCandPlusOne += snp2pq[i]*snp2pqCand;
        }
    }
        
    float logCurr = -0.5f*(curr*sumLog2pq + snpEffWtdSumSqCurr/sigmaSq + curr*curr/var);
    float logCand = -0.5f*(cand*sumLog2pq + snpEffWtdSumSqCand/sigmaSq + cand*cand/var);
    
    //cout << "curr " << curr << " logCurr " << logCurr << " cand " << cand << " logCand " << logCand << " sigmaSq " << sigmaSq << endl;

    float scaleCurr = scale;
    float scaleCand = 0.5f*vg/sum2pqCandPlusOne; // based on the mean of scaled inverse chisq distribution

    // terms due to scale factor of scaled-inverse chi-square distribution
//    float logChisqCurr = 2.0f*log(scaleCurr) - 2.0f*scaleCurr/sigmaSq;
//    float logChisqCand = 2.0f*log(scaleCand) - 2.0f*scaleCand/sigmaSq;
    
    //cout << "curr " << curr << " logChisqCurr " << logChisqCurr << " cand " << cand << " logChisqCand " <<  logChisqCand << endl;
    //cout << "scaleCurr " << scaleCurr << " scaleCand " << scaleCand << endl;
    
//    if (abs(logCand-logCurr) > abs(logChisqCand-logChisqCurr)*10) {  // to avoid the prior of variance dominating the posterior when number of nonzeros are very small
//        logCurr += logChisqCurr;
//        logCand += logChisqCand;
//    }
    
    //cout << "prob " << exp(logCand-logCurr) << endl;
    
    if (Stat::ranf() < exp(logCand-logCurr)) {  // accept
        value = cand;
        scale = scaleCand;
        snp2pqPowS = snp2pq.array().pow(cand);
        sum2pqSplusOne = sum2pqCandPlusOne;
        ar.count(1, 0.1, 0.5);
    } else {
        ar.count(0, 0.1, 0.5);
    }
    
    if (!(ar.cnt % 10)) {
        if      (ar.value < 0.2) varProp *= 0.8;
        else if (ar.value > 0.5) varProp *= 1.2;
    }
    
    tuner.value = varProp;
}

void BayesS::Sp::hmcSampler(const unsigned numNonZeros, const float sigmaSq, const VectorXf &snpEffects,
                            const VectorXf &snp2pq, ArrayXf &snp2pqPowS, const ArrayXf &logSnp2pq,
                            const float vg, float &scale, float &sum2pqSplusOne){
    // Hamiltonian Monte Carlo
    // note that the scale factor of sigmaSq will be simultaneously updated
    
    // Cautious:
    // The sampled value of SNP effect can be exactly zero even it is in the model. In this case, the numNonZeros will be inflated and cause zero element at the end of snp2pqDelta1 vector.
    // To get around this, recalculate numNonZeros here.
    
    unsigned nnz = 0;
    for (unsigned i=0; i<numSnps; ++i)
        if (snpEffects[i]) ++nnz;
    
    // Prepare
    ArrayXf snpEffectDelta1(nnz);
    ArrayXf snp2pqDelta1(nnz);
    ArrayXf logSnp2pqDelta1(nnz);
    
    for (unsigned i=0, j=0; i<numSnps; ++i) {
        if (snpEffects[i]) {
            snpEffectDelta1[j] = snpEffects[i];
            snp2pqDelta1[j] = snp2pq[i];
            logSnp2pqDelta1[j] = logSnp2pq[i];
            ++j;
        }
    }
    
    float snp2pqLogSumDelta1 = logSnp2pqDelta1.sum();
    
    float curr = value;
    float curr_p = Stat::snorm();
    
    float cand = curr;
    // Make a half step for momentum at the beginning
    float cand_p = curr_p - 0.5*stepSize * gradientU(curr,  snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, sigmaSq, vg);
    
    for (unsigned i=0; i<numSteps; ++i) {
        // Make a full step for the position
        cand += stepSize * cand_p;
        if (i < numSteps-1) {
            // Make a full step for the momentum, except at end of trajectory
            cand_p -= stepSize * gradientU(cand, snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, sigmaSq, vg);
        } else {
            // Make a half step for momentum at the end
            cand_p -= 0.5*stepSize * gradientU(cand, snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, sigmaSq, vg);
        }
        //cout << i << " " << cand << endl;
    }

    // Evaluate potential (negative log posterior) and kinetic energies at start and end of trajectory
    float scaleCurr, scaleCand;
    float curr_U_chisq, cand_U_chisq;
    float curr_H = computeU(curr, snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, sigmaSq, vg, scaleCurr, curr_U_chisq) + 0.5*curr_p*curr_p;
    float cand_H = computeU(cand, snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, sigmaSq, vg, scaleCand, cand_U_chisq) + 0.5*cand_p*cand_p;
    
//    if (abs(curr_H-cand_H) > abs(curr_U_chisq-cand_U_chisq)*10) { // temporary fix to avoid the prior of variance dominating the posterior (especially when number of nonzeros are very small)
//        curr_H += curr_U_chisq;
//        cand_H += cand_U_chisq;
//    }
    
    //cout << " curr " << curr << " curr_H " << curr_H << " curr_U " << curr_H - 0.5*curr_p*curr_p << " curr_p " << 0.5*curr_p*curr_p << " curr_scale " << scaleCurr << " sigmaSq " << sigmaSq << endl;
    //cout << " cand " << cand << " cand_H " << cand_H << " cand_U " << cand_H - 0.5*cand_p*cand_p << " curr_p " << 0.5*cand_p*cand_p << " cand_scale " << scaleCand << endl;
    //cout << "curr_H-cand_H " << curr_H-cand_H << endl;
    
    if (Stat::ranf() < exp(curr_H-cand_H)) {  // accept
        value = cand;
        scale = scaleCand;
        snp2pqPowS = snp2pq.array().pow(cand);
        sum2pqSplusOne = snp2pqDelta1.pow(1.0+value).sum();
        ar.count(1, 0.5, 0.9);
    } else {
        ar.count(0, 0.5, 0.9);
    }
    
    if (!(ar.cnt % 10)) {
        if      (ar.value < 0.6) stepSize *= 0.8;
        else if (ar.value > 0.8) stepSize *= 1.2;
    }

    if (ar.consecRej > 20) stepSize *= 0.8;

    tuner.value = stepSize;
}

float BayesS::Sp::gradientU(const float S, const ArrayXf &snpEffects, const float snp2pqLogSum, const ArrayXf &snp2pq, const ArrayXf &logSnp2pq, const float sigmaSq, const float vg){
    // compute the first derivative of the negative log posterior    
    long size = snp2pq.size();
    long chunkSize = size/omp_get_max_threads();
    ArrayXf snp2pqPowS(size);
#pragma omp parallel for schedule(dynamic, chunkSize)
    for (unsigned i=0; i<size; ++i) {
        snp2pqPowS[i] = powf(snp2pq[i], S);
    }
//    ArrayXf snp2pqPowS = snp2pq.pow(S);
    float constantA = snp2pqLogSum;
    float constantB = (snpEffects.square()*logSnp2pq/snp2pqPowS).sum();
    //float constantC = (snp2pq/snp2pqPowS).sum();
    //float constantD = (logSnp2pq*snp2pq/snp2pqPowS).sum();
    float ret = 0.5*constantA - 0.5/sigmaSq*constantB + S/var;
    //float dchisq = - 2.0/constantC*constantD + vg/(sigmaSq*constantC*constantC)*constantD;
    //ret += dchisq;
    //cout << ret << " " << dchisq << endl;
    return ret;
}

float BayesS::Sp::computeU(const float S, const ArrayXf &snpEffects, const float snp2pqLogSum, const ArrayXf &snp2pq, const ArrayXf &logSnp2pq, const float sigmaSq, const float vg, float &scale, float &U_chisq){
    // compute negative log posterior and scale
    ArrayXf snp2pqPowS = snp2pq.pow(S);
    float constantA = snp2pqLogSum;
    float constantB = (snpEffects.square()/snp2pqPowS).sum();
    float constantC = (snp2pq*snp2pqPowS).sum();
    // cout << "Hello I'm in computeU" << endl;
    scale = 0.5*vg/constantC;
    float ret = 0.5*S*constantA + 0.5/sigmaSq*constantB + 0.5*S*S/var;
    U_chisq = 2.0*logf(constantC) + scale/sigmaSq;
    //cout << abs(ret) << " " << dchisq << endl;
    //if (abs(ret) > abs(dchisq)) ret += dchisq;
    return ret;
}

void BayesS::Sp::regression(const VectorXf &snpEffects, const ArrayXf &logSnp2pq, ArrayXf &snp2pqPowS, float &sigmaSq){
    unsigned nnz = 0;
    for (unsigned i=0; i<numSnps; ++i)
        if (snpEffects[i]) ++nnz;
    
    VectorXf y(nnz);
    MatrixXf X(nnz, 2);
    X.col(0) = VectorXf::Ones(nnz);

    for (unsigned i=0, j=0; i<numSnps; ++i) {
        if (snpEffects[i]) {
            y[j]  = snpEffects[i];
            X(j,1) = logSnp2pq[i];
            ++j;
        }
    }

    VectorXf b = X.householderQr().solve(y);
    value = b[1];
    sigmaSq = expf(b[0]);
    snp2pqPowS = (b[1]*logSnp2pq).exp();
}

void BayesS::Sp::sampleFromFC2(const unsigned numNonZeros, const VectorXf &snpEffects,
                               const VectorXf &snp2pq, ArrayXf &snp2pqPowS, const ArrayXf &logSnp2pq,
                               const float varg, float &sum2pqSplusOne) {
    // Hamiltonian Monte Carlo
    // This is for the robust parameterisation where sigmaSq is determined by heritability and sum of 2pq to the power of (S+1)
    
    // Cautious:
    // The sampled value of SNP effect can be exactly zero even it is in the model. In this case, the numNonZeros will be inflated and cause zero element at the end of snp2pqDelta1 vector.
    // To get around this, recalculate numNonZeros here.
    
    unsigned nnz = 0;
    for (unsigned i=0; i<numSnps; ++i)
        if (snpEffects[i]) ++nnz;
    
    // Prepare
    ArrayXf snpEffectDelta1(nnz);
    ArrayXf snp2pqDelta1(nnz);
    ArrayXf logSnp2pqDelta1(nnz);
    
    for (unsigned i=0, j=0; i<numSnps; ++i) {
        if (snpEffects[i]) {
            snpEffectDelta1[j] = snpEffects[i];
            snp2pqDelta1[j] = snp2pq[i];
            logSnp2pqDelta1[j] = logSnp2pq[i];
            ++j;
        }
    }
    
    float snp2pqLogSumDelta1 = logSnp2pqDelta1.sum();
    
    float curr = value;
    float curr_p = Stat::snorm();
    
    float cand = curr;
    // Make a half step for momentum at the beginning
    float cand_p = curr_p - 0.5*stepSize * gradientU2(curr, snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, varg);
    
    for (unsigned i=0; i<numSteps; ++i) {
        // Make a full step for the position
        cand += stepSize * cand_p;
        if (i < numSteps-1) {
            // Make a full step for the momentum, except at end of trajectory
            cand_p -= stepSize * gradientU2(cand, snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, varg);
        } else {
            // Make a half step for momentum at the end
            cand_p -= 0.5*stepSize * gradientU2(cand, snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, varg);
        }
        //cout << i << " " << cand << endl;
    }

    // Evaluate potential (negative log posterior) and kinetic energies at start and end of trajectory
    float scaleCurr, scaleCand;
    float curr_U_chisq, cand_U_chisq;
    float curr_H = computeU2(curr, snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, varg) + 0.5*curr_p*curr_p;
    float cand_H = computeU2(cand, snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, varg) + 0.5*cand_p*cand_p;
    
    if (Stat::ranf() < exp(curr_H-cand_H)) {  // accept
        value = cand;
        snp2pqPowS = snp2pq.array().pow(cand);
        sum2pqSplusOne = snp2pqDelta1.pow(1.0+value).sum();
        ar.count(1, 0.5, 0.9);
    } else {
        ar.count(0, 0.5, 0.9);
    }
    
    if (!(ar.cnt % 10)) {
        if      (ar.value < 0.6) stepSize *= 0.8;
        else if (ar.value > 0.8) stepSize *= 1.2;
    }

    if (ar.consecRej > 20) stepSize *= 0.8;

    tuner.value = stepSize;
}

float BayesS::Sp::gradientU2(const float S, const ArrayXf &snpEffects, const float snp2pqLogSum, const ArrayXf &snp2pq, const ArrayXf &logSnp2pq, const float varg){
    // compute the first derivative of the negative log posterior
    long size = snp2pq.size();
    long chunkSize = size/omp_get_max_threads();
    ArrayXf snp2pqPowS(size);
#pragma omp parallel for schedule(dynamic, chunkSize)
    for (unsigned i=0; i<size; ++i) {
        snp2pqPowS[i] = powf(snp2pq[i], S);
    }
    ArrayXf snp2pqPowSplusOne = snp2pqPowS*snp2pq;
    ArrayXf snpEffectSqOverSnp2pqPowS = snpEffects.square()/snp2pqPowS;
    float constantA = snp2pqLogSum;
    float constantB = snp2pqPowSplusOne.sum();
    float constantC = (logSnp2pq*snp2pqPowSplusOne).sum();
    float constantD = snpEffectSqOverSnp2pqPowS.sum();
    float constantE = (snpEffectSqOverSnp2pqPowS*logSnp2pq).sum();
    float ret = 0.5*constantA - 0.5*float(size)*constantC/constantB + 0.5*constantC*constantD/varg - 0.5*constantB*constantE/varg + S/var;
    return ret;
}

float BayesS::Sp::computeU2(const float S, const ArrayXf &snpEffects, const float snp2pqLogSum, const ArrayXf &snp2pq, const ArrayXf &logSnp2pq, const float varg){
    // compute negative log posterior
    ArrayXf snp2pqPowS = snp2pq.pow(S);
    ArrayXf snp2pqPowSplusOne = snp2pqPowS*snp2pq;
    float m = snp2pq.size();
    float constantA = snp2pqLogSum;
    float constantB = snp2pqPowSplusOne.sum();
    float constantC = (snpEffects.square()/snp2pqPowS).sum();
    // cout << "Hello I'm in computeU" << endl;
    float ret = 0.5*S*constantA - 0.5*m*log(constantB) + 0.5*constantB*constantC/varg + 0.5*S*S/var;
    return ret;
}



void BayesS::SnpEffects::sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes,
                                      const float sigmaSq, const float pi, const float vare,
                                      const ArrayXf &snp2pqPowS, const VectorXf &snp2pq,
                                      const float vg, float &scale, VectorXf &ghat){
    wtdSumSq = 0.0;
    numNonZeros = 0;
    
    ghat.setZero(ycorr.size());

    float oldSample;
    float rhs, invLhs, uhat;
    float logDelta0, logDelta1, probDelta1;
    float logPi = log(pi);
    float logPiComp = log(1.0-pi);
    float invVare = 1.0f/vare;
    float invSigmaSq = 1.0f/sigmaSq;
    
    for (unsigned i=0; i<size; ++i) {
        if (!ZPZdiag[i]) continue;
        
        oldSample = values[i];
        rhs = Z.col(i).dot(ycorr);
        rhs += ZPZdiag[i]*oldSample;
        rhs *= invVare;
        invLhs = 1.0f/(ZPZdiag[i]*invVare + invSigmaSq/snp2pqPowS[i]);
        uhat = invLhs*rhs;
        logDelta1 = 0.5*(logf(invLhs) - logf(snp2pqPowS[i]*sigmaSq) + uhat*rhs) + logPi;
        logDelta0 = logPiComp;
        
        probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
        
        if (bernoulli.sample(probDelta1)) {
            values[i] = normal.sample(uhat, invLhs);
            ycorr += Z.col(i) * (oldSample - values[i]);
            if (weightedRes) ghat += Z.col(i).cwiseProduct(Rsqrt) * values[i];
            else ghat  += Z.col(i) * values[i];
            wtdSumSq += values[i]*values[i]/snp2pqPowS[i];
            ++numNonZeros;
        } else {
            if (oldSample) ycorr += Z.col(i) * oldSample;
            values[i] = 0.0;
        }
    }
}

void BayesS::sampleUnknowns(){
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }

    unsigned cnt=0;
    do {
        snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, data.Rsqrt, data.weightedRes, sigmaSq.value, pi.value, vare.value, snp2pqPowS, data.snp2pq, genVarPrior, sigmaSq.scale, ghat);
        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros == 0);
    
    sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
    
    if (estimatePi) pi.sampleFromFC(snpEffects.size, snpEffects.numNonZeros);
    vare.sampleFromFC(ycorr);
    
//    float sse = ycorr.squaredNorm();
//    float ypy = data.y.squaredNorm();
//    float bpr = snpEffects.values.dot(data.Z.transpose()*data.y);
//    cout << "sse " << sse << " ypy " << ypy << " b'r " << bpr << " b'X'Xb " << sse - ypy + 2*bpr << endl;
    

    S.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros, sigmaSq.value, snpEffects.values, data.snp2pq, snp2pqPowS, logSnp2pq, genVarPrior, sigmaSq.scale, snpEffects.sum2pqSplusOne);
    if (iter >= 2000) sigmaSq.scale = scalePrior;
    scale.getValue(sigmaSq.scale);
    
    varg.compute(ghat);
    hsq.compute(varg.value, vare.value);
    
    rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);
    nnzSnp.getValue(snpEffects.numNonZeros);
    
    if (++iter < 2000) {
        genVarPrior += (varg.value - genVarPrior)/iter;
        scalePrior  += (sigmaSq.scale - scalePrior)/iter;
    }
}

void BayesS::sampleStartVal(){
    sigmaSq.sampleFromPrior();
    if (estimatePi) pi.sampleFromPrior();
    S.sampleFromPrior();
    cout << "  Starting value for " << sigmaSq.label << ": " << sigmaSq.value << endl;
    if (estimatePi) cout << "  Starting value for " << pi.label << ": " << pi.value << endl;
    cout << "  Starting value for " << S.label << ": " << S.value << endl;
    cout << endl;
}


void BayesS::findStartValueForS(const vector<float> &val){
    long size = val.size();
    float start;
    if (size == 1) start = val[0];
    else {
        cout << "Finding the optimal starting value for S ..." << endl;
        float loglike=0, topLoglike=0, optimal=0;
        unsigned idx = 0;
        for (unsigned i=0; i<size; ++i) {
            vector<float> cand = {val[i]};
            BayesS *model = new BayesS(data, varg.value, vare.value, sigmaSqRand.value, pi.value, pi.alpha, pi.beta, estimatePi, S.var, cand, "", false);
            unsigned numiter = 100;
            for (unsigned iter=0; iter<numiter; ++iter) {
                model->sampleUnknownsWarmup();
            }
            loglike = model->computeLogLikelihood();
            if (i==0) {
                topLoglike = loglike;
                optimal = model->S.value;
            }
            if (loglike > topLoglike) {
                idx = i;
                topLoglike = loglike;
                optimal = model->S.value;
            }
            //cout << val[i] <<" " << loglike << " " << model->S.value << endl;
            delete model;
        }
        start = optimal;
        cout << "The optimal starting value for S is " << start << endl;
    }
    S.value = start;
}

float BayesS::computeLogLikelihood(){
    float sse = ycorr.squaredNorm();
    return -0.5f*data.numKeptInds*log(vare.value) - 0.5f*sse/vare.value;
}

void BayesS::sampleUnknownsWarmup(){
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, data.Rsqrt, data.weightedRes, sigmaSq.value, pi.value, vare.value, snp2pqPowS, data.snp2pq, varg.value, sigmaSq.scale, ghat);
    sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
    if (estimatePi) pi.sampleFromFC(snpEffects.size, snpEffects.numNonZeros);
    vare.sampleFromFC(ycorr);
    S.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros, sigmaSq.value, snpEffects.values, data.snp2pq, snp2pqPowS, logSnp2pq, varg.value, sigmaSq.scale, snpEffects.sum2pqSplusOne);
    scale.getValue(sigmaSq.scale);
    varg.compute(ghat);
}


void BayesNS::SnpEffects::sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes,
                                       const float sigmaSq, const float pi, const float vare,
                                       const ArrayXf &snp2pqPowS, const VectorXf &snp2pq,
                                       const float vg, float &scale, VectorXf &ghat){
    wtdSumSq = 0.0;
    numNonZeros = 0;
    numNonZeroWind = 0;
    
    ghat.setZero(ycorr.size());
    
    float oldSample;
    float rhs, invLhs, uhat;
    float logDelta0, logDelta1, probDelta1;
    float logPi = log(pi);
    float logPiComp = log(1.0-pi);
    float invVare = 1.0f/vare;
    float invSigmaSq = 1.0f/sigmaSq;
    float diffQuadSum;
    float logDelta0MinusLogDelta1;
    float snp2pqOneMinusS;
    
    unsigned start, end;
    
    for (unsigned i=0; i<numWindows; ++i) {
        start = windStart[i];
        end = i+1 < numWindows ? windStart[i] + windSize[i] : size;
        
        // sample window delta
        diffQuadSum = 0.0;
        snp2pqOneMinusS = 0.0;
        if (windDelta[i]) {
            for (unsigned j=start; j<end; ++j) {
                if (snpDelta[j]) {
                    rhs = Z.col(j).dot(ycorr);
                    diffQuadSum += 2.0f*beta[j]*rhs + beta[j]*beta[j]*ZPZdiag[j];
                }
            }
        } else {
            for (unsigned j=start; j<end; ++j) {
                if (snpDelta[j]) {
                    rhs = Z.col(j).dot(ycorr);
                    diffQuadSum += 2.0f*beta[j]*rhs - beta[j]*beta[j]*ZPZdiag[j];
                }
            }
        }
        
        diffQuadSum *= invVare;
        logDelta0MinusLogDelta1 = -0.5f*diffQuadSum + logPiComp - logPi;
        probDelta1 = 1.0f/(1.0f + expf(logDelta0MinusLogDelta1));
        
        if (bernoulli.sample(probDelta1)) {
            if (!windDelta[i]) {
                for (unsigned j=start; j<end; ++j) {
                    if (snpDelta[j]) {
                        ycorr -= Z.col(j) * beta[j];
                    }
                }
            }
            windDelta[i] = 1.0;
            ++numNonZeroWind;

            for (unsigned j=start; j<end; ++j) {
                oldSample = beta[j]*snpDelta[j];
                rhs = Z.col(j).dot(ycorr);
                rhs += ZPZdiag[j]*oldSample;
                rhs *= invVare;
                invLhs = 1.0f/(ZPZdiag[j]*invVare + invSigmaSq/snp2pqPowS[j]);
                uhat = invLhs*rhs;
                logDelta1 = 0.5*(logf(invLhs) - logf(snp2pqPowS[j]*sigmaSq) + uhat*rhs) + logLocalPi[i];
                logDelta0 = logLocalPiComp[i];
                
                probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
                
                if (bernoulli.sample(probDelta1)) {
                    values[j] = beta[j] = normal.sample(uhat, invLhs);
                    ycorr += Z.col(j) * (oldSample - values[j]);
                    if (weightedRes) ghat += Z.col(j).cwiseProduct(Rsqrt) * values[j];
                    else ghat  += Z.col(j) * values[j];
                    wtdSumSq += values[j]*values[j]/snp2pqPowS[j];
                    snpDelta[j] = 1.0;
                    ++cumDelta[j];
                    ++numNonZeros;
                } else {
                    if (oldSample) ycorr += Z.col(j) * oldSample;
                    beta[j] = normal.sample(0.0, snp2pqPowS[j]*sigmaSq);
                    snpDelta[j] = 0.0;
                    values[j] = 0.0;
                }
            }
        }
        else {
            for (unsigned j=start; j<end; ++j) {
                beta[j] = normal.sample(0.0, snp2pqPowS[j]*sigmaSq);
                snpDelta[j] = bernoulli.sample(localPi[i]);
                if (values[j]) ycorr += Z.col(j) * values[j];
                values[j] = 0.0;
            }
            windDelta[i] = 0.0;
        }
    }
}

void BayesNS::Sp::sampleFromFC(const unsigned numNonZeros, const float sigmaSq, const VectorXf &snpEffects, const VectorXf &snp2pq, ArrayXf &snp2pqPowS, const ArrayXf &logSnp2pq){
    // do not update scale factor of sigmaSq
    
    // Prepare
    ArrayXf snpEffectDelta1(numNonZeros);
    ArrayXf snp2pqDelta1(numNonZeros);
    ArrayXf logSnp2pqDelta1(numNonZeros);
    
    for (unsigned i=0, j=0; i<numSnps; ++i) {
        if (snpEffects[i]) {
            snpEffectDelta1[j] = snpEffects[i];
            snp2pqDelta1[j] = snp2pq[i];
            logSnp2pqDelta1[j] = logSnp2pq[i];
            ++j;
        }
    }
    
    float snp2pqLogSumDelta1 = logSnp2pqDelta1.sum();
    
    float curr = value;
    float curr_p = Stat::snorm();
    
    float cand = curr;
    // Make a half step for momentum at the beginning
    float cand_p = curr_p - 0.5*stepSize * gradientU(curr,  snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, sigmaSq);
    
    for (unsigned i=0; i<numSteps; ++i) {
        // Make a full step for the position
        cand += stepSize * cand_p;
        if (i < numSteps-1) {
            // Make a full step for the momentum, except at end of trajectory
            cand_p -= stepSize * gradientU(cand, snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, sigmaSq);
        } else {
            // Make a half step for momentum at the end
            cand_p -= 0.5*stepSize * gradientU(cand, snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, sigmaSq);
        }
    }

    // Evaluate potential (negative log posterior) and kinetic energies at start and end of trajectory
    float curr_H = computeU(curr, snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, sigmaSq) + 0.5*curr_p*curr_p;
    float cand_H = computeU(cand, snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, sigmaSq) + 0.5*cand_p*cand_p;
    
    if (Stat::ranf() < exp(curr_H-cand_H)) {  // accept
        value = cand;
        snp2pqPowS = snp2pq.array().pow(cand);
        ar.count(1, 0.5, 0.9);
    } else {
        ar.count(0, 0.5, 0.9);
    }
    
    if      (ar.value < 0.5) stepSize *= 0.8;
    else if (ar.value > 0.9) stepSize *= 1.2;
    
    tuner.value = stepSize;
}

float BayesNS::Sp::gradientU(const float S, const ArrayXf &snpEffects, const float snp2pqLogSum,
                             const ArrayXf &snp2pq, const ArrayXf &logSnp2pq, const float sigmaSq){
    // compute the first derivative of the negative log posterior
    return 0.5*snp2pqLogSum - 0.5/sigmaSq*(snpEffects.square()*logSnp2pq/snp2pq.pow(S)).sum() + S/var;
}

float BayesNS::Sp::computeU(const float S, const ArrayXf &snpEffects, const float snp2pqLogSum,
                            const ArrayXf &snp2pq, const ArrayXf &logSnp2pq, const float sigmaSq){
    // compute negative log posterior
    return 0.5*S*snp2pqLogSum + 0.5/sigmaSq*(snpEffects.square()/snp2pq.pow(S)).sum() + 0.5*S*S/var;
}

void BayesNS::sampleUnknowns(){
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }

    unsigned cnt=0;
    do {
        snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, data.Rsqrt, data.weightedRes, sigmaSq.value, pi.value, vare.value, snp2pqPowS, data.snp2pq, genVarPrior, sigmaSq.scale, ghat);
        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros == 0);
    
    sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
    
    //scale.sampleFromFC(sigmaSq.value, sigmaSq.df, sigmaSq.scale);
    
    if (estimatePi) pi.sampleFromFC(snpEffects.numWindows, snpEffects.numNonZeroWind);
    vare.sampleFromFC(ycorr);
    
    S.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros, sigmaSq.value, snpEffects.values, data.snp2pq, snp2pqPowS, logSnp2pq, genVarPrior, sigmaSq.scale, snpEffects.sum2pqSplusOne);
    if (iter >= 2000) sigmaSq.scale = scalePrior;
    scale.getValue(sigmaSq.scale);
    
    varg.compute(ghat);
    hsq.compute(varg.value, vare.value);
    
    rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);
    nnzSnp.getValue(snpEffects.numNonZeros);
    nnzWind.getValue(snpEffects.numNonZeroWind);
    windDelta.getValues(snpEffects.windDelta);
    
    if (++iter < 2000) {
        genVarPrior += (varg.value - genVarPrior)/iter;
        scalePrior += (sigmaSq.scale - scalePrior)/iter;
    }
}


void BayesRS::SnpEffects::sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes, const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare, const ArrayXf &snp2pqPowS, const VectorXf &snp2pq, const float varg, float &scale, VectorXf &ghat, const bool originalModel) {
    
    wtdSumSq = 0.0;
    numNonZeros = 0.0;

    ghat.setZero(ycorr.size());
    
    ArrayXf wtdSigmaSq(ndist);
    ArrayXf invWtdSigmaSq(ndist);
    ArrayXf logWtdSigmaSq(ndist);
    ArrayXf logPis = pis.array().log();
    ArrayXf log2pqPowS = snp2pqPowS.log();
    
    if (originalModel) {
        wtdSigmaSq = gamma * 0.01 * varg;
    } else {
        wtdSigmaSq = gamma * sigmaSq;
    }
    
    invWtdSigmaSq = wtdSigmaSq.inverse();
    logWtdSigmaSq = wtdSigmaSq.log();
    
    numSnpMix.setZero(ndist);
    snpset.resize(ndist);
    
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
    }
    
    float oldSample;
    float rhs;
    float invVare = 1.0f/vare;

    ArrayXf invLhs(ndist);
    ArrayXf uhat(ndist);
    ArrayXf logDelta(ndist);
    ArrayXf probDelta(ndist);
    
    unsigned delta;

    for (unsigned i=0; i<size; ++i) {
        
        oldSample = values[i];
        rhs = Z.col(i).dot(ycorr);
        rhs += ZPZdiag[i] * oldSample;
        rhs *= invVare;
        
        invLhs = (ZPZdiag[i]*invVare + invWtdSigmaSq/snp2pqPowS[i]).inverse();
        uhat = invLhs*rhs;
        
        logDelta = 0.5*(invLhs.log() - log2pqPowS[i] - logWtdSigmaSq + uhat*rhs) + logPis;
        logDelta[0] = logPis[0];
        
        for (unsigned k=0; k<ndist; ++k) {
            probDelta[k] = 1.0f/(logDelta-logDelta[k]).exp().sum();
        }
        
        delta = bernoulli.sample(probDelta);
        
        snpset[delta].push_back(i);
        numSnpMix[delta]++;
        
        if (delta) {
            values[i] = normal.sample(uhat[delta], invLhs[delta]);
            ycorr += Z.col(i) * (oldSample - values[i]);
            if (weightedRes) ghat += Z.col(i).cwiseProduct(Rsqrt) * values[i];
            else ghat  += Z.col(i) * values[i];
            wtdSumSq += (values[i] * values[i]) / (gamma[delta]*snp2pqPowS[i]);
            ++numNonZeros;
        }
        else {
            if (oldSample) ycorr += Z.col(i) * oldSample;
            values[i] = 0.0;
        }
    }
}

void BayesRS::Sp::sampleFromFC(vector<vector<unsigned> > &snpset, const VectorXf &snpEffects,
                                     float &sigmaSq, const VectorXf &gamma,
                                     const VectorXf &snp2pq, ArrayXf &snp2pqPowS, const ArrayXf &logSnp2pq,
                                     const float vg, float &scale, float &sum2pqSplusOne) {
    // Hamiltonian Monte Carlo
    // note that the scale factor of sigmaSq will be simultaneously updated
    
    unsigned nnzMix = snpset.size() - 1; // nonzero component
    
    // Prepare
    vector<ArrayXf> snpEffectMix(nnzMix);
    vector<ArrayXf> snp2pqMix(nnzMix);
    vector<ArrayXf> logSnp2pqMix(nnzMix);
    
    float snp2pqLogSumNZ = 0.0;
    
    for (unsigned i=0; i<nnzMix; ++i) {
        unsigned k=i+1;
        long isize = snpset[k].size();
        snpEffectMix[i].resize(isize);
        snp2pqMix[i].resize(isize);
        logSnp2pqMix[i].resize(isize);
        for (unsigned j=0; j<isize; ++j) {
            snpEffectMix[i][j] = snpEffects[snpset[k][j]];
            snp2pqMix[i][j] = snp2pq[snpset[k][j]];
            logSnp2pqMix[i][j] = logSnp2pq[snpset[k][j]];
        }
        snp2pqLogSumNZ += logSnp2pqMix[i].sum();
    }
    
    float curr = value;
    float curr_p = Stat::snorm();
    
    float cand = curr;
    // Make a half step for momentum at the beginning
    float cand_p = curr_p - 0.5*stepSize * gradientU(curr, nnzMix, snpEffectMix, snp2pqLogSumNZ, snp2pqMix, logSnp2pqMix, sigmaSq, gamma, vg);

    for (unsigned i=0; i<numSteps; ++i) {
        // Make a full step for the position
        cand += stepSize * cand_p;
        if (i < numSteps-1) {
            // Make a full step for the momentum, except at end of trajectory
            cand_p -= stepSize * gradientU(cand, nnzMix, snpEffectMix, snp2pqLogSumNZ, snp2pqMix, logSnp2pqMix, sigmaSq, gamma, vg);
        } else {
            // Make a half step for momentum at the end
            cand_p -= 0.5*stepSize * gradientU(cand, nnzMix, snpEffectMix, snp2pqLogSumNZ, snp2pqMix, logSnp2pqMix, sigmaSq, gamma, vg);
        }
        //cout << i << " " << cand << endl;
    }

    // Evaluate potential (negative log posterior) and kinetic energies at start and end of trajectory
    float scaleCurr, scaleCand;
    float curr_H = computeU(curr, nnzMix, snpEffectMix, snp2pqLogSumNZ, snp2pqMix, logSnp2pqMix, sigmaSq, gamma, vg, scaleCurr) + 0.5*curr_p*curr_p;
    float cand_H = computeU(cand, nnzMix, snpEffectMix, snp2pqLogSumNZ, snp2pqMix, logSnp2pqMix, sigmaSq, gamma, vg, scaleCand) + 0.5*cand_p*cand_p;
    
    if (Stat::ranf() < exp(curr_H-cand_H)) {  // accept
        value = cand;
        scale = scaleCand;
        snp2pqPowS = snp2pq.array().pow(cand);
        sum2pqSplusOne = 0.0;
        for (unsigned i=0; i<nnzMix; ++i) sum2pqSplusOne += snp2pqMix[i].pow(1.0+value).sum();
        ar.count(1, 0.5, 0.9);
    } else {
        ar.count(0, 0.5, 0.9);
    }
    
    if (!(ar.cnt % 10)) {
        if      (ar.value < 0.6) stepSize *= 0.8;
        else if (ar.value > 0.8) stepSize *= 1.2;
    }
    
    if (ar.consecRej > 20) stepSize *= 0.8;
    
    tuner.value = stepSize;
}

float BayesRS::Sp::gradientU(const float S, const unsigned nnzMix, const vector<ArrayXf> &snpEffectMix, const float snp2pqLogSum, const vector<ArrayXf> &snp2pqMix, const vector<ArrayXf> &logSnp2pqMix, const float sigmaSq, const VectorXf &gamma, const float vg){
    float constantA = snp2pqLogSum;
    ArrayXf constantB(nnzMix);
    for (unsigned i=0; i<nnzMix; ++i) {
        constantB[i] = (snpEffectMix[i].square()*logSnp2pqMix[i]/snp2pqMix[i].pow(S)).sum()/gamma[i+1];
    }
    return 0.5*constantA - 0.5/sigmaSq*constantB.sum() + S/var;
}

float BayesRS::Sp::computeU(const float S, const unsigned nnzMix, const vector<ArrayXf> &snpEffectMix, const float snp2pqLogSum, const vector<ArrayXf> &snp2pqMix, const vector<ArrayXf> &logSnp2pqMix, const float sigmaSq, const VectorXf &gamma, const float vg, float &scale) {
    vector<ArrayXf> snp2pqPowSMix(nnzMix);
    float constantA = snp2pqLogSum;
    ArrayXf constantB(nnzMix);
    ArrayXf constantC(nnzMix);
    for (unsigned i=0; i<nnzMix; ++i) {
        snp2pqPowSMix[i] = snp2pqMix[i].pow(S);
        constantB[i] = (snpEffectMix[i].square()/snp2pqPowSMix[i]).sum()/gamma[i+1];
        constantC[i] = (snp2pqMix[i]*snp2pqPowSMix[i]).sum();
    }
    scale = 0.5*vg/constantC.sum();
    return 0.5*S*constantA + 0.5/sigmaSq*constantB.sum() + 0.5*S*S/var;
}

void BayesRS::sampleUnknowns() {
    static int iter = 0;
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }

    unsigned cnt=0;
    do {
        snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, data.Rsqrt, data.weightedRes, sigmaSq.value, Pis.values, gamma.values, vare.value, snp2pqPowS, data.snp2pq, varg.value, sigmaSq.scale, ghat, originalModel);
        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros == 0);
    sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
    if (estimatePi) Pis.sampleFromFC(snpEffects.numSnpMix);
    numSnps.getValues(snpEffects.numSnpMix);
    nnzSnp.getValue(snpEffects.numNonZeros);
    vare.sampleFromFC(ycorr);
    varg.compute(ghat);
    hsq.compute(varg.value, vare.value);
    
    S.sampleFromFC(snpEffects.snpset, snpEffects.values, sigmaSq.value, gamma.values, data.snp2pq, snp2pqPowS, logSnp2pq, genVarPrior, sigmaSq.scale, snpEffects.sum2pqSplusOne);

    if (iter >= 2000) sigmaSq.scale = scalePrior;
    scale.getValue(sigmaSq.scale);
    // cout << "iter " << iter << " scalePrior " << scalePrior << "sigmaSq.scale " << sigmaSq.scale << endl;
    
    rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);

    nnzSnp.getValue(snpEffects.numNonZeros);
    
    //    numSnpVg.compute(snpEffects.values, data.ZPZdiag, varg.value, vare.nobs);
    if (originalModel) Vgs.compute(snpEffects.values, data.Z, snpEffects.snpset, varg.value);
    
    if (++iter < 2000) {
        genVarPrior += (varg.value - genVarPrior)/iter;
        scalePrior  += (sigmaSq.scale - scalePrior)/iter;
    }
}



void ApproxBayesC::FixedEffects::sampleFromFC(const MatrixXf &XPX, const VectorXf &XPXdiag,
                                              const MatrixXf &ZPX, const VectorXf &XPy,
                                              const VectorXf &snpEffects, const float vare,
                                              VectorXf &rcorr){
    for (unsigned i=0; i<size; ++i) {
        float oldSample = values[i];
        float XPZa = ZPX.col(i).dot(snpEffects);
        float rhs = XPy[i] - XPZa - XPX.row(i).dot(values) + XPXdiag[i]*values[i];
        float invLhs = 1.0f/XPXdiag[i];
        float bhat = invLhs*rhs;
        values[i] = Normal::sample(bhat, invLhs*vare);
        //rcorr += ZPX.col(i) * (oldSample - values[i]);
    }

}

void ApproxBayesC::SnpEffects::sampleFromFC(VectorXf &rcorr,const vector<SparseVector<float> > &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const VectorXf &snp2pq, const VectorXf &LDsamplVar,
                                            const float sigmaSq, const float pi, const float vare, const float varg, const float ps, const float overdispersion){
    
    long numChr = chromInfoVec.size();

    float ssq[numChr], s2pq[numChr], nnz[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(s2pq,0,sizeof(float)*numChr);
    memset(nnz,0, sizeof(float)*numChr);

//    for (unsigned chr=0; chr<numChr; ++chr) {
//        ChromInfo *chromInfo = chromInfoVec[chr];
//        unsigned chrStart = chromInfo->startSnpIdx;
//        unsigned chrEnd   = chromInfo->endSnpIdx;
//        if (iter==0) {
//            cout << "chr " << chr+1 << " start " << chrStart << " end " << chrEnd << endl;
//        }
//    }
//    if (iter==0) cout << endl; 

    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking precents the writing from other threads

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
#pragma omp parallel for
    for (unsigned chr=0; chr<numChr; ++chr) {
        //cout << " thread " << omp_get_thread_num() << " chr " << chr << endl;
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned windEnd, j;
        
        float oldSample;
        float rhs, invLhs, uhat;
        float logDelta0, logDelta1, probDelta1;
        float logPi = log(pi);
        float logPiComp = log(1.0-pi);
        float logSigmaSq = log(sigmaSq);
        float invSigmaSq = 1.0f/sigmaSq;
        float varei;
        
        for (unsigned i=chrStart; i<=chrEnd; ++i) {

            oldSample = valuesPtr[i];

            if (leaveout[i]) {
                probDelta1 = 0;
            }
            else {

            //cout << i << " " << chrStart << " " << chrEnd << " " << ssq << " " << nnz << endl;
//            if (!(iter % 100)) {
//                //float varei = (sse[i] - values.segment(windStart[i], windSize[i]).dot(ZPy.segment(windStart[i], windSize[i]) + rcorr.segment(windStart[i], windSize[i])))/n[i];
//                windEnd = windStart[i] + windSize[i];
//                varei[i] = tss[i];
//                for (j=windStart[i]; j<windEnd; ++j) {
//                    if (valuesPtr[j]) varei[i] -= valuesPtr[j]*(ZPy[j] + rcorr[j]);
//                }
//                varei[i] /= n[i];
//            }
                
            varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
                
            rhs = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            invLhs = 1.0f/(ZPZdiag[i]/varei + invSigmaSq);
            uhat = invLhs*rhs;
            logDelta1 = 0.5*(logf(invLhs) - logSigmaSq + uhat*rhs) + logPi;
            logDelta0 = logPiComp;
            probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
                
            }
            
            //cout << rhs << " " << invLhs << " " << logDelta1 << " " << logSigmaSq << " " << sigmaSq << endl;
//            if (bernoulli.sample(probDelta1)) {
            if (urnd[i] < probDelta1) {
//                valuesPtr[i] = normal.sample(uhat, invLhs);
                valuesPtr[i] = uhat + nrnd[i]*sqrtf(invLhs);
//                rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*(oldSample - values[i]);
                float sampleDiff = oldSample - valuesPtr[i];
                for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {
                    rcorr[it.index()] += it.value() * sampleDiff;
                }
                ssq[chr]  += valuesPtr[i]*valuesPtr[i];
                s2pq[chr] += snp2pq[i];
                ++nnz[chr];
            } else {
                if (oldSample) {
//                    rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*oldSample;
                    for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {
                        rcorr[it.index()] += it.value() * oldSample;
                    }
                }
                valuesPtr[i] = 0.0;
            }
        }
    }
    
    //cout << ssq << " " << nnz << endl;

    sumSq = 0.0;
    sum2pq = 0.0;
    numNonZeros = 0;
    nnzPerChr.setZero(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        sumSq += ssq[i];
        sum2pq += s2pq[i];
        numNonZeros += nnz[i];
        nnzPerChr[i] = nnz[i];
    }

    values = VectorXf::Map(valuesPtr, size);
}

void ApproxBayesC::SnpEffects::sampleFromFC(VectorXf &rcorr,const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const VectorXf &snp2pq, const VectorXf &LDsamplVar,
                                            const float sigmaSq, const float pi, const float vare, const float varg, const float ps, const float overdispersion){
    
    long numChr = chromInfoVec.size();
    
    float ssq[numChr], nnz[numChr], s2pq[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(nnz,0,sizeof(float)*numChr);
    memset(s2pq,0,sizeof(float)*numChr);

    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking precents the writing from other threads
  
    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
#pragma omp parallel for
    for (unsigned chr=0; chr<numChr; ++chr) {
        //cout << " thread " << omp_get_thread_num() << " chr " << chr << endl;
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned windEnd, j;
        
        float oldSample;
        float rhs, invLhs, uhat;
        float logDelta0, logDelta1, probDelta1;
        float logPi = log(pi);
        float logPiComp = log(1.0-pi);
        float logSigmaSq = log(sigmaSq);
        float invSigmaSq = 1.0f/sigmaSq;
        float varei;
        
        for (unsigned i=chrStart; i<=chrEnd; ++i) {

            oldSample = valuesPtr[i];

            if (leaveout[i]) {
                probDelta1 = 0;
            }
            else {
            
            //cout << i << " " << chrStart << " " << chrEnd << " " << ssq << " " << nnz << endl;
            // cout << "Varei i " << varei[i] << endl;
//            if (!(iter % 100)) {
//                //float varei = (sse[i] - values.segment(windStart[i], windSize[i]).dot(ZPy.segment(windStart[i], windSize[i]) + rcorr.segment(windStart[i], windSize[i])))/n[i];
//                windEnd = windStart[i] + windSize[i];
//                varei[i] = tss[i];
//                for (j=windStart[i]; j<windEnd; ++j) {
//                    if (valuesPtr[j]) varei[i] -= valuesPtr[j]*(ZPy[j] + rcorr[j]);
//                }
//                varei[i] /= n[i];
//                // cout << "Varei 100 " << varei[i] << endl;
//            }
            // varei[i] = tss[i] / n[i];
            //            varei = se[i]*se[i]*ZPZdiag[i];
            
            varei = LDsamplVar[i]*varg + vare + ps + overdispersion;

            rhs = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            invLhs = 1.0f/(ZPZdiag[i]/varei + invSigmaSq);
            uhat = invLhs*rhs;
            logDelta1 = 0.5*(logf(invLhs) - logSigmaSq + uhat*rhs) + logPi;
            logDelta0 = logPiComp;
            probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
            //cout << rhs << " " << invLhs << " " << logDelta1 << " " << logSigmaSq << " " << sigmaSq << endl;
                
            }
            
//            if (bernoulli.sample(probDelta1)) {
            if (urnd[i] < probDelta1) {
//                valuesPtr[i] = normal.sample(uhat, invLhs);
                valuesPtr[i] = uhat + nrnd[i]*sqrtf(invLhs);
                rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*(oldSample - valuesPtr[i]);
                ssq[chr] += valuesPtr[i]*valuesPtr[i];
                s2pq[chr] += snp2pq[i];
                ++nnz[chr];
            } else {
                if (oldSample) rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*oldSample;
                valuesPtr[i] = 0.0;
            }
        }
    }
    // cout << "Varei 1 max" << varei.maxCoeff() << endl;
    //cout << ssq << " " << nnz << endl;
    
    sumSq = 0.0;
    sum2pq = 0.0;
    numNonZeros = 0.0;
    nnzPerChr.setZero(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        sumSq += ssq[i];
        sum2pq += s2pq[i];
        numNonZeros += nnz[i];
        nnzPerChr[i] = nnz[i];
    }
    
    values = VectorXf::Map(valuesPtr, size);
}

void ApproxBayesC::SnpEffects::hmcSampler(VectorXf &rcorr, const VectorXf &ZPy, const vector<VectorXf> &ZPZ,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const float sigmaSq, const float pi, const float vare){
    
    float stepSize = 0.001;
    unsigned numSteps = 1;
    
    
    //#pragma omp parallel for   // this multi-thread may not work due to vector locking when write to the vector
    for (unsigned chr=0; chr<chromInfoVec.size(); ++chr) {
        //cout << " thread " << omp_get_thread_num() << " chr " << chr << endl;
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned chrSize  = chromInfo->size;
        
        VectorXf chrZPy = ZPy.segment(chrStart, chrSize);
        VectorXi chrWindStart = windStart.segment(chrStart, chrSize);
        VectorXi chrWindSize = windSize.segment(chrStart, chrSize);
        chrWindStart.array() -= chrStart;
        

        VectorXf delta;
        delta.setZero(chrSize);
        for (unsigned i=chrStart, j=0; i<=chrEnd; ++i) {
            if (values[i]) {
                delta[j++] = 1;
            }
        }
        
        
        VectorXf curr = values.segment(chrStart, chrSize);
        VectorXf curr_p(chrSize);
        
        for (unsigned i=0; i<chrSize; ++i) {
            curr_p[i] = Stat::snorm();
        }
        
        VectorXf cand = curr.cwiseProduct(delta);
        // Make a half step for momentum at the beginning
        VectorXf rc = chrZPy;
        VectorXf cand_p = curr_p.cwiseProduct(delta) - 0.5*stepSize * gradientU(curr, rc, chrZPy, ZPZ, chrWindStart, chrWindSize, chrStart, chrSize, sigmaSq, vare).cwiseProduct(delta);
        
        for (unsigned i=0; i<numSteps; ++i) {
            cand += stepSize * cand_p.cwiseProduct(delta);
            if (i < numSteps-1) {
                cand_p -= stepSize * gradientU(cand, rc, chrZPy, ZPZ, chrWindStart, chrWindSize, chrStart, chrSize, sigmaSq, vare).cwiseProduct(delta);
            } else {
                cand_p -= 0.5* stepSize * gradientU(cand, rc, chrZPy, ZPZ, chrWindStart, chrWindSize, chrStart, chrSize, sigmaSq, vare).cwiseProduct(delta);
            }
        }
        
        float curr_H = computeU(curr, rcorr.segment(chrStart, chrSize), chrZPy, sigmaSq, vare) + 0.5*curr_p.squaredNorm();
        float cand_H = computeU(cand, rc, chrZPy, sigmaSq, vare) + 0.5*cand_p.squaredNorm();
        
        if (Stat::ranf() < exp(curr_H-cand_H)) {  // accept
            values.segment(chrStart, chrSize) = cand;
            rcorr.segment(chrStart, chrSize) = rc;
            ++mhr;
        }
    }
    
    sumSq = values.squaredNorm();
    //numNonZeros = size;
    
    for (unsigned i=0; i<size; ++i) {
        if(values[i]) ++numNonZeros;
    }
    //cout << sumSq << " " << nnz << " " << numNonZeros << endl;
    
    //cout << values.head(10).transpose() << endl;
    
//    if (!(++cnt % 100) && myMPI::rank==0) {
//        float ar = mhr/float(cnt*22);
//        if      (ar < 0.5) cout << "Warning: acceptance rate for SNP effects is too low "  << ar << endl;
//        else if (ar > 0.9) cout << "Warning: acceptance rate for SNP effects is too high " << ar << endl;
//    }

}

VectorXf ApproxBayesC::SnpEffects::gradientU(const VectorXf &effects, VectorXf &rcorr, const VectorXf &ZPy, const vector<VectorXf> &ZPZ, 
                                             const VectorXi &windStart, const VectorXi &windSize, const unsigned chrStart, const unsigned chrSize,
                                             const float sigmaSq, const float vare){
    rcorr = ZPy;
    for (unsigned i=0; i<chrSize; ++i) {
        if (effects[i]) {
            rcorr.segment(windStart[i], windSize[i]) -= ZPZ[chrStart+i]*effects[i];
        }
    }
    return -rcorr/vare + effects/sigmaSq;
}

float ApproxBayesC::SnpEffects::computeU(const VectorXf &effects, const VectorXf &rcorr, const VectorXf &ZPy,                                             const float sigmaSq, const float vare){
    return -0.5f/vare*effects.dot(ZPy+rcorr) + 0.5/sigmaSq*effects.squaredNorm();
}

//void ApproxBayesC::ResidualVar::sampleFromFC(VectorXf &rcorr, const SpMat &ZPZinv){ // this would not work if ZPZ has no full rank, e.g. exist SNPs in complete LD.
//    float sse = rcorr.transpose()*ZPZinv*rcorr;
//    float dfTilde = df + nobs;
//    float scaleTilde = sse + df*scale;
//    value = InvChiSq::sample(dfTilde, scaleTilde);
//}

//void ApproxBayesC::ResidualVar::sampleFromFC(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr, const float varg, const float nnz){
//    float sse = ypy - effects.dot(ZPy) - effects.dot(rcorr) + nobs*varg*nnz*icrsq;
////    if (sse < 0) sse = 0.0;
////    if (sse > ypy) sse = ypy;
//    float dfTilde = df + nobs;
//    float scaleTilde = sse + df*scale;
//    value = InvChiSq::sample(dfTilde, scaleTilde);
//}

void ApproxBayesC::ResidualVar::sampleFromFC(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr, const float covg) {
    float sse = ypy - effects.dot(ZPy) - effects.dot(rcorr) + nobs*covg;
    if (sse < 0) {
        string vare_str = to_string(static_cast<long double>(sse/nobs));
        throw("\nError: Residual variance is negative (" + vare_str + "). This may indicate that effect sizes are \"blowing up\" likely due to a convergence problem. If SigmaSq variable is increasing with MCMC iterations, then this further indicates MCMC may not converge.");
    }
//    if (sse < 0) sse = 0; //-sse;
//    if (sse > ypy) sse = ypy;
    float dfTilde = df + nobs;
    float scaleTilde = sse + df*scale;
    value = InvChiSq::sample(dfTilde, scaleTilde);
}


//void ApproxBayesC::ResidualVar::sampleFromFCshrink(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr, const float hsq, const float phi){
//    float sse = ypy - effects.dot(ZPy) - effects.dot(rcorr);
////    if (sse < 0) sse = 0.0;
////    if (sse > ypy) sse = ypy;
//    float df = nobs*hsq*phi;   // shrink the residual variance more toward the prior mean when hsq is higher to compensate the bias
//    float dfTilde = df + nobs;
//    float scaleTilde = sse + df*ypy/nobs;    // the prior value of the scale parameter is set to be zero
//    value = InvChiSq::sample(dfTilde, scaleTilde);
//}

//void ApproxBayesC::ResidualVar::sampleFromFC2(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &ghat){
//    float sse = ypy - effects.dot(ZPy)*2.0f + ghat.dot(ghat);
////    cout << "sse " << sse << " ypy " << ypy << " b'r " << effects.dot(ZPy) << " b'X'Xb " << Gadget::calcVariance(ghat)*nobs << endl;
////    if (sse < 0) sse = 0.0;
////    if (sse > ypy) sse = ypy;
//    float dfTilde = df + nobs;
//    float scaleTilde = sse + df*scale;
//    value = InvChiSq::sample(dfTilde, scaleTilde);
//}

//void ApproxBayesC::ResidualVar::randomWalkMHsampler(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr, const VectorXf &ZPZrss, const float sigmaSq, const float pi){
//    // Random walk Mentroplis-Hastings taking into account the variability of sse due to the reduced LD matrix
//    // Assume sse has a gamma distribution
//    
//    float mean = ypy - effects.dot(ZPy) - effects.dot(rcorr);
//    float variance = effects.dot((pi*ZPZrss).cwiseProduct(effects))*sigmaSq;
//    
//    float scale = variance/mean;
//    float shape = mean/scale;
//    
//    float varProp = 0.01*value;
//    
//    float curr = value;
//    float cand = value + Stat::snorm()*sqrtf(varProp);
//
//    float logCurr = (-0.5f*(df+nobs+1)+shape)*logf(curr) - shape*logf(scale+2.0f*curr) - df*scale/(2.0f*curr);
//    float logCand = (-0.5f*(df+nobs+1)+shape)*logf(cand) - shape*logf(scale+2.0f*cand) - df*scale/(2.0f*cand);
//    
//    if (Stat::ranf() < exp(logCand-logCurr)) {  // accept
//        value = cand;
//    }

//    cout << "sse " << mean << " ypy " << ypy << " b'r " << effects.dot(ZPy) << " b'X'Xb " << effects.dot(ZPy) - effects.dot(rcorr) << endl;
//    cout << " mean " << mean << " variance " << variance << " shape " << shape << " scale " << scale << endl;
    
//}

//void ApproxBayesC::Overdispersion::sampleFromFC(const VectorXf &y, const VectorXf &ghat){
//    float sse = (y-ghat).squaredNorm();
//    float dfTilde = df + nobs;
//    float scaleTilde = sse + df*scale;
//    value = InvChiSq::sample(dfTilde, scaleTilde);
//}


//void ApproxBayesC::GenotypicVar::compute(const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr){
//    float modelSS = effects.dot(ZPy) - effects.dot(rcorr);
////    if (modelSS < 0) modelSS = 0;
//    value = modelSS/nobs;
//}

void ApproxBayesC::GenotypicVar::compute(const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr, const float covg){
    float modelSS = effects.dot(ZPy) - effects.dot(rcorr) + nobs*covg;
    if (modelSS < 0) modelSS = 0; // -modelSS;
    value = modelSS/nobs;
}

void ApproxBayesC::Rounding::computeRcorr(const VectorXf &ZPy, const vector<SparseVector<float> > &ZPZ,
                                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                          const VectorXf &snpEffects, VectorXf &rcorr){
    if (count++ % 100) return;
    VectorXf rcorrOld = rcorr;
    rcorr = ZPy;
#pragma omp parallel for
    for (unsigned chr=0; chr<chromInfoVec.size(); ++chr) {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {
                //rcorr[windStart[i]+it.index()] -= it.value() * snpEffects[i];
                rcorr[it.index()] -= it.value() * snpEffects[i];
            }
//            rcorr.segment(windStart[i], windSize[i]) -= ZPZ[i]*snpEffects[i];
        }
    }
    value = sqrt(Gadget::calcVariance(rcorrOld-rcorr));
}

void ApproxBayesC::Rounding::computeRcorr(const VectorXf &ZPy, const vector<VectorXf> &ZPZ,
                                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                          const VectorXf &snpEffects, VectorXf &rcorr){
    if (count++ % 100) return;
    VectorXf rcorrOld = rcorr;
    rcorr = ZPy;
#pragma omp parallel for
    for (unsigned chr=0; chr<chromInfoVec.size(); ++chr) {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            rcorr.segment(windStart[i], windSize[i]) -= ZPZ[i]*snpEffects[i];
        }
    }
    value = sqrt(Gadget::calcVariance(rcorrOld-rcorr));
}

void ApproxBayesC::Rounding::computeGhat(const MatrixXf &Z, const VectorXf &snpEffects, VectorXf &ghat){
    if (count++ % 100) return;
    VectorXf ghatOld = ghat;
    ghat.setZero(ghat.size());
    for (unsigned i=0; i<snpEffects.size(); ++i) {
        if (snpEffects[i]) ghat += Z.col(i)*snpEffects[i];
    }
    value = sqrt(Gadget::calcVariance(ghatOld-ghat));
}

void ApproxBayesC::PopulationStratification::compute(const VectorXf &rcorr, const VectorXf &ZPZdiag, const VectorXf &LDsamplVar, const float varg, const float vare, const VectorXf &chisq){
    
    value = (rcorr.array().square()/(ZPZdiag.array() * (LDsamplVar.array()*varg + value + vare))).mean() - 1.0;
    value = value < -0.01 ? -0.01 : value;
    
//    VectorXf varEta = ZPZdiag.array() * (LDsamplVar.array()*varg + value + vare);
////    VectorXf wt = varEta.array().square().inverse();
//    VectorXf wt = (rcorr.array().square()/ZPZdiag.array().square()).square().inverse();
////    VectorXf zsq = rcorr.array().square()/varEta.array();
//    VectorXf zsq = rcorr.array().square()/ZPZdiag.array() - LDsamplVar.array()*varg - vare;
//    value = zsq.cwiseProduct(wt).sum()/wt.sum();

    
//    VectorXf tmp = rcorr.array().square()/ZPZdiag.array() - LDsamplVar.array()*varg - vare;
//    float ssq = 0.0;
//    long size = rcorr.size();
//    long cnt = 0;
//    for (unsigned i=0; i<size; ++i) {
////        if (chisq[i] < 20 && !(i%20)) {
//            ssq += tmp[i];
//            ++cnt;
////        }
//    }
////    ssq /= float(cnt);
//    if (ssq < 0) ssq = 0.0;
//    float dfTilde = df + cnt;
//    float scaleTilde = ssq + df*scale;
//    value = InvChiSq::sample(dfTilde, scaleTilde);

    
//        ofstream out("rcorr.txt");
//        out << rcorr.array().square()/(ZPZdiag.array() * (LDsamplVar.array()*varg + value + vare)) << endl;
//        out.close();
    
}

void ApproxBayesC::PopulationStratification::compute(const VectorXf &rcorr, const VectorXf &ZPZdiag, const VectorXf &LDsamplVar, const float varg, const float vare, const vector<ChromInfo*> chromInfoVec){
    
    for (unsigned i=0; i<22; ++i) {
        unsigned start = chromInfoVec[i]->startSnpIdx;
        unsigned end = chromInfoVec[i]->endSnpIdx;
        unsigned size = end - start + 1;
        chrSpecific[i] = (rcorr.segment(start,size).array().square()/(ZPZdiag.segment(start,size).array() * (LDsamplVar.segment(start,size).array()*varg + chrSpecific[i] + vare))).mean() - 1.0;
    }    
}

void ApproxBayesC::NumResidualOutlier::compute(const VectorXf &rcorr, const VectorXf &ZPZdiag, const VectorXf &LDsamplVar, const float varg, const float vare, const vector<string> &snpName, VectorXi &leaveout, const vector<SparseVector<float> > &ZPZ, const VectorXf &ZPy, const VectorXf &snpEffects) {
    iter++;
    //if (iter<10) return;
    
    VectorXf tss = ZPy.array().square()/ZPZdiag.array();
    VectorXf sse = rcorr.array().square()/ZPZdiag.array();
    value = 0;
    long size = tss.size();
    for (unsigned i=0; i<size; ++i) {
        if (sse[i] > 10 && tss[i] < 10) ++value;
        //if (sse[i] > 30 && sse[i] > tss[i]) ++value;
    }
    
//    MatrixXf X(tss.size(), 2);
//    X.col(0) = VectorXf::Ones(tss.size());
//    X.col(1) = snpEffects;
//    VectorXf bhat = ZPy.array()/ZPZdiag.array();
//    VectorXf b = X.householderQr().solve(bhat);
//    value = b[1];
    
//    cout << sse.array().maxCoeff() << endl;
    
//    cout << tss.mean() << " " << sse.mean() << endl;
    
//    value = 0;
//    
//    VectorXf tmp = rcorr.array().square()/(ZPZdiag.array() * (LDsamplVar.array()*varg + value + vare));
//    //VectorXf tmp = rcorr.array().square()/ZPZdiag.array() - LDsamplVar.array()*varg;
//    
//    long size = tmp.size();
//    stringstream ss;
//    for (unsigned i=0; i<size; ++i) {
//        if (tmp[i]>20) {
//            ++value;
//            //leaveout[i] = 1;
//            ss << " " << snpName[i];
//            
////            for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {
////                leaveout[it.index()] = 1;
////            }
//            
//        }
//    }
//    if (value) out << iter << ss.str() << endl;
}

void ApproxBayesC::ldScoreReg(const VectorXf &chisq, const VectorXf &LDscore, const VectorXf &LDsamplVar,
                              const float varg, const float vare, float &ps){
    
    long nrow = chisq.size();
    
//    ps = (chisq - LDscore - LDsamplVar*varg - VectorXf::Ones(nrow)*vare).mean();
//    ps = (chisq - LDscore*0.17/float(nrow) - VectorXf::Ones(nrow)).mean();

//    return;
    
    VectorXf y = chisq - LDsamplVar*varg - VectorXf::Ones(nrow)*vare;
//    VectorXf y = chisq;
//    VectorXf weight = 2.0*(LDscore*vargj + LDsamplVar*varg + VectorXf::Ones(nrow)*vare).array().square();
//    VectorXf weight = 2.0*(LDscore*vargj + VectorXf::Ones(nrow)).array().square();
//    VectorXf weightInv = weight.cwiseInverse();
    
    MatrixXf X(nrow, 2);
//    X.col(0) = weight;
//    X.col(1) = LDscore.cwiseProduct(weight);
//    y.array() *= weight.array();

    X.col(0) = VectorXf::Ones(nrow);
    X.col(1) = LDscore;
    
//    unsigned m = 0;
//    for (unsigned i=0; i<nrow; ++i) {
//        if (chisq[i] < 30) ++m;
//    }
//    VectorXf ysub(m);
//    MatrixXf Xsub(m, 2);
//    unsigned j=0;
//    for (unsigned i=0; i<nrow; ++i) {
//        if (chisq[i] < 30) {
//            ysub[j] = y[i];
//            Xsub.row(j) = X.row(i);
//            ++j;
//        }
//    }
//    VectorXf b = Xsub.householderQr().solve(ysub);

    
    VectorXf b = X.householderQr().solve(y);
//    VectorXf b = (X.transpose()*weightInv.asDiagonal()*X).inverse()*X.transpose()*weightInv.asDiagonal()*y;

    ps = b[0];
    
//    cout << b.transpose() << endl;
    
}

void ApproxBayesC::InterChrGenetCov::compute(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr) {
    if (!spouseCorrelation) return;
    float bZPy = effects.dot(ZPy);
    float brcorr = effects.dot(rcorr);
    float varg = (bZPy - brcorr)/nobs;
//    float vare = (ypy - bZPy - brcorr)/nobs;
    float varp = ypy/nobs;
    float hsq = varg/varp;
    float R = spouseCorrelation*hsq / (1 - spouseCorrelation*hsq);
    value = varg * R; // * 0.95;
}

void ApproxBayesC::NnzGwas::compute(const VectorXf &effects, const vector<SparseVector<float> > &ZPZ, const VectorXf &ZPZdiag) {
    if (iter++ % 100) return;
    value = 0;
    long numSnps = effects.size();
    unsigned i, j;
    for (i=0; i<numSnps; ++i) {
        for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {
            j = it.index();
            if (effects[j]){
                if (it.value()*it.value() > 0.1*ZPZdiag[i]*ZPZdiag[j]) {
                    ++value;
                    break;
                }
            }
        }
    }
}

void ApproxBayesC::PiGwas::compute(const float nnzGwas, const unsigned int numSnps) {
    if (iter++ % 100) return;
    value = nnzGwas/float(numSnps);
}

void ApproxBayesC::checkHsq(vector<float> &hsqMCMC) {
    long niter = hsqMCMC.size();
    VectorXf y(niter);
    MatrixXf X(niter, 2);
    for (unsigned i=0; i<niter; ++i) {
        y[i] = hsqMCMC[i];
        X(i,0) = 1;
        X(i,1) = i;
    }
    VectorXf b = X.householderQr().solve(y);
    float slope = b[1];
    float vare = (y.squaredNorm() - (X*b).squaredNorm())/float(niter);
    float se = sqrt(vare/X.col(1).squaredNorm());
    if ((slope/se) > 3) {   // 3 corresponds to P = 0.001 at one-way test
        string slope_str = to_string(static_cast<float>(slope));
        string se_str = to_string(static_cast<float>(se));
        throw("\nError: The SNP-heritability is increasing over MCMC iterations (slope: " + slope_str + "; se: " + se_str + "). This may indicate that effect sizes are \"blowing up\" likely due to a convergence problem.");
    }

}

void ApproxBayesC::sampleUnknowns(){
    static int iter = 0;
//    fixedEffects.sampleFromFC(data.XPX, data.XPXdiag, data.ZPX, data.XPy, snpEffects.values, vare.value, rcorr);
    unsigned cnt=0;
    //do {
        //snpEffects.sampleFromFC(rcorr, data.ZPZ, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, data.se, tss, data.n, data.snp2pq, sigmaSq.value, pi.value, vare.value);
        //snpEffects.hmcSampler(rcorr, data.ZPy, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, sigmaSq.value, pi.value, vare.value);
        if (sparse)
            snpEffects.sampleFromFC(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, data.se, data.tss, varei,
                                    data.n, data.snp2pq, data.LDsamplVar, sigmaSq.value, pi.value, vare.value, varg.value, ps.value, overdispersion);
        else
            snpEffects.sampleFromFC(rcorr, data.ZPZ, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, data.se, data.tss, varei,
                                    data.n, data.snp2pq, data.LDsamplVar, sigmaSq.value, pi.value, vare.value, varg.value, ps.value, overdispersion);
    //    if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    //} while (snpEffects.numNonZeros == 0);
    if (diagnose) nro.compute(rcorr, data.ZPZdiag, data.LDsamplVar, varg.value, vare.value, snpEffects.header, snpEffects.leaveout, data.ZPZsp, data.ZPy, snpEffects.values);
    
    if (robustMode) {
        if (noscale) {
            sigmaSq.value = varg.value/(data.snp2pq.array().sum()*pi.value);
        } else {
            sigmaSq.value = varg.value/(data.numIncdSnps*pi.value);  // LDpred2's parameterisation
        }
    } else {
        sigmaSq.sampleFromFC(snpEffects.sumSq, snpEffects.numNonZeros);
    }

    if (estimatePi) pi.sampleFromFC(data.numIncdSnps, snpEffects.numNonZeros);
    nnzSnp.getValue(snpEffects.numNonZeros);
    sigmaSqG.compute(sigmaSq.value, snpEffects.sum2pq);

    covg.compute(data.ypy, snpEffects.values, data.ZPy, rcorr);
    varg.compute(snpEffects.values, data.ZPy, rcorr, covg.value);
//    varg.value = sigmaSqG.value;
    vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, covg.value);
    hsq.compute(varg.value, vare.value);

    if (iter >= 2000) sigmaSq.scale = scalePrior;
    //sigmaSq.scale = scalePrior;
    scale.getValue(sigmaSq.scale);
        
    if (sparse)
        rounding.computeRcorr(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    else
        rounding.computeRcorr(data.ZPy, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    if (modelPS) ps.compute(rcorr, data.ZPZdiag, data.LDsamplVar, varg.value, vare.value, data.chisq);
//    if (modelPS) ldScoreReg(data.chisq, data.LDscore, data.LDsamplVar, varg.value, vare.value, ps.value);
    
//    if (sparse) {
//        nnzgwas.compute(snpEffects.values, data.ZPZsp, data.ZPZdiag);
//        pigwas.compute(nnzgwas.value, data.numIncdSnps);
//    }
    
    float scaleIteri = 0;
    if (++iter < 2000) {
        if (noscale)
        {
            scaleIteri = 0.5f * varg.value / (data.snp2pq.array().sum()*(pi.value));
        } else
        {
            scaleIteri = 0.5f * varg.value / (data.snp2pq.size()*(pi.value));
        }
        genVarPrior += (varg.value - genVarPrior)/iter;
        scalePrior += (scaleIteri - scalePrior)/iter;
    }
    //if(iter>1990 && iter < 2010){
    //    cout << "iter " << iter << " , scalePrior " << scalePrior << " , sigmaSq.scale " << sigmaSq.scale << " , genVarPrior " << genVarPrior << " , varg " << varg.value << " , pi " << pi.value << endl;
    //}
    
//    if (iter > 100 & !(iter % 10)) {
//        hsqMCMC.push_back(hsq.value);
//        if (!(iter % 1000)) checkHsq(hsqMCMC);
//    }
}


void ApproxBayesB::SnpEffects::sampleFromFC(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const VectorXf &snp2pq, const VectorXf &LDsamplVar,
                                            const VectorXf &sigmaSq, const float pi, const float vare, const float varg, const float ps, const float overdispersion) {
    long numChr = chromInfoVec.size();
    
    float ssq[numChr], s2pq[numChr], nnz[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(s2pq,0,sizeof(float)*numChr);
    memset(nnz,0, sizeof(float)*numChr);
    
    //    for (unsigned chr=0; chr<numChr; ++chr) {
    //        ChromInfo *chromInfo = chromInfoVec[chr];
    //        unsigned chrStart = chromInfo->startSnpIdx;
    //        unsigned chrEnd   = chromInfo->endSnpIdx;
    //        if (iter==0) {
    //            cout << "chr " << chr+1 << " start " << chrStart << " end " << chrEnd << endl;
    //        }
    //    }
    //    if (iter==0) cout << endl;
    
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking precents the writing from other threads
    
    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
#pragma omp parallel for
    for (unsigned chr=0; chr<numChr; ++chr) {
        //cout << " thread " << omp_get_thread_num() << " chr " << chr << endl;
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned windEnd, j;
        
        float oldSample;
        float rhs, invLhs, uhat;
        float logDelta0, logDelta1, probDelta1;
        float logPi = log(pi);
        float logPiComp = log(1.0-pi);
//        float logSigmaSq = log(sigmaSq);
//        float invSigmaSq = 1.0f/sigmaSq;
        float varei;
        
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            
            oldSample = valuesPtr[i];
            
            if (leaveout[i]) {
                probDelta1 = 0;
            }
            else {
                
                //cout << i << " " << chrStart << " " << chrEnd << " " << ssq << " " << nnz << endl;
                //            if (!(iter % 100)) {
                //                //float varei = (sse[i] - values.segment(windStart[i], windSize[i]).dot(ZPy.segment(windStart[i], windSize[i]) + rcorr.segment(windStart[i], windSize[i])))/n[i];
                //                windEnd = windStart[i] + windSize[i];
                //                varei[i] = tss[i];
                //                for (j=windStart[i]; j<windEnd; ++j) {
                //                    if (valuesPtr[j]) varei[i] -= valuesPtr[j]*(ZPy[j] + rcorr[j]);
                //                }
                //                varei[i] /= n[i];
                //            }
                
                varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
                
                rhs = rcorr[i] + ZPZdiag[i]*oldSample;
                rhs /= varei;
                invLhs = 1.0f/(ZPZdiag[i]/varei + 1.0f/sigmaSq[i]);
                uhat = invLhs*rhs;
                logDelta1 = 0.5*(logf(invLhs) - logf(sigmaSq[i]) + uhat*rhs) + logPi;
                logDelta0 = logPiComp;
                probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
                
            }
            
            //cout << rhs << " " << invLhs << " " << logDelta1 << " " << logSigmaSq << " " << sigmaSq << endl;
            //            if (bernoulli.sample(probDelta1)) {
            if (urnd[i] < probDelta1) {
                //                valuesPtr[i] = normal.sample(uhat, invLhs);
                valuesPtr[i] = uhat + nrnd[i]*sqrtf(invLhs);
                //                rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*(oldSample - values[i]);
                float sampleDiff = oldSample - valuesPtr[i];
                for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {
                    rcorr[it.index()] += it.value() * sampleDiff;
                }
                betaSq[i] = valuesPtr[i]*valuesPtr[i];
                ssq[chr]  += valuesPtr[i]*valuesPtr[i];
                s2pq[chr] += snp2pq[i];
                ++nnz[chr];
            } else {
                if (oldSample) {
                    //                    rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*oldSample;
                    for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {
                        rcorr[it.index()] += it.value() * oldSample;
                    }
                }
                valuesPtr[i] = 0.0;
            }
        }
    }
    
    //cout << ssq << " " << nnz << endl;
    
    sumSq = 0.0;
    sum2pq = 0.0;
    numNonZeros = 0;
    nnzPerChr.setZero(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        sumSq += ssq[i];
        sum2pq += s2pq[i];
        numNonZeros += nnz[i];
        nnzPerChr[i] = nnz[i];
    }
    
    values = VectorXf::Map(valuesPtr, size);
}

void ApproxBayesB::SnpEffects::sampleFromFC(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const VectorXf &snp2pq, const VectorXf &LDsamplVar,
                                            const VectorXf &sigmaSq, const float pi, const float vare, const float varg, const float ps, const float overdispersion) {
    long numChr = chromInfoVec.size();
    
    float ssq[numChr], nnz[numChr], s2pq[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(nnz,0,sizeof(float)*numChr);
    memset(s2pq,0,sizeof(float)*numChr);
    
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking precents the writing from other threads
    
    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
#pragma omp parallel for
    for (unsigned chr=0; chr<numChr; ++chr) {
        //cout << " thread " << omp_get_thread_num() << " chr " << chr << endl;
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned windEnd, j;
        
        float oldSample;
        float rhs, invLhs, uhat;
        float logDelta0, logDelta1, probDelta1;
        float logPi = log(pi);
        float logPiComp = log(1.0-pi);
//        float logSigmaSq = log(sigmaSq);
//        float invSigmaSq = 1.0f/sigmaSq;
        float varei;
        
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            
            oldSample = valuesPtr[i];
            
            if (leaveout[i]) {
                probDelta1 = 0;
            }
            else {
                
                //cout << i << " " << chrStart << " " << chrEnd << " " << ssq << " " << nnz << endl;
                // cout << "Varei i " << varei[i] << endl;
                //            if (!(iter % 100)) {
                //                //float varei = (sse[i] - values.segment(windStart[i], windSize[i]).dot(ZPy.segment(windStart[i], windSize[i]) + rcorr.segment(windStart[i], windSize[i])))/n[i];
                //                windEnd = windStart[i] + windSize[i];
                //                varei[i] = tss[i];
                //                for (j=windStart[i]; j<windEnd; ++j) {
                //                    if (valuesPtr[j]) varei[i] -= valuesPtr[j]*(ZPy[j] + rcorr[j]);
                //                }
                //                varei[i] /= n[i];
                //                // cout << "Varei 100 " << varei[i] << endl;
                //            }
                // varei[i] = tss[i] / n[i];
                //            varei = se[i]*se[i]*ZPZdiag[i];
                
                varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
                
                rhs = rcorr[i] + ZPZdiag[i]*oldSample;
                rhs /= varei;
                invLhs = 1.0f/(ZPZdiag[i]/varei + 1.0f/sigmaSq[i]);
                uhat = invLhs*rhs;
                logDelta1 = 0.5*(logf(invLhs) - logf(sigmaSq[i]) + uhat*rhs) + logPi;
                logDelta0 = logPiComp;
                probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
                //cout << rhs << " " << invLhs << " " << logDelta1 << " " << logSigmaSq << " " << sigmaSq << endl;
                
            }
            
            //            if (bernoulli.sample(probDelta1)) {
            if (urnd[i] < probDelta1) {
                //                valuesPtr[i] = normal.sample(uhat, invLhs);
                valuesPtr[i] = uhat + nrnd[i]*sqrtf(invLhs);
                rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*(oldSample - valuesPtr[i]);
                betaSq[i] = valuesPtr[i]*valuesPtr[i];
                ssq[chr] += valuesPtr[i]*valuesPtr[i];
                s2pq[chr] += snp2pq[i];
                ++nnz[chr];
            } else {
                if (oldSample) rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*oldSample;
                valuesPtr[i] = 0.0;
            }
        }
    }
    // cout << "Varei 1 max" << varei.maxCoeff() << endl;
    //cout << ssq << " " << nnz << endl;
    
    sumSq = 0.0;
    sum2pq = 0.0;
    numNonZeros = 0.0;
    nnzPerChr.setZero(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        sumSq += ssq[i];
        sum2pq += s2pq[i];
        numNonZeros += nnz[i];
        nnzPerChr[i] = nnz[i];
    }
    
    values = VectorXf::Map(valuesPtr, size);
}

void ApproxBayesB::sampleUnknowns() {
    static int iter = 0;
    unsigned cnt=0;
    do {
        if (sparse)
        snpEffects.sampleFromFC(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, data.se, data.tss, varei,
                                data.n, data.snp2pq, data.LDsamplVar, sigmaSq.values, pi.value, vare.value, varg.value, ps.value, overdispersion);
        else
        snpEffects.sampleFromFC(rcorr, data.ZPZ, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, data.se, data.tss, varei,
                                data.n, data.snp2pq, data.LDsamplVar, sigmaSq.values, pi.value, vare.value, varg.value, ps.value, overdispersion);
        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros == 0);
    if (diagnose) nro.compute(rcorr, data.ZPZdiag, data.LDsamplVar, varg.value, vare.value, snpEffects.header, snpEffects.leaveout, data.ZPZsp, data.ZPy, snpEffects.values);

    if (robustMode) {
        if (noscale) {
            sigmaSq.value = varg.value/(data.snp2pq.array().sum()*pi.value);
        } else {
            sigmaSq.value = varg.value/(data.numIncdSnps*pi.value);  // LDpred2's parameterisation
        }
    } else {
        sigmaSq.sampleFromFC(snpEffects.betaSq);
    }

    if (estimatePi) pi.sampleFromFC(data.numIncdSnps, snpEffects.numNonZeros);
    nnzSnp.getValue(snpEffects.numNonZeros);
    sigmaSqG.compute(sigmaSq.value, snpEffects.sum2pq);
    
    covg.compute(data.ypy, snpEffects.values, data.ZPy, rcorr);
    varg.compute(snpEffects.values, data.ZPy, rcorr, covg.value);
    //    varg.value = sigmaSqG.value;
    vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, covg.value);
    hsq.compute(varg.value, vare.value);
    
    if (sparse)
    rounding.computeRcorr(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    else
    rounding.computeRcorr(data.ZPy, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    if (modelPS) ps.compute(rcorr, data.ZPZdiag, data.LDsamplVar, varg.value, vare.value, data.chisq);
    //    if (modelPS) ldScoreReg(data.chisq, data.LDscore, data.LDsamplVar, varg.value, vare.value, ps.value);
    
    //    if (sparse) {
    //        nnzgwas.compute(snpEffects.values, data.ZPZsp, data.ZPZdiag);
    //        pigwas.compute(nnzgwas.value, data.numIncdSnps);
    //    }
    
//    if (iter > 100 & !(iter % 10)) {
//        hsqMCMC.push_back(hsq.value);
//        if (!(iter % 1000)) checkHsq(hsqMCMC);
//    }

    ++iter;
}


void ApproxBayesS::SnpEffects::sampleFromFC(VectorXf &rcorr,const vector<SparseVector<float> > &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const float sigmaSq, const float pi, const float vare,
                                            const VectorXf &snp2pqPowS, const VectorXf &snp2pq, const VectorXf &LDsamplVar,
                                            const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n,
                                            const float varg, const float ps, const float overdispersion){
    long numChr = chromInfoVec.size();
    
    float ssq[numChr], nnz[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(nnz,0,sizeof(float)*numChr);
    //ssq.setZero(numChr);
    //nnz.setZero(numChr);
    
//    for (unsigned chr=0; chr<numChr; ++chr) {
//        ChromInfo *chromInfo = chromInfoVec[chr];
//        unsigned chrStart = chromInfo->startSnpIdx;
//        unsigned chrEnd   = chromInfo->endSnpIdx;
//        if (iter==0) {
//            cout << "chr " << chr+1 << " start " << chrStart << " end " << chrEnd << endl;
//        }
//    }
//    if (iter==0) cout << endl;
    
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    

#pragma omp parallel for
    for (unsigned chr=0; chr<numChr; ++chr) {
        //cout << " thread " << omp_get_thread_num() << " chr " << chr << endl;
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        
        float oldSample;
        float rhs, invLhs, uhat;
        float logDelta0, logDelta1, probDelta1;
        float logPi = log(pi);
        float logPiComp = log(1.0-pi);
        float invSigmaSq = 1.0f/sigmaSq;
        float varei;
        
//        float vare = tss.segment(chrStart, chrSize).mean() - values.segment(chrStart, chrSize).dot(ZPy.segment(chrStart, chrSize)) - values.segment(chrStart, chrSize).dot(rcorr.segment(chrStart, chrSize));
//        vare /= n.segment(chrStart, chrSize).mean();

//        float vare = tss.mean() - values.dot(ZPy) - values.dot(rcorr);
//        vare /= n.mean();

        for (unsigned i=chrStart; i<=chrEnd; ++i) {

            oldSample = valuesPtr[i];

            if (leaveout[i]) {
                probDelta1 = 0;
            }
            else {
                
                varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
                
//                float a = rcorr[i]*rcorr[i]/(ZPZdiag[i]*varei);
//                if (a > 15) varei += a-15;
                
                //float varei = se[i]*se[i]*ZPZdiag[i];
                
                rhs  = rcorr[i] + ZPZdiag[i]*oldSample;
                rhs /= varei;
                invLhs = 1.0f/(ZPZdiag[i]/varei + invSigmaSq/snp2pqPowS[i]);
                uhat = invLhs*rhs;
                
                logDelta1 = 0.5*(logf(invLhs) - logf(snp2pqPowS[i]*sigmaSq) + uhat*rhs) + logPi;
                logDelta0 = logPiComp;
                
                probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
                
            }
            
//            if (bernoulli.sample(probDelta1)) {
            if (urnd[i] < probDelta1) {
//                valuesPtr[i] = normal.sample(uhat, invLhs);

                valuesPtr[i] = uhat + nrnd[i]*sqrtf(invLhs);

//                if (iter < 10) {
//                valuesPtr[i] = uhat + nrnd[i]*sqrtf(invLhs);
//                }
//                else {
//                unsigned cnt = 0;
//                do {
//                valuesPtr[i] = uhat + nrnd[i]*sqrtf(invLhs);
////                    if(++cnt == 100) throw("Error: The effect of SNP " + header[i] + " is larger than 10 times the SD of random effect distribution!");
//                    if(++cnt == 100) {valuesPtr[i] -= uhat; leaveout[i] = 1; break;}
//                } while (valuesPtr[i] < 10*sqrt(sigmaSq));
//                }
                
                float sampleDiff = oldSample - valuesPtr[i];
                for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {
                    rcorr[it.index()] += it.value() * sampleDiff;
                }
                ssq[chr] += valuesPtr[i]*valuesPtr[i]/snp2pqPowS[i];
                ++nnz[chr];
                
            } else {
                if (oldSample) {
                    for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {
                        rcorr[it.index()] += it.value() * oldSample;
                    }
                }
                valuesPtr[i] = 0.0;
            }
        }
    }
    
    //wtdSumSq = ssq.sum();
    //numNonZeros = nnz.sum();
    wtdSumSq = 0.0;
    numNonZeros = 0;
    nnzPerChr.setZero(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        wtdSumSq += ssq[i];
        numNonZeros += nnz[i];
        nnzPerChr[i] = nnz[i];
    }

    values = VectorXf::Map(valuesPtr, size);
}

void ApproxBayesS::SnpEffects::sampleFromFC(VectorXf &rcorr,const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const float sigmaSq, const float pi, const float vare,
                                            const VectorXf &snp2pqPowS, const VectorXf &snp2pq, const VectorXf &LDsamplVar,
                                            const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n,
                                            const float varg, const float ps, const float overdispersion){
    long numChr = chromInfoVec.size();
    
    float ssq[numChr], nnz[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(nnz,0,sizeof(float)*numChr);
    //ssq.setZero(numChr);
    //nnz.setZero(numChr);
    
//    for (unsigned chr=0; chr<numChr; ++chr) {
//        ChromInfo *chromInfo = chromInfoVec[chr];
//        unsigned chrStart = chromInfo->startSnpIdx;
//        unsigned chrEnd   = chromInfo->endSnpIdx;
//        if (iter==0) {
//            cout << "chr " << chr+1 << " start " << chrStart << " end " << chrEnd << endl;
//        }
//    }
//    if (iter==0) cout << endl;
    
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking precents the writing from other threads
    
    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }

#pragma omp parallel for
    for (unsigned chr=0; chr<numChr; ++chr) {
        //cout << " thread " << omp_get_thread_num() << " chr " << chr << endl;
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        
        float oldSample;
        float rhs, invLhs, uhat;
        float logDelta0, logDelta1, probDelta1;
        float logPi = log(pi);
        float logPiComp = log(1.0-pi);
        float invSigmaSq = 1.0f/sigmaSq;
        float varei;
        
        for (unsigned i=chrStart; i<=chrEnd; ++i) {

            oldSample = valuesPtr[i];

            if (leaveout[i]) {
                probDelta1 = 0;
            }
            else {
                
                varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
                
                //float varei = se[i]*se[i]*ZPZdiag[i];
                
                rhs  = rcorr[i] + ZPZdiag[i]*oldSample;
                rhs /= varei;
                invLhs = 1.0f/(ZPZdiag[i]/varei + invSigmaSq/snp2pqPowS[i]);
                uhat = invLhs*rhs;
                
                logDelta1 = 0.5*(logf(invLhs) - logf(snp2pqPowS[i]*sigmaSq) + uhat*rhs) + logPi;
                logDelta0 = logPiComp;
                
                probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
                
            }
            
//            if (bernoulli.sample(probDelta1)) {
            if (urnd[i] < probDelta1) {
//                valuesPtr[i] = normal.sample(uhat, invLhs);
                valuesPtr[i] = uhat + nrnd[i]*sqrtf(invLhs);
                rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*(oldSample - valuesPtr[i]);
                ssq[chr] += valuesPtr[i]*valuesPtr[i]/snp2pqPowS[i];
                ++nnz[chr];
            } else {
                if (oldSample) {
                    rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*oldSample;
                }
                valuesPtr[i] = 0.0;
            }
        }
    }
    
    //wtdSumSq = ssq.sum();
    //numNonZeros = nnz.sum();
    wtdSumSq = 0.0;
    numNonZeros = 0;
    nnzPerChr.setZero(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        wtdSumSq += ssq[i];
        numNonZeros += nnz[i];
        nnzPerChr[i] = nnz[i];
    }

    values = VectorXf::Map(valuesPtr, size);
}

void ApproxBayesS::SnpEffects::sampleFromFC(const VectorXf &ZPy,const MatrixXf &Z, const VectorXf &ZPZdiag,
                                            const float sigmaSq, const float pi, const float vare,
                                            const VectorXf &snp2pqPowS, const VectorXf &snp2pq, VectorXf &ghat){
    wtdSumSq = 0.0;
    numNonZeros = 0;
    
//    cout << "size " << size << " ZPy " << ZPy.size() << " Z " << Z.rows() << " " << Z.cols() << " ZPZdiag " << ZPZdiag.size() << " 2pq " << snp2pq.size() << " ghat " << ghat.size() << endl;
    
//    cout << ZPZdiag.transpose() << endl;
//    cout << ZPy.transpose() << endl;
    
    float oldSample;
    float rhs, invLhs, uhat;
    float logDelta0, logDelta1, probDelta1;
    float logPi = log(pi);
    float logPiComp = log(1.0-pi);
    float invVare = 1.0f/vare;
    float invSigmaSq = 1.0f/sigmaSq;
    
    for (unsigned i=0; i<size; ++i) {
        oldSample = values[i];
        
        rhs = ZPy[i] + ZPZdiag[i]*oldSample - Z.col(i).dot(ghat);
        
        //rhs = ZPy[i] + ZPZdiag[i]*oldSample - (Z.col(i).transpose()*Z).dot(values);
        
//        cout << ZPy[i] << " " << ZPZdiag[i]*oldSample << " " << Z.col(i).dot(ghat) << endl;
        rhs *= invVare;
        invLhs = 1.0f/(ZPZdiag[i]*invVare + invSigmaSq/snp2pqPowS[i]);
        uhat = invLhs*rhs;
        logDelta1 = 0.5*(logf(invLhs) - logf(snp2pqPowS[i]*sigmaSq) + uhat*rhs) + logPi;
        logDelta0 = logPiComp;
        probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
        
        if (bernoulli.sample(probDelta1)) {
            values[i] = normal.sample(uhat, invLhs);
            ghat  += Z.col(i) * (values[i] - oldSample);
            wtdSumSq += values[i]*values[i]/snp2pqPowS[i];
            ++numNonZeros;
        } else {
            if (oldSample) ghat -= Z.col(i) * oldSample;
            values[i] = 0.0;
        }
    }
}


void ApproxBayesS::SnpEffects::hmcSampler(VectorXf &rcorr, const VectorXf &ZPy, const vector<VectorXf> &ZPZ,
                                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                          const float sigmaSq, const float pi, const float vare, const VectorXf &snp2pqPowS){
    
    float stepSize = 0.001;
    unsigned numSteps = 1;
    
    
//#pragma omp parallel for   // this multi-thread may not work due to vector locking when write to the vector
    for (unsigned chr=0; chr<chromInfoVec.size(); ++chr) {
        //cout << " thread " << omp_get_thread_num() << " chr " << chr << endl;
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned chrSize  = chromInfo->size;
        
        VectorXf chrZPy = ZPy.segment(chrStart, chrSize);
        VectorXf chrSnp2pqPowS = snp2pqPowS.segment(chrStart, chrSize);
        VectorXi chrWindStart = windStart.segment(chrStart, chrSize);
        VectorXi chrWindSize = windSize.segment(chrStart, chrSize);
        chrWindStart.array() -= chrStart;
        
        
        VectorXf delta;
        delta.setZero(chrSize);
        for (unsigned i=chrStart, j=0; i<=chrEnd; ++i) {
            if (values[i]) {
                delta[j++] = 1;
            }
        }
        
        
        VectorXf curr = values.segment(chrStart, chrSize);
        VectorXf curr_p(chrSize);
        
        for (unsigned i=0; i<chrSize; ++i) {
            curr_p[i] = Stat::snorm();
        }
        
        VectorXf cand = curr.cwiseProduct(delta);
        // Make a half step for momentum at the beginning
        VectorXf rc = chrZPy;
        VectorXf cand_p = curr_p.cwiseProduct(delta) - 0.5*stepSize * gradientU(curr, rc, chrZPy, ZPZ, chrWindStart, chrWindSize, chrStart, chrSize, sigmaSq, vare, chrSnp2pqPowS).cwiseProduct(delta);
        
        for (unsigned i=0; i<numSteps; ++i) {
            cand += stepSize * cand_p.cwiseProduct(delta);
            if (i < numSteps-1) {
                cand_p -= stepSize * gradientU(cand, rc, chrZPy, ZPZ, chrWindStart, chrWindSize, chrStart, chrSize, sigmaSq, vare, chrSnp2pqPowS).cwiseProduct(delta);
            } else {
                cand_p -= 0.5* stepSize * gradientU(cand, rc, chrZPy, ZPZ, chrWindStart, chrWindSize, chrStart, chrSize, sigmaSq, vare, chrSnp2pqPowS).cwiseProduct(delta);
            }
        }
        
        float curr_H = computeU(curr, rcorr.segment(chrStart, chrSize), chrZPy, sigmaSq, vare, chrSnp2pqPowS) + 0.5*curr_p.squaredNorm();
        float cand_H = computeU(cand, rc, chrZPy, sigmaSq, vare, chrSnp2pqPowS) + 0.5*cand_p.squaredNorm();
        
        if (Stat::ranf() < exp(curr_H-cand_H)) {  // accept
            values.segment(chrStart, chrSize) = cand;
            rcorr.segment(chrStart, chrSize) = rc;
            ++mhr;
            //cout << "accept " << curr_H << " " << cand_H << " " << exp(curr_H-cand_H) << endl;
        } else {
            //cout << "reject!!" << endl;
        }
    }
    
    sumSq = values.squaredNorm();
    //numNonZeros = size;
    
    for (unsigned i=0; i<size; ++i) {
        if(values[i]) ++numNonZeros;
    }
    //cout << sumSq << " " << nnz << " " << numNonZeros << endl;
    
    //cout << values.head(10).transpose() << endl;
    
//    if (!(++cnt % 100) && myMPI::rank==0) {
//        float ar = mhr/float(cnt);
//        if      (ar < 0.5) cout << "Warning: acceptance rate for SNP effects is too low "  << ar << endl;
//        else if (ar > 0.9) cout << "Warning: acceptance rate for SNP effects is too high " << ar << endl;
//    }
    
}

VectorXf ApproxBayesS::SnpEffects::gradientU(const VectorXf &effects, VectorXf &rcorr, const VectorXf &ZPy, const vector<VectorXf> &ZPZ,
                                             const VectorXi &windStart, const VectorXi &windSize, const unsigned chrStart, const unsigned chrSize,
                                             const float sigmaSq, const float vare, const VectorXf &snp2pqPowS){
    rcorr = ZPy;
    for (unsigned i=0; i<chrSize; ++i) {
        if (effects[i]) {
            rcorr.segment(windStart[i], windSize[i]) -= ZPZ[chrStart+i]*effects[i];
        }
    }
    return -rcorr/vare + effects.cwiseProduct(snp2pqPowS.cwiseInverse())/sigmaSq;
}

float ApproxBayesS::SnpEffects::computeU(const VectorXf &effects, const VectorXf &rcorr, const VectorXf &ZPy,                                             const float sigmaSq, const float vare, const VectorXf &snp2pqPowS){
    return -0.5f/vare*effects.dot(ZPy+rcorr) + 0.5/sigmaSq*effects.cwiseProduct(snp2pqPowS.cwiseInverse()).squaredNorm();
}

void ApproxBayesS::MeanEffects::sampleFromFC(const vector<SparseVector<float> > &ZPZ, const VectorXf &snpEffects, const VectorXf &snp2pq, const float vare, VectorXf &rcorr) {
    long numSnps = snpEffects.size();
    VectorXf snp2pqPowSmuDelta = snp2pqPowSmu;
    for (unsigned i=0; i<numSnps; ++i) {
        if (!snpEffects[i]) snp2pqPowSmuDelta[i] = 0;
    }
    float rhs = snp2pqPowSmuDelta.dot(rcorr);
    float lhs = 0.0;
    float oldSample = value;
    for (unsigned i=0; i<numSnps; ++i) {
        lhs += snp2pqPowSmuDelta[i]*(ZPZ[i].dot(snp2pqPowSmuDelta));
    }
    float invLhs = 1.0f/lhs;
    float bhat = invLhs*rhs;
    value = Normal::sample(bhat, invLhs*vare);
    snp2pqPowSmu = snp2pq.array().pow(value);
    rcorr += snp2pqPowSmu * (oldSample - value);
}

void ApproxBayesS::Smu::sampleFromFC(const vector<SparseVector<float> > &ZPZ, const VectorXf &snpEffects, const VectorXf &snp2pq, const float vare, VectorXf &snp2pqPowSmu, VectorXf &rcorr) {
    // random walk MH algorithm
    long numSnps = snpEffects.size();

    float curr = value;
    float cand = Normal::sample(value, varProp);
    
    VectorXf snp2pqPowSmuDeltaCurr = snp2pqPowSmu;
    VectorXf snp2pqPowSmuCand = snp2pq.array().pow(cand);
    VectorXf snp2pqPowSmuDeltaCand = snp2pqPowSmuCand;
    
    for (unsigned i=0; i<numSnps; ++i) {
        if (!snpEffects[i]) {
            snp2pqPowSmuDeltaCurr[i] = 0;
            snp2pqPowSmuDeltaCand[i] = 0;
        }
    }
    
    float rhsCurr = snp2pqPowSmu.dot(rcorr);
    float lhsCurr = 0.0;
    for (unsigned i=0; i<numSnps; ++i) {
        lhsCurr += snp2pqPowSmuDeltaCurr[i]*(ZPZ[i].dot(snp2pqPowSmuDeltaCurr));
    }

    float rhsCand = snp2pqPowSmuDeltaCand.dot(rcorr);
    float lhsCand = 0.0;
    for (unsigned i=0; i<numSnps; ++i) {
        lhsCand += snp2pqPowSmuDeltaCand[i]*(ZPZ[i].dot(snp2pqPowSmuDeltaCand));
    }
    
    float logCurr = -0.5f*(-log(lhsCurr) + rhsCurr/(lhsCurr*vare)) + curr*curr;
    float logCand = -0.5f*(-log(lhsCand) + rhsCand/(lhsCand*vare)) + cand*cand;
    
    if (Stat::ranf() < exp(logCand-logCurr)) {  // accept
        value = cand;
        snp2pqPowSmu = snp2pqPowSmuCand;
        ar.count(1, 0.1, 0.5);
    } else {
        ar.count(0, 0.1, 0.5);
    }
    
    if (!(ar.cnt % 10)) {
        if      (ar.value < 0.2) varProp *= 0.8;
        else if (ar.value > 0.5) varProp *= 1.2;
    }
    
    tuner.value = varProp;

}


void ApproxBayesS::sampleUnknowns(){
    
//    if (iter==0) vareiMean.setZero(data.numIncdSnps); ///TMP

//    fixedEffects.sampleFromFC(data.XPX, data.XPXdiag, data.ZPX, data.XPy, snpEffects.values, vare.value, rcorr);
    
//    tauSq.sampleFromFC(data.y, ghat);
    
    unsigned cnt=0;
    do {
        if (sparse) {
            snpEffects.sampleFromFC(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, sigmaSq.value, pi.value, vare.value,
                                    snp2pqPowS, data.snp2pq, data.LDsamplVar, data.se, data.tss, varei, data.n, varg.value, ps.value, overdispersion);
        } else {
            snpEffects.sampleFromFC(rcorr, data.ZPZ, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, sigmaSq.value, pi.value, vare.value,
                                    snp2pqPowS, data.snp2pq, data.LDsamplVar, data.se, data.tss, varei, data.n, varg.value, ps.value, overdispersion);
        }
        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros == 0);

    if (diagnose) nro.compute(rcorr, data.ZPZdiag, data.LDsamplVar, varg.value, vare.value, snpEffects.header, snpEffects.leaveout, data.ZPZsp, data.ZPy, snpEffects.values);

    if (estimatePi) pi.sampleFromFC(data.numIncdSnps, snpEffects.numNonZeros);
    
    if (estimateEffectMean) {
        mu.sampleFromFC(data.ZPZsp, snpEffects.values, data.snp2pq, vare.value, rcorr);
        Su.sampleFromFC(data.ZPZsp, snpEffects.values, data.snp2pq, vare.value, mu.snp2pqPowSmu, rcorr);
    }
    
    nnzSnp.getValue(snpEffects.numNonZeros);

    covg.compute(data.ypy, snpEffects.values, data.ZPy, rcorr);
    varg.compute(snpEffects.values, data.ZPy, rcorr, covg.value);
//    varg.value = sigmaSqG.value;
    vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, covg.value);
    //vare.value = data.ypy/data.numKeptInds;
    hsq.compute(varg.value, vare.value);
    
    if (robustMode) {
        S.sampleFromFC2(snpEffects.numNonZeros, snpEffects.values, data.snp2pq, snp2pqPowS, logSnp2pq, varg.value, snpEffects.sum2pqSplusOne);
        sigmaSq.value = varg.value/snpEffects.sum2pqSplusOne;
    } else {
        sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
        S.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros, sigmaSq.value, snpEffects.values, data.snp2pq, snp2pqPowS, logSnp2pq, genVarPrior, sigmaSq.scale, snpEffects.sum2pqSplusOne);
    }
    
    sigmaSqG.compute(sigmaSq.value, snpEffects.sum2pqSplusOne);

    if (iter >= 2000) sigmaSq.scale = scalePrior;
    scale.getValue(sigmaSq.scale);
    // cout << "sigmaSq.scale " << sigmaSq.scale << endl;

    if (sparse)
        rounding.computeRcorr(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    else
        rounding.computeRcorr(data.ZPy, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    
    if (modelPS) ps.compute(rcorr, data.ZPZdiag, data.LDsamplVar, varg.value, vare.value, data.chisq);
    
//    if (!(iter % 100)) vareiMean += varei;  ///TMP

    if (++iter < 2000) {
        genVarPrior += (varg.value - genVarPrior)/iter;
        scalePrior  += (sigmaSq.scale - scalePrior)/iter;
    }
    
//    if (sparse) {
//        nnzgwas.compute(snpEffects.values, data.ZPZsp, data.ZPZdiag);
//        pigwas.compute(nnzgwas.value, data.numIncdSnps);
//    }
    
}



void ApproxBayesST::SnpEffects::sampleFromFC(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                             const vector<ChromInfo *> &chromInfoVec, const VectorXf &LDsamplVar, const ArrayXf &hSlT, const VectorXf &snp2pq,
                                             const float sigmaSq, const float pi, const float vare, const float varg,
                                             const float ps, const float overdispersion) {
    long numChr = chromInfoVec.size();
    
    float ssq[numChr], nnz[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(nnz,0,sizeof(float)*numChr);
    
    sum2pqhSlT = 0.0;
    
//    for (unsigned chr=0; chr<numChr; ++chr) {
//        ChromInfo *chromInfo = chromInfoVec[chr];
//        unsigned chrStart = chromInfo->startSnpIdx;
//        unsigned chrEnd   = chromInfo->endSnpIdx;
//        if (iter==0) {
//            cout << "chr " << chr+1 << " start " << chrStart << " end " << chrEnd << endl;
//        }
//    }
//    if (iter==0) cout << endl;
    
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads
    
    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
#pragma omp parallel for
    for (unsigned chr=0; chr<numChr; ++chr) {
        //cout << " thread " << omp_get_thread_num() << " chr " << chr << endl;
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        
        float oldSample;
        float rhs, invLhs, uhat;
        float logDelta0, logDelta1, probDelta1;
        float logPi = log(pi);
        float logPiComp = log(1.0-pi);
        float invSigmaSq = 1.0f/sigmaSq;
        float varei;
        
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            oldSample = valuesPtr[i];
            varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
            rhs  = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            invLhs = 1.0f/(ZPZdiag[i]/varei + invSigmaSq/hSlT[i]);
            uhat = invLhs*rhs;
            
            logDelta1 = 0.5*(logf(invLhs) - logf(hSlT[i]*sigmaSq) + uhat*rhs) + logPi;
            logDelta0 = logPiComp;
            
            probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
            
            if (urnd[i] < probDelta1) {
                valuesPtr[i] = uhat + nrnd[i]*sqrtf(invLhs);
                float sampleDiff = oldSample - valuesPtr[i];
                for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {
                    rcorr[it.index()] += it.value() * sampleDiff;
                }
                ssq[chr] += valuesPtr[i]*valuesPtr[i]/hSlT[i];
                sum2pqhSlT += snp2pq[i]*hSlT[i];
                ++nnz[chr];
                
            } else {
                if (oldSample) {
                    for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {
                        rcorr[it.index()] += it.value() * oldSample;
                    }
                }
                valuesPtr[i] = 0.0;
            }
        }
    }
    
    wtdSumSq = 0.0;
    numNonZeros = 0;
    nnzPerChr.setZero(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        wtdSumSq += ssq[i];
        numNonZeros += nnz[i];
        nnzPerChr[i] = nnz[i];
    }
    
    values = VectorXf::Map(valuesPtr, size);
}

void ApproxBayesST::Sp::sampleFromFC(const unsigned int numNonZeros, const float sigmaSq, const VectorXf &snpEffects,
                                     const VectorXf &snp2pq, const ArrayXf &logSnp2pq,
                                     const VectorXf &ldsc, const ArrayXf &logLdsc,
                                     const float varg, float &scale, float &T, ArrayXf &hSlT) {
    unsigned nnz = 0;
    for (unsigned i=0; i<numSnps; ++i)
        if (snpEffects[i]) ++nnz;
    
    // Prepare
    ArrayXf snpEffectSqDelta1(nnz);
    ArrayXf snp2pqDelta1(nnz);
    ArrayXf ldscDelta1(nnz);
    ArrayXf logSnp2pqDelta1(nnz);
    ArrayXf logLdscDelta1(nnz);
    
    for (unsigned i=0, j=0; i<numSnps; ++i) {
        if (snpEffects[i]) {
            snpEffectSqDelta1[j] = snpEffects[i]*snpEffects[i];
            snp2pqDelta1[j] = snp2pq[i];
            ldscDelta1[j] = ldsc[i];
            logSnp2pqDelta1[j] = logSnp2pq[i];
            logLdscDelta1[j] = logLdsc[i];
            ++j;
        }
    }
    
    float snp2pqLogSumDelta1 = logSnp2pqDelta1.sum();
    float ldscLogSumDelta1 = logLdscDelta1.sum();
    
    Vector2f curr; curr << value, T;
    Vector2f curr_p; curr_p << Stat::snorm(), Stat::snorm();
    Vector2f cand = curr;
    
    // Make a half step for momentum at the beginning
    Vector2f cand_p = curr_p - 0.5*stepSize * gradientU(curr,  snpEffectSqDelta1, sigmaSq, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, ldscLogSumDelta1, ldscDelta1, logLdscDelta1);
    
    for (unsigned i=0; i<numSteps; ++i) {
        // Make a full step for the position
        cand += stepSize * cand_p;
        if (i < numSteps-1) {
            // Make a full step for the momentum, except at end of trajectory
            cand_p -= stepSize * gradientU(cand, snpEffectSqDelta1, sigmaSq, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, ldscLogSumDelta1, ldscDelta1, logLdscDelta1);
        } else {
            // Make a half step for momentum at the end
            cand_p -= 0.5*stepSize * gradientU(cand, snpEffectSqDelta1, sigmaSq, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, ldscLogSumDelta1, ldscDelta1, logLdscDelta1);
        }
        //cout << i << " " << cand << endl;
    }
    
    // Evaluate potential (negative log posterior) and kinetic energies at start and end of trajectory
    float curr_H = computeU(curr, snpEffectSqDelta1, sigmaSq, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, ldscLogSumDelta1, ldscDelta1, logLdscDelta1) + 0.5*curr_p.squaredNorm();
    float cand_H = computeU(cand, snpEffectSqDelta1, sigmaSq, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, ldscLogSumDelta1, ldscDelta1, logLdscDelta1) + 0.5*cand_p.squaredNorm();
        
    if (Stat::ranf() < exp(curr_H-cand_H)) {  // accept
        value = cand[0];
        T = cand[1];
        scale = varg/(snp2pqDelta1.pow(1.0+value)*(ldscDelta1.pow(T))).sum();
        hSlT = snp2pq.array().pow(value) * ldsc.array().pow(T);
        ar.count(1, 0.5, 0.9);
    } else {
        ar.count(0, 0.5, 0.9);
    }
    
    if (!(ar.cnt % 10)) {
        if      (ar.value < 0.6) stepSize *= 0.8;
        else if (ar.value > 0.8) stepSize *= 1.2;
    }
    
    if (ar.consecRej > 20) stepSize *= 0.8;
    
    tuner.value = stepSize;
}

Vector2f ApproxBayesST::Sp::gradientU(const Vector2f &ST, const ArrayXf &snpEffectSq, const float sigmaSq,
                                      const float snp2pqLogSum, const ArrayXf &snp2pq, const ArrayXf &logSnp2pq,
                                      const float ldscLogSum, const ArrayXf &ldsc, const ArrayXf &logLdsc) {
    float S = ST[0];
    float T = ST[1];
    long size = snp2pq.size();
    long chunkSize = size/omp_get_max_threads();
    ArrayXf snp2pqPowS(size);
    ArrayXf ldscPowT(size);
#pragma omp parallel for schedule(dynamic, chunkSize)
    for (unsigned i=0; i<size; ++i) {
        snp2pqPowS[i] = powf(snp2pq[i], S);
        ldscPowT[i] = powf(ldsc[i], T);
    }
    ArrayXf hSlT = snp2pqPowS*ldscPowT;
    Vector2f ret;
    ret[0] = 0.5*snp2pqLogSum - 0.5/sigmaSq*(snpEffectSq*logSnp2pq/hSlT).sum() + S;
    ret[1] = 0.5*ldscLogSum - 0.5/sigmaSq*(snpEffectSq*logLdsc/hSlT).sum() + T;
    return ret;
}

float ApproxBayesST::Sp::computeU(const Vector2f &ST, const ArrayXf &snpEffectSq, const float sigmaSq,
                                  const float snp2pqLogSum, const ArrayXf &snp2pq, const ArrayXf &logSnp2pq,
                                  const float ldscLogSum, const ArrayXf &ldsc, const ArrayXf &logLdsc) {
    float S = ST[0];
    float T = ST[1];
    ArrayXf snp2pqPowS = snp2pq.pow(S);
    ArrayXf ldscPowT = ldsc.pow(T);
    ArrayXf hSlT = snp2pqPowS*ldscPowT;
    return 0.5*S*snp2pqLogSum + 0.5*T*ldscLogSum + 0.5/sigmaSq*(snpEffectSq/hSlT).sum() + 0.5*S*S + 0.5*T*T;
}

void ApproxBayesST::Tp::sampleFromFC(const unsigned int numNonZeros, const float sigmaSq, const VectorXf &snpEffects,
                                     const VectorXf &snp2pq, const VectorXf &ldsc, const ArrayXf &logLdsc,
                                     const float varg, float &scale, ArrayXf &hSlT) {
    unsigned nnz = 0;
    for (unsigned i=0; i<numSnps; ++i)
        if (snpEffects[i]) ++nnz;
    
    // Prepare
    ArrayXf snpEffectSqDelta1(nnz);
    ArrayXf snp2pqDelta1(nnz);
    ArrayXf ldscDelta1(nnz);
    ArrayXf logLdscDelta1(nnz);
    
    for (unsigned i=0, j=0; i<numSnps; ++i) {
        if (snpEffects[i]) {
            snpEffectSqDelta1[j] = snpEffects[i]*snpEffects[i];
            snp2pqDelta1[j] = snp2pq[i];
            ldscDelta1[j] = ldsc[i];
            logLdscDelta1[j] = logLdsc[i];
            ++j;
        }
    }
    
    float ldscLogSumDelta1 = logLdscDelta1.sum();
    
    float curr = value;
    float curr_p = Stat::snorm();
    float cand = curr;
    
    // Make a half step for momentum at the beginning
    float cand_p = curr_p - 0.5*stepSize * gradientU(curr,  snpEffectSqDelta1, sigmaSq, ldscLogSumDelta1, ldscDelta1, logLdscDelta1);
    
    for (unsigned i=0; i<numSteps; ++i) {
        // Make a full step for the position
        cand += stepSize * cand_p;
        if (i < numSteps-1) {
            // Make a full step for the momentum, except at end of trajectory
            cand_p -= stepSize * gradientU(cand, snpEffectSqDelta1, sigmaSq, ldscLogSumDelta1, ldscDelta1, logLdscDelta1);
        } else {
            // Make a half step for momentum at the end
            cand_p -= 0.5*stepSize * gradientU(cand, snpEffectSqDelta1, sigmaSq, ldscLogSumDelta1, ldscDelta1, logLdscDelta1);
        }
        //cout << i << " " << cand << endl;
    }
    
    // Evaluate potential (negative log posterior) and kinetic energies at start and end of trajectory
    float curr_H = computeU(curr, snpEffectSqDelta1, sigmaSq, ldscLogSumDelta1, ldscDelta1, logLdscDelta1) + 0.5*curr_p*curr_p;
    float cand_H = computeU(cand, snpEffectSqDelta1, sigmaSq, ldscLogSumDelta1, ldscDelta1, logLdscDelta1) + 0.5*cand_p*cand_p;
    
    if (Stat::ranf() < exp(curr_H-cand_H)) {  // accept
        value = cand;
        scale = varg/(snp2pqDelta1.array()*ldscDelta1.pow(value)).sum();
        hSlT = ldsc.array().pow(value);
        ar.count(1, 0.5, 0.9);
    } else {
        ar.count(0, 0.5, 0.9);
    }
    
    if (!(ar.cnt % 10)) {
        if      (ar.value < 0.6) stepSize *= 0.8;
        else if (ar.value > 0.8) stepSize *= 1.2;
    }
    
    if (ar.consecRej > 20) stepSize *= 0.8;
    
    tuner.value = stepSize;
}

float ApproxBayesST::Tp::gradientU(const float &T, const ArrayXf &snpEffectSq, const float sigmaSq,
                                  const float ldscLogSum, const ArrayXf &ldsc, const ArrayXf &logLdsc) {
    long size = ldsc.size();
    long chunkSize = size/omp_get_max_threads();
    ArrayXf ldscPowT(size);
#pragma omp parallel for schedule(dynamic, chunkSize)
    for (unsigned i=0; i<size; ++i) {
        ldscPowT[i] = powf(ldsc[i], T);
    }
    return 0.5*ldscLogSum - 0.5/sigmaSq*(snpEffectSq*logLdsc/ldscPowT).sum() + T;
}

float ApproxBayesST::Tp::computeU(const float &T, const ArrayXf &snpEffectSq, const float sigmaSq,
                                 const float ldscLogSum, const ArrayXf &ldsc, const ArrayXf &logLdsc) {
    return 0.5*T*ldscLogSum + 0.5/sigmaSq*(snpEffectSq/ldsc.pow(T)).sum() + 0.5*T*T;

}


void ApproxBayesST::sampleUnknowns(){
    unsigned cnt=0;
    do {
        snpEffects.sampleFromFC(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.chromInfoVec,
                                data.LDsamplVar, hSlT, data.snp2pq, sigmaSq.value, pi.value, vare.value,
                                varg.value, ps.value, overdispersion);
        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros == 0);
    sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
    if (estimatePi) pi.sampleFromFC(data.numIncdSnps, snpEffects.numNonZeros);
    nnzSnp.getValue(snpEffects.numNonZeros);
    sigmaSqG.compute(sigmaSq.value, snpEffects.sum2pqhSlT);
    varg.compute(snpEffects.values, data.ZPy, rcorr, covg.value);
    vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, covg.value);
    hsq.compute(varg.value, vare.value);
    if (estimateS)
        S.sampleFromFC(snpEffects.numNonZeros, sigmaSq.value, snpEffects.values, data.snp2pq, logSnp2pq,
                       data.LDscore, logLdsc, varg.value, sigmaSq.scale, T.value, hSlT);
    else
        T.sampleFromFC(snpEffects.numNonZeros, sigmaSq.value, snpEffects.values, data.snp2pq, data.LDscore,
                       logLdsc, varg.value, sigmaSq.scale, hSlT);
    scale.getValue(sigmaSq.scale);
    rounding.computeRcorr(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    if (modelPS) ps.compute(rcorr, data.ZPZdiag, data.LDsamplVar, varg.value, vare.value, data.chisq);
//    nnzgwas.compute(snpEffects.values, data.ZPZsp, data.ZPZdiag);
//    pigwas.compute(nnzgwas.value, data.numIncdSnps);
}

void ApproxBayesST::sampleStartVal(){
    sigmaSq.sampleFromPrior();
    if (estimatePi) pi.sampleFromPrior();
    S.sampleFromPrior();
    T.sampleFromPrior();
    cout << "  Starting value for " << sigmaSq.label << ": " << sigmaSq.value << endl;
    if (estimatePi) cout << "  Starting value for " << pi.label << ": " << pi.value << endl;
    cout << "  Starting value for " << S.label << ": " << S.value << endl;
    cout << "  Starting value for " << T.label << ": " << T.value << endl;
    cout << endl;
}


// *******************************************************
// Bayes R - Approximate
// *******************************************************

void ApproxBayesR::sampleUnknowns(){
    //sigmaSq.value = 0.000275;   // TMP_JZ
    static int iter = 0;
//    fixedEffects.sampleFromFC(data.XPX, data.XPXdiag, data.ZPX, data.XPy, snpEffects.values, vare.value, rcorr);
    unsigned cnt=0;
    do {
        if (data.Z.size()) {
            snpEffects.sampleFromFC(data.ZPy, data.ZPZdiag, data.Z, data.Z.rows(), data.numKeptInds, sigmaSq.value, Pis.values, gamma.values, vare.value, snpStore, ghat, varg.value, originalModel);
        }
        else {
            if (lowRankModel) {
                snpEffects.sampleFromFC(wcorrBlocks, data.Qblocks, whatBlocks, data.keptLdBlockInfoVec, data.nGWASblock, vareBlk.values, sigmaSq.value, Pis.values, gamma.values, snpStore, varg.value, originalModel);
            }
            else if (sparse)
                snpEffects.sampleFromFC(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, data.se, data.tss, varei, data.n, data.snp2pq, data.LDsamplVar, sigmaSq.value, Pis.values, gamma.values, vare.value, snpStore,
                                        varg.value, ps.value, overdispersion, originalModel);
            else
                snpEffects.sampleFromFC(rcorr, data.ZPZ, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, data.se, data.tss, varei, data.n, data.snp2pq, data.LDsamplVar, sigmaSq.value, Pis.values, gamma.values, vare.value, snpStore,
                                        varg.value, ps.value, overdispersion, originalModel);
        }
        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros == 0);
        
    if (algorithm == cg) {
        snpEffects.adjustByCG(data.ZPy, data.ZPZsp, rcorr);
    }
    
    if (diagnose) nro.compute(rcorr, data.ZPZdiag, data.LDsamplVar, varg.value, vare.value, snpEffects.header, snpEffects.leaveout, data.ZPZsp, data.ZPy, snpEffects.values);
    
    if (robustMode) {
        if (noscale) {
            sigmaSq.value = varg.value/(data.snp2pq.array().sum()*gamma.values.dot(Pis.values));
        } else {
            sigmaSq.value = varg.value/(data.numIncdSnps*gamma.values.dot(Pis.values));  // LDpred2's parameterisation
        }
    } else if (originalModel) {
        sigmaSq.sampleFromFC(snpEffects.sumSq, snpEffects.numNonZeros);
    } else {
        if (estimateSigmaSq) sigmaSq.sampleFromFC(snpEffects.sumSq, snpEffects.numNonZeros);
    }
        
    if (estimatePi) Pis.sampleFromFC(snpStore);
    numSnps.getValues(snpStore);
    nnzSnp.getValue(snpEffects.numNonZeros);
    sigmaSqG.compute(sigmaSq.value, snpEffects.sum2pq);

    if (estimateHsq) {
        if (data.Z.size()) {   // TMP_JZ
            //        ghat.setZero(data.Z.rows());
            //        for (unsigned i=0; i<snpEffects.size; ++i) {
            //            if (snpEffects.values[i]) ghat += data.Z.col(i)*snpEffects.values[i];
            //        }
            varg.value = Gadget::calcVariance(ghat);
            float n_ref = data.Z.rows();
            float n_ratio = n_ref/float(data.numKeptInds);
            vare.value = (data.ypy*n_ratio - 2.0f*snpEffects.values.dot(data.ZPy)*n_ratio + ghat.dot(ghat))/n_ref;
            //cout << "varg " << varg.value << " vare " << vare.value << endl;
            //cout << "n_ratio " << n_ratio << " ypy " << data.ypy*n_ratio << " ypg " << 2.0f*snpEffects.values.dot(data.ZPy)*n_ratio << " gpg " << ghat.dot(ghat) << " n_ref " << n_ref << endl;
        }
        else if (lowRankModel) {
            vargBlk.compute(whatBlocks);
            vareBlk.sampleFromFC(wcorrBlocks, snpEffects.ssqBlocks, data.nGWASblock, data.numEigenvalBlock);
            varg.value = vargBlk.total;
            vare.value = vareBlk.mean;
        }
        else {
            covg.compute(data.ypy, snpEffects.values, data.ZPy, rcorr);
            varg.compute(snpEffects.values, data.ZPy, rcorr, covg.value);
            vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, covg.value);
        }
        
        //hsq.compute(varg.value, vare.value);
        hsq.value = varg.value / data.varPhenotypic;
    }
    
    if (iter >= 2000) sigmaSq.scale = scalePrior;
    scale.getValue(sigmaSq.scale);
    // cout << "iter " << iter << " scalePrior " << scalePrior << "sigmaSq.scale " << sigmaSq.scale << endl;

//    if (sparse)
//        rounding.computeRcorr(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
//    else
//        rounding.computeRcorr(data.ZPy, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    if (modelPS) ps.compute(rcorr, data.ZPZdiag, data.LDsamplVar, varg.value, vare.value, data.chisq);

    nnzSnp.getValue(snpEffects.numNonZeros);
    sigmaSqG.compute(sigmaSq.value, snpEffects.sum2pq);

//    numSnpVg.compute(snpEffects.values, data.ZPZdiag, varg.value, vare.nobs);
    if (originalModel) {
        Vgs.compute(snpEffects.values, data.ZPy, rcorr, snpEffects.snpset, varg.value, vare.nobs);
//        if (sparse)
//            Vgs.compute(snpEffects.values, data.ZPZsp, snpEffects.snpset, varg.value, vare.nobs);
//        else
//            Vgs.compute(snpEffects.values, data.ZPZ, snpEffects.snpset, varg.value, vare.nobs);
    }

    float scaleIteri = 0;
    if (++iter < 2000) {
        if (noscale)
        {
            scaleIteri = 0.5f * varg.value / (data.snp2pq.array().sum()*gamma.values.dot(Pis.values));
        } else
        {
            scaleIteri = 0.5f * varg.value / (data.snp2pq.size()*gamma.values.dot(Pis.values));
        }
        genVarPrior += (varg.value - genVarPrior)/iter;
        scalePrior += (scaleIteri - scalePrior)/iter;
    }
    
//    if (iter > 100 & !(iter % 10)) {
//        hsqMCMC.push_back(hsq.value);
//        if (!(iter % 1000)) checkHsq(hsqMCMC);
//    }
}

void ApproxBayesR::VgMixComps::compute(const VectorXf &snpEffects, const VectorXf &ZPy, const VectorXf &rcorr, const vector<vector<unsigned> > snpset, const float varg, const float nobs) {
    values.setZero(ndist);
    for (unsigned k=1; k<ndist; ++k) {
        unsigned size = snpset[k].size();
        for (unsigned j=0; j<size; ++j) {
            unsigned snpIdx = snpset[k][j];
            float varj = snpEffects[snpIdx] * (ZPy[snpIdx] - rcorr[snpIdx]);
            values[k] += varj;
        }
        values[k] /= varg * nobs;
        (*this)[k]->value = values[k];
    }
}


//void ApproxBayesR::VgMixComps::compute(const VectorXf &snpEffects, const vector<SparseVector<float> > &ZPZsp, const vector<vector<unsigned> > snpset, const float varg, const float nobs) {
//    values.setZero(ndist);
//    for (unsigned k=0; k<ndist; ++k) {
//        if (k!=zeroIdx && k!=minIdx) {
//            long numSnps = snpset[k].size();
//            float vargk = 0.0;
//            for (unsigned i=0, j=0; i<numSnps; ++i) {
//                unsigned ki = snpset[k][i];
//                unsigned kj = snpset[k][j];
//                for (SparseVector<float>::InnerIterator it(ZPZsp[ki]); it; ++it) {
//                    if (it.index() == kj) {
//                        vargk += snpEffects[ki]*snpEffects[kj]*it.value();
//                        kj = snpset[k][++j];
//                        if (j==numSnps) break;
//                    }
//                }
//            }
//            (*this)[k]->value = values[k] = vargk/(varg*nobs);
//        }
//    }
//    float sum = values.sum();
//    (*this)[minIdx]->value = values[minIdx] = 1.0 - sum;
//}
//
//void ApproxBayesR::VgMixComps::compute(const VectorXf &snpEffects, const vector<VectorXf> &ZPZ, const vector<vector<unsigned> > snpset, const float varg, const float nobs) {
//    values.setZero(ndist);
//    for (unsigned k=0; k<ndist; ++k) {
//        if (k!=zeroIdx && k!=minIdx) {
//            long numSnps = snpset[k].size();
//            float vargk = 0.0;
//            for (unsigned i=0; i<numSnps; ++i) {
//                unsigned ki = snpset[k][i];
//                for (unsigned j=0; j<numSnps; ++j) {
//                    unsigned kj = snpset[k][j];
//                    vargk += snpEffects[ki]*snpEffects[kj]*ZPZ[ki][kj];
//                }
//            }
//            (*this)[k]->value = values[k] = vargk/(varg*nobs);
//        }
//    }
//    float sum = values.sum();
//    (*this)[minIdx]->value = values[minIdx] = 1.0 - sum;
//}

// ==============================================================
// Sparse vector version
// ==============================================================

void ApproxBayesR::SnpEffects::sampleFromFC(VectorXf &rcorr, const vector<SparseVector<float>> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const VectorXf &snp2pq, const VectorXf &LDsamplVar,
                                            const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare, VectorXf &snpStore,
                                            const float varg, const float ps, const float overdispersion,
                                            const bool originalModel){
    // -----------------------------------------
    // Initialise the parameters in MCMC sampler
    // -----------------------------------------
    static unsigned iter = 0;
    long numChr = chromInfoVec.size();

    float ssq[numChr], s2pq[numChr], nnz[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(s2pq,0,sizeof(float)*numChr);
    memset(nnz,0, sizeof(float)*numChr);

    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
    // R specific parameters
    int ndist;
    VectorXf gp;
    snpStore.setZero(pis.size());
    // --------------------------------------------------------------------------------
    // Scale the variances in each of the normal distributions by the genetic variance
    // and initialise the class membership probabilities
    // --------------------------------------------------------------------------------
    ndist = pis.size();
    if (originalModel)
        gp = gamma * 0.01 * varg;
    else
        gp = gamma * sigmaSq;
    snpset.resize(ndist);
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
    }
    
    VectorXf invGamma = gamma.array().inverse();
    invGamma[0] = 0.0;
    
    lambdaVec.setZero(size);
    uhatVec.setZero(size);
    invGammaVec.setZero(size);
    deltaNZ.setZero(size);
    
    deltaNzIdx.clear();
    deltaNzIdx.reserve(size);
    
    // --------------------------------------------------------------------------------
    // Cycle over all variants in the window and sample the genetics effects
    // --------------------------------------------------------------------------------

//#pragma omp parallel for  // openmp is not working for SBayesR
    for (unsigned chr=0; chr<numChr; ++chr) 
    {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned windEnd, j;
        
        // R specific parameters
        int indistflag;
        double rhs, v1,  b_ls, ssculm, r;
        VectorXf ll, pll, snpindist, var_b_ls;
        ll.setZero(pis.size());
        pll.setZero(pis.size());

        float oldSample, varei;
        
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            oldSample = valuesPtr[i]; 
            varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
                        
//            varei = (tss[i] + oldSample*oldSample*ZPZdiag[i])/n[i];
//            float ssei = rcorr[i]*rcorr[i]/ZPZdiag[i];

//            float dfTilde = 10 + 1;
//            float scaleTilde = ssei + 10*(LDsamplVar[i]*varg + vare + ps + overdispersion);
//            Stat::InvChiSq invchisq;
//            varei = invchisq.sample(dfTilde, scaleTilde);
//
//            cout << i << " " << varei << endl;
            
            // ------------------------------
            // Derived Bayes R implementation
            // ------------------------------
            // ----------------------------------------------------
            // Add back the content for the corrected rhs for SNP k
            // ----------------------------------------------------
            rhs = rcorr[i] + ZPZdiag[i] * oldSample;
            // ------------------------------------------------------
            // Calculate the beta least squares updates and variances
            // ------------------------------------------------------
            b_ls = rhs / ZPZdiag[i];
            var_b_ls = gp.array() + varei / ZPZdiag[i];
            // ------------------------------------------------------
            // Calculate the likelihoods for each distribution
            // ------------------------------------------------------
            ll = (-1.0 / 2.0) * var_b_ls.array().log()  - (b_ls * b_ls)  / (2 * var_b_ls.array()) + pis.array().log();
            // --------------------------------------------------------------
            // Calculate probability that snp is in each of the distributions
            // in this iteration
            // --------------------------------------------------------------
            // pll = (ll.array().exp().cwiseProduct(pis.array())) / ((ll.array().exp()).cwiseProduct(pis.array())).sum();
            for (unsigned k=0; k<pis.size(); ++k) {
              pll[k] = 1.0 / (exp(ll.array() - ll[k])).sum();
            }
            // --------------------------------------------------------------
            // Sample the group based on the calculated probabilities
            // --------------------------------------------------------------
            ssculm = 0.0;
            r = urnd[i];
            indistflag = 1;
            for (int kk = 0; kk < ndist; kk++)
            {
                ssculm += pll(kk);
                if (r < ssculm)
                {
                    indistflag = kk + 1;
                    snpStore(kk) = snpStore(kk) + 1; 
                    break;
                }
            }
            snpset[indistflag-1].push_back(i);
            // --------------------------------------------------------------
            // Sample the effect given the group and adjust the rhs                                                                                                                 
            // --------------------------------------------------------------                                                                                                       
            if (indistflag != 1)                                                                                                                                                    
            {                                                                                                                                                                       
                v1 = ZPZdiag[i] + varei / gp((indistflag - 1));                                                                                                                     
//                valuesPtr[i] = normal.sample(rhs / v1, varei / v1);                                                                                                                 
                valuesPtr[i] = rhs / v1 + nrnd[i]*sqrtf(varei / v1);
                float sampleDiff = oldSample - valuesPtr[i];
                for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {                                                                                                     
                    rcorr[it.index()] += it.value() * sampleDiff;                                                                                                                   
                }                                                                                                                                                                   
                ssq[chr]  += (valuesPtr[i]*valuesPtr[i]) / gamma[indistflag - 1];
                s2pq[chr] += snp2pq[i];
                deltaNZ[i] = 1;
                ++nnz[chr];
                deltaNzIdx.push_back(i);
            } else {                                                                                                                                                                
                if (oldSample) {                                                                                                                                                    
                    for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {                                                                                                 
                        rcorr[it.index()] += it.value() * oldSample;                                                                                                                
                    }                                                                                                                                                               
                }                                                                                                                                                                   
                valuesPtr[i] = 0.0;                                                                                                                                                 
            }
            
            uhatVec[i] = rhs/v1;
            lambdaVec[i] = vare/gp[indistflag-1];
            invGammaVec[i] = invGamma[indistflag-1];
        }
    }
    // ---------------------------------------------------------------------
    // Tally up the effect sum of squares and the number of non-zero effects
    // ---------------------------------------------------------------------
    sumSq = 0.0;                                                                                                                                                                    
    sum2pq = 0.0;                                                                                                                                                                   
    numNonZeros = 0;                                                                                                                                                                
    nnzPerChr.setZero(numChr);                                                                                                                                                      
    for (unsigned i=0; i<numChr; ++i) {
        sumSq += ssq[i];
        sum2pq += s2pq[i];                                                                                                                                                          
        numNonZeros += nnz[i];                                                                                                                                                      
        nnzPerChr[i] = nnz[i];                                                                                                                                                      
    }
    ++iter;

    values = VectorXf::Map(valuesPtr, size);
}

// ==============================================================
// Vector of vectors version
// ==============================================================

void ApproxBayesR::SnpEffects::sampleFromFC(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const VectorXf &snp2pq, const VectorXf &LDsamplVar,
                                            const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare, VectorXf &snpStore,
                                            const float varg, const float ps, const float overdispersion,
                                            const bool originalModel){
    // -----------------------------------------
    // Initialise the parameters in MCMC sampler
    // -----------------------------------------
    static unsigned iter = 0;                                                                                                                                                       
    long numChr = chromInfoVec.size();                                                                                                                                              
    
    float ssq[numChr], nnz[numChr], s2pq[numChr];                                                                                                                                   
    memset(ssq,0,sizeof(float)*numChr);                                                                                                                                             
    memset(nnz,0,sizeof(float)*numChr);                                                                                                                                             
    memset(s2pq,0,sizeof(float)*numChr);                                                                                                                                            

    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking precents the writing from other threads                     
  
    vector<float> urnd(size), nrnd(size);                                                                                                                                           
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work                                                                                                               
        urnd[i] = Stat::ranf();                                                                                                                                                     
        nrnd[i] = Stat::snorm();                                                                                                                                                    
    }

    // R specific parameters
    int ndist;
    VectorXf gp;
    snpStore.setZero(pis.size());
    // --------------------------------------------------------------------------------
    // Scale the variances in each of the normal distributions by the genetic variance
    // and initialise the class membership probabilities
    // --------------------------------------------------------------------------------
    ndist = pis.size();
    if (originalModel)
        gp = gamma * 0.01 * varg;
    else
        gp = gamma * sigmaSq;
    snpset.resize(ndist);
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
    }
    // --------------------------------------------------------------------------------
    // Cycle over all variants in the window and sample the genetics effects
    // --------------------------------------------------------------------------------

#pragma omp parallel for
    for (unsigned chr=0; chr<numChr; ++chr) 
    {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned windEnd, j;
        float oldSample, varei;
        double rhs, invLhs, uhat;
        
        int indistflag;
        double v1,  b_ls, ssculm, r;
        VectorXf ll, pll, snpindist, var_b_ls;
        ll.setZero(pis.size());
        pll.setZero(pis.size());

        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            oldSample = valuesPtr[i];
            // ---------------------------------------------
            // Calculate residual variance including a 
            // correction for the sampling variation and
            // LD ignored
            // ---------------------------------------------
            varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
            // ------------------------------
            // Derived Bayes R implementation
            // ------------------------------
            // ----------------------------------------------------
            // Add back the content for the corrected rhs for SNP k
            // ----------------------------------------------------
            rhs = rcorr[i] + ZPZdiag[i] * oldSample;
            // ------------------------------------------------------
            // Calculate the beta least squares updates and variances
            // ------------------------------------------------------
            b_ls = rhs / ZPZdiag[i];
            var_b_ls = gp.array() + varei / ZPZdiag[i];
            // ------------------------------------------------------
            // Calculate the likelihoods for each distribution
            // ------------------------------------------------------
            // ll  = (-1.0 / 2.0) * var_b_ls.array().log()  - (b_ls * b_ls)  / (2 * var_b_ls.array());
            ll = (-1.0 / 2.0) * var_b_ls.array().log()  - (b_ls * b_ls)  / (2 * var_b_ls.array()) + pis.array().log();
            // --------------------------------------------------------------
            // Calculate probability that snp is in each of the distributions
            // in this iteration
            // --------------------------------------------------------------
            // pll = (ll.array().exp().cwiseProduct(pis.array())) / ((ll.array().exp()).cwiseProduct(pis.array())).sum();
            for (unsigned k=0; k<pis.size(); ++k) {
              pll[k] = 1.0 / (exp(ll.array() - ll[k])).sum();
            }
            // if (i < 10) {
            //   cout << "P likelihood 1 " << pll << endl;
            //   cout << "P likelihood 2 " << pll2 << endl;
            // }
            // --------------------------------------------------------------
            // Sample the group based on the calculated probabilities
            // --------------------------------------------------------------
            ssculm = 0.0;
            r = urnd[i];
            indistflag = 1;
            for (int kk = 0; kk < ndist; kk++)
            {
                ssculm += pll(kk);
                if (r < ssculm)
                {
                    indistflag = kk + 1;
                    snpStore(kk) = snpStore(kk) + 1; 
                    break;
                }
            }
            snpset[indistflag-1].push_back(i);
            // --------------------------------------------------------------
            // Sample the effect given the group and adjust the rhs
            // --------------------------------------------------------------
            if (indistflag != 1)                                                                                                                                                    
            {                                                                                                                                                                       
                v1 = ZPZdiag[i] + varei / gp((indistflag - 1));                                                                                                                     
//                valuesPtr[i] = normal.sample(rhs / v1, varei / v1);                                                                                                                 
                valuesPtr[i] = rhs / v1 + nrnd[i]*sqrtf(varei / v1);
                rcorr.segment(windStart[i], windSize[i]) += ZPZ[i] * (oldSample - valuesPtr[i]);
                ssq[chr] += (valuesPtr[i] * valuesPtr[i]) / gamma[indistflag - 1];                                                                                                  
                s2pq[chr] += snp2pq[i];                                                                                                                                             
                ++nnz[chr];                                                                                                                                                         
            } else {                                                                                                                                                                
                if (oldSample) rcorr.segment(windStart[i], windSize[i]) += ZPZ[i] * oldSample;                                                                                      
                valuesPtr[i] = 0.0;                                                                                                                                                 
            }  
        }
    }
    // ---------------------------------------------------------------------                                                                                                        
    // Tally up the effect sum of squares and the number of non-zero effects                                                                                                        
    // ---------------------------------------------------------------------                                                                                                        
    sumSq = 0.0;                                                                                                                                                                    
    sum2pq = 0.0;                                                                                                                                                                   
    numNonZeros = 0.0;                                                                                                                                                              
    nnzPerChr.setZero(numChr);                                                                                                                                                      
    for (unsigned i=0; i<numChr; ++i) {                                                                                                                                             
        sumSq += ssq[i];                                                                                                                                                            
        sum2pq += s2pq[i];                                                                                                                                                          
        numNonZeros += nnz[i];                                                                                                                                                      
        nnzPerChr[i] = nnz[i];                                                                                                                                                      
    }
    ++iter;
    
    values = VectorXf::Map(valuesPtr, size); 
}

void ApproxBayesR::SnpEffects::sampleFromFC(const VectorXf &ZPy, const VectorXf &ZPZdiag, const MatrixXf &Z, const float n_ref, const float n_gwas,
                                            const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare,
                                            VectorXf &snpStore, VectorXf &ghat, const float varg, const bool originalModel) {
        sumSq = 0.0;
        numNonZeros = 0;
            
        ghat.setZero(n_ref);
        float oldSample;
        float my_rhs, rhs;
        // -----------------------------------------
        // Initialise the parameters in MCMC sampler
        // -----------------------------------------
        // ----------------
        // Bayes R specific
        // ----------------
        int ndist, indistflag;
        double v1,  b_ls, ssculm, r;
        VectorXf gp, ll, ll2, pll, snpindist, var_b_ls;
        ndist = pis.size();
        snpStore.setZero(pis.size());
        pll.setZero(pis.size());
        // --------------------------------------------------------------------------------
        // Scale the variances in each of the normal distributions by the genetic variance
        // and initialise the class membership probabilities
        // --------------------------------------------------------------------------------
        if (originalModel)
            gp = gamma * 0.01 * varg;
        else
            gp = gamma * sigmaSq;
    //    cout << varg << " " << gp.transpose() << endl;
        snpset.resize(ndist);
        for (unsigned k=0; k<ndist; ++k) {
            snpset[k].resize(0);
        }
        
        for (unsigned i=0; i<size; ++i) {
            // ------------------------------
            // Derived Bayes R implementation
            // ------------------------------
            // ----------------------------------------------------
            // Add back the content for the corrected rhs for SNP k
            // ----------------------------------------------------
            //my_rhs = Z.col(i).dot(ycorr);
            oldSample = values[i];
            rhs = ZPy[i] - n_gwas/n_ref*Z.col(i).dot(ghat) + ZPZdiag[i]*oldSample;
            // ------------------------------------------------------
            // Calculate the beta least squares updates and variances
            // ------------------------------------------------------
            b_ls = rhs / ZPZdiag[i];
            var_b_ls = gp.array() + vare / ZPZdiag[i];
            // ------------------------------------------------------
            // Calculate the likelihoods for each distribution
            // ------------------------------------------------------
            // ll  = (-1.0 / 2.0) * var_b_ls.array().log()  - (b_ls * b_ls)  / (2 * var_b_ls.array());
            ll = (-1.0 / 2.0) * var_b_ls.array().log()  - (b_ls * b_ls)  / (2 * var_b_ls.array()) + pis.array().log();
            // --------------------------------------------------------------
            // Calculate probability that snp is in each of the distributions
            // in this iteration
            // --------------------------------------------------------------
            // pll = (ll.array().exp().cwiseProduct(pis.array())) / ((ll.array().exp()).cwiseProduct(pis.array())).sum();
            for (unsigned k=0; k<pis.size(); ++k) {
                pll[k] = 1.0 / (exp(ll.array() - ll[k])).sum();
            }
            // --------------------------------------------------------------
            // Sample the group based on the calculated probabilities
            // --------------------------------------------------------------
            ssculm = 0.0;
            r = Stat::ranf();
            indistflag = 1;
            for (int kk = 0; kk < ndist; kk++)
            {
                ssculm += pll(kk);
                if (r < ssculm)
                {
                    indistflag = kk + 1;
                    snpStore(kk) = snpStore(kk) + 1;
                    break;
                }
            }
            snpset[indistflag-1].push_back(i);
            // --------------------------------------------------------------
            // Sample the effect given the group and adjust the rhs
            // --------------------------------------------------------------
            if (indistflag != 1)
            {
                v1 = ZPZdiag[i] + vare / gp((indistflag - 1));
                values[i] = normal.sample(rhs / v1, vare / v1);
                ghat  += Z.col(i) * (values[i] - oldSample);
                sumSq += (values[i] * values[i]) / gamma[indistflag - 1];
                ++numNonZeros;
            } else {
                if (oldSample) ghat -= Z.col(i) * oldSample;
                values[i] = 0.0;
            }
        }
}



void ApproxBayesR::SnpEffects::adjustByCG(const VectorXf &ZPy, const vector<SparseVector<float> > &ZPZsp, VectorXf &rcorr) {
    // construct mixed model equations for those SNPs with nonzero effects and solve the equations using conjugate gradient method
    // then adjust the Gibbs samples with the CG solutions
    
    VectorXf ZPyNZ(numNonZeros);

    vector<Triplet<float> > tripletList;
    tripletList.reserve(numNonZeros);
    
    for (unsigned i=0; i<numNonZeros; ++i) {
        unsigned row = deltaNzIdx[i];
        VectorXf val;
        val.setZero(size);
        for (SparseVector<float>::InnerIterator it(ZPZsp[row]); it; ++it) {
            val[it.index()] = it.value();
        }
        val[row] += lambdaVec[row];
//        cout << "val " << val.transpose() << endl;
        for (unsigned j=0; j<numNonZeros; ++j) {
            unsigned col = deltaNzIdx[j];
//            cout << i << " " << j << " " << row << " " << col << endl;
            tripletList.push_back(Triplet<float>(i, j, val[col]));
//            cout << i << " " << j << " " << val[col] << endl;
        }
        ZPyNZ[i] = ZPy[row];
    }
    
    SpMat C(numNonZeros, numNonZeros);
    C.setFromTriplets(tripletList.begin(), tripletList.end());
    C.makeCompressed();
    tripletList.clear();

//    cout << "C \n" << C.block(0,0,10,10) << endl;
    
    SimplicialLLT<SpMat> solverC;
    solverC.compute(C);
    
    if(solverC.info()!=Success) {
        cout << "Oh: Very bad" << endl;
    }
    
    SpMat eye(numNonZeros, numNonZeros);
    eye.setIdentity();
    
    SpMat Cinv = solverC.solve(eye);
    
    LLT<MatrixXf> llt;
    llt.compute(Cinv); // cholesky decomposition
    VectorXf nrnd(numNonZeros);
    for (unsigned i=0; i<numNonZeros; ++i) {
        nrnd[i] = Stat::snorm();
    }

//    ConjugateGradient<SpMat, Lower|Upper> cg;
//    cg.compute(C);
    VectorXf sol(numNonZeros);
//    sol = cg.solve(ZPyNZ + llt.matrixL()*nrnd);
    
    sol = Cinv * ZPyNZ + llt.matrixL()*nrnd;

//    cout << "numNonZeros " << numNonZeros << endl;
//    cout << "size C " << C.size() << endl;
////
//    cout << "ZPyNZ " << ZPyNZ << endl;
//        cout << "sol " << sol << endl;
//    cout << "uhatVec " << uhatVec << endl;
//
//    cout << "#nonZero:        " << numNonZeros << endl;
//    cout << "#iterations:     " << cg.iterations() << endl;
//    cout << "estimated error: " << cg.error()      << endl;

    float oldSample;
    for (unsigned i=0, j=0; i<size; ++i) {
        if (deltaNZ[i]) {
            oldSample = values[i];
            values[i] = sol[j];
            //values[i] *= sol[j]/uhatVec[i];
//            cout << i << " " << sol[j]/uhatVec[i] << endl;
            for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                rcorr[it.index()] += it.value() * (oldSample - values[i]);
            }
            ++j;
        }
    }
    
//    cout << "old sumsq " << sumSq << endl;
    sumSq = values.cwiseProduct(invGammaVec).dot(values);
//    cout << "new sumsq " << sumSq << endl;

}

void ApproxBayesR::SnpEffects::sampleFromFC(const VectorXf &ZPy, const SpMat &ZPZsp, const VectorXf &ZPZdiag,
                                            VectorXf &rcorr, const VectorXf &LDsamplVar,
                                            const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, VectorXf &snpStore,
                                            const float varg, const float vare, const float ps, const float overdispersion, const bool originalModel) {
    // CG-accelerated Gibbs sampling algorithm
    // first sample delta conditional on beta for all SNPs
    // then construct mixed model equations for which the solutions are samples from the Gibbs sampling
    // and solve the equations by conjugate gradient method
    
    VectorXf lambdaVec(size);
    VectorXf invGammaVec(size);
    
    unsigned ndist = gamma.size();
    snpStore.setZero(ndist);
    
    float varei;
    float rhs;
    
    ArrayXf wtdSigmaSq(ndist);
    ArrayXf invWtdSigmaSq(ndist);
    ArrayXf logWtdSigmaSq(ndist);
    ArrayXf logPis = pis.array().log();
    ArrayXf invLhs(ndist);
    ArrayXf uhat(ndist);
    ArrayXf logDelta(ndist);
    ArrayXf probDelta(ndist);
    
    unsigned delta;
    
    if (originalModel) {
        wtdSigmaSq = gamma * 0.01 * varg;
    } else {
        wtdSigmaSq = gamma * sigmaSq;
    }
    
    invWtdSigmaSq = wtdSigmaSq.inverse();
    logWtdSigmaSq = wtdSigmaSq.log();
    
    VectorXf invGamma = gamma.inverse();
    invGamma[0] = 0;


    for (unsigned i=0; i<size; ++i) {
        
        varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
        
        rhs  = rcorr[i] + ZPZdiag[i]*values[i];
        
        invLhs = (ZPZdiag[i] + varei*invWtdSigmaSq).inverse();
        uhat = invLhs*rhs;
        
        logDelta = 0.5*(invLhs.log() - logWtdSigmaSq + uhat*rhs) + logPis;
        logDelta[0] = logPis[0];
        
        for (unsigned k=0; k<ndist; ++k) {
            probDelta[k] = 1.0f/(logDelta-logDelta[k]).exp().sum();
        }
        
        delta = bernoulli.sample(probDelta);
        
        deltaNZ[i] = delta ? 1:0;
        
        snpset[delta].push_back(i);
        snpStore[delta]++;

        lambdaVec[i] = varei*invWtdSigmaSq[delta];
        invGammaVec[i] = invGamma[delta];
    }
    
    numNonZeros = deltaNZ.sum();
    
    VectorXf lambdaNZ(numNonZeros);
    VectorXf RHS(numNonZeros);
    SpMat eye(numNonZeros, numNonZeros);
    vector<Triplet<float> > tripletList;
    tripletList.reserve(numNonZeros);
    for (unsigned i=0, j=0; i<size; ++i) {
        if (deltaNZ[i]) {
            tripletList.push_back(Triplet<float>(i,i,1));
            lambdaNZ[j] = lambdaVec[i];
            RHS[j] = ZPy[i] + normal.sample(0.0, ZPZdiag[i] + lambdaVec[i]);
            ++j;
        }
    }
    eye.setFromTriplets(tripletList.begin(), tripletList.end());
    eye.makeCompressed();
    tripletList.clear();

    SpMat LHS = eye * ZPZsp * eye;
    LHS.diagonal() += lambdaNZ;
    
    ConjugateGradient<SpMat, Lower|Upper> cg;
    cg.compute(LHS);
    values = cg.solve(RHS);
    
    sumSq = values.cwiseProduct(invGamma).dot(values);
    
    rcorr = ZPy - ZPZsp * values;
}

void ApproxBayesR::SnpEffects::sampleFromFC(vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
                                            const vector<LDBlockInfo*> keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks,
                                            const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, VectorXf &snpStore, const float varg,
                                            const bool originalModel) {
    // -----------------------------------------
    // This method uses low-rank model with eigen-decomposition of LD matrices
    // -----------------------------------------
    long nBlocks = keptLdBlockInfoVec.size();
    
    whatBlocks.resize(nBlocks);
    ssqBlocks.resize(nBlocks);
    for (unsigned i=0; i<nBlocks; ++i) {
        whatBlocks[i].resize(wcorrBlocks[i].size());
    }

    float ssq[nBlocks], s2pq[nBlocks], nnz[nBlocks];
    memset(ssq,0, sizeof(float)*nBlocks);
    memset(s2pq,0,sizeof(float)*nBlocks);
    memset(nnz,0, sizeof(float)*nBlocks);

    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
    // R specific parameters
    int ndist = pis.size();
    ArrayXf logPis = pis.array().log();
    ArrayXf wtdSigmaSq(ndist);
    ArrayXf invWtdSigmaSq(ndist);
    ArrayXf logWtdSigmaSq(ndist);

    if (originalModel) {
        wtdSigmaSq = gamma * 0.01 * varg;
    } else {
        wtdSigmaSq = gamma * sigmaSq;
    }

    snpStore.setZero(pis.size());
    
    invWtdSigmaSq = wtdSigmaSq.inverse();
    logWtdSigmaSq = wtdSigmaSq.log();

    snpset.resize(ndist);
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
    }
    
    // --------------------------------------------------------------------------------
    // Cycle over all variants in the window and sample the genetics effects
    // --------------------------------------------------------------------------------

    #pragma omp parallel for schedule(dynamic)
    for(unsigned blk = 0; blk < nBlocks; blk++){
        Ref<const MatrixXf> Q = Qblocks[blk];
        Ref<VectorXf> wcorr = wcorrBlocks[blk];
        Ref<VectorXf> what = whatBlocks[blk];

        what.setZero();
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[blk];
        
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd   = blockInfo->endSnpIdx;
        
        float vareDn = nGWASblocks[blk] / vareBlocks[blk];

        ArrayXf invLhs = 1.0/(vareDn + invWtdSigmaSq);
        ArrayXf logInvLhsMsigma = invLhs.log() - logWtdSigmaSq;

        for(unsigned i = blockStart; i <= blockEnd; i++){
            float oldSample = valuesPtr[i];
            Ref<const VectorXf> Qi = Q.col(i - blockStart);
            float rhs = (Qi.dot(wcorr) + oldSample)*vareDn;
            ArrayXf uhat = invLhs * rhs;
            ArrayXf logDelta = 0.5*(logInvLhsMsigma + uhat*rhs) + logPis;
            logDelta[0] = logPis[0];
            
            ArrayXf probDelta(ndist);
            for (unsigned k=0; k<ndist; ++k) {
                probDelta[k] = 1.0f/(logDelta-logDelta[k]).exp().sum();
            }
                        

            unsigned delta;
            #pragma omp critical
            {
                delta = bernoulli.sample(probDelta);

                snpset[delta].push_back(i);
                snpStore[delta]++;
            }
            
            if (delta) {
                valuesPtr[i] = uhat[delta] + nrnd[i]*sqrtf(invLhs[delta]);
                wcorr += Qi*(oldSample - valuesPtr[i]);
                what  += Qi* valuesPtr[i];
                ssq[blk] += (valuesPtr[i] * valuesPtr[i]) / gamma[delta];
                ++nnz[blk];
            }
            else {
                if (oldSample) wcorr += Qi * oldSample;
                valuesPtr[i] = 0.0;
            }
        }

    }

    // ---------------------------------------------------------------------
    // Tally up the effect sum of squares and the number of non-zero effects
    // ---------------------------------------------------------------------
    sumSq = 0.0;
    numNonZeros = 0;
    nnzPerBlk.setZero(nBlocks);
    for (unsigned blk=0; blk<nBlocks; ++blk) {
        sumSq += ssq[blk];
        numNonZeros += nnz[blk];
        nnzPerBlk[blk] = nnz[blk];
        ssqBlocks[blk] = ssq[blk];
    }
    values = VectorXf::Map(valuesPtr, size);
 
}


void ApproxBayesR::BlockGenotypicVar::compute(const vector<VectorXf> &whatBlocks){
    for (unsigned i=0; i<numBlocks; ++i) {
        values[i] = whatBlocks[i].squaredNorm();
        //cout << "varg " << i << " " << values[i] << endl;
    }
    total = values.sum();
}

void ApproxBayesR::BlockResidualVar::sampleFromFC(vector<VectorXf> &wcorrBlocks, VectorXf &ssqBlocks, const VectorXf &nGWASblocks, const VectorXf &numEigenvalBlock){
    for (unsigned i=0; i<numBlocks; ++i) {
        float sse = wcorrBlocks[i].squaredNorm() * nGWASblocks[i];
        float dfTilde = df + numEigenvalBlock[i];
        float scaleTilde = sse + df*scale;
        float sample = InvChiSq::sample(dfTilde, scaleTilde);
        if (ssqBlocks[i]/sample > threshold) {
            values[i] = sample;
        } else {
            values[i] = vary;
        }
        //cout << "vare " << i << " " << values[i] << endl;
    }
    mean = values.mean();
}




// *******************************************************
// Approximate Bayes RS
// *******************************************************

void ApproxBayesRS::SnpEffects::sampleFromFC(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                             const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo *> &chromInfoVec,
                                             const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare,
                                             const ArrayXf &snp2pqPowS, const VectorXf &snp2pq,
                                             const VectorXf &LDsamplVar, const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const float varg,
                                             const float ps, const float overdispersion, const bool originalModel) {
    // sample SNP effects with a sparse LD matrix
    
    long numChr = chromInfoVec.size();
    
    float ssq[numChr], s2pq[numChr], nnz[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(s2pq,0,sizeof(float)*numChr);
    memset(nnz,0, sizeof(float)*numChr);
    
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
    ArrayXf wtdSigmaSq(ndist);
    ArrayXf invWtdSigmaSq(ndist);
    ArrayXf logWtdSigmaSq(ndist);
    ArrayXf logPis = pis.array().log();
    ArrayXf log2pqPowS = snp2pqPowS.log();
    
    if (originalModel) {
        wtdSigmaSq = gamma * 0.01 * varg;
    } else {
        wtdSigmaSq = gamma * sigmaSq;
    }
    
    invWtdSigmaSq = wtdSigmaSq.inverse();
    logWtdSigmaSq = wtdSigmaSq.log();
    
    numSnpMix.setZero(ndist);
    snpset.resize(ndist);
    
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
    }
    
    //#pragma omp parallel for  // openmp is not working for SBayesR
    for (unsigned chr=0; chr<numChr; ++chr) {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;

        float oldSample;
        float sampleDiff;
        float rhs;
        float varei;
        
        ArrayXf invLhs(ndist);
        ArrayXf uhat(ndist);
        ArrayXf logDelta(ndist);
        ArrayXf probDelta(ndist);
        
        unsigned delta;
        
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            oldSample = valuesPtr[i];

            varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
            
            rhs  = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            
            invLhs = (ZPZdiag[i]/varei + invWtdSigmaSq/snp2pqPowS[i]).inverse();
            uhat = invLhs*rhs;
            
            logDelta = 0.5*(invLhs.log() - log2pqPowS[i] - logWtdSigmaSq + uhat*rhs) + logPis;
            logDelta[0] = logPis[0];
            
            for (unsigned k=0; k<ndist; ++k) {
                probDelta[k] = 1.0f/(logDelta-logDelta[k]).exp().sum();
            }
            
            delta = bernoulli.sample(probDelta);
            
            snpset[delta].push_back(i);
            numSnpMix[delta]++;
            
            if (delta) {
                valuesPtr[i] = uhat[delta] + nrnd[i]*sqrtf(invLhs[delta]);
                sampleDiff = oldSample - valuesPtr[i];
                for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                    rcorr[it.index()] += it.value() * sampleDiff;
                }
                ssq[chr] += (valuesPtr[i] * valuesPtr[i]) / (gamma[delta]*snp2pqPowS[i]);
                ++nnz[chr];
            }
            else {
                if (oldSample) {
                    for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                        rcorr[it.index()] += it.value() * oldSample;
                    }
                }
                valuesPtr[i] = 0.0;
            }
        }
    }
    
    wtdSumSq = 0.0;
    numNonZeros = 0.0;
    nnzPerChr.setZero(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        wtdSumSq += ssq[i];
        numNonZeros += nnz[i];
        nnzPerChr[i] = nnz[i];
    }
    values = VectorXf::Map(valuesPtr, size);
}

void ApproxBayesRS::SnpEffects::sampleFromFC(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                             const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo *> &chromInfoVec,
                                             const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare,
                                             const ArrayXf &snp2pqPowS, const VectorXf &snp2pq,
                                             const VectorXf &LDsamplVar, const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const float varg,
                                             const float ps, const float overdispersion, const bool originalModel) {
    // sample SNP effects with a full LD matrix
    
    long numChr = chromInfoVec.size();
    
    float ssq[numChr], s2pq[numChr], nnz[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(s2pq,0,sizeof(float)*numChr);
    memset(nnz,0, sizeof(float)*numChr);
    
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
    ArrayXf wtdSigmaSq(ndist);
    ArrayXf invWtdSigmaSq(ndist);
    ArrayXf logWtdSigmaSq(ndist);
    ArrayXf logPis = pis.array().log();
    ArrayXf log2pqPowS = snp2pqPowS.log();
    
    if (originalModel) {
        wtdSigmaSq = gamma * 0.01 * varg;
    } else {
        wtdSigmaSq = gamma * sigmaSq;
    }
    
    invWtdSigmaSq = wtdSigmaSq.inverse();
    logWtdSigmaSq = wtdSigmaSq.log();
    
    numSnpMix.setZero(ndist);
    snpset.resize(ndist);
    
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
    }
    
    //#pragma omp parallel for  // openmp is not working for SBayesR
    for (unsigned chr=0; chr<numChr; ++chr) {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;

        float oldSample;
        float sampleDiff;
        float rhs;
        float varei;
        
        ArrayXf invLhs(ndist);
        ArrayXf uhat(ndist);
        ArrayXf logDelta(ndist);
        ArrayXf probDelta(ndist);
        
        unsigned delta;
        
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            oldSample = valuesPtr[i];

            varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
            
            rhs  = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            
            invLhs = (ZPZdiag[i]/varei + invWtdSigmaSq/snp2pqPowS[i]).inverse();
            uhat = invLhs*rhs;
            
            logDelta = 0.5*(invLhs.log() - log2pqPowS[i] - logWtdSigmaSq + uhat*rhs) + logPis;
            logDelta[0] = logPis[0];
            
            for (unsigned k=0; k<ndist; ++k) {
                probDelta[k] = 1.0f/(logDelta-logDelta[k]).exp().sum();
            }
            
            delta = bernoulli.sample(probDelta);
            
            snpset[delta].push_back(i);
            numSnpMix[delta]++;
            
            if (delta) {
                valuesPtr[i] = uhat[delta] + nrnd[i]*sqrtf(invLhs[delta]);
                rcorr.segment(windStart[i], windSize[i]) += ZPZ[i] * (oldSample - valuesPtr[i]);
                ssq[chr] += (valuesPtr[i] * valuesPtr[i]) / (gamma[delta]*snp2pqPowS[i]);
                ++nnz[chr];
            }
            else {
                if (oldSample) {
                    if (oldSample) rcorr.segment(windStart[i], windSize[i]) += ZPZ[i] * oldSample;
                }
                valuesPtr[i] = 0.0;
            }
        }
    }
    
    wtdSumSq = 0.0;
    numNonZeros = 0.0;
    nnzPerChr.setZero(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        wtdSumSq += ssq[i];
        numNonZeros += nnz[i];
        nnzPerChr[i] = nnz[i];
    }
    values = VectorXf::Map(valuesPtr, size);
}

void ApproxBayesRS::Sp::sampleFromFC(vector<vector<unsigned> > &snpset, const VectorXf &snpEffects,
                                     float &sigmaSq, const VectorXf &gamma,
                                     const VectorXf &snp2pq, ArrayXf &snp2pqPowS, const ArrayXf &logSnp2pq,
                                     const float vg, float &scale, float &sum2pqSplusOne, const bool originalModel) {
    // Hamiltonian Monte Carlo
    // note that the scale factor of sigmaSq will be simultaneously updated
    
    unsigned nnzMix = snpset.size() - 1; // nonzero component
    
    // Prepare
    vector<ArrayXf> snpEffectMix(nnzMix);
    vector<ArrayXf> snp2pqMix(nnzMix);
    vector<ArrayXf> logSnp2pqMix(nnzMix);
    
    float snp2pqLogSumNZ = 0.0;
    
    for (unsigned i=0; i<nnzMix; ++i) {
        unsigned k=i+1;
        long isize = snpset[k].size();
        snpEffectMix[i].resize(isize);
        snp2pqMix[i].resize(isize);
        logSnp2pqMix[i].resize(isize);
        for (unsigned j=0; j<isize; ++j) {
            snpEffectMix[i][j] = snpEffects[snpset[k][j]];
            snp2pqMix[i][j] = snp2pq[snpset[k][j]];
            logSnp2pqMix[i][j] = logSnp2pq[snpset[k][j]];
        }
        snp2pqLogSumNZ += logSnp2pqMix[i].sum();
    }
    
    float curr = value;
    float curr_p = Stat::snorm();
    
    float cand = curr;
    // Make a half step for momentum at the beginning
    float cand_p = curr_p - 0.5*stepSize * gradientU(curr, nnzMix, snpEffectMix, snp2pqLogSumNZ, snp2pqMix, logSnp2pqMix, sigmaSq, gamma, vg);

    for (unsigned i=0; i<numSteps; ++i) {
        // Make a full step for the position
        cand += stepSize * cand_p;
        if (i < numSteps-1) {
            // Make a full step for the momentum, except at end of trajectory
            cand_p -= stepSize * gradientU(cand, nnzMix, snpEffectMix, snp2pqLogSumNZ, snp2pqMix, logSnp2pqMix, sigmaSq, gamma, vg);
        } else {
            // Make a half step for momentum at the end
            cand_p -= 0.5*stepSize * gradientU(cand, nnzMix, snpEffectMix, snp2pqLogSumNZ, snp2pqMix, logSnp2pqMix, sigmaSq, gamma, vg);
        }
        //cout << i << " " << cand << endl;
    }

    // Evaluate potential (negative log posterior) and kinetic energies at start and end of trajectory
    float scaleCurr, scaleCand;
    float curr_H = computeU(curr, nnzMix, snpEffectMix, snp2pqLogSumNZ, snp2pqMix, logSnp2pqMix, sigmaSq, gamma, vg, scaleCurr, originalModel) + 0.5*curr_p*curr_p;
    float cand_H = computeU(cand, nnzMix, snpEffectMix, snp2pqLogSumNZ, snp2pqMix, logSnp2pqMix, sigmaSq, gamma, vg, scaleCand, originalModel) + 0.5*cand_p*cand_p;
    
    if (Stat::ranf() < exp(curr_H-cand_H)) {  // accept
        value = cand;
        scale = scaleCand;
        snp2pqPowS = snp2pq.array().pow(cand);
        sum2pqSplusOne = 0.0;
        for (unsigned i=0; i<nnzMix; ++i) sum2pqSplusOne += snp2pqMix[i].pow(1.0+value).sum();
        ar.count(1, 0.5, 0.9);
    } else {
        ar.count(0, 0.5, 0.9);
    }
    
    if (!(ar.cnt % 10)) {
        if      (ar.value < 0.6) stepSize *= 0.8;
        else if (ar.value > 0.8) stepSize *= 1.2;
    }
    
    if (ar.consecRej > 20) stepSize *= 0.8;
    
    tuner.value = stepSize;
}

float ApproxBayesRS::Sp::gradientU(const float S, const unsigned nnzMix, const vector<ArrayXf> &snpEffectMix, const float snp2pqLogSum, const vector<ArrayXf> &snp2pqMix, const vector<ArrayXf> &logSnp2pqMix, const float sigmaSq, const VectorXf &gamma, const float vg){
    float constantA = snp2pqLogSum;
    ArrayXf constantB(nnzMix);
    for (unsigned i=0; i<nnzMix; ++i) {
        constantB[i] = (snpEffectMix[i].square()*logSnp2pqMix[i]/snp2pqMix[i].pow(S)).sum()/gamma[i+1];
    }
    return 0.5*constantA - 0.5/sigmaSq*constantB.sum() + S/var;
}

float ApproxBayesRS::Sp::computeU(const float S, const unsigned nnzMix, const vector<ArrayXf> &snpEffectMix, const float snp2pqLogSum, const vector<ArrayXf> &snp2pqMix, const vector<ArrayXf> &logSnp2pqMix, const float sigmaSq, const VectorXf &gamma, const float vg, float &scale, const bool originalModel) {
    vector<ArrayXf> snp2pqPowSMix(nnzMix);
    float constantA = snp2pqLogSum;
    ArrayXf constantB(nnzMix);
    ArrayXf constantC(nnzMix);
    for (unsigned i=0; i<nnzMix; ++i) {
        snp2pqPowSMix[i] = snp2pqMix[i].pow(S);
        constantB[i] = (snpEffectMix[i].square()/snp2pqPowSMix[i]).sum()/gamma[i+1];
        if (originalModel) constantC[i] = snp2pqPowSMix[i].sum()*gamma[i+1];
        else constantC[i] = (snp2pqMix[i]*snp2pqPowSMix[i]).sum()*gamma[i+1];
    }
    scale = 0.5*vg/constantC.sum();
    return 0.5*S*constantA + 0.5/sigmaSq*constantB.sum() + 0.5*S*S/var;
}

void ApproxBayesRS::sampleUnknowns() {
    static int iter = 0;
    unsigned cnt=0;
    do {
        if (sparse)
            snpEffects.sampleFromFC(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, sigmaSq.value, Pis.values, gamma.values, vare.value, snp2pqPowS, data.snp2pq, data.LDsamplVar, data.se, data.tss, varei, data.n, varg.value, ps.value, overdispersion, originalModel);
        else
            snpEffects.sampleFromFC(rcorr, data.ZPZ, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, sigmaSq.value, Pis.values, gamma.values, vare.value, snp2pqPowS, data.snp2pq, data.LDsamplVar, data.se, data.tss, varei, data.n, varg.value, ps.value, overdispersion, originalModel);
        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros == 0);
    if (diagnose) nro.compute(rcorr, data.ZPZdiag, data.LDsamplVar, varg.value, vare.value, snpEffects.header, snpEffects.leaveout, data.ZPZsp, data.ZPy, snpEffects.values);
    sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
    if (estimatePi) Pis.sampleFromFC(snpEffects.numSnpMix);
    numSnps.getValues(snpEffects.numSnpMix);
    nnzSnp.getValue(snpEffects.numNonZeros);
    sigmaSqG.compute(sigmaSq.value, snpEffects.sum2pq);
    
    covg.compute(data.ypy, snpEffects.values, data.ZPy, rcorr);
    varg.compute(snpEffects.values, data.ZPy, rcorr, covg.value);
    vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, covg.value);
    
    hsq.compute(varg.value, vare.value);
    
    S.sampleFromFC(snpEffects.snpset, snpEffects.values, sigmaSq.value, gamma.values, data.snp2pq, snp2pqPowS, logSnp2pq, genVarPrior, sigmaSq.scale, snpEffects.sum2pqSplusOne, originalModel);

    if (iter >= 2000) sigmaSq.scale = scalePrior;
    scale.getValue(sigmaSq.scale);
    // cout << "iter " << iter << " scalePrior " << scalePrior << "sigmaSq.scale " << sigmaSq.scale << endl;
    
    if (sparse)
    rounding.computeRcorr(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    else
    rounding.computeRcorr(data.ZPy, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    if (modelPS) ps.compute(rcorr, data.ZPZdiag, data.LDsamplVar, varg.value, vare.value, data.chisq);
    
    nnzSnp.getValue(snpEffects.numNonZeros);
    sigmaSqG.compute(sigmaSq.value, snpEffects.sum2pq);
    
    //    numSnpVg.compute(snpEffects.values, data.ZPZdiag, varg.value, vare.nobs);
    if (originalModel) {
        Vgs.compute(snpEffects.values, data.ZPy, rcorr, snpEffects.snpset, varg.value, vare.nobs);
//        if (sparse)
//        Vgs.compute(snpEffects.values, data.ZPZsp, snpEffects.snpset, varg.value, vare.nobs);
//        else
//        Vgs.compute(snpEffects.values, data.ZPZ, snpEffects.snpset, varg.value, vare.nobs);
    }
    
    if (++iter < 2000) {
        genVarPrior += (varg.value - genVarPrior)/iter;
        scalePrior  += (sigmaSq.scale - scalePrior)/iter;
    }

}


// *******************************************************
// Bayes Kappa - Approximate
// *******************************************************

void ApproxBayesKappa::AcceptanceRate::count(const bool state, const float lower, const float upper){
    accepted += state;
    value = accepted/float(++cnt);
    if (!state) ++consecRej;
    else consecRej = 0;
}

void ApproxBayesKappa::Kappa::randomWalkMHsampler(const float sigmaSq, const VectorXf &snpEffects, const VectorXf &snpindist){
    // Random walk Mentroplis-Hastings
    // note that the scale factor of sigmaSq will be simultaneously updated
    // cout << "snpindist " << snpindist.segment(1,100) << endl;
    // cout << "snpEffects " << snpEffects.segment(1,10) << endl;
    // cout << "sigmaSq   " << sigmaSq << endl;
    boost::math::normal_distribution <> d(0 ,1);
    float curr = value;
    float cand = 0.0;
    do {
      cand = sample(value, varProp);
      //cout << cdf(d, 1) << endl;
      //cout << cand << endl;
    } while(!(cand >= 0)); // This part implements in the truncated sampler as the value can be negative, which does make any sense for kappa. Taken from # darrenjw.wordpress.com
    //cout << "cand at top " << cand << endl;
    float gamKap = boost::math::gamma_p_derivative(k0, curr / theta0) / theta0;
    float betasSqrCurr = 0;
    float betasSqrCand = 0;
    const double PiVal = boost::math::constants::pi<double>();
    //cout << "gamKap " << gamKap << " Pi " << PiVal << endl; 
    int numSnps = snpEffects.size();
    VectorXf scalesCurr;
    VectorXf scalesCand;
    scalesCurr = -curr * (snpindist.array() - 1);
    scalesCand = -cand * (snpindist.array() - 1);
    scalesCurr = scalesCurr.array().exp();
    scalesCand = scalesCand.array().exp();
    //cout << "Hello from Random Walk" << " current value is " << value << endl;  // " scalesCurr " << scalesCurr << "  scalesCand " <<  scalesCand << endl;
    float detValCurr = 0;
    float detValCand = 0;
    for (unsigned i=0; i<numSnps; ++i) {
        if (snpEffects[i]) {
          detValCurr += -0.5f * logf(2.0f * PiVal * sigmaSq * scalesCurr[i]);
          detValCand += -0.5f * logf(2.0f * PiVal * sigmaSq * scalesCand[i]); 
          betasSqrCurr += snpEffects[i] * snpEffects[i] / scalesCurr[i];
          betasSqrCand += snpEffects[i] * snpEffects[i] / scalesCand[i];
        }
    }
  
    float logCurr = detValCurr - (1 / (2.0f * sigmaSq)) * betasSqrCurr + logf(gamKap);
    float logCand = detValCand - (1 / (2.0f * sigmaSq)) * betasSqrCand + logf(gamKap);
    
    //cout << "curr " << curr << " logCurr " << logCurr << " cand " << cand << " logCand " << logCand << " sigmaSq " << sigmaSq << endl;
    float A;
    A = exp(logCand-logCurr)*cdf(d, curr)/cdf(d, cand); // This accounts for the non-symmetry of the truncated sampler
    //cout << " exp(logCand-logCurr) " << exp(logCand-logCurr) << " cdf(d, curr)/cdf(d, cand) " << cdf(d, curr)/cdf(d, cand) << "A " << A << endl;
    if (Stat::ranf() < A) {  // accept
        value = cand;
        // cout << "value " << value << endl;
        ar.count(1, 0.1, 0.5);
    } else {
        ar.count(0, 0.1, 0.5);
    }
    float delta = min(exp(0.05), exp(1/sqrt(ar.cnt)));
    
    if (!(ar.cnt % 50)) {
        if      (ar.value < 0.44) varProp /= delta;
        else if (ar.value >= 0.44) varProp *= delta;
         cout << "ar.count " << ar.cnt  << " Delta " << delta << " Acc rate " << ar.value <<  endl;
    }
    
    tuner.value = varProp;
}


void ApproxBayesKappa::sampleUnknowns(){
    static int iter = 0;
//    fixedEffects.sampleFromFC(data.XPX, data.XPXdiag, data.ZPX, data.XPy, snpEffects.values, vare.value, rcorr);
    unsigned cnt=0;
    do {
        if (sparse)
            snpEffects.sampleFromFC(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, data.se, data.tss, varei, data.n, data.snp2pq, sigmaSq.value, Pis.values, gamma.values, vare.value, snpStore, kappa.value, snpindist.values);
        else
            snpEffects.sampleFromFC(rcorr, data.ZPZ, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, data.se, data.tss, varei, data.n, data.snp2pq, sigmaSq.value, Pis.values, gamma.values, vare.value, snpStore, kappa.value, snpindist.values);
        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros == 0);
    sigmaSq.sampleFromFC(snpEffects.sumSq, snpEffects.numNonZeros);
    kappa.randomWalkMHsampler(sigmaSq.value, snpEffects.values, snpindist.values);
    Pis.sampleFromFC(snpStore);
    nnzSnp.getValue(snpEffects.numNonZeros);
    sigmaSqG.compute(sigmaSq.value, snpEffects.sum2pq);
    varg.compute(snpEffects.values, data.ZPy, rcorr, 0);
    vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, 0);
    hsq.compute(varg.value, vare.value);

    if (iter >= 2000) sigmaSq.scale = scalePrior;
    scale.getValue(sigmaSq.scale);
    // cout << "iter " << iter << " scalePrior " << scalePrior << "sigmaSq.scale " << sigmaSq.scale << endl;

    if (sparse)
        rounding.computeRcorr(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    else
        rounding.computeRcorr(data.ZPy, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    nnzSnp.getValue(snpEffects.numNonZeros);
    sigmaSqG.compute(sigmaSq.value, snpEffects.sum2pq);

    if (++iter < 2000) {
        if (noscale)
        {
            scalePrior  = 0.5f * varg.value / (data.snp2pq.array().sum()*(1-Pis.values[(Pis.values.size()-1)]));
        } else
        {
            scalePrior  = 0.5f * varg.value / (data.snp2pq.size()*(1-Pis.values[(Pis.values.size()-1)]));
        }
        genVarPrior += (varg.value - genVarPrior)/iter;
        scalePrior  += (sigmaSq.scale - scalePrior)/iter;
    }
}

// ==============================================================
// Kappa Sparse vector version
// ==============================================================

void ApproxBayesKappa::SnpEffects::sampleFromFC(VectorXf &rcorr, const vector<SparseVector<float>> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const VectorXf &snp2pq,
                                            const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare, VectorXf &snpStore, 
                                            const float kappa, VectorXf &snpindist){
    // -----------------------------------------
    // Initialise the parameters in MCMC sampler
    // -----------------------------------------
    static unsigned iter = 0;
    long numChr = chromInfoVec.size();

    float ssq[numChr], s2pq[numChr], nnz[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(s2pq,0,sizeof(float)*numChr);
    memset(nnz,0, sizeof(float)*numChr);

    for (unsigned chr=0; chr<numChr; ++chr) {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        if (iter==0) {
            cout << "chr " << chr+1 << " start " << chrStart << " end " << chrEnd << endl;
        }
    }
    if (iter==0) cout << endl;

    // ----------------
    // Bayes Kappa specific
    // ----------------
    int ndist, indistflag;
    double rhs, v1,  b_ls, ssculm, r;
    VectorXf gp, ll, gamgam, pll, var_b_ls;
    gamgam.setZero(pis.size());
    snpStore.setZero(pis.size());
    ll.setZero(pis.size());
    pll.setZero(pis.size());
    //snpindist.setZero(tss.size());
    // kappa=2.302585;
    // --------------------------------------------------------------------------------
    // Scale the variances in each of the normal distributions by the genetic variance
    // and initialise the class membership probabilities
    // --------------------------------------------------------------------------------
    ndist = pis.size();
    for (int dstInd=0; dstInd<(ndist-1); ++dstInd)
    {
      gamgam[dstInd] = exp(-kappa * dstInd);
    }
    // cout << "gamgam " << gamgam << endl;
    gp = gamgam * sigmaSq;
    // --------------------------------------------------------------------------------
    // Cycle over all variants in the window and sample the genetics effects
    // --------------------------------------------------------------------------------
#pragma omp parallel for
    for (unsigned chr=0; chr<numChr; ++chr) 
    {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned windEnd, j;

        float oldSample;
        double rhs, invLhs, uhat;
        
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            // ---------------------------------------------
            // Calculate residual variance from local region
            // ---------------------------------------------
            //cout << "Varei i " << varei[i] << endl;
             if (!(iter % 100)) {
                //float varei = (sse[i] - values.segment(windStart[i], windSize[i]).dot(ZPy.segment(windStart[i], windSize[i]) + rcorr.segment(windStart[i], windSize[i])))/n[i];
                windEnd = windStart[i] + windSize[i];
                varei[i] = tss[i];
                //cout << "Varei 100 " << varei[i] << endl;
                for (j=windStart[i]; j<windEnd; ++j) {
                    if (values[j]) varei[i] -= values[j]*(ZPy[j] + rcorr[j]);
                }
                varei[i] /= n[i];
            }  
            // ------------------------------
            // Derived Bayes R implementation
            // ------------------------------
            // ----------------------------------------------------
            // Add back the content for the corrected rhs for SNP k
            // ----------------------------------------------------
            oldSample = values[i];
            rhs = rcorr[i] + ZPZdiag[i] * oldSample;
            // ------------------------------------------------------
            // Calculate the beta least squares updates and variances
            // ------------------------------------------------------
            b_ls = rhs / ZPZdiag[i];
            var_b_ls = gp.array() + varei[i] / ZPZdiag[i];
            // ------------------------------------------------------
            // Calculate the likelihoods for each distribution
            // ------------------------------------------------------
            ll = (-1.0 / 2.0) * var_b_ls.array().log()  - (b_ls * b_ls)  / (2 * var_b_ls.array()) + pis.array().log();
            // --------------------------------------------------------------
            // Calculate probability that snp is in each of the distributions
            // in this iteration
            // --------------------------------------------------------------
            // pll = (ll.array().exp().cwiseProduct(pis.array())) / ((ll.array().exp()).cwiseProduct(pis.array())).sum();
            for (unsigned k=0; k<pis.size(); ++k) {
              pll[k] = 1.0 / (exp(ll.array() - ll[k])).sum();
            }
            // --------------------------------------------------------------
            // Sample the group based on the calculated probabilities
            // --------------------------------------------------------------
            ssculm = 0.0;
            r = Stat::ranf();
            indistflag = 1;
            for (int kk = 0; kk < ndist; kk++)
            {
                ssculm += pll(kk);
                if (r < ssculm)
                {
                    indistflag = kk + 1;
                    snpStore(kk) = snpStore(kk) + 1;
                    snpindist[i] = kk + 1; 
                    break;
                }
            }
            // --------------------------------------------------------------
            // Sample the effect given the group and adjust the rhs
            // --------------------------------------------------------------
            if (indistflag != ndist)
            {
                v1 = ZPZdiag[i] + varei[i] / gp((indistflag - 1));
                values[i] = normal.sample(rhs / v1, varei[i] / v1);
                float sampleDiff = oldSample - values[i];
                for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {
                    rcorr[it.index()] += it.value() * sampleDiff;
                }
                ssq[chr] += (values[i] * values[i]) / gamgam[indistflag - 1];
                s2pq[chr] += snp2pq[i];
                ++nnz[chr];
            } else {
                if (oldSample) {
                    for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {
                        rcorr[it.index()] += it.value() * oldSample;
                    }
                }
                values[i] = 0.0;
            }
        }
    }
    // ---------------------------------------------------------------------
    // Tally up the effect sum of squares and the number of non-zero effects
    // ---------------------------------------------------------------------
    sumSq = 0.0;                                                                                                                                                                                
    sum2pq = 0.0;                                                                                                                                                                               
    numNonZeros = 0;                                                                                                                                                                            
    for (unsigned i=0; i<numChr; ++i) {                                                                                                                                                         
        sumSq += ssq[i];                                                                                                                                                                        
        sum2pq += s2pq[i];                                                                                                                                                                      
        numNonZeros += nnz[i];                                                                                                                                                                  
    }                                                                                                                                                                                           
    ++iter; 
}

// ==============================================================
// Kappa Vector of vectors version
// ==============================================================

void ApproxBayesKappa::SnpEffects::sampleFromFC(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const VectorXf &snp2pq,
                                            const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare, VectorXf &snpStore, 
                                            const float kappa, VectorXf &snpindist){
    // -----------------------------------------
    // Initialise the parameters in MCMC sampler
    // -----------------------------------------
    static unsigned iter = 0;
    long numChr = chromInfoVec.size();
    VectorXf ssq, s2pq, nnz;
    ssq.setZero(numChr);
    s2pq.setZero(numChr);
    nnz.setZero(numChr);
    // ----------------
    // Bayes Kappa specific
    // ----------------
    int ndist, indistflag;
    double rhs, v1,  b_ls, ssculm, r;
    VectorXf gp, gamgam, ll, pll, var_b_ls;
    gamgam.setZero(pis.size());
    snpStore.setZero(pis.size());
    ll.setZero(pis.size());
    pll.setZero(pis.size());
    //snpindist.setZero(tss.size());
    // --------------------------------------------------------------------------------
    // Scale the variances in each of the normal distributions by the genetic variance
    // and initialise the class membership probabilities
    // --------------------------------------------------------------------------------
    ndist = pis.size();
    for (int dstInd=0; dstInd<(ndist-1); ++dstInd)
    {
      gamgam[dstInd] = exp(-kappa * dstInd);
    }
    //cout << "Gammas after Kappa " << gamgam << endl;
    //gp = gamma * sigmaSq;
    gp = gamgam * sigmaSq;
    // --------------------------------------------------------------------------------
    // Cycle over all variants in the window and sample the genetics effects
    // --------------------------------------------------------------------------------
    for (unsigned chr=0; chr<numChr; ++chr) 
    {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned windEnd, j;
        float oldSample;
        double rhs, invLhs, uhat;
        
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            // ---------------------------------------------
            // Calculate residual variance from local region
            // ---------------------------------------------
            // cout << "Varei i " << varei[i] << endl;
            if (!(iter % 100)) {
                //float varei = (sse[i] - values.segment(windStart[i], windSize[i]).dot(ZPy.segment(windStart[i], windSize[i]) + rcorr.segment(windStart[i], windSize[i])))/n[i];
                windEnd = windStart[i] + windSize[i];
                varei[i] = tss[i];
                for (j=windStart[i]; j<windEnd; ++j) {
                    if (values[j]) varei[i] -= values[j]*(ZPy[j] + rcorr[j]);
                }
                varei[i] /= n[i];
                // cout << "Varei 1 max" << varei.maxCoeff() << endl;
            }   
            // varei[i] = tss[i] / n[i];
            // ------------------------------
            // Derived Bayes R implementation
            // ------------------------------
            // ----------------------------------------------------
            // Add back the content for the corrected rhs for SNP k
            // ----------------------------------------------------
            oldSample = values[i];
            rhs = rcorr[i] + ZPZdiag[i] * oldSample;
            // ------------------------------------------------------
            // Calculate the beta least squares updates and variances
            // ------------------------------------------------------
            b_ls = rhs / ZPZdiag[i];
            var_b_ls = gp.array() + varei[i] / ZPZdiag[i];
            // ------------------------------------------------------
            // Calculate the likelihoods for each distribution
            // ------------------------------------------------------
            ll = (-1.0 / 2.0) * var_b_ls.array().log()  - (b_ls * b_ls)  / (2 * var_b_ls.array()) + pis.array().log();
            // --------------------------------------------------------------
            // Calculate probability that snp is in each of the distributions
            // in this iteration
            // --------------------------------------------------------------
            // pll = (ll.array().exp().cwiseProduct(pis.array())) / ((ll.array().exp()).cwiseProduct(pis.array())).sum();
            for (unsigned k=0; k<pis.size(); ++k) {
              pll[k] = 1.0 / (exp(ll.array() - ll[k])).sum();
            }
            // cout << "Likelihoods " << pll << endl;
            // --------------------------------------------------------------
            // Sample the group based on the calculated probabilities
            // --------------------------------------------------------------
            ssculm = 0.0;
            r = Stat::ranf();
            indistflag = 1;
            for (int kk = 0; kk < ndist; kk++)
            {
                ssculm += pll(kk);
                if (r < ssculm)
                {
                    indistflag = kk + 1;
                    snpStore(kk) = snpStore(kk) + 1; 
                    snpindist[i] = kk + 1;
                    break;
                }
            }
            // --------------------------------------------------------------
            // Sample the effect given the group and adjust the rhs
            // --------------------------------------------------------------
            //cout << "Indistflag " << indistflag << " ndist " << ndist << endl;
            if (indistflag != ndist)
            {
                v1 = ZPZdiag[i] + varei[i] / gp((indistflag - 1));
                values[i] = normal.sample(rhs / v1, varei[i] / v1);
                rcorr.segment(windStart[i], windSize[i]) += ZPZ[i] * (oldSample - values[i]);
                ssq[chr] += (values[i] * values[i]) / gamgam[indistflag - 1];
                s2pq[chr] += snp2pq[i];
                ++nnz[chr];
            } else {
                if (oldSample) rcorr.segment(windStart[i], windSize[i]) += ZPZ[i] * oldSample;
                values[i] = 0.0;
            }
        }
    }
    // ---------------------------------------------------------------------
    // Tally up the effect sum of squares and the number of non-zero effects
    // ---------------------------------------------------------------------
    sumSq = ssq.sum();
    sum2pq = s2pq.sum();
    numNonZeros = nnz.sum();
    ++iter;
}



void ApproxBayesSMix::SnpEffects::sampleFromFC(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy, const vector<ChromInfo *> &chromInfoVec, const VectorXf &LDsamplVar, const ArrayXf &snp2pqPowS, const VectorXf &snp2pq, const Vector2f &sigmaSq, const Vector3f &pi, const float vare, const float varg, const float ps, const float overdispersion, VectorXf &deltaS) {
    long numChr = chromInfoVec.size();
    
    wtdSum2pq.setZero();
    wtdSumSq.setZero();
    numNonZeros.setZero();
    numSnpMixComp.setZero();
    
    valuesMixCompS.setZero(size);
    deltaS.setZero(size);
    
//    for (unsigned chr=0; chr<numChr; ++chr) {
//        ChromInfo *chromInfo = chromInfoVec[chr];
//        unsigned chrStart = chromInfo->startSnpIdx;
//        unsigned chrEnd   = chromInfo->endSnpIdx;
//        if (iter==0) {
//            cout << "chr " << chr+1 << " start " << chrStart << " end " << chrEnd << endl;
//        }
//    }
//    if (iter==0) cout << endl;
    
    for (unsigned chr=0; chr<numChr; ++chr) {
        //cout << " thread " << omp_get_thread_num() << " chr " << chr << endl;
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        
        float oldSample;
        float sampleDiff;
        float rhs;
        float varei;
        float logSigmaSqC = log(sigmaSq[0]);
        
        Array3f logPi = pi.array().log();  // zero, C, S
        Array2f invSigmaSq = sigmaSq.cwiseInverse();
        Array3f invLhs;
        Array3f uhat;
        Array3f logDelta;
        Array3f probDelta;
        Array3f weight; weight << 0, 1, 1;
        
        unsigned delta;
        
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            
            oldSample = values[i];
            weight[2] = snp2pqPowS[i];
            
            varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
            
            rhs  = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            
            invLhs[0] = 0.0;
            invLhs[1] = 1.0f/(ZPZdiag[i]/varei + invSigmaSq[0]);
            invLhs[2] = 1.0f/(ZPZdiag[i]/varei + invSigmaSq[1]/snp2pqPowS[i]);
            uhat = invLhs*rhs;
            
            logDelta[0] = logPi[0];
            logDelta[1] = 0.5*(logf(invLhs[1]) - logSigmaSqC + uhat[1]*rhs) + logPi[1];
            logDelta[2] = 0.5*(logf(invLhs[2]) - logf(snp2pqPowS[i]*sigmaSq[1]) + uhat[2]*rhs) + logPi[2];
            
            for (unsigned j=0; j<3; ++j) {
                probDelta[j] = 1.0f/(logDelta-logDelta[j]).exp().sum();
            }
            
            delta = bernoulli.sample(probDelta);
            numSnpMixComp[delta]++;
            
            if (delta) {
                values[i] = normal.sample(uhat[delta], invLhs[delta]);
                sampleDiff = oldSample - values[i];
                for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                    rcorr[it.index()] += it.value() * sampleDiff;
                }
                wtdSum2pq[delta-1] += snp2pq[i]*weight[delta];
                wtdSumSq[delta-1]  += values[i]*values[i]/weight[delta];
                if (delta == 2) {
                    valuesMixCompS[i] = values[i];
                    deltaS[i] = 1;
                }
                ++numNonZeros[delta-1];
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

void ApproxBayesSMix::SnpEffects::sampleFromFC(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy, const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo *> &chromInfoVec, const VectorXf &LDsamplVar, const ArrayXf &snp2pqPowS, const VectorXf &snp2pq, const Vector2f &sigmaSq, const Vector3f &pi, const float vare, const float varg, const float ps, const float overdispersion, VectorXf &deltaS) {
    long numChr = chromInfoVec.size();
    
    wtdSum2pq.setZero();
    wtdSumSq.setZero();
    numNonZeros.setZero();
    numSnpMixComp.setZero();
    
    valuesMixCompS.setZero(size);
    deltaS.setZero(size);
    
//    for (unsigned chr=0; chr<numChr; ++chr) {
//        ChromInfo *chromInfo = chromInfoVec[chr];
//        unsigned chrStart = chromInfo->startSnpIdx;
//        unsigned chrEnd   = chromInfo->endSnpIdx;
//        if (iter==0) {
//            cout << "chr " << chr+1 << " start " << chrStart << " end " << chrEnd << endl;
//        }
//    }
//    if (iter==0) cout << endl;
    
    for (unsigned chr=0; chr<numChr; ++chr) {
        //cout << " thread " << omp_get_thread_num() << " chr " << chr << endl;
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        
        float oldSample;
        float sampleDiff;
        float rhs;
        float varei;
        float logSigmaSqC = log(sigmaSq[0]);
        
        Array3f logPi = pi.array().log();  // zero, C, S
        Array2f invSigmaSq = sigmaSq.cwiseInverse();
        Array3f invLhs;
        Array3f uhat;
        Array3f logDelta;
        Array3f probDelta;
        Array3f weight; weight << 0, 1, 1;
        
        unsigned delta;
        
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            
            oldSample = values[i];
            weight[2] = snp2pqPowS[i];
            
            varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
            
            rhs  = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            
            invLhs[0] = 0.0;
            invLhs[1] = 1.0f/(ZPZdiag[i]/varei + invSigmaSq[0]);
            invLhs[2] = 1.0f/(ZPZdiag[i]/varei + invSigmaSq[1]/snp2pqPowS[i]);
            uhat = invLhs*rhs;
            
            logDelta[0] = logPi[0];
            logDelta[1] = 0.5*(logf(invLhs[1]) - logSigmaSqC + uhat[1]*rhs) + logPi[1];
            logDelta[2] = 0.5*(logf(invLhs[2]) - logf(snp2pqPowS[i]*sigmaSq[1]) + uhat[2]*rhs) + logPi[2];
            
            for (unsigned j=0; j<3; ++j) {
                probDelta[j] = 1.0f/(logDelta-logDelta[j]).exp().sum();
            }
            
            delta = bernoulli.sample(probDelta);
            
            numSnpMixComp[delta]++;
            
            if (delta) {
                values[i] = normal.sample(uhat[delta], invLhs[delta]);
                sampleDiff = oldSample - values[i];
                rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*(oldSample - values[i]);
                wtdSum2pq[delta-1] += snp2pq[i]*weight[delta];
                wtdSumSq[delta-1]  += values[i]*values[i]/weight[delta];
                if (delta == 2) {
                    valuesMixCompS[i] = values[i];
                    deltaS[i] = 1;
                }
                ++numNonZeros[delta-1];
            } else {
                if (oldSample) {
                    rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*oldSample;
                }
                values[i] = 0.0;
            }
        }
    }
}


void ApproxBayesSMix::PiMixComp::sampleFromFC(const VectorXf &numSnpMixComp) {
    VectorXf alphaTilde;
    alphaTilde = numSnpMixComp + alpha;
    values = Dirichlet::sample(ndist, alphaTilde);
    for (unsigned i=0; i<ndist; ++i) {
        (*this)[i]->value = values[i];
    }
}

void ApproxBayesSMix::VarEffects::sampleFromFC(const Vector2f &snpEffSumSq, const Vector2f &numSnpEff) {
    for (unsigned i=0; i<2; ++i) {
        (*this)[i]->sampleFromFC(snpEffSumSq[i], numSnpEff[i]);
        values[i] = (*this)[i]->value;
    }
}

void ApproxBayesSMix::VarEffects::computeScale(const Vector2f &varg, const Vector2f &wtdSum2pq) {
    for (unsigned i=0; i<2; ++i) {
        (*this)[i]->computeScale(varg[i], wtdSum2pq[i]);
    }
}

void ApproxBayesSMix::GenotypicVarMixComp::compute(const Vector2f &sigmaSq, const Vector2f &wtdSum2pq) {
    values = sigmaSq.cwiseProduct(wtdSum2pq);
    for (unsigned i=0; i<2; ++i) {
        (*this)[i]->value = values[i];
    }
}

void ApproxBayesSMix::HeritabilityMixComp::compute(const Vector2f &vargMixComp, const float varg, const float vare) {
    for (unsigned i=0; i<2; ++i) {
        (*this)[i]->value = values[i] = vargMixComp[i]/(varg+vare);
    }
}

void ApproxBayesSMix::sampleUnknowns() {
    unsigned cnt=0;
    do {
        if (sparse) {
            snpEffects.sampleFromFC(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.chromInfoVec,
                                    data.LDsamplVar, snp2pqPowS, data.snp2pq, sigmaSq.values,
                                    piMixComp.values, vare.value, varg.value, ps.value, overdispersion, deltaS.values);
        } else {
            snpEffects.sampleFromFC(rcorr, data.ZPZ, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec,
                                    data.LDsamplVar, snp2pqPowS, data.snp2pq, sigmaSq.values,
                                    piMixComp.values, vare.value, varg.value, ps.value, overdispersion, deltaS.values);
        }
        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros.sum() == 0);
    piMixComp.sampleFromFC(snpEffects.numSnpMixComp);
    pi.sampleFromFC(data.numIncdSnps, snpEffects.numNonZeros.sum());
    sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
    nnzSnp.getValue(snpEffects.numNonZeros.sum());
    varg.compute(snpEffects.values, data.ZPy, rcorr, covg.value);
    vargMixComp.compute(sigmaSq.values, snpEffects.wtdSum2pq);
    vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, covg.value);
    hsq.compute(varg.value, vare.value);
    hsqMixComp.compute(vargMixComp.values, varg.value, vare.value);
    
//    if (++iter < 2000) {
//        sigmaSq.computeScale(hsqMixComp.values, snpEffects.wtdSum2pq);
//        scalePrior += (sigmaSq[0]->scale - scalePrior)/iter;
//    } else {
//        sigmaSq[0]->scale = scalePrior;
//    }
    
    S.sampleFromFC(snpEffects.wtdSumSq[1], snpEffects.numNonZeros[1], sigmaSq.values[1], snpEffects.valuesMixCompS, data.snp2pq, snp2pqPowS, logSnp2pq, vargMixComp.values[1], sigmaSq[1]->scale, snpEffects.sum2pqSplusOne);
    
    if (sparse)
        rounding.computeRcorr(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    else
        rounding.computeRcorr(data.ZPy, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    
    if (modelPS) ps.compute(rcorr, data.ZPZdiag, data.LDsamplVar, varg.value, vare.value, data.chisq);
}


void BayesSMix::SnpEffects::sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const ArrayXf &snp2pqPowS, const VectorXf &snp2pq, const Vector2f &sigmaSq, const Vector3f &pi, const float vare, VectorXf &deltaS, VectorXf &ghat, vector<VectorXf> &ghatMixComp){
    
    wtdSum2pq.setZero();
    wtdSumSq.setZero();
    numNonZeros.setZero();
    numSnpMixComp.setZero();
    
    valuesMixCompS.setZero(size);
    deltaS.setZero(size);
    
    ghat.setZero(ycorr.size());
    ghatMixComp[0].setZero(ycorr.size());
    ghatMixComp[1].setZero(ycorr.size());
    
    float oldSample;
    float rhs;
    float invVare = 1.0f/vare;
    float logSigmaSqC = log(sigmaSq[0]);
    
    Array3f logPi = pi.array().log();  // zero, C, S
    Array2f invSigmaSq = sigmaSq.cwiseInverse();
    Array3f invLhs;
    Array3f uhat;
    Array3f logDelta;
    Array3f probDelta;
    Array3f weight; weight << 0, 1, 1;
    
    unsigned delta;
    
    for (unsigned i=0; i<size; ++i) {
        if (!ZPZdiag[i]) continue;
        
        oldSample = values[i];
        weight[2] = snp2pqPowS[i];
        
        rhs = Z.col(i).dot(ycorr);
        rhs += ZPZdiag[i]*oldSample;
        rhs *= invVare;
        
        invLhs[0] = 0.0;
        invLhs[1] = 1.0f/(ZPZdiag[i]*invVare + invSigmaSq[0]);
        invLhs[2] = 1.0f/(ZPZdiag[i]*invVare + invSigmaSq[1]/snp2pqPowS[i]);
        uhat = invLhs*rhs;
        
        logDelta[0] = logPi[0];
        logDelta[1] = 0.5*(logf(invLhs[1]) - logSigmaSqC + uhat[1]*rhs) + logPi[1];
        logDelta[2] = 0.5*(logf(invLhs[2]) - logf(snp2pqPowS[i]*sigmaSq[1]) + uhat[2]*rhs) + logPi[2];
        
        for (unsigned j=0; j<3; ++j) {
            probDelta[j] = 1.0f/(logDelta-logDelta[j]).exp().sum();
        }
        
//        if (iter==837) cout << i << " " << probDelta.transpose() << endl;
        
        delta = bernoulli.sample(probDelta);
        
        numSnpMixComp[delta]++;
        
        if (delta) {
            values[i] = normal.sample(uhat[delta], invLhs[delta]);
            ycorr += Z.col(i) * (oldSample - values[i]);
            ghat  += Z.col(i) * values[i];
            ghatMixComp[delta-1] += Z.col(i) * values[i];
            wtdSum2pq[delta-1] += snp2pq[i]*weight[delta];
            wtdSumSq[delta-1]  += values[i]*values[i]/weight[delta];
            if (delta == 2) {
                valuesMixCompS[i] = values[i];
                deltaS[i] = 1;
            }
            ++numNonZeros[delta-1];
        } else {
            if (oldSample) ycorr += Z.col(i) * oldSample;
            values[i] = 0.0;
        }
    }
    
}

void BayesSMix::GenotypicVarMixComp::compute(const vector<VectorXf> &ghatMixComp){
    for (unsigned i=0; i<2; ++i) {
        (*this)[i]->value = values[i] = Gadget::calcVariance(ghatMixComp[i]);
    }
}

void BayesSMix::sampleUnknowns(){    
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }

    unsigned cnt=0;
    do {
        snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, snp2pqPowS, data.snp2pq, sigmaSq.values, piMixComp.values, vare.value, deltaS.values, ghat, ghatMixComp);
        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros.sum() == 0);
    
    sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
    
    if (estimatePi) {
        pi.sampleFromFC(snpEffects.size, snpEffects.numNonZeros.sum());
        piMixComp.sampleFromFC(snpEffects.numSnpMixComp);
    }
    
    nnzSnp.getValue(snpEffects.numNonZeros.sum());

    varg.compute(ghat);
    vare.sampleFromFC(ycorr);
    hsq.compute(varg.value, vare.value);
    
//    vargMixComp.compute(sigmaSq.values, snpEffects.wtdSum2pq);
    vargMixComp.compute(ghatMixComp);
    
    hsqMixComp.compute(vargMixComp.values, varg.value, vare.value);
    
    S.sampleFromFC(snpEffects.wtdSumSq[1], snpEffects.numNonZeros[1], sigmaSq.values[1], snpEffects.valuesMixCompS, data.snp2pq, snp2pqPowS, logSnp2pq, vargMixComp.values[1], sigmaSq[1]->scale, snpEffects.sum2pqSplusOne);
    
    rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);
    
}


void ApproxBayesRC::SnpEffects::sampleFromFC(VectorXf &rcorr, const vector<SparseVector<float>> &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const VectorXf &LDsamplVar,
                                            const float sigmaSq, const MatrixXf &snpPi, const VectorXf &gamma, const float vare,
                                            const float varg, const float ps, const float overdispersion,
                                            const bool originalModel, DeltaPi &deltaPi){


    long numChr = chromInfoVec.size();

    float ssq[numChr], s2pq[numChr], nnz[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(s2pq,0,sizeof(float)*numChr);
    memset(nnz,0, sizeof(float)*numChr);

    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
    z.setZero(size, ndist-1);   // indicator variables for conditional membership
    
    // R specific parameters
    ArrayXf wtdSigmaSq(ndist);
    ArrayXf invWtdSigmaSq(ndist);
    ArrayXf logWtdSigmaSq(ndist);
    MatrixXf logPi = snpPi.array().log().matrix();
        
    if (originalModel) {
        wtdSigmaSq = gamma * 0.01 * varg;
    } else {
        wtdSigmaSq = gamma * sigmaSq;
    }
    
    invWtdSigmaSq = wtdSigmaSq.inverse();
    logWtdSigmaSq = wtdSigmaSq.log();
    
    numSnpMix.setZero(ndist);
    snpset.resize(ndist);
    
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
        deltaPi[k]->values.setZero(size);
    }


    #pragma omp parallel for  // openmp is not working for SBayesR
    for (unsigned chr=0; chr<numChr; ++chr)
    {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned windEnd, j;
        
        float oldSample;
        float sampleDiff;
        float rhs;
        float varei;
        
        ArrayXf invLhs(ndist);
        ArrayXf uhat(ndist);
        ArrayXf logDelta(ndist);
        ArrayXf probDelta(ndist);
        
        unsigned delta;
                
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            oldSample = valuesPtr[i];
            
            varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
                        
            rhs  = rcorr[i] + ZPZdiag[i] * oldSample;
            rhs /= varei;
            
            invLhs = (ZPZdiag[i]/varei + invWtdSigmaSq).inverse();
            uhat = invLhs*rhs;
            
            logDelta = 0.5*(invLhs.log() - logWtdSigmaSq + uhat*rhs) + logPi.row(i).transpose().array();
            logDelta[0] = logPi(i,0);
            
            for (unsigned k=0; k<ndist; ++k) {
                probDelta[k] = 1.0f/(logDelta-logDelta[k]).exp().sum();
                deltaPi[k]->values[i] = probDelta[k];
            }
            delta = bernoulli.sample(probDelta);
            
            snpset[delta].push_back(i);
            numSnpMix[delta]++;
            
            if (delta) {
                valuesPtr[i] = uhat[delta] + nrnd[i]*sqrtf(invLhs[delta]);
                sampleDiff = oldSample - valuesPtr[i];
                for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                    rcorr[it.index()] += it.value() * sampleDiff;
                }
                ssq[chr] += (valuesPtr[i] * valuesPtr[i]) / gamma[delta];
                ++nnz[chr];
                z(i,0) = 1;
                if (delta > 1) z(i,1) = 1;
                if (delta > 2) z(i,2) = 1;
            }
            else {
                if (oldSample) {
                    for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                        rcorr[it.index()] += it.value() * oldSample;
                    }
                }
                valuesPtr[i] = 0.0;
            }
        }
    }
    
    sumSq = 0.0;
    numNonZeros = 0;
    nnzPerChr.setZero(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        sumSq += ssq[i];
        numNonZeros += nnz[i];
        nnzPerChr[i] = nnz[i];
    }
    values = VectorXf::Map(valuesPtr, size);
}

void ApproxBayesRC::SnpEffects::sampleFromFC(VectorXf &rcorr, const vector<SparseVector<float>> &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                  const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                  const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const VectorXf &LDsamplVar,
                  const float sigmaSq, const MatrixXf &snpPi, const VectorXf &gamma, const float vare,
                  const VectorXf &snpVarg, const float vary, const float ps, const float overdispersion,
                   const bool originalModel, DeltaPi &deltaPi){
    // -----------------------------------------
    // This method uses per-SNP genetic variance
    // -----------------------------------------
    long numChr = chromInfoVec.size();

    float ssq[numChr], s2pq[numChr], nnz[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(s2pq,0,sizeof(float)*numChr);
    memset(nnz,0, sizeof(float)*numChr);

    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
    z.setZero(size, ndist-1);   // indicator variables for conditional membership
    
    // R specific parameters
    MatrixXf wtdSigmaSq(size, ndist);
    MatrixXf invWtdSigmaSq(size, ndist);
    MatrixXf logWtdSigmaSq(size, ndist);
    MatrixXf logPi = snpPi.array().log().matrix();
        
    if (originalModel) {
        for (unsigned i=0; i<size; ++i) {
            wtdSigmaSq.row(i).array() = gamma.array().transpose() * 0.01 * snpVarg[i];
            invWtdSigmaSq.row(i) = wtdSigmaSq.row(i).array().inverse();
            logWtdSigmaSq.row(i) = wtdSigmaSq.row(i).array().log();
        }
    } else {
        throw("Per-SNP-GV model is only available with --original-model!");
    }
    
//    cout << "wtdSigmaSq " << wtdSigmaSq.row(0) << endl;
//    cout << "invWtdSigmaSq " << invWtdSigmaSq.row(0) << endl;
//    cout << "logWtdSigmaSq " << logWtdSigmaSq.row(0) << endl;

    
    numSnpMix.setZero(ndist);
    snpset.resize(ndist);
    
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
        deltaPi[k]->values.setZero(size);
    }

    // --------------------------------------------------------------------------------
    // Cycle over all variants in the window and sample the genetics effects
    // --------------------------------------------------------------------------------

    #pragma omp parallel for  // openmp is not working for SBayesR
    for (unsigned chr=0; chr<numChr; ++chr)
    {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned windEnd, j;
        
        float oldSample;
        float sampleDiff;
        float rhs;
        float varei;
        
        ArrayXf invLhs(ndist);
        ArrayXf uhat(ndist);
        ArrayXf logDelta(ndist);
        ArrayXf probDelta(ndist);
        
        unsigned delta;
        
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            oldSample = valuesPtr[i];
            
            varei = vary;
//            varei = LDsamplVar[i]*varg[] + vare + ps + overdispersion;
                        
            rhs  = rcorr[i] + ZPZdiag[i] * oldSample;
            rhs /= varei;
            
            invLhs = (ZPZdiag[i]/varei + invWtdSigmaSq.row(i).transpose().array()).inverse();
            uhat = invLhs*rhs;
            
            logDelta = 0.5*(invLhs.log() - logWtdSigmaSq.row(i).transpose().array() + uhat*rhs) + logPi.row(i).transpose().array();
            logDelta[0] = logPi(i,0);
//            cout << "i " << i << endl;
//            cout << invWtdSigmaSq.row(i).array() << endl;
//            cout << logPi.row(i).array() << endl;
//            cout << logDelta.transpose() << endl;
            
            for (unsigned k=0; k<ndist; ++k) {
                probDelta[k] = 1.0f/(logDelta-logDelta[k]).exp().sum();
                deltaPi[k]->values[i] = probDelta[k];
            }
                        
            delta = bernoulli.sample(probDelta);
            
            snpset[delta].push_back(i);
            numSnpMix[delta]++;
            
            if (delta) {
                valuesPtr[i] = uhat[delta] + nrnd[i]*sqrtf(invLhs[delta]);
                sampleDiff = oldSample - valuesPtr[i];
                for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                    rcorr[it.index()] += it.value() * sampleDiff;
                }
                ssq[chr] += (valuesPtr[i] * valuesPtr[i]) / gamma[delta];
                ++nnz[chr];
                z(i,0) = 1;
                if (delta > 1) z(i,1) = 1;
                if (delta > 2) z(i,2) = 1;
            }
            else {
                if (oldSample) {
                    for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                        rcorr[it.index()] += it.value() * oldSample;
                    }
                }
                valuesPtr[i] = 0.0;
            }
        }
    }
    
    // ---------------------------------------------------------------------
    // Tally up the effect sum of squares and the number of non-zero effects
    // ---------------------------------------------------------------------
    sumSq = 0.0;
    numNonZeros = 0;
    nnzPerChr.setZero(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        sumSq += ssq[i];
        numNonZeros += nnz[i];
        nnzPerChr[i] = nnz[i];
    }
    values = VectorXf::Map(valuesPtr, size);
}

void ApproxBayesRC::SnpEffects::sampleFromFC(vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
                                             const vector<LDBlockInfo*> keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks,
                                             const MatrixXf &snpPi, const VectorXf &gamma, const float varg,
                                             DeltaPi &deltaPi){
    // -----------------------------------------
    // This method uses low-rank model with eigen-decomposition of LD matrices
    // -----------------------------------------
    long nBlocks = keptLdBlockInfoVec.size();
    
    whatBlocks.resize(nBlocks);
    ssqBlocks.resize(nBlocks);
    for (unsigned i=0; i<nBlocks; ++i) {
        whatBlocks[i].resize(wcorrBlocks[i].size());
    }

    float ssq[nBlocks], s2pq[nBlocks], nnz[nBlocks];
    memset(ssq,0, sizeof(float)*nBlocks);
    memset(s2pq,0,sizeof(float)*nBlocks);
    memset(nnz,0, sizeof(float)*nBlocks);

    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
    z.setZero(size, ndist-1);   // indicator variables for conditional membership
    
    ArrayXf wtdSigmaSq = gamma.array().transpose() * 0.01 * varg;
    ArrayXf invWtdSigmaSq = wtdSigmaSq.inverse();
    ArrayXf logWtdSigmaSq = wtdSigmaSq.log();
    
    MatrixXf logPi = snpPi.array().log().matrix();
        
    
//    cout << "wtdSigmaSq " << wtdSigmaSq.row(0) << endl;
//    cout << "invWtdSigmaSq " << invWtdSigmaSq.row(0) << endl;
//    cout << "logWtdSigmaSq " << logWtdSigmaSq.row(0) << endl;

    
    numSnpMix.setZero(ndist);
    snpset.resize(ndist);
    
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
        deltaPi[k]->values.setZero(size);
    }

    // --------------------------------------------------------------------------------
    // Cycle over all variants in the window and sample the genetics effects
    // --------------------------------------------------------------------------------

    //cout << "Run 1.1" << std::endl;
    #pragma omp parallel for schedule(dynamic)
    for(unsigned blk = 0; blk < nBlocks; blk++){
        Ref<const MatrixXf> Q = Qblocks[blk];
        Ref<VectorXf> wcorr = wcorrBlocks[blk];
        Ref<VectorXf> what = whatBlocks[blk];

        what.setZero();
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[blk];
        
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd   = blockInfo->endSnpIdx;
        
        float vareDn = nGWASblocks[blk] / vareBlocks[blk];

        ArrayXf invLhs = 1.0/(vareDn + invWtdSigmaSq);
        ArrayXf logInvLhsMsigma = invLhs.log() - logWtdSigmaSq;

        for(unsigned i = blockStart; i <= blockEnd; i++){
            float oldSample = valuesPtr[i];
            Ref<const VectorXf> Qi = Q.col(i - blockStart);
            float rhs = (Qi.dot(wcorr) + oldSample)*vareDn;
            ArrayXf uhat = invLhs * rhs;
            ArrayXf logDelta = 0.5*(logInvLhsMsigma + uhat*rhs) + logPi.row(i).transpose().array();
            logDelta[0] = logPi(i,0);
            
//            cout << i << " rhs " << rhs << " vareDn " << vareDn << " invWtdSigmaSq " << invWtdSigmaSq.transpose() << " uhat " << uhat.transpose() << endl;
            
            ArrayXf probDelta(ndist);
            for (unsigned k=0; k<ndist; ++k) {
                probDelta[k] = 1.0f/(logDelta-logDelta[k]).exp().sum();
                if(isnan(probDelta[k])) probDelta[k] = 0;
                deltaPi[k]->values[i] = probDelta[k];
            }

            unsigned delta;
            #pragma omp critical
            {
                delta = bernoulli.sample(probDelta);
                snpset[delta].push_back(i);
                numSnpMix[delta]++;
            }
            
            if (delta) {
                valuesPtr[i] = uhat[delta] + nrnd[i]*sqrtf(invLhs[delta]);
                wcorr += Qi*(oldSample - valuesPtr[i]);                
                what  += Qi* valuesPtr[i];
                ssq[blk] += (valuesPtr[i] * valuesPtr[i]) / gamma[delta];
                ++nnz[blk];
                z(i,0) = 1;
                if (delta > 1) z(i,1) = 1;
                if (delta > 2) z(i,2) = 1;
            }
            else {
                if (oldSample) wcorr += Qi * oldSample;
                valuesPtr[i] = 0.0;
            }
        }

    }
    
    // ---------------------------------------------------------------------
    // Tally up the effect sum of squares and the number of non-zero effects
    // ---------------------------------------------------------------------
    sumSq = 0.0;
    numNonZeros = 0;
    nnzPerBlk.setZero(nBlocks);
    for (unsigned blk=0; blk<nBlocks; ++blk) {
        sumSq += ssq[blk];
        numNonZeros += nnz[blk];
        nnzPerBlk[blk] = nnz[blk];
        ssqBlocks[blk] = ssq[blk];
    }
    values = VectorXf::Map(valuesPtr, size);
 
}

void ApproxBayesRC::AnnoEffects::sampleFromFC_Gibbs(MatrixXf &z, const MatrixXf &annoMat, const VectorXf &sigmaSq, MatrixXf &snpP) {
//    cout << "sampling anno effects..." << endl;
    
//    static unsigned iter=0;
    
    VectorXf numOnes(numComp);
    #pragma omp parallel for
    for (unsigned i=0; i<numComp; ++i) {
        numOnes[i] = z.col(i).sum();
    }
    
    //cout << numOnes.transpose() << endl;
    
    unsigned numSnps = z.rows();
    for (unsigned i=0; i<numComp; ++i) {
        VectorXf &alphai = (*this)[i]->values;
        VectorXf y, zi;
        unsigned numDP;  // number of data points for each component
        if (i==0) numDP = numSnps;
        else numDP = numOnes[i-1];

        if(numDP == 0){
            alphai.setZero();
            alphai[0] = -10.0;
            ssq[i] = 0;
        }else{
            y.setZero(numDP);
            zi.setZero(numDP);
            const MatrixXf *annotMatP;
            MatrixXf annoMatPO;
            // get annotation coefficient matrix for component i
            if (i==0) {
                annotMatP = &annoMat;
                zi = z.col(i);
            } else {
                annoMatPO.setZero(numDP, numAnno);
                for (unsigned j=0, idx=0; j<numSnps; ++j) {
                    if (z(j,i-1)) {
                        annoMatPO.row(idx) = annoMat.row(j);
                        zi[idx] = z(j,i);
                        ++idx;
                    }
                }
                annotMatP = &annoMatPO;
            }
            const MatrixXf &annoMati = (*annotMatP);


            VectorXf annoDiagi(numAnno);
            //        for (unsigned k=1; k<numAnno; ++k) {   // skip the first annotation because the first annotation is the intercept
            //            annoMean[i][k] = annoMati.col(k).mean();
            //            annoMati.col(k).array() -= annoMean[i][k];
            //        }
            if (i==0) {
                annoDiagi = annoDiag;
            } else {
                annoDiagi[0] = numOnes[i-1];
                #pragma omp parallel for
                for (unsigned k=1; k<numAnno; ++k) {
                    annoDiagi[k] = annoMati.col(k).squaredNorm();
                }
            }

            // compute the mean of truncated normal distribution
            VectorXf mean = annoMati * alphai;

            // sample latent variables
            for (unsigned j=0; j<numDP; ++j) {
                //            cout << j << " mean[j] " << mean[j] << " anno " << annoMati.row(j) << endl;
                if (zi[j]) y[j] = TruncatedNormal::sample_lower_truncated(mean[j], 1.0, 0.0);
                else y[j] = TruncatedNormal::sample_upper_truncated(mean[j], 1.0, 0.0);
            }

            // adjust the latent variable by all annotation effects;
            y -= mean;

            // intercept is fitted with a flat prior
            float oldSample = alphai[0];
            float rhs = y.sum() + annoDiagi[0]*oldSample;
            float invLhs = 1.0/annoDiagi[0];
            float ahat = invLhs*rhs;
            alphai[0] = Normal::sample(ahat, invLhs);
            y.array() += oldSample - alphai[0];
            //        cout << i << " alphai[0] " << alphai[0] << endl;

            // annotations are fitted with a normal prior
            ssq[i] = 0;
            for (unsigned k=1; k<numAnno; ++k) {
                oldSample = alphai[k];
                rhs = annoMati.col(k).dot(y) + annoDiagi[k]*oldSample;
                invLhs = 1.0/(annoDiagi[k] + 1.0/sigmaSq[i]);
                ahat = invLhs*rhs;
                alphai[k] = Normal::sample(ahat, invLhs);
                y += annoMati.col(k) * (oldSample - alphai[k]);
                ssq[i] += alphai[k] * alphai[k];
                //            cout << i << " " << k << " " << alphai[k] << " " << ahat << " " << invLhs << " " << annoDiagi[k] << " " << sigmaSq[i] << endl;
            }
        }
        //cout << i << " " << alphai.transpose() << endl;
        
        #pragma omp parallel for
        for (unsigned j=0; j<numSnps; ++j) {
            snpP(j,i) = Normal::cdf_01(annoMat.row(j).dot(alphai));
        }
    }
//    ++iter;
    
//    cout << "sampling anno effects finished." << endl;

}

void ApproxBayesRC::AnnoEffects::sampleFromFC_MH(MatrixXf &z, const MatrixXf &annoMat, const VectorXf &sigmaSq, MatrixXf &snpP) {
    // random-walk Mentropolis-Hastings sampling
    
    //    cout << "sampling anno effects..." << endl;
    
//    static unsigned iter=0;
    
    VectorXf numOnes(numComp);
    for (unsigned i=0; i<numComp; ++i) {
        numOnes[i] = z.col(i).sum();
    }
    
    //cout << numOnes.transpose() << endl;
        
    unsigned numSnps = z.rows();
    for (unsigned i=0; i<numComp; ++i) {
        VectorXf curr_alpha = (*this)[i]->values;
        VectorXf cand_alpha(numAnno);
        for (unsigned k=0; k<numAnno; ++k) {
            cand_alpha[k] = Normal::sample(curr_alpha[k], varProp[i]);
        }

        float logPriorCurr = -0.5f * ((curr_alpha.squaredNorm() - curr_alpha[0]*curr_alpha[0])/sigmaSq[i]);  // first annotation is intercept which has a flat prior
        float logPriorCand = -0.5f * ((cand_alpha.squaredNorm() - cand_alpha[0]*cand_alpha[0])/sigmaSq[i]);
        
        float logLikCurr = 0.0;
        float logLikCand = 0.0;
        
        for (unsigned j=0; j<numSnps; ++j) {
            if (i==0) {
                float curr_p = snpP(j,i);
                float cand_p = 1.0/(1.0 + expf(- annoMat.row(j).dot(cand_alpha)));
                if (z(j,i)) {
                    logLikCurr += logf(curr_p);
                    logLikCand += logf(cand_p);
                } else {
                    logLikCurr += logf(1.0 - curr_p);
                    logLikCand += logf(1.0 - cand_p);
                }
            }
            else {
                if (z(j,i-1)) {
                    float curr_p = snpP(j,i);
                    float cand_p = 1.0/(1.0 + expf(- annoMat.row(j).dot(cand_alpha)));
                    if (z(j,i)) {
                        logLikCurr += logf(curr_p);
                        logLikCand += logf(cand_p);
                    } else {
                        logLikCurr += logf(1.0 - curr_p);
                        logLikCand += logf(1.0 - cand_p);
                    }
                }
            }
        }
        
        float logPostCurr = logLikCurr + logPriorCurr;
        float logPostCand = logLikCand + logPriorCand;
        
        if (Stat::ranf() < exp(logPostCand-logPostCurr)) {  // accept
            (*this)[i]->values = cand_alpha;
            snpP.col(i) = 1.0/(1.0 + (- (annoMat * cand_alpha).array()).exp());
            ar[i]->count(1, 0.1, 0.5);
        } else {
            ar[i]->count(0, 0.1, 0.5);
        }

        if (!(ar[i]->cnt % 10)) {
            if      (ar[i]->value < 0.2) varProp[i] *= 0.8;
            else if (ar[i]->value > 0.5) varProp[i] *= 1.2;
        }
        
        ssq[i] = (*this)[i]->values.squaredNorm() - (*this)[i]->values[0]*(*this)[i]->values[0];
       //cout << i << " " << alphai.transpose() << endl;
    }
//    ++iter;
}

void ApproxBayesRC::VarAnnoEffects::sampleFromFC(const VectorXf &ssq){
//    cout << "sampling anno effects variance..." << endl;
    for (unsigned i=0; i<size; ++i) {
        float dfTilde = df + (numAnno-1); // exclude the intercept
        float scaleTilde = ssq[i] + df*scale;
        values[i] = InvChiSq::sample(dfTilde, scaleTilde);
    }
}

//void ApproxBayesRC::AnnoEffects::sampleFromFC(MatrixXf &snpP, const MatrixXf &annoMat) {
////    cout << "\nanno sample " << endl;
//    wcorr = (snpP.array()/(1.0-snpP.array())).log().matrix() - wcorr;
//    for (unsigned i=0; i<numComp; ++i) {
////        cout << endl;
//        varwcorr[i] = wcorr.col(i).squaredNorm()/float(wcorr.rows());
////        cout << "wcorr.col(i) " << wcorr.col(i).head(5).transpose() << endl;
//        for (unsigned j=0; j<numAnno; ++j) {
//            float oldSample = (*this)[i]->values[j];
//            float rhs = annoMat.col(j).dot(wcorr.col(i));
//            rhs += annoDiag[j]*oldSample;
//            float invLhs = 1.0f/annoDiag[j];
//            float ahat = invLhs*rhs;
//            float sample = Normal::sample(ahat, invLhs*varwcorr[i]);
//            wcorr.col(i) += annoMat.col(i) * (oldSample - sample);
//            (*this)[i]->values[j] = sample;
////            cout << j << " oldSample " << oldSample << " annoDiag " << annoDiag[j] << " newSample " << sample << endl;
//        }
////        cout << (*this)[i]->values.transpose() << endl;
//        snpP.col(i) = 1.0/(1.0 + (- annoMat * (*this)[i]->values).array().exp());
//        wcorr.col(i) = (snpP.col(i).array()/(1.0-snpP.col(i).array())).log();
//    }
//}

//void ApproxBayesRC::computePfromPi(const MatrixXf &snpPi, MatrixXf &snpP) {
//    snpP.col(0) = snpPi.col(1) + snpPi.col(2) + snpPi.col(3);
//    snpP.col(1) = (snpPi.col(2).array() + snpPi.col(3).array()) / snpP.col(0).array();
//    snpP.col(2) = snpPi.col(3).array() / (snpPi.col(2).array() + snpPi.col(3).array());
//    unsigned nrow = snpP.rows();
//    unsigned ncol = snpP.cols();
//    float lb = 0.000001;
//    float ub = 0.999999;
//    for (unsigned i=0; i<nrow; ++i) {
//        for (unsigned j=0; j<ncol; ++j) {
//            if (snpP(i,j) < lb) snpP(i,j) = lb;
//            if (snpP(i,j) > ub) snpP(i,j) = ub;
//        }
//    }
//}

void ApproxBayesRC::computePiFromP(const MatrixXf &snpP, MatrixXf &snpPi) {
//    cout << "computing Pi from p ..." << endl;
//    cout << "snpP" << endl;
//    cout << snpP << endl;
    unsigned numDist = snpPi.cols();
    unsigned numSnps = snpPi.rows();
    
    for (unsigned i=0; i<numDist; ++i) {
        if (i < numDist-1) snpPi.col(i) = (1.0 - snpP.col(i).array());
        else snpPi.col(i).setOnes();
        if (i) {
            for (unsigned j=0; j<i; ++j) {
                snpPi.col(i).array() *= snpP.col(j).array();
            }
        }
    }
//    snpPi.col(0) = 1.0 - snpP.col(0).array();
//    snpPi.col(1) = (1.0 - snpP.col(1).array()) * snpP.col(0).array();
//    snpPi.col(2) = (1.0 - snpP.col(2).array()) * snpP.col(0).array() * snpP.col(1).array();
//    snpPi.col(3) = snpP.col(0).array() * snpP.col(1).array() * snpP.col(2).array();
    
//    cout << snpPi.row(0) << endl;
//    cout << snpPi.row(1) << endl;
//    cout << snpPi.row(2) << endl;
}

void ApproxBayesRC::initSnpPandPi(const VectorXf &pis, const unsigned numSnps, MatrixXf &snpP, MatrixXf &snpPi) {
    unsigned ndist = pis.size();
    snpP.setZero(numSnps, ndist-1);
    snpPi.setZero(numSnps, ndist);
    VectorXf p(ndist-1);
    
    for (unsigned i=1; i<ndist; ++i) {
        p[i-1] = pis.tail(ndist-i).sum();
        if (i>1) p[i-1] /= pis.tail(ndist-i+1).sum();
    }
//    p[0] = pis[1] + pis[2] + pis[3];
//    p[1] = (pis[2] + pis[3]) / p[0];
//    p[2] = pis[3] / (pis[2] + pis[3]);
    for (unsigned i=0; i<numSnps; ++i) {
        snpPi.row(i) = pis;
        snpP.row(i) = p;
    }
}

void ApproxBayesRC::AnnoEffects::initIntercept_probit(const VectorXf &pis){
    VectorXf p(numComp);
    unsigned ndist = pis.size();
    for (unsigned i=1; i<ndist; ++i) {
        p[i-1] = pis.tail(ndist-i).sum();
        if (i>1) p[i-1] /= pis.tail(ndist-i+1).sum();
    }
//    p[0] = pis[1] + pis[2] + pis[3];
//    p[1] = (pis[2] + pis[3]) / p[0];
//    p[2] = pis[3] / (pis[2] + pis[3]);
    for (unsigned i = 0; i<numComp; ++i) {
        (*this)[i]->values[0] = Normal::quantile_01(p[i]);
    }
}

void ApproxBayesRC::AnnoEffects::initIntercept_logistic(const VectorXf &pis){
    VectorXf p(numComp);
    unsigned ndist = pis.size();
    for (unsigned i=1; i<ndist; ++i) {
        p[i-1] = pis.tail(ndist-i).sum();
        if (i>1) p[i-1] /= pis.tail(ndist-i+1).sum();
    }
//    p[0] = pis[1] + pis[2] + pis[3];
//    p[1] = (pis[2] + pis[3]) / p[0];
//    p[2] = pis[3] / (pis[2] + pis[3]);
    for (unsigned i = 0; i<numComp; ++i) {
        (*this)[i]->values[0] = log(p[i]/(1-p[i]));
    }
}

void ApproxBayesRC::AnnoCondProb::compute_probit(const AnnoEffects &annoEffects, const VectorXf &annoSD){
    for (unsigned i=0; i<annoEffects.numComp; ++i) {
        for (unsigned j=0; j<annoEffects.numAnno; ++j) {
            VectorXf &alpha = annoEffects[i]->values;
            if (j==0) (*this)[i]->values[j] = Normal::cdf_01(alpha[j]);
            else (*this)[i]->values[j] = Normal::cdf_01(alpha[0] + annoSD[j]*alpha[j]);  // NEW
        }
    }
}

void ApproxBayesRC::AnnoCondProb::compute_logistic(const AnnoEffects &annoEffects, const VectorXf &annoSD){
//    cout << "computing conditional prob... " << endl;
    for (unsigned i=0; i<numComp; ++i) {
        for (unsigned j=0; j<annoEffects.numAnno; ++j) {
            VectorXf &alpha = annoEffects[i]->values;
            if (j==0) (*this)[i]->values[j] = 1.0/(1.0 + exp(-alpha[j]));
            else (*this)[i]->values[j] = 1.0/(1.0 + exp(- alpha[0] - annoSD[j]*alpha[j]));
        }
    }
}

void ApproxBayesRC::AnnoJointProb::compute(const AnnoCondProb &annoCondProb){
//    cout << "computing joint prob... " << endl;
    for (unsigned k=0; k<annoCondProb.numAnno; ++k) {
        for (unsigned i=0; i<numDist; ++i) {
           if (i < numDist-1) (*this)[i]->values[k] = 1.0 - annoCondProb[i]->values[k];
            else (*this)[i]->values[k] = 1.0;
            if (i) {
                for (unsigned j=0; j<i; ++j) {
                    (*this)[i]->values[k] *= annoCondProb[j]->values[k];
                }
            }
        }
//        (*this)[0]->values[j] =  1.0 - annoCondProb[0]->values[j];
//        (*this)[1]->values[j] = (1.0 - annoCondProb[1]->values[j]) * annoCondProb[0]->values[j];
//        (*this)[2]->values[j] = (1.0 - annoCondProb[2]->values[j]) * annoCondProb[0]->values[j] * annoCondProb[1]->values[j];
//        (*this)[3]->values[j] = annoCondProb[0]->values[j] * annoCondProb[1]->values[j] * annoCondProb[2]->values[j];
    }
//    cout << "computing joint prob finished." << endl;
}

void ApproxBayesRC::AnnoGenVar::compute(const VectorXf &snpEffects, const vector<vector<unsigned> > &snpset, const VectorXf &ZPy, const VectorXf &rcorr, const MatrixXf &annoMat){
    for (unsigned i=0; i<numComp; ++i) {
        (*this)[i]->values.setZero(numAnno);
        unsigned size = snpset[i+1].size();
        for (unsigned j=0; j<size; ++j) {
            unsigned snpIdx = snpset[i+1][j];
            float varj = snpEffects[snpIdx] * (ZPy[snpIdx] - rcorr[snpIdx]);
            for (unsigned k=0; k<numAnno; ++k) {
                if (annoMat(snpIdx,k) > 0) {  // for centered annotations
                    (*this)[i]->values[k] += varj;
                }
            }
        }
        (*this)[i]->values.array() /= nobs;
    }
}

void ApproxBayesRC::AnnoGenVar::compute(const VectorXf &snpEffects, const vector<vector<unsigned> > &snpset, const MatrixXf &annoMat){
    for (unsigned i=0; i<numComp; ++i) {
        (*this)[i]->values.setZero(numAnno);
        unsigned size = snpset[i+1].size();
        for (unsigned j=0; j<size; ++j) {
            unsigned snpIdx = snpset[i+1][j];
            float varj = snpEffects[snpIdx] * snpEffects[snpIdx];
            for (unsigned k=0; k<numAnno; ++k) {
                if (annoMat(snpIdx,k) > 0) {  // for centered annotations
                    (*this)[i]->values[k] += varj;
                }
            }
        }
    }
}

void ApproxBayesRC::AnnoTotalGenVar::compute(const AnnoGenVar &annoGenVar){
    values.setZero(size);
    for (unsigned i=0; i<annoGenVar.numComp; ++i) {
        values.array() += annoGenVar[i]->values.array();
    }
}

void ApproxBayesRC::AnnoPerSnpHsqEnrichment::compute(const VectorXf &annoTotalGenVar, const float varg){
    values = annoTotalGenVar.array()/varg * invSnpProp.array();
//    cout << "AnnoPerSnpHsqEnrichment " << values.transpose() << endl;
}

void ApproxBayesRC::AnnoPerSnpHsqEnrichment::compute(const VectorXf &snpEffects, const MatrixXf &annoMat, const unsigned nnz){
    unsigned numAnnos = annoMat.cols();
    unsigned numSnps = snpEffects.size();
    VectorXf y(nnz);
    MatrixXf X(nnz, numAnnos);
    for (unsigned i=0, j=0; i<numSnps; ++i) {
        if (snpEffects[i]) {
            y[j] = snpEffects[i] * snpEffects[i];
            X.row(j) = annoMat.row(i);
            ++j;
        }
    }
    MatrixXf XPX = X.transpose() * X;
    for (unsigned k=0; k<numAnnos; ++k) {
        XPX(k,k) += 0.01;
    }
    VectorXf XPy = X.transpose() * y;
    values = XPX.householderQr().solve(XPy);
}
    
void ApproxBayesRC::computeSnpVarg(const MatrixXf &annoMat, const VectorXf &annoPerSnpHsqEnrich, const float varg, const unsigned numSnps){
    VectorXf tau = annoPerSnpHsqEnrich.array() - 1.0;
    for (unsigned i=0; i<numSnps; ++i) {
        snpVarg[i] = (1.0 + annoMat.row(i).dot(tau))*varg;
    }
}

void ApproxBayesRC::AnnoDistribution::compute(const MatrixXf &z, const MatrixXf &annoMat, const ArrayXf &numSnpMix){
    unsigned numSnps = z.rows();
    VectorXi delta = z.rowwise().sum().cast<int>();
    for (unsigned i=0; i<numDist; ++i) {
//        cout << "i " << i << endl;
        (*this)[i]->values.setZero(numAnno);
        unsigned nsnpDisti = numSnpMix[i];
//        cout << "i " << i << " " << nsnpDisti << endl;
        if (nsnpDisti == 0) continue;
        MatrixXf annoMatCompi(nsnpDisti, numAnno);
        unsigned idx = 0;
        for (unsigned j=0; j<numSnps; ++j) {
            if (delta[j] == i) {
                annoMatCompi.row(idx) = annoMat.row(j);
                ++idx;
            }
        }
        for (unsigned k=0; k<numAnno; ++k) {
            VectorXf annoSrt = annoMatCompi.col(k);
            std::sort(annoSrt.data(), annoSrt.data() + annoSrt.size());
//            cout << "k " << k << " " << annoSrt.size() << endl;
            (*this)[i]->values[k] = annoSrt[annoSrt.size()/2];   // median value
        }
    }
}


void ApproxBayesRC::sampleUnknowns(){
    static int iter = 0;    
//    fixedEffects.sampleFromFC(data.XPX, data.XPXdiag, data.ZPX, data.XPy, snpEffects.values, vare.value, rcorr);
    unsigned cnt=0;
    do {
        if (lowRankModel) {
            snpEffects.sampleFromFC(wcorrBlocks, data.Qblocks, whatBlocks,
                                    data.keptLdBlockInfoVec, data.nGWASblock, vareBlk.values,
                                    snpPi, gamma.values, varg.value, deltaPi);
            //cout << snpEffects.values << endl;
        }
        else if (allowPerSnpGV) {
            computeSnpVarg(data.annoMat, annoPerSnpHsqEnrich.values, varg.value, data.numIncdSnps);
//            cout << "iter " << iter << endl;
//            cout << snpVarg.head(5).transpose() << endl;
//            cout << snpVarg.tail(5).transpose() << endl;
            snpEffects.sampleFromFC(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, data.se, data.tss, varei, data.n, data.LDsamplVar, sigmaSq.value, snpPi, gamma.values, vare.value,
                snpVarg, data.varPhenotypic, ps.value, overdispersion, originalModel, deltaPi);
        } else {
            snpEffects.sampleFromFC(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, data.se, data.tss, varei, data.n, data.LDsamplVar, sigmaSq.value, snpPi, gamma.values, vare.value,
                varg.value, ps.value, overdispersion, originalModel, deltaPi);
        }
        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros == 0);
                   
    if (robustMode) {
        if (noscale) {
            sigmaSq.value = varg.value/(data.snp2pq.array().sum()*gamma.values.dot(Pis.values));
        } else {
            sigmaSq.value = varg.value/(data.numIncdSnps*gamma.values.dot(Pis.values));  // LDpred2's parameterisation
        }
    } else {
        if (estimateSigmaSq) sigmaSq.sampleFromFC(snpEffects.sumSq, snpEffects.numNonZeros);
    }
    
    //cout << "snpPi " << snpPi.row(0) << endl;
    //cout << "snpP " << snpP.row(0) << endl;

    if (estimatePi) {
        //Pis.sampleFromFC(snpEffects.numSnpMix);
        //initSnpPandPi(Pis.values, data.numIncdSnps, snpP, snpPi);
        //computePfromPi(snpPi, snpP);
        if (algorithm == gibbs) {
            annoEffects.sampleFromFC_Gibbs(snpEffects.z, data.annoMat, sigmaSqAnno.values, snpP);
            annoCondProb.compute_probit(annoEffects, data.annoSD);
        } else {
            annoEffects.sampleFromFC_MH(snpEffects.z, data.annoMat, sigmaSqAnno.values, snpP);
            annoCondProb.compute_logistic(annoEffects, data.annoSD);
        }
        sigmaSqAnno.sampleFromFC(annoEffects.ssq);
        computePiFromP(snpP, snpPi);
        annoJointProb.compute(annoCondProb);
    }
            
    numSnps.getValues(snpEffects.numSnpMix);
    nnzSnp.getValue(snpEffects.numNonZeros);
    sigmaSqG.compute(sigmaSq.value, snpEffects.sum2pq);

//    cout << "check 2 " << endl;

    if (lowRankModel) {
        vargBlk.compute(whatBlocks);
        vareBlk.sampleFromFC(wcorrBlocks, snpEffects.ssqBlocks, data.nGWASblock, data.numEigenvalBlock);
        varg.value = vargBlk.total;
        vare.value = vareBlk.mean;
    }
    else {
        covg.compute(data.ypy, snpEffects.values, data.ZPy, rcorr);
        varg.compute(snpEffects.values, data.ZPy, rcorr, covg.value);
        vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, covg.value);
    }
        
    //hsq.compute(varg.value, vare.value);
    hsq.value = varg.value / data.varPhenotypic;
    
//    cout << "check 3 " << endl;

    annoGenVar.compute(snpEffects.values, snpEffects.snpset, data.ZPy, rcorr, data.annoMat);
//    cout << "check 4 " << endl;
    annoTotalGenVar.compute(annoGenVar);
//    cout << "check 5 " << endl;
    annoPerSnpHsqEnrich.compute(annoTotalGenVar.values, varg.value);
    //annoPerSnpHsqEnrich.compute(snpEffects.values, data.annoMat, snpEffects.numNonZeros);

//    cout << "check 6 " << endl;

//    annoDist.compute(snpEffects.z, data.annoMat, snpEffects.numSnpMix);
    
//    cout << "check 7 " << endl;

    if (iter >= 2000) sigmaSq.scale = scalePrior;
    scale.getValue(sigmaSq.scale);
    // cout << "iter " << iter << " scalePrior " << scalePrior << "sigmaSq.scale " << sigmaSq.scale << endl;

//    if (sparse)
//        rounding.computeRcorr(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
//    else
//        rounding.computeRcorr(data.ZPy, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    if (modelPS) ps.compute(rcorr, data.ZPZdiag, data.LDsamplVar, varg.value, vare.value, data.chisq);

    nnzSnp.getValue(snpEffects.numNonZeros);
    sigmaSqG.compute(sigmaSq.value, snpEffects.sum2pq);

//    numSnpVg.compute(snpEffects.values, data.ZPZdiag, varg.value, vare.nobs);
    if (originalModel) {
        Vgs.compute(snpEffects.values, data.ZPy, rcorr, snpEffects.snpset, varg.value, vare.nobs);
//        if (sparse)
//            Vgs.compute(snpEffects.values, data.ZPZsp, snpEffects.snpset, varg.value, vare.nobs);
//        else
//            Vgs.compute(snpEffects.values, data.ZPZ, snpEffects.snpset, varg.value, vare.nobs);
    }

    float scaleIteri = 0;
    if (++iter < 2000) {
        if (noscale)
        {
            scaleIteri = 0.5f * varg.value / (data.snp2pq.array().sum()*gamma.values.dot(Pis.values));
        } else
        {
            scaleIteri = 0.5f * varg.value / (data.snp2pq.size()*gamma.values.dot(Pis.values));
        }
        genVarPrior += (varg.value - genVarPrior)/iter;
        scalePrior += (scaleIteri - scalePrior)/iter;
    }
}


void BayesRC::SnpEffects::sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes, const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare, VectorXf &ghat, const MatrixXf &snpPi, const float varg, const bool originalModel, DeltaPi &deltaPi){
    sumSq = 0.0;
    numNonZeros = 0;
    
    ghat.setZero(ycorr.size());

    z.setZero(size, ndist-1);   // indicator variables for conditional membership
    
    // R specific parameters
    ArrayXf wtdSigmaSq(ndist);
    ArrayXf invWtdSigmaSq(ndist);
    ArrayXf logWtdSigmaSq(ndist);
    MatrixXf logPi = snpPi.array().log().matrix();
        
    if (originalModel) {
        wtdSigmaSq = gamma * 0.01 * varg;
    } else {
        wtdSigmaSq = gamma * sigmaSq;
    }
    
    invWtdSigmaSq = wtdSigmaSq.inverse();
    logWtdSigmaSq = wtdSigmaSq.log();
    
    numSnpMix.setZero(ndist);
    snpset.resize(ndist);
    
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
        deltaPi[k]->values.setZero(size);
    }
    
    float invVare = 1.0/vare;

    for (unsigned i=0; i<size; ++i) {
        float oldSample;
        float sampleDiff;
        float rhs;
        
        ArrayXf invLhs(ndist);
        ArrayXf uhat(ndist);
        ArrayXf logDelta(ndist);
        ArrayXf probDelta(ndist);
        
        unsigned delta;
        
        oldSample = values[i];
        rhs  = Z.col(i).dot(ycorr);
        rhs += ZPZdiag[i] * oldSample;
        rhs *= invVare;

        invLhs = (ZPZdiag[i]*invVare + invWtdSigmaSq).inverse();
        uhat = invLhs*rhs;
        
        logDelta = 0.5*(invLhs.log() - logWtdSigmaSq + uhat*rhs) + logPi.row(i).transpose().array();
        logDelta[0] = logPi(i,0);
        
        for (unsigned k=0; k<ndist; ++k) {
            probDelta[k] = 1.0f/(logDelta-logDelta[k]).exp().sum();
            deltaPi[k]->values[i] = probDelta[k];
        }
                    
        delta = bernoulli.sample(probDelta);
        
        snpset[delta].push_back(i);
        numSnpMix[delta]++;
        
        if (delta) {
            values[i] = normal.sample(uhat[delta], invLhs[delta]);
            ycorr += Z.col(i) * (oldSample - values[i]);
            if (weightedRes) ghat += Z.col(i).cwiseProduct(Rsqrt) * values[i];
            else ghat  += Z.col(i) * values[i];
            sumSq += (values[i] * values[i]) / gamma[delta];
            ++numNonZeros;
            z(i,0) = 1;
            if (delta > 1) z(i,1) = 1;
            if (delta > 2) z(i,2) = 1;
        }
        else {
            if (oldSample) ycorr += Z.col(i) * oldSample;
            values[i] = 0.0;
        }
    }
}

void BayesRC::computePiFromP(const MatrixXf &snpP, MatrixXf &snpPi) {
//    cout << "computing Pi from p ..." << endl;
//    cout << "snpP" << endl;
//    cout << snpP << endl;
    snpPi.col(0) = 1.0 - snpP.col(0).array();
    snpPi.col(1) = (1.0 - snpP.col(1).array()) * snpP.col(0).array();
    snpPi.col(2) = (1.0 - snpP.col(2).array()) * snpP.col(0).array() * snpP.col(1).array();
    snpPi.col(3) = snpP.col(0).array() * snpP.col(1).array() * snpP.col(2).array();
}

void BayesRC::initSnpPandPi(const VectorXf &pis, const unsigned numSnps, MatrixXf &snpP, MatrixXf &snpPi) {
    unsigned ndist = pis.size();
    snpP.setZero(numSnps, ndist-1);
    snpPi.setZero(numSnps, ndist);
    VectorXf p(ndist-1);
    p[0] = pis[1] + pis[2] + pis[3];
    p[1] = (pis[2] + pis[3]) / p[0];
    p[2] = pis[3] / (pis[2] + pis[3]);
    for (unsigned i=0; i<numSnps; ++i) {
        snpPi.row(i) = pis;
        snpP.row(i) = p;
    }
}

void BayesRC::sampleUnknowns(){
    static unsigned iter=0;
    
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }

    unsigned cnt=0;
    do {
        snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, data.Rsqrt, data.weightedRes, sigmaSq.value, Pis.values, gamma.values, vare.value, ghat, snpPi, varg.value, originalModel, deltaPi);
        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros == 0);
                
    sigmaSq.sampleFromFC(snpEffects.sumSq, snpEffects.numNonZeros);
    
    //cout << "snpPi " << snpPi.row(0) << endl;
    //cout << "snpP " << snpP.row(0) << endl;

    if (estimatePi) {
        //Pis.sampleFromFC(snpEffects.numSnpMix);
        //initSnpPandPi(Pis.values, data.numIncdSnps, snpP, snpPi);
        //computePfromPi(snpPi, snpP);
        if (algorithm == gibbs) {
            annoEffects.sampleFromFC_Gibbs(snpEffects.z, data.annoMat, sigmaSqAnno.values, snpP);
            annoCondProb.compute_probit(annoEffects, data.annoSD);
        } else {
            annoEffects.sampleFromFC_MH(snpEffects.z, data.annoMat, sigmaSqAnno.values, snpP);
            annoCondProb.compute_logistic(annoEffects, data.annoSD);
        }
        sigmaSqAnno.sampleFromFC(annoEffects.ssq);
        computePiFromP(snpP, snpPi);
        annoJointProb.compute(annoCondProb);
    }
    
    numSnps.getValues(snpEffects.numSnpMix);
    nnzSnp.getValue(snpEffects.numNonZeros);

    vare.sampleFromFC(ycorr);
    varg.compute(ghat);
    hsq.compute(varg.value, vare.value);
    
    annoGenVar.compute(snpEffects.values, snpEffects.snpset, data.annoMat);
    annoTotalGenVar.compute(annoGenVar);
    annoPerSnpHsqEnrich.compute(annoTotalGenVar.values, varg.value);

    if (iter >= 2000) sigmaSq.scale = scalePrior;
    scale.getValue(sigmaSq.scale);
    // cout << "iter " << iter << " scalePrior " << scalePrior << "sigmaSq.scale " << sigmaSq.scale << endl;

    if (originalModel) Vgs.compute(snpEffects.values, data.Z, snpEffects.snpset, varg.value);
    rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);
    nnzSnp.getValue(snpEffects.numNonZeros);

    float scaleIteri = 0;
    if (++iter < 2000) {
        if (noscale)
        {
            scaleIteri = 0.5f * varg.value / (data.snp2pq.array().sum()*gamma.values.dot(Pis.values));
        } else
        {
            scaleIteri = 0.5f * varg.value / (data.snp2pq.size()*gamma.values.dot(Pis.values));
        }
        genVarPrior += (varg.value - genVarPrior)/iter;
        scalePrior += (scaleIteri - scalePrior)/iter;
    }
}

