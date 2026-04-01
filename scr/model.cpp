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
    //if (algorithm == gibbs) {
        gibbsSampler(ycorr, Z, ZPZdiag, Rsqrt, weightedRes, sigmaSq, pi, vare, ghat);
    //} else if (algorithm == hmc) {
    //    hmcSampler(ycorr, Z, ZPZdiag, sigmaSq, pi, vare, ghat);
    //}
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
    
    // shuffling the SNP index for faster convergence
    if (shuffle) Gadget::shuffle_vector(snpIndexVec);

    unsigned i;
    for (unsigned t=0; t<size; ++t) {
        i = snpIndexVec[t];
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

void BayesC::SnpEffects::computePosteriorMean(const unsigned int iter){
    posteriorMean.array() += (values - posteriorMean).array()/(iter+1);
    posteriorMeanPIP.array() += (pip - posteriorMeanPIP).array()/(iter+1);
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
    else {
        unsigned n = 0;
        for (unsigned i = 0; i < (unsigned)snp2pq.size(); ++i)
            if (snp2pq[i] > 0.f) ++n;
        if (!n) n = 1;
        scale = (df-2)/df * varg/(n*pi);
    }
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

void BayesC::HeritabilityFixedEffects::compute(const VectorXf &snpEffects, const Data &data){
    float sumSq = 0;
    for (unsigned i = 0; i < data.numIncdSnps; ++i) {
        if (data.incdSnpInfoVec[i]->fitAsFixedEffect) {
            float a = snpEffects[i];
            sumSq += a * a;
        }
    }
    value = data.varPhenotypic > 0 ? sumSq / data.varPhenotypic : 0;
}

void BayesC::RandomVar::compute(const VectorXf &rhat){
    //value = Gadget::calcVariance(ghat);
    float sum = rhat.sum();
    float ssq = rhat.squaredNorm();
    unsigned size = (unsigned)rhat.size();
    float mean = sum/size;
    value = ssq/size - mean*mean;
}

void BayesC::SnpHsqPEP::compute(const VectorXf &snpEffects, const float varg){
    float rndSnpHsq = varg/float(size);
    values.setZero(size);
    for (unsigned i=0; i<size; ++i) {
        if (snpEffects[i]*snpEffects[i] > rndSnpHsq) values[i] = 1.0;
    }
}

void BayesC::Rounding::computeYcorr(const VectorXf &y, const MatrixXf &X, const MatrixXf &W, const MatrixXf &Z,
                                    const VectorXf &fixedEffects, const VectorXf &randomEffects, const VectorXf &snpEffects,
                                    VectorXf &ycorr){
    VectorXf oldYcorr = ycorr;
    ycorr = y - X*fixedEffects;
    if (randomEffects.size()) ycorr -= W*randomEffects;
    for (unsigned i=0; i<snpEffects.size(); ++i) {
        if (snpEffects[i]) ycorr -= Z.col(i)*snpEffects[i];
    }
    float ss = (ycorr - oldYcorr).squaredNorm();
    value = sqrt(ss);
}

void BayesC::sampleUnknowns(const unsigned iter){
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }

    snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, data.Rsqrt, data.weightedRes, sigmaSq.value, pi.value, vare.value, ghat);
    snpPip.getValues(snpEffects.pip);
    nnzSnp.getValue(snpEffects.numNonZeros);

    sigmaSq.sampleFromFC(snpEffects.sumSq, snpEffects.numNonZeros);

    if (estimatePi) pi.sampleFromFC(snpEffects.size, snpEffects.numNonZeros);
    
    varg.compute(ghat);
    vare.sampleFromFC(ycorr);
    hsq.compute(varg.value, vare.value);
    
    if (!(iter % 100)) rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);
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
    
    pip.setZero(size);
    ghat.setZero(ycorr.size());
    
    float oldSample;
    float rhs, invLhs, uhat;
    float logDelta0, logDelta1, probDelta1;
    float logPi = log(pi);
    float logPiComp = log(1.0-pi);
    float invVare = 1.0f/vare;
    float beta;
    
    if (shuffle) Gadget::shuffle_vector(snpIndexVec);

    unsigned i;
    for (unsigned t=0; t<size; ++t) {
        i = snpIndexVec[t];
        oldSample = values[i];
        rhs = Z.col(i).dot(ycorr);
        rhs += ZPZdiag[i]*oldSample;
        rhs *= invVare;
        invLhs = 1.0f/(ZPZdiag[i]*invVare + 1.0f/sigmaSq[i]);
        uhat = invLhs*rhs;
        logDelta1 = 0.5*(logf(invLhs) - logf(sigmaSq[i]) + uhat*rhs) + logPi;
        logDelta0 = logPiComp;
        probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
        pip[i] = probDelta1;
        
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

void BayesB::sampleUnknowns(const unsigned iter){
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }

    snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, data.Rsqrt, data.weightedRes, sigmaSq.values, pi.value, vare.value, ghat);
    snpPip.getValues(snpEffects.pip);
    nnzSnp.getValue(snpEffects.numNonZeros);

    sigmaSq.sampleFromFC(snpEffects.betaSq);
    
    if (estimatePi) pi.sampleFromFC(snpEffects.size, snpEffects.numNonZeros);
    
    varg.compute(ghat);
    vare.sampleFromFC(ycorr);
    hsq.compute(varg.value, vare.value);
    
    if (!(iter % 100)) rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);
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

void BayesN::sampleUnknowns(const unsigned iter){
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }
    
    snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, data.Rsqrt, data.weightedRes, sigmaSq.value, pi.value, vare.value, ghat);
    snpPip.getValues(snpEffects.pip);
    nnzSnp.getValue(snpEffects.numNonZeros);
    nnzWind.getValue(snpEffects.numNonZeroWind);
    windDelta.getValues(snpEffects.windDelta);

    sigmaSq.sampleFromFC(snpEffects.sumSq, snpEffects.numNonZeros);

    if (estimatePi) pi.sampleFromFC(snpEffects.numWindows, snpEffects.numNonZeroWind);
    
    varg.compute(ghat);
    vare.sampleFromFC(ycorr);
    hsq.compute(varg.value, vare.value);
    
    if (!(iter % 100)) rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);
}

// ----------------------------------------------------------------------------------------
// Bayes R
// ----------------------------------------------------------------------------------------

void BayesR::ProbMixComps::sampleFromFC(const VectorXf &snpStore) {
	VectorXf dirx;
	dirx = snpStore + alphaVec;
    values = Dirichlet::sample(size, dirx);
    for (unsigned i=0; i<size; ++i) {
      (*this)[i]->value = values[i];
    }
}

void BayesR::NumSnpMixComps::getValues(const VectorXf &snpStore) {
    values = snpStore;
    for (unsigned i=0; i<size; ++i) {
        (*this)[i]->value = values[i];
    }
}

void BayesR::VgMixComps::compute(const VectorXf &snpEffects, const MatrixXf &Z, const vector<vector<unsigned> > snpset, const float varg) {
    values.setZero(size);
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

    for (unsigned k=0; k<size; ++k) {
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
    for (unsigned k=0; k<size; ++k) {
        (*this)[k]->value = values[k] = values[k]/sum;
    }

}

void BayesR::SnpEffects::sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes,
                                      const float sigmaSq, const VectorXf &pis, const VectorXf &gamma,
                                      const float vare, VectorXf &ghat, VectorXf &snpStore,
                                      const float varg, const bool hsqPercModel, DeltaPi &deltaPi){
    sumSq = 0.0;
    wtdSumSq = 0.0;
    numNonZeros = 0;
    
    pip.setZero(size);
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
    if (hsqPercModel && varg)
        gp = gamma * 0.01 * varg;
    else
        gp = gamma * sigmaSq;
//    cout << varg << " " << gp.transpose() << endl;
    snpset.resize(ndist);
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
    }
    
    for (unsigned k=0; k<ndist; ++k) {
        deltaPi[k]->values.setZero(size);
    }
    
    if (shuffle) Gadget::shuffle_vector(snpIndexVec);

    unsigned i;
    for (unsigned t = 0; t < size; t++) {
        i = snpIndexVec[t];
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
            deltaPi[k]->values[i] = pll[k];
        }
        pip[i] = 1.0f - pll[0];
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
            sumSq += values[i] * values[i];
            wtdSumSq += (values[i] * values[i]) / gamma[indistflag - 1];
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
    else {
        unsigned n = 0;
        for (unsigned i = 0; i < (unsigned)snp2pq.size(); ++i)
            if (snp2pq[i] > 0.f) ++n;
        if (!n) n = 1;
        scale = (df-2)/df * varg/(n*gamma.dot(pi));
    }
}

void BayesR::sampleUnknowns(const unsigned iter){
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }

    snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, data.Rsqrt, data.weightedRes, sigmaSq.value, Pis.values, gamma.values, vare.value, ghat, snpStore, varg.value, hsqPercModel, deltaPi);
    snpPip.getValues(snpEffects.pip);
    nnzSnp.getValue(snpEffects.numNonZeros);
    numSnps.getValues(snpStore);

    sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
    
    if (estimatePi) Pis.sampleFromFC(snpStore);

    varg.compute(ghat);
    vare.sampleFromFC(ycorr);
    hsq.compute(varg.value, vare.value);
    Vgs.compute(snpEffects.values, data.Z, snpEffects.snpset, varg.value);
    
    if (!(iter % 100)) rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);
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
                              const float vg, float &scale, float &sum2pqSplusOne, bool scaledGeno){
    //if (algorithm == random_walk) {
    //    randomWalkMHsampler(snpEffWtdSumSq, numNonZeros, sigmaSq, snpEffects, snp2pq, snp2pqPowS, logSnp2pq, vg, scale, sum2pqSplusOne);
    //} else if (algorithm == hmc) {
        hmcSampler(numNonZeros, sigmaSq, snpEffects, snp2pq, snp2pqPowS, logSnp2pq, vg, scale, sum2pqSplusOne, scaledGeno);
    //} else if (algorithm == reg) {
    //    regression(snpEffects, logSnp2pq, snp2pqPowS, sigmaSq);
    //}
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
                            const float vg, float &scale, float &sum2pqSplusOne, bool scaledGeno){
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
    float cand_p = curr_p - 0.5*stepSize * gradientU(curr,  snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, sigmaSq, vg, scaledGeno);
    
    for (unsigned i=0; i<numSteps; ++i) {
        // Make a full step for the position
        cand += stepSize * cand_p;
        if (i < numSteps-1) {
            // Make a full step for the momentum, except at end of trajectory
            cand_p -= stepSize * gradientU(cand, snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, sigmaSq, vg, scaledGeno);
        } else {
            // Make a half step for momentum at the end
            cand_p -= 0.5*stepSize * gradientU(cand, snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, sigmaSq, vg, scaledGeno);
        }
        //cout << i << " " << cand << endl;
    }

    // Evaluate potential (negative log posterior) and kinetic energies at start and end of trajectory
    float scaleCurr, scaleCand;
    float curr_U_chisq, cand_U_chisq;
    float curr_H = computeU(curr, snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, sigmaSq, vg, scaleCurr, curr_U_chisq, scaledGeno) + 0.5*curr_p*curr_p;
    float cand_H = computeU(cand, snpEffectDelta1, snp2pqLogSumDelta1, snp2pqDelta1, logSnp2pqDelta1, sigmaSq, vg, scaleCand, cand_U_chisq, scaledGeno) + 0.5*cand_p*cand_p;
    
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
        snp2pqPowS = scaledGeno ? snp2pq.array().pow(cand + 1.0f) : snp2pq.array().pow(cand);
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

float BayesS::Sp::gradientU(const float S, const ArrayXf &snpEffects, const float snp2pqLogSum, const ArrayXf &snp2pq, const ArrayXf &logSnp2pq, const float sigmaSq, const float vg, bool scaledGeno){
    // compute the first derivative of the negative log posterior
    ArrayXf snp2pqPowS = scaledGeno ? snp2pq.pow(S + 1.0f) : snp2pq.pow(S);
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

float BayesS::Sp::computeU(const float S, const ArrayXf &snpEffects, const float snp2pqLogSum, const ArrayXf &snp2pq, const ArrayXf &logSnp2pq, const float sigmaSq, const float vg, float &scale, float &U_chisq, bool scaledGeno){
    // compute negative log posterior and scale
    ArrayXf snp2pqPowS = scaledGeno ? snp2pq.pow(S + 1.0f) : snp2pq.pow(S);
    float constantA = snp2pqLogSum;
    float constantB = (snpEffects.square()/snp2pqPowS).sum();
    float constantC = scaledGeno ? snp2pqPowS.sum() : (snp2pq*snp2pqPowS).sum();
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


void BayesS::SnpEffects::sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes,
                                      const float sigmaSq, const float pi, const float vare,
                                      const ArrayXf &snp2pqPowS, const VectorXf &snp2pq,
                                      const float vg, float &scale, VectorXf &ghat){
    wtdSumSq = 0.0;
    numNonZeros = 0;
    
    pip.setZero(size);
    ghat.setZero(ycorr.size());

    float oldSample;
    float rhs, invLhs, uhat;
    float logDelta0, logDelta1, probDelta1;
    float logPi = log(pi);
    float logPiComp = log(1.0-pi);
    float invVare = 1.0f/vare;
    float invSigmaSq = 1.0f/sigmaSq;
    
    if (shuffle) Gadget::shuffle_vector(snpIndexVec);

    unsigned i;
    for (unsigned t=0; t<size; ++t) {
        i = snpIndexVec[t];
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
        pip[i] = probDelta1;

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

void BayesS::sampleUnknowns(const unsigned iter){
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }

    snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, data.Rsqrt, data.weightedRes, sigmaSq.value, pi.value, vare.value, snp2pqPowS, data.snp2pq, genVarPrior, sigmaSq.scale, ghat);
    snpPip.getValues(snpEffects.pip);
    nnzSnp.getValue(snpEffects.numNonZeros);

    sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
    
    if (estimatePi) pi.sampleFromFC(snpEffects.size, snpEffects.numNonZeros);

    S.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros, sigmaSq.value, snpEffects.values, data.snp2pq, snp2pqPowS, logSnp2pq, genVarPrior, sigmaSq.scale, snpEffects.sum2pqSplusOne, scaledGeno);   

    varg.compute(ghat);
    vare.sampleFromFC(ycorr);
    hsq.compute(varg.value, vare.value);
        
    if (iter < 1000) {
        genVarPrior += (varg.value - genVarPrior)/iter;
        scalePrior  += (sigmaSq.scale - scalePrior)/iter;
        sigmaSq.scale = scalePrior;
    }
    scale.getValue(sigmaSq.scale);
    
    if (!(iter % 100)) rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);
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
    S.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros, sigmaSq.value, snpEffects.values, data.snp2pq, snp2pqPowS, logSnp2pq, varg.value, sigmaSq.scale, snpEffects.sum2pqSplusOne, scaledGeno);
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
    
    pip.setZero(size);
    windPip.setZero(numWindows);

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
                invLhs = 1.0f/(ZPZdiag[j]*invVare + invSigmaSq/snp2pqPowS[j]);
                uhat = invLhs*rhs;
                logDelta1 = 0.5*(logf(invLhs) - logf(snp2pqPowS[j]*sigmaSq) + uhat*rhs) + logLocalPi[i];
                logDelta0 = logLocalPiComp[i];
                
                probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
                pip[j] = probDelta1;

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
                pip[j] = localPi[i];
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

void BayesNS::sampleUnknowns(const unsigned iter){
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }

    snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, data.Rsqrt, data.weightedRes, sigmaSq.value, pi.value, vare.value, snp2pqPowS, data.snp2pq, genVarPrior, sigmaSq.scale, ghat);
    snpPip.getValues(snpEffects.pip);
    nnzSnp.getValue(snpEffects.numNonZeros);
    nnzWind.getValue(snpEffects.numNonZeroWind);
    windDelta.getValues(snpEffects.windDelta);
    
    sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
        
    if (estimatePi) pi.sampleFromFC(snpEffects.numWindows, snpEffects.numNonZeroWind);
    
    S.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros, sigmaSq.value, snpEffects.values, data.snp2pq, snp2pqPowS, logSnp2pq, genVarPrior, sigmaSq.scale, snpEffects.sum2pqSplusOne, scaledGeno);

    varg.compute(ghat);
    vare.sampleFromFC(ycorr);
    hsq.compute(varg.value, vare.value);
    
    if (iter < 1000) {
        genVarPrior += (varg.value - genVarPrior)/iter;
        scalePrior += (sigmaSq.scale - scalePrior)/iter;
        sigmaSq.scale = scalePrior;
    }
    scale.getValue(sigmaSq.scale);

    if (!(iter % 100)) rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);
}


void BayesRS::SnpEffects::sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes, const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare, const ArrayXf &snp2pqPowS, const VectorXf &snp2pq, const float varg, float &scale, VectorXf &ghat, const bool hsqPercModel) {
    
    wtdSumSq = 0.0;
    numNonZeros = 0.0;

    pip.setZero(size);
    ghat.setZero(ycorr.size());
    
    ArrayXf wtdSigmaSq(ndist);
    ArrayXf invWtdSigmaSq(ndist);
    ArrayXf logWtdSigmaSq(ndist);
    ArrayXf logPis = pis.array().log();
    ArrayXf log2pqPowS = snp2pqPowS.log();
    
    if (hsqPercModel && varg) {
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

    if (shuffle) Gadget::shuffle_vector(snpIndexVec);

    unsigned i;
    for (unsigned t=0; t<size; ++t) {
        i = snpIndexVec[t];
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
        pip[i] = 1.0f - probDelta[0];

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
        for (unsigned i=0; i<nnzMix; ++i) sum2pqSplusOne += snp2pqMix[i].pow(value+1.0).sum();
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

void BayesRS::sampleUnknowns(const unsigned iter) {
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }

    snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, data.Rsqrt, data.weightedRes, sigmaSq.value, Pis.values, gamma.values, vare.value, snp2pqPowS, data.snp2pq, varg.value, sigmaSq.scale, ghat, hsqPercModel);
    snpPip.getValues(snpEffects.pip);
    nnzSnp.getValue(snpEffects.numNonZeros);
    numSnps.getValues(snpEffects.numSnpMix);

    sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
    
    if (estimatePi) Pis.sampleFromFC(snpEffects.numSnpMix);
    
    // BayesRS::Sp::sampleFromFC has different signature, doesn't use scaledGeno
    S.sampleFromFC(snpEffects.snpset, snpEffects.values, sigmaSq.value, gamma.values, data.snp2pq, snp2pqPowS, logSnp2pq, genVarPrior, sigmaSq.scale, snpEffects.sum2pqSplusOne);

    varg.compute(ghat);
    vare.sampleFromFC(ycorr);
    hsq.compute(varg.value, vare.value);
    Vgs.compute(snpEffects.values, data.Z, snpEffects.snpset, varg.value);
    
    if (iter < 1000) {
        genVarPrior += (varg.value - genVarPrior)/iter;
        scalePrior  += (sigmaSq.scale - scalePrior)/iter;
        sigmaSq.scale = scalePrior;
    }
    scale.getValue(sigmaSq.scale);
    
    if (!(iter % 100)) rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);
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

void ApproxBayesC::FixedEffects::sampleFromFC_eigen(SnpEffects &snpEffects, vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
    const vector<LDBlockInfo*> &keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks, const unsigned numFixedEffectSnps){
    if (!numFixedEffectSnps) return;
    vector<float> nrnd(snpEffects.size);
    for (unsigned i=0; i<snpEffects.size; ++i) nrnd[i] = Stat::snorm();
    unsigned nBlocks = keptLdBlockInfoVec.size();
    whatBlocks.resize(nBlocks);
    for (unsigned blk=0; blk<nBlocks; ++blk) {
        whatBlocks[blk].resize(wcorrBlocks[blk].size());
    }
    for (unsigned blk=0; blk<nBlocks; ++blk) {
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[blk];
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd = blockInfo->endSnpIdx;
        Ref<const MatrixXf> Q = Qblocks[blk];
        Ref<VectorXf> wcorr = wcorrBlocks[blk];
        Ref<VectorXf> what = whatBlocks[blk];
        float invVareDn = nGWASblocks[blk] / vareBlocks[blk];
        for (unsigned si=blockStart; si<=blockEnd; ++si) {
            SnpInfo *snp = blockInfo->snpInfoVec[si - blockStart];
            if (snp->skip) continue;
            if (!snp->fitAsFixedEffect) continue;
            float *valuesPtr = snpEffects.values.data();
            float *pipPtr = snpEffects.pip.data();
            float oldSample = valuesPtr[si];
            int qcol = eigenQColIndex(blockInfo, si);
            if (qcol < 0) continue;
            Ref<const VectorXf> Qi = Q.col(qcol);
            float rhs = (Qi.dot(wcorr) + oldSample) * invVareDn;
            float invLhsFe = 1.0f / invVareDn;
            float uhat = invLhsFe * rhs;
            valuesPtr[si] = uhat + nrnd[si] * sqrtf(invLhsFe);
            wcorr += Qi * (oldSample - valuesPtr[si]);
            what += Qi * valuesPtr[si];
            pipPtr[si] = 1.0f;
        }
    }
}

void ApproxBayesC::FixedEffects::sampleFromFC_sparse(SnpEffects &snpEffects, VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag,
    const VectorXf &snp2pq, const vector<SnpInfo*> &incdSnpInfoVec, float varg, float vare, const unsigned numFixedEffectSnps){
    if (!numFixedEffectSnps) return;
    float varei = varg + vare;
    vector<float> nrnd(snpEffects.size);
    for (unsigned i=0; i<snpEffects.size; ++i) nrnd[i] = Stat::snorm();
    for (unsigned i=0; i<incdSnpInfoVec.size(); ++i) {
        if (!incdSnpInfoVec[i]->fitAsFixedEffect) continue;
        float *valuesPtr = snpEffects.values.data();
        float *pipPtr = snpEffects.pip.data();
        float oldSample = valuesPtr[i];
        float rhs = rcorr[i] + ZPZdiag[i] * oldSample;
        float v1 = ZPZdiag[i];
        valuesPtr[i] = rhs / v1 + nrnd[i] * sqrtf(varei / v1);
        float sampleDiff = oldSample - valuesPtr[i];
        for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
            rcorr[it.index()] += it.value() * sampleDiff;
        }
        pipPtr[i] = 1.0f;
    }
}

void ApproxBayesC::FixedEffects::sampleFromFC_full(SnpEffects &snpEffects, VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag,
    const VectorXi &windStart, const VectorXi &windSize, const VectorXf &snp2pq, const vector<SnpInfo*> &incdSnpInfoVec, float varg, float vare, const unsigned numFixedEffectSnps){
    if (!numFixedEffectSnps) return;
    float varei = varg + vare;
    vector<float> nrnd(snpEffects.size);
    for (unsigned i=0; i<snpEffects.size; ++i) nrnd[i] = Stat::snorm();
    for (unsigned i=0; i<incdSnpInfoVec.size(); ++i) {
        if (!incdSnpInfoVec[i]->fitAsFixedEffect) continue;
        float *valuesPtr = snpEffects.values.data();
        float *pipPtr = snpEffects.pip.data();
        float oldSample = valuesPtr[i];
        double rhs = rcorr[i] + ZPZdiag[i] * oldSample;
        double v1 = ZPZdiag[i];
        valuesPtr[i] = rhs / v1 + nrnd[i] * sqrtf(varei / v1);
        rcorr.segment(windStart[i], windSize[i]) += ZPZ[i] * (oldSample - valuesPtr[i]);
        pipPtr[i] = 1.0f;
    }
}

void ApproxBayesC::SnpEffects::sampleFromFC_sparse(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                                   const vector<ChromInfo*> &chromInfoVec, const VectorXf &snp2pq,
                                                   const float sigmaSq, const float pi, const float vare, const float varg,
                                                   const vector<SnpInfo*> &incdSnpInfoVec){
    
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

    pip.setZero(size);
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads
    float *pipPtr = pip.data();

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
#pragma omp parallel for schedule(dynamic)
    for (unsigned chr=0; chr<numChr; ++chr) {
        //cout << " thread " << omp_get_thread_num() << " chr " << chr << endl;
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned chrSize = chrEnd - chrStart + 1;

        float oldSample;
        float rhs, invLhs, uhat;
        float logDelta0, logDelta1, probDelta1;
        float logPi = log(pi);
        float logPiComp = log(1.0-pi);
        float logSigmaSq = log(sigmaSq);
        float invSigmaSq = 1.0f/sigmaSq;
        float varei = varg + vare;
        float sampleDiff;
        
        // shuffling the SNP index for faster convergence
        vector<int> snpIndexVec = Gadget::shuffle_index(chrStart, chrEnd);

        //for (unsigned i=chrStart; i<=chrEnd; ++i) {
        unsigned i;
        for (unsigned t = 0; t < chrSize; t++) {
            i = snpIndexVec[t];
            if (badSnps[i]) {
                valuesPtr[i] = 0.0;
                continue;
            }
            if (incdSnpInfoVec[i]->fitAsFixedEffect) {
                continue;
            }
            oldSample = valuesPtr[i];
            rhs = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            invLhs = 1.0f/(ZPZdiag[i]/varei + invSigmaSq);
            uhat = invLhs*rhs;
            logDelta1 = 0.5*(logf(invLhs) - logSigmaSq + uhat*rhs) + logPi;
            logDelta0 = logPiComp;
            probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
            pipPtr[i] = probDelta1;
            
            //            cout << i << " " << rhs << " " << invLhs << " " << uhat << " " << probDelta1 << " " << sigmaSq << endl;
            //            int tmp;
            //            cin >> tmp;
            
            //if (bernoulli.sample(probDelta1)) {
            if (urnd[i] < probDelta1) {
                //valuesPtr[i] = normal.sample(uhat, invLhs);
                valuesPtr[i] = uhat + nrnd[i]*sqrtf(invLhs);
                sampleDiff = oldSample - valuesPtr[i];
                for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                    rcorr[it.index()] += it.value() * sampleDiff;
                }
                ssq[chr]  += valuesPtr[i]*valuesPtr[i];
                s2pq[chr] += snp2pq[i];
                ++nnz[chr];
            } else {
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
    pip = VectorXf::Map(pipPtr, size);
}

void ApproxBayesC::SnpEffects::sampleFromFC_full(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                                 const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec, const VectorXf &snp2pq,
                                                 const float sigmaSq, const float pi, const float vare, const float varg,
                                                 const vector<SnpInfo*> &incdSnpInfoVec){
    
    long numChr = chromInfoVec.size();
    
    float ssq[numChr], nnz[numChr], s2pq[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(nnz,0,sizeof(float)*numChr);
    memset(s2pq,0,sizeof(float)*numChr);

    pip.setZero(size);
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads
    float *pipPtr = pip.data();

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
#pragma omp parallel for schedule(dynamic)
    for (unsigned chr=0; chr<numChr; ++chr) {
        //cout << " thread " << omp_get_thread_num() << " chr " << chr << endl;

        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned chrSize = chrEnd - chrStart + 1;

        float oldSample;
        float rhs, invLhs, uhat;
        float logDelta0, logDelta1, probDelta1;
        float logPi = log(pi);
        float logPiComp = log(1.0-pi);
        float logSigmaSq = log(sigmaSq);
        float invSigmaSq = 1.0f/sigmaSq;
        float varei = varg + vare;
        float sampleDiff;

        // shuffling the SNP index for faster convergence
        vector<int> snpIndexVec = Gadget::shuffle_index(chrStart, chrEnd);

        unsigned i;
        for (unsigned t = 0; t < chrSize; t++) {
            i = snpIndexVec[t];
            if (badSnps[i]) {
                valuesPtr[i] = 0.0;
                continue;
            }
            if (incdSnpInfoVec[i]->fitAsFixedEffect) {
                continue;
            }
            oldSample = valuesPtr[i];
            rhs = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            invLhs = 1.0f/(ZPZdiag[i]/varei + invSigmaSq);
            uhat = invLhs*rhs;
            logDelta1 = 0.5*(logf(invLhs) - logSigmaSq + uhat*rhs) + logPi;
            logDelta0 = logPiComp;
            probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
            pipPtr[i] = probDelta1;
            
            if (urnd[i] < probDelta1) {
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
    pip = VectorXf::Map(pipPtr, size);
}

void ApproxBayesC::SnpEffects::sampleFromFC_eigen(vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
                                            const vector<LDBlockInfo*> keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks,
                                            const float sigmaSq, const float pi, const float varg, const VectorXf &snp2pq){
    
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

    pip.setZero(size);
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads
    float *pipPtr = pip.data();

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
    float logPi = log(pi);
    float logPiComp = log(1.0-pi);
    float logSigmaSq = log(sigmaSq);
    float invSigmaSq = 1.0f/sigmaSq;

#pragma omp parallel for schedule(dynamic)
    for(unsigned blk = 0; blk < nBlocks; blk++){
        Ref<const MatrixXf> Q = Qblocks[blk];
        Ref<VectorXf> wcorr = wcorrBlocks[blk];
        Ref<VectorXf> what = whatBlocks[blk];

        what.setZero();
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[blk];
        
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd   = blockInfo->endSnpIdx;
        unsigned blockSize  = blockEnd - blockStart + 1;

        float invVareDn = nGWASblocks[blk] / vareBlocks[blk];

        float invLhs = 1.0/(invVareDn + invSigmaSq);
        float logInvLhsMsigma = logf(invLhs) - logSigmaSq;

        // shuffling the SNP index for faster convergence
        vector<int> snpIndexVec = Gadget::shuffle_index(blockStart, blockEnd);

        //for(unsigned i = blockStart; i <= blockEnd; i++){
        for (unsigned t = 0; t < blockSize; t++) {
            unsigned i = snpIndexVec[t];
            SnpInfo *snp = blockInfo->snpInfoVec[i - blockStart];
            if (snp->skip) {
                valuesPtr[i] = 0.0;
                continue;
            }
            if (badSnps[i]) {
                valuesPtr[i] = 0.0;
                continue;
            }
            if (snp->fitAsFixedEffect) {
                continue;
            }
            float oldSample = valuesPtr[i];
            int qcol = eigenQColIndex(blockInfo, i);
            if (qcol < 0) continue;
            Ref<const VectorXf> Qi = Q.col(qcol);
            float rhs = (Qi.dot(wcorr) + oldSample)*invVareDn;
            float uhat = invLhs * rhs;
            float logDelta1 = 0.5*(logInvLhsMsigma + uhat*rhs) + logPi;
            float logDelta0 = logPiComp;
            float probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
            pipPtr[i] = probDelta1;
            
            if (urnd[i] < probDelta1) {
                valuesPtr[i] = uhat + nrnd[i]*sqrtf(invLhs);
                wcorr += Qi*(oldSample - valuesPtr[i]);
                what  += Qi* valuesPtr[i];
                ssq[blk] += (valuesPtr[i] * valuesPtr[i]);
                s2pq[blk] += snp2pq[i];
                ++nnz[blk];
            } else {
                if (oldSample) wcorr += Qi * oldSample;
                valuesPtr[i] = 0.0;
            }
        }
    }
    
    sumSq = 0.0;
    sum2pq = 0.0;
    numNonZeros = 0.0;
    nnzPerBlk.setZero(nBlocks);
    for (unsigned blk=0; blk<nBlocks; ++blk) {
        sumSq += ssq[blk];
        sum2pq += s2pq[blk];
        numNonZeros += nnz[blk];
        nnzPerBlk[blk] = nnz[blk];
        ssqBlocks[blk] = ssq[blk];
    }
    
    values = VectorXf::Map(valuesPtr, size);
    pip = VectorXf::Map(pipPtr, size);
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

void ApproxBayesC::SnpEffects::computeFromBLUP(vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
                                            const vector<LDBlockInfo*> keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks,
                                            const float sigmaSq, const float pi, const float varg, const VectorXf &snp2pq){
    
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
      
    float logPi = log(pi);
    float logPiComp = log(1.0-pi);
    float logSigmaSq = log(sigmaSq);
    float invSigmaSq = 1.0f/sigmaSq;

    
#pragma omp parallel for schedule(dynamic)
    for(unsigned blk = 0; blk < nBlocks; blk++){
        Ref<const MatrixXf> Q = Qblocks[blk];
        Ref<VectorXf> wcorr = wcorrBlocks[blk];
        Ref<VectorXf> what = whatBlocks[blk];

        what.setZero();
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[blk];
        
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd   = blockInfo->endSnpIdx;
        
        float invVareDn = nGWASblocks[blk] / vareBlocks[blk];

        float invLhs = 1.0/(invVareDn + invSigmaSq);
        float logInvLhsMsigma = logf(invLhs) - logSigmaSq;

        for(unsigned i = blockStart; i <= blockEnd; i++){
            if (badSnps[i]) {
                valuesPtr[i] = 0.0;
                continue;
            }
            int qcol = eigenQColIndex(blockInfo, i);
            if (qcol < 0) continue;
            float oldSample = valuesPtr[i];
            Ref<const VectorXf> Qi = Q.col(qcol);
            float rhs = (Qi.dot(wcorr) + oldSample)*invVareDn;
            float uhat = invLhs * rhs;
            
//            if (i < 5) cout << "i " << i << " rhs " << rhs << " invLhs " << invLhs << " invVareDn " << invVareDn << " invSigmaSq " << invSigmaSq << " Qi.dot(wcorr) " << Qi.dot(wcorr) << " nGWASblocks[blk] " << nGWASblocks[blk] << endl;
            
//            if (i==0) {
//                cout << "\nQi\n" << Qi.transpose() << endl;
//                cout << "\nwcorr\n" << wcorr.transpose() << endl;
//            }
            
            float logDelta1 = 0.5*(logInvLhsMsigma + uhat*rhs) + logPi;
            float logDelta0 = logPiComp;
            float probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
            
            valuesPtr[i] = uhat;
            wcorr += Qi*(oldSample - valuesPtr[i]);
            what  += Qi* valuesPtr[i];
            ssq[blk] += (valuesPtr[i] * valuesPtr[i]);
            s2pq[blk] += snp2pq[i];
            ++nnz[blk];
        }
    }
    
    sumSq = 0.0;
    sum2pq = 0.0;
    numNonZeros = 0.0;
    nnzPerBlk.setZero(nBlocks);
    for (unsigned blk=0; blk<nBlocks; ++blk) {
        sumSq += ssq[blk];
        sum2pq += s2pq[blk];
        numNonZeros += nnz[blk];
        nnzPerBlk[blk] = nnz[blk];
        ssqBlocks[blk] = ssq[blk];
    }
    
    values = VectorXf::Map(valuesPtr, size);
}

void ApproxBayesC::VarEffects::computeRobustMode(const float varg, const VectorXf &snp2pq, const float pi, const bool noscale){
    if (noscale) {
        value = varg/(snp2pq.array().sum()*pi);
    } else {
        unsigned n = 0;
        for (unsigned i = 0; i < (unsigned)snp2pq.size(); ++i)
            if (snp2pq[i] > 0.f) ++n;
        if (!n) n = 1;
        value = varg/(n*pi);  // LDpred2's parameterisation (masked 2pq: count random SNPs)
    }
}

void ApproxBayesC::ResidualVar::sampleFromFC(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr, string &message) {
    float sse = ypy - effects.dot(ZPy) - effects.dot(rcorr);
    if (sse < 0) {
        sse = ypy;
        ++nNegVal;
        if (nNegVal > 10) {
            value = sse/nobs;
            message = "Negative residual variance";
            cout << message << ": " << value << endl;
            return;
        }
    }
    float dfTilde = df + nobs;
    float scaleTilde = sse + df*scale;
    value = InvChiSq::sample(dfTilde, scaleTilde);
}

void ApproxBayesC::GenotypicVar::compute(const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr){
    float modelSS = effects.dot(ZPy) - effects.dot(rcorr);
    if (modelSS < 0) modelSS = 0; // -modelSS;
    value = modelSS/nobs;
}

void ApproxBayesC::BlockGenotypicVar::compute(const vector<VectorXf> &whatBlocks){
    for (unsigned i=0; i<numBlocks; ++i) {
        values[i] = whatBlocks[i].squaredNorm();
        //cout << "varg " << i << " " << values[i] << endl;
    }
    total = values.sum();
}

void ApproxBayesC::BlockResidualVar::sampleFromFC(vector<VectorXf> &wcorrBlocks, VectorXf &vargBlocks, VectorXf &ssqBlocks, const VectorXf &nGWASblocks, const VectorXf &numEigenvalBlock){
    for (unsigned i=0; i<numBlocks; ++i) {
        float sse = wcorrBlocks[i].squaredNorm() * nGWASblocks[i];
        float dfTilde = df + numEigenvalBlock[i];
        float scaleTilde = sse + df*scale;
        float sample = InvChiSq::sample(dfTilde, scaleTilde);
        if (ssqBlocks[i]/vargBlocks[i] > threshold & sample/vary > 0.9) {
            values[i] = sample;
        } else {
            values[i] = vary;
        }
        //cout << "vare " << i << " " << nGWASblocks[i] << " " << sse << " " << dfTilde << " " << values[i] << " " << ssqBlocks[i] << " " << vargBlocks[i] << endl;
    }
    mean = values.mean();
}

void ApproxBayesC::BlockResidualVar::sampleFromFC(const vector<VectorXf> &wcorrBlocks, const VectorXf &beta, const VectorXf &b, const VectorXf &nGWASblocks, const vector<LDBlockInfo*> keptLdBlockInfoVec){
    for (unsigned i=0; i<numBlocks; ++i) {
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[i];
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd   = blockInfo->endSnpIdx;
        unsigned blockSize  = blockEnd - blockStart + 1;
        VectorXf betai = beta.segment(blockStart, blockSize);
        VectorXf bi = b.segment(blockStart, blockSize);
        float sse = nGWASblocks[i] - betai.dot(bi) - betai.dot(wcorrBlocks[i]);
        float dfTilde = df + nGWASblocks[i];
        float scaleTilde = sse + df*scale;
        values[i] = InvChiSq::sample(dfTilde, scaleTilde);
    }
    mean = values.mean();
}

void ApproxBayesC::Rounding::computeRcorr_sparse(const VectorXf &ZPy, const vector<SparseVector<float> > &ZPZ,
                                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                          const VectorXf &snpEffects, VectorXf &rcorr){
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

void ApproxBayesC::Rounding::computeRcorr_full(const VectorXf &ZPy, const vector<VectorXf> &ZPZ,
                                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                          const VectorXf &snpEffects, VectorXf &rcorr){
    VectorXf rcorrOld = rcorr;
    rcorr = ZPy;
#pragma omp parallel for
    for (unsigned chr=0; chr<chromInfoVec.size(); ++chr) {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            if (snpEffects[i]) rcorr.segment(windStart[i], windSize[i]) -= ZPZ[i]*snpEffects[i];
        }
    }
    value = sqrt(Gadget::calcVariance(rcorrOld-rcorr));
}

void ApproxBayesC::Rounding::computeWcorr_eigen(const vector<VectorXf> &wBlocks, const vector<MatrixXf> &Qblocks,
                                                const vector<LDBlockInfo*> keptLdBlockInfoVec,
                                                const VectorXf &snpEffects, vector<VectorXf> &wcorrBlocks){
    long nBlocks = keptLdBlockInfoVec.size();
    VectorXf res(nBlocks);
    
#pragma omp parallel for
    for(unsigned blk = 0; blk < nBlocks; blk++){
        Ref<const MatrixXf> Q = Qblocks[blk];
        Ref<VectorXf> wcorr = wcorrBlocks[blk];
        
        VectorXf wcorrOld = wcorr;
        wcorr = wBlocks[blk];

        LDBlockInfo *blockInfo = keptLdBlockInfoVec[blk];
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd   = blockInfo->endSnpIdx;

        for(unsigned i = blockStart; i <= blockEnd; i++){
            int qcol = eigenQColIndex(blockInfo, i);
            if (qcol < 0) continue;
            Ref<const VectorXf> Qi = Q.col(qcol);
            if (snpEffects[i]) wcorr -= Qi*snpEffects[i];
        }
        res[blk] = sqrt(Gadget::calcVariance(wcorrOld-wcorr));
    }
    value = res.sum();
}

void ApproxBayesC::Rounding::computeGhat(const MatrixXf &Z, const VectorXf &snpEffects, VectorXf &ghat){
    VectorXf ghatOld = ghat;
    ghat.setZero(ghat.size());
    for (unsigned i=0; i<snpEffects.size(); ++i) {
        if (snpEffects[i]) ghat += Z.col(i)*snpEffects[i];
    }
    value = sqrt(Gadget::calcVariance(ghatOld-ghat));
}

//void ApproxBayesC::PopulationStratification::compute(const VectorXf &rcorr, const VectorXf &ZPZdiag, const VectorXf &LDsamplVar, const float varg, const float vare, const VectorXf &chisq){
//    
//    value = (rcorr.array().square()/(ZPZdiag.array() * (LDsamplVar.array()*varg + value + vare))).mean() - 1.0;
//    value = value < -0.01 ? -0.01 : value;
//    
////    VectorXf varEta = ZPZdiag.array() * (LDsamplVar.array()*varg + value + vare);
//////    VectorXf wt = varEta.array().square().inverse();
////    VectorXf wt = (rcorr.array().square()/ZPZdiag.array().square()).square().inverse();
//////    VectorXf zsq = rcorr.array().square()/varEta.array();
////    VectorXf zsq = rcorr.array().square()/ZPZdiag.array() - LDsamplVar.array()*varg - vare;
////    value = zsq.cwiseProduct(wt).sum()/wt.sum();
//
//    
////    VectorXf tmp = rcorr.array().square()/ZPZdiag.array() - LDsamplVar.array()*varg - vare;
////    float ssq = 0.0;
////    long size = rcorr.size();
////    long cnt = 0;
////    for (unsigned i=0; i<size; ++i) {
//////        if (chisq[i] < 20 && !(i%20)) {
////            ssq += tmp[i];
////            ++cnt;
//////        }
////    }
//////    ssq /= float(cnt);
////    if (ssq < 0) ssq = 0.0;
////    float dfTilde = df + cnt;
////    float scaleTilde = ssq + df*scale;
////    value = InvChiSq::sample(dfTilde, scaleTilde);
//
//    
////        ofstream out("rcorr.txt");
////        out << rcorr.array().square()/(ZPZdiag.array() * (LDsamplVar.array()*varg + value + vare)) << endl;
////        out.close();
//    
//}
//
//void ApproxBayesC::PopulationStratification::compute(const VectorXf &rcorr, const VectorXf &ZPZdiag, const VectorXf &LDsamplVar, const float varg, const float vare, const vector<ChromInfo*> chromInfoVec){
//    
//    for (unsigned i=0; i<22; ++i) {
//        unsigned start = chromInfoVec[i]->startSnpIdx;
//        unsigned end = chromInfoVec[i]->endSnpIdx;
//        unsigned size = end - start + 1;
//        chrSpecific[i] = (rcorr.segment(start,size).array().square()/(ZPZdiag.segment(start,size).array() * (LDsamplVar.segment(start,size).array()*varg + chrSpecific[i] + vare))).mean() - 1.0;
//    }    
//}
//
//void ApproxBayesC::NumResidualOutlier::compute(const VectorXf &rcorr, const VectorXf &ZPZdiag, const VectorXf &LDsamplVar, const float varg, const float vare, const vector<string> &snpName, VectorXi &leaveout, const vector<SparseVector<float> > &ZPZ, const VectorXf &ZPy, const VectorXf &snpEffects) {
//    //if (iter<10) return;
//    
//    VectorXf tss = ZPy.array().square()/ZPZdiag.array();
//    VectorXf sse = rcorr.array().square()/ZPZdiag.array();
//    value = 0;
//    long size = tss.size();
//    for (unsigned i=0; i<size; ++i) {
//        if (sse[i] > 10 && tss[i] < 10) ++value;
//        //if (sse[i] > 30 && sse[i] > tss[i]) ++value;
//    }
//    
////    MatrixXf X(tss.size(), 2);
////    X.col(0) = VectorXf::Ones(tss.size());
////    X.col(1) = snpEffects;
////    VectorXf bhat = ZPy.array()/ZPZdiag.array();
////    VectorXf b = X.householderQr().solve(bhat);
////    value = b[1];
//    
////    cout << sse.array().maxCoeff() << endl;
//    
////    cout << tss.mean() << " " << sse.mean() << endl;
//    
////    value = 0;
////    
////    VectorXf tmp = rcorr.array().square()/(ZPZdiag.array() * (LDsamplVar.array()*varg + value + vare));
////    //VectorXf tmp = rcorr.array().square()/ZPZdiag.array() - LDsamplVar.array()*varg;
////    
////    long size = tmp.size();
////    stringstream ss;
////    for (unsigned i=0; i<size; ++i) {
////        if (tmp[i]>20) {
////            ++value;
////            //leaveout[i] = 1;
////            ss << " " << snpName[i];
////            
//////            for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {
//////                leaveout[it.index()] = 1;
//////            }
////            
////        }
////    }
////    if (value) out << iter << ss.str() << endl;
//}
//
//void ApproxBayesC::ldScoreReg(const VectorXf &chisq, const VectorXf &LDscore, const VectorXf &LDsamplVar,
//                              const float varg, const float vare, float &ps){
//    
//    long nrow = chisq.size();
//    
////    ps = (chisq - LDscore - LDsamplVar*varg - VectorXf::Ones(nrow)*vare).mean();
////    ps = (chisq - LDscore*0.17/float(nrow) - VectorXf::Ones(nrow)).mean();
//
////    return;
//    
//    VectorXf y = chisq - LDsamplVar*varg - VectorXf::Ones(nrow)*vare;
////    VectorXf y = chisq;
////    VectorXf weight = 2.0*(LDscore*vargj + LDsamplVar*varg + VectorXf::Ones(nrow)*vare).array().square();
////    VectorXf weight = 2.0*(LDscore*vargj + VectorXf::Ones(nrow)).array().square();
////    VectorXf weightInv = weight.cwiseInverse();
//    
//    MatrixXf X(nrow, 2);
////    X.col(0) = weight;
////    X.col(1) = LDscore.cwiseProduct(weight);
////    y.array() *= weight.array();
//
//    X.col(0) = VectorXf::Ones(nrow);
//    X.col(1) = LDscore;
//    
////    unsigned m = 0;
////    for (unsigned i=0; i<nrow; ++i) {
////        if (chisq[i] < 30) ++m;
////    }
////    VectorXf ysub(m);
////    MatrixXf Xsub(m, 2);
////    unsigned j=0;
////    for (unsigned i=0; i<nrow; ++i) {
////        if (chisq[i] < 30) {
////            ysub[j] = y[i];
////            Xsub.row(j) = X.row(i);
////            ++j;
////        }
////    }
////    VectorXf b = Xsub.householderQr().solve(ysub);
//
//    
//    VectorXf b = X.householderQr().solve(y);
////    VectorXf b = (X.transpose()*weightInv.asDiagonal()*X).inverse()*X.transpose()*weightInv.asDiagonal()*y;
//
//    ps = b[0];
//    
////    cout << b.transpose() << endl;
//    
//}
//
//void ApproxBayesC::InterChrGenetCov::compute(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr) {
//    if (!spouseCorrelation) return;
//    float bZPy = effects.dot(ZPy);
//    float brcorr = effects.dot(rcorr);
//    float varg = (bZPy - brcorr)/nobs;
////    float vare = (ypy - bZPy - brcorr)/nobs;
//    float varp = ypy/nobs;
//    float hsq = varg/varp;
//    float R = spouseCorrelation*hsq / (1 - spouseCorrelation*hsq);
//    value = varg * R; // * 0.95;
//}
//
//void ApproxBayesC::NnzGwas::compute(const VectorXf &effects, const vector<SparseVector<float> > &ZPZ, const VectorXf &ZPZdiag) {
//    value = 0;
//    long numSnps = effects.size();
//    unsigned i, j;
//    for (i=0; i<numSnps; ++i) {
//        for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {
//            j = it.index();
//            if (effects[j]){
//                if (it.value()*it.value() > 0.1*ZPZdiag[i]*ZPZdiag[j]) {
//                    ++value;
//                    break;
//                }
//            }
//        }
//    }
//}
//
//void ApproxBayesC::PiGwas::compute(const float nnzGwas, const unsigned int numSnps) {
//    value = nnzGwas/float(numSnps);
//}
//
//void ApproxBayesC::checkHsq(vector<float> &hsqMCMC) {
//    long niter = hsqMCMC.size();
//    VectorXf y(niter);
//    MatrixXf X(niter, 2);
//    for (unsigned i=0; i<niter; ++i) {
//        y[i] = hsqMCMC[i];
//        X(i,0) = 1;
//        X(i,1) = i;
//    }
//    VectorXf b = X.householderQr().solve(y);
//    float slope = b[1];
//    float vare = (y.squaredNorm() - (X*b).squaredNorm())/float(niter);
//    float se = sqrt(vare/X.col(1).squaredNorm());
//    if ((slope/se) > 3) {   // 3 corresponds to P = 0.001 at one-way test
//        string slope_str = to_string(static_cast<float>(slope));
//        string se_str = to_string(static_cast<float>(se));
//        throw("\nError: The SNP-heritability is increasing over MCMC iterations (slope: " + slope_str + "; se: " + se_str + "). This may indicate that effect sizes are \"blowing up\" likely due to a convergence problem.");
//    }
//}

void ApproxBayesC::NumBadSnps::compute_sparse(VectorXi &badSnps, VectorXf &effects, VectorXf &effectMean, const VectorXf &b, VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const vector<ChromInfo*> &chromInfoVec, const int iter, const Data &data) {
    //cout << "computing NumBadSnps..." << endl;
    
    value = 0;
    
    float rate_thresh1 = 0, rate_thresh2 = 0;
    if(iter < 300){
        rate_thresh1 = 4.0;
        rate_thresh2 = 2.0;
    }else if(iter < 600){
        rate_thresh1 = 3.0;
        rate_thresh2 = 1.5;
    }else if(iter < 900){
        rate_thresh1 = 2.0;
        rate_thresh2 = 1.3;
    }else{
        rate_thresh1 = 1.5;
        rate_thresh2 = 1.1;
    }
    
    long numChr = chromInfoVec.size();
    for (unsigned chr=0; chr<numChr; ++chr) {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;

        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            if (badSnps[i]) {
                continue;
            }
            if (data.incdSnpInfoVec[i]->fitAsFixedEffect) {
                continue;
            }

            float rate_b = abs((effectMean[i] - b[i])/b[i]);
            bool sameSign = (b[i] >= 0.0f) == (effectMean[i] >= 0.0f);
            float compare_rate = sameSign ? rate_thresh1 : rate_thresh2;

            //cout << "SNP: " << i << "\t" << effectMean[i] << "\t" << b[i] << "\t" << rate_b << "\t" << compare_rate << endl;

            if(abs(effectMean[i]) > betaThresh && rate_b > compare_rate){
                //cout << "DEL:" << "\t" << betaVal << "\t" << b[idx] << "\t" << rate_b << "\t" << compare_rate << std::endl;
                for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                    rcorr[it.index()] += it.value() * effects[i];
                }
                
                effects[i] = 0.0;
                effectMean[i] = 0.0;
                badSnps[i] = 1;
                badSnpIdx.push_back(i);
                badSnpName.push_back(snpNames[i]);
                if (writeTxt) out << i+1 << "\t" << snpNames[i] << endl;
                ++value;
            }
        }
    }
}

void ApproxBayesC::NumBadSnps::compute_full(VectorXi &badSnps, VectorXf &effects, VectorXf &effectMean, const VectorXf &b, VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec, const int iter, const Data &data) {
    //cout << "computing NumBadSnps..." << endl;
    
    value = 0;
    
    float rate_thresh1 = 0, rate_thresh2 = 0;
    if(iter < 300){
        rate_thresh1 = 4.0;
        rate_thresh2 = 2.0;
    }else if(iter < 600){
        rate_thresh1 = 3.0;
        rate_thresh2 = 1.5;
    }else if(iter < 900){
        rate_thresh1 = 2.0;
        rate_thresh2 = 1.3;
    }else{
        rate_thresh1 = 1.5;
        rate_thresh2 = 1.1;
    }
    
    long numChr = chromInfoVec.size();
    for (unsigned chr=0; chr<numChr; ++chr) {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;

        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            if (badSnps[i]) {
                continue;
            }
            if (data.incdSnpInfoVec[i]->fitAsFixedEffect) {
                continue;
            }

            float rate_b = abs((effectMean[i] - b[i])/b[i]);
            bool sameSign = (b[i] >= 0.0f) == (effectMean[i] >= 0.0f);
            float compare_rate = sameSign ? rate_thresh1 : rate_thresh2;

            //cout << "SNP: " << i << "\t" << effectMean[i] << "\t" << b[i] << "\t" << rate_b << "\t" << compare_rate << endl;

            if(abs(effectMean[i]) > betaThresh && rate_b > compare_rate){
                //cout << "DEL:" << "\t" << betaVal << "\t" << b[idx] << "\t" << rate_b << "\t" << compare_rate << std::endl;
                rcorr.segment(windStart[i], windSize[i]) += ZPZ[i] * effects[i];
                
                effects[i] = 0.0;
                effectMean[i] = 0.0;
                badSnps[i] = 1;
                badSnpIdx.push_back(i);
                badSnpName.push_back(snpNames[i]);
                if (writeTxt) out << i+1 << "\t" << snpNames[i] << endl;
                ++value;
            }
        }
    }
}

void ApproxBayesC::NumBadSnps::compute_eigen(VectorXi &badSnps, VectorXf &effects, VectorXf &effectMean, const VectorXf &b, vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, const vector<LDBlockInfo*> keptLdBlockInfoVec, const int iter, const Data &data) {
    //cout << "computing NumBadSnps..." << endl;
    
    value = 0;
    
    float rate_thresh1 = 0, rate_thresh2 = 0;
    if(iter < 300){
        rate_thresh1 = 4.0;
        rate_thresh2 = 2.0;
    }else if(iter < 600){
        rate_thresh1 = 3.0;
        rate_thresh2 = 1.5;
    }else if(iter < 900){
        rate_thresh1 = 2.0;
        rate_thresh2 = 1.3;
    }else{
        rate_thresh1 = 1.5;
        rate_thresh2 = 1.1;
    }
    
    unsigned nBlocks = Qblocks.size();
    for(unsigned blk = 0; blk < nBlocks; blk++){
        Ref<const MatrixXf> Q = Qblocks[blk];
        Ref<VectorXf> wcorr = wcorrBlocks[blk];
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[blk];
        
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd   = blockInfo->endSnpIdx;
        
        for(unsigned i = blockStart; i <= blockEnd; i++){
            if (badSnps[i]) {
                continue;
            }
            if (data.incdSnpInfoVec[i]->fitAsFixedEffect) {
                continue;
            }

            float rate_b = abs((effectMean[i] - b[i])/b[i]);
            bool sameSign = (b[i] >= 0.0f) == (effectMean[i] >= 0.0f);
            float compare_rate = sameSign ? rate_thresh1 : rate_thresh2;

            //cout << "SNP: " << i << "\t" << effectMean[i] << "\t" << b[i] << "\t" << rate_b << "\t" << compare_rate << endl;

            if(abs(effectMean[i]) > betaThresh && rate_b > compare_rate){
                //cout << "DEL:" << "\t" << betaVal << "\t" << b[idx] << "\t" << rate_b << "\t" << compare_rate << std::endl;
                int qcol = eigenQColIndex(blockInfo, i);
                if (qcol < 0) continue;
                Ref<const VectorXf> Qi = Q.col(qcol);
                wcorr = wcorr + Qi * effects[i];
                effects[i] = 0.0;
                effectMean[i] = 0.0;
                badSnps[i] = 1;
                badSnpIdx.push_back(i);
                badSnpName.push_back(snpNames[i]);
                if (writeTxt) out << i+1 << "\t" << snpNames[i] << endl;
                ++value;
            }
        }
    }
}

void ApproxBayesC::sampleUnknowns(const unsigned iter){
    if (lowRankModel) {
        fixedSnpEffects.sampleFromFC_eigen(snpEffects, wcorrBlocks, data.Qblocks, whatBlocks, data.keptLdBlockInfoVec, data.nGWASblock, vareBlk.values, data.numFixedEffectSnps);
        snpEffects.sampleFromFC_eigen(wcorrBlocks, data.Qblocks, whatBlocks, data.keptLdBlockInfoVec, data.nGWASblock, vareBlk.values, sigmaSq.value, pi.value, vargRandom.value, data.snp2pq);
    } else if (sparse) {
        fixedSnpEffects.sampleFromFC_sparse(snpEffects, rcorr, data.ZPZsp, data.ZPZdiag, data.snp2pq, data.incdSnpInfoVec, varg.value, vare.value, data.numFixedEffectSnps);
        snpEffects.sampleFromFC_sparse(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.chromInfoVec, data.snp2pq, sigmaSq.value, pi.value, vare.value, vargRandom.value, data.incdSnpInfoVec);
    } else {
        fixedSnpEffects.sampleFromFC_full(snpEffects, rcorr, data.ZPZ, data.ZPZdiag, data.windStart, data.windSize, data.snp2pq, data.incdSnpInfoVec, varg.value, vare.value, data.numFixedEffectSnps);
        snpEffects.sampleFromFC_full(rcorr, data.ZPZ, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, data.snp2pq, sigmaSq.value, pi.value, vare.value, vargRandom.value, data.incdSnpInfoVec);
    }
    snpEffects.computePosteriorMean(iter);
    snpPip.getValues(snpEffects.pip);
    nnzSnp.getValue(snpEffects.numNonZeros);
    
    if (robustMode) {
        sigmaSq.computeRobustMode(vargRandom.value, data.numFixedEffectSnps ? data.snp2pqForRandomSnpEffects : data.snp2pq, pi.value, noscale);
    } else {
        sigmaSq.sampleFromFC(snpEffects.sumSq, snpEffects.numNonZeros);
    }

    if (estimatePi) pi.sampleFromFC(data.numIncdSnps, snpEffects.numNonZeros);
    
    if (lowRankModel) {
        vargBlk.compute(whatBlocks);
        vareBlk.sampleFromFC(wcorrBlocks, vargBlk.values, snpEffects.ssqBlocks, data.nGWASblock, data.numEigenvalBlock);
        //vareBlk.sampleFromFC(wcorrBlocks, snpEffects.values, data.b, data.nGWASblock, data.keptLdBlockInfoVec);
        vargRandom.value = vargBlk.total;
        vare.value = vareBlk.mean;
    } else {
        varg.compute(snpEffects.values, data.ZPy, rcorr);
        vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, status);
        if (data.numFixedEffectSnps && data.varGenotypic > 0.f)
            vargRandom.value = varg.value * (data.vargRandomSnpEffects / data.varGenotypic);
        else
            vargRandom.value = varg.value;
    }
    hsqRandom.compute(vargRandom.value, data.varPhenotypic);
    if (data.numFixedEffectSnps) {
        hsqFixed.compute(snpEffects.values, data);
        hsq.value = hsqRandom.value + hsqFixed.value;
        if (lowRankModel)
            varg.value = vargRandom.value + hsqFixed.value * data.varPhenotypic;
    } else {
        hsq.value = hsqRandom.value;
        if (lowRankModel)
            varg.value = vargRandom.value;
    }
    
    if (!(iter % 10)) {
        if (lowRankModel) {
            nBadSnps.compute_eigen(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, wcorrBlocks, data.Qblocks, data.keptLdBlockInfoVec, iter, data);
        } else if (sparse) {
            nBadSnps.compute_sparse(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, rcorr, data.ZPZsp, data.chromInfoVec, iter, data);
        } else {
            nBadSnps.compute_full(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, rcorr, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, iter, data);
        }
    }

    if (!(iter % 100)) {
        if (lowRankModel) {
            rounding.computeWcorr_eigen(data.wcorrBlocks, data.Qblocks, data.keptLdBlockInfoVec, snpEffects.values, wcorrBlocks);
        } else if (sparse) {
            rounding.computeRcorr_sparse(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
        } else {
            rounding.computeRcorr_full(data.ZPy, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
        }
    }
}


void ApproxBayesB::SnpEffects::sampleFromFC_sparse(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                            const vector<ChromInfo*> &chromInfoVec,
                                            const VectorXf &snp2pq,
                                            const VectorXf &sigmaSq, const float pi, const float vare, const float varg) {
    long numChr = chromInfoVec.size();
    
    float ssq[numChr], s2pq[numChr], nnz[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(s2pq,0,sizeof(float)*numChr);
    memset(nnz,0, sizeof(float)*numChr);
        
    pip.setZero(size);
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads
    float *pipPtr = pip.data();

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
#pragma omp parallel for schedule(dynamic)
    for (unsigned chr=0; chr<numChr; ++chr) {
        //cout << " thread " << omp_get_thread_num() << " chr " << chr << endl;
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned chrSize = chrEnd - chrStart + 1;

        float oldSample;
        float rhs, invLhs, uhat;
        float logDelta0, logDelta1, probDelta1;
        float logPi = log(pi);
        float logPiComp = log(1.0-pi);
        float varei = varg + vare;
        float sampleDiff;

        // shuffling the SNP index for faster convergence
        vector<int> snpIndexVec = Gadget::shuffle_index(chrStart, chrEnd);

        unsigned i;
        for (unsigned t = 0; t < chrSize; t++) {
            i = snpIndexVec[t];
            if (badSnps[i]) {
                valuesPtr[i] = 0.0;
                continue;
            }
            oldSample = valuesPtr[i];
            rhs = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            invLhs = 1.0f/(ZPZdiag[i]/varei + 1.0f/sigmaSq[i]);
            uhat = invLhs*rhs;
            logDelta1 = 0.5*(logf(invLhs) - logf(sigmaSq[i]) + uhat*rhs) + logPi;
            logDelta0 = logPiComp;
            probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
            pipPtr[i] = probDelta1;
            
            //if (bernoulli.sample(probDelta1)) {
            if (urnd[i] < probDelta1) {
                //valuesPtr[i] = normal.sample(uhat, invLhs);
                valuesPtr[i] = uhat + nrnd[i]*sqrtf(invLhs);
                sampleDiff = oldSample - valuesPtr[i];
                for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {
                    rcorr[it.index()] += it.value() * sampleDiff;
                }
                betaSq[i] = valuesPtr[i]*valuesPtr[i];
                ssq[chr]  += valuesPtr[i]*valuesPtr[i];
                s2pq[chr] += snp2pq[i];
                ++nnz[chr];
            } else {
                if (oldSample) {
                    //rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*oldSample;
                    for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {
                        rcorr[it.index()] += it.value() * oldSample;
                    }
                }
                valuesPtr[i] = 0.0;
            }
        }
    }
        
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
    pip = VectorXf::Map(pipPtr, size);
}

void ApproxBayesB::SnpEffects::sampleFromFC_full(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const VectorXf &snp2pq,
                                            const VectorXf &sigmaSq, const float pi, const float vare, const float varg) {
    long numChr = chromInfoVec.size();
    
    float ssq[numChr], nnz[numChr], s2pq[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(nnz,0,sizeof(float)*numChr);
    memset(s2pq,0,sizeof(float)*numChr);
    
    pip.setZero(size);
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads
    float *pipPtr = pip.data();

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
#pragma omp parallel for
    for (unsigned chr=0; chr<numChr; ++chr) {
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned chrSize = chrEnd - chrStart + 1;

        float oldSample;
        float rhs, invLhs, uhat;
        float logDelta0, logDelta1, probDelta1;
        float logPi = log(pi);
        float logPiComp = log(1.0-pi);
//        float logSigmaSq = log(sigmaSq);
//        float invSigmaSq = 1.0f/sigmaSq;
        float varei = varg + vare;
        float sampleDiff;

        // shuffling the SNP index for faster convergence
        vector<int> snpIndexVec = Gadget::shuffle_index(chrStart, chrEnd);

        unsigned i;
        for (unsigned t = 0; t < chrSize; t++) {
            i = snpIndexVec[t];
            if (badSnps[i]) {
                valuesPtr[i] = 0.0;
                continue;
            }
            oldSample = valuesPtr[i];
            rhs = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            invLhs = 1.0f/(ZPZdiag[i]/varei + 1.0f/sigmaSq[i]);
            uhat = invLhs*rhs;
            logDelta1 = 0.5*(logf(invLhs) - logf(sigmaSq[i]) + uhat*rhs) + logPi;
            logDelta0 = logPiComp;
            probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
            pipPtr[i] = probDelta1;
            
            //if (bernoulli.sample(probDelta1)) {
            if (urnd[i] < probDelta1) {
                //valuesPtr[i] = normal.sample(uhat, invLhs);
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
    pip = VectorXf::Map(pipPtr, size);
}

void ApproxBayesB::SnpEffects::sampleFromFC_eigen(vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
                                            const vector<LDBlockInfo*> keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks,
                                            const VectorXf &sigmaSq, const float pi, const float varg, const VectorXf &snp2pq){
    
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

    pip.setZero(size);
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads
    float *pipPtr = pip.data();

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
    float logPi = log(pi);
    float logPiComp = log(1.0-pi);
    
#pragma omp parallel for schedule(dynamic)
    for(unsigned blk = 0; blk < nBlocks; blk++){
        Ref<const MatrixXf> Q = Qblocks[blk];
        Ref<VectorXf> wcorr = wcorrBlocks[blk];
        Ref<VectorXf> what = whatBlocks[blk];

        what.setZero();
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[blk];
        
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd   = blockInfo->endSnpIdx;
        unsigned blockSize  = blockEnd - blockStart + 1;

        float invVareDn = nGWASblocks[blk] / vareBlocks[blk];

        // shuffling the SNP index for faster convergence
        vector<int> snpIndexVec = Gadget::shuffle_index(blockStart, blockEnd);

        //for(unsigned i = blockStart; i <= blockEnd; i++){
        for (unsigned t = 0; t < blockSize; t++) {
            unsigned i = snpIndexVec[t];
            if (badSnps[i]) {
                valuesPtr[i] = 0.0;
                continue;
            }
            int qcolB = eigenQColIndex(blockInfo, i);
            if (qcolB < 0) continue;
            float oldSample = valuesPtr[i];
            Ref<const VectorXf> Qi = Q.col(qcolB);
            float rhs = (Qi.dot(wcorr) + oldSample)*invVareDn;
            float invLhs = 1.0/(invVareDn + 1.0f/sigmaSq[i]);
            float uhat = invLhs * rhs;
            float logInvLhsMsigma = logf(invLhs) - logf(sigmaSq[i]);
            float logDelta1 = 0.5*(logInvLhsMsigma + uhat*rhs) + logPi;
            float logDelta0 = logPiComp;
            float probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
            pipPtr[i] = probDelta1;

//            if (bernoulli.sample(probDelta1)) {
            if (urnd[i] < probDelta1) {
//                valuesPtr[i] = normal.sample(uhat, invLhs);
                valuesPtr[i] = uhat + nrnd[i]*sqrtf(invLhs);
                wcorr += Qi*(oldSample - valuesPtr[i]);
                what  += Qi* valuesPtr[i];
                ssq[blk] += (valuesPtr[i] * valuesPtr[i]);
                s2pq[blk] += snp2pq[i];
                ++nnz[blk];
            } else {
                if (oldSample) wcorr += Qi * oldSample;
                valuesPtr[i] = 0.0;
            }
        }
    }
    
    sumSq = 0.0;
    sum2pq = 0.0;
    numNonZeros = 0.0;
    nnzPerBlk.setZero(nBlocks);
    for (unsigned blk=0; blk<nBlocks; ++blk) {
        sumSq += ssq[blk];
        sum2pq += s2pq[blk];
        numNonZeros += nnz[blk];
        nnzPerBlk[blk] = nnz[blk];
        ssqBlocks[blk] = ssq[blk];
    }
    
    values = VectorXf::Map(valuesPtr, size);
    pip = VectorXf::Map(pipPtr, size);
}


void ApproxBayesB::sampleUnknowns(const unsigned iter) {
    if (lowRankModel) {
        snpEffects.sampleFromFC_eigen(wcorrBlocks, data.Qblocks, whatBlocks, data.keptLdBlockInfoVec, data.nGWASblock, vareBlk.values, sigmaSq.values, pi.value, varg.value, data.snp2pq);
    } else if (sparse) {
        snpEffects.sampleFromFC_sparse(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.chromInfoVec, data.snp2pq, sigmaSq.values, pi.value, vare.value, varg.value);
    } else {
        snpEffects.sampleFromFC_full(rcorr, data.ZPZ, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, data.snp2pq, sigmaSq.values, pi.value, vare.value, varg.value);
    }
    snpEffects.computePosteriorMean(iter);
    snpPip.getValues(snpEffects.pip);
    nnzSnp.getValue(snpEffects.numNonZeros);

    sigmaSq.sampleFromFC(snpEffects.betaSq);

    if (estimatePi) pi.sampleFromFC(data.numIncdSnps, snpEffects.numNonZeros);
        
    if (lowRankModel) {
        vargBlk.compute(whatBlocks);
        vareBlk.sampleFromFC(wcorrBlocks, vargBlk.values, snpEffects.ssqBlocks, data.nGWASblock, data.numEigenvalBlock);
        //vareBlk.sampleFromFC(wcorrBlocks, snpEffects.values, data.b, data.nGWASblock, data.keptLdBlockInfoVec);
        varg.value = vargBlk.total;
        vare.value = vareBlk.mean;
    }
    else {
        varg.compute(snpEffects.values, data.ZPy, rcorr);
        vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, status);
    }
    //hsq.compute(varg.value, vare.value);
    hsq.value = varg.value / data.varPhenotypic;

    if (!(iter % 10)) {
        if (lowRankModel) {
            nBadSnps.compute_eigen(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, wcorrBlocks, data.Qblocks, data.keptLdBlockInfoVec, iter, data);
        } else if (sparse) {
            nBadSnps.compute_sparse(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, rcorr, data.ZPZsp, data.chromInfoVec, iter, data);
        } else {
            nBadSnps.compute_full(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, rcorr, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, iter, data);
        }
    }

    if (!(iter % 100)) {
        if (lowRankModel) {
            rounding.computeWcorr_eigen(data.wcorrBlocks, data.Qblocks, data.keptLdBlockInfoVec, snpEffects.values, wcorrBlocks);
        } else if (sparse) {
            rounding.computeRcorr_sparse(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
        } else {
            rounding.computeRcorr_full(data.ZPy, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
        }
    }
}


void ApproxBayesS::SnpEffects::sampleFromFC_sparse(VectorXf &rcorr,const vector<SparseVector<float> > &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                            const vector<ChromInfo*> &chromInfoVec,
                                            const float sigmaSq, const float pi, const float vare,
                                            const VectorXf &snp2pqPowS, const VectorXf &snp2pq, const float varg){
    long numChr = chromInfoVec.size();
    
    float ssq[numChr], nnz[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(nnz,0,sizeof(float)*numChr);
    
//    for (unsigned chr=0; chr<numChr; ++chr) {
//        ChromInfo *chromInfo = chromInfoVec[chr];
//        unsigned chrStart = chromInfo->startSnpIdx;
//        unsigned chrEnd   = chromInfo->endSnpIdx;
//        if (iter==0) {
//            cout << "chr " << chr+1 << " start " << chrStart << " end " << chrEnd << endl;
//        }
//    }
    
    pip.setZero(size);
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads
    float *pipPtr = pip.data();

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
#pragma omp parallel for schedule(dynamic)
    for (unsigned chr=0; chr<numChr; ++chr) {
        //cout << " thread " << omp_get_thread_num() << " chr " << chr << endl;
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned chrSize = chrEnd - chrStart + 1;

        float oldSample;
        float rhs, invLhs, uhat;
        float logDelta0, logDelta1, probDelta1;
        float logPi = log(pi);
        float logPiComp = log(1.0-pi);
        float invSigmaSq = 1.0f/sigmaSq;
        float varei = varg + vare;
        float sampleDiff;

        // shuffling the SNP index for faster convergence
        vector<int> snpIndexVec = Gadget::shuffle_index(chrStart, chrEnd);

        unsigned i;
        for (unsigned t = 0; t < chrSize; t++) {
            i = snpIndexVec[t];
            if (badSnps[i]) {
                valuesPtr[i] = 0.0;
                continue;
            }
            oldSample = valuesPtr[i];
            rhs  = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            invLhs = 1.0f/(ZPZdiag[i]/varei + invSigmaSq/snp2pqPowS[i]);
            uhat = invLhs*rhs;
            
            logDelta1 = 0.5*(logf(invLhs) - logf(snp2pqPowS[i]*sigmaSq) + uhat*rhs) + logPi;
            logDelta0 = logPiComp;
            
            probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
            pipPtr[i] = probDelta1;
            
            //if (bernoulli.sample(probDelta1)) {
            if (urnd[i] < probDelta1) {
                //valuesPtr[i] = normal.sample(uhat, invLhs);
                valuesPtr[i] = uhat + nrnd[i]*sqrtf(invLhs);
                sampleDiff = oldSample - valuesPtr[i];
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
    
    wtdSumSq = 0.0;
    numNonZeros = 0;
    nnzPerChr.setZero(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        wtdSumSq += ssq[i];
        numNonZeros += nnz[i];
        nnzPerChr[i] = nnz[i];
    }

    values = VectorXf::Map(valuesPtr, size);
    pip = VectorXf::Map(pipPtr, size);
}

void ApproxBayesS::SnpEffects::sampleFromFC_full(VectorXf &rcorr,const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const float sigmaSq, const float pi, const float vare,
                                            const VectorXf &snp2pqPowS, const VectorXf &snp2pq, const float varg){
    long numChr = chromInfoVec.size();
    
    float ssq[numChr], nnz[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(nnz,0,sizeof(float)*numChr);
    
    pip.setZero(size);
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads
    float *pipPtr = pip.data();

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }

#pragma omp parallel for schedule(dynamic)
    for (unsigned chr=0; chr<numChr; ++chr) {
        //cout << " thread " << omp_get_thread_num() << " chr " << chr << endl;
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned chrSize = chrEnd - chrStart + 1;

        float oldSample;
        float rhs, invLhs, uhat;
        float logDelta0, logDelta1, probDelta1;
        float logPi = log(pi);
        float logPiComp = log(1.0-pi);
        float invSigmaSq = 1.0f/sigmaSq;
        float varei = varg + vare;
        float sampleDiff;

        // shuffling the SNP index for faster convergence
        vector<int> snpIndexVec = Gadget::shuffle_index(chrStart, chrEnd);

        unsigned i;
        for (unsigned t = 0; t < chrSize; t++) {
            i = snpIndexVec[t];
            if (badSnps[i]) {
                valuesPtr[i] = 0.0;
                continue;
            }
            oldSample = valuesPtr[i];
            rhs  = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            invLhs = 1.0f/(ZPZdiag[i]/varei + invSigmaSq/snp2pqPowS[i]);
            uhat = invLhs*rhs;
            
            logDelta1 = 0.5*(logf(invLhs) - logf(snp2pqPowS[i]*sigmaSq) + uhat*rhs) + logPi;
            logDelta0 = logPiComp;
            
            probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
            pipPtr[i] = probDelta1;
            
            
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
    
    wtdSumSq = 0.0;
    numNonZeros = 0;
    nnzPerChr.setZero(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        wtdSumSq += ssq[i];
        numNonZeros += nnz[i];
        nnzPerChr[i] = nnz[i];
    }

    values = VectorXf::Map(valuesPtr, size);
    pip = VectorXf::Map(pipPtr, size);
}


void ApproxBayesS::SnpEffects::sampleFromFC_ind(const VectorXf &ZPy,const MatrixXf &Z, const VectorXf &ZPZdiag,
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

void ApproxBayesS::SnpEffects::sampleFromFC_eigen(vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
                               const vector<LDBlockInfo*> keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks,
                               const float sigmaSq, const float pi, const float varg,
                                const VectorXf &snp2pqPowS, const VectorXf &snp2pq){
    
    long nBlocks = keptLdBlockInfoVec.size();

    whatBlocks.resize(nBlocks);
    ssqBlocks.resize(nBlocks);
    for (unsigned i=0; i<nBlocks; ++i) {
        whatBlocks[i].resize(wcorrBlocks[i].size());
    }

    float ssq[nBlocks], nnz[nBlocks];
    memset(ssq,0, sizeof(float)*nBlocks);
    memset(nnz,0, sizeof(float)*nBlocks);

    pip.setZero(size);
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads
    float *pipPtr = pip.data();

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
    float logPi = log(pi);
    float logPiComp = log(1.0-pi);
    float invSigmaSq = 1.0f/sigmaSq;

    
#pragma omp parallel for schedule(dynamic)
    for(unsigned blk = 0; blk < nBlocks; blk++){
        Ref<const MatrixXf> Q = Qblocks[blk];
        Ref<VectorXf> wcorr = wcorrBlocks[blk];
        Ref<VectorXf> what = whatBlocks[blk];

        what.setZero();
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[blk];
        
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd   = blockInfo->endSnpIdx;
        unsigned blockSize  = blockEnd - blockStart + 1;

        float invVareDn = nGWASblocks[blk] / vareBlocks[blk];

        // shuffling the SNP index for faster convergence
        vector<int> snpIndexVec = Gadget::shuffle_index(blockStart, blockEnd);

        //for(unsigned i = blockStart; i <= blockEnd; i++){
        for (unsigned t = 0; t < blockSize; t++) {
            unsigned i = snpIndexVec[t];
            SnpInfo *snp = blockInfo->snpInfoVec[i - blockStart];
            if (snp->skip) {
                valuesPtr[i] = 0.0;
                continue;
            }
            if (badSnps[i]) {
                valuesPtr[i] = 0.0;
                continue;
            }
            if (snp->fitAsFixedEffect) { // flat prior: no spike-and-slab / N(0, snp2pq^S sigma^2) (match ApproxBayesC/R eigen)
                float oldSample = valuesPtr[i];
                int qcolS = eigenQColIndex(blockInfo, i);
                if (qcolS < 0) continue;
                Ref<const VectorXf> Qi = Q.col(qcolS);
                float rhs = (Qi.dot(wcorr) + oldSample) * invVareDn;
                float invLhsFe = 1.0f / invVareDn;
                float uhat = invLhsFe * rhs;
                valuesPtr[i] = uhat + nrnd[i] * sqrtf(invLhsFe);
                wcorr += Qi * (oldSample - valuesPtr[i]);
                what += Qi * valuesPtr[i];
                //ssq[blk] += valuesPtr[i] * valuesPtr[i] / snp2pqPowS[i];
                //++nnz[blk];
                pipPtr[i] = 1.0f;
                continue;
            }
            float oldSample = valuesPtr[i];
            int qcolS2 = eigenQColIndex(blockInfo, i);
            if (qcolS2 < 0) continue;
            Ref<const VectorXf> Qi = Q.col(qcolS2);
            float rhs = (Qi.dot(wcorr) + oldSample)*invVareDn; // times 2pq because the diagonal of ZPZ is not 1 but the variance of genptypes in the case of unstandardised genotypes
            float invLhs = 1.0/(invVareDn + invSigmaSq/snp2pqPowS[i]); // times 2pq because the diagonal of ZPZ is not 1 but the variance of genptypes in the case of unstandardised genotypes
            float uhat = invLhs * rhs;
            float logInvLhsMsigma = logf(invLhs) - logf(snp2pqPowS[i]*sigmaSq);
            float logDelta1 = 0.5*(logInvLhsMsigma + uhat*rhs) + logPi;
            float logDelta0 = logPiComp;
            float probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
            pipPtr[i] = probDelta1;

//            if (bernoulli.sample(probDelta1)) {
            if (urnd[i] < probDelta1) {
//                valuesPtr[i] = normal.sample(uhat, invLhs);
                valuesPtr[i] = uhat + nrnd[i]*sqrtf(invLhs);
                wcorr += Qi*(oldSample - valuesPtr[i]);
                what  += Qi* valuesPtr[i];
                ssq[blk] += valuesPtr[i] * valuesPtr[i] /snp2pqPowS[i];
                ++nnz[blk];
            } else {
                if (oldSample) wcorr += Qi * oldSample;
                valuesPtr[i] = 0.0;
            }
        }
    }
    // cout << "Varei 1 max" << varei.maxCoeff() << endl;
    //cout << ssq << " " << nnz << endl;
    
    wtdSumSq = 0.0;
    numNonZeros = 0;
    nnzPerBlk.setZero(nBlocks);
    for (unsigned blk=0; blk<nBlocks; ++blk) {
        wtdSumSq += ssq[blk];
        numNonZeros += nnz[blk];
        nnzPerBlk[blk] = nnz[blk];
        ssqBlocks[blk] = ssq[blk];
    }
    
    values = VectorXf::Map(valuesPtr, size);
    pip = VectorXf::Map(pipPtr, size);
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


void ApproxBayesS::sampleUnknowns(const unsigned iter){
        
    if (lowRankModel) {
        snpEffects.sampleFromFC_eigen(wcorrBlocks, data.Qblocks, whatBlocks, data.keptLdBlockInfoVec, data.nGWASblock, vareBlk.values, sigmaSq.value, pi.value, varg.value, snp2pqPowS, data.snp2pq);
    } else if (sparse) {
        snpEffects.sampleFromFC_sparse(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.chromInfoVec, sigmaSq.value, pi.value, vare.value, snp2pqPowS, data.snp2pq, varg.value);
    } else {
        snpEffects.sampleFromFC_full(rcorr, data.ZPZ, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, sigmaSq.value, pi.value, vare.value, snp2pqPowS, data.snp2pq, varg.value);
    }
    
    snpEffects.computePosteriorMean(iter);
    snpPip.getValues(snpEffects.pip);
    nnzSnp.getValue(snpEffects.numNonZeros);

    if (estimatePi) pi.sampleFromFC(data.numIncdSnps, snpEffects.numNonZeros);
    
    if (estimateEffectMean) {
        mu.sampleFromFC(data.ZPZsp, snpEffects.values, data.snp2pq, vare.value, rcorr);
        Su.sampleFromFC(data.ZPZsp, snpEffects.values, data.snp2pq, vare.value, mu.snp2pqPowSmu, rcorr);
    }
    
    if (lowRankModel) {
        vargBlk.compute(whatBlocks);
        vareBlk.sampleFromFC(wcorrBlocks, vargBlk.values, snpEffects.ssqBlocks, data.nGWASblock, data.numEigenvalBlock);
        //vareBlk.sampleFromFC(wcorrBlocks, snpEffects.values, data.b, data.nGWASblock, data.keptLdBlockInfoVec);
        varg.value = vargBlk.total;
        vare.value = vareBlk.mean;
    }
    else {
        varg.compute(snpEffects.values, data.ZPy, rcorr);
        vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, status);
        //vare.value = data.varPhenotypic;
    }
//    hsq.compute(varg.value, vare.value);
    hsq.value = varg.value / data.varPhenotypic;

    if (snpEffects.numNonZeros) {
        sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
        S.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros, sigmaSq.value, snpEffects.values, data.snp2pq, snp2pqPowS, logSnp2pq, genVarPrior, sigmaSq.scale, snpEffects.sum2pqSplusOne, scaledGeno);
    }
    if (iter < 1000) {
        genVarPrior += (varg.value - genVarPrior)/iter;
        scalePrior  += (sigmaSq.scale - scalePrior)/iter;
        sigmaSq.scale = scalePrior;
    }
    scale.getValue(sigmaSq.scale);

    if (!(iter % 10)) {
        if (lowRankModel) {
            nBadSnps.compute_eigen(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, wcorrBlocks, data.Qblocks, data.keptLdBlockInfoVec, iter, data);
        } else if (sparse) {
            nBadSnps.compute_sparse(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, rcorr, data.ZPZsp, data.chromInfoVec, iter, data);
        } else {
            nBadSnps.compute_full(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, rcorr, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, iter, data);
        }
    }

    if (!(iter % 100)) {
        if (lowRankModel) {
            rounding.computeWcorr_eigen(data.wcorrBlocks, data.Qblocks, data.keptLdBlockInfoVec, snpEffects.values, wcorrBlocks);
        } else if (sparse) {
            rounding.computeRcorr_sparse(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
        } else {
            rounding.computeRcorr_full(data.ZPy, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
        }
    }
}



void ApproxBayesST::SnpEffects::sampleFromFC(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                             const vector<ChromInfo *> &chromInfoVec, const ArrayXf &hSlT, const VectorXf &snp2pq,
                                             const float sigmaSq, const float pi, const float vare, const float varg) {
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
    
    pip.setZero(size);
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads
    float *pipPtr = pip.data();

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
#pragma omp parallel for schedule(dynamic)
    for (unsigned chr=0; chr<numChr; ++chr) {
        //cout << " thread " << omp_get_thread_num() << " chr " << chr << endl;
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned chrSize = chrEnd - chrStart + 1;

        float oldSample;
        float rhs, invLhs, uhat;
        float logDelta0, logDelta1, probDelta1;
        float logPi = log(pi);
        float logPiComp = log(1.0-pi);
        float invSigmaSq = 1.0f/sigmaSq;
        float varei = varg + vare;
        float sampleDiff;
        
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            oldSample = valuesPtr[i];
            rhs  = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            invLhs = 1.0f/(ZPZdiag[i]/varei + invSigmaSq/hSlT[i]);
            uhat = invLhs*rhs;
            
            logDelta1 = 0.5*(logf(invLhs) - logf(hSlT[i]*sigmaSq) + uhat*rhs) + logPi;
            logDelta0 = logPiComp;
            
            probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
            
            if (urnd[i] < probDelta1) {
                valuesPtr[i] = uhat + nrnd[i]*sqrtf(invLhs);
                sampleDiff = oldSample - valuesPtr[i];
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
    pip = VectorXf::Map(pipPtr, size);
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


void ApproxBayesST::sampleUnknowns(const unsigned iter){
    unsigned cnt=0;
//    do {
        snpEffects.sampleFromFC(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.chromInfoVec, hSlT, data.snp2pq, sigmaSq.value, pi.value, vare.value,
                                varg.value);
//        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
//    } while (snpEffects.numNonZeros == 0);
    sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
    if (estimatePi) pi.sampleFromFC(data.numIncdSnps, snpEffects.numNonZeros);
    nnzSnp.getValue(snpEffects.numNonZeros);
    varg.compute(snpEffects.values, data.ZPy, rcorr);
    vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, status);
    hsq.compute(varg.value, vare.value);
    if (estimateS)
        S.sampleFromFC(snpEffects.numNonZeros, sigmaSq.value, snpEffects.values, data.snp2pq, logSnp2pq,
                       data.LDscore, logLdsc, varg.value, sigmaSq.scale, T.value, hSlT);
    else
        T.sampleFromFC(snpEffects.numNonZeros, sigmaSq.value, snpEffects.values, data.snp2pq, data.LDscore,
                       logLdsc, varg.value, sigmaSq.scale, hSlT);
    scale.getValue(sigmaSq.scale);
    if (!(iter % 100)) rounding.computeRcorr_sparse(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
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

void ApproxBayesR::VgMixComps::compute(const VectorXf &snpEffects, const VectorXf &ZPy, const VectorXf &rcorr, const vector<vector<unsigned> > &snpset, const float varg, const float nobs) {
    values.setZero(size);
    if (varg) {
        for (unsigned k=1; k<size; ++k) {
            unsigned snpSetSize = snpset[k].size();
            for (unsigned j=0; j<snpSetSize; ++j) {
                unsigned snpIdx = snpset[k][j];
                float varj = snpEffects[snpIdx] * (ZPy[snpIdx] - rcorr[snpIdx]);
                values[k] += varj;
            }
            values[k] /= varg * nobs;
            (*this)[k]->value = values[k];
        }
    }
}

void ApproxBayesR::VgMixComps::compute(const VectorXf &snpEffects, const vector<vector<unsigned> > &snpset) {
    values.setZero(size);
    float totalVar = snpEffects.squaredNorm();
    if (totalVar) {
        for (unsigned k=1; k<size; ++k) {
            unsigned snpSetSize = snpset[k].size();
            for (unsigned j=0; j<snpSetSize; ++j) {
                unsigned snpIdx = snpset[k][j];
                float varj = snpEffects[snpIdx] * snpEffects[snpIdx];
                values[k] += varj;
            }
            values[k] /= totalVar;
            (*this)[k]->value = values[k];
        }
    }
}

void ApproxBayesR::VgMixComps::compute(const VectorXf &snpEffects, const vector<unsigned> &membership, const vector<MatrixXf> &Qblocks, const vector<LDBlockInfo*> &keptLdBlockInfoVec) {
    values.setZero(size);
    unsigned nBlocks = Qblocks.size();
    vector<vector<float> > vgBlocks;
    vgBlocks.resize(nBlocks);
    for (unsigned i=0; i<nBlocks; ++i) {
        vgBlocks[i].resize(size-1);
    }
    
#pragma omp parallel for schedule(dynamic)
    for(unsigned blk = 0; blk < nBlocks; blk++){
        Ref<const MatrixXf> Q = Qblocks[blk];
        vector<VectorXf> whatBlock(size-1);
        for (unsigned k=0; k<size-1; ++k) {
            whatBlock[k].resize(Q.rows());
            whatBlock[k].setZero();
        }
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[blk];
        
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd   = blockInfo->endSnpIdx;
        
        for (unsigned j = blockStart; j <= blockEnd; j++){
            int qcol = eigenQColIndex(blockInfo, j);
            if (qcol < 0) continue;
            Ref<const VectorXf> Qj = Q.col(qcol);
            unsigned delta = membership[j];
            if (delta) whatBlock[delta-1] += Qj*snpEffects[j];
        }
        
        for (unsigned k=0; k<size-1; ++k) {
            vgBlocks[blk][k] = whatBlock[k].squaredNorm();
        }
    }
    
    for (unsigned k=1; k<size; ++k) {
        for (unsigned i=0; i<nBlocks; ++i) {
            values[k] += vgBlocks[i][k-1];
        }
        (*this)[k]->value = values[k];
    }
}

void ApproxBayesR::VarEffects::computeRobustMode(const float varg, const VectorXf &snp2pq, const VectorXf &gamma, const VectorXf &pi, const bool noscale){
    if (noscale) {
        value = varg/(snp2pq.array().sum()*gamma.dot(pi));
    } else {
        unsigned n = 0;
        for (unsigned i = 0; i < (unsigned)snp2pq.size(); ++i)
            if (snp2pq[i] > 0.f) ++n;
        if (!n) n = 1;
        value = varg/(n*gamma.dot(pi));  // LDpred2's parameterisation (masked 2pq: count random SNPs)
    }
}

// ==============================================================
// Sparse vector version
// ==============================================================

void ApproxBayesR::SnpEffects::sampleFromFC_sparse(VectorXf &rcorr, const vector<SparseVector<float>> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                            const vector<ChromInfo*> &chromInfoVec,
                                            const VectorXf &snp2pq,
                                            const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare, VectorXf &snpStore,
                                            const float varg,
                                            const bool hsqPercModel, DeltaPi &deltaPi,
                                            const vector<SnpInfo*> &incdSnpInfoVec){
    // -----------------------------------------
    // Initialise the parameters in MCMC sampler
    // -----------------------------------------
    long numChr = chromInfoVec.size();

    float wtdssq[numChr], ssq[numChr], s2pq[numChr], nnz[numChr];
    memset(wtdssq,0,sizeof(float)*numChr);
    memset(ssq,0,sizeof(float)*numChr);
    memset(s2pq,0,sizeof(float)*numChr);
    memset(nnz,0, sizeof(float)*numChr);

    pip.setZero(size);
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads
    float *pipPtr = pip.data();

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
    if (hsqPercModel && varg)
        gp = gamma * 0.01 * varg;
    else
        gp = gamma * sigmaSq;
    
    vector<vector<vector<unsigned> > > snpsetChr(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        snpsetChr[i].resize(ndist);
        for (unsigned k=0; k<ndist; ++k) {
            snpsetChr[i][k].resize(0);
        }
    }

    for (unsigned k=0; k<ndist; ++k) {
        deltaPi[k]->values.setZero(size);
    }

    VectorXf invGamma = gamma.array().inverse();
    invGamma[0] = 0.0;
    
    lambdaVec.setZero(size);
    uhatVec.setZero(size);
    invGammaVec.setZero(size);
    deltaNZ.setZero(size);
    
    // --------------------------------------------------------------------------------
    // Cycle over all variants in the window and sample the genetics effects
    // --------------------------------------------------------------------------------

#pragma omp parallel for schedule(dynamic)
    for (unsigned chr=0; chr<numChr; ++chr) 
    {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned chrSize = chrEnd - chrStart + 1;

        // R specific parameters
        int indistflag;
        double rhs, v1,  b_ls, ssculm, r;
        VectorXf ll, pll, snpindist, var_b_ls;
        ll.setZero(pis.size());
        pll.setZero(pis.size());

        float oldSample;
        float sampleDiff;
        float varei = varg + vare;
        
        // shuffling the SNP index for faster convergence
        vector<int> snpIndexVec = Gadget::shuffle_index(chrStart, chrEnd);

        //for (unsigned i=chrStart; i<=chrEnd; ++i) {
        unsigned i;
        for (unsigned t = 0; t < chrSize; t++) {
            i = snpIndexVec[t];
            if (badSnps[i]) {
                valuesPtr[i] = 0.0;
                continue;
            }
            if (incdSnpInfoVec[i]->fitAsFixedEffect) {
                continue;
            }
            oldSample = valuesPtr[i];
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
              deltaPi[k]->values[i] = pll[k];
           }
            pipPtr[i] = 1.0f - pll[0];

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
                    break;
                }
            }
            snpsetChr[chr][indistflag-1].push_back(i);
            // --------------------------------------------------------------
            // Sample the effect given the group and adjust the rhs
            // --------------------------------------------------------------
            if (indistflag != 1)
            {
                v1 = ZPZdiag[i] + varei / gp((indistflag - 1));
//                valuesPtr[i] = normal.sample(rhs / v1, varei / v1);
                valuesPtr[i] = rhs / v1 + nrnd[i]*sqrtf(varei / v1);
                sampleDiff = oldSample - valuesPtr[i];
                for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {
                    rcorr[it.index()] += it.value() * sampleDiff;
                }
                ssq[chr]  += valuesPtr[i]*valuesPtr[i];
                wtdssq[chr]  += (valuesPtr[i]*valuesPtr[i]) / gamma[indistflag - 1];
                s2pq[chr] += snp2pq[i];
                deltaNZ[i] = 1;
                ++nnz[chr];
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
    wtdSumSq = 0.0;
    sum2pq = 0.0;
    numNonZeros = 0;                                                                                                                                                                
    nnzPerChr.setZero(numChr);                                                                                                                                                      
    snpStore.setZero(ndist);
    snpset.resize(ndist);
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
    }
    for (unsigned i=0; i<numChr; ++i) {
        sumSq += ssq[i];
        wtdSumSq += wtdssq[i];
        sum2pq += s2pq[i];
        numNonZeros += nnz[i];
        nnzPerChr[i] = nnz[i];                                                                                                                                                      
        for (unsigned k=0; k<ndist; ++k) {
            for (unsigned j=0; j<snpsetChr[i][k].size(); ++j) {
                snpset[k].push_back(snpsetChr[i][k][j]);
                snpStore[k]++;
            }
        }
    }

    values = VectorXf::Map(valuesPtr, size);
    pip = VectorXf::Map(pipPtr, size);
}

// ==============================================================
// Vector of vectors version
// ==============================================================

void ApproxBayesR::SnpEffects::sampleFromFC_full(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const VectorXf &snp2pq,
                                            const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare, VectorXf &snpStore,
                                            const float varg,
                                            const bool hsqPercModel, DeltaPi &deltaPi,
                                            const vector<SnpInfo*> *incdSnpInfoVec){
    // -----------------------------------------
    // Initialise the parameters in MCMC sampler
    // -----------------------------------------
    long numChr = chromInfoVec.size();
    
    float wtdssq[numChr], ssq[numChr], nnz[numChr], s2pq[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(wtdssq,0,sizeof(float)*numChr);
    memset(nnz,0,sizeof(float)*numChr);
    memset(s2pq,0,sizeof(float)*numChr);

    pip.setZero(size);
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads
    float *pipPtr = pip.data();

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
    if (hsqPercModel && varg)
        gp = gamma * 0.01 * varg;
    else
        gp = gamma * sigmaSq;

    vector<vector<vector<unsigned> > > snpsetChr(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        snpsetChr[i].resize(ndist);
        for (unsigned k=0; k<ndist; ++k) {
            snpsetChr[i][k].resize(0);
        }
    }

    for (unsigned k=0; k<ndist; ++k) {
        deltaPi[k]->values.setZero(size);
    }

    // --------------------------------------------------------------------------------
    // Cycle over all variants in the window and sample the genetics effects
    // --------------------------------------------------------------------------------

#pragma omp parallel for schedule(dynamic)
    for (unsigned chr=0; chr<numChr; ++chr) 
    {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned chrSize = chrEnd - chrStart + 1;
        
        float oldSample;
        float varei = varg + vare;
        double rhs, invLhs, uhat;
        
        int indistflag;
        double v1,  b_ls, ssculm, r;
        VectorXf ll, pll, snpindist, var_b_ls;
        ll.setZero(pis.size());
        pll.setZero(pis.size());

        // shuffling the SNP index for faster convergence
        vector<int> snpIndexVec = Gadget::shuffle_index(chrStart, chrEnd);

        unsigned i;
        for (unsigned t = 0; t < chrSize; t++) {
            i = snpIndexVec[t];
            if (badSnps[i]) {
                valuesPtr[i] = 0.0;
                continue;
            }
            oldSample = valuesPtr[i];
            if (incdSnpInfoVec && (*incdSnpInfoVec)[i]->fitAsFixedEffect) {
                continue;
            }
            // ---------------------------------------------
            // Calculate residual variance including a 
            // correction for the sampling variation and
            // LD ignored
            // ---------------------------------------------
            //varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
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
              deltaPi[k]->values[i] = pll[k];
           }
            pipPtr[i] = 1.0f - pll[0];
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
                    break;
                }
            }
            snpsetChr[chr][indistflag-1].push_back(i);
            // --------------------------------------------------------------
            // Sample the effect given the group and adjust the rhs
            // --------------------------------------------------------------
            if (indistflag != 1)
            {
                v1 = ZPZdiag[i] + varei / gp((indistflag - 1));
//                valuesPtr[i] = normal.sample(rhs / v1, varei / v1);
                valuesPtr[i] = rhs / v1 + nrnd[i]*sqrtf(varei / v1);
                rcorr.segment(windStart[i], windSize[i]) += ZPZ[i] * (oldSample - valuesPtr[i]);
                ssq[chr] += valuesPtr[i] * valuesPtr[i];
                wtdssq[chr] += (valuesPtr[i] * valuesPtr[i]) / gamma[indistflag - 1];
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
    wtdSumSq = 0.0;
    sum2pq = 0.0;
    numNonZeros = 0.0;
    
    nnzPerChr.setZero(numChr);
    snpStore.setZero(ndist);
    snpset.resize(ndist);
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
    }
    for (unsigned i=0; i<numChr; ++i) {
        sumSq += ssq[i]; 
        wtdSumSq += wtdssq[i];
        sum2pq += s2pq[i];
        numNonZeros += nnz[i];
        nnzPerChr[i] = nnz[i];
        for (unsigned k=0; k<ndist; ++k) {
            for (unsigned j=0; j<snpsetChr[i][k].size(); ++j) {
                snpset[k].push_back(snpsetChr[i][k][j]);
                snpStore[k]++;
            }
        }
    }
    
    values = VectorXf::Map(valuesPtr, size); 
    pip = VectorXf::Map(pipPtr, size);
}

void ApproxBayesR::SnpEffects::sampleFromFC_eigen(vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
                                            const vector<LDBlockInfo*> &keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks,
                                            const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, VectorXf &snpStore, const float varg,
                                            const bool hsqPercModel, DeltaPi &deltaPi) {
    // -----------------------------------------
    // This method uses low-rank model with eigen-decomposition of LD matrices
    // -----------------------------------------
    long nBlocks = keptLdBlockInfoVec.size();
    
    whatBlocks.resize(nBlocks);
    for (unsigned i=0; i<nBlocks; ++i) {
        whatBlocks[i].resize(wcorrBlocks[i].size());
    }

    float ssq[nBlocks], wtdssq[nBlocks], s2pq[nBlocks], nnz[nBlocks];
    memset(ssq,0, sizeof(float)*nBlocks);
    memset(wtdssq,0, sizeof(float)*nBlocks);
    memset(s2pq,0,sizeof(float)*nBlocks);
    memset(nnz,0, sizeof(float)*nBlocks);
    
    pip.setZero(size);
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads
    float *pipPtr = pip.data();

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
    membership.setZero(size);
    
    // R specific parameters
    int ndist = pis.size();
    ArrayXf logPis = pis.array().log();
    ArrayXf wtdSigmaSq(ndist);
    ArrayXf invWtdSigmaSq(ndist);
    ArrayXf logWtdSigmaSq(ndist);

    if (hsqPercModel && varg) {
        wtdSigmaSq = gamma * 0.01 * varg;
    } else {
        wtdSigmaSq = gamma * sigmaSq;
    }
    
    invWtdSigmaSq = wtdSigmaSq.inverse();
    logWtdSigmaSq = wtdSigmaSq.log();

    vector<vector<vector<unsigned> > > snpsetBlocks(nBlocks);
    for (unsigned i=0; i<nBlocks; ++i) {
        snpsetBlocks[i].resize(ndist);
        for (unsigned k=0; k<ndist; ++k) {
            snpsetBlocks[i][k].resize(0);
        }
    }

    for (unsigned k=0; k<ndist; ++k) {
        deltaPi[k]->values.setZero(size);
    }

    // --------------------------------------------------------------------------------
    // Cycle over all variants in the window and sample the genetics effects
    // --------------------------------------------------------------------------------
    
    #pragma omp parallel for schedule(dynamic)
    for(unsigned blk = 0; blk < nBlocks; blk++){
        //printf("  Inner: Thread %d of %d\n", omp_get_thread_num(), omp_get_num_threads());
        Ref<const MatrixXf> Q = Qblocks[blk];
        Ref<VectorXf> wcorr = wcorrBlocks[blk];
        Ref<VectorXf> what = whatBlocks[blk];

        what.setZero();
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[blk];
        
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd   = blockInfo->endSnpIdx;
        unsigned blockSize  = blockEnd - blockStart + 1;

        float invVareDn = nGWASblocks[blk] / vareBlocks[blk];

        ArrayXf invLhs = 1.0/(invVareDn + invWtdSigmaSq);
        ArrayXf logInvLhsMsigma = invLhs.log() - logWtdSigmaSq;
        
        // shuffling the SNP index for faster convergence
        vector<int> snpIndexVec = Gadget::shuffle_index(blockStart, blockEnd);

        //for(unsigned i = blockStart; i <= blockEnd; i++){
        for (unsigned t = 0; t < blockSize; t++) {
            unsigned i = snpIndexVec[t];
            SnpInfo *snp = blockInfo->snpInfoVec[i - blockStart];
            if (snp->skip) {
                valuesPtr[i] = 0.0;
                // Skipped SNPs are not sampled; leave mixture probs at null-only so output PIP (1 - DeltaPi0) is not spuriously 1 (see gctb.cpp output for R).
                pipPtr[i] = 0.0f;
                for (unsigned k = 0; k < (unsigned)ndist; ++k)
                    deltaPi[k]->values[i] = (k == 0) ? 1.0f : 0.0f;
                continue;
            }
            if (badSnps[i]) {
                valuesPtr[i] = 0.0;
                continue;
            }
            if (snp->fitAsFixedEffect) {
                continue;
            }
            float oldSample = valuesPtr[i];
            int qcolR = eigenQColIndex(blockInfo, i);
            if (qcolR < 0) continue;
            Ref<const VectorXf> Qi = Q.col(qcolR);
            float rhs = (Qi.dot(wcorr) + oldSample) * invVareDn;
            ArrayXf uhat = invLhs * rhs;
            ArrayXf logDelta = 0.5 * (logInvLhsMsigma + uhat * rhs) + logPis;
            logDelta[0] = logPis[0];

            ArrayXf probDelta(ndist);
            for (unsigned k = 0; k < ndist; ++k) {
                probDelta[k] = 1.0f / (logDelta - logDelta[k]).exp().sum();
                deltaPi[k]->values[i] = probDelta[k];
            }
            pipPtr[i] = 1.0f - probDelta[0];

            unsigned delta = bernoulli.sample(probDelta, urnd[i]);
            membership[i] = delta;
            snpsetBlocks[blk][delta].push_back(i);

            if (delta) {
                valuesPtr[i] = uhat[delta] + nrnd[i] * sqrtf(invLhs[delta]);
                wcorr += Qi * (oldSample - valuesPtr[i]);
                what += Qi * valuesPtr[i];
                ssq[blk] += valuesPtr[i] * valuesPtr[i];
                wtdssq[blk] += (valuesPtr[i] * valuesPtr[i]) / gamma[delta];
                ++nnz[blk];
            } else {
                if (oldSample) wcorr += Qi * oldSample;
                valuesPtr[i] = 0.0;
            }
        }

    }

    // ---------------------------------------------------------------------
    // Tally up the effect sum of squares and the number of non-zero effects
    // ---------------------------------------------------------------------
    sumSq = 0.0;
    wtdSumSq = 0.0;
    numNonZeros = 0;
    nnzPerBlk.setZero(nBlocks);
    ssqBlocks.setZero(nBlocks);
    snpStore.setZero(ndist);
    snpset.resize(ndist);
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
    }
    for (unsigned blk=0; blk<nBlocks; ++blk) {
        sumSq += ssq[blk];
        wtdSumSq += wtdssq[blk];
        numNonZeros += nnz[blk];
        nnzPerBlk[blk] = nnz[blk];
        ssqBlocks[blk] = ssq[blk];
        for (unsigned k=0; k<ndist; ++k) {
            for (unsigned j=0; j<snpsetBlocks[blk][k].size(); ++j) {
                snpset[k].push_back(snpsetBlocks[blk][k][j]);
                snpStore[k]++;
            }
        }
    }
    values = VectorXf::Map(valuesPtr, size);
    pip = VectorXf::Map(pipPtr, size);
}

void ApproxBayesR::SnpEffects::sampleFromFC(const VectorXf &ZPy, const VectorXf &ZPZdiag, const MatrixXf &Z, const float n_ref, const float n_gwas,
                                            const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare,
                                            VectorXf &snpStore, VectorXf &ghat, const float varg, const bool hsqPercModel, DeltaPi &deltaPi) {
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
        if (hsqPercModel && varg)
            gp = gamma * 0.01 * varg;
        else
            gp = gamma * sigmaSq;
    //    cout << varg << " " << gp.transpose() << endl;
        snpset.resize(ndist);
        for (unsigned k=0; k<ndist; ++k) {
            snpset[k].resize(0);
        }
        
    for (unsigned k=0; k<ndist; ++k) {
        deltaPi[k]->values.setZero(size);
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
                deltaPi[k]->values[i] = pll[k];
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
                sumSq += values[i] * values[i];
                wtdSumSq += (values[i] * values[i]) / gamma[indistflag - 1];
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
                                            const float varg, const float vare, const float ps, const float overdispersion, const bool hsqPercModel, DeltaPi &deltaPi) {
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
    
    if (hsqPercModel && varg) {
        wtdSigmaSq = gamma * 0.01 * varg;
    } else {
        wtdSigmaSq = gamma * sigmaSq;
    }
    
    invWtdSigmaSq = wtdSigmaSq.inverse();
    logWtdSigmaSq = wtdSigmaSq.log();
    
    VectorXf invGamma = gamma.inverse();
    invGamma[0] = 0;

    for (unsigned k=0; k<ndist; ++k) {
        deltaPi[k]->values.setZero(size);
    }


    for (unsigned i=0; i<size; ++i) {
        
        varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
        
        rhs  = rcorr[i] + ZPZdiag[i]*values[i];
        
        invLhs = (ZPZdiag[i] + varei*invWtdSigmaSq).inverse();
        uhat = invLhs*rhs;
        
        logDelta = 0.5*(invLhs.log() - logWtdSigmaSq + uhat*rhs) + logPis;
        logDelta[0] = logPis[0];
        
        for (unsigned k=0; k<ndist; ++k) {
            probDelta[k] = 1.0f/(logDelta-logDelta[k]).exp().sum();
            deltaPi[k]->values[i] = probDelta[k];
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
    
    sumSq = values.squaredNorm();
    wtdSumSq = values.cwiseProduct(invGamma).dot(values);
    
    rcorr = ZPy - ZPZsp * values;
}

void ApproxBayesR::SnpEffects::sampleFromPrior(const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float varg, const bool hsqPercModel){
    int ndist = pis.size();
    ArrayXf wtdSigmaSq(ndist);

    if (hsqPercModel && varg) {
        wtdSigmaSq = gamma * 0.01 * varg;
    } else {
        wtdSigmaSq = gamma * sigmaSq;
    }
    
    for (unsigned i=0; i<size; ++i) {
        unsigned delta = bernoulli.sample(pis);
        if (delta) {
            values[i] = Stat::snorm()*sqrtf(wtdSigmaSq[delta]);
        }
        else {
            values[i] = 0.0;
        }
    }
}

void ApproxBayesR::updateRHSfull(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec, const VectorXf &snpEffects){
#pragma omp parallel for
    for (unsigned chr=0; chr<chromInfoVec.size(); ++chr) {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            rcorr.segment(windStart[i], windSize[i]) -= ZPZ[i]*snpEffects[i];
        }
    }
}

void ApproxBayesR::updateRHSsparse(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZ, const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec, const VectorXf &snpEffects){
#pragma omp parallel for
    for (unsigned chr=0; chr<chromInfoVec.size(); ++chr) {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            for (SparseVector<float>::InnerIterator it(ZPZ[i]); it; ++it) {
                rcorr[it.index()] -= it.value() * snpEffects[i];
            }
        }
    }
}

void ApproxBayesR::updateRHSlowRankModel(vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, const vector<LDBlockInfo*> &keptLdBlockInfoVec, const VectorXf &snpEffects){
    long nBlocks = keptLdBlockInfoVec.size();
    
#pragma omp parallel for schedule(dynamic)
    for(unsigned blk = 0; blk < nBlocks; blk++){
        Ref<const MatrixXf> Q = Qblocks[blk];
        Ref<VectorXf> wcorr = wcorrBlocks[blk];
                
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[blk];
        
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd   = blockInfo->endSnpIdx;
                
        for(unsigned i = blockStart; i <= blockEnd; i++){
            int qcol = eigenQColIndex(blockInfo, i);
            if (qcol < 0) continue;
            Ref<const VectorXf> Qi = Q.col(qcol);
            if (snpEffects[i]) {
                wcorr -= Qi*snpEffects[i];
            }
        }
    }

}

void ApproxBayesR::sampleUnknowns(const unsigned iter){
    if (lowRankModel) {
        fixedSnpEffects.sampleFromFC_eigen(snpEffects, wcorrBlocks, data.Qblocks, whatBlocks, data.keptLdBlockInfoVec, data.nGWASblock, vareBlk.values, data.numFixedEffectSnps);
        snpEffects.sampleFromFC_eigen(wcorrBlocks, data.Qblocks, whatBlocks, data.keptLdBlockInfoVec, data.nGWASblock, vareBlk.values, sigmaSq.value, Pis.values, gamma.values, snpStore, vargRandom.value, hsqPercModel, deltaPi);
    } else if (sparse) {
        fixedSnpEffects.sampleFromFC_sparse(snpEffects, rcorr, data.ZPZsp, data.ZPZdiag, data.snp2pq, data.incdSnpInfoVec, varg.value, vare.value, data.numFixedEffectSnps);
        snpEffects.sampleFromFC_sparse(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.chromInfoVec, data.snp2pq, sigmaSq.value, Pis.values, gamma.values, vare.value, snpStore, vargRandom.value, hsqPercModel, deltaPi, data.incdSnpInfoVec);
    } else {
        fixedSnpEffects.sampleFromFC_full(snpEffects, rcorr, data.ZPZ, data.ZPZdiag, data.windStart, data.windSize, data.snp2pq, data.incdSnpInfoVec, varg.value, vare.value, data.numFixedEffectSnps);
        snpEffects.sampleFromFC_full(rcorr, data.ZPZ, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, data.snp2pq, sigmaSq.value, Pis.values, gamma.values, vare.value, snpStore, vargRandom.value, hsqPercModel, deltaPi, &data.incdSnpInfoVec);
    }
    
//    if (!(iter % 10)) {   // To improve mixing, apply tempered Gibbs sampling on high-LD SNPs in every 10 iterations
//        snpEffects.sampleFromTGS_eigen(highLDsnpSet, wcorrBlocks, data.Qblocks, whatBlocks,
//                                      data.keptLdBlockInfoVec, data.nGWASblock, vareBlk.values,
//                                      Pis.values, gamma.values, vargRandom.value, hsqPercModel, sigmaSq.value);
//    }
    
    snpEffects.computePosteriorMean(iter);
    snpPip.getValues(snpEffects.pip);
    nnzSnp.getValue(snpEffects.numNonZeros);
    numSnps.getValues(snpStore);
    
    if (robustMode) {
        sigmaSq.computeRobustMode(varg.value, data.numFixedEffectSnps ? data.snp2pqForRandomSnpEffects : data.snp2pq, Pis.values, gamma.values, noscale);
    } else {
        sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
    }
    
    if (estimatePi) Pis.sampleFromFC(snpStore);
    
    if (lowRankModel) {
        vargBlk.compute(whatBlocks);
        vareBlk.sampleFromFC(wcorrBlocks, vargBlk.values, snpEffects.ssqBlocks, data.nGWASblock, data.numEigenvalBlock);
        //vareBlk.sampleFromFC(wcorrBlocks, snpEffects.values, data.b, data.nGWASblock, data.keptLdBlockInfoVec);
        vargRandom.value = vargBlk.total;
        vare.value = vareBlk.mean;
    } else {
        varg.compute(snpEffects.values, data.ZPy, rcorr);
        vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, status);
        if (data.numFixedEffectSnps && data.varGenotypic > 0.f)
            vargRandom.value = varg.value * (data.vargRandomSnpEffects / data.varGenotypic);
        else
            vargRandom.value = varg.value;
    }
    hsqRandom.compute(vargRandom.value, data.varPhenotypic);
    if (data.numFixedEffectSnps) {
        hsqFixed.compute(snpEffects.values, data);
        hsq.value = hsqRandom.value + hsqFixed.value;
        if (lowRankModel)
            varg.value = vargRandom.value + hsqFixed.value * data.varPhenotypic;
    } else {
        hsq.value = hsqRandom.value;
        if (lowRankModel)
            varg.value = vargRandom.value;
    }
    
    Vgs.compute(snpEffects.values, snpEffects.snpset);
    
    snpHsqPep.compute(snpEffects.values, vargRandom.value);

    if (!(iter % 10)) {
        if (lowRankModel) {
            nBadSnps.compute_eigen(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, wcorrBlocks, data.Qblocks, data.keptLdBlockInfoVec, iter, data);
        } else if (sparse) {
            nBadSnps.compute_sparse(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, rcorr, data.ZPZsp, data.chromInfoVec, iter, data);
        } else {
            nBadSnps.compute_full(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, rcorr, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, iter, data);
        }
    }

    if (!(iter % 100)) {
        if (lowRankModel) {
            rounding.computeWcorr_eigen(data.wcorrBlocks, data.Qblocks, data.keptLdBlockInfoVec, snpEffects.values, wcorrBlocks);
        } else if (sparse) {
            rounding.computeRcorr_sparse(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
        } else {
            rounding.computeRcorr_full(data.ZPy, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
        }
    }
}


// *******************************************************
// Approximate Bayes RS
// *******************************************************

void ApproxBayesRS::SnpEffects::sampleFromFC_sparse(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                             const vector<ChromInfo *> &chromInfoVec,
                                             const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare,
                                             const ArrayXf &snp2pqPowS, const VectorXf &snp2pq,
                                             const float varg,
                                             const bool hsqPercModel) {
    // sample SNP effects with a sparse LD matrix
    
    long numChr = chromInfoVec.size();
    
    float ssq[numChr], s2pq[numChr], nnz[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(s2pq,0,sizeof(float)*numChr);
    memset(nnz,0, sizeof(float)*numChr);
    
    pip.setZero(size);
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads
    float *pipPtr = pip.data();

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
    
    if (hsqPercModel && varg) {
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
        unsigned chrSize = chrEnd - chrStart + 1;

        float oldSample;
        float sampleDiff;
        float rhs;
        float varei = varg + vare;
        
        ArrayXf invLhs(ndist);
        ArrayXf uhat(ndist);
        ArrayXf logDelta(ndist);
        ArrayXf probDelta(ndist);
        
        unsigned delta;
        
        // shuffling the SNP index for faster convergence
        vector<int> snpIndexVec = Gadget::shuffle_index(chrStart, chrEnd);

        unsigned i;
        for (unsigned t = 0; t < chrSize; t++) {
            i = snpIndexVec[t];
            if (badSnps[i]) {
                valuesPtr[i] = 0.0;
                continue;
            }
            oldSample = valuesPtr[i];
            
            rhs  = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            
            invLhs = (ZPZdiag[i]/varei + invWtdSigmaSq/snp2pqPowS[i]).inverse();
            uhat = invLhs*rhs;
            
            logDelta = 0.5*(invLhs.log() - log2pqPowS[i] - logWtdSigmaSq + uhat*rhs) + logPis;
            logDelta[0] = logPis[0];
            
            for (unsigned k=0; k<ndist; ++k) {
                probDelta[k] = 1.0f/(logDelta-logDelta[k]).exp().sum();
            }
            pipPtr[i] = 1.0f - probDelta[0];

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
    pip = VectorXf::Map(pipPtr, size);
}

void ApproxBayesRS::SnpEffects::sampleFromFC_full(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                             const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo *> &chromInfoVec,
                                             const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare,
                                             const ArrayXf &snp2pqPowS, const VectorXf &snp2pq,
                                             const float varg,
                                             const bool hsqPercModel) {
    // sample SNP effects with a full LD matrix
    
    long numChr = chromInfoVec.size();
    
    float ssq[numChr], s2pq[numChr], nnz[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(s2pq,0,sizeof(float)*numChr);
    memset(nnz,0, sizeof(float)*numChr);
    
    pip.setZero(size);
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads
    float *pipPtr = pip.data();

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
    
    if (hsqPercModel && varg) {
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
        unsigned chrSize = chrEnd - chrStart + 1;

        float oldSample;
        float sampleDiff;
        float rhs;
        float varei = varg + vare;
        
        ArrayXf invLhs(ndist);
        ArrayXf uhat(ndist);
        ArrayXf logDelta(ndist);
        ArrayXf probDelta(ndist);
        
        unsigned delta;
        
        // shuffling the SNP index for faster convergence
        vector<int> snpIndexVec = Gadget::shuffle_index(chrStart, chrEnd);

        unsigned i;
        for (unsigned t = 0; t < chrSize; t++) {
            i = snpIndexVec[t];
            if (badSnps[i]) {
                valuesPtr[i] = 0.0;
                continue;
            }
            oldSample = valuesPtr[i];
            
            rhs  = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            
            invLhs = (ZPZdiag[i]/varei + invWtdSigmaSq/snp2pqPowS[i]).inverse();
            uhat = invLhs*rhs;
            
            logDelta = 0.5*(invLhs.log() - log2pqPowS[i] - logWtdSigmaSq + uhat*rhs) + logPis;
            logDelta[0] = logPis[0];
            
            for (unsigned k=0; k<ndist; ++k) {
                probDelta[k] = 1.0f/(logDelta-logDelta[k]).exp().sum();
            }
            pipPtr[i] = 1.0f - probDelta[0];

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
    pip = VectorXf::Map(pipPtr, size);
}

void ApproxBayesRS::Sp::sampleFromFC(vector<vector<unsigned> > &snpset, const VectorXf &snpEffects,
                                     float &sigmaSq, const VectorXf &gamma,
                                     const VectorXf &snp2pq, ArrayXf &snp2pqPowS, const ArrayXf &logSnp2pq,
                                     const float vg, float &scale, float &sum2pqSplusOne, const bool hsqPercModel) {
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
    float curr_H = computeU(curr, nnzMix, snpEffectMix, snp2pqLogSumNZ, snp2pqMix, logSnp2pqMix, sigmaSq, gamma, vg, scaleCurr, hsqPercModel) + 0.5*curr_p*curr_p;
    float cand_H = computeU(cand, nnzMix, snpEffectMix, snp2pqLogSumNZ, snp2pqMix, logSnp2pqMix, sigmaSq, gamma, vg, scaleCand, hsqPercModel) + 0.5*cand_p*cand_p;
    
    if (Stat::ranf() < exp(curr_H-cand_H)) {  // accept
        value = cand;
        scale = scaleCand;
        snp2pqPowS = snp2pq.array().pow(cand);
        sum2pqSplusOne = 0.0;
        for (unsigned i=0; i<nnzMix; ++i) sum2pqSplusOne += snp2pqMix[i].pow(value+1.0).sum();
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

float ApproxBayesRS::Sp::computeU(const float S, const unsigned nnzMix, const vector<ArrayXf> &snpEffectMix, const float snp2pqLogSum, const vector<ArrayXf> &snp2pqMix, const vector<ArrayXf> &logSnp2pqMix, const float sigmaSq, const VectorXf &gamma, const float vg, float &scale, const bool hsqPercModel) {
    vector<ArrayXf> snp2pqPowSMix(nnzMix);
    float constantA = snp2pqLogSum;
    ArrayXf constantB(nnzMix);
    ArrayXf constantC(nnzMix);
    for (unsigned i=0; i<nnzMix; ++i) {
        snp2pqPowSMix[i] = snp2pqMix[i].pow(S);
        constantB[i] = (snpEffectMix[i].square()/snp2pqPowSMix[i]).sum()/gamma[i+1];
        if (hsqPercModel) constantC[i] = snp2pqPowSMix[i].sum()*gamma[i+1];
        else constantC[i] = (snp2pqMix[i]*snp2pqPowSMix[i]).sum()*gamma[i+1];
    }
    scale = 0.5*vg/constantC.sum();
    return 0.5*S*constantA + 0.5/sigmaSq*constantB.sum() + 0.5*S*S/var;
}

void ApproxBayesRS::sampleUnknowns(const unsigned iter) {
    if (sparse) {
        snpEffects.sampleFromFC_sparse(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.chromInfoVec, sigmaSq.value, Pis.values, gamma.values, vare.value, snp2pqPowS, data.snp2pq, varg.value, hsqPercModel);
    } else {
        snpEffects.sampleFromFC_full(rcorr, data.ZPZ, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, sigmaSq.value, Pis.values, gamma.values, vare.value, snp2pqPowS, data.snp2pq, varg.value, hsqPercModel);
    }
    snpEffects.computePosteriorMean(iter);
    snpPip.getValues(snpEffects.pip);
    nnzSnp.getValue(snpEffects.numNonZeros);
    numSnps.getValues(snpStore);

    sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
    
    if (estimatePi) Pis.sampleFromFC(snpEffects.numSnpMix);
    
    varg.compute(snpEffects.values, data.ZPy, rcorr);
    vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, status);
    hsq.compute(varg.value, vare.value);
    
    Vgs.compute(snpEffects.values, snpEffects.snpset);

    S.sampleFromFC(snpEffects.snpset, snpEffects.values, sigmaSq.value, gamma.values, data.snp2pq, snp2pqPowS, logSnp2pq, genVarPrior, sigmaSq.scale, snpEffects.sum2pqSplusOne, hsqPercModel);
    
    if (iter < 1000) {
        genVarPrior += (varg.value - genVarPrior)/iter;
        scalePrior  += (sigmaSq.scale - scalePrior)/iter;
        sigmaSq.scale = scalePrior;
    }
    scale.getValue(sigmaSq.scale);

    if (!(iter % 10)) {
        if (sparse) {
            nBadSnps.compute_sparse(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, rcorr, data.ZPZsp, data.chromInfoVec, iter, data);
        } else {
            nBadSnps.compute_full(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, rcorr, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, iter, data);
        }
    }

    if (!(iter % 100)) {
        if (sparse)
            rounding.computeRcorr_sparse(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
        else
            rounding.computeRcorr_full(data.ZPy, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    }
}


void ApproxBayesSMix::SnpEffects::sampleFromFC_sparse(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy, const vector<ChromInfo *> &chromInfoVec, const ArrayXf &snp2pqPowS, const VectorXf &snp2pq, const Vector2f &sigmaSq, const Vector3f &pi, const float vare, const float varg, VectorXf &deltaS) {
    long numChr = chromInfoVec.size();
    
    wtdSum2pq.setZero();
    wtdSumSq.setZero();
    numNonZeros.setZero();
    numSnpMixComp.setZero();
    
    valuesMixCompS.setZero(size);
    deltaS.setZero(size);
    pip.setZero(size);
    
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
        unsigned chrSize = chrEnd - chrStart + 1;

        float oldSample;
        float sampleDiff;
        float rhs;
        float varei = varg + vare;
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
            pip[i] = 1.0f - probDelta[0];
            
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

void ApproxBayesSMix::SnpEffects::sampleFromFC_full(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy, const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo *> &chromInfoVec, const ArrayXf &snp2pqPowS, const VectorXf &snp2pq, const Vector2f &sigmaSq, const Vector3f &pi, const float vare, const float varg, VectorXf &deltaS) {
    long numChr = chromInfoVec.size();
    
    wtdSum2pq.setZero();
    wtdSumSq.setZero();
    numNonZeros.setZero();
    numSnpMixComp.setZero();
    
    valuesMixCompS.setZero(size);
    deltaS.setZero(size);
    pip.setZero(size);
    
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
        unsigned chrSize = chrEnd - chrStart + 1;

        float oldSample;
        float sampleDiff;
        float rhs;
        float varei = varg + vare;
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
            pip[i] = 1.0f - probDelta[0];

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

void ApproxBayesSMix::sampleUnknowns(const unsigned iter) {
    unsigned cnt=0;
//    do {
        if (sparse) {
            snpEffects.sampleFromFC_sparse(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.chromInfoVec,
                                    snp2pqPowS, data.snp2pq, sigmaSq.values,
                                    piMixComp.values, vare.value, varg.value, deltaS.values);
        } else {
            snpEffects.sampleFromFC_full(rcorr, data.ZPZ, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec,
                                    snp2pqPowS, data.snp2pq, sigmaSq.values,
                                    piMixComp.values, vare.value, varg.value, deltaS.values);
        }
//        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
//    } while (snpEffects.numNonZeros.sum() == 0);
    piMixComp.sampleFromFC(snpEffects.numSnpMixComp);
    pi.sampleFromFC(data.numIncdSnps, snpEffects.numNonZeros.sum());
    sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
    nnzSnp.getValue(snpEffects.numNonZeros.sum());
    varg.compute(snpEffects.values, data.ZPy, rcorr);
    vargMixComp.compute(sigmaSq.values, snpEffects.wtdSum2pq);
    vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, status);
    hsq.compute(varg.value, vare.value);
    hsqMixComp.compute(vargMixComp.values, varg.value, vare.value);
    
//    if (iter < 1000) {
//        sigmaSq.computeScale(hsqMixComp.values, snpEffects.wtdSum2pq);
//        scalePrior += (sigmaSq[0]->scale - scalePrior)/iter;
//        sigmaSq[0]->scale = scalePrior;
//    }
    
    S.sampleFromFC(snpEffects.wtdSumSq[1], snpEffects.numNonZeros[1], sigmaSq.values[1], snpEffects.valuesMixCompS, data.snp2pq, snp2pqPowS, logSnp2pq, vargMixComp.values[1], sigmaSq[1]->scale, snpEffects.sum2pqSplusOne, true);
    
    if (!(iter % 10)) {
        if (sparse) {
            nBadSnps.compute_sparse(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, rcorr, data.ZPZsp, data.chromInfoVec, iter, data);
        } else {
            nBadSnps.compute_full(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, rcorr, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, iter, data);
        }
    }

    if (!(iter % 100)) {
        if (sparse)
            rounding.computeRcorr_sparse(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
        else
            rounding.computeRcorr_full(data.ZPy, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    }
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

void BayesSMix::sampleUnknowns(const unsigned iter){
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }

    unsigned cnt=0;
//    do {
        snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, snp2pqPowS, data.snp2pq, sigmaSq.values, piMixComp.values, vare.value, deltaS.values, ghat, ghatMixComp);
//        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
//    } while (snpEffects.numNonZeros.sum() == 0);
    
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
    
    // BayesSMix inherits from BayesS, so it has scaledGeno member (defaults to true since noscale=false)
    S.sampleFromFC(snpEffects.wtdSumSq[1], snpEffects.numNonZeros[1], sigmaSq.values[1], snpEffects.valuesMixCompS, data.snp2pq, snp2pqPowS, logSnp2pq, vargMixComp.values[1], sigmaSq[1]->scale, snpEffects.sum2pqSplusOne, scaledGeno);
    
    if (!(iter % 100)) rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);
    
}

void ApproxBayesRC::SnpEffects::sampleFromFC_sparse(VectorXf &rcorr, const vector<SparseVector<float>> &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                            const vector<ChromInfo*> &chromInfoVec,
                                            const float sigmaSq, const MatrixXf &snpPi, const VectorXf &gamma, const float vare, const float varg,
                                            const bool hsqPercModel, DeltaPi &deltaPi, const vector<SnpInfo*> &incdSnpInfoVec){


    long numChr = chromInfoVec.size();

    float ssq[numChr], wtdssq[numChr], s2pq[numChr], nnz[numChr];
    memset(ssq,0,sizeof(float)*numChr);
    memset(wtdssq,0,sizeof(float)*numChr);
    memset(s2pq,0,sizeof(float)*numChr);
    memset(nnz,0, sizeof(float)*numChr);

    pip.setZero(size);
    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads
    float *pipPtr = pip.data();

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

    if (hsqPercModel && varg) {
        wtdSigmaSq = gamma * 0.01 * varg;
    } else {
        wtdSigmaSq = gamma * sigmaSq;
    }

    invWtdSigmaSq = wtdSigmaSq.inverse();
    logWtdSigmaSq = wtdSigmaSq.log();

    vector<vector<vector<unsigned> > > snpsetChr(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        snpsetChr[i].resize(ndist);
        for (unsigned k=0; k<ndist; ++k) {
            snpsetChr[i][k].resize(0);
        }
    }

    for (unsigned k=0; k<ndist; ++k) {
        deltaPi[k]->values.setZero(size);
    }


    #pragma omp parallel for schedule(dynamic)
    for (unsigned chr=0; chr<numChr; ++chr)
    {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned chrSize = chrEnd - chrStart + 1;

        float oldSample;
        float sampleDiff;
        float rhs;
        float varei = varg + vare;

        ArrayXf invLhs(ndist);
        ArrayXf uhat(ndist);
        ArrayXf logDelta(ndist);
        ArrayXf probDelta(ndist);

        unsigned delta;

        // shuffling the SNP index for faster convergence
        vector<int> snpIndexVec = Gadget::shuffle_index(chrStart, chrEnd);

        //for (unsigned i=chrStart; i<=chrEnd; ++i) {
        for (unsigned t = 0; t < chrSize; t++) {
            unsigned i = snpIndexVec[t];
            if (badSnps[i]) {
                valuesPtr[i] = 0.0;
                continue;
            }
            if (incdSnpInfoVec[i]->fitAsFixedEffect) continue;
            oldSample = valuesPtr[i];

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
            pipPtr[i] = 1.0f - probDelta[0];

            delta = bernoulli.sample(probDelta, urnd[i]);
            snpsetChr[chr][delta].push_back(i);

            if (delta) {
                valuesPtr[i] = uhat[delta] + nrnd[i]*sqrtf(invLhs[delta]);
                sampleDiff = oldSample - valuesPtr[i];
                for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                    rcorr[it.index()] += it.value() * sampleDiff;
                }
                ssq[chr] += valuesPtr[i] * valuesPtr[i];
                wtdssq[chr] += (valuesPtr[i] * valuesPtr[i]) / gamma[delta];
                ++nnz[chr];
                for(unsigned k2 = 0; k2 < delta ; k2++){
                    z(i, k2) = 1;
                }
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
    wtdSumSq = 0.0;
    numNonZeros = 0;
    nnzPerChr.setZero(numChr);
    numSnpMix.setZero(ndist);
    snpset.resize(ndist);
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
    }
    for (unsigned i=0; i<numChr; ++i) {
        sumSq += ssq[i];
        wtdSumSq += wtdssq[i];
        numNonZeros += nnz[i];
        nnzPerChr[i] = nnz[i];
        for (unsigned k=0; k<ndist; ++k) {
            for (unsigned j=0; j<snpsetChr[i][k].size(); ++j) {
                snpset[k].push_back(snpsetChr[i][k][j]);
                numSnpMix[k]++;
            }
        }
    }
    values = VectorXf::Map(valuesPtr, size);
    pip = VectorXf::Map(pipPtr, size);
}


void ApproxBayesRC::SnpEffects::sampleFromFC_eigen(vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
                                             const vector<LDBlockInfo*> &keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks,
                                             const MatrixXf &snpPi, const VectorXf &gamma, const float varg,
                                             DeltaPi &deltaPi, const bool hsqPercModel, const float sigmaSq, const vector<SnpInfo*> &incdSnpInfoVec){
    // -----------------------------------------
    // This method uses low-rank model with eigen-decomposition of LD matrices
    // -----------------------------------------
    long nBlocks = keptLdBlockInfoVec.size();
    
    whatBlocks.resize(nBlocks);
    for (unsigned i=0; i<nBlocks; ++i) {
        whatBlocks[i].resize(wcorrBlocks[i].size());
    }

    float ssq[nBlocks], wtdssq[nBlocks], s2pq[nBlocks], nnz[nBlocks];
    memset(ssq,0, sizeof(float)*nBlocks);
    memset(wtdssq,0, sizeof(float)*nBlocks);
    memset(s2pq,0,sizeof(float)*nBlocks);
    memset(nnz,0, sizeof(float)*nBlocks);

    float *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads
    float *fcMeanPtr = fcMean.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads

    vector<float> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
    pip.setZero(size);
    membership.setZero(size);
    z.setZero(size, ndist-1);   // indicator variables for conditional membership
    
    ArrayXf wtdSigmaSq(ndist);
    
    if (hsqPercModel && varg) {
        wtdSigmaSq = gamma * 0.01 * varg;
    } else {
        wtdSigmaSq = gamma * sigmaSq;
    }

    ArrayXf invWtdSigmaSq = wtdSigmaSq.inverse();
    ArrayXf logWtdSigmaSq = wtdSigmaSq.log();
    
    MatrixXf logPi = snpPi.array().log().matrix();
        
//    cout << "varg " << varg << endl;
//    cout << "gamma " << gamma.transpose() << endl;
//    cout << "wtdSigmaSq " << wtdSigmaSq.transpose() << endl;
//    cout << "invWtdSigmaSq " << invWtdSigmaSq.transpose() << endl;
//    cout << "logWtdSigmaSq " << logWtdSigmaSq.row(0) << endl;
    
    vector<vector<vector<unsigned> > > snpsetBlocks(nBlocks);
    for (unsigned i=0; i<nBlocks; ++i) {
        snpsetBlocks[i].resize(ndist);
        for (unsigned k=0; k<ndist; ++k) {
            snpsetBlocks[i][k].resize(0);
        }
    }

    for (unsigned k=0; k<ndist; ++k) {
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
        unsigned blockSize  = blockEnd - blockStart + 1;

        float invVareDn = nGWASblocks[blk] / (vareBlocks[blk]);

        ArrayXf invLhs = 1.0/(invVareDn + invWtdSigmaSq);
        ArrayXf logInvLhsMsigma = invLhs.log() - logWtdSigmaSq;

        // shuffling the SNP index for faster convergence
        vector<int> snpIndexVec = Gadget::shuffle_index(blockStart, blockEnd);

        //for(unsigned i = blockStart; i <= blockEnd; i++){
        for (unsigned t = 0; t < blockSize; t++) {
            unsigned i = snpIndexVec[t];
            SnpInfo *snp = blockInfo->snpInfoVec[i-blockStart];
            if (snp->skip) {
                valuesPtr[i] = 0.0;
                pip[i] = 0.0f;
                for (unsigned k = 0; k < (unsigned)ndist; ++k)
                    deltaPi[k]->values[i] = (k == 0) ? 1.0f : 0.0f;
                continue;
            }
            if (badSnps[i]) {
                valuesPtr[i] = 0.0;
                continue;
            }
            if (snp->fitAsFixedEffect) continue;
            float oldSample = valuesPtr[i];
            int qcolRC = eigenQColIndex(blockInfo, i);
            if (qcolRC < 0) continue;
            Ref<const VectorXf> Qi = Q.col(qcolRC);
            float rhs = (Qi.dot(wcorr) + oldSample)*invVareDn;
            ArrayXf uhat = invLhs * rhs;
            ArrayXf logDelta = 0.5*(logInvLhsMsigma + uhat*rhs) + logPi.row(i).transpose().array();
            logDelta[0] = logPi(i,0);
            
//            if (i==0) cout << i << " rhs " << rhs << " invVareDn " << invVareDn << " invWtdSigmaSq " << invWtdSigmaSq.transpose() << " uhat " << uhat.transpose() << endl;
            
            ArrayXf probDelta(ndist);
            for (unsigned k=0; k<ndist; ++k) {
                probDelta[k] = 1.0f/(logDelta-logDelta[k]).exp().sum();
                if(isnan(probDelta[k])) probDelta[k] = 0;
                deltaPi[k]->values[i] = probDelta[k];
            }
            pip[i] = 1.0f - probDelta[0];
            
            //            unsigned delta;
            //            #pragma omp critical
            //            {
            unsigned delta = bernoulli.sample(probDelta, urnd[i]);
            snpsetBlocks[blk][delta].push_back(i);
            membership[i] = delta;
            //            }
            
            if (delta) {
                valuesPtr[i] = uhat[delta] + nrnd[i]*sqrtf(invLhs[delta]);
                wcorr += Qi*(oldSample - valuesPtr[i]);
                what  += Qi* valuesPtr[i];
                ssq[blk] += valuesPtr[i] * valuesPtr[i];
                wtdssq[blk] += (valuesPtr[i] * valuesPtr[i]) / gamma[delta];
                ++nnz[blk];
                for(unsigned k2 = 0; k2 < delta ; k2++){
                    z(i, k2) = 1;
                }
            }
            else {
                if (oldSample) wcorr += Qi * oldSample;
                valuesPtr[i] = 0.0;
            }
            
            uhat[0] = 0.0;
            fcMeanPtr[i] = (uhat * probDelta).sum();  // full conditional mean
        }

    }
    for (unsigned i = 0; i < size; ++i) {
        if (incdSnpInfoVec[i]->fitAsFixedEffect)
            fcMeanPtr[i] = valuesPtr[i];
    }
    
    // ---------------------------------------------------------------------
    // Tally up the effect sum of squares and the number of non-zero effects
    // ---------------------------------------------------------------------
    sumSq = 0.0;
    wtdSumSq = 0.0;
    numNonZeros = 0;
    nnzPerBlk.setZero(nBlocks);
    ssqBlocks.setZero(nBlocks);
    numSnpMix.setZero(ndist);
    snpset.resize(ndist);
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
    }
    for (unsigned blk=0; blk<nBlocks; ++blk) {
        sumSq += ssq[blk];
        wtdSumSq += wtdssq[blk];
        numNonZeros += nnz[blk];
        nnzPerBlk[blk] = nnz[blk];
        ssqBlocks[blk] = ssq[blk];
        for (unsigned k=0; k<ndist; ++k) {
            for (unsigned j=0; j<snpsetBlocks[blk][k].size(); ++j) {
                snpset[k].push_back(snpsetBlocks[blk][k][j]);
                numSnpMix[k]++;
            }
        }
    }
    values = VectorXf::Map(valuesPtr, size);
    fcMean = VectorXf::Map(fcMeanPtr, size);
}

//void ApproxBayesRC::SnpEffects::sampleFromTGS_eigen(const vector<vector<int> > &selectedSnps, vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
//                                             const vector<LDBlockInfo*> &keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks,
//                                             const MatrixXf &snpPi, const VectorXf &gamma, const float varg,
//                                                   const bool hsqPercModel, const float sigmaSq){
//    // -----------------------------------------
//    // This method uses tempered Gibbs sampler to improve mixing
//    // -----------------------------------------
//    
//    unsigned numSelSnp = selectedSnps.size();
//    
//    unsigned niter = numSelSnp;
//    
//    ArrayXf wtdSigmaSq(ndist);
//    if (hsqPercModel && varg) {
//        wtdSigmaSq = gamma * 0.01 * varg;
//    } else {
//        wtdSigmaSq = gamma * sigmaSq;
//    }
//    ArrayXf invWtdSigmaSq = wtdSigmaSq.inverse();
//    ArrayXf logWtdSigmaSq = wtdSigmaSq.log();
//    MatrixXf logPi = snpPi.array().log().matrix();
//    
//    MatrixXf probDelta(numSelSnp, ndist);
//    VectorXi delta(numSelSnp);
//    VectorXf probDelta_current(numSelSnp);
//    VectorXf p_delta(numSelSnp);
//    VectorXf selSnpIndices(numSelSnp);
//    VectorXf selSnpBlkIdx(numSelSnp);
//    
//    // compute full conditional probabilities
//    for (unsigned i=0; i<numSelSnp; ++i) {
//        unsigned chr = selectedSnps[i][0];
//        unsigned blk = selectedSnps[i][1];
//        unsigned snp = selectedSnps[i][2];
//        unsigned blockStart = keptLdBlockInfoVec[blk]->startSnpIdx;
//        float invVareDn = nGWASblocks[blk] / (vareBlocks[blk]);
//        ArrayXf invLhs = 1.0/(invVareDn + invWtdSigmaSq);
//        ArrayXf logInvLhsMsigma = invLhs.log() - logWtdSigmaSq;
//
//        float oldSample = values[snp];
//        Ref<const VectorXf> Qi = Qblocks[blk].col(snp - blockStart);
//        Ref<VectorXf> wcorr = wcorrBlocks[blk];
//        float rhs = (Qi.dot(wcorr) + oldSample)*invVareDn;
//        ArrayXf uhat = invLhs * rhs;
//        ArrayXf logDelta = 0.5*(logInvLhsMsigma + uhat*rhs) + logPi.row(snp).transpose().array();
//        logDelta[0] = logPi(snp,0);
//        for (unsigned k=0; k<ndist; ++k) {
//            probDelta(i,k) = 1.0f/(logDelta-logDelta[k]).exp().sum();
//            if(isnan(probDelta(i,k))) probDelta(i,k) = 0;
//        }
//        
//        delta[i] = membership[snp];
//        selSnpIndices[i] = snp;
//        selSnpBlkIdx[i] = blk;
//    }
//    
//    VectorXf weight(niter);
//    VectorXf probDelta0_sum;
//    probDelta0_sum.setZero(numSelSnp);
//
//    for (unsigned t=0; t<niter; ++t) {
//        // get full conditional probabilities for current delta
//        for (unsigned i=0; i<numSelSnp; ++i) {
//            probDelta_current[i] = probDelta(i,delta[i]);
//        }
//        
//        // compute p_delta = g(delta|else)/f(delta|else)
//        p_delta = 1.0/ndist/probDelta_current.array();
//        
//        // sample focal SNP
//        unsigned focal_snp;
//        vector<unsigned> zeroIndices;
//        for (unsigned i=0; i<numSelSnp; ++i) {
//            if (probDelta_current[i] == 0) {
//                zeroIndices.push_back(i);
//            }
//        }
//        unsigned numZero = zeroIndices.size();
//        VectorXf probVec;
//        probVec.setZero(numSelSnp);
//        if (numZero){
//            for (unsigned i=0; i<numZero; ++i) {
//                probVec[zeroIndices[i]] = 1/float(numZero);
//            }
//        } else {
//            probVec = p_delta/p_delta.sum();
//        }
//        focal_snp = bernoulli.sample(probVec);
//        unsigned focal_snp_idx = selSnpIndices[focal_snp];
//        unsigned focal_snp_blk = selSnpBlkIdx[focal_snp];
//        
//        // sample delta_focal from uniform distribution
//        unsigned delta_focal_old = delta[focal_snp];
//        VectorXf otherDeltaStates(ndist-1);
//        probVec.resize(ndist-1);
//        for (unsigned k=0, idx=0; k<ndist; ++k) {
//            if (k != delta_focal_old) {
//                otherDeltaStates[idx] = k;
//                probVec[idx] = 1/float(ndist-1);
//                ++idx;
//            }
//        }
//        delta[focal_snp] = otherDeltaStates[bernoulli.sample(probVec)];
//        
//        // update probDelta
//        if (0 == delta_focal_old || 0 == delta[focal_snp]) {  // update probDelta for SNPs in the same block when the old or new delta = 0
//            unsigned blockStart = keptLdBlockInfoVec[focal_snp_blk]->startSnpIdx;
//            Ref<const VectorXf> Q_focal = Qblocks[focal_snp_blk].col(focal_snp_idx - blockStart);
//            float beta_focal = values[focal_snp_idx];
//            
//            for (unsigned i=0; i<numSelSnp; ++i) {
//                unsigned chr = selectedSnps[i][0];
//                unsigned blk = selectedSnps[i][1];
//                unsigned snp = selectedSnps[i][2];
//                
//                if (blk != focal_snp_blk) continue;
//                
//                blockStart = keptLdBlockInfoVec[blk]->startSnpIdx;
//                float invVareDn = nGWASblocks[blk] / (vareBlocks[blk]);
//                ArrayXf invLhs = 1.0/(invVareDn + invWtdSigmaSq);
//                ArrayXf logInvLhsMsigma = invLhs.log() - logWtdSigmaSq;
//
//                float oldSample = values[snp];
//                Ref<const VectorXf> Qi = Qblocks[blk].col(snp - blockStart);
//                Ref<VectorXf> wcorr = wcorrBlocks[blk];
//                float rhs = Qi.dot(wcorr) + oldSample;
//                if (0 == delta_focal_old) {
//                    rhs -= Qi.dot(Q_focal)*beta_focal;
//                } else {
//                    rhs += Qi.dot(Q_focal)*beta_focal;
//                }
//                rhs *= invVareDn;
//                ArrayXf uhat = invLhs * rhs;
//                ArrayXf logDelta = 0.5*(logInvLhsMsigma + uhat*rhs) + logPi.row(snp).transpose().array();
//                logDelta[0] = logPi(snp,0);
//                for (unsigned k=0; k<ndist; ++k) {
//                    probDelta(i,k) = 1.0f/(logDelta-logDelta[k]).exp().sum();
//                    if(isnan(probDelta(i,k))) probDelta(i,k) = 0;
//                }
//            }
//        }
//        else { // change state from one nonzero component to another nonzero component only affect the focal SNP
//            unsigned blockStart = keptLdBlockInfoVec[focal_snp_blk]->startSnpIdx;
//            Ref<const VectorXf> Q_focal = Qblocks[focal_snp_blk].col(focal_snp_idx - blockStart);
//            Ref<VectorXf> wcorr = wcorrBlocks[focal_snp_blk];
//            float invVareDn = nGWASblocks[focal_snp_blk] / (vareBlocks[focal_snp_blk]);
//            ArrayXf invLhs = 1.0/(invVareDn + invWtdSigmaSq);
//            ArrayXf logInvLhsMsigma = invLhs.log() - logWtdSigmaSq;
//
//            float oldSample = values[focal_snp_idx];
//            float rhs = (Q_focal.dot(wcorr) + oldSample)*invVareDn;
//            ArrayXf uhat = invLhs * rhs;
//            ArrayXf logDelta = 0.5*(logInvLhsMsigma + uhat*rhs) + logPi.row(focal_snp_idx).transpose().array();
//            logDelta[0] = logPi(focal_snp_idx,0);
//            for (unsigned k=0; k<ndist; ++k) {
//                probDelta(focal_snp,k) = 1.0f/(logDelta-logDelta[k]).exp().sum();
//                if(isnan(probDelta(focal_snp,k))) probDelta(focal_snp,k) = 0;
//            }
//        }
//        
//        // get full conditional probabilities for current delta
//        for (unsigned i=0; i<numSelSnp; ++i) {
//            probDelta_current[i] = probDelta(i,delta[i]);
//        }
//        
//        // update p_delta = g(delta|else)/f(delta|else)
//        p_delta = 1.0/ndist/probDelta_current.array();
//        
//        // compute weight
//        float sum_p_delta = 0.0;
//        for (unsigned i=0; i<numSelSnp; ++i) {
//            if (std::isfinite(p_delta[i])) sum_p_delta += p_delta[i];
//        }
//        weight[t] = float(numSelSnp)/sum_p_delta;
//        
//        // update pi0 = 1 - PIP
//        probDelta0_sum += weight[t]*probDelta.col(0);
//        
//        // update focal SNP effect
//        unsigned blockStart = keptLdBlockInfoVec[focal_snp_blk]->startSnpIdx;
//        Ref<const VectorXf> Q_focal = Qblocks[focal_snp_blk].col(focal_snp_idx - blockStart);
//        Ref<VectorXf> wcorr = wcorrBlocks[focal_snp_blk];
//        float invVareDn = nGWASblocks[focal_snp_blk] / (vareBlocks[focal_snp_blk]);
//        ArrayXf invLhs = 1.0/(invVareDn + invWtdSigmaSq);
//        ArrayXf logInvLhsMsigma = invLhs.log() - logWtdSigmaSq;
//        float oldSample = values[focal_snp_idx];
//        float rhs = (Q_focal.dot(wcorr) + oldSample)*invVareDn;
//        ArrayXf uhat = invLhs * rhs;
//        
//        if (delta[focal_snp]) {
//            values[focal_snp_idx] = normal.sample(uhat[delta[focal_snp]], invLhs[delta[focal_snp]]);
//            wcorr += Q_focal*(oldSample - values[focal_snp_idx]);
//            for(unsigned k2 = 0; k2 < delta[focal_snp]; k2++){
//                z(focal_snp_idx, k2) = 1;
//            }
//        }
//        else {
//            if (oldSample) wcorr += Q_focal * oldSample;
//            values[focal_snp_idx] = 0.0;
//        }
//    }
//    
//    
//    // update PIPs for the selected SNPs
//    float weight_sum = weight.sum();
//    for (unsigned i=0; i<numSelSnp; ++i) {
//        unsigned chr = selectedSnps[i][0];
//        unsigned blk = selectedSnps[i][1];
//        unsigned snp = selectedSnps[i][2];
//        pip[snp] = 1.0 - probDelta0_sum[i]/weight_sum;
//    }
//    
//}

void ApproxBayesRC::SnpEffects::sampleFromTGS_eigen(vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
                                                    const map<SnpInfo*, vector<SnpInfo*> > &LDmap, const vector<LDBlockInfo*> &keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks,
                                                    const MatrixXf &snpPi, const VectorXf &gamma, const float varg,
                                                    DeltaPi &deltaPi, const bool hsqPercModel, const float sigmaSq){
    // -----------------------------------------
    // This method uses tempered Gibbs sampler to improve mixing
    // Apply to the LD friend set of each SNP with sampled nonzero effect and running PIP > 0.01
    // -----------------------------------------
    
    ArrayXf wtdSigmaSq(ndist);
    if (hsqPercModel && varg) {
        wtdSigmaSq = gamma * 0.01 * varg;
    } else {
        wtdSigmaSq = gamma * sigmaSq;
    }
    ArrayXf invWtdSigmaSq = wtdSigmaSq.inverse();
    ArrayXf logWtdSigmaSq = wtdSigmaSq.log();
    
//    cout << "LDmap size " << LDmap.size() << endl;
    unsigned cnt_tmp = 0;
    
    map<SnpInfo*, vector<SnpInfo*> >::const_iterator it, end = LDmap.end();
    for (it = LDmap.begin(); it != end; ++it) {
        SnpInfo *snpInfo = it->first;
        if (snpInfo->fitAsFixedEffect) continue;
        unsigned snpIdx = snpInfo->index;
                
//        cout << cnt_tmp++ << " " << snpIdx << " " << snpInfo->ID << endl;
//        cout << "membership " << membership[snpIdx] << " membership_szie " << membership.size() << " posteriorMeanPIP " << posteriorMeanPIP[snpIdx] << " posteriorMeanPIP_size " << posteriorMeanPIP.size()<< endl;

        if (!membership[snpIdx]) continue;   // skip zero effect SNPs
        if (posteriorMeanPIP[snpIdx] < 0.05) continue;  // skip uninformative SNPs
        
        unsigned blockIdx = snpInfo->blockIdx;
        
        vector<int> selectedSnps;
        selectedSnps.push_back(snpIdx);
        
        for (unsigned j=0; j<it->second.size(); ++j) {
            selectedSnps.push_back(it->second[j]->index);
        }
        
        unsigned numSelSnp = selectedSnps.size();
        
//        cout << snpIdx << endl;
//        for (unsigned i=0; i<numSelSnp; ++i) {
//            cout << selectedSnps[i] << " ";
//        }
//        cout << endl;
        
        unsigned niter = numSelSnp;
        
        MatrixXf logPi(numSelSnp, ndist);
        for (unsigned j=0; j<numSelSnp; ++j) {
            unsigned j_idx = selectedSnps[j];
            logPi.row(j) = snpPi.row(j_idx).array().log();
//            cout << j << " logPi.row(j) " << logPi.row(j) << endl;
        }
        
        MatrixXf probDelta(numSelSnp, ndist);
        VectorXi delta(numSelSnp);
        VectorXf probDelta_current(numSelSnp);
        VectorXf p_delta(numSelSnp);
        VectorXf weight(niter);
        MatrixXf probDelta_sum;
        probDelta_sum.setZero(numSelSnp, ndist);
        
//        cout << "blockIdx " << blockIdx << endl;
        
        Ref<const MatrixXf> Q = Qblocks[blockIdx];
        Ref<VectorXf> wcorr = wcorrBlocks[blockIdx];
        unsigned blockStart = keptLdBlockInfoVec[blockIdx]->startSnpIdx;
        const LDBlockInfo *blockInfo = keptLdBlockInfoVec[blockIdx];
        float invVareDn = nGWASblocks[blockIdx] / vareBlocks[blockIdx];
        ArrayXf invLhs = 1.0/(invVareDn + invWtdSigmaSq);
        ArrayXf logInvLhsMsigma = invLhs.log() - logWtdSigmaSq;
        
        // compute full conditional probabilities
        for (unsigned i=0; i<numSelSnp; ++i) {
            unsigned snpIdx = selectedSnps[i];
//            cout << i << " " << snpIdx << endl;
            for (unsigned k=0; k<ndist; ++k) {
//                cout << k << " " << deltaPi[k]->values[snpIdx] << endl;
                probDelta(i,k) = deltaPi[k]->values[snpIdx];
            }
//            cout << i << " probDelta.row(i) " << probDelta.row(i) << endl;

            delta[i] = membership[snpIdx];
        }
        
        
        // TGS begins
        for (unsigned t=0; t<niter; ++t) {
            // get full conditional probabilities for current delta
            for (unsigned i=0; i<numSelSnp; ++i) {
                probDelta_current[i] = probDelta(i,delta[i]);
            }

            // compute p_delta = g(delta|else)/f(delta|else)
            p_delta = 1.0/ndist/probDelta_current.array();
            
            // sample focal SNP
            unsigned focal_snp;
            vector<unsigned> zeroIndices;
            for (unsigned i=0; i<numSelSnp; ++i) {
                if (probDelta_current[i] == 0) {
                    zeroIndices.push_back(i);
                }
            }
            unsigned numZero = zeroIndices.size();
            VectorXf probVec;
            probVec.setZero(numSelSnp);
            if (numZero){
                for (unsigned i=0; i<numZero; ++i) {
                    probVec[zeroIndices[i]] = 1/float(numZero);
                }
            } else {
                probVec = p_delta/p_delta.sum();
            }
            focal_snp = bernoulli.sample(probVec);
            unsigned focal_snp_idx = selectedSnps[focal_snp];
            
            // sample delta_focal from uniform distribution
            unsigned delta_focal_old = delta[focal_snp];
            VectorXf otherDeltaStates(ndist-1);
            probVec.resize(ndist-1);
            for (unsigned k=0, idx=0; k<ndist; ++k) {
                if (k != delta_focal_old) {
                    otherDeltaStates[idx] = k;
                    probVec[idx] = 1/float(ndist-1);
                    ++idx;
                }
            }
            delta[focal_snp] = otherDeltaStates[bernoulli.sample(probVec)];
            
            // update probDelta
            int qcolF = eigenQColIndex(blockInfo, focal_snp_idx);
            if (qcolF < 0) continue;
            Ref<const VectorXf> Q_focal = Q.col(qcolF);
            float old_beta_focal = values[focal_snp_idx];
            float rhs_focal = (Q_focal.dot(wcorr) + old_beta_focal)*invVareDn;
            ArrayXf uhat_focal = invLhs * rhs_focal;
            ArrayXf logDelta = 0.5*(logInvLhsMsigma + uhat_focal*rhs_focal) + logPi.row(focal_snp).transpose().array();
            logDelta[0] = logPi(focal_snp,0);
            for (unsigned k=0; k<ndist; ++k) {
                probDelta(focal_snp,k) = 1.0f/(logDelta-logDelta[k]).exp().sum();
                if(isnan(probDelta(focal_snp,k))) probDelta(focal_snp,k) = 0;
            }
            if (0 == delta_focal_old || 0 == delta[focal_snp]) {  // update probDelta for other SNPs in the set when the old or new delta = 0
                for (unsigned i=0; i<numSelSnp; ++i) {
                    if (i == focal_snp) continue;
                    unsigned snpIdx = selectedSnps[i];
                    float oldSample = values[snpIdx];
                    int qcolT = eigenQColIndex(blockInfo, snpIdx);
                    if (qcolT < 0) continue;
                    Ref<const VectorXf> Qi = Q.col(qcolT);
                    float rhs = Qi.dot(wcorr) + oldSample;
                    if (0 == delta_focal_old) {
                        rhs -= Qi.dot(Q_focal)*old_beta_focal;
                    } else {
                        rhs += Qi.dot(Q_focal)*old_beta_focal;
                    }
                    rhs *= invVareDn;
                    ArrayXf uhat = invLhs * rhs;
                    ArrayXf logDelta = 0.5*(logInvLhsMsigma + uhat*rhs) + logPi.row(i).transpose().array();
                    logDelta[0] = logPi(i,0);
                    for (unsigned k=0; k<ndist; ++k) {
                        probDelta(i,k) = 1.0f/(logDelta-logDelta[k]).exp().sum();
                        if(isnan(probDelta(i,k))) probDelta(i,k) = 0;
                    }
                }
            }
            
            // get full conditional probabilities for current delta
            for (unsigned i=0; i<numSelSnp; ++i) {
                probDelta_current[i] = probDelta(i,delta[i]);
            }
            
            // update p_delta = g(delta|else)/f(delta|else)  and  compute weight
            if ((probDelta_current.array() > 0).all()) {
                p_delta = 1.0/ndist/probDelta_current.array();
                weight[t] = float(numSelSnp)/p_delta.sum();
            } else {
                weight[t] = 0;
            }
//            if (!std::isfinite(weight[t])) {
//                cout << "focal_snp_idx " << focal_snp_idx << endl;
//                cout << "delta_focal_old " << delta_focal_old << " delta_focal_new " << delta[focal_snp] << endl;
//                ArrayXf logDelta = 0.5*(logInvLhsMsigma + uhat_focal*rhs_focal) + logPi.row(focal_snp).transpose().array();
//                logDelta[0] = logPi(focal_snp,0);
//                for (unsigned k=0; k<ndist; ++k) {
//                    probDelta(focal_snp,k) = 1.0f/(logDelta-logDelta[k]).exp().sum();
//                    if(isnan(probDelta(focal_snp,k))) probDelta(focal_snp,k) = 0;
//                    cout << "k " << k << " probDelta(focal_snp,k) " << probDelta(focal_snp,k) << " logDelta[k] " << logDelta[k] << endl;
//                }
//
//                for (unsigned j=0; j<numSelSnp; ++j) {
//                    cout << selectedSnps[j] << " ";
//                }
//                cout << endl;
//                cout << "delta " << delta.transpose() << endl;
//                cout << "p_delta " << p_delta.transpose() << endl;
//                cout << "probDelta_current " << probDelta_current.transpose() << endl;
//            }
            
            // update full conditional values of pi
            probDelta_sum += weight[t] * probDelta;
            
            // update focal SNP effect
            if (delta[focal_snp]) {
                values[focal_snp_idx] = normal.sample(uhat_focal[delta[focal_snp]], invLhs[delta[focal_snp]]);
                wcorr += Q_focal*(old_beta_focal - values[focal_snp_idx]);
            }
            else {
                if (old_beta_focal) wcorr += Q_focal * old_beta_focal;
                values[focal_snp_idx] = 0.0;
            }
        }
        
        // update PIPs for the selected SNPs
        float weight_sum = weight.sum();
        if (weight_sum) probDelta_sum /= weight_sum;
        for (unsigned i=0; i<numSelSnp; ++i) {
            unsigned snpIdx = selectedSnps[i];
            membership[snpIdx] = delta[i];
            if (weight_sum) pip[snpIdx] = 1.0 - probDelta_sum(i,0);
            for (unsigned k=0; k<ndist; ++k) {
                if (weight_sum) deltaPi[k]->values[snpIdx] = probDelta_sum(i,k);
            }
            if (delta[i]) {
                for(unsigned k2 = 0; k2 < delta[i]; k2++){
                    z(snpIdx, k2) = 1;
                }
            } else {
                z.row(snpIdx) *= 0;
            }
        }
    }
}


void ApproxBayesRC::AnnoEffects::sampleFromFC_Gibbs(MatrixXf &z, const MatrixXf &annoMat, const VectorXf &sigmaSq, MatrixXf &snpP,
                                                    const vector<SnpInfo*> *incdSnpInfoVec, const VectorXf *baselineSnpP) {
    auto inAnno = [&](unsigned j) -> bool {
        return !incdSnpInfoVec || !(*incdSnpInfoVec)[j]->fitAsFixedEffect;
    };

    unsigned numSnps = z.rows();
    if (incdSnpInfoVec && baselineSnpP && baselineSnpP->size() == (Eigen::Index)numComp) {
        unsigned nAnno = 0;
        for (unsigned j = 0; j < numSnps; ++j) {
            if (inAnno(j)) ++nAnno;
        }
        if (nAnno == 0) {
            for (unsigned j = 0; j < numSnps; ++j) {
                for (unsigned t = 0; t < numComp; ++t) {
                    snpP(j, t) = (*baselineSnpP)[t];
                }
            }
            return;
        }
    }

    VectorXf numOnes(numComp);
    #pragma omp parallel for schedule(dynamic)
    for (unsigned i=0; i<numComp; ++i) {
        float s = 0.f;
        for (unsigned j=0; j<numSnps; ++j) {
            if (!inAnno(j)) continue;
            s += z(j,i);
        }
        numOnes[i] = s;
    }

    MatrixXf annoMatPO;
    for (unsigned i=0; i<numComp; ++i) {
        VectorXf &alphai = (*this)[i]->values;
        VectorXf y, zi;
        unsigned numDP;
        if (i==0) {
            numDP = 0;
            for (unsigned j=0; j<numSnps; ++j) {
                if (inAnno(j)) ++numDP;
            }
        } else {
            numDP = (unsigned)numOnes[i-1];
        }

        if(numDP == 0){
            alphai.setZero();
            alphai[0] = -10.0;
            ssq[i] = 0;
        }else{
            y.setZero(numDP);
            zi.setZero(numDP);
            const MatrixXf *annotMatP;
            if (i==0) {
                annoMatPO.resize(numDP, numAnno);
                for (unsigned j=0, idx=0; j<numSnps; ++j) {
                    if (!inAnno(j)) continue;
                    annoMatPO.row(idx) = annoMat.row(j);
                    zi[idx] = z(j,i);
                    ++idx;
                }
                annotMatP = &annoMatPO;
            } else {
                annoMatPO.setZero(numDP, numAnno);
                for (unsigned j=0, idx=0; j<numSnps; ++j) {
                    if (!inAnno(j)) continue;
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
            if (i==0) {
                annoDiagi[0] = float(annoMati.rows());
                #pragma omp parallel for schedule(dynamic)
                for (unsigned k=1; k<numAnno; ++k) {
                    annoDiagi[k] = annoMati.col(k).squaredNorm();
                }
            } else {
                annoDiagi[0] = numOnes[i-1];
                #pragma omp parallel for schedule(dynamic)
                for (unsigned k=1; k<numAnno; ++k) {
                    annoDiagi[k] = annoMati.col(k).squaredNorm();
                }
            }

            VectorXf mean = annoMati * alphai;

            for (unsigned j=0; j<numDP; ++j) {
                if (zi[j]) y[j] = TruncatedNormal::sample_lower_truncated(mean[j], 1.0, 0.0);
                else y[j] = TruncatedNormal::sample_upper_truncated(mean[j], 1.0, 0.0);
            }

            y -= mean;

            float oldSample = alphai[0];
            float rhs = y.sum() + annoDiagi[0]*oldSample;
            float invLhs = 1.0/annoDiagi[0];
            float ahat = invLhs*rhs;
            alphai[0] = Normal::sample(ahat, invLhs);
            y.array() += oldSample - alphai[0];

            vector<int> shuffled_index = Gadget::shuffle_index(1, numAnno-1);

            ssq[i] = 0;
            for (unsigned t=0; t<shuffled_index.size(); ++t) {
                unsigned k = shuffled_index[t];
                oldSample = alphai[k];
                rhs = annoMati.col(k).dot(y) + annoDiagi[k]*oldSample;
                invLhs = 1.0/(annoDiagi[k] + 1.0/sigmaSq[i]);
                ahat = invLhs*rhs;
                alphai[k] = Normal::sample(ahat, invLhs);
                y += annoMati.col(k) * (oldSample - alphai[k]);
                ssq[i] += alphai[k] * alphai[k];
            }
        }

        #pragma omp parallel for schedule(dynamic)
        for (unsigned j=0; j<numSnps; ++j) {
            if (!inAnno(j)) continue;
            snpP(j,i) = Normal::cdf_01(annoMat.row(j).dot(alphai));
        }
    }

    if (incdSnpInfoVec && baselineSnpP && baselineSnpP->size() == (Eigen::Index)numComp) {
        for (unsigned j = 0; j < numSnps; ++j) {
            if (!(*incdSnpInfoVec)[j]->fitAsFixedEffect) continue;
            for (unsigned t = 0; t < numComp; ++t) {
                snpP(j, t) = (*baselineSnpP)[t];
            }
        }
    }
}

void ApproxBayesRC::AnnoEffects::sampleFromFC_MH(MatrixXf &z, const MatrixXf &annoMat, const VectorXf &sigmaSq, MatrixXf &snpP,
                                                 const vector<SnpInfo*> *incdSnpInfoVec, const VectorXf *baselineSnpP) {
    auto inAnno = [&](unsigned j) -> bool {
        return !incdSnpInfoVec || !(*incdSnpInfoVec)[j]->fitAsFixedEffect;
    };

    unsigned numSnps = z.rows();
    for (unsigned i=0; i<numComp; ++i) {
        VectorXf curr_alpha = (*this)[i]->values;
        VectorXf cand_alpha(numAnno);
        for (unsigned k=0; k<numAnno; ++k) {
            cand_alpha[k] = Normal::sample(curr_alpha[k], varProp[i]);
        }

        float logPriorCurr = -0.5f * ((curr_alpha.squaredNorm() - curr_alpha[0]*curr_alpha[0])/sigmaSq[i]);
        float logPriorCand = -0.5f * ((cand_alpha.squaredNorm() - cand_alpha[0]*cand_alpha[0])/sigmaSq[i]);
        
        float logLikCurr = 0.0;
        float logLikCand = 0.0;
        
        for (unsigned j=0; j<numSnps; ++j) {
            if (!inAnno(j)) continue;
            if (i==0) {
                float curr_p = snpP(j,i);
                float cand_p = 1.0f/(1.0f + expf(- annoMat.row(j).dot(cand_alpha)));
                if (z(j,i)) {
                    logLikCurr += logf(curr_p);
                    logLikCand += logf(cand_p);
                } else {
                    logLikCurr += logf(1.0f - curr_p);
                    logLikCand += logf(1.0f - cand_p);
                }
            }
            else {
                if (z(j,i-1)) {
                    float curr_p = snpP(j,i);
                    float cand_p = 1.0f/(1.0f + expf(- annoMat.row(j).dot(cand_alpha)));
                    if (z(j,i)) {
                        logLikCurr += logf(curr_p);
                        logLikCand += logf(cand_p);
                    } else {
                        logLikCurr += logf(1.0f - curr_p);
                        logLikCand += logf(1.0f - cand_p);
                    }
                }
            }
        }
        
        float logPostCurr = logLikCurr + logPriorCurr;
        float logPostCand = logLikCand + logPriorCand;
        
        if (Stat::ranf() < expf(logPostCand-logPostCurr)) {
            (*this)[i]->values = cand_alpha;
            VectorXf lin = annoMat * cand_alpha;
            for (unsigned j=0; j<numSnps; ++j) {
                if (!inAnno(j)) {
                    if (baselineSnpP && baselineSnpP->size() == (Eigen::Index)numComp)
                        snpP(j,i) = (*baselineSnpP)[i];
                    continue;
                }
                snpP(j,i) = 1.0f/(1.0f + expf(-lin[j]));
            }
            ar[i]->count(1, 0.1, 0.5);
        } else {
            ar[i]->count(0, 0.1, 0.5);
        }

        if (!(ar[i]->cnt % 10)) {
            if      (ar[i]->value < 0.2) varProp[i] *= 0.8f;
            else if (ar[i]->value > 0.5) varProp[i] *= 1.2f;
        }
        
        ssq[i] = (*this)[i]->values.squaredNorm() - (*this)[i]->values[0]*(*this)[i]->values[0];
    }

    if (incdSnpInfoVec && baselineSnpP && baselineSnpP->size() == (Eigen::Index)numComp) {
        for (unsigned j = 0; j < numSnps; ++j) {
            if (!(*incdSnpInfoVec)[j]->fitAsFixedEffect) continue;
            for (unsigned t = 0; t < numComp; ++t) {
                snpP(j, t) = (*baselineSnpP)[t];
            }
        }
    }
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

void ApproxBayesRC::AnnoCondProb::compute_probit(const AnnoEffects &annoEffects, const vector<AnnoInfo*> &annoInfoVec){
    for (unsigned i=0; i<annoEffects.numComp; ++i) {
        for (unsigned j=0; j<annoEffects.numAnno; ++j) {
            VectorXf &alpha = annoEffects[i]->values;
            if (j==0) (*this)[i]->values[j] = Normal::cdf_01(alpha[j]);
            else (*this)[i]->values[j] = Normal::cdf_01(alpha[0] + annoInfoVec[j]->sd*alpha[j]);  // NEW
        }
    }
}

void ApproxBayesRC::AnnoCondProb::compute_logistic(const AnnoEffects &annoEffects, const vector<AnnoInfo*> &annoInfoVec){
//    cout << "computing conditional prob... " << endl;
    for (unsigned i=0; i<numComp; ++i) {
        for (unsigned j=0; j<annoEffects.numAnno; ++j) {
            VectorXf &alpha = annoEffects[i]->values;
            if (j==0) (*this)[i]->values[j] = 1.0/(1.0 + exp(-alpha[j]));
            else (*this)[i]->values[j] = 1.0/(1.0 + exp(- alpha[0] - annoInfoVec[j]->sd*alpha[j]));
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
                    //cout << i << " " << j << " " << k << " " << varj << " " << snpEffects[snpIdx] << " " << ZPy[snpIdx] << " " << rcorr[snpIdx] << endl;
                }
            }
        }
        (*this)[i]->values.array() /= nobs;
    }
}

void ApproxBayesRC::AnnoGenVar::compute(const VectorXf &snpEffects, const vector<unsigned> &membership, const vector<MatrixXf> &Qblocks, const vector<LDBlockInfo*> &keptLdBlockInfoVec, const MatrixXf &annoMat, const vector<AnnoInfo*> &annoInfoVec){
    
    unsigned nBlocks = Qblocks.size();
    vector<vector<vector<float> > > vgBlocks;
    vgBlocks.resize(nBlocks);
    for (unsigned i=0; i<nBlocks; ++i) {
        vgBlocks[i].resize(numComp);
        for (unsigned k=0; k<numComp; ++k) {
            vgBlocks[i][k].resize(numAnno);
        }
    }

#pragma omp parallel for schedule(dynamic)
    for(unsigned blk = 0; blk < nBlocks; blk++){
        Ref<const MatrixXf> Q = Qblocks[blk];
        vector<vector<VectorXf> > whatBlock(numComp);
        for (unsigned k=0; k<numComp; ++k) {
            whatBlock[k].resize(numAnno);
            for (unsigned c=0; c<numAnno; ++c) {
                whatBlock[k][c].resize(Q.rows());
                whatBlock[k][c].setZero();
            }
        }
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[blk];
        
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd   = blockInfo->endSnpIdx;
        
        for (unsigned j = blockStart; j <= blockEnd; j++){
            int qcol = eigenQColIndex(blockInfo, j);
            if (qcol < 0) continue;
            Ref<const VectorXf> Qj = Q.col(qcol);
            unsigned delta = membership[j];
            if (delta) {
                for (unsigned c=0; c<numAnno; ++c) {
                    if (annoInfoVec[c]->isBinary && annoMat(j,c)) {  // for binary annotations
                        whatBlock[delta-1][c] += Qj*snpEffects[j];
                    }
                }
            }
        }
        
        for (unsigned k=0; k<numComp; ++k) {
            for (unsigned c=0; c<numAnno; ++c) {
                if (annoInfoVec[c]->isBinary) {  // for binary annotations
                    vgBlocks[blk][k][c] = whatBlock[k][c].squaredNorm();
                }
            }
        }
    }

    
    for (unsigned k=0; k<numComp; ++k) {
        (*this)[k]->values.setZero(numAnno);
        for (unsigned c=0; c<numAnno; ++c) {
            for(unsigned blk = 0; blk < nBlocks; blk++){
                if (annoInfoVec[c]->isBinary) {  // for binary annotations
                    (*this)[k]->values[c] += vgBlocks[blk][k][c];
                }
            }
        }
    }
}

void ApproxBayesRC::AnnoGenVar::compute(const VectorXf &snpEffects, const vector<vector<unsigned> > &snpset, const MatrixXf &annoMat, const vector<SnpInfo*> *incdSnpInfoVec){
    for (unsigned i=0; i<numComp; ++i) {
        (*this)[i]->values.setZero(numAnno);
        unsigned size = snpset[i+1].size();
        for (unsigned j=0; j<size; ++j) {
            unsigned snpIdx = snpset[i+1][j];
            if (incdSnpInfoVec && (*incdSnpInfoVec)[snpIdx]->fitAsFixedEffect) continue;
            float varj = snpEffects[snpIdx] * snpEffects[snpIdx];
            for (unsigned k=0; k<numAnno; ++k) {
                if (annoMat(snpIdx,k)) {
                    (*this)[i]->values[k] += annoMat(snpIdx,k) * varj;
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

void ApproxBayesRC::AnnoPerSnpHsqEnrichment::compute(const VectorXf &annoTotalGenVar, const vector<AnnoInfo*> &annoInfoVec){
    for (unsigned i=0; i<size; ++i) {
        AnnoInfo *anno = annoInfoVec[i];
        if (anno->isBinary) {
            values[i] = annoTotalGenVar[i]/annoTotalGenVar[0] * invSnpProp[i];
        } else {
            MatrixXf XPX(2,2);
            XPX(0,0) = size;
            XPX(0,1) = XPX(1,0) = anno->sum;
            XPX(1,1) = anno->ssq;
            VectorXf XPy(2);
            XPy(0) = annoTotalGenVar[0];
            XPy(1) = annoTotalGenVar[i];
            VectorXf coef = XPX.householderQr().solve(XPy);
            values[i] = 1.0 + coef(1) / annoTotalGenVar[0] * invSnpProp[0];
        }
    }
    //values = annoTotalGenVar.array()/annoTotalGenVar[0] * invSnpProp.array();  // first is the intercept
//    cout << "AnnoPerSnpHsqEnrichment " << values.transpose() << endl;
}

void ApproxBayesRC::AnnoPerSnpRsqEnrichment::compute(const VectorXf &snpEffectMeans, const MatrixXf &annoMat, const vector<AnnoInfo*> &annoInfoVec){
    VectorXf rsqVec;
    rsqVec.setZero(size);
    unsigned numSnps = snpEffectMeans.size();
    
    for (unsigned k=0; k<size; ++k) {
        AnnoInfo *anno = annoInfoVec[k];
        for (unsigned j=0; j<anno->size; ++j) {
            SnpInfo *snp = anno->memberSnpVec[j];
            rsqVec[k] += annoMat(snp->index,k) * snpEffectMeans[snp->index]*snpEffectMeans[snp->index];
        }
        if (anno->isBinary) {
            values[k] = rsqVec[k]/rsqVec[0] * invSnpProp[k];
        } else {
            MatrixXf XPX(2,2);
            XPX(0,0) = size;
            XPX(0,1) = XPX(1,0) = anno->sum;
            XPX(1,1) = anno->ssq;
            VectorXf XPy(2);
            XPy(0) = rsqVec[0];
            XPy(1) = rsqVec[k];
            VectorXf coef = XPX.householderQr().solve(XPy);
            values[k] = 1.0 + coef(1) / rsqVec[0] * invSnpProp[0];
        }
    }
}

//void ApproxBayesRC::AnnoPerSnpHsqEnrichment::compute(const VectorXf &snpEffects, const VectorXf &annoTotalGenVar, const MatrixXf &annoMat, const MatrixXf &APA, const vector<AnnoInfo*> &annoInfoVec){
//    VectorXf y = snpEffects.array().square();
//    VectorXf APy = annoMat.transpose() * y;
//    VectorXf coef = APA.householderQr().solve(APy);
//    values = coef.array() * annoTotalGenVar.array().inverse * invSnpProp.array();
////    for (unsigned i=0; i<size; ++i) {
////        AnnoInfo *anno = annoInfoVec[i];
////        if (anno->isBinary) {
////            values[i] = ;
////        } else {
////            values[i]
////        }
////    }
//}

void ApproxBayesRC::AnnoJointPerSnpHsqEnrichment::compute(const AnnoJointProb &annoJointProb, const vector<AnnoInfo*> &annoInfoVec, const VectorXf &gamma, const float varg, const bool hsqPercModel, const float sigmaSq){
    // assuming annotation has been mean centred and variance standardised
    
    // compute alpha values for joint probability pi
    unsigned numDist = annoJointProb.numDist;
    MatrixXf alphaForPi(numDist, size);
    for (unsigned i=0; i<size; ++i) {
        AnnoInfo *anno = annoInfoVec[i];
        for (unsigned k=0; k<numDist; ++k) {
            if (i == 0) {
                alphaForPi(k,i) =  Normal::quantile_01(annoJointProb[k]->values[i]);
            } else {
                alphaForPi(k,i) = (Normal::quantile_01(annoJointProb[k]->values[i]) - Normal::quantile_01(annoJointProb[k]->values[0])) / anno->sd;
            }
        }
    }
    
    // compute parititioned per-SNP hsq components
    VectorXf hsqPartition;
    hsqPartition.setZero(size);
    for (unsigned i=0; i<size; ++i) {
        AnnoInfo *anno = annoInfoVec[i];
        if (i == 0) {
            for (unsigned k=0; k<numDist; ++k) {
                if (hsqPercModel) {
                    //hsqPartition[i] += Normal::cdf_01(alphaForPi(k,0)) * gamma[k] * 0.01*varg;
                    hsqPartition[i] += annoJointProb[k]->values[i] * gamma[k] * 0.01*varg;
                } else {
                    //hsqPartition[i] += Normal::cdf_01(alphaForPi(k,0)) * gamma[k] * sigmaSq;
                    hsqPartition[i] += annoJointProb[k]->values[i] * gamma[k] * sigmaSq;
                }
            }
        } else {
            for (unsigned k=0; k<numDist; ++k) {
                if (hsqPercModel) {
                    //float dev = anno->mean * Normal::pdf_01(alphaForPi(k,0))*alphaForPi(k,i) * gamma[k] * 0.01*varg;
                    //dev += anno->sd*anno->sd * (-0.5)*alphaForPi(k,0)*Normal::pdf_01(alphaForPi(k,0))*alphaForPi(k,i)*alphaForPi(k,i) * gamma[k] * 0.01*varg;
                    //if (std::isfinite(dev)) hsqPartition[i] += dev;
                    hsqPartition[i] += annoJointProb[k]->values[i] * gamma[k] * 0.01*varg;
                } else {
                    //float dev = anno->mean * Normal::pdf_01(alphaForPi(k,0))*alphaForPi(k,i) * gamma[k] * sigmaSq;
                    //dev += anno->sd*anno->sd * (-0.5)*alphaForPi(k,0)*Normal::pdf_01(alphaForPi(k,0))*alphaForPi(k,i)*alphaForPi(k,i) * gamma[k] * sigmaSq;
                    //if (std::isfinite(dev)) hsqPartition[i] += dev;
                    hsqPartition[i] += annoJointProb[k]->values[i] * gamma[k] * sigmaSq;
                }
            }
        }
    }
    
//    cout << "alphaForPi\n" << alphaForPi << endl;
//    cout << "hsqPartition\n" << hsqPartition.transpose() << endl;
//    cout << "enrich\n" << values.transpose() << endl;
//    cout << "anno0 " << annoInfoVec[0]->sum << " " << annoInfoVec[0]->sd << " " << annoInfoVec[0]->ssq << endl;
//    cout << "anno1 " << annoInfoVec[1]->sum << " " << annoInfoVec[1]->sd << " " << annoInfoVec[1]->ssq << endl;

    // compute enrichment
    for (unsigned i=0; i<size; ++i) {
        if (i == 0) {
            values[i] = 1.0;
        } else {
            if (hsqPartition[0] == 0) values[i] = 1.0;
//            else values[i] = 1.0 + hsqPartition[i]/hsqPartition[0];
            else values[i] = hsqPartition[i]/hsqPartition[0];
        }
    }
}

void ApproxBayesRC::AnnoJointPerSnpHsqEnrichment::compute(const VectorXf &snpEffects, const MatrixXf &annoMat, const AnnoJointProb &annoJointProb, const vector<AnnoInfo*> &annoInfoVec, const VectorXf &gamma, const VectorXf &snpAnnoCntInv, const vector<SnpInfo*> *incdSnpInfoVec){
    unsigned numAnno = annoMat.cols();
    unsigned numSnps = annoMat.rows();
    unsigned numDist = annoJointProb.numDist;
    VectorXf hsqAnnoObs;
    VectorXf hsqAnnoNull;
    hsqAnnoObs.setZero(numAnno);
    hsqAnnoNull.setZero(numAnno);
    
    float betaSqNull = snpEffects.squaredNorm()/float(numSnps);
    VectorXf probAnnoObs(numAnno);
    for (unsigned j=0; j<numSnps; ++j) {
        if (incdSnpInfoVec && (*incdSnpInfoVec)[j]->fitAsFixedEffect) continue;
        if (snpEffects[j]) {
            probAnnoObs.setZero(numAnno);
            for (unsigned i=0; i<numAnno; ++i) {
                if (annoMat(j,i)) {
                    for (unsigned k=1; k<numDist; ++k) {  // loop over nonzero distributions
                        probAnnoObs[i]  += annoJointProb[k]->values[i] * gamma[k];
                    }
                }
            }
            probAnnoObs  /= probAnnoObs.sum();
            float betajSq = snpEffects[j]*snpEffects[j];
            for (unsigned i=0; i<numAnno; ++i) {
                if (annoMat(j,i)) {
                    hsqAnnoObs[i]  += probAnnoObs[i]  * betajSq;
                    hsqAnnoNull[i] += snpAnnoCntInv[j] * betaSqNull;
                }
            }
        } else {
            for (unsigned i=0; i<numAnno; ++i) {
                if (annoMat(j,i)) {
                    hsqAnnoNull[i] += snpAnnoCntInv[j] * betaSqNull;
                }
            }
        }
    }
    
    values = hsqAnnoObs.array()/hsqAnnoNull.array();
    
//    cout << "Intercept " << values[0] << endl;

}

void ApproxBayesRC::AnnoJointPerSnpHsqEnrichment::compute(const VectorXf &snpEffects, const MatrixXf &annoMat, const AnnoCondProb &annoCondProb, const vector<AnnoInfo*> &annoInfoVec, const VectorXf &snpAnnoCntInv, const vector<SnpInfo*> *incdSnpInfoVec){
    unsigned numAnno = annoMat.cols();
    unsigned numSnps = annoMat.rows();
    VectorXf hsqAnnoObs;
    VectorXf hsqAnnoNull;
    hsqAnnoObs.setZero(numAnno);
    hsqAnnoNull.setZero(numAnno);
    
    float betaSqNull = snpEffects.squaredNorm()/float(numSnps);
    VectorXf probAnnoObs(numAnno);
    for (unsigned j=0; j<numSnps; ++j) {
        if (incdSnpInfoVec && (*incdSnpInfoVec)[j]->fitAsFixedEffect) continue;
        if (snpEffects[j]) {
            //cout << j << endl;
            probAnnoObs.setZero(numAnno);
            for (unsigned i=0; i<numAnno; ++i) {
                if (annoMat(j,i)) probAnnoObs[i] = annoCondProb[0]->values[i];
            }
            probAnnoObs /= probAnnoObs.sum();
//            cout << j << " annoCondProb " << annoCondProb[0]->values.head(5).transpose() << endl;
//            cout << j << " annoMat " << annoMat.row(j).head(5) << endl;
//            cout << j << " probAnno " << probAnno.head(5).transpose() << endl;
            float betajSq = snpEffects[j]*snpEffects[j];
            for (unsigned i=0; i<numAnno; ++i) {
                if (annoMat(j,i)) {
                    hsqAnnoObs[i]  += probAnnoObs[i] * betajSq;
                    hsqAnnoNull[i] += snpAnnoCntInv[j] * betaSqNull;
                }
            }
        } else {
            for (unsigned i=0; i<numAnno; ++i) {
                if (annoMat(j,i)) {
                    hsqAnnoNull[i] += snpAnnoCntInv[j] * betaSqNull;
                }
            }
        }
    }
    
    values = hsqAnnoObs.array()/hsqAnnoNull.array();
    
//    cout << "annoCondProb " << annoCondProb[0]->values.head(5).transpose() << endl;
//    cout << "hsqAnnoObs " << hsqAnnoObs.head(5).transpose() << endl;
//    cout << "hsqAnnoNull " << hsqAnnoNull.head(5).transpose() << endl;
//    cout << "Intercept " << values[0] << " Coding_UCSC " << values[1] << endl;

}

//void ApproxBayesRC::AnnoJointPerSnpHsqEnrichment::compute(const VectorXf &snpEffects, const MatrixXf &annoMat, const AnnoEffects &alpha, const vector<AnnoInfo*> &annoInfoVec, const MatrixXf &snpPi){
//    // partition the SNP effect into annotation-specific components based on SBayesRC model
//    unsigned numAnno = annoMat.cols();
//    unsigned numSnps = annoMat.rows();
//    VectorXf hsqAnnoVec;
//    hsqAnnoVec.setZero(size);
//    
//    float cdfMu = Normal::cdf_01(alpha[0]->values[0]);
//    float pdfMu = Normal::pdf_01(alpha[0]->values[0]);
//    VectorXf weights = pdfMu * alpha[0]->values;
//    weights[0] = cdfMu;
//    
//    VectorXf betajc;
//    for (unsigned j=0; j<numSnps; ++j) {
//        if (snpEffects[j]) {
//            float p2 = annoMat.row(j).dot(weights);
//            betajc = weights/p2 * snpEffects[j];
////            cout << j << " " << snpEffects[j] << " " << betajc[0] << " " << p2 << " " << 1.0-snpPi(j,0) << endl;
//            for (unsigned i=0; i<numAnno; ++i) {
//                if (0==i) {
//                    hsqAnnoVec[i] += betajc[0]*betajc[0];
//                } else {
//                    float tmp = betajc[0] + annoMat(j,i)*betajc[i];
//                    hsqAnnoVec[i] += tmp*tmp;
//                }
//            }
//        }
//    }
//    
//    float gwPerSnpHsq = snpEffects.squaredNorm()/(float)numSnps;
//    for (unsigned i=0; i<numAnno; ++i) {
//        AnnoInfo *anno = annoInfoVec[i];
//        if (anno->isBinary) {
//            values[i] = hsqAnnoVec[i]/(float)anno->size /gwPerSnpHsq;
//        } else {
//            values[i] = hsqAnnoVec[i]/anno->ssq /gwPerSnpHsq;
//        }
//        if(0==i) cout << i << " " << values[i] << endl;
//    }
//    
//    if (values[0] > 20) {
//        for (unsigned j=0; j<numSnps; ++j) {
//            if (snpEffects[j]) {
//                float p2 = annoMat.row(j).dot(weights);
//                betajc = weights/p2 * snpEffects[j];
//                cout << j << " " << snpEffects[j] << " " << betajc[0] << " " << p2 << " " << 1.0-snpPi(j,0) << endl;
//            }
//        }
//    }
//}

//void ApproxBayesRC::AnnoJointPerSnpHsqEnrichment::compute(const VectorXf &snpEffects, const MatrixXf &annoMat, const AnnoEffects &alpha, const vector<AnnoInfo*> &annoInfoVec, const MatrixXf &snpPi){
//    // partition the SNP effect into annotation-specific components based on SBayesRC model
//    unsigned numAnno = annoMat.cols();
//    unsigned numSnps = annoMat.rows();
//    VectorXf hsqAnnoVec;
//    hsqAnnoVec.setZero(size);
//    
//    float CDFmu = Normal::cdf_01(alpha[0]->values[0]);
//    float PDFmu = Normal::pdf_01(alpha[0]->values[0]);
//    float logCDFmu = log(CDFmu);
//    VectorXf weights = PDFmu/CDFmu * alpha[0]->values;  // Taylor expansion of logCDF(mu+delta)
//    weights[0] = logCDFmu;
//    
//    for (unsigned j=0; j<numSnps; ++j) {
//        if (snpEffects[j]) {
//            float pj2 = expf(annoMat.row(j).dot(weights));
//            float wj0 = 2.0f*weights[0] + logf(snpEffects[j]*snpEffects[j]/(pj2*pj2));
////            cout << j << " " << snpEffects[j] << " " << betajc[0] << " " << pj2 << " " << 1.0-snpPi(j,0) << endl;
//            for (unsigned i=0; i<numAnno; ++i) {
//                if (0==i) {
//                    hsqAnnoVec[i] += expf(wj0);
//                } else {
//                    float wjc = 2.0f*weights[i];
//                    hsqAnnoVec[i] += expf(wj0 + annoMat(j,i)*wjc);
//                }
//            }
//        }
//    }
//    
//    float gwPerSnpHsq = snpEffects.squaredNorm()/(float)numSnps;
//    for (unsigned i=0; i<numAnno; ++i) {
//        AnnoInfo *anno = annoInfoVec[i];
//        if (anno->isBinary) {
//            values[i] = hsqAnnoVec[i]/(float)anno->size /gwPerSnpHsq;
//        } else {
//            values[i] = hsqAnnoVec[i]/anno->ssq /gwPerSnpHsq;
//        }
//        if(0==i) cout << i << " " << values[i] << endl;
//    }
//    
//    if (values[0] > 20) {
//        for (unsigned j=0; j<numSnps; ++j) {
//            if (snpEffects[j]) {
//                float pj2 = expf(annoMat.row(j).dot(weights));
//                float wj0 = 2.0f*weights[0] + logf(snpEffects[j]*snpEffects[j]/(pj2*pj2));
//                cout << j << " " << snpEffects[j]*snpEffects[j] << " " << expf(wj0) << " " << pj2 << " " << 1.0-snpPi(j,0) << endl;
//            }
//        }
//    }
//}
//
//void ApproxBayesRC::AnnoJointPerSnpRsqEnrichment::compute(const VectorXf &snpEffectMeans, const MatrixXf &annoMat, const AnnoEffects &alpha, const vector<AnnoInfo*> &annoInfoVec){
//    // partition the SNP effect into annotation-specific components based on SBayesRC model
//    unsigned numAnno = annoMat.cols();
//    unsigned numSnps = annoMat.rows();
//    VectorXf hsqAnnoVec;
//    hsqAnnoVec.setZero(size);
//    
//    float cdfMu = Normal::cdf_01(alpha[0]->values[0]);
//    float pdfMu = Normal::pdf_01(alpha[0]->values[0]);
//    VectorXf weights = pdfMu * alpha[0]->values;
//    weights[0] = cdfMu;
//    
//    VectorXf betajc;
//    for (unsigned j=0; j<numSnps; ++j) {
//        float p2 = annoMat.row(j).dot(weights);
//        betajc = weights/p2 * snpEffectMeans[j];
//        for (unsigned i=0; i<numAnno; ++i) {
//            if (0==i) {
//                hsqAnnoVec[i] += betajc[0]*betajc[0];
//            } else {
//                float tmp = betajc[0] + annoMat(j,i)*betajc[i];
//                hsqAnnoVec[i] += tmp*tmp;
//            }
//        }
//    }
//    
//    float gwPerSnpHsq = snpEffectMeans.squaredNorm()/(float)numSnps;
//    for (unsigned i=0; i<numAnno; ++i) {
//        AnnoInfo *anno = annoInfoVec[i];
//        if (anno->isBinary) {
//            values[i] = hsqAnnoVec[i]/(float)anno->size /gwPerSnpHsq;
//        } else {
//            values[i] = hsqAnnoVec[i]/anno->ssq /gwPerSnpHsq;
//        }
//    }
//
//}

void ApproxBayesRC::AnnoJointPerSnpRsqEnrichment::compute(const VectorXf &snpEffectMeans, const MatrixXf &annoMat, const AnnoJointProb &annoJointProb, const vector<AnnoInfo*> &annoInfoVec, const VectorXf &gamma, const VectorXf &snpAnnoCntInv){
    unsigned numAnno = annoMat.cols();
    unsigned numSnps = annoMat.rows();
    unsigned numDist = annoJointProb.numDist;
    VectorXf rsqAnnoObs;
    VectorXf rsqAnnoNull;
    rsqAnnoObs.setZero(numAnno);
    rsqAnnoNull.setZero(numAnno);
    
    float betaSqNull = snpEffectMeans.squaredNorm()/float(numSnps);
    VectorXf probAnnoObs(numAnno);
    for (unsigned j=0; j<numSnps; ++j) {
        probAnnoObs.setZero(numAnno);
        for (unsigned i=0; i<numAnno; ++i) {
            if (annoMat(j,i)) {
                for (unsigned k=1; k<numDist; ++k) {  // loop over nonzero distributions
                    probAnnoObs[i]  += annoJointProb[k]->values[i] * gamma[k];
                }
            }
        }
        probAnnoObs  /= probAnnoObs.sum();
        float betajSq = snpEffectMeans[j]*snpEffectMeans[j];
        for (unsigned i=0; i<numAnno; ++i) {
            if (annoMat(j,i)) {
                rsqAnnoObs[i]  += probAnnoObs[i]  * betajSq;
                rsqAnnoNull[i] += snpAnnoCntInv[j] * betaSqNull;
            }
        }
    }
    
    values = rsqAnnoObs.array()/rsqAnnoNull.array();
    
//    cout << "Intercept " << values[0] << endl;

}

void ApproxBayesRC::AnnoJointPerSnpRsqEnrichment::compute(const VectorXf &snpEffectMeans, const MatrixXf &annoMat, const AnnoCondProb &annoCondProb, const vector<AnnoInfo*> &annoInfoVec, const VectorXf &snpAnnoCntInv){
    unsigned numAnno = annoMat.cols();
    unsigned numSnps = annoMat.rows();
    VectorXf rsqAnnoObs;
    VectorXf rsqAnnoNull;
    rsqAnnoObs.setZero(numAnno);
    rsqAnnoNull.setZero(numAnno);
    
    float betaSqNull = snpEffectMeans.squaredNorm()/float(numSnps);
    VectorXf probAnnoObs(numAnno);
    for (unsigned j=0; j<numSnps; ++j) {
        probAnnoObs.setZero(numAnno);
        for (unsigned i=0; i<numAnno; ++i) {
            if (annoMat(j,i)) probAnnoObs[i] = annoCondProb[0]->values[i];
        }
        probAnnoObs /= probAnnoObs.sum();
        float betajSq = snpEffectMeans[j]*snpEffectMeans[j];
        for (unsigned i=0; i<numAnno; ++i) {
            if (annoMat(j,i)) {
                rsqAnnoObs[i]  += probAnnoObs[i] * betajSq;
                rsqAnnoNull[i] += snpAnnoCntInv[j] * betaSqNull;
            }
        }
    }
    
    values = rsqAnnoObs.array()/rsqAnnoNull.array();
}


//void ApproxBayesRC::computeSnpVarg(const MatrixXf &annoMat, const VectorXf &annoPerSnpHsqEnrich, const float varg, const unsigned numSnps){
//    VectorXf tau = annoPerSnpHsqEnrich.array() - 1.0;
//    for (unsigned i=0; i<numSnps; ++i) {
//        snpVarg[i] = (1.0 + annoMat.row(i).dot(tau))*varg;
//    }
//}

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

void ApproxBayesRC::getSnpAnnoCntInv(const MatrixXf &annoMat, VectorXf &snpAnnoCntInv){
    // get number of nonzero annotations for each SNP
    unsigned numAnno = annoMat.cols();
    unsigned numSnps = annoMat.rows();
    VectorXf snpAnnoCnt;
    snpAnnoCnt.setZero(numSnps);
    for (unsigned j=0; j<numSnps; ++j) {
        for (unsigned i=0; i<numAnno; ++i) {
            if (annoMat(j,i)) ++snpAnnoCnt[j];
        }
    }
    snpAnnoCntInv = snpAnnoCnt.array().inverse();
}

void ApproxBayesRC::sampleUnknowns(const unsigned iter){
    if (lowRankModel) {
        fixedSnpEffects.sampleFromFC_eigen(snpEffects, wcorrBlocks, data.Qblocks, whatBlocks, data.keptLdBlockInfoVec, data.nGWASblock, vareBlk.values, data.numFixedEffectSnps);
        snpEffects.sampleFromFC_eigen(wcorrBlocks, data.Qblocks, whatBlocks,
                                      data.keptLdBlockInfoVec, data.nGWASblock, vareBlk.values,
                                      snpPi, gamma.values, vargRandom.value, deltaPi, hsqPercModel, sigmaSq.value, data.incdSnpInfoVec);
    } else {
        fixedSnpEffects.sampleFromFC_sparse(snpEffects, rcorr, data.ZPZsp, data.ZPZdiag, data.snp2pq, data.incdSnpInfoVec, varg.value, vare.value, data.numFixedEffectSnps);
        snpEffects.sampleFromFC_sparse(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.chromInfoVec, sigmaSq.value, snpPi, gamma.values, vare.value,
                                       vargRandom.value, hsqPercModel, deltaPi, data.incdSnpInfoVec);
    }
    
    if (algorithm == tgs) {  // To improve mixing, apply tempered Gibbs sampling on high-LD SNP
        snpEffects.sampleFromTGS_eigen(wcorrBlocks, data.Qblocks, whatBlocks,
                                       data.LDmap, data.keptLdBlockInfoVec, data.nGWASblock, vareBlk.values,
                                       snpPi, gamma.values, vargRandom.value, deltaPi, hsqPercModel, sigmaSq.value);
    } else if (algorithm == tgs_thin) {
        if (!(iter % 10)) snpEffects.sampleFromTGS_eigen(wcorrBlocks, data.Qblocks, whatBlocks,
                                                         data.LDmap, data.keptLdBlockInfoVec, data.nGWASblock, vareBlk.values,
                                                         snpPi, gamma.values, vargRandom.value, deltaPi, hsqPercModel, sigmaSq.value);
    }
    
    snpEffects.computePosteriorMean(iter);
    snpPip.getValues(snpEffects.pip);
    nnzSnp.getValue(snpEffects.numNonZeros);
    numSnps.getValues(snpEffects.numSnpMix);

    if (robustMode) {
        sigmaSq.value = vargRandom.value/(data.numIncdSnps*gamma.values.dot(Pis.values));  // LDpred2's parameterisation
    } else {
        sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
    }
    
    if (estimatePi) {
        //if (algorithm == gibbs) {
            annoEffects.sampleFromFC_Gibbs(snpEffects.z, data.annoMat, sigmaSqAnno.values, snpP,
                data.numFixedEffectSnps ? &data.incdSnpInfoVec : nullptr,
                data.numFixedEffectSnps ? &baselineSnpP : nullptr);
            annoCondProb.compute_probit(annoEffects, data.annoInfoVec);
        //} else {
        //    annoEffects.sampleFromFC_MH(snpEffects.z, data.annoMat, sigmaSqAnno.values, snpP, ...);
        //    annoCondProb.compute_logistic(annoEffects, data.annoInfoVec);
        //}
        sigmaSqAnno.sampleFromFC(annoEffects.ssq);
        computePiFromP(snpP, snpPi);
        annoJointProb.compute(annoCondProb);
    }
            
    if (lowRankModel) {
        vargBlk.compute(whatBlocks);
        vareBlk.sampleFromFC(wcorrBlocks, vargBlk.values, snpEffects.ssqBlocks, data.nGWASblock, data.numEigenvalBlock);
        //vareBlk.sampleFromFC(wcorrBlocks, snpEffects.values, data.b, data.nGWASblock, data.keptLdBlockInfoVec);
        vargRandom.value = vargBlk.total;
        vare.value = vareBlk.mean;
    } else {
        varg.compute(snpEffects.values, data.ZPy, rcorr);
        vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, status);
        if (data.numFixedEffectSnps && data.varGenotypic > 0.f)
            vargRandom.value = varg.value * (data.vargRandomSnpEffects / data.varGenotypic);
        else
            vargRandom.value = varg.value;
    }
    hsqRandom.compute(vargRandom.value, data.varPhenotypic);
    if (data.numFixedEffectSnps) {
        hsqFixed.compute(snpEffects.values, data);
        hsq.value = hsqRandom.value + hsqFixed.value;
        if (lowRankModel)
            varg.value = vargRandom.value + hsqFixed.value * data.varPhenotypic;
    } else {
        hsq.value = hsqRandom.value;
        if (lowRankModel)
            varg.value = vargRandom.value;
    }
    
    annoGenVar.compute(snpEffects.values, snpEffects.snpset, data.annoMat, data.numFixedEffectSnps ? &data.incdSnpInfoVec : nullptr);
    annoTotalGenVar.compute(annoGenVar);
    annoPerSnpHsqEnrich.compute(annoTotalGenVar.values, data.annoInfoVec);
    annoJointPerSnpHsqEnrich.compute(snpEffects.values, data.annoMat, annoJointProb, data.annoInfoVec, gamma.values, snpAnnoCntInv, data.numFixedEffectSnps ? &data.incdSnpInfoVec : nullptr);
    if (estimateRsqEnrich) {
        annoPerSnpRsqEnrich.compute(snpEffects.fcMean, data.annoMat, data.annoInfoVec);
        annoJointPerSnpRsqEnrich.compute(snpEffects.fcMean, data.annoMat, annoJointProb, data.annoInfoVec, gamma.values, snpAnnoCntInv);
    }

    Vgs.compute(snpEffects.values, snpEffects.snpset);
    
    snpHsqPep.compute(snpEffects.values, vargRandom.value);

    if (!(iter % 10)) {
        if (lowRankModel) {
            nBadSnps.compute_eigen(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, wcorrBlocks, data.Qblocks, data.keptLdBlockInfoVec, iter, data);
        } else {
            nBadSnps.compute_sparse(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, rcorr, data.ZPZsp, data.chromInfoVec, iter, data);
        }
    }

    if (!(iter % 100)) {
        if (lowRankModel) {
            rounding.computeWcorr_eigen(data.wcorrBlocks, data.Qblocks, data.keptLdBlockInfoVec, snpEffects.values, wcorrBlocks);
        } else {
            rounding.computeRcorr_sparse(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
        }
    }
}

//void ApproxBayesRC::sampleUnknownsTGS(vector<vector<int> > &selectedSnps){
//    snpEffects.sampleFromTGS_eigen(selectedSnps, wcorrBlocks, data.Qblocks, whatBlocks,
//                                  data.keptLdBlockInfoVec, data.nGWASblock, vareBlk.values,
//                                  snpPi, gamma.values, varg.value, deltaPi, hsqPercModel, sigmaSq.value);
//}


void BayesRC::SnpEffects::sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes, const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare, VectorXf &ghat, const MatrixXf &snpPi, const float varg, const bool hsqPercModel, DeltaPi &deltaPi){
    sumSq = 0.0;
    wtdSumSq = 0.0;
    numNonZeros = 0;
    
    ghat.setZero(ycorr.size());

    pip.setZero(size);
    z.setZero(size, ndist-1);   // indicator variables for conditional membership
    
    // R specific parameters
    ArrayXf wtdSigmaSq(ndist);
    ArrayXf invWtdSigmaSq(ndist);
    ArrayXf logWtdSigmaSq(ndist);
    MatrixXf logPi = snpPi.array().log().matrix();
        
    if (hsqPercModel && varg) {
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

    if (shuffle) Gadget::shuffle_vector(snpIndexVec);

    unsigned i;
    for (unsigned t = 0; t < size; t++) {
        i = snpIndexVec[t];
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
        pip[i] = 1.0f - probDelta[0];

        delta = bernoulli.sample(probDelta);
        
        snpset[delta].push_back(i);
        numSnpMix[delta]++;
        
        if (delta) {
            values[i] = normal.sample(uhat[delta], invLhs[delta]);
            ycorr += Z.col(i) * (oldSample - values[i]);
            if (weightedRes) ghat += Z.col(i).cwiseProduct(Rsqrt) * values[i];
            else ghat  += Z.col(i) * values[i];
            sumSq += values[i] * values[i];
            wtdSumSq += (values[i] * values[i]) / gamma[delta];
            ++numNonZeros;
            for(unsigned k2 = 0; k2 < delta ; k2++){
                z(i, k2) = 1;
            }
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
}

void BayesRC::initSnpPandPi(const VectorXf &pis, const unsigned numSnps, MatrixXf &snpP, MatrixXf &snpPi) {
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

void BayesRC::getSnpAnnoCntInv(const MatrixXf &annoMat, VectorXf &snpAnnoCntInv){
    // get number of nonzero annotations for each SNP
    unsigned numAnno = annoMat.cols();
    unsigned numSnps = annoMat.rows();
    VectorXf snpAnnoCnt;
    snpAnnoCnt.setZero(numSnps);
    for (unsigned j=0; j<numSnps; ++j) {
        for (unsigned i=0; i<numAnno; ++i) {
            if (annoMat(j,i)) ++snpAnnoCnt[j];
        }
    }
    snpAnnoCntInv = snpAnnoCnt.array().inverse();
}

void BayesRC::sampleUnknowns(const unsigned iter){
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }

    snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, data.Rsqrt, data.weightedRes, sigmaSq.value, Pis.values, gamma.values, vare.value, ghat, snpPi, varg.value, hsqPercModel, deltaPi);
    snpPip.getValues(snpEffects.pip);
    nnzSnp.getValue(snpEffects.numNonZeros);
    numSnps.getValues(snpEffects.numSnpMix);

    sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
    
    if (estimatePi) {
        //if (algorithm == gibbs) {
            annoEffects.sampleFromFC_Gibbs(snpEffects.z, data.annoMat, sigmaSqAnno.values, snpP,
                data.numFixedEffectSnps ? &data.incdSnpInfoVec : nullptr,
                data.numFixedEffectSnps ? &baselineSnpP : nullptr);
            annoCondProb.compute_probit(annoEffects, data.annoInfoVec);
        //} else {
        //    annoEffects.sampleFromFC_MH(snpEffects.z, data.annoMat, sigmaSqAnno.values, snpP, ...);
        //    annoCondProb.compute_logistic(annoEffects, data.annoInfoVec);
        //}
        sigmaSqAnno.sampleFromFC(annoEffects.ssq);
        computePiFromP(snpP, snpPi);
        annoJointProb.compute(annoCondProb);
    }
    
    varg.compute(ghat);
    vare.sampleFromFC(ycorr);
    hsq.compute(varg.value, vare.value);
    
    annoGenVar.compute(snpEffects.values, snpEffects.snpset, data.annoMat, data.numFixedEffectSnps ? &data.incdSnpInfoVec : nullptr);
    annoTotalGenVar.compute(annoGenVar);
    annoPerSnpHsqEnrich.compute(annoTotalGenVar.values, data.annoInfoVec);
//    annoJointPerSnpHsqEnrich.compute(annoJointProb, data.annoInfoVec, gamma.values, varg.value, hsqPercModel, sigmaSq.value);
    annoJointPerSnpHsqEnrich.compute(snpEffects.values, data.annoMat, annoJointProb, data.annoInfoVec, gamma.values, snpAnnoCntInv, data.numFixedEffectSnps ? &data.incdSnpInfoVec : nullptr);

    Vgs.compute(snpEffects.values, data.Z, snpEffects.snpset, varg.value);
    
    if (!(iter % 100)) rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);
}


void ApproxBayesRD::AnnoEffects::initIntercept(const VectorXf &pis){
    VectorXf p(numComp);
    unsigned ndist = pis.size();
    for (unsigned i=1; i<ndist; ++i) {
        p[i-1] = pis.tail(ndist-i).sum();
        if (i>1) p[i-1] /= pis.tail(ndist-i+1).sum();
    }
    for (unsigned i = 0; i<numComp; ++i) {
        (*this)[i]->values[0] = Normal::quantile_01(p[i]);
    }
}

void ApproxBayesRD::AnnoEffects::sampleFromFC_indep(MatrixXf &z, const MatrixXf &annoMat, const VectorXf &sigmaSq, const float pi, MatrixXf &snpP){
    VectorXf numOnes(numComp);
    #pragma omp parallel for schedule(dynamic)
    for (unsigned i=0; i<numComp; ++i) {
        numOnes[i] = z.col(i).sum();
    }
        
    unsigned numSnps = z.rows();
    numNonZeros.setZero(numComp);
    
    float logPi = log(pi);
    float logPiComp = log(1.0-pi);
    float logDelta0, logDelta1, probDelta1;

    for (unsigned i=0; i<numComp; ++i) {
        VectorXf &alphai = (*this)[i]->values;
        VectorXf &pip = (*this)[i]->pip;
        Stat::Bernoulli &bernoulli = (*this)[i]->bernoulli;
        
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
            if (i==0) {
                annoDiagi = annoDiag;
            } else {
                annoDiagi[0] = numOnes[i-1];
                #pragma omp parallel for schedule(dynamic)
                for (unsigned k=1; k<numAnno; ++k) {
                    annoDiagi[k] = annoMati.col(k).squaredNorm();
                }
            }

            // compute the mean of truncated normal distribution
            VectorXf mean = annoMati * alphai;

            // sample latent variables
            for (unsigned j=0; j<numDP; ++j) {
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
            float logSigmaSq = log(sigmaSq[i]);
            alphai[0] = Normal::sample(ahat, invLhs);
            y.array() += oldSample - alphai[0];
            pip[0] = 1.0;

            // annotations are fitted with a BayesC-type mixture prior
            // shuffle the annotations
            vector<int> shuffled_index = Gadget::shuffle_index(1, numAnno-1);

            ssq[i] = 0;
            for (unsigned t=0; t<shuffled_index.size(); ++t) {
                unsigned k = shuffled_index[t];
                oldSample = alphai[k];
                rhs = annoMati.col(k).dot(y) + annoDiagi[k]*oldSample;
                invLhs = 1.0/(annoDiagi[k] + 1.0/sigmaSq[i]);
                ahat = invLhs*rhs;
                logDelta1 = 0.5*(logf(invLhs) - logSigmaSq + ahat*rhs) + logPi;
                logDelta0 = logPiComp;
                probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
                pip[k] = probDelta1;
                
                if (bernoulli.sample(probDelta1)) {
                    alphai[k] = Normal::sample(ahat, invLhs);
                    y += annoMati.col(k) * (oldSample - alphai[k]);
                    ssq[i] += alphai[k] * alphai[k];
                    ++numNonZeros[i];
                    //            cout << i << " " << k << " " << alphai[k] << " " << ahat << " " << invLhs << " " << annoDiagi[k] << " " << sigmaSq[i] << endl;
                } else {
                    if (oldSample) y += annoMati.col(k) * oldSample;
                    alphai[k] = 0.0;
                }
            }
        }
        //cout << i << " " << alphai.transpose() << endl;
        
        values.col(i) = alphai;
        
        #pragma omp parallel for schedule(dynamic)
        for (unsigned j=0; j<numSnps; ++j) {
            snpP(j,i) = Normal::cdf_01(annoMat.row(j).dot(alphai));
        }
    }
    
    nnz = numNonZeros.sum();
    numAnnoTotal = (numAnno - 1) * numComp;
    
    pip.setZero(numAnno);
    for (unsigned i=0; i<numComp; ++i) {
        pip += (*this)[i]->pip;
    }
    pip /= float(numComp);
}

//void ApproxBayesRD::AnnoEffects::sampleFromFC_joint(MatrixXf &z, const MatrixXf &annoMat, const VectorXf &sigmaSq, const float pi, MatrixXf &snpP){
//    VectorXf numOnes(numComp);
//#pragma omp parallel for schedule(dynamic)
//    for (unsigned i=0; i<numComp; ++i) {
//        numOnes[i] = z.col(i).sum();
//    }
//    
//    unsigned numSnps = z.rows();
//    numNonZeros.setZero(numComp);
//    
//    pip.setZero(numAnno);
//    pip[0] = 1.0;  // always fit intercept in the model
//    nnz = 0;
//    
//    VectorXi numDP(numComp);  // number of data points for each component
//    unsigned numNZComp = 0;   // number of components with at least one nonzero SNP
//    for (unsigned i=0; i<numComp; ++i) {
//        VectorXf &alphai = (*this)[i]->values;
//        if (i==0) numDP[i] = numSnps;
//        else numDP[i] = numOnes[i-1];
//        if(numDP[i] == 0){
//            alphai.setZero();
//            alphai[0] = -10.0;
//            ssq[i] = 0;
//        }else{
//            ++numNZComp;
//        }
//    }
//    
//    vector<VectorXf> APY(numNZComp);
//    vector<MatrixXf> APA(numNZComp);
//    MatrixXf Alpha(numAnno, numNZComp);
//    MatrixXf APAdiag(numAnno, numNZComp);
//    
//    Stat::Bernoulli &bernoulli = (*this)[0]->bernoulli;
//    
//    for (unsigned i=0; i<numNZComp; ++i) {
//        VectorXf &alphai = (*this)[i]->values;
//        VectorXf y, zi;
//        y.setZero(numDP[i]);
//        zi.setZero(numDP[i]);
//        const MatrixXf *annotMatP;
//        MatrixXf annoMatPO;
//        // get annotation coefficient matrix for component i
//        if (i==0) {
//            annotMatP = &annoMat;
//            zi = z.col(i);
//        } else {
//            annoMatPO.setZero(numDP[i], numAnno);
//            for (unsigned j=0, idx=0; j<numSnps; ++j) {
//                if (z(j,i-1)) {
//                    annoMatPO.row(idx) = annoMat.row(j);
//                    zi[idx] = z(j,i);
//                    ++idx;
//                }
//            }
//            annotMatP = &annoMatPO;
//        }
//        const MatrixXf &annoMati = (*annotMatP);
//        
//        
//        VectorXf annoDiagi(numAnno);
//        if (i==0) {
//            annoDiagi = annoDiag;
//        } else {
//            annoDiagi[0] = numOnes[i-1];
//#pragma omp parallel for schedule(dynamic)
//            for (unsigned k=1; k<numAnno; ++k) {
//                annoDiagi[k] = annoMati.col(k).squaredNorm();
//            }
//        }
//        
//        // compute the mean of truncated normal distribution
//        VectorXf mean = annoMati * alphai;
//        
//        // sample latent variables
//        for (unsigned j=0; j<numDP[i]; ++j) {
//            if (zi[j]) y[j] = TruncatedNormal::sample_lower_truncated(mean[j], 1.0, 0.0);
//            else y[j] = TruncatedNormal::sample_upper_truncated(mean[j], 1.0, 0.0);
//        }
//        
//        // adjust the latent variable by all annotation effects;
//        y -= mean;
//        
//        // intercept is fitted with a flat prior
//        float oldSample = alphai[0];
//        float rhs = y.sum() + annoDiagi[0]*oldSample;
//        float invLhs = 1.0/annoDiagi[0];
//        float ahat = invLhs*rhs;
//        float logSigmaSq = log(sigmaSq[i]);
//        alphai[0] = Normal::sample(ahat, invLhs);
//        y.array() += oldSample - alphai[0];
//        
//        APY[i] = annoMati.transpose()*y;
//        APA[i] = annoMati.transpose()*annoMati;
//        APAdiag.col(i) = annoDiagi;
//        
//        ssq[i] = 0;
//        Alpha.col(i) = alphai;
//    }
//    
//    // annotations are fitted with a BayesC-type mixture prior
//    float logPi = log(pi);
//    float logPiComp = log(1.0-pi);
//    VectorXf invSigmaSq = sigmaSq.array().inverse();
//    VectorXf logSigmaSq = sigmaSq.array().log();
//    invSigmaSq.conservativeResize(numNZComp);
//    logSigmaSq.conservativeResize(numNZComp);
//    
//    // shuffle the annotations
//    vector<int> shuffled_index = Gadget::shuffle_index(1, numAnno-1);
//    for (unsigned t=0; t<shuffled_index.size(); ++t) {
//        unsigned k = shuffled_index[t];
//        VectorXf oldSample = Alpha.row(k);
//        VectorXf rhs(numNZComp);
//        for (unsigned i=0; i<numNZComp; ++i) {
//            rhs[i] = APY[i][k] + APAdiag(k,i)*oldSample[i];
//        }
//        VectorXf invLhs = (APAdiag.row(k).transpose() + invSigmaSq).array().inverse();
//        VectorXf ahat = invLhs.cwiseProduct(rhs);
//        VectorXf logDelta1Vec = 0.5*(invLhs.array().log() - logSigmaSq.array() + ahat.array()*rhs.array());
//        float logDelta1 = logDelta1Vec.sum() + logPi;
//        float logDelta0 = logPiComp;
//        float probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
//        pip[k] = probDelta1;
////        cout << k << " probDelta1 " << probDelta1 << " logDelta1 " << logDelta1 << " logDelta1Vec " << logDelta1Vec.transpose() << " logPi " << logPi << " logPiComp " << logPiComp << endl;
//        
//        if (bernoulli.sample(probDelta1)) {
//            for (unsigned i=0; i<numNZComp; ++i) {
//                Alpha(k,i) = Normal::sample(ahat[i], invLhs[i]);
//                APY[i] += APA[i].col(k) * (oldSample[i] - Alpha(k,i));
//                ssq[i] += Alpha(k,i) * Alpha(k,i);
//            }
//            ++nnz;
//        } else {
//            for (unsigned i=0; i<numNZComp; ++i) {
//                if (oldSample[i]) APY[i] += APA[i].col(k) * oldSample[i];
//                Alpha(k,i) = 0.0;
//            }
//        }
//    }
//    
//    //    cout << "Alpha\n" << Alpha << endl;
//    
//    for (unsigned i=0; i<numNZComp; ++i) {
//        (*this)[i]->values = Alpha.col(i);
//    }
//        
//#pragma omp parallel for schedule(dynamic)
//    for (unsigned i=0; i<numComp; ++i) {
//        VectorXf &alphai = (*this)[i]->values;
//        values.col(i) = alphai;
//        for (unsigned j=0; j<numSnps; ++j) {
//            float val = annoMat.row(j).dot(alphai);
//            snpP(j,i) = Normal::cdf_01(annoMat.row(j).dot(alphai));
//        }
//    }
//    
//    numAnnoTotal = numAnno - 1;
//}

void ApproxBayesRD::AnnoEffects::sampleFromFC_joint(MatrixXf &z, const MatrixXf &annoMat, const VectorXf &sigmaSq, const float pi, MatrixXf &snpP){
    VectorXf numOnes(numComp);
#pragma omp parallel for schedule(dynamic)
    for (unsigned i=0; i<numComp; ++i) {
        numOnes[i] = z.col(i).sum();
    }
    
    unsigned numSnps = z.rows();
    numNonZeros.setZero(numComp);
    
    pip.setZero(numAnno);
    pip[0] = 1.0;  // always fit intercept in the model
    nnz = 0;
    
    VectorXi numDP(numComp);  // number of data points for each component
    unsigned numNZComp = 0;   // number of components with at least one nonzero SNP
    for (unsigned i=0; i<numComp; ++i) {
        VectorXf &alphai = (*this)[i]->values;
        if (i==0) numDP[i] = numSnps;
        else numDP[i] = numOnes[i-1];
        if(numDP[i] == 0){
            alphai.setZero();
            alphai[0] = -10.0;
            ssq[i] = 0;
        }else{
            ++numNZComp;
        }
    }
    
    vector<VectorXf> Y(numNZComp);
    vector<MatrixXf> Amat(numNZComp);
    MatrixXf Alpha(numAnno, numNZComp);
    MatrixXf APAdiag(numAnno, numNZComp);
    
    Stat::Bernoulli &bernoulli = (*this)[0]->bernoulli;
    
    for (unsigned i=0; i<numNZComp; ++i) {
        VectorXf &alphai = (*this)[i]->values;
        VectorXf y, zi;
        y.setZero(numDP[i]);
        zi.setZero(numDP[i]);
        const MatrixXf *annotMatP;
        MatrixXf annoMatPO;
        // get annotation coefficient matrix for component i
        if (i==0) {
            annotMatP = &annoMat;
            zi = z.col(i);
        } else {
            annoMatPO.setZero(numDP[i], numAnno);
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
        if (i==0) {
            annoDiagi = annoDiag;
        } else {
            annoDiagi[0] = numOnes[i-1];
#pragma omp parallel for schedule(dynamic)
            for (unsigned k=1; k<numAnno; ++k) {
                annoDiagi[k] = annoMati.col(k).squaredNorm();
            }
        }
        
        // compute the mean of truncated normal distribution
        VectorXf mean = annoMati * alphai;
        
        // sample latent variables
        for (unsigned j=0; j<numDP[i]; ++j) {
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
        float logSigmaSq = log(sigmaSq[i]);
        alphai[0] = Normal::sample(ahat, invLhs);
        y.array() += oldSample - alphai[0];
        
        Y[i] = y;
        Amat[i] = annoMati;
        APAdiag.col(i) = annoDiagi;
        
        ssq[i] = 0;
        Alpha.col(i) = alphai;
    }
    
    // annotations are fitted with a BayesC-type mixture prior
    float logPi = log(pi);
    float logPiComp = log(1.0-pi);
    VectorXf invSigmaSq = sigmaSq.array().inverse();
    VectorXf logSigmaSq = sigmaSq.array().log();
    invSigmaSq.conservativeResize(numNZComp);
    logSigmaSq.conservativeResize(numNZComp);
    
    // shuffle the annotations
    vector<int> shuffled_index = Gadget::shuffle_index(1, numAnno-1);
    for (unsigned t=0; t<shuffled_index.size(); ++t) {
        unsigned k = shuffled_index[t];
        VectorXf oldSample = Alpha.row(k);
        VectorXf rhs(numNZComp);
        for (unsigned i=0; i<numNZComp; ++i) {
            rhs[i] = Amat[i].col(k).dot(Y[i]) + APAdiag(k,i)*oldSample[i];
        }
        VectorXf invLhs = (APAdiag.row(k).transpose() + invSigmaSq).array().inverse();
        VectorXf ahat = invLhs.cwiseProduct(rhs);
        VectorXf logDelta1Vec = 0.5*(invLhs.array().log() - logSigmaSq.array() + ahat.array()*rhs.array());
        float logDelta1 = logDelta1Vec.sum() + logPi;
        float logDelta0 = logPiComp;
        float probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
        pip[k] = probDelta1;
//        cout << k << " probDelta1 " << probDelta1 << " logDelta1 " << logDelta1 << " logDelta1Vec " << logDelta1Vec.transpose() << " logPi " << logPi << " logPiComp " << logPiComp << endl;
        
        if (bernoulli.sample(probDelta1)) {
            for (unsigned i=0; i<numNZComp; ++i) {
                Alpha(k,i) = Normal::sample(ahat[i], invLhs[i]);
                Y[i] += Amat[i].col(k) * (oldSample[i] - Alpha(k,i));
                ssq[i] += Alpha(k,i) * Alpha(k,i);
            }
            ++nnz;
        } else {
            for (unsigned i=0; i<numNZComp; ++i) {
                if (oldSample[i]) Y[i] += Amat[i].col(k) * oldSample[i];
                Alpha(k,i) = 0.0;
            }
        }
    }
    
    //    cout << "Alpha\n" << Alpha << endl;
    
    for (unsigned i=0; i<numNZComp; ++i) {
        (*this)[i]->values = Alpha.col(i);
    }
        
#pragma omp parallel for schedule(dynamic)
    for (unsigned i=0; i<numComp; ++i) {
        VectorXf &alphai = (*this)[i]->values;
        values.col(i) = alphai;
        for (unsigned j=0; j<numSnps; ++j) {
            float val = annoMat.row(j).dot(alphai);
            snpP(j,i) = Normal::cdf_01(annoMat.row(j).dot(alphai));
        }
    }
    
    numAnnoTotal = numAnno - 1;
}

void ApproxBayesRD::AnnoCondProb::compute(const AnnoEffects &annoEffects, const vector<AnnoInfo*> &annoInfoVec){
    for (unsigned i=0; i<annoEffects.numComp; ++i) {
        for (unsigned j=0; j<annoEffects.numAnno; ++j) {
            VectorXf &alpha = annoEffects[i]->values;
            if (j==0) (*this)[i]->values[j] = Normal::cdf_01(alpha[j]);
            else (*this)[i]->values[j] = Normal::cdf_01(alpha[0] + annoInfoVec[j]->sd*alpha[j]);  // NEW
        }
    }
}

void ApproxBayesRD::sampleUnknowns(const unsigned iter){
    if (lowRankModel) {
        fixedSnpEffects.sampleFromFC_eigen(snpEffects, wcorrBlocks, data.Qblocks, whatBlocks, data.keptLdBlockInfoVec, data.nGWASblock, vareBlk.values, data.numFixedEffectSnps);
        snpEffects.sampleFromFC_eigen(wcorrBlocks, data.Qblocks, whatBlocks,
                                data.keptLdBlockInfoVec, data.nGWASblock, vareBlk.values,
                                snpPi, gamma.values, vargRandom.value, deltaPi, hsqPercModel, sigmaSq.value, data.incdSnpInfoVec);
    } else {
        fixedSnpEffects.sampleFromFC_sparse(snpEffects, rcorr, data.ZPZsp, data.ZPZdiag, data.snp2pq, data.incdSnpInfoVec, varg.value, vare.value, data.numFixedEffectSnps);
        snpEffects.sampleFromFC_sparse(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.chromInfoVec, sigmaSq.value, snpPi, gamma.values, vare.value,
                                vargRandom.value, hsqPercModel, deltaPi, data.incdSnpInfoVec);
    }
    snpEffects.computePosteriorMean(iter);
    snpPip.getValues(snpEffects.pip);
    nnzSnp.getValue(snpEffects.numNonZeros);
    numSnps.getValues(snpEffects.numSnpMix);

    if (robustMode) {
        sigmaSq.value = vargRandom.value/(data.numIncdSnps*gamma.values.dot(Pis.values));
    } else {
        sigmaSq.sampleFromFC(snpEffects.wtdSumSq, snpEffects.numNonZeros);
    }
    
    if (estimatePi) {
        //annoEffects.sampleFromFC_indep(snpEffects.z, data.annoMat, sigmaSqAnno.values, piAnno.value, snpP);
        annoEffects.sampleFromFC_joint(snpEffects.z, data.annoMat, sigmaSqAnno.values, piAnno.value, snpP);
        sigmaSqAnno.sampleFromFC(annoEffects.ssq);
        piAnno.sampleFromFC(annoEffects.numAnnoTotal, annoEffects.nnz);
        computePiFromP(snpP, snpPi);
        annoCondProb.compute(annoEffects, data.annoInfoVec);
        annoJointProb.compute(annoCondProb);
        annoPip.getValues(annoEffects.pip);
    }
            
    if (lowRankModel) {
        vargBlk.compute(whatBlocks);
        vareBlk.sampleFromFC(wcorrBlocks, vargBlk.values, snpEffects.ssqBlocks, data.nGWASblock, data.numEigenvalBlock);
        //vareBlk.sampleFromFC(wcorrBlocks, snpEffects.values, data.b, data.nGWASblock, data.keptLdBlockInfoVec);
        vargRandom.value = vargBlk.total;
        vare.value = vareBlk.mean;
    } else {
        varg.compute(snpEffects.values, data.ZPy, rcorr);
        vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, status);
        if (data.numFixedEffectSnps && data.varGenotypic > 0.f)
            vargRandom.value = varg.value * (data.vargRandomSnpEffects / data.varGenotypic);
        else
            vargRandom.value = varg.value;
    }
    hsqRandom.compute(vargRandom.value, data.varPhenotypic);
    if (data.numFixedEffectSnps) {
        hsqFixed.compute(snpEffects.values, data);
        hsq.value = hsqRandom.value + hsqFixed.value;
        if (lowRankModel)
            varg.value = vargRandom.value + hsqFixed.value * data.varPhenotypic;
    } else {
        hsq.value = hsqRandom.value;
        if (lowRankModel)
            varg.value = vargRandom.value;
    }
    
    annoGenVar.compute(snpEffects.values, snpEffects.snpset, data.annoMat, data.numFixedEffectSnps ? &data.incdSnpInfoVec : nullptr);
    annoTotalGenVar.compute(annoGenVar);
    annoPerSnpHsqEnrich.compute(annoTotalGenVar.values, data.annoInfoVec);
    if (estimateRsqEnrich) {
        annoJointPerSnpHsqEnrich.compute(annoJointProb, data.annoInfoVec, gamma.values, vargRandom.value, hsqPercModel, sigmaSq.value);
    }
    
    Vgs.compute(snpEffects.values, snpEffects.snpset);
    
    snpHsqPep.compute(snpEffects.values, vargRandom.value);

    if (!(iter % 10)) {
        if (lowRankModel) {
            nBadSnps.compute_eigen(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, wcorrBlocks, data.Qblocks, data.keptLdBlockInfoVec, iter, data);
        } else {
            nBadSnps.compute_sparse(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, rcorr, data.ZPZsp, data.chromInfoVec, iter, data);
        }
    }

    if (!(iter % 100)) {
        if (lowRankModel) {
            rounding.computeWcorr_eigen(data.wcorrBlocks, data.Qblocks, data.keptLdBlockInfoVec, snpEffects.values, wcorrBlocks);
        } else {
            rounding.computeRcorr_sparse(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
        }
    }
}


// ApproxBayesAPP: Bivariate Bayesian Analysis with Annotation and Pleiotropy
// ============================================================================

// Helper function to sample from Inverse-Wishart distribution (2x2)
// X ~ InvWishart(df, scale) means: X^-1 ~ Wishart(df, scale^-1)
static Matrix2f sampleInverseWishart2x2(float df, const Matrix2f &scale) {
    // Sample from Wishart(df, scale^-1) then invert
    Matrix2f scaleInv = scale.inverse();
    
    // Sample from Wishart(df, scaleInv) using Bartlett decomposition
    // For 2x2 Wishart: W = L * L' where L is lower triangular
    Matrix2f L;
    L.setZero();
    
    // L[0,0] ~ sqrt(Chi^2(df)) where Chi^2(df) = Gamma(df/2, 2)
    Stat::Gamma gammaSampler;
    float chi2_0 = gammaSampler.sample(0.5f * df, 2.0f);
    L(0,0) = sqrtf(chi2_0);
    
    // L[1,0] ~ N(0, 1) (standard normal)
    L(1,0) = Stat::snorm();
    
    // L[1,1] ~ sqrt(Chi^2(df-1))
    float chi2_1 = gammaSampler.sample(0.5f * (df - 1.0f), 2.0f);
    L(1,1) = sqrtf(chi2_1);
    
    // Scale by Cholesky of scaleInv
    LLT<Matrix2f> chol(scaleInv);
    if (chol.info() != Success) {
        // If Cholesky fails, use identity
        Matrix2f W = L * L.transpose();
        return W.inverse();
    }
    Matrix2f L_scale = chol.matrixL();
    L = L_scale * L;
    
    Matrix2f W = L * L.transpose();
    return W.inverse();
}


void ApproxBayesAPP::AnnoGenotypicVarBivariate::accumulate(const VarBlocks &varBlocks) {
    const unsigned nBlocks = varBlocks.totalVargBlk.size();
    for (unsigned cat = 0; cat < numCategories; ++cat) {
        G_vec[cat].setZero();
        for (unsigned b = 0; b < nBlocks; ++b) {
            G_vec[cat](0,0) += varBlocks.vargBlkCat[b](0, cat);
            G_vec[cat](1,1) += varBlocks.vargBlkCat[b](1, cat);
            G_vec[cat](0,1) += varBlocks.vargCovBlkCat[b][cat];
        }
        G_vec[cat](1,0) = G_vec[cat](0,1);
    }
}

unsigned ApproxBayesAPP::DeltaBivariate::countNonZero(const unsigned numCategories) const {
    unsigned nnz = 0;
    for (unsigned cat = 0; cat < numCategories; ++cat) {
        for (unsigned i = 0; i < size; ++i) {
            unsigned markerIndex = cat * size + i;
            if (deltaMatrix(markerIndex, 0) > 0.5f || deltaMatrix(markerIndex, 1) > 0.5f)
                ++nnz;
        }
    }
    return nnz;
}

// Map [delta1, delta2] to Pi index matching Julia order: 0=[1,1], 1=[1,0], 2=[0,1], 3=[0,0]
static inline unsigned deltaToPiIdx(float d0, float d1) {
    return (d0 < 0.5f ? 2u : 0u) + (d1 < 0.5f ? 1u : 0u);
}

void ApproxBayesAPP::DeltaBivariate::countLociByState(const unsigned numCategories) {
    for (unsigned cat = 0; cat < numCategories; ++cat) {
        nLociCounts[cat].setZero();
        for (unsigned i = 0; i < size; ++i) {
            unsigned markerIndex = cat * size + i;
            unsigned idx = deltaToPiIdx(deltaMatrix(markerIndex, 0), deltaMatrix(markerIndex, 1));
            nLociCounts[cat][idx] += 1.0f;
        }
    }
}

void ApproxBayesAPP::VarBlocks::resize(const unsigned nBlocks, const unsigned nCategories) {
    vargBlkCat.resize(nBlocks);
    ssqBlkCat.resize(nBlocks);
    vargCovBlkCat.resize(nBlocks);
    totalVargBlk.resize(nBlocks);
    for (unsigned b = 0; b < nBlocks; ++b) {
        vargBlkCat[b]    = MatrixXf::Zero(2, nCategories);
        ssqBlkCat[b]     = VectorXf::Zero(3 * nCategories);
        vargCovBlkCat[b] = VectorXf::Zero(nCategories);
        totalVargBlk[b]  = Matrix2f::Zero();
    }
}

void ApproxBayesAPP::VarBlocks::setZero(const unsigned nBlocks) {
    for (unsigned b = 0; b < nBlocks; ++b) {
        vargBlkCat[b].setZero();
        ssqBlkCat[b].setZero();
        vargCovBlkCat[b].setZero();
        totalVargBlk[b].setZero();
    }
}

Matrix2f ApproxBayesAPP::VarBlocks::computeTotal() const {
    Matrix2f total = Matrix2f::Zero();
    for (const auto &m : totalVargBlk) total += m;
    return total;
}

void ApproxBayesAPP::VarBlocks::compute(const MatrixXf &alphaMatrix, const unsigned nSNPs,
                                        const vector<LDBlockInfo*> &keptLdBlockInfoVec,
                                        const vector<MatrixXf> &Qblocks, const unsigned numCategories) {
    unsigned nBlocks = keptLdBlockInfoVec.size();
    setZero(nBlocks);
    for (unsigned b = 0; b < nBlocks; ++b) {
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[b];
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd = blockInfo->endSnpIdx;
        Ref<const MatrixXf> Q = Qblocks[b];

        VectorXf what_array_total1 = VectorXf::Zero(Q.rows());
        VectorXf what_array_total2 = VectorXf::Zero(Q.rows());

        for (unsigned cat = 0; cat < numCategories; ++cat) {
            VectorXf alphaArray_c1 = VectorXf::Zero(blockEnd - blockStart + 1);
            VectorXf alphaArray_c2 = VectorXf::Zero(blockEnd - blockStart + 1);
            for (unsigned i = blockStart; i <= blockEnd; ++i) {
                unsigned markerIndex = cat * nSNPs + i;
                alphaArray_c1[i - blockStart] = alphaMatrix(markerIndex, 0);
                alphaArray_c2[i - blockStart] = alphaMatrix(markerIndex, 1);
            }

            VectorXf what_array_c1 = Q * alphaArray_c1;
            VectorXf what_array_c2 = Q * alphaArray_c2;

            vargBlkCat[b](0, cat) = what_array_c1.dot(what_array_c1);
            vargBlkCat[b](1, cat) = what_array_c2.dot(what_array_c2);
            ssqBlkCat[b][0 * numCategories + cat] = alphaArray_c1.dot(alphaArray_c1);
            ssqBlkCat[b][1 * numCategories + cat] = alphaArray_c2.dot(alphaArray_c2);
            vargCovBlkCat[b][cat] = what_array_c1.dot(what_array_c2);
            ssqBlkCat[b][2 * numCategories + cat] = alphaArray_c1.dot(alphaArray_c2);

            what_array_total1 += what_array_c1;
            what_array_total2 += what_array_c2;
        }

        totalVargBlk[b](0,0) = what_array_total1.dot(what_array_total1);
        totalVargBlk[b](1,1) = what_array_total2.dot(what_array_total2);
        totalVargBlk[b](0,1) = totalVargBlk[b](1,0) = what_array_total1.dot(what_array_total2);
    }
}

void ApproxBayesAPP::SnpPIPBivariate::compute(const MatrixXf &deltaMatrix, const unsigned numCategories) {
    unsigned nSNPs = pip1().size;
    VectorXf p1 = VectorXf::Zero(nSNPs);
    VectorXf p2 = VectorXf::Zero(nSNPs);
    for (unsigned cat = 0; cat < numCategories; ++cat) {
        for (unsigned i = 0; i < nSNPs; ++i) {
            unsigned markerIndex = cat * nSNPs + i;
            if (deltaMatrix(markerIndex, 0) > 0.5f) p1[i] = 1.0f;
            if (deltaMatrix(markerIndex, 1) > 0.5f) p2[i] = 1.0f;
        }
    }
    pip1().getValues(p1);
    pip2().getValues(p2);
}

void ApproxBayesAPP::AnnoPiBivariate::sampleFromFC(const vector<VectorXf> &nLociCounts) {
    // Sample from Dirichlet distribution for each category
    for (unsigned c = 0; c < numCategories; ++c) {
        VectorXf counts = nLociCounts[c];
        VectorXf alpha = counts.array() + alphaVec.array();
        piVec[c] = Stat::Dirichlet().sample(4, alpha);
    }
}

void ApproxBayesAPP::PiBivariate::compute(const AnnoPiBivariate &piAnno, const VectorXf &nLociAnno) {
    VectorXf margPi = VectorXf::Zero(4);
    float totalLoci = nLociAnno.sum();
    if (totalLoci > 0.0f) {
        for (unsigned c = 0; c < piAnno.numCategories; ++c)
            margPi += piAnno.piVec[c] * nLociAnno[c];
        margPi /= totalLoci;
    }
    // margPi order: 0=[1,1], 1=[1,0], 2=[0,1], 3=[0,0]. Output: Pi_00, Pi_10, Pi_01, Pi_11
    (*this)[0]->value = margPi[3];  // Pi_00 = P([0,0])
    (*this)[1]->value = margPi[1];  // Pi_10 = P([1,0])
    (*this)[2]->value = margPi[2];  // Pi_01 = P([0,1])
    (*this)[3]->value = margPi[0];  // Pi_11 = P([1,1])
}

void ApproxBayesAPP::HsqBivariate::compute(const Matrix2f &vargTotal, const vector<Matrix2f> &R_blk) {
    float vare1 = 0.0f, vare2 = 0.0f, cove = 0.0f;
    for (const auto &R : R_blk) { vare1 += R(0,0); vare2 += R(1,1); cove += R(0,1); }
    vare1 /= R_blk.size();
    vare2 /= R_blk.size();
    cove /= R_blk.size();
    (*this)[0]->value = vargTotal(0, 0);  // varg is 1 for both traits
    (*this)[1]->value = vargTotal(1, 1);
    // Co-heritability = Covg_12 / Covp_12
    (*this)[2]->value = vargTotal(0, 1) / cove;
}

void ApproxBayesAPP::GeneticCorrelation::compute(const Matrix2f &vargTotal) {
    float vg1 = vargTotal(0, 0);
    float vg2 = vargTotal(1, 1);
    float cvg = vargTotal(0, 1);
    float denom = sqrtf(vg1 * vg2);
    (*this)[0]->value = (denom > 0.0f) ? (cvg / denom) : 0.0f;
}

void ApproxBayesAPP::AnnoVarEffectsBivariate::sampleFromFC(const MatrixXf &betaMatrix, const VectorXf &nLociAnno) {
    // Compute SSE per category from betaMatrix
    unsigned nSNPs = betaMatrix.rows() / numCategories;
    vector<Matrix2f> SSE_vec(numCategories);
    for (unsigned c = 0; c < numCategories; ++c) {
        SSE_vec[c].setZero();
        for (unsigned i = 0; i < nSNPs; ++i) {
            unsigned idx = c * nSNPs + i;
            float b1 = betaMatrix(idx, 0), b2 = betaMatrix(idx, 1);
            SSE_vec[c](0,0) += b1 * b1;
            SSE_vec[c](1,1) += b2 * b2;
            SSE_vec[c](0,1) += b1 * b2;
        }
        SSE_vec[c](1,0) = SSE_vec[c](0,1);
    }

    // Sample from Inverse-Wishart for each category
    const float df_G = 6.0f;  // df = 4 + nTraits

    for (unsigned c = 0; c < numCategories; ++c) {
        Matrix2f scale_G = A_vec[c] * (df_G - 3.0f);  // df_G - nTraits - 1
        Matrix2f SSE = SSE_vec[c];
        Matrix2f scale_tilde = scale_G + SSE;
        
        // Ensure scale_tilde is positive definite
        LLT<Matrix2f> chol(scale_tilde);
        if (chol.info() != Success) {
            // If not positive definite, add small value to diagonal
            scale_tilde(0,0) += 0.01f;
            scale_tilde(1,1) += 0.01f;
        }
        
        A_vec[c] = sampleInverseWishart2x2(df_G + nLociAnno[c], scale_tilde);
        Ainv_vec[c] = A_vec[c].inverse();
    }
}

void ApproxBayesAPP::ResidualVarBivariate::normalise(const vector<Vector2f> &nGWASblocks) {
    for (unsigned b = 0; b < R_blk.size(); ++b) {
        Matrix2f Rnorm = R_blk[b];
        const Vector2f &n = nGWASblocks[b];
        for (unsigned i = 0; i < 2; ++i)
            for (unsigned j = i; j < 2; ++j) {
                Rnorm(i,j) /= sqrtf(n[i] * n[j]);
                Rnorm(j,i)  = Rnorm(i,j);
            }
        Rinv_blk[b] = Rnorm.inverse();
    }
}

void ApproxBayesAPP::ResidualVarBivariate::sampleFromFC(vector<VectorXf> &wcorrBlocks, const vector<Vector2f> &nGWASblocks,
                                                    const VectorXf &numEigenvalBlock, const VarBlocks &varBlocks) {
    const vector<MatrixXf> &vargBlkCat    = varBlocks.vargBlkCat;
    const vector<VectorXf> &ssqBlkCat     = varBlocks.ssqBlkCat;
    const vector<Matrix2f> &totalVargBlk  = varBlocks.totalVargBlk;
    // Sample residual variance per block, matching Julia implementation
    unsigned nBlocks = wcorrBlocks.size();
    
    for (unsigned b = 0; b < nBlocks; ++b) {
        unsigned nEigen = wcorrBlocks[b].size() / 2;  // wcorr is [w1; w2]
        VectorXf w1 = wcorrBlocks[b].segment(0, nEigen);
        VectorXf w2 = wcorrBlocks[b].segment(nEigen, nEigen);
        
        // Get sample sizes per trait (Julia: nInd = [nInd1, nInd2])
        float n1 = nGWASblocks[b][0];
        float n2 = nGWASblocks[b][1];
        
        // Compute SSE as in Julia: dot(ycorri, ycorrj) * sqrt(nInd[i] * nInd[j])
        Matrix2f SSE;
        SSE(0,0) = w1.dot(w1) * sqrtf(n1 * n1);
        SSE(1,1) = w2.dot(w2) * sqrtf(n2 * n2);
        SSE(0,1) = SSE(1,0) = w1.dot(w2) * sqrtf(n1 * n2);
        
        // Sample R (Julia: sample_variance_sumstats)
        Matrix2f scale_tilde = scale_R + SSE;
        
        // Ensure positive definite
        LLT<Matrix2f> chol(scale_tilde);
        if (chol.info() != Success) {
            scale_tilde(0,0) += 0.01f;
            scale_tilde(1,1) += 0.01f;
        }
        
        Matrix2f sampled_R = sampleInverseWishart2x2(df_R + numEigenvalBlock[b], scale_tilde);
        
        // R tuning logic from Julia (lines 639-650)
        float Rcor = sampled_R(0,1) / sqrtf(sampled_R(0,0) * sampled_R(1,1));  // correlation
        
        // Check threshold for each trait (Julia lines 639-648)
        for (unsigned traiti = 0; traiti < 2; ++traiti) {
            float totalVarg = totalVargBlk[b](traiti, traiti);
            if (totalVarg <= 0.0f) {
                R_blk[b](traiti, traiti) = 1.0f;
                continue;
            }
            // Compute thres = sum(ssq_blk_cat[b][traiti, :]) / totalvarg_blk[b][traiti, traiti]
            float ssq_sum = 0.0f;
            for (unsigned c = 0; c < vargBlkCat[b].cols(); ++c) {
                ssq_sum += ssqBlkCat[b][traiti * vargBlkCat[b].cols() + c];
            }
            float thres = ssq_sum / totalVarg;
            
            if (thres > 1.1f) {
                R_blk[b](traiti, traiti) = sampled_R(traiti, traiti);
            } else {
                R_blk[b](traiti, traiti) = 1.0f;
            }
        }
        
        // Tune covariance in R_blk
        float Rcov = Rcor * sqrtf(R_blk[b](0,0) * R_blk[b](1,1));
        R_blk[b](0,1) = R_blk[b](1,0) = Rcov;
    }
    
    // Average across blocks for next iteration (as in Julia code line 820-822)
    Matrix2f R_blk_sum = Matrix2f::Zero();
    for (unsigned b = 0; b < nBlocks; ++b) {
        R_blk_sum += R_blk[b];
    }
    Matrix2f R_blkmean = R_blk_sum / nBlocks;
    
    // Update all blocks with mean
    for (unsigned b = 0; b < nBlocks; ++b) {
        R_blk[b] = R_blkmean;
    }
    normalise(nGWASblocks);
}

void ApproxBayesAPP::AnnoGenotypicVarBivariate::compute(const vector<VectorXf> &whatBlocks, const vector<MatrixXf> &Qblocks,
                                                const vector<LDBlockInfo*> &keptLdBlockInfoVec,
                                                const MatrixXf &alphaMatrix, const MatrixXf &annoMatrix) {
    // Compute genetic variance per category (alphaMatrix is (nSNPs*nCategory) x 2)
    unsigned nBlocks = keptLdBlockInfoVec.size();
    unsigned nSNPs = alphaMatrix.rows() / numCategories;
    
    // Reset
    for (unsigned c = 0; c < numCategories; ++c) {
        G_vec[c].setZero();
    }
    
    for (unsigned b = 0; b < nBlocks; ++b) {
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[b];
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd = blockInfo->endSnpIdx;
        
        Ref<const MatrixXf> Q = Qblocks[b];
        
        for (unsigned c = 0; c < numCategories; ++c) {
            VectorXf what1_c = VectorXf::Zero(Q.rows());
            VectorXf what2_c = VectorXf::Zero(Q.rows());
            
            for (unsigned i = blockStart; i <= blockEnd; ++i) {
                if (annoMatrix(i, c) != 0.0) {  // SNP in this category
                    int qcol = eigenQColIndex(blockInfo, i);
                    if (qcol < 0) continue;
                    unsigned markerIndex = c * nSNPs + i;
                    Ref<const VectorXf> Qi = Q.col(qcol);
                    what1_c += Qi * alphaMatrix(markerIndex, 0);
                    what2_c += Qi * alphaMatrix(markerIndex, 1);
                }
            }
            
            G_vec[c](0,0) += what1_c.dot(what1_c);
            G_vec[c](1,1) += what2_c.dot(what2_c);
            G_vec[c](0,1) = G_vec[c](1,0) += what1_c.dot(what2_c);
        }
    }
}

void ApproxBayesAPP::SnpEffectsBivariate::sampleFromFC(MatrixXf &deltaMatrix, vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks,
                                         vector<VectorXf> &whatBlocks, VectorXf &vareBlocks,
                                         const vector<Matrix2f> &A_vec, const AnnoPiBivariate &Pi,
                                         const VectorXf &snp2pq, const vector<Matrix2f> &Rinv_blk,
                                         const vector<LDBlockInfo*> &keptLdBlockInfoVec,
                                         const vector<Vector2f> &nGWASblocks, const MatrixXf &annoMat) {

    unsigned nBlocks = keptLdBlockInfoVec.size();
    unsigned nCon = 0;  // TODO: Support continuous annotations

    // Initialise whatBlocks (persist across iterations; only resize when needed)
    if (whatBlocks.size() != nBlocks) {
        whatBlocks.resize(nBlocks);
        for (unsigned b = 0; b < nBlocks; ++b)
            whatBlocks[b].setZero(wcorrBlocks[b].size());
    }
    for (unsigned b = 0; b < nBlocks; ++b)
        whatBlocks[b].setZero();
    
    // Pre-compute xpx and xArray for each block and category (matching Julia lines 248-265)
    vector<vector<VectorXf> > xpx_dict(nBlocks);  // blk -> vector of xpx per category
    vector<vector<MatrixXf> > xArray_dict(nBlocks);  // blk -> vector of xArray per category
    
    for (unsigned blk = 0; blk < nBlocks; ++blk) {
        Ref<const MatrixXf> Q = Qblocks[blk];
        Ref<VectorXf> wcorr = wcorrBlocks[blk];
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[blk];
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd = blockInfo->endSnpIdx;
        unsigned nMarkerb = blockEnd - blockStart + 1;
        
        if (blockInfo->eigenColRemap.size() > 0 && (unsigned)Q.cols() != nMarkerb) {
            throw("Error: ApproxBayesAPP does not support sub-LD eigen (many SNPs --skip); use SBayesC/R/RC or omit --skip for this analysis.");
        }
        
        // Verify wcorr size matches expected bivariate size
        if (wcorr.size() != 2 * Q.rows()) {
            throw("Error in ApproxBayesAPP::SnpEffects::sampleFromFC: Block " + to_string(blk) + 
                  " - wcorr.size()=" + to_string(wcorr.size()) + " but expected " + to_string(2*Q.rows()) + 
                  " (2 * Q.rows()). Q.rows()=" + to_string(Q.rows()));
        }
        
        xpx_dict[blk].resize(nCon + 1);
        xArray_dict[blk].resize(nCon + 1);
        
        // For continuous annotations (if any)
        for (unsigned c = 0; c < nCon; ++c) {
            xpx_dict[blk][c] = VectorXf(nMarkerb);
            for (unsigned i = 0; i < nMarkerb; ++i) {
                unsigned snpIdx = blockStart + i;
                float annoVal = annoMat(snpIdx, c);
                Ref<const VectorXf> Qcol = Q.col(i);
                xpx_dict[blk][c][i] = annoVal * annoVal * Qcol.dot(Qcol);
            }
            // xArray = Q * diag(annoMat[:, c])
            xArray_dict[blk][c] = Q;
            for (unsigned i = 0; i < nMarkerb; ++i) {
                unsigned snpIdx = blockStart + i;
                xArray_dict[blk][c].col(i) *= annoMat(snpIdx, c);
            }
        }
        
        // For categorical annotations (last element)
        xpx_dict[blk][nCon] = VectorXf(nMarkerb);
        
        if (Q.cols() != nMarkerb) {
            throw("Error: Q.cols()=" + to_string(Q.cols()) + " != nMarkerb=" + to_string(nMarkerb) + " in block " + to_string(blk));
        }
        
        for (unsigned i = 0; i < nMarkerb; ++i) {
            Ref<const VectorXf> Qcol = Q.col(i);
            xpx_dict[blk][nCon][i] = Qcol.dot(Qcol);
        }
        xArray_dict[blk][nCon] = Q;
    }
    
    // Process each block (matching Julia lines 491-652)
    for (unsigned blk = 0; blk < nBlocks; ++blk) {
        if (blk >= Qblocks.size() || blk >= wcorrBlocks.size() || blk >= whatBlocks.size()) {
            throw("Error: blk=" + to_string(blk) + " out of bounds. Qblocks.size()=" + to_string(Qblocks.size()) + ", wcorrBlocks.size()=" + to_string(wcorrBlocks.size()) + ", whatBlocks.size()=" + to_string(whatBlocks.size()));
        }
        
        // Check keptLdBlockInfoVec first
        if (blk >= keptLdBlockInfoVec.size()) {
            throw("Error: blk=" + to_string(blk) + " >= keptLdBlockInfoVec.size()=" + to_string(keptLdBlockInfoVec.size()));
        }
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[blk];
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd = blockInfo->endSnpIdx;
        unsigned nMarkerb = blockEnd - blockStart + 1;
        
        // Create Ref objects after we have nMarkerb
        Ref<const MatrixXf> Q = Qblocks[blk];
        Ref<VectorXf> wcorr = wcorrBlocks[blk];
        Ref<VectorXf> what = whatBlocks[blk];
        
        // Get sample sizes per trait (Julia: nInd = [nInd1, nInd2])
        if (blk >= nGWASblocks.size()) {
            throw("Error: blk=" + to_string(blk) + " >= nGWASblocks.size()=" + to_string(nGWASblocks.size()));
        }
        Vector2f nInd;
        // Access Vector2f elements directly to avoid potential Ref issues
        Vector2f nGWASblk = nGWASblocks[blk];  // Make a copy
        nInd[0] = nGWASblk[0];
        nInd[1] = nGWASblk[1];
        
        // Get pre-normalised Rinv for this block (computed by ResidualVarBivariate::sampleFromFC)
        Matrix2f Rinv = Rinv_blk[blk];
        
        // Shuffle marker order (Julia line 514)
        vector<unsigned> markerOrder;
        for (unsigned i = 0; i < nMarkerb; ++i) {
            markerOrder.push_back(i);
        }
        random_shuffle(markerOrder.begin(), markerOrder.end());
        
        // Initialize xArrayc/xpxc (Julia lines 502-503)
        MatrixXf xArrayc = xArray_dict[blk][0];
        VectorXf xpxc = xpx_dict[blk][0];
        
        // Loop through markers (Julia lines 515-599)
        for (unsigned marker_idx = 0; marker_idx < markerOrder.size(); ++marker_idx) {
            unsigned marker = markerOrder[marker_idx];
            unsigned true_marker_num = blockStart + marker;  // marker position across all SNPs
            
            // Get annotation categories for this marker (Julia line 517)
            VectorXf annoindexm = annoMat.row(true_marker_num).transpose();
            
            // Loop through categories (Julia lines 519-598)
            for (unsigned cat = 0; cat < numCategories; ++cat) {
                // Update xArrayc/xpxc for continuous groups (Julia lines 520-524)
                if (nCon > 0) {
                    if (cat <= nCon) {
                        xArrayc = xArray_dict[blk][cat];
                        xpxc = xpx_dict[blk][cat];
                    } else {
                        // For categorical groups after nCon, use the last element
                        xArrayc = xArray_dict[blk][nCon];
                        xpxc = xpx_dict[blk][nCon];
                    }
                } else {
                    // No continuous annotations, use categorical (index nCon = 0)
                    if (xArray_dict[blk].size() == 0 || xArray_dict[blk][0].size() == 0) {
                        throw("Error: xArray_dict not properly initialized for block " + to_string(blk));
                    }
                    xArrayc = xArray_dict[blk][0];
                    xpxc = xpx_dict[blk][0];
                }
                
                Matrix2f Ginv = A_vec[cat].inverse();
                VectorXf piVec = Pi.piVec[cat];
                
                // Check if marker belongs to this category (Julia line 529)
                if (annoindexm[cat] != 0.0) {
                    // markerIndex = (cat - 1) * my_nsnp + true_marker_num (Julia line 530)
                    unsigned markerIndex = cat * this->betaTotal[0].size + true_marker_num;
                    
                    VectorXf x = xArrayc.col(marker);
                    float xpx_marker = xpxc[marker];
                    
                    // Get current values (Julia lines 532-536)
                    Vector2f beta, oldAlpha, delta;
                    beta[0] = betaMatrix(markerIndex, 0);
                    beta[1] = betaMatrix(markerIndex, 1);
                    delta[0] = deltaMatrix(markerIndex, 0);
                    delta[1] = deltaMatrix(markerIndex, 1);
                    oldAlpha[0] = alphaMatrix(markerIndex, 0);
                    oldAlpha[1] = alphaMatrix(markerIndex, 1);
                    
                    // Compute w (Julia line 536)
                    // For bivariate, wcorr should be size 2 * numEigenvalBlock[blk]
                    unsigned nEigen = Q.rows();  // Number of eigenvalues for this block
                    
                    // Check sizes before dot product
                    if (wcorr.size() != 2 * nEigen) {
                        throw("Error: wcorr size mismatch in ApproxBayesAPP block " + to_string(blk) + ". Expected " + to_string(2*nEigen) + ", got " + to_string(wcorr.size()) + ". Q.rows()=" + to_string(Q.rows()));
                    }
                    if (x.size() != nEigen) {
                        throw("Error: x size mismatch in ApproxBayesAPP block " + to_string(blk) + " marker " + to_string(marker) + ". Expected " + to_string(nEigen) + ", got " + to_string(x.size()) + ". xArrayc.rows()=" + to_string(xArrayc.rows()) + ", xArrayc.cols()=" + to_string(xArrayc.cols()));
                    }
                    
                    Vector2f w;
                    VectorXf wcorr1 = wcorr.segment(0, nEigen);
                    VectorXf wcorr2 = wcorr.segment(nEigen, nEigen);
                    
                    // Final size check before dot product
                    if (x.size() != wcorr1.size()) {
                        throw("Error: Final size check failed - x.size()=" + to_string(x.size()) + ", wcorr1.size()=" + to_string(wcorr1.size()) + ", nEigen=" + to_string(nEigen));
                    }
                    
                    w[0] = x.dot(wcorr1) + xpx_marker * oldAlpha[0];
                    w[1] = x.dot(wcorr2) + xpx_marker * oldAlpha[1];
                    
                    // Sample delta and beta for each trait (Julia lines 539-578)
                    vector<unsigned> traitOrder = {0, 1};
                    random_shuffle(traitOrder.begin(), traitOrder.end());
                    
                    Vector2f newAlpha = oldAlpha;
                    
                    for (unsigned k_idx = 0; k_idx < 2; ++k_idx) {
                        unsigned k = traitOrder[k_idx];
                        unsigned nok = 1 - k;
                        
                        float Ginv11 = Ginv(k, k);
                        float Ginv12 = Ginv(k, nok);
                        
                        // When delta[k] = 0 (Julia lines 546-549)
                        float invLhs0 = 1.0f / Ginv11;
                        float rhs0 = -Ginv12 * beta[nok];
                        float gHat0 = rhs0 * invLhs0;
                        
                        // When delta[k] = 1 (Julia lines 550-553)
                        float C11 = Ginv11 + Rinv(k, k) * xpx_marker;
                        float C12 = Ginv12 + xpx_marker * delta[nok] * Rinv(k, nok);
                        float invLhs1 = 1.0f / C11;
                        Vector2f w_vec;
                        w_vec << w[0], w[1];
                        float rhs1 = w_vec.dot(Rinv.col(k)) - C12 * beta[nok];
                        float gHat1 = rhs1 * invLhs1;
                        
                        // Compute log probabilities (Julia lines 555-563)
                        Vector2f d0 = delta;
                        Vector2f d1 = delta;
                        d0[k] = 0.0;
                        d1[k] = 1.0;
                        
                        // Get Pi values; states: 0=[1,1], 1=[1,0], 2=[0,1], 3=[0,0] (matching Julia)
                        unsigned idx0 = deltaToPiIdx(d0[0], d0[1]);
                        unsigned idx1 = deltaToPiIdx(d1[0], d1[1]);
                        
                        float logDelta0 = -0.5f * (logf(Ginv11) - gHat0*gHat0*Ginv11) + logf(piVec[idx0]);
                        float logDelta1 = -0.5f * (logf(C11) - gHat1*gHat1*C11) + logf(piVec[idx1]);
                        float probDelta1 = 1.0f / (1.0f + expf(logDelta0 - logDelta1));
                        
                        // Sample delta[k] (Julia lines 565-577)
                        if (Stat::ranf() < probDelta1) {
                            delta[k] = 1.0f;
                            beta[k] = gHat1 + Stat::snorm() * sqrtf(invLhs1);
                            newAlpha[k] = beta[k];
                            wcorr.segment(k * nEigen, nEigen) += x * (oldAlpha[k] - newAlpha[k]);
                        } else {
                            beta[k] = gHat0 + Stat::snorm() * sqrtf(invLhs0);
                            delta[k] = 0.0f;
                            newAlpha[k] = 0.0f;
                            if (oldAlpha[k] != 0.0f) {
                                wcorr.segment(k * nEigen, nEigen) += x * oldAlpha[k];
                            }
                        }
                    }
                    
                    // Store final values (Julia lines 589-592)
                    deltaMatrix(markerIndex, 0) = delta[0];
                    deltaMatrix(markerIndex, 1) = delta[1];
                    betaMatrix(markerIndex, 0) = beta[0];
                    betaMatrix(markerIndex, 1) = beta[1];
                    alphaMatrix(markerIndex, 0) = newAlpha[0];
                    alphaMatrix(markerIndex, 1) = newAlpha[1];
                }
            }
        }
        
        // Update what for genetic variance computation (will be done separately)
        // This is computed in sampleUnknowns after all blocks are processed
    }
    updateValues();
}

void ApproxBayesAPP::sampleUnknowns(const unsigned iter) {
    // Main MCMC sampling loop matching Julia implementation
    
    // 1. Sample SNP effects (delta and beta) - matching Julia lines 461-599
    snpEffects.sampleFromFC(delta.deltaMatrix, wcorrBlocks, data.Qblocks, whatBlocks,
                            vareBlk.values, sigmaSq.A_vec, piAnno,
                            data.snp2pq, vare.Rinv_blk,
                            data.keptLdBlockInfoVec, nGWASblocks, annoMatIncd);
    
    // Compute per-SNP PIP from delta (marginalize over categories)
    snpPip.compute(delta.deltaMatrix, numCategories);
    
    // 2. Compute genetic variance per block per category (Julia lines 603-633)
    varBlocks.compute(snpEffects.alphaMatrix, snpEffects.betaTotal[0].size,
                      data.keptLdBlockInfoVec, data.Qblocks, numCategories);
    
    // 3. Sample residual variance (Julia lines 635-651, 820-822)
    if (estimateVare) {
        vare.sampleFromFC(wcorrBlocks, nGWASblocks, data.numEigenvalBlock, varBlocks);
    }
    
    // 4. Count loci in each state per category (Julia lines 580-587)
    delta.countLociByState(numCategories);

    // 5. Sample annotation-stratified Pi, then compute marginal Pi
    if (estimatePi) {
        piAnno.sampleFromFC(delta.nLociCounts);
    }
    // Compute nLociAnno (number of loci per annotation) for the marginal pi
    pi.compute(piAnno, nLociPerAnno);

    // 6 & 7. Sample marker effect variance (Julia lines 751-760, 789-808)
    if (estimateVara)
        sigmaSq.sampleFromFC(snpEffects.betaMatrix, nLociPerAnno);

    // 8. Accumulate genetic variance per category across blocks (Julia lines 726-747)
    varg.accumulate(varBlocks);

    // 9. Compute total genetic variance, heritability, and global genetic correlation
    vargTotal = varBlocks.computeTotal();
    hsq.compute(vargTotal, vare.R_blk);
    geneticCorrelation.compute(vargTotal);

    // 10. Update nnz (count non-zero SNPs across all categories)
    nnzSnp.getValue(delta.countNonZero(numCategories));
}

void ApproxBayesAPP::SnpEffectsBivariate::updateValues(void) {
    // Compute total SNP effect (betaTotal = sum over categories) for each trait
    // Structure: alphaMatrix is (nSNPs * nCategory) × 2, stores category-specific alpha = delta * beta
    // We sum across categories to get total effect per SNP, stored in betaTotal[0].values and betaTotal[1].values
    unsigned nSNPs = betaTotal[0].size;
    
    // Initialize values to zero
    betaTotal[0].values.setZero();
    betaTotal[1].values.setZero();
    
    // Sum across categories for each SNP
    for (unsigned cat = 0; cat < numCategories; ++cat) {
        for (unsigned i = 0; i < nSNPs; ++i) {
            unsigned markerIndex = cat * nSNPs + i;
            // Use alpha (delta * beta) for values, matching the actual effects
            betaTotal[0].values[i] += alphaMatrix(markerIndex, 0);
            betaTotal[1].values[i] += alphaMatrix(markerIndex, 1);
        }
    }
}

