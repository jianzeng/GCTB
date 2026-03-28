//
//  model.hpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#ifndef model_hpp
#define model_hpp

#include <iostream>
#include <math.h>
#include "stat.hpp"
#include "data.hpp"

using namespace std;


class Parameter {
    // base class for a single parameter
public:
    string label;
    float value;   // sampled value
    
    unsigned numChains;      // for multiple chains
    VectorXf perChainValue;  // for multiple chains

    Parameter(const string &label): label(label){
        value = 0.0;
        numChains = 1;
    }
};

class ParamSet {
    // base class for a set of parameters of same kind, e.g. fixed effects, snp effects ...
public:
    string label;
    const vector<string> &header;
    unsigned size;
    VectorXf values;
        
    unsigned numChains;      // for multiple chains
    MatrixXf perChainValues; // for multiple chains

    ParamSet(const string &label, const vector<string> &header)
    : label(label), header(header), size(int(header.size())){
        values.setZero(size);
        numChains = 1;
    }
};

class ParamVec : public vector<Parameter*> {
    // a vector of parameters, e.g., pi values in BayesR model
public:
    string label;
    unsigned size;  // number of parameters
    VectorXf values;
    
    ParamVec(const string &label, const VectorXf &values):
    label(label), size(values.size()), values(values){
        for (unsigned i=0; i<size; ++i) {
            this->push_back(new Parameter(label + to_string(static_cast<long long>(i + 1))));
            (*this)[i]->value = values[i];
        }
    }
    
    ParamVec(const string &label, const unsigned size):
    label(label), size(size){
        values.setZero(size);
        for (unsigned i=0; i<size; ++i) {
            this->push_back(new Parameter(label + to_string(static_cast<long long>(i + 1))));
            (*this)[i]->value = 0;
        }
    }
};

class ParamSetVec : public vector<ParamSet*> {
    // a vector of sets of parameters, e.g, delta values for distribution component membership in BayesR model
public:
    string label;
    unsigned size;  // number of parameter sets
    
    ParamSetVec(const string &label, const vector<string> &header, const unsigned numSets):
    label(label), size(numSets){
        for (unsigned i=0; i<size; ++i) {
            this->push_back(new ParamSet(label + to_string(static_cast<long long>(i + 1)), header));
        }
    }
};


class Model {
public:
    string bayesType;
    string status;
        
    vector<ParamSet*> paramSetVec;
    vector<Parameter*> paramVec;
    vector<Parameter*> paramToPrint;
    vector<ParamSet*> paramSetToPrint;
    
    virtual void sampleUnknowns(const unsigned iter) = 0;
    virtual void sampleStartVal(void) = 0;
};


class BayesC : public Model {
    // model settings and prior specifications in class constructors
public:
    
    class FixedEffects : public ParamSet, public Stat::Flat {
        // all fixed effects has flat prior
    public:
        FixedEffects(const vector<string> &header, const string &lab = "CovEffects")
        : ParamSet(lab, header){}
        
        void sampleFromFC(VectorXf &ycorr, const MatrixXf &X, const VectorXf &XPXdiag, const float vare);
    };
    
    class RandomEffects : public ParamSet, public Stat::Normal {
        // random covariate effects
    public:
        float ssq;  // sum of squares

        RandomEffects(const vector<string> &header, const string &lab = "RandCovEffects")
        : ParamSet(lab, header){}
        
        void sampleFromFC(VectorXf &ycorr, const MatrixXf &W, const VectorXf &WPWdiag, const VectorXf &Rsqrt, const bool weightedRes, const float sigmaSqRand, const float vare, VectorXf &rhat);
    };
    
    class VarRandomEffects : public Parameter, public Stat::InvChiSq {
        // variance of random covariate effects has a scaled-inverse chi-square prior
    public:
        const float df;  // hyperparameter
        float scale;     // hyperparameter

        VarRandomEffects(const float varRandom, const float numRandomEffects, const string &lab = "SigmaSqRand")
        : Parameter(lab), df(4)
        {
            //value = varRandom/numRandomEffects;
            value = varRandom;
            scale = 0.5f*value/numRandomEffects;  // due to df = 4
        }
        
        void sampleFromFC(const float randEffSumSq, const unsigned numRandEff);
    };

    
    class SnpEffects : public ParamSet, public Stat::NormalZeroMixture {
        // all snp effects has a mixture prior of a nomral distribution and a point mass at zero
    public:
        float sumSq;
        unsigned numNonZeros;
        
        VectorXf posteriorMean;
        VectorXf posteriorMeanPIP;
        VectorXf pip;
        
        enum {gibbs, hmc} algorithm;
        
        unsigned cnt;
        float mhr;
        
        bool shuffle;
        vector<int> snpIndexVec;
        
        SnpEffects(const vector<string> &header, const string &alg, const string &lab = "SnpEffects")
        : ParamSet(lab, header){
            sumSq = 0.0;
            numNonZeros = 0;
            posteriorMean.setZero(size);
            posteriorMeanPIP.setZero(size);
            pip.setZero(size);
            if (alg=="HMC") algorithm = hmc;
            else algorithm = gibbs;
            cnt = 0;
            mhr = 0.0;
            
            shuffle = true;  // shuffle SNP order
            snpIndexVec.resize(size);
            for (unsigned i=0; i<size; i++) snpIndexVec[i] = i;
        }
        
        void sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes,
                          const float sigmaSq, const float pi, const float vare, VectorXf &ghat);
        void gibbsSampler(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes,
                          const float sigmaSq, const float pi, const float vare, VectorXf &ghat);
        void hmcSampler(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag,
                        const float sigmaSq, const float pi, const float vare, VectorXf &ghat);
        ArrayXf gradientU(const VectorXf &alpha, const MatrixXf &ZPZ, const VectorXf &ypZ,
                        const float sigmaSq, const float vare);
        float computeU(const VectorXf &alpha, const MatrixXf &ZPZ, const VectorXf &ypZ,
                       const float sigmaSq, const float vare);
        
        void sampleFromFC_omp(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag,
                              const float sigmaSq, const float pi, const float vare, VectorXf &ghat);
        void computePosteriorMean(const unsigned iter);

    };
    
    class SnpPIP : public ParamSet {
    public:
        SnpPIP(const vector<string> &header, const string &lab = "PIP") : ParamSet(lab, header){}
        
        void getValues(const VectorXf &pip){values = pip;}
    };
    
    class VarEffects : public Parameter, public Stat::InvChiSq {
        // variance of snp effects has a scaled-inverse chi-square prior
    public:
        const float df;  // hyperparameter
        float scale;     // hyperparameter
        bool noscale;  // no scaling on the genotypes

        VarEffects(const float vg, const VectorXf &snp2pq, const float pi, const bool noscale, const string &lab = "SigmaSq")
        : Parameter(lab), df(4), noscale(noscale)
        {
            // cout << "To scale or not to scale " << noscale << endl;
            // cout << "Scale value 1 " << value << endl;
            if (noscale) {
                value = vg / (snp2pq.sum() * pi);  // derived from prior knowledge on Vg and pi
            } else {
                value = vg / (snp2pq.size() * pi);  // derived from prior knowledge on Vg and pi
            }
            
            scale = 0.5f*value;  // due to df = 4
            
            //cout << value << " " << vg << " " << snp2pq.sum() << " " << pi << " " << noscale << endl;
        }
        
        void sampleFromFC(const float snpEffSumSq, const unsigned numSnpEff);
        void sampleFromPrior(void);
        void computeScale(const float varg, const VectorXf &snp2pq, const float pi);
        void computeScale(const float varg, const float sum2pq);
        void compute(const float snpEffSumSq, const float numSnpEff);

    };
    
    class ScaleVar : public Parameter, public Stat::Gamma {
        // scale factor of variance variable
    public:
        const float shape;
        const float scale;
        
        ScaleVar(const float val, const string &lab = "Scale"): shape(1.0), scale(1.0), Parameter(lab){
            value = val;  // starting value
        }
        
        void sampleFromFC(const float sigmaSq, const float df, float &scaleVar);
        void getValue(const float val){ value = val; };
    };
    
    class Pi : public Parameter, public Stat::Beta {
        // prior probability of a snp with a non-zero effect has a beta prior
    public:
        const float alpha;  // hyperparameter
        const float beta;   // hyperparameter
        
        Pi(const float pi, const float alpha, const float beta, const string &lab = "Pi"): Parameter(lab), alpha(alpha), beta(beta){  // informative prior
            value = pi;
        }
        
        void sampleFromFC(const unsigned numSnps, const unsigned numSnpEff);
        void sampleFromPrior(void);
        void compute(const float numSnps, const float numSnpEff);
    };
    
    
    class ResidualVar : public Parameter, public Stat::InvChiSq {
        // residual variance has a scaled-inverse chi-square prior
    public:
        const float df;      // hyperparameter
        const float scale;   // hyperparameter
        unsigned nobs;
        
        ResidualVar(const float vare, const unsigned n, const string &lab = "ResVar")
        : Parameter(lab), df(4)
        , scale(0.5f*vare){
            nobs = n;
            value = vare;  // due to df = 4
        }
        
        void sampleFromFC(VectorXf &ycorr);
    };
    
    class GenotypicVar : public Parameter {
        // compute genotypic variance from the sampled SNP effects
        // strictly speaking, this is not a model parameter
    public:
        GenotypicVar(const float varg, const string &lab = "GenVar"): Parameter(lab){
            value = varg;
        };
        void compute(const VectorXf &ghat);
    };
    
    class RandomVar : public Parameter {
        // compute variance explained due to random covariate effects
    public:
        RandomVar(const float varRandom, const string &lab = "RanVar"): Parameter (lab){
            value = varRandom;
        }
        
        void compute(const VectorXf &rhat);
    };
    
    class Heritability : public Parameter {
        // compute heritability based on sampled values of genotypic and residual variances
        // strictly speaking, this is not a model parameter
    public:
        Heritability(const string &lab = "hsq"): Parameter(lab){};
        void compute(const float genVar, const float resVar){
            value = genVar/(genVar+resVar);
        }
    };

    /** Marginal h^2 from --fixed-effect SNPs: sum_i alpha_i^2 / varPhenotypic (LD ignored). Effects are on SNP genotype SD scale (Var(X)=1), so no 2pq weighting. */
    class HeritabilityFixedEffects : public Parameter {
    public:
        HeritabilityFixedEffects(const string &lab = "hsqFixed"): Parameter(lab){}
        void compute(const VectorXf &snpEffects, const Data &data);
    };
    
    class Rounding : public Parameter {
        // re-compute ycorr to eliminate rounding errors
    public:
        Rounding(const string &lab = "Rounding"): Parameter(lab){}
        void computeYcorr(const VectorXf &y, const MatrixXf &X, const MatrixXf &W, const MatrixXf &Z,
                          const VectorXf &fixedEffects, const VectorXf &randomEffects, const VectorXf &snpEffects,
                          VectorXf &ycorr);
    };
    
    class NumNonZeroSnp : public Parameter {
        // number of non-zero SNP effects
    public:
        NumNonZeroSnp(const string &lab = "NnzSnp"): Parameter(lab){};
        void getValue(const unsigned nnz){ value = nnz; };
    };

    class varEffectScaled : public Parameter {
        // Alternative way to estimate genetic variance: sum 2pq sigmaSq
    public:
        varEffectScaled(const string &lab = "SigmaSqG"): Parameter(lab){};
        void compute(const float sigmaSq, const float sum2pq){value = sigmaSq*sum2pq;};
    };

    class SnpHsqPEP : public ParamSet {
        // per-SNP heritability posterior enrichment probability
    public:
        SnpHsqPEP(const vector<string> &header, const string &lab = "PEP") : ParamSet(lab, header){}
        
        void compute(const VectorXf &snpEffects, const float varg);
    };
    
public:
    const Data &data;
    
    VectorXf ycorr;   // corrected y for mcmc sampling
    VectorXf ghat;    // predicted total genotypic values
    VectorXf rhat;    // predicted total random covariate values
    
    bool estimatePi;
        
    FixedEffects fixedEffects;
    RandomEffects randomEffects;
    SnpEffects snpEffects;
    SnpPIP snpPip;
    VarEffects sigmaSq;
    VarRandomEffects sigmaSqRand;
    ScaleVar scale;
    Pi pi;
    ResidualVar vare;
    
    GenotypicVar varg;
    Heritability hsq;
    RandomVar varRand;
    Rounding rounding;
    NumNonZeroSnp nnzSnp;
    
    BayesC(const Data &data, const float varGenotypic, const float varResidual, const float varRandom, const float pival, const float piAlpha, const float piBeta, const bool estimatePi, const bool noscale,
           const string &algorithm = "Gibbs", const bool message = true):
    data(data),
    ycorr(data.y),
    fixedEffects(data.fixedEffectNames),
    randomEffects(data.randomEffectNames),
    sigmaSqRand(varRandom, data.numRandomEffects),
    snpEffects(data.snpEffectNames, algorithm),
    snpPip(data.snpEffectNames),
    sigmaSq(varGenotypic, data.snp2pq, pival, noscale),
    scale(sigmaSq.scale),
    pi(pival, piAlpha, piBeta),
    vare(varResidual, data.numKeptInds),
    varg(varGenotypic),
    varRand(varRandom),
    estimatePi(estimatePi)
    {
        bayesType = "C";
        paramSetVec = {&snpEffects, &fixedEffects, &snpPip};           // for which collect mcmc samples
        paramVec = {&pi, &nnzSnp, &sigmaSq, &varg, &vare, &hsq};       // for which collect mcmc samples
        paramToPrint = {&pi, &nnzSnp, &sigmaSq, &varg, &vare, &hsq};   // print in order
        if (data.numRandomEffects) {
            paramSetVec.push_back(&randomEffects);
            paramVec.push_back(&sigmaSqRand);
            paramVec.push_back(&varRand);
            paramToPrint.push_back(&varRand);
        }
        //paramToPrint.push_back(&rounding);
        if (message) {
            string alg = algorithm;
            if (alg!="HMC") alg = "Gibbs (default)";
            cout << "\nBayesC model fitted. Algorithm: " << alg << "." << endl;
            cout << "scale factor: " << sigmaSq.scale << endl;
            if (noscale) {
                cout << "Fitting model assuming unscaled genotypes " << endl;
            } else {
                cout << "Fitting model assuming scaled genotypes " << endl;
            }
        }
    }
    
    void sampleUnknowns(const unsigned iter);
    void sampleStartVal(void);
};

class BayesB : public BayesC {
public:
    
    class SnpEffects : public BayesC::SnpEffects {
    // for the ease of sampling, we model the SNP effect to be alpha_j = beta_j * delta_j where beta_j has a univariate normal prior.
    public:
        VectorXf betaSq;     // save sample squres of full conditional normal distribution regardless of delta values
        
        SnpEffects(const vector<string> &header): BayesC::SnpEffects(header, "Gibbs"){
            betaSq.setZero(size);
        }
        
        void sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes,
                          const VectorXf &sigmaSq, const float pi, const float vare, VectorXf &ghat);
    };

    class VarEffects : public ParamSet, public BayesC::VarEffects {
    public:
        VarEffects(const float vg, const VectorXf &snp2pq, const float pi, const bool noscale):
        ParamSet("SigmaSqs", vector<string>(snp2pq.size())),
        BayesC::VarEffects(vg, snp2pq, pi, noscale){
            values.setConstant(size, value);
        }
        
        void sampleFromFC(const VectorXf &betaSq);
    };
    
    SnpEffects snpEffects;
    VarEffects sigmaSq;

    BayesB(const Data &data, const float varGenotypic, const float varResidual, const float varRandom, const float pival, const float piAlpha, const float piBeta,
           const bool estimatePi, const bool noscale, const bool message = true):
    BayesC(data, varGenotypic, varResidual, varRandom, pival, piAlpha, piBeta, estimatePi, noscale, "Gibbs", false),
    snpEffects(data.snpEffectNames),
    sigmaSq(varGenotypic, data.snp2pq, pival, noscale)
    {
        bayesType = "B";
        paramSetVec = {&snpEffects, &fixedEffects, &snpPip};           // for which collect mcmc samples
        paramVec = {&pi, &nnzSnp, &varg, &vare, &hsq};       // for which collect mcmc samples
        paramToPrint = {&pi, &nnzSnp, &varg, &vare, &hsq};   // print in order
        if (data.numRandomEffects) {
            paramSetVec.push_back(&randomEffects);
            paramVec.push_back(&sigmaSqRand);
            paramVec.push_back(&varRand);
            paramToPrint.push_back(&varRand);
        }
        //paramToPrint.push_back(&rounding);
        if (message) {
            cout << "\nBayesB model fitted." << endl;
            cout << "scale factor: " << sigmaSq.scale << endl;
            if (noscale) {
                cout << "Fitting model assuming unscaled genotypes " << endl;
            } else {
                cout << "Fitting model assuming scaled genotypes " << endl;
            }
        }
    }
    
    void sampleUnknowns(const unsigned iter);

};

class BayesN : public BayesC {
    // Nested model
public:
    
    class WindowDelta : public ParamSet {
    public:
        WindowDelta(const vector<string> &header, const string &lab = "WindowDelta"): ParamSet(lab, header){}
        void getValues(const VectorXf &val){ values = val; };
    };
    
    class SnpEffects : public BayesC::SnpEffects {
    public:
        unsigned numWindows;
        unsigned numNonZeroWind;
        
        const VectorXi &windStart;
        const VectorXi &windSize;
        
        VectorXf localPi, logLocalPi, logLocalPiComp;
        VectorXf windDelta;
        VectorXf snpDelta;
        VectorXf beta;     // save samples of full conditional normal distribution regardless of delta values
        ArrayXf cumDelta;  // for Polya urn proposal
        
        VectorXf windPip;
        
        SnpEffects(const vector<string> &header, const VectorXi &windStart, const VectorXi &windSize, const unsigned snpFittedPerWindow):
        BayesC::SnpEffects(header, "Gibbs"), windStart(windStart), windSize(windSize){
            numWindows = (unsigned) windStart.size();
            windDelta.setZero(numWindows);
            localPi.setOnes(numWindows);
            windPip.setZero(numWindows);
            snpDelta.setZero(size);
            beta.setZero(size);
            cumDelta.setZero(size);
            for (unsigned i=0; i<numWindows; ++i) {
                if (snpFittedPerWindow < windSize[i])
                    localPi[i] = snpFittedPerWindow/float(windSize[i]);
            }
            logLocalPi = localPi.array().log().matrix();
            logLocalPiComp = (1.0f-localPi.array()).log().matrix();
        }

        void sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes,
                          const float sigmaSq, const float pi, const float vare, VectorXf &ghat);
    };
    
    class VarEffects : public BayesC::VarEffects {
    public:
        VarEffects(const float vg, const VectorXf &snp2pq, const float pi,
                   const VectorXf &localPi, const unsigned snpFittedPerWindow):
        BayesC::VarEffects(vg, snp2pq, pi, false){
            value /= localPi.mean();
            scale = 0.5*value;
        }
    };
    
    class NumNonZeroWind : public Parameter {
        // number of non-zero window effects
    public:
        NumNonZeroWind(const string &lab = "NNZwind"): Parameter(lab){};
        void getValue(const unsigned nnz){ value = nnz; };
    };
    
    
    SnpEffects snpEffects;
    VarEffects sigmaSq;
    NumNonZeroWind nnzWind;
    WindowDelta windDelta;
    
    BayesN(const Data &data, const float varGenotypic, const float varResidual, const float varRandom, const float pival, const float piAlpha, const float piBeta,
           const bool estimatePi, const bool noscale, const unsigned snpFittedPerWindow, const bool message = true):
    BayesC(data, varGenotypic, varResidual, varRandom, pival, piAlpha, piBeta, estimatePi, noscale, "Gibbs", false),
    snpEffects(data.snpEffectNames, data.windStart, data.windSize, snpFittedPerWindow),
    sigmaSq(varGenotypic, data.snp2pq, pival, snpEffects.localPi, snpFittedPerWindow),
    windDelta(vector<string>(snpEffects.numWindows))
    {
        bayesType = "N";
        paramSetVec = {&snpEffects, &fixedEffects, &snpPip, &windDelta};           // for which collect mcmc samples
        paramVec = {&pi, &nnzWind, &nnzSnp, &sigmaSq, &varg, &vare, &hsq};       // for which collect mcmc samples
        paramToPrint = {&pi, &nnzWind, &nnzSnp, &sigmaSq, &varg, &vare, &hsq};   // print in order
        if (data.numRandomEffects) {
            paramSetVec.push_back(&randomEffects);
            paramVec.push_back(&sigmaSqRand);
            paramVec.push_back(&varRand);
            paramToPrint.push_back(&varRand);
        }
        //paramToPrint.push_back(&rounding);
        if (message) {
            cout << "\nBayesN model fitted." << endl;
            cout << "scale factor: " << sigmaSq.scale << endl;
            if (noscale) {
                cout << "Fitting model assuming unscaled genotypes " << endl;
            } else {
                cout << "Fitting model assuming scaled genotypes " << endl;
            }
        }
    }

    void sampleUnknowns(const unsigned iter);
};

// -----------------------------------------------------------------------------------------------
// Bayes R
// -----------------------------------------------------------------------------------------------

class BayesR : public BayesC {
    // Prior for snp efect pi_1 * N(0, 0) + pi_2 * N(0, sigma^2_beta * gamma_2) + pi_3 * N(0, sigma^2_beta * gamma_3) + pi_3 * N(0, sigma^2_beta * gamma_4)
public:

    class DeltaPi : public ParamSetVec {
    public:
        DeltaPi(const vector<string> &header, const unsigned numDist, const string &lab = "DeltaPi"):
        ParamSetVec(lab, header, numDist) {}
    };
    
    class SnpEffects : public BayesC::SnpEffects {
    public:
        vector<vector<unsigned> > snpset;
        float sum2pq;
        float wtdSumSq;
        
        SnpEffects(const vector<string> &header, const string &alg): BayesC::SnpEffects(header, "Gibbs"){
            sum2pq = 0.0;
            wtdSumSq = 0.0;
        }
        
        void sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes,
                          const float sigmaSq, const VectorXf &pis,  const VectorXf &gamma,
                          const float vare, VectorXf &ghat, VectorXf &snpStore,
                          const float varg, const bool hsqPercModel, DeltaPi &deltaPi);
    };

    class ProbMixComps : public ParamVec, public Stat::Dirichlet {
        // prior probability of a snp being in any of the distributions effect has a dirichlet prior
    public:
        VectorXf alphaVec;  // hyperparameter

        ProbMixComps(const VectorXf &pis, const VectorXf &alphas, const string &lab = "Pi"):
        ParamVec(lab, pis){
            if (alphas.size() != size) alphaVec.setOnes(size);
            else alphaVec = alphas;
        }
        
        void sampleFromFC(const VectorXf &snpStore);
    };
    
    class VarEffects : public BayesC::VarEffects {
    public:
        VarEffects(const float vg, const VectorXf &snp2pq, const VectorXf &gamma, const VectorXf &pi, const bool noscale, const string &lab = "SigmaSq"):
        BayesC::VarEffects(vg, snp2pq, 1-pi[0], noscale, lab) {
            if (noscale) {
                value = vg / (snp2pq.sum() * gamma.dot(pi));  // derived from prior knowledge on Vg and pi
            } else {
                value = vg / (snp2pq.size() * gamma.dot(pi));  // derived from prior knowledge on Vg and pi
            }
            
            scale = (df-2)/df*value;
            
//            cout << "noscale " << noscale << " vg " << vg << " snp2pq.sum() " << snp2pq.sum() << " snp2pq.size() " << snp2pq.size() << " gamma.dot(pi) " << gamma.dot(pi) << " scale " << scale << endl;
        }
        
        void computeScale(const float varg, const VectorXf &snp2pq, const VectorXf &gamma, const VectorXf &pi);
   };
    
    class VgMixComps : public ParamVec {
    public:
        unsigned zeroIdx, minIdx;
        
        VgMixComps(const VectorXf &gamma, const string &lab = "Vg"):
        ParamVec(lab, gamma.size()){
            float min = 1.0;
            minIdx = 0;
            for (unsigned i=0; i<size; ++i) {
                if (gamma[i] == 0) zeroIdx = i;
                else if (gamma[i] < min) {
                    min = gamma[i];
                    minIdx = i;
                }
            }
        }
        
        void compute(const VectorXf &snpEffects, const MatrixXf &Z, const vector<vector<unsigned> > snpset, const float varg);
    };

    class NumSnpMixComps : public ParamVec {
    public:
        NumSnpMixComps(const VectorXf &pis, const string &lab = "NumSnp"):
        ParamVec(lab, pis.size()){}
        
        void getValues(const VectorXf &snpStore);
    };

    class Gammas : public ParamSet {
        // Set of scaling factors for each of the distributions
    public:
        Gammas(const VectorXf &gamma, const vector<string> &header, const string &lab = "gamma"): ParamSet(lab, header){
            values = gamma;
        }
    };
    
public:
    VectorXf snpStore;   
    SnpEffects snpEffects;
    VarEffects sigmaSq;
    ProbMixComps Pis;
    VgMixComps Vgs;
    NumSnpMixComps numSnps;
    Gammas gamma;
    DeltaPi deltaPi;

    bool hsqPercModel;

    BayesR(const Data &data, const float varGenotypic, const float varResidual, const float varRandom, const VectorXf pis, const VectorXf &piPar, const VectorXf gamma, const bool estimatePi, const bool noscale, const bool hsqPercModel,
           const string &algorithm, const bool message = true):
    BayesC(data, varGenotypic, varResidual, varRandom, 1-pis[0], piPar[0], piPar[1], estimatePi, noscale, "Gibbs", false),
    Pis(pis, piPar),
    numSnps(pis),
    Vgs(gamma),
    gamma(gamma, vector<string>(gamma.size())),
    snpEffects(data.snpEffectNames, algorithm),
    sigmaSq(varGenotypic, data.snp2pq, gamma, pis, noscale),
    deltaPi(data.snpEffectNames, pis.size()),
    hsqPercModel(hsqPercModel)
    {
        bayesType = "R";

        if (data.numKeptInds < 1000) this->hsqPercModel = false;  // when sample size is small, sample the common variance variable as hsq estimate may be instable
        
        paramSetVec  = {&snpEffects, &fixedEffects, &snpPip};
        paramSetVec.insert(paramSetVec.end(), deltaPi.begin(), deltaPi.end());
        paramVec     = {&nnzSnp, &sigmaSq, &vare, &varg, &hsq};
        paramVec.insert(paramVec.end(), numSnps.begin(), numSnps.end());
        paramVec.insert(paramVec.end(), Vgs.begin(), Vgs.end());
        
        paramToPrint = {&sigmaSq, &varg, &vare, &hsq};
        paramToPrint.insert(paramToPrint.begin(), Vgs.begin(), Vgs.end());
        paramToPrint.insert(paramToPrint.begin(), numSnps.begin(), numSnps.end());
        
        if (data.numRandomEffects) {
            paramSetVec.push_back(&randomEffects);
            paramVec.push_back(&sigmaSqRand);
            paramVec.push_back(&varRand);
            paramToPrint.push_back(&varRand);
        }
        
        if (message) {
            string alg = algorithm;
            if (alg!="HMC") alg = "Gibbs (default)";
            cout << "\nBayesR model fitted. Algorithm: " << alg << "." << endl;
            cout << "scale factor: " << sigmaSq.scale << endl;
            cout << "Gamma: " << gamma.transpose() << endl;
            if (noscale) {
                cout << "Fitting model assuming unscaled genotypes " << endl;
            } else {
                cout << "Fitting model assuming scaled genotypes " << endl;
            }
        }
    }
    
    void sampleUnknowns(const unsigned iter);
};
    

class BayesS : public BayesC {
    // Prior for snp efect alpha_j ~ N(0, sigma^2_a / (2p_j q_j)^S)
    // consider S as unknown to make inference on the relationship between MAF and effect size
public:
    
    class AcceptanceRate : public Parameter {
    public:
        unsigned cnt;
        unsigned accepted;
        unsigned consecRej;
        
        AcceptanceRate(): Parameter("AR"){
            cnt = 0;
            accepted = 0;
            value = 0.0;
            consecRej = 0;
        };
        
        void count(const bool state, const float lower, const float upper);
    };
    
    class Sp : public Parameter, public Stat::Normal {
        // S parameter for genotypes or equivalently for the variance of snp effects
        
        // random-walk MH and HMC algorithms implemented
        
    public:
        const float mean;  // prior
        const float var;   // prior
        const unsigned numSnps;
        
        float varProp;     // variance of proposal normal for random walk MH
        
        float stepSize;     // for HMC
        unsigned numSteps;  // for HMC
        
        enum {random_walk, hmc, reg} algorithm;
        
        AcceptanceRate ar;
        Parameter tuner;
        
        Sp(const unsigned m, const float var, const float start, const string &alg, const string &lab = "S"): Parameter(lab), mean(0), var(var), numSnps(m)
        , tuner(alg=="RWMH" ? "varProp" : "Stepsize"){
            value = start;  // starting value
            varProp = 0.01;
            stepSize = 0.001;
            numSteps = 100;
            if (alg=="RWMH") algorithm = random_walk;
            else if (alg=="Reg") algorithm = reg;
            else algorithm = hmc;
            //else throw("Error: Invalid algorithm for sampling S: " + alg + " (the available are RWMH, HMC, Reg)!");
        }
        
        // note that the scale factor of sigmaSq will be simultaneously updated
        void sampleFromFC(const float snpEffWtdSumSq, const unsigned numNonZeros, float &sigmaSq, const VectorXf &snpEffects,
                          const VectorXf &snp2pq, ArrayXf &snp2pqPowS, const ArrayXf &logSnp2pq,
                          const float vg, float &scale, float &sum2pqSplusOne, bool scaledGeno);
        void sampleFromPrior(void);
        void randomWalkMHsampler(const float snpEffWtdSumSq, const unsigned numNonZeros, const float sigmaSq, const VectorXf &snpEffects,
                                 const VectorXf &snp2pq, ArrayXf &snp2pqPowS, const ArrayXf &logSnp2pq,
                                 const float vg, float &scale, float &sum2pqSplusOne);
        void hmcSampler(const unsigned numNonZeros, const float sigmaSq, const VectorXf &snpEffects,
                        const VectorXf &snp2pq, ArrayXf &snp2pqPowS, const ArrayXf &logSnp2pq,
                        const float vg, float &scale, float &sum2pqSplusOne, bool scaledGeno);
        float gradientU(const float S, const ArrayXf &snpEffects, const float snp2pqLogSum, const ArrayXf &snp2pq, const ArrayXf &logSnp2pq, const float sigmaSq, const float vg, bool scaledGeno);
        float computeU(const float S, const ArrayXf &snpEffects, const float snp2pqLogSum, const ArrayXf &snp2pq, const ArrayXf &logSnp2pq, const float sigmaSq, const float vg, float &scale, float &U_chisq, bool scaledGeno);
        void regression(const VectorXf &snpEffects, const ArrayXf &logSnp2pq, ArrayXf &snp2pqPowS, float &sigmaSq);
 
    };
    
    class SnpEffects : public BayesC::SnpEffects {
    public:
        float wtdSumSq;  // weighted sum of squares by 2pq^S
        float sum2pqSplusOne;  // sum of delta_j* (2p_j q_j)^{1+S}
        
        SnpEffects(const vector<string> &header, const VectorXf &snp2pq, const float pi): BayesC::SnpEffects(header, "Gibbs") {
            wtdSumSq = 0.0;
            sum2pqSplusOne = snp2pq.sum()*pi;  // starting value of S is 0
        }
        
        void sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes,
                          const float sigmaSq, const float pi, const float vare,
                          const ArrayXf &snp2pqPowS, const VectorXf &snp2pq,
                          const float vg, float &scale, VectorXf &ghat);
    };
    
    
public:
    float genVarPrior;
    float scalePrior;

    ArrayXf snp2pqPowS;
    ArrayXf snp2pqPowSplusOne;  // snp2pq^S-1
    const ArrayXf logSnp2pq;
    
    Sp S;
    SnpEffects snpEffects;
    
    bool scaledGeno;  // whether genotypes are scaled
    
    BayesS(const Data &data, const float varGenotypic, const float varResidual, const float varRandom, const float pival, const float piAlpha, const float piBeta, const bool estimatePi, const float varS, const vector<float> &svalue,
           const string &algorithm, const bool noscale = false, const bool message = true):
    BayesC(data, varGenotypic, varResidual, varRandom, pival, piAlpha, piBeta, estimatePi, noscale, "Gibbs", false),
    logSnp2pq(data.snp2pq.array().log()),
    S(data.numIncdSnps, varS, svalue[0], algorithm),
    snpEffects(data.snpEffectNames, data.snp2pq, pival),
    genVarPrior(varGenotypic),
    scalePrior(sigmaSq.scale)
    {
        scaledGeno = !noscale;
        
        bayesType = "S";
        
        findStartValueForS(svalue);

        snp2pqPowS = scaledGeno ? data.snp2pq.array().pow(S.value + 1.0f) : data.snp2pq.array().pow(S.value);
        sigmaSq.value = scaledGeno ? varGenotypic/(snp2pqPowS.sum()*pival) : varGenotypic/((snp2pqPowS*data.snp2pq.array()).sum()*pival);        
        scale.value = sigmaSq.scale = 0.5*sigmaSq.value;
 
        paramSetVec = {&snpEffects, &fixedEffects, &snpPip};
        paramVec = {&pi, &nnzSnp, &sigmaSq, &S, &varg, &vare, &hsq};
        paramToPrint = {&pi, &nnzSnp, &sigmaSq, &scale, &S, &varg, &vare, &hsq, &S.ar, &S.tuner};
        if (data.numRandomEffects) {
            paramSetVec.push_back(&randomEffects);
            paramVec.push_back(&sigmaSqRand);
            paramVec.push_back(&varRand);
            paramToPrint.push_back(&varRand);
        }
        //paramToPrint.push_back(&rounding);
        if (message) {
            string alg = algorithm;
            if (alg!="RWMH" && alg!="Reg") alg = "HMC";
            cout << "\nBayesS model fitted. Algorithm: " << alg << "." << endl;
            cout << "scale factor: " << sigmaSq.scale << endl;
            if (scaledGeno) {
                cout << "Fitting model assuming scaled genotypes " << endl;
            } else {
                cout << "Fitting model assuming unscaled genotypes " << endl;
            }
        }
    }
    
    void sampleUnknowns(const unsigned iter);
    void sampleStartVal(void);
    void findStartValueForS(const vector<float> &val);
    float computeLogLikelihood(void);
    void sampleUnknownsWarmup(void);
};


class BayesNS : public BayesS {
    // combine BayesN and BayesS primarily for speed
public:
    
    class SnpEffects : public BayesN::SnpEffects {
    public:
        unsigned burnin;
        
        float wtdSumSq;  // weighted sum of squares by 2pq^S
        float sum2pqSplusOne;  // sum of delta_j* (2p_j q_j)^{1+S}
        
        ArrayXf varPseudoPrior;
        
        SnpEffects(const vector<string> &header, const VectorXi &windStart, const VectorXi &windSize,
                   const unsigned snpFittedPerWindow, const VectorXf &snp2pq, const float pi):
        BayesN::SnpEffects(header, windStart, windSize, snpFittedPerWindow){
            burnin = 2000;
            wtdSumSq = 0.0;
            sum2pqSplusOne = 0.0;
            //sum2pqSplusOne = snp2pq.sum()*(1.0f-pi)*(1.0f-snpFittedPerWindow/float(windSize));  // starting value of S is 0
            varPseudoPrior.setZero(size);
        }
        
        void sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes,
                          const float sigmaSq, const float pi, const float vare,
                          const ArrayXf &snp2pqPowS, const VectorXf &snp2pq,
                          const float vg, float &scale, VectorXf &ghat);

    };
    
    class Sp : public BayesS::Sp {  //**** NOT working ****
        // difference to BayesS::Sp is that since a gamma prior is given to the scale factor of sigmaSq,
        // S parameter is no longer present in the density function of sigmaSq
    public:
        Sp(const unsigned numSnps, const float var, const float start, const string &alg): BayesS::Sp(numSnps, var, start, alg){}
        
        void sampleFromFC(const unsigned numNonZeros, const float sigmaSq, const VectorXf &snpEffects,
                          const VectorXf &snp2pq, ArrayXf &snp2pqPowS, const ArrayXf &logSnp2pq);
        float gradientU(const float S, const ArrayXf &snpEffects, const float snp2pqLogSum,
                        const ArrayXf &snp2pq, const ArrayXf &logSnp2pq, const float sigmaSq);
        float computeU(const float S, const ArrayXf &snpEffects, const float snp2pqLogSum,
                       const ArrayXf &snp2pq, const ArrayXf &logSnp2pq, const float sigmaSq);
    };
    
    SnpEffects snpEffects;
    //Sp S;
    BayesN::VarEffects sigmaSq;
    BayesC::ScaleVar scale;
    BayesN::NumNonZeroWind nnzWind;
    BayesN::WindowDelta windDelta;
    
    BayesNS(const Data &data, const float varGenotypic, const float varResidual, const float varRandom, const float pival, const float piAlpha, const float piBeta,
            const bool estimatePi, const float varS, const vector<float> &svalue, const unsigned snpFittedPerWindow,
            const string &algorithm, const bool message = true):
    BayesS(data, varGenotypic, varResidual, varRandom, pival, piAlpha, piBeta, estimatePi, varS, svalue, algorithm, false),
    snpEffects(data.snpEffectNames, data.windStart, data.windSize, snpFittedPerWindow, data.snp2pq, pival),
    //S(data.numIncdSnps, "HMC"),
    sigmaSq(varGenotypic, data.snp2pq, pival, snpEffects.localPi, snpFittedPerWindow),
    scale(sigmaSq.scale),
    windDelta(vector<string>(snpEffects.numWindows))
    {
        bayesType = "NS";
        paramSetVec = {&snpEffects, &fixedEffects, &snpPip, &windDelta};
        paramVec = {&pi, &nnzWind, &nnzSnp, &sigmaSq, &S, &varg, &vare, &hsq};
        paramToPrint = {&pi, &nnzWind, &nnzSnp, &sigmaSq, &scale, &S, &varg, &vare, &hsq, &S.ar, &S.tuner};
        if (data.numRandomEffects) {
            paramSetVec.push_back(&randomEffects);
            paramVec.push_back(&sigmaSqRand);
            paramVec.push_back(&varRand);
            paramToPrint.push_back(&varRand);
        }
        //paramToPrint.push_back(&rounding);
        if (message) {
            string alg = algorithm;
            if (alg!="RWMH" && alg!="Reg") alg = "HMC";
            cout << "\nBayesNS model fitted. Algorithm: " << alg << "." << endl;
            cout << "scale factor: " << sigmaSq.scale << endl;
        }
    }
    
    void sampleUnknowns(const unsigned iter);
};

class BayesRS : public BayesR {
public:
    
    class SnpEffects : public BayesS::SnpEffects {
    public:
        unsigned ndist;
        ArrayXf numSnpMix;
        vector<vector<unsigned> > snpset;
        
        SnpEffects(const vector<string> &header, const VectorXf &snp2pq, const VectorXf &pis): BayesS::SnpEffects(header, snp2pq, 1.0-pis[0]) {
            ndist = pis.size();
            numSnpMix.setZero(ndist);
        }
        
        void sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes,
                          const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare,
                          const ArrayXf &snp2pqPowS, const VectorXf &snp2pq,
                          const float varg, float &scale, VectorXf &ghat, const bool hsqPercModel);
    };
    
    class Sp : public BayesS::Sp {
    public:
        Sp(const unsigned m, const float var, const float start): BayesS::Sp(m, var, start, "HMC"){}
        
        void sampleFromFC(vector<vector<unsigned> > &snpset, const VectorXf &snpEffects,
                          float &sigmaSq, const VectorXf &gamma,
                          const VectorXf &snp2pq, ArrayXf &snp2pqPowS, const ArrayXf &logSnp2pq,
                          const float vg, float &scale, float &sum2pqSplusOne);
        float gradientU(const float S, const unsigned nnzMix, const vector<ArrayXf> &snpEffectMix, const float snp2pqLogSum, const vector<ArrayXf> &snp2pqMix, const vector<ArrayXf> &logSnp2pqMix, const float sigmaSq, const VectorXf &gamma, const float vg);
        float computeU(const float S, const unsigned nnzMix, const vector<ArrayXf> &snpEffectMix, const float snp2pqLogSum, const vector<ArrayXf> &snp2pqMix, const vector<ArrayXf> &logSnp2pqMix, const float sigmaSq, const VectorXf &gamma, const float vg, float &scale);
    };
    
    SnpEffects snpEffects;
    Sp S;
    
    ArrayXf logSnp2pq;
    ArrayXf snp2pqPowS;

    float genVarPrior;
    float scalePrior;

    bool scaledGeno;  // whether genotypes are scaled
    
    BayesRS(const Data &data, const float varGenotypic, const float varResidual, const float varRandom, const VectorXf pis, const VectorXf &piPar, const VectorXf gamma, const bool estimatePi, const float varS, const vector<float> &svalue, const bool noscale, const bool hsqPercModel, const string &algorithm, const bool message = true):
    BayesR(data, varGenotypic, varResidual, varRandom, pis, piPar, gamma, estimatePi, noscale, hsqPercModel, algorithm, false),
    snpEffects(data.snpEffectNames, data.snp2pq, pis),
    S(data.numIncdSnps, varS, svalue[0]),
    genVarPrior(varGenotypic),
    scalePrior(sigmaSq.scale)
    {
        scaledGeno = !noscale;  // scaled genotypes by default
        
        bayesType = "RS";
        
        logSnp2pq = data.snp2pq.array().log();
        snp2pqPowS = data.snp2pq.array().pow(S.value);
        sigmaSq.value = varGenotypic/((snp2pqPowS*data.snp2pq.array()).sum()*(1.0-pis[0]));
        scale.value = sigmaSq.scale = 0.5*sigmaSq.value;

        paramSetVec = {&snpEffects, &fixedEffects, &snpPip};
        paramSetVec.insert(paramSetVec.end(), deltaPi.begin(), deltaPi.end());
        paramVec    = {&nnzSnp, &sigmaSq, &S, &varg, &vare, &hsq};
        paramVec.insert(paramVec.end(), numSnps.begin(), numSnps.end());
        paramVec.insert(paramVec.end(), Vgs.begin(), Vgs.end());
        
        paramToPrint = {&sigmaSq, &S, &varg, &vare, &hsq};
        paramToPrint.insert(paramToPrint.begin(), Vgs.begin(), Vgs.end());
        paramToPrint.insert(paramToPrint.begin(), numSnps.begin(), numSnps.end());
        
        if (data.numRandomEffects) {
            paramSetVec.push_back(&randomEffects);
            paramVec.push_back(&sigmaSqRand);
            paramVec.push_back(&varRand);
            paramToPrint.push_back(&varRand);
        }
        
        if (message) {
            cout << "\nBayesRS model fitted." << endl;
            cout << "scale factor: " << sigmaSq.scale << endl;
            cout << "Gamma: " << gamma.transpose() << endl;
            if (noscale)
            {
                cout << "Fitting model assuming unscaled genotypes " << endl;
            } else
            {
                cout << "Fitting model assuming scaled genotypes "  << endl;
            }
        }
    }
    
    void sampleUnknowns(const unsigned iter);
};


class ApproxBayesC : public BayesC {
public:
    class SnpEffects; // forward decl for FixedEffects Gibbs helpers

    class FixedEffects : public BayesC::FixedEffects {
    public:
        FixedEffects(const vector<string> &header): BayesC::FixedEffects(header){}
        
        void sampleFromFC(const MatrixXf &XPX, const VectorXf &XPXdiag,
                          const MatrixXf &ZPX, const VectorXf &XPy,
                          const VectorXf &snpEffects, const float vare,
                          VectorXf &rcorr);

        /** Gibbs sweep over all --fixed-effect SNPs (call before mixture SnpEffects::sampleFromFC_*). */
        void sampleFromFC_eigen(SnpEffects &snpEffects, vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
            const vector<LDBlockInfo*> &keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks, const unsigned numFixedEffectSnps);
        void sampleFromFC_sparse(SnpEffects &snpEffects, VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag,
            const VectorXf &snp2pq, const vector<SnpInfo*> &incdSnpInfoVec, float varg, float vare, const unsigned numFixedEffectSnps);
        void sampleFromFC_full(SnpEffects &snpEffects, VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag,
            const VectorXi &windStart, const VectorXi &windSize, const VectorXf &snp2pq, const vector<SnpInfo*> &incdSnpInfoVec, float varg, float vare, const unsigned numFixedEffectSnps);
    };
    
    class SnpEffects : public BayesC::SnpEffects {
    public:
        float sum2pq;
        
        VectorXf nnzPerChr;
        VectorXf nnzPerBlk;
        VectorXi leaveout;
        VectorXf ssqBlocks;
        VectorXi badSnps;

        SnpEffects(const vector<string> &header): BayesC::SnpEffects(header, "Gibbs"){
            sum2pq = 0.0;
            leaveout.setZero(size);
            badSnps.setZero(size);
        }
        
        void sampleFromFC_sparse(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                 const vector<ChromInfo*> &chromInfoVec, const VectorXf &snp2pq,
                                 const float sigmaSq, const float pi, const float vare, const float varg,
                                 const vector<SnpInfo*> &incdSnpInfoVec);
        void sampleFromFC_full(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXf &snp2pq,
                          const float sigmaSq, const float pi, const float vare, const float varg,
                          const vector<SnpInfo*> &incdSnpInfoVec);
        void sampleFromFC_eigen(vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
                          const vector<LDBlockInfo*> keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks,
                          const float sigmaSq, const float pi, const float varg, const VectorXf &snp2pq);
        
        void hmcSampler(VectorXf &rcorr, const VectorXf &ZPy, const vector<VectorXf> &ZPZ,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const float sigmaSq, const float pi, const float vare);
        VectorXf gradientU(const VectorXf &effects, VectorXf &rcorr, const VectorXf &ZPy, const vector<VectorXf> &ZPZ,
                                                     const VectorXi &windStart, const VectorXi &windSize, const unsigned chrStart, const unsigned chrSize,
                                                     const float sigmaSq, const float vare);
        float computeU(const VectorXf &effects, const VectorXf &rcorr, const VectorXf &ZPy, const float sigmaSq, const float vare);
        
        void computeFromBLUP(vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
                          const vector<LDBlockInfo*> keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks,
                          const float sigmaSq, const float pi, const float varg, const VectorXf &snp2pq);

    };
    
    class VarEffects : public BayesC::VarEffects {
    public:
        VarEffects(const float varg, const VectorXf &snp2pq, const float pi, const bool noscale):
        BayesC::VarEffects(varg, snp2pq, pi, noscale){}
        
        void computeRobustMode(const float varg, const VectorXf &snp2pq, const float pi, const bool noscale);
    };
    
    class ResidualVar : public BayesC::ResidualVar {
    public:
        unsigned nNegVal;
        ResidualVar(const float vare, const unsigned nobs): BayesC::ResidualVar(vare, nobs){
            nNegVal = 0;
        }
        
        void sampleFromFC(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr, string &message);
    };
    
    class GenotypicVar : public BayesC::GenotypicVar {
    public:
        const unsigned nobs;
        
        GenotypicVar(const float varg, const unsigned n): BayesC::GenotypicVar(varg), nobs(n){}
        void compute(const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr);
    };
    
    class BlockGenotypicVar : public ParamSet, public BayesC::GenotypicVar {
    public:
        const unsigned nobs;
        unsigned numBlocks;
        float total;
        
        BlockGenotypicVar(const vector<string> &header, const float varg, const unsigned n, const string &lab = "BlockGenVar"):
        ParamSet(lab, header), BayesC::GenotypicVar(varg), nobs(n){
            numBlocks = header.size();
            total = 0.0;
        }

        void compute(const vector<VectorXf> &whatBlocks);
    };
    
    class BlockResidualVar : public ParamSet, public Stat::InvChiSq {
    public:
        const float df;      // hyperparameter
        const float scale;   // hyperparameter
        
        const float vary;

        unsigned numBlocks;
        float threshold;
        float mean;
        
        BlockResidualVar(const vector<string> &header, const float varPhenotypic, const string &lab = "BlockResVar"):
        ParamSet(lab, header), df(4), scale(0.5f*varPhenotypic), vary(varPhenotypic) {
            values.setConstant(size, varPhenotypic);
            numBlocks = header.size();
            threshold = 1.1;
            mean = varPhenotypic;
        }
        
        void sampleFromFC(vector<VectorXf> &wcorrBlocks, VectorXf &vargBlocks, VectorXf &ssqBlocks, const VectorXf &nGWASblocks, const VectorXf &numEigenvalBlock);
        void sampleFromFC(const vector<VectorXf> &wcorrBlocks, const VectorXf &beta, const VectorXf &b, const VectorXf &nGWASblocks, const vector<LDBlockInfo*> keptLdBlockInfoVec);
    };


    class Rounding : public BayesC::Rounding {
    public:
        void computeRcorr_sparse(const VectorXf &ZPy, const vector<SparseVector<float> > &ZPZsp,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXf &snpEffects, VectorXf &rcorr);
        void computeRcorr_full(const VectorXf &ZPy, const vector<VectorXf> &ZPZ,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXf &snpEffects, VectorXf &rcorr);
        void computeWcorr_eigen(const vector<VectorXf> &wBlocks, const vector<MatrixXf> &Qblocks, const vector<LDBlockInfo*> keptLdBlockInfoVec,
                                const VectorXf &snpEffects, vector<VectorXf> &wcorrBlocks);
        void computeGhat(const MatrixXf &Z, const VectorXf &snpEffects, VectorXf &ghat);
    };
        
//    class PopulationStratification : public Parameter, public Stat::InvChiSq {
//    public:
//        const float df;
//        const float scale;
//        
//        VectorXf chrSpecific;
//        
//        PopulationStratification(): Parameter("PS"), df(4), scale(0.5){
//            chrSpecific.setZero(22);
//        }
//        
//        void compute(const VectorXf &rcorr, const VectorXf &ZPZdiag, const VectorXf &LDsamplVar, const float varg, const float vare, const VectorXf &chisq);
//        void compute(const VectorXf &rcorr, const VectorXf &ZPZdiag, const VectorXf &LDsamplVar, const float varg, const float vare, const vector<ChromInfo*> chromInfoVec);
//    };
//    
//    class NumResidualOutlier : public Parameter {
//    public:
//        ofstream out;
//        
//        NumResidualOutlier(): Parameter("Nro"){}
//        
//        void compute(const VectorXf &rcorr, const VectorXf &ZPZdiag, const VectorXf &LDsamplVar, const float varg, const float vare, const vector<string> &snpName, VectorXi &leaveout, const vector<SparseVector<float> > &ZPZ, const VectorXf &ZPy, const VectorXf &snpEffects);
//    };
//    
//    class InterChrGenetCov : public Parameter {
//    public:
//        const float spouseCorrelation;
//        const unsigned nobs;
//        
//        InterChrGenetCov(const float corr, const unsigned nobs): Parameter("GenCov"), spouseCorrelation(corr), nobs(nobs) {}
//        
//        void compute(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr);
//    };
//    
//    class NnzGwas : public Parameter {
//    public:
//        
//        NnzGwas(): Parameter("NnzGwas"){}
//        
//        void compute(const VectorXf &effects, const vector<SparseVector<float> > &ZPZ, const VectorXf &ZPZdiag);
//    };
//    
//    class PiGwas : public Parameter {
//    public:
//        PiGwas(): Parameter("PiGwas"){}
//        
//        void compute(const float nnzGwas, const unsigned numSnps);
//    };
    
    class NumBadSnps : public Parameter {
    public:
        float betaThresh;
        vector<string> snpNames;
        vector<string> badSnpName;
        vector<unsigned> badSnpIdx;
        ofstream out;
        
        bool writeTxt;
        
        NumBadSnps(const string &title, const VectorXf &b, const vector<string> &snpNames): Parameter("NumSkeptSnp"), snpNames(snpNames){
            VectorXf abs_b = b.array().abs();
            std::sort(abs_b.data(), abs_b.data() + abs_b.size());
            int index8 = 0.8 * (abs_b.size() - 1);
            betaThresh = abs_b[index8];
            //cout << "Set beta cutoff threshold: " << betaThresh << endl;
            //cout << b.head(10) << endl;
            
            string filename = title + ".skepticalSNPs";
            out.open(filename.c_str());
            writeTxt = true;
        }
        
        void compute_sparse(VectorXi &badSnps, VectorXf &effects, VectorXf &effectMean, const VectorXf &b, VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const vector<ChromInfo*> &chromInfoVec, const int iter, const Data &data);
        void compute_full(VectorXi &badSnps, VectorXf &effects, VectorXf &effectMean, const VectorXf &b, VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec, const int iter, const Data &data);
        void compute_eigen(VectorXi &badSnps, VectorXf &effects, VectorXf &effectMean, const VectorXf &b, vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, const vector<LDBlockInfo*> keptLdBlockInfoVec, const int iter, const Data &data);
    };
    

    VectorXf rcorr;
        
    bool sparse;
    bool robustMode;
    bool noscale;
    bool lowRankModel;

    SnpEffects snpEffects;
    VarEffects sigmaSq;
    ResidualVar vare;
    GenotypicVar varg;
    Rounding rounding;
    NumBadSnps nBadSnps;
    BlockGenotypicVar vargBlk;
    BlockResidualVar vareBlk;
    
    vector<VectorXf> wcorrBlocks;
    vector<VectorXf> whatBlocks;

    /** Used only for --fixed-effect SNP Gibbs (distinct from BayesC::fixedEffects intercept/covariates). */
    FixedEffects fixedSnpEffects;

    HeritabilityFixedEffects hsqFixed;

    ApproxBayesC(const Data &data, const bool lowRank, const float varGenotypic, const float varResidual, const float varRandom, const float pival, const float piAlpha, const float piBeta, const bool estimatePi, const bool noscale, const bool robustMode, const bool message = true)
    : BayesC(data, varGenotypic, varResidual, 0.0, pival, piAlpha, piBeta, estimatePi, noscale, "Gibbs", false)
    , rcorr(data.ZPy)
    , wcorrBlocks(data.wcorrBlocks)
    , fixedSnpEffects(data.fixedEffectNames)
    , hsqFixed("hsqFixed")
    , snpEffects(data.snpEffectNames)
    , sigmaSq(varGenotypic, data.snp2pq, pival, noscale)
    , vare(varResidual, data.numKeptInds)
    , varg(varGenotypic, data.numKeptInds)
    , vargBlk(data.ldblockNames, varGenotypic, data.numKeptInds)
    , vareBlk(data.ldblockNames, data.varPhenotypic)
    , nBadSnps(data.title, data.b, data.snpEffectNames)
    , noscale(noscale)
    , sparse(data.sparseLDM)
    , robustMode(robustMode)
    , lowRankModel(lowRank)
    {
        
        paramSetVec = {&snpEffects, &snpPip};
        paramVec = {&pi, &nnzSnp, &sigmaSq, &hsq, &vare};
        paramToPrint = {&pi, &nnzSnp, &sigmaSq, &hsq, &vare};
        if (data.numFixedEffectSnps) {
            paramVec.insert(paramVec.begin() + 4, &hsqFixed);
            paramToPrint.insert(paramToPrint.begin() + 4, &hsqFixed);
        }

        if (lowRankModel) {
            paramSetVec.push_back(&vargBlk);
            paramSetVec.push_back(&vareBlk);
            paramToPrint.push_back(&nBadSnps);
        }

        if (message) {
            cout << "\nSBayesC" << endl;
            if (lowRankModel) {
                cout << "Using the low-rank model" << endl;
            }
            cout << "sigmaSq: " << sigmaSq.value << " scale factor: " << sigmaSq.scale << endl;
            if (noscale){
               cout << "Fitting model assuming unscaled genotypes " << endl; 
            } else {
               cout << "Fitting model assuming scaled genotypes "  << endl;
            }
            if (robustMode) cout << "Using a more robust parameterisation " << endl;
            if (data.numFixedEffectSnps)
                cout << "Fitting " << data.numFixedEffectSnps << " SNP(s) as fixed effects (flat prior)." << endl;
        }
    }
    
    void sampleUnknowns(const unsigned iter);
    static void ldScoreReg(const VectorXf &chisq, const VectorXf &LDscore, const VectorXf &LDsamplVar,
                           const float varg, const float vare, float &ps);
    void checkHsq(vector<float> &hsqMCMC);
};


class ApproxBayesB : public ApproxBayesC {
public:

    class SnpEffects : public ApproxBayesC::SnpEffects {
    public:
        VectorXf betaSq;     // save sample squres of full conditional normal distribution regardless of delta values
        
        SnpEffects(const vector<string> &header): ApproxBayesC::SnpEffects(header){
            betaSq.setZero(size);
        }
        
        void sampleFromFC_sparse(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const vector<ChromInfo*> &chromInfoVec,
                        const VectorXf &snp2pq,
                          const VectorXf &sigmaSq, const float pi, const float vare, const float varg);
        void sampleFromFC_full(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXf &snp2pq,
                          const VectorXf &sigmaSq, const float pi, const float vare, const float varg);
        void sampleFromFC_eigen(vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
                          const vector<LDBlockInfo*> keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks,
                          const VectorXf &sigmaSq, const float pi, const float varg, const VectorXf &snp2pq);
    };
    
    SnpEffects snpEffects;
    BayesB::VarEffects sigmaSq;

    ApproxBayesB(const Data &data, const float lowRank, const float varGenotypic, const float varResidual, const float pival, const float piAlpha, const float piBeta, const bool estimatePi, const bool noscale, const bool message = true)
    : ApproxBayesC(data, lowRank, varGenotypic, varResidual, pival, piAlpha, piBeta, estimatePi, noscale, false, false),
    snpEffects(data.snpEffectNames),
    sigmaSq(varGenotypic, data.snp2pq, pival, noscale)
    {
        paramSetVec = {&snpEffects, &snpPip};
        paramVec = {&pi, &nnzSnp, &hsq, &vare};
        paramToPrint = {&pi, &nnzSnp, &hsq, &vare};
        
        if (lowRankModel) {
            paramSetVec.push_back(&vargBlk);
            paramSetVec.push_back(&vareBlk);
            paramToPrint.push_back(&nBadSnps);
        }
        
        if (message) {
            cout << "\nSBayesB" << endl;
            if (lowRankModel) {
                cout << "Using the low-rank model" << endl;
            }
            cout << "scale factor: " << sigmaSq.scale << endl;
            if (noscale)
            {
                cout << "Fitting model assuming unscaled genotypes " << endl;
            } else
            {
                cout << "Fitting model assuming scaled genotypes "  << endl;
            }
        }
    }
    
    void sampleUnknowns(const unsigned iter);

};


class ApproxBayesS : public BayesS {
public:
    
    class SnpEffects : public ApproxBayesC::SnpEffects {
    public:
        float wtdSumSq;  // weighted sum of squares by 2pq^S
        float sum2pqSplusOne;  // sum of delta_j* (2p_j q_j)^{1+S}
        
        SnpEffects(const vector<string> &header, const VectorXf &snp2pq, const float pi): ApproxBayesC::SnpEffects(header) {
            wtdSumSq = 0.0;
            sum2pqSplusOne = snp2pq.sum()*pi;  // starting value of S is 0
        }
        
        void sampleFromFC_sparse(VectorXf &rcorr,const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const vector<ChromInfo*> &chromInfoVec,
                          const float sigmaSq, const float pi, const float vare,
                          const VectorXf &snp2pqPowS, const VectorXf &snp2pq, const float varg);
        void sampleFromFC_full(VectorXf &rcorr,const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const float sigmaSq, const float pi, const float vare,
                          const VectorXf &snp2pqPowS, const VectorXf &snp2pq, const float varg);

        void sampleFromFC_ind(const VectorXf &ZPy,const MatrixXf &Z, const VectorXf &ZPZdiag,
                          const float sigmaSq, const float pi, const float vare,
                          const VectorXf &snp2pqPowS, const VectorXf &snp2pq, VectorXf &ghat);

        void hmcSampler(VectorXf &rcorr, const VectorXf &ZPy, const vector<VectorXf> &ZPZ,
                        const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                        const float sigmaSq, const float pi, const float vare, const VectorXf &snp2pqPowS);
        VectorXf gradientU(const VectorXf &effects, VectorXf &rcorr, const VectorXf &ZPy, const vector<VectorXf> &ZPZ,
                           const VectorXi &windStart, const VectorXi &windSize, const unsigned chrStart, const unsigned chrSize,
                           const float sigmaSq, const float vare, const VectorXf &snp2pqPowS);
        float computeU(const VectorXf &effects, const VectorXf &rcorr, const VectorXf &ZPy,
                       const float sigmaSq, const float vare, const VectorXf &snp2pqPowS);
        
        void sampleFromFC_eigen(vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
                          const vector<LDBlockInfo*> keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks,
                          const float sigmaSq, const float pi, const float varg,
                          const VectorXf &snp2pqPowS, const VectorXf &snp2pq);

    };

    class MeanEffects : public Parameter, public Stat::Normal {
    public:
        VectorXf snp2pqPowSmu;
        
        MeanEffects(const unsigned numSnps, const string &lab = "Mu"): Parameter(lab){
            snp2pqPowSmu.setOnes(numSnps);
        }
        
        void sampleFromFC(const vector<SparseVector<float> > &ZPZ, const VectorXf &snpEffects, const VectorXf &snp2pq, const float vare, VectorXf &rcorr);
        
        void sampleFromFC(const VectorXf &snpEffects, const VectorXf &snp2pq);
    };
    
    class Smu : public Parameter, public Stat::Normal {
    public:
        float varProp;

        AcceptanceRate ar;
        Parameter tuner;

        Smu(const string &lab = "Smu"): Parameter(lab), tuner("varProp"){
            varProp = 0.01;
        }
    
        void sampleFromFC(const vector<SparseVector<float> > &ZPZ, const VectorXf &snpEffects, const VectorXf &snp2pq, const float vare, VectorXf &snp2pqPowSmu, VectorXf &rcorr);
    };
    

public:
    VectorXf rcorr;
        
    bool sparse;
    bool lowRankModel;
    bool estimateEffectMean;

    SnpEffects snpEffects;
    ApproxBayesC::ResidualVar vare;
    ApproxBayesC::GenotypicVar varg;
    ApproxBayesC::Rounding rounding;
    ApproxBayesC::NumBadSnps nBadSnps;
    ApproxBayesC::BlockGenotypicVar vargBlk;
    ApproxBayesC::BlockResidualVar vareBlk;
    
    vector<VectorXf> wcorrBlocks;
    vector<VectorXf> whatBlocks;

    MeanEffects mu;
    Smu Su;

    ApproxBayesS(const Data &data, const bool lowRank, const float varGenotypic, const float varResidual, const float pival, const float piAlpha, const float piBeta, const bool estimatePi,
                 const float varS, const vector<float> &svalue,
                 const string &algorithm, const bool noscale = false, const bool message = true)
    : BayesS(data, varGenotypic, varResidual, 0.0, pival, piAlpha, piBeta, estimatePi, varS, svalue, algorithm, !(!noscale || lowRank), false)
    , rcorr(data.ZPy)
    , wcorrBlocks(data.wcorrBlocks)
    , snpEffects(data.snpEffectNames, data.snp2pq, pival)
    , vare(varResidual, data.numKeptInds)
    , varg(varGenotypic, data.numKeptInds)
    , mu(data.numIncdSnps)
    , vargBlk(data.ldblockNames, varGenotypic, data.numKeptInds)
    , vareBlk(data.ldblockNames, data.varPhenotypic)
    , nBadSnps(data.title, data.b, data.snpEffectNames)
    , sparse(data.sparseLDM)
    , lowRankModel(lowRank)
   {

        estimateEffectMean = false;
    
        // Override paramVec to use ApproxBayesS's own varg and vare members
        paramSetVec = {&snpEffects, &fixedEffects, &snpPip};
        paramVec = {&pi, &nnzSnp, &sigmaSq, &S, &varg, &vare, &hsq};
        paramToPrint = {&pi, &nnzSnp, &sigmaSq, &scale, &S, &varg, &vare, &hsq, &S.ar, &S.tuner};
        if (data.numRandomEffects) {
            paramSetVec.push_back(&randomEffects);
            paramVec.push_back(&sigmaSqRand);
            paramVec.push_back(&varRand);
            paramToPrint.push_back(&varRand);
        }

       if (message) {
            string alg = algorithm;
            if (alg!="RWMH" && alg!="Reg") alg = "HMC";
            cout << "\nSBayesS" << endl;
            if (lowRankModel) {
                cout << "Using the low-rank model" << endl;
            }
            if (scaledGeno) {
                cout << "Fitting model assuming scaled genotypes " << endl;
            } else {
                cout << "Fitting model assuming unscaled genotypes " << endl;
            }
            cout << "Algorithm: " << alg << "." << endl;
            cout << "scale factor: " << sigmaSq.scale << endl;
        }

        //if (randomStart) sampleStartVal();
    }
    
    void sampleUnknowns(const unsigned iter);
};


class ApproxBayesST : public ApproxBayesS {
    // Approximate BayesST is a model to account for both MAF- and LD-dependent architecture
public:
    
    class SnpEffects : public ApproxBayesS::SnpEffects {
    public:
        float wtdSumSq;  // weighted sum of squares by 2pq^S * ldsc^T
        float sum2pqhSlT;  // sum of delta_j* (2p_j q_j)* h_j^S * l_j^T
        
        SnpEffects(const vector<string> &header, const VectorXf &snp2pq, const float pi): ApproxBayesS::SnpEffects(header, snp2pq, pi) {
            wtdSumSq = 0.0;
            sum2pqhSlT = snp2pq.sum()*pi;  // starting value of S,T is 0
        }
        
        void sampleFromFC(VectorXf &rcorr,const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const vector<ChromInfo*> &chromInfoVec, const ArrayXf &hSlT, const VectorXf &snp2pq,
                          const float sigmaSq, const float pi, const float vare, const float varg);
    };
    
    class Sp : public BayesS::Sp {
    public:
        Sp(const unsigned m): BayesS::Sp(m, 1, 0, "HMC", "S"){}
        
        //sample S and T jointly
        void sampleFromFC(const unsigned numNonZeros, const float sigmaSq, const VectorXf &snpEffects,
                          const VectorXf &snp2pq, const ArrayXf &logSnp2pq,
                          const VectorXf &ldsc, const ArrayXf &logLdsc,
                          const float varg, float &scale, float &T, ArrayXf &hSlT);
        Vector2f gradientU(const Vector2f &ST, const ArrayXf &snpEffectSq, const float sigmaSq,
                           const float snp2pqLogSum, const ArrayXf &snp2pq, const ArrayXf &logSnp2pq,
                           const float ldscLogSum, const ArrayXf &ldsc, const ArrayXf &logLdsc);
        float computeU(const Vector2f &ST, const ArrayXf &snpEffectSq, const float sigmaSq,
                       const float snp2pqLogSum, const ArrayXf &snp2pq, const ArrayXf &logSnp2pq,
                       const float ldscLogSum, const ArrayXf &ldsc, const ArrayXf &logLdsc);
        
    };

    class Tp : public BayesS::Sp {
    public:
        Tp(const unsigned m): BayesS::Sp(m, 1, 0, "HMC", "T"){}
        
        void sampleFromFC(const unsigned numNonZeros, const float sigmaSq, const VectorXf &snpEffects,
                          const VectorXf &snp2pq, const VectorXf &ldsc, const ArrayXf &logLdsc, const float varg, float &scale, ArrayXf &hSlT);
        float gradientU(const float &T, const ArrayXf &snpEffectSq, const float sigmaSq,
                        const float ldscLogSum, const ArrayXf &ldsc, const ArrayXf &logLdsc);
        float computeU(const float &T, const ArrayXf &snpEffectSq, const float sigmaSq,
                       const float ldscLogSum, const ArrayXf &ldsc, const ArrayXf &logLdsc);

        
    };
    
    const bool estimateS;
    const ArrayXf logLdsc;
    ArrayXf hSlT;
    
    SnpEffects snpEffects;
    Sp S;
    Tp T;
    
    ApproxBayesST(const Data &data, const bool lowRank, const float varGenotypic, const float varResidual, const float pival,
                  const float piAlpha, const float piBeta, const bool estimatePi, const float varS, const vector<float> &svalue, const bool estimateS, const bool noscale = false, const bool message = true):
    ApproxBayesS(data, lowRank, varGenotypic, varResidual, pival, piAlpha, piBeta, estimatePi, varS, svalue, "HMC", noscale, false),
    estimateS(estimateS),
    logLdsc(data.LDscore.array().log()),
    hSlT(snp2pqPowS),
    snpEffects(data.snpEffectNames, data.snp2pq, pival),
    S(data.numIncdSnps),
    T(data.numIncdSnps)
    {
        paramSetVec = {&snpEffects, &snpPip};
        paramVec = {&pi, &nnzSnp, &sigmaSq, &S, &T, &varg, &vare, &hsq};
        paramToPrint = {&pi, &nnzSnp, &sigmaSq, &S, &T, &varg, &vare, &hsq};
        if (lowRankModel) {
            paramSetVec.push_back(&vargBlk);
            paramSetVec.push_back(&vareBlk);
        }
        if (message) {
            cout << "\nSBayesST" << endl;
            if (lowRankModel) {
                cout << "Using the low-rank model" << endl;
            }
            if (scaledGeno) {
                cout << "Fitting model assuming scaled genotypes " << endl;
            } else {
                cout << "Fitting model assuming unscaled genotypes " << endl;
            }
        }
        
        //if (randomStart) sampleStartVal();
    }

    void sampleUnknowns(const unsigned iter);
    void sampleStartVal(void);
};


// -----------------------------------------------------------------------------------------------
// Approximate Bayes R
// -----------------------------------------------------------------------------------------------

class ApproxBayesR : public BayesR {
    
public:
        
    class SnpEffects : public ApproxBayesC::SnpEffects {
    public:
        vector<vector<unsigned> > snpset;
        VectorXi membership;
        VectorXi deltaNzIdx;
        VectorXf deltaNZ;
        VectorXf lambdaVec;
        VectorXf uhatVec;
        VectorXf invGammaVec;
        float sum2pq;
        float wtdSumSq;  // weighted inverse of gamma


        SnpEffects(const vector<string> &header): ApproxBayesC::SnpEffects(header){
            sum2pq = 0.0;
            wtdSumSq = 0.0;
            deltaNZ.setZero(size);
        }
        
        void sampleFromFC_sparse(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const vector<ChromInfo*> &chromInfoVec,
                          const VectorXf &snp2pq,
                          const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare, VectorXf &snpStore, 
                          const float varg,
                          const bool hsqPercModel, DeltaPi &deltaPi,
                          const vector<SnpInfo*> &incdSnpInfoVec);
        void sampleFromFC_full(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXf &snp2pq,
                          const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare, VectorXf &snpStore,
                          const float varg,
                          const bool hsqPercModel, DeltaPi &deltaPi,
                          const vector<SnpInfo*> *incdSnpInfoVec = nullptr);
        
        void sampleFromFC_eigen(vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
                          const vector<LDBlockInfo*> &keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks,
                          const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, VectorXf &snpStore, const float varg,
                          const bool hsqPercModel, DeltaPi &deltaPi);
        
        // tempered Gibbs sampler
//        void sampleFromTGS_eigen(const vector<vector<int> > &selectedSnps, vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
//                                 const vector<LDBlockInfo*> &keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks,
//                                 const VectorXf &pis, const VectorXf &gamma, const float varg, const bool hsqPercModel, const float sigmaSq);
        

        // obsoleted
        void sampleFromFC(const VectorXf &ZPy, const SpMat &ZPZsp, const VectorXf &ZPZdiag,
                          VectorXf &rcorr, const VectorXf &LDsamplVar,
                          const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, VectorXf &snpStore,
                          const float varg, const float vare, const float ps, const float overdispersion, const bool hsqPercModel, DeltaPi &deltaPi);
        
        void sampleFromFC(const VectorXf &ZPy, const VectorXf &ZPZdiag, const MatrixXf &Z, const float n_ref, const float n_gwas,
                          const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare,
                          VectorXf &snpStore, VectorXf &ghat, const float varg, const bool hsqPercModel, DeltaPi &deltaPi);

        void adjustByCG(const VectorXf &ZPy, const vector<SparseVector<float> > &ZPZsp, VectorXf &rcorr);
        
        void sampleFromPrior(const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float varg, const bool hsqPercModel);
    };
    
    class VarEffects : public BayesR::VarEffects {
    public:
        VarEffects(const float varg, const VectorXf &snp2pq, const VectorXf &gamma, const VectorXf &pi, const bool noscale):
        BayesR::VarEffects(varg, snp2pq, gamma, pi, noscale){}
        
        void computeRobustMode(const float varg, const VectorXf &snp2pq, const VectorXf &gamma, const VectorXf &pi, const bool noscale);
    };
    
    class VgMixComps : public BayesR::VgMixComps {
    public:
        VgMixComps(const VectorXf &gamma): BayesR::VgMixComps(gamma){}
        
        void compute(const VectorXf &snpEffects, const VectorXf &ZPy, const VectorXf &rcorr, const vector<vector<unsigned> > &snpset, const float varg, const float nobs);
        void compute(const VectorXf &snpEffects, const vector<vector<unsigned> > &snpset);
        void compute(const VectorXf &snpEffects, const vector<unsigned> &membership, const vector<MatrixXf> &Qblocks, const vector<LDBlockInfo*> &keptLdBlockInfoVec);
    };
    
    
    VectorXf rcorr;

    bool sparse;
    bool robustMode;
    bool noscale;
    bool lowRankModel;

    SnpEffects snpEffects;
    VarEffects sigmaSq;
    VgMixComps Vgs;
    ApproxBayesC::ResidualVar vare;
    ApproxBayesC::GenotypicVar varg;
    ApproxBayesC::Rounding rounding;
    ApproxBayesC::NumBadSnps nBadSnps;
    ApproxBayesC::BlockGenotypicVar vargBlk;
    ApproxBayesC::BlockResidualVar vareBlk;
    BayesC::SnpHsqPEP snpHsqPep;

    vector<VectorXf> wcorrBlocks;
    vector<VectorXf> whatBlocks;

    ApproxBayesC::FixedEffects fixedSnpEffects;
    BayesC::HeritabilityFixedEffects hsqFixed;
    
    enum {gibbs, cg, mh, tgs, tgs_thin} algorithm;
    
    vector<vector<int> > highLDsnpSet;
        
    ApproxBayesR(const Data &data, const bool lowRank, const float varGenotypic, const float varResidual, const VectorXf pis, const VectorXf &piPar, const VectorXf gamma, const bool estimatePi, const bool noscale, const bool hsqPercModel, const bool robustMode, const string &alg, const bool message = true):
    BayesR(data, varGenotypic, varResidual, 0.0, pis, piPar, gamma, estimatePi, noscale, hsqPercModel, alg, false)
    , rcorr(data.ZPy)
    , snpEffects(data.snpEffectNames)
    , sigmaSq(varGenotypic, data.snp2pq, gamma, pis, noscale)
    , Vgs(gamma)
    , vare(varResidual, data.numKeptInds)
    , varg(varGenotypic, data.numKeptInds)
    , vargBlk(data.ldblockNames, varGenotypic, data.numKeptInds)
    , vareBlk(data.ldblockNames, data.varPhenotypic)
    , nBadSnps(data.title, data.b, data.snpEffectNames)
    , snpHsqPep(data.snpEffectNames)
    , wcorrBlocks(data.wcorrBlocks)
    , fixedSnpEffects(data.fixedEffectNames)
    , hsqFixed("hsqFixed")
    , noscale(noscale)
    , sparse(data.sparseLDM)
    , robustMode(robustMode)
    , lowRankModel(lowRank)
    {

        if (alg == "cg") algorithm = cg;
        else if (alg == "MH") algorithm = mh;
        else if (alg == "TGS") algorithm = tgs;
        else if (alg == "TGS_thin") algorithm = tgs_thin;
        else algorithm = gibbs;

        paramSetVec = {&snpEffects, &snpPip, &snpHsqPep};
        paramSetVec.insert(paramSetVec.end(), deltaPi.begin(), deltaPi.end());
        paramVec     = {&nnzSnp, &sigmaSq, &hsq, &vare};
        if (data.numFixedEffectSnps)
            paramVec.insert(paramVec.begin() + 3, &hsqFixed);
        paramVec.insert(paramVec.end(), numSnps.begin(), numSnps.end());
        paramVec.insert(paramVec.end(), Vgs.begin(), Vgs.end());
        
        paramToPrint = {&sigmaSq, &hsq, &vare};
        if (data.numFixedEffectSnps)
            paramToPrint.insert(paramToPrint.begin() + 2, &hsqFixed);
        paramToPrint.insert(paramToPrint.begin(), Vgs.begin(), Vgs.end());
        paramToPrint.insert(paramToPrint.begin(), numSnps.begin(), numSnps.end());
        
       if (lowRankModel) {
            paramSetVec.push_back(&vargBlk);
            paramSetVec.push_back(&vareBlk);
            paramToPrint.push_back(&nBadSnps);
        }
                        
        if (message) {
            cout << "\nSBayesR" << endl;
            if (lowRankModel) {
                cout << "Using the low-rank model" << endl;
            }
//            cout << "scale factor: " << sigmaSq.scale << endl;
            cout << "Gamma: " << gamma.transpose() << endl;
            if (noscale){
               cout << "Fitting model assuming unscaled genotypes " << endl; 
            } else {
               cout << "Fitting model assuming scaled genotypes "  << endl;
            }
            if (robustMode) cout << "Using a more robust parameterisation " << endl;
            if (algorithm == cg) cout << "Conjugate gradient-adjusted Gibbs sampling" << endl;
            if (algorithm == tgs_thin) cout << "Using tempered Gibbs sampling (TGS)" << endl;
            if (!hsqPercModel) cout << "The SNP effect prior is a mixture distribution with an unknown variance variable." << endl;
            if (data.numFixedEffectSnps)
                cout << "Fitting " << data.numFixedEffectSnps << " SNP(s) as fixed effects (flat prior)." << endl;
        }
    }
    
    void sampleUnknowns(const unsigned iter);
    void updateRHSfull(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec, const VectorXf &snpEffects);
   void updateRHSsparse(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZ, const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec, const VectorXf &snpEffects);
    void updateRHSlowRankModel(vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, const vector<LDBlockInfo*> &keptLdBlockInfoVec, const VectorXf &snpEffects);
};

// -----------------------------------------------------------------------------------------------
// Approximate Bayes RS
// -----------------------------------------------------------------------------------------------

class ApproxBayesRS : public ApproxBayesR {
public:
    
    class SnpEffects : public ApproxBayesS::SnpEffects {
    public:
        unsigned ndist;
        ArrayXf numSnpMix;
        vector<vector<unsigned> > snpset;
        
        SnpEffects(const vector<string> &header, const VectorXf &snp2pq, const VectorXf &pis): ApproxBayesS::SnpEffects(header, snp2pq, 1.0-pis[0]) {
            ndist = pis.size();
            numSnpMix.setZero(ndist);
        }
        
        void sampleFromFC_sparse(VectorXf &rcorr,const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const vector<ChromInfo*> &chromInfoVec,
                          const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare,
                          const ArrayXf &snp2pqPowS, const VectorXf &snp2pq,
                          const float varg, const bool hsqPercModel);
        
        void sampleFromFC_full(VectorXf &rcorr,const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare,
                          const ArrayXf &snp2pqPowS, const VectorXf &snp2pq,
                          const float varg, const bool hsqPercModel);
    };
        
    class Sp : public BayesS::Sp {
    public:
    
        Sp(const unsigned m, const float var, const float start): BayesS::Sp(m, var, start, "HMC"){}
        
        void sampleFromFC(vector<vector<unsigned> > &snpset, const VectorXf &snpEffects,
                          float &sigmaSq, const VectorXf &gamma,
                          const VectorXf &snp2pq, ArrayXf &snp2pqPowS, const ArrayXf &logSnp2pq,
                          const float vg, float &scale, float &sum2pqSplusOne, const bool hsqPercModel);
        float gradientU(const float S, const unsigned nnzMix, const vector<ArrayXf> &snpEffectMix, const float snp2pqLogSum, const vector<ArrayXf> &snp2pqMix, const vector<ArrayXf> &logSnp2pqMix, const float sigmaSq, const VectorXf &gamma, const float vg);
        float computeU(const float S, const unsigned nnzMix, const vector<ArrayXf> &snpEffectMix, const float snp2pqLogSum, const vector<ArrayXf> &snp2pqMix, const vector<ArrayXf> &logSnp2pqMix, const float sigmaSq, const VectorXf &gamma, const float vg, float &scale, const bool hsqPercModel);

    };
    
    float genVarPrior;
    float scalePrior;

    SnpEffects snpEffects;
    Sp S;
    
    ArrayXf logSnp2pq;
    ArrayXf snp2pqPowS;
    
    ApproxBayesRS(const Data &data, const bool lowRank, const float varGenotypic, const float varResidual, const VectorXf pis, const VectorXf &piPar, const VectorXf gamma, const bool estimatePi, const float varS, const vector<float> &svalue, const bool noscale, const bool hsqPercModel, const bool robustMode, const string &alg, const bool message = true):
    ApproxBayesR(data, lowRank, varGenotypic, varResidual, pis, piPar, gamma, estimatePi, noscale, hsqPercModel, robustMode, alg, false),
    snpEffects(data.snpEffectNames, data.snp2pq, pis),
    S(data.numIncdSnps, varS, svalue[0]),
    genVarPrior(varGenotypic),
    scalePrior(sigmaSq.scale)
    {
        logSnp2pq = data.snp2pq.array().log();
        snp2pqPowS = data.snp2pq.array().pow(S.value);
                
        sigmaSq.value = varGenotypic/((snp2pqPowS*data.snp2pq.array()).sum()*(1.0-pis[0]));
        scale.value = sigmaSq.scale = 0.5*sigmaSq.value;

        paramSetVec = {&snpEffects, &snpPip};
        paramVec     = {&nnzSnp, &sigmaSq, &S, &hsq, &vare};
        paramVec.insert(paramVec.end(), numSnps.begin(), numSnps.end());
        paramVec.insert(paramVec.end(), Vgs.begin(), Vgs.end());
        
        paramToPrint = {&sigmaSq, &S, &hsq, &vare};
        paramToPrint.insert(paramToPrint.begin(), Vgs.begin(), Vgs.end());
        paramToPrint.insert(paramToPrint.begin(), numSnps.begin(), numSnps.end());
        
        if (message) {
            cout << "\nApproximate BayesRS model fitted." << endl;
            cout << "scale factor: " << sigmaSq.scale << endl;
            cout << "Gamma: " << gamma.transpose() << endl;
            if (noscale){
                cout << "Fitting model assuming unscaled genotypes " << endl;
            } else {
                cout << "Fitting model assuming scaled genotypes "  << endl;
            }
            if (algorithm == tgs || algorithm == tgs_thin) {
                cout << "Using tempered Gibbs sampling (TGS) algorithm for high-LD SNPs" << endl;
                if (!data.LDmap.size()) {
                    cout << "\nError: To use tempered GIbbs sampling, you need to give pairwise LD file by --pwld-file " << endl;
                }
            }
        }
    }

    void sampleUnknowns(const unsigned iter);
};


class ApproxBayesSMix : public ApproxBayesS {
    // Approximate BayesSMix is a mixture model of zero (component 1), BayesC (component 2) and BayesS (component 3) prior
public:
    
    class SnpEffects : public ApproxBayesS::SnpEffects {
    public:
        Vector2f wtdSum2pq;
        Vector2f wtdSumSq;
        Vector2f numNonZeros;
        Vector3f numSnpMixComp;
        VectorXf valuesMixCompS;
        
        SnpEffects(const vector<string> &header, const VectorXf &snp2pq, const float pi): ApproxBayesS::SnpEffects(header, snp2pq, pi) {}
        
        void sampleFromFC_sparse(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const vector<ChromInfo*> &chromInfoVec, const ArrayXf &snp2pqPowS, const VectorXf &snp2pq,
                          const Vector2f &sigmaSq, const Vector3f &pi, const float vare, const float varg,
                          VectorXf &deltaS);
        
        void sampleFromFC_full(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec, const ArrayXf &snp2pqPowS, const VectorXf &snp2pq,
                          const Vector2f &sigmaSq, const Vector3f &pi, const float vare, const float varg,
                           VectorXf &deltaS);

    };
    
    class DeltaS : public ParamSet {  // indicator variable for S model component
    public:
        DeltaS(const vector<string> &header): ParamSet("DeltaS", header){};
    };
    
    class PiMixComp : public vector<Parameter*>, public Stat::Dirichlet {
    public:
        const unsigned ndist;
        Vector3f alpha;
        Vector3f values;
        
        PiMixComp(const float pival): ndist(3) {
            vector<string> label = {"0", "C", "S"};
            for (unsigned i=0; i<ndist; ++i) {
                this->push_back(new Parameter("Pi" + label[i]));
            }
            alpha.setOnes();
            values << 1.0-pival, 0.5*pival, 0.5*pival;
        }
        
        void sampleFromFC(const VectorXf &numSnpMixComp);
    };
    
    class VarEffects : public vector<BayesC::VarEffects*> {
    public:
        Vector2f values;
        
        VarEffects(const float vg, const float pi, const VectorXf &snp2pq) {
            vector<string> label = {"C", "S"};
            for (unsigned i=0; i<2; ++i) {
                this->push_back(new BayesC::VarEffects(vg, snp2pq, 0.5*pi, true, "SigmaSq" + label[i]));
                values[i] = (*this)[i]->value;
            }
        }
        
        void sampleFromFC(const Vector2f &snpEffSumSq, const Vector2f &numSnpEff);
        void computeScale(const Vector2f &varg, const Vector2f &wtdSum2pq);
    };
    
    
    class GenotypicVarMixComp : public vector<Parameter*> {
    public:
        Vector2f values;
        
        GenotypicVarMixComp() {
            vector<string> label = {"C", "S"};
            for (unsigned i=0; i<2; ++i) {
                this->push_back(new Parameter("GenVar" + label[i]));
            }
        }
        
        void compute(const Vector2f &sigmaSq, const Vector2f &wtdSum2pq);

    };
    
    class HeritabilityMixComp : public vector<Parameter*> {
    public:
        Vector2f values;
        
        HeritabilityMixComp() {
            vector<string> label = {"C", "S"};
            for (unsigned i=0; i<2; ++i) {
                this->push_back(new Parameter("hsq" + label[i]));
            }
        }
        
        void compute(const Vector2f &vargMixComp, const float varg, const float vare);
    };
    
    SnpEffects snpEffects;
    DeltaS deltaS;
    PiMixComp piMixComp;
    VarEffects sigmaSq;
    GenotypicVarMixComp vargMixComp;
    HeritabilityMixComp hsqMixComp;
    
    ApproxBayesSMix(const Data &data, const bool lowRank, const float varGenotypic, const float varResidual, const float pival, const float varS, const vector<float> &svalue,
                  const bool noscale = false, const bool message = true):
    ApproxBayesS(data, lowRank, varGenotypic, varResidual, pival, 1, 1, true, varS, svalue, "HMC", noscale, false),
    snpEffects(data.snpEffectNames, data.snp2pq, 0.5*pival),
    deltaS(data.snpEffectNames),
    piMixComp(pival),
    sigmaSq(varGenotypic, pival, data.snp2pq)
    {
        paramSetVec = {&snpEffects, &snpPip, &deltaS};
        paramVec = {piMixComp[2], piMixComp[1], &pi, &nnzSnp, sigmaSq[1], &S, sigmaSq[0], hsqMixComp[1], hsqMixComp[0], &hsq};
        paramToPrint = {piMixComp[2], piMixComp[1], &pi, &nnzSnp, sigmaSq[1], &S, sigmaSq[0], hsqMixComp[1], hsqMixComp[0], &hsq};
        if (lowRankModel) {
            paramSetVec.push_back(&vargBlk);
            paramSetVec.push_back(&vareBlk);
        }
        if (message) {
            cout << "\nSBayesSMix" << endl;
            if (lowRankModel) {
                cout << "Using the low-rank model" << endl;
            }
            if (scaledGeno) {
                cout << "Fitting model assuming scaled genotypes " << endl;
            } else {
                cout << "Fitting model assuming unscaled genotypes " << endl;
            }
        }
    }
    
    void sampleUnknowns(const unsigned iter);

};

class BayesSMix : public BayesS {
    // BayesSMix is a mixture model of zero (component 1), BayesC (component 2) and BayesS (component 3) prior
public:

    class SnpEffects : public BayesS::SnpEffects {
    public:
        Vector2f wtdSum2pq;
        Vector2f wtdSumSq;
        Vector2f numNonZeros;
        Vector3f numSnpMixComp;
        VectorXf valuesMixCompS;
        
        SnpEffects(const vector<string> &header, const VectorXf &snp2pq, const float pi): BayesS::SnpEffects(header, snp2pq, pi) {}
        
        void sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const ArrayXf &snp2pqPowS, const VectorXf &snp2pq,
                          const Vector2f &sigmaSq, const Vector3f &pi, const float vare, VectorXf &deltaS, VectorXf &ghat, vector<VectorXf> &ghatMixComp);
    };
    
    class GenotypicVarMixComp : public ApproxBayesSMix::GenotypicVarMixComp {
    public:
        
        void compute(const vector<VectorXf> &ghatMixComp);
    };

    
    vector<VectorXf> ghatMixComp;
    
    SnpEffects snpEffects;
    GenotypicVarMixComp vargMixComp;
    ApproxBayesSMix::DeltaS deltaS;
    ApproxBayesSMix::PiMixComp piMixComp;
    ApproxBayesSMix::VarEffects sigmaSq;
    ApproxBayesSMix::HeritabilityMixComp hsqMixComp;

    BayesSMix(const Data &data, const float varGenotypic, const float varResidual, const float varRandom, const float pival, const float piAlpha, const float piBeta, const bool estimatePi, const float varS, const vector<float> &svalue, const string &algorithm, const bool message = true):
    BayesS(data, varGenotypic, varResidual, varRandom, pival, piAlpha, piBeta, estimatePi, varS, svalue, "HMC", false),
    snpEffects(data.snpEffectNames, data.snp2pq, 0.5*pival),
    deltaS(data.snpEffectNames),
    piMixComp(pival),
    sigmaSq(varGenotypic, pival, data.snp2pq)
    {
        ghatMixComp.resize(2);
        paramSetVec = {&snpEffects, &deltaS, &fixedEffects};
        paramVec = {piMixComp[2], piMixComp[1], &pi, &nnzSnp, sigmaSq[1], &S, sigmaSq[0], hsqMixComp[1], hsqMixComp[0], &hsq};
        paramToPrint = {piMixComp[2], piMixComp[1], &pi, &nnzSnp, sigmaSq[1], &S, sigmaSq[0], hsqMixComp[1], hsqMixComp[0], &hsq, &rounding};
        if (message) {
            cout << "\nBayesSMix model fitted." << endl;
        }
    }
    
    void sampleUnknowns(const unsigned iter);
};


// -----------------------------------------------------------------------------------------------
// Approximate Bayes RC: fitting functional annotations
// -----------------------------------------------------------------------------------------------

class ApproxBayesRC : public ApproxBayesR {
public:
        
    class SnpEffects : public ApproxBayesR::SnpEffects {
    public:
        unsigned ndist;
        ArrayXf numSnpMix;
        MatrixXf z;
        vector<vector<unsigned> > snpset;
        VectorXf fcMean;
        
        SnpEffects(const vector<string> &header, const VectorXf &pis): ApproxBayesR::SnpEffects(header){
            ndist = pis.size();
            numSnpMix.setZero(ndist);
            z.setZero(size, ndist-1);
            fcMean.setZero(size);
        }
        
        void sampleFromFC_sparse(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                                 const vector<ChromInfo*> &chromInfoVec,
                                 const float sigmaSq, const MatrixXf &snpPi, const VectorXf &gamma,
                                 const float vare, const float varg,
                                 const bool hsqPercModel, DeltaPi &deltaPi, const vector<SnpInfo*> &incdSnpInfoVec);
        void sampleFromFC_eigen(vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
                                const vector<LDBlockInfo*> &keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks,
                                const MatrixXf &snpPi, const VectorXf &gamma, const float varg,
                                DeltaPi &deltaPi, const bool hsqPercModel, const float sigmaSq, const vector<SnpInfo*> &incdSnpInfoVec);
        // tempered Gibbs sampler
        void sampleFromTGS_eigen(vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks, vector<VectorXf> &whatBlocks,
                                 const map<SnpInfo*, vector<SnpInfo*> > &LDmap, const vector<LDBlockInfo*> &keptLdBlockInfoVec, const VectorXf &nGWASblocks, const VectorXf &vareBlocks,
                                 const MatrixXf &snpPi, const VectorXf &gamma, const float varg,
                                 DeltaPi &deltaPi, const bool hsqPercModel, const float sigmaSq);
        
    };
    
    class AnnoEffects : public vector<BayesC::FixedEffects*>, public Stat::TruncatedNormal  {
    public:
        unsigned numComp;  // number of components = number of mixture components - 1
        unsigned numAnno;  // number of annotations
        MatrixXf wcorr;
        //VectorXf varwcorr;
        VectorXf annoDiag;
        VectorXf ssq;
        vector<string> colnames;

        VectorXf varProp;
        vector<BayesS::AcceptanceRate*> ar;

        AnnoEffects(const vector<string> &header, const unsigned ndist, const MatrixXf &annoMat) {
            numComp = ndist - 1;
            colnames.resize(numComp);
            ar.resize(numComp);
            varProp.setZero(numComp);
            numAnno = header.size();
            unsigned numSnps = annoMat.rows();
            for (unsigned i = 0; i<numComp; ++i) {
                colnames[i] = "AnnoEffects_p" + to_string(static_cast<long long>(i + 2));
                this->push_back(new BayesC::FixedEffects(header, colnames[i]));
                ar[i] = new BayesS::AcceptanceRate;
                varProp[i] = 0.01;
            }
            wcorr.setZero(numSnps, numComp);
            //varwcorr.setZero(numComp);
            annoDiag.setZero(numAnno);
            annoDiag[0] = numSnps;  // first annotation is intercept
            for (unsigned j=1; j<numAnno; ++j) {
                annoDiag[j] = annoMat.col(j).squaredNorm();
            }
            ssq.setZero(numComp);
        }
        
//        void sampleFromFC(MatrixXf &snpP, const MatrixXf &annoMat);
        void sampleFromFC_Gibbs(MatrixXf &z, const MatrixXf &annoMat, const VectorXf &sigmaSq, MatrixXf &snpP);
        void sampleFromFC_MH(MatrixXf &z, const MatrixXf &annoMat, const VectorXf &sigmaSq, MatrixXf &snpP);
        void initIntercept_probit(const VectorXf &pis);
        void initIntercept_logistic(const VectorXf &pis);
    };
    
    class VarAnnoEffects : public ParamSet, public Stat::InvChiSq {
    public:
        const float df;
        const float scale;
        unsigned numAnno;
        
        VarAnnoEffects(const vector<string> &header, const unsigned numAnno, const string &lab = "SigmaSqAnno"):
        ParamSet(lab, header), df(4), scale(1), numAnno(numAnno) {
            values.setOnes(size);
        }
        
        void sampleFromFC(const VectorXf &ssq);
    };

    class AnnoCondProb : public vector<ParamSet*>, public Stat::Normal {
    public:
        vector<string> colnames;
        unsigned numComp;
        unsigned numAnno;

        AnnoCondProb(const vector<string> &header, const unsigned numComp, const string &lab = "AnnoCondProb"):
        numComp(numComp) {
            colnames.resize(numComp);
            numAnno = header.size();
            for (unsigned i = 0; i<numComp; ++i) {
                colnames[i] = "AnnoCondProb_p" + to_string(static_cast<long long>(i + 2));
                this->push_back(new ParamSet(colnames[i], header));
            }
        }
        
        void compute_probit(const AnnoEffects &annoEffects, const vector<AnnoInfo*> &annoInfoVec);
        void compute_logistic(const AnnoEffects &annoEffects, const vector<AnnoInfo*> &annoInfoVec);
    };
    
    class AnnoJointProb : public vector<ParamSet*> {
    public:
        vector<string> colnames;
        unsigned numDist;

        AnnoJointProb(const vector<string> &header, const unsigned numDist, const string &lab = "AnnoJointProb"):
        numDist(numDist) {
            colnames.resize(numDist);
            for (unsigned i = 0; i<numDist; ++i) {
                colnames[i] = "AnnoJointProb_pi" + to_string(static_cast<long long>(i + 1));
                this->push_back(new ParamSet(colnames[i], header));
            }
        }
        
        void compute(const AnnoCondProb &annoCondProb);
    };
    
    class AnnoGenVar : public vector<ParamSet*> {
    public:
        vector<string> colnames;
        unsigned numComp;
        unsigned numAnno;
        float nobs;
        
        AnnoGenVar(const vector<string> &header, const unsigned numDist, const unsigned nobs, const string &lab = "AnnoGenVar"):
        numComp(numDist-1), nobs(nobs) {
            numAnno = header.size();
            colnames.resize(numComp);
            for (unsigned i = 0; i<numComp; ++i) {
                colnames[i] = "AnnoGenVar_pi" + to_string(static_cast<long long>(i + 2));
                this->push_back(new ParamSet(colnames[i], header));
            }
        }
        
        void compute(const VectorXf &snpEffects, const vector<vector<unsigned> > &snpset, const VectorXf &ZPy, const VectorXf &rcorr, const MatrixXf &annoMat);
        void compute(const VectorXf &snpEffects, const vector<unsigned> &membership, const vector<MatrixXf> &Qblocks, const vector<LDBlockInfo*> &keptLdBlockInfoVec, const MatrixXf &annoMat, const vector<AnnoInfo*> &annoInfoVec);
        void compute(const VectorXf &snpEffects, const vector<vector<unsigned> > &snpset, const MatrixXf &annoMat);
    };
    
    class AnnoTotalGenVar : public ParamSet {
    public:
        
        AnnoTotalGenVar(const vector<string> &header, const string &lab = "AnnoTotalGenVar"):
        ParamSet(lab, header) {}
        
        void compute(const AnnoGenVar &annoGenVar);
    };
    
    class AnnoPerSnpHsqEnrichment : public ParamSet {
    public:
        VectorXf invSnpProp;
        
        AnnoPerSnpHsqEnrichment(const vector<string> &header, const vector<AnnoInfo*> &annoVec, const string &lab = "Marginal_Heritability_Enrichment"):
        ParamSet(lab, header) {
            values.setOnes(size);
            invSnpProp.setZero(size);
            for (unsigned i=0; i<size; ++i) {
                invSnpProp[i] = 1.0/annoVec[i]->fraction;
            }
        }
        
        void compute(const VectorXf &annoTotalGenVar, const vector<AnnoInfo*> &annoInfoVec);
        void compute(const VectorXf &snpEffects, const VectorXf &annoTotalGenVar, const MatrixXf &annoMat, const MatrixXf &APA, const vector<AnnoInfo*> &annoInfoVec);
    };
    
    class AnnoPerSnpRsqEnrichment : public ParamSet {
    public:
        VectorXf invSnpProp;
        
        AnnoPerSnpRsqEnrichment(const vector<string> &header, const vector<AnnoInfo*> &annoVec, const string &lab = "Marginal_Predictability_Enrichment"):
        ParamSet(lab, header) {
            values.setOnes(size);
            invSnpProp.setZero(size);
            for (unsigned i=0; i<size; ++i) {
                invSnpProp[i] = 1.0/annoVec[i]->fraction;
            }
        }
        
        void compute(const VectorXf &snpEffectMeans, const MatrixXf &annoMat, const vector<AnnoInfo*> &annoInfoVec);
    };

    class AnnoJointPerSnpHsqEnrichment : public AnnoPerSnpHsqEnrichment, public Stat::Normal {
    public:
        
        AnnoJointPerSnpHsqEnrichment(const vector<string> &header, const vector<AnnoInfo*> &annoVec): AnnoPerSnpHsqEnrichment(header, annoVec, "Joint_Heritability_Enrichment") {}
        
        void compute(const AnnoJointProb &annoJointProb, const vector<AnnoInfo*> &annoInfoVec, const VectorXf &gamma, const float varg, const bool hsqPercModel, const float sigmaSq);
        void compute(const VectorXf &snpEffects, const MatrixXf &annoMat, const AnnoCondProb &annoCondProb, const vector<AnnoInfo*> &annoInfoVec, const VectorXf &snpAnnoCntInv);
        void compute(const VectorXf &snpEffects, const MatrixXf &annoMat, const AnnoJointProb &AnnoJointProb, const vector<AnnoInfo*> &annoInfoVec, const VectorXf &gamma, const VectorXf &snpAnnoCntInv);
//        void compute(const VectorXf &snpEffectMeans, const MatrixXf &annoMat, const AnnoEffects &alpha, const vector<AnnoInfo*> &annoInfoVec, const MatrixXf &snpPi);
    };
    
    class AnnoJointPerSnpRsqEnrichment : public AnnoPerSnpRsqEnrichment, public Stat::Normal {
    public:
        
        AnnoJointPerSnpRsqEnrichment(const vector<string> &header, const vector<AnnoInfo*> &annoVec): AnnoPerSnpRsqEnrichment(header, annoVec, "Joint_Predictability_Enrichment") {}
        
//        void compute(const VectorXf &snpEffectMeans, const MatrixXf &annoMat, const AnnoEffects &alpha, const vector<AnnoInfo*> &annoInfoVec);
        void compute(const VectorXf &snpEffectMeans, const MatrixXf &annoMat, const AnnoCondProb &annoCondProb, const vector<AnnoInfo*> &annoInfoVec, const VectorXf &snpAnnoCntInv);
        void compute(const VectorXf &snpEffectMeans, const MatrixXf &annoMat, const AnnoJointProb &AnnoJointProb, const vector<AnnoInfo*> &annoInfoVec, const VectorXf &gamma, const VectorXf &snpAnnoCntInv);
    };

    class AnnoDistribution : public vector<ParamSet*> {
    public:
        vector<string> colnames;
        unsigned numDist;
        unsigned numAnno;

        AnnoDistribution(const vector<string> &header, const unsigned numDist, const string &lab = "AnnoDistribution"):
        numDist(numDist) {
            colnames.resize(numDist);
            numAnno = header.size();
            for (unsigned i = 0; i<numDist; ++i) {
                colnames[i] = "AnnoDistribution_k" + to_string(static_cast<long long>(i + 1));
                this->push_back(new ParamSet(colnames[i], header));
            }
        }
        
        void compute(const MatrixXf &z, const MatrixXf &annoMat, const ArrayXf &numSnpMix);
    };
            
    SnpEffects snpEffects;
    AnnoEffects annoEffects;
    VarAnnoEffects sigmaSqAnno;
    AnnoCondProb annoCondProb;
    AnnoJointProb annoJointProb;
    AnnoGenVar annoGenVar;
    AnnoTotalGenVar annoTotalGenVar;
    AnnoPerSnpHsqEnrichment annoPerSnpHsqEnrich;
    AnnoPerSnpRsqEnrichment annoPerSnpRsqEnrich;
    AnnoJointPerSnpHsqEnrichment annoJointPerSnpHsqEnrich;
    AnnoJointPerSnpRsqEnrichment annoJointPerSnpRsqEnrich;
    AnnoDistribution annoDist;
        
    MatrixXf snpP;    // p = Pr(k>i | k>i-1); p2 = pi2+pi3+pi4; p3 = (pi3+pi4)/(pi2+pi3+pi4); p4 = pi4/(pi3+pi4)
    MatrixXf snpPi;   // pi1 = 1-p2; pi2 = (1-p3)*p2; pi3 = (1-p4)*p2*p3; pi4 = p2*p3*p4
    
    VectorXf snpAnnoCntInv;
    
    bool estimateRsqEnrich;
            
    ApproxBayesRC(const Data &data, const bool lowRank, const float varGenotypic, const float varResidual, const VectorXf pis, const VectorXf &piPar, const VectorXf gamma, const bool estimatePi, const bool noscale, const bool hsqPercModel, const bool robustMode, const bool estimateRsqEnrich, const string &alg, const bool message = true):
    ApproxBayesR(data, lowRank, varGenotypic, varResidual, pis, piPar, gamma, estimatePi, noscale, hsqPercModel, robustMode, alg, false),
    snpEffects(data.snpEffectNames, pis),
    annoEffects(data.annoNames, pis.size(), data.annoMat),
    sigmaSqAnno(annoEffects.colnames, annoEffects.numAnno),
    annoCondProb(data.annoNames, annoEffects.numComp),
    annoJointProb(data.annoNames, pis.size()),
    annoGenVar(data.annoNames, pis.size(), data.numKeptInds),
    annoTotalGenVar(data.annoNames),
    annoPerSnpHsqEnrich(data.annoNames, data.annoInfoVec),
    annoPerSnpRsqEnrich(data.annoNames, data.annoInfoVec),
    annoJointPerSnpHsqEnrich(data.annoNames, data.annoInfoVec),
    annoJointPerSnpRsqEnrich(data.annoNames, data.annoInfoVec),
    annoDist(data.annoNames, pis.size()),
    estimateRsqEnrich(estimateRsqEnrich)
    {
                
        initSnpPandPi(pis, data.numIncdSnps, snpP, snpPi);
        
        annoEffects.initIntercept_probit(pis);
//        if (algorithm == gibbs) annoEffects.initIntercept_probit(pis);
//        else if (algorithm == mh) annoEffects.initIntercept_logistic(pis);
//        else cout << "ERROR: unknown algorithm " << algorithm << endl;
        
        paramSetVec = {&snpEffects, &snpPip, &snpHsqPep};
        paramSetVec.insert(paramSetVec.end(), deltaPi.begin(), deltaPi.end());
        paramSetVec.insert(paramSetVec.end(), annoEffects.begin(), annoEffects.end());
        paramSetVec.insert(paramSetVec.end(), annoCondProb.begin(), annoCondProb.end());
        paramSetVec.insert(paramSetVec.end(), annoJointProb.begin(), annoJointProb.end());
        paramSetVec.insert(paramSetVec.end(), annoGenVar.begin(), annoGenVar.end());
        paramSetVec.push_back(&annoTotalGenVar);
        paramSetVec.push_back(&annoPerSnpHsqEnrich);
        paramSetVec.push_back(&annoJointPerSnpHsqEnrich);
        if (estimateRsqEnrich) {
            paramSetVec.push_back(&annoPerSnpRsqEnrich);
            paramSetVec.push_back(&annoJointPerSnpRsqEnrich);
        }

        paramVec    = {&nnzSnp, &sigmaSq, &hsq, &vare};
        paramVec.insert(paramVec.end(), numSnps.begin(), numSnps.end());
        paramVec.insert(paramVec.end(), Vgs.begin(), Vgs.end());
        if (data.numFixedEffectSnps)
            paramVec.insert(paramVec.begin() + 3, &hsqFixed);
        
        paramSetToPrint.resize(0);
        paramSetToPrint.insert(paramSetToPrint.end(), annoEffects.begin(), annoEffects.end());
        paramSetToPrint.insert(paramSetToPrint.end(), annoCondProb.begin(), annoCondProb.end());
        paramSetToPrint.insert(paramSetToPrint.end(), annoJointProb.begin(), annoJointProb.end());
        paramSetToPrint.insert(paramSetToPrint.end(), annoGenVar.begin(), annoGenVar.end());
        paramSetToPrint.push_back(&annoTotalGenVar);
        paramSetToPrint.push_back(&annoPerSnpHsqEnrich);
        paramSetToPrint.push_back(&annoJointPerSnpHsqEnrich);
        if (estimateRsqEnrich) {
            paramSetToPrint.push_back(&annoPerSnpRsqEnrich);
            paramSetToPrint.push_back(&annoJointPerSnpRsqEnrich);
        }

        paramToPrint = {&sigmaSq, &hsq, &vare};
        if (data.numFixedEffectSnps)
            paramToPrint.insert(paramToPrint.begin() + 2, &hsqFixed);
        paramToPrint.insert(paramToPrint.begin(), Vgs.begin(), Vgs.end());
        paramToPrint.insert(paramToPrint.begin(), numSnps.begin(), numSnps.end());

        if (lowRankModel) {
            paramSetVec.push_back(&vargBlk);
            paramSetVec.push_back(&vareBlk);
            paramToPrint.push_back(&nBadSnps);
        }

        if (message) {
            cout << "\nSBayesRC" << endl;
            if (lowRankModel) {
                cout << "Using the low-rank model" << endl;
            }
            cout << "scale factor: " << sigmaSq.scale << endl;
            cout << "Gamma: " << gamma.transpose() << endl;
            if (noscale) {
               cout << "Fitting model assuming unscaled genotypes " << endl;
            } else {
               cout << "Fitting model assuming scaled genotypes "  << endl;
            }
            if (!hsqPercModel) cout << "The SNP effect prior is a mixture distribution with an unknown variance variable." << endl;
            if (robustMode) cout << "Using a more robust parameterisation " << endl;
            if (algorithm == tgs_thin) cout << "Using tempered Gibbs sampling (TGS)" << endl;
            if (data.numFixedEffectSnps)
                cout << "Fitting " << data.numFixedEffectSnps << " SNP(s) as fixed effects (flat prior)." << endl;
        }
        
        getSnpAnnoCntInv(data.annoMat, snpAnnoCntInv);
    }

    void sampleUnknowns(const unsigned iter);
//    void sampleUnknownsTGS(vector<vector<int> > &selectedSnps);
    void computePfromPi(const MatrixXf &snpPi, MatrixXf &snpP);
    void computePiFromP(const MatrixXf &snpP, MatrixXf &snpPi);
    void initSnpPandPi(const VectorXf &pis, const unsigned numSnps, MatrixXf &snpP, MatrixXf &snpPi);
    void getSnpAnnoCntInv(const MatrixXf &annoMat, VectorXf &snpAnnoCntInv);
};


class BayesRC : public BayesR {
public:

    class SnpEffects : public BayesR::SnpEffects {
    public:
        unsigned ndist;
        ArrayXf numSnpMix;
        MatrixXf z;
        vector<vector<unsigned> > snpset;
        
        SnpEffects(const vector<string> &header, const VectorXf &pis): BayesR::SnpEffects(header, "Gibbs"){
            ndist = pis.size();
            numSnpMix.setZero(ndist);
            z.setZero(size, ndist-1);
        }
        
        void sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &Rsqrt, const bool weightedRes,
                          const float sigmaSq, const VectorXf &pis,  const VectorXf &gamma,
                          const float vare, VectorXf &ghat, const MatrixXf &snpPi,
                          const float varg, const bool hsqPercModel, DeltaPi &deltaPi);
    };
    
    SnpEffects snpEffects;
    ApproxBayesRC::AnnoEffects annoEffects;
    ApproxBayesRC::VarAnnoEffects sigmaSqAnno;
    ApproxBayesRC::AnnoCondProb annoCondProb;
    ApproxBayesRC::AnnoJointProb annoJointProb;
    ApproxBayesRC::AnnoGenVar annoGenVar;
    ApproxBayesRC::AnnoTotalGenVar annoTotalGenVar;
    ApproxBayesRC::AnnoPerSnpHsqEnrichment annoPerSnpHsqEnrich;
    ApproxBayesRC::AnnoPerSnpRsqEnrichment annoPerSnpRsqEnrich;
    ApproxBayesRC::AnnoJointPerSnpHsqEnrichment annoJointPerSnpHsqEnrich;
    ApproxBayesRC::AnnoJointPerSnpRsqEnrichment annoJointPerSnpRsqEnrich;

    MatrixXf snpP;    // p = Pr(k>i | k>i-1); p2 = pi2+pi3+pi4; p3 = (pi3+pi4)/(pi2+pi3+pi4); p4 = pi4/(pi3+pi4)
    MatrixXf snpPi;   // pi1 = 1-p2; pi2 = (1-p3)*p2; pi3 = (1-p4)*p2*p3; pi4 = p2*p3*p4

    enum {gibbs, mh} algorithm;
    
    bool estimateRsqEnrich;

    VectorXf snpAnnoCntInv;

    BayesRC(const Data &data, const float varGenotypic, const float varResidual, const float varRandom, const VectorXf pis, const VectorXf &piPar, const VectorXf gamma, const bool estimatePi, const bool noscale, const bool hsqPercModel, const bool estimateRsqEnrich,
            const string &alg, const bool message = true):
    BayesR(data, varGenotypic, varResidual, varRandom, pis, piPar, gamma, estimatePi, noscale, hsqPercModel, alg, false),
    snpEffects(data.snpEffectNames, pis),
    annoEffects(data.annoNames, pis.size(), data.annoMat),
    sigmaSqAnno(annoEffects.colnames, annoEffects.numAnno),
    annoCondProb(data.annoNames, annoEffects.numComp),
    annoJointProb(data.annoNames, pis.size()),
    annoGenVar(data.annoNames, pis.size(), data.numKeptInds),
    annoTotalGenVar(data.annoNames),
    annoPerSnpHsqEnrich(data.annoNames, data.annoInfoVec),
    annoPerSnpRsqEnrich(data.annoNames, data.annoInfoVec),
    annoJointPerSnpHsqEnrich(data.annoNames, data.annoInfoVec),
    annoJointPerSnpRsqEnrich(data.annoNames, data.annoInfoVec),
    estimateRsqEnrich(estimateRsqEnrich)
    {
        initSnpPandPi(pis, data.numIncdSnps, snpP, snpPi);
        
        //if (alg == "Gibbs") {
        //    algorithm = gibbs;
            annoEffects.initIntercept_probit(pis);
        //} else if (alg == "MH") {
        //    algorithm = mh;
        //    annoEffects.initIntercept_logistic(pis);
        //} else cout << "ERROR: unknown algorithm " << alg << endl;
        
        paramSetVec  = {&snpEffects, &fixedEffects, &snpPip};
        paramSetVec.insert(paramSetVec.end(), deltaPi.begin(), deltaPi.end());
        paramSetVec.insert(paramSetVec.end(), annoEffects.begin(), annoEffects.end());
        paramSetVec.insert(paramSetVec.end(), annoCondProb.begin(), annoCondProb.end());
        paramSetVec.insert(paramSetVec.end(), annoJointProb.begin(), annoJointProb.end());
        paramSetVec.insert(paramSetVec.end(), annoGenVar.begin(), annoGenVar.end());
        paramSetVec.push_back(&annoTotalGenVar);
        paramSetVec.push_back(&annoPerSnpHsqEnrich);
        paramSetVec.push_back(&annoJointPerSnpHsqEnrich);
        if (estimateRsqEnrich) {
            paramSetVec.push_back(&annoPerSnpRsqEnrich);
            paramSetVec.push_back(&annoJointPerSnpRsqEnrich);
        }

        paramVec     = {&nnzSnp, &sigmaSq, &varg, &vare, &hsq};
        paramVec.insert(paramVec.end(), numSnps.begin(), numSnps.end());
        paramVec.insert(paramVec.end(), Vgs.begin(), Vgs.end());
        
        paramSetToPrint.resize(0);
        paramSetToPrint.insert(paramSetToPrint.end(), annoEffects.begin(), annoEffects.end());
        paramSetToPrint.insert(paramSetToPrint.end(), annoCondProb.begin(), annoCondProb.end());
        paramSetToPrint.insert(paramSetToPrint.end(), annoJointProb.begin(), annoJointProb.end());
        paramSetToPrint.insert(paramSetToPrint.end(), annoGenVar.begin(), annoGenVar.end());
        paramSetToPrint.push_back(&annoTotalGenVar);
        paramSetToPrint.push_back(&annoPerSnpHsqEnrich);
        paramSetToPrint.push_back(&annoJointPerSnpHsqEnrich);
        if (estimateRsqEnrich) {
            paramSetToPrint.push_back(&annoPerSnpRsqEnrich);
            paramSetToPrint.push_back(&annoJointPerSnpRsqEnrich);
        }

        paramToPrint = {&sigmaSq, &varg, &vare, &hsq};
        paramToPrint.insert(paramToPrint.begin(), Vgs.begin(), Vgs.end());
        paramToPrint.insert(paramToPrint.begin(), numSnps.begin(), numSnps.end());
        
        if (data.numRandomEffects) {
            paramSetVec.push_back(&randomEffects);
            paramVec.push_back(&sigmaSqRand);
            paramVec.push_back(&varRand);
            paramToPrint.push_back(&varRand);
        }
        //paramToPrint.push_back(&rounding);
        
        if (message) {
            cout << "\nBayesRC model fitted. Algorithm: " << alg << "." << endl;
            cout << "scale factor: " << sigmaSq.scale << endl;
            cout << "Gamma: " << gamma.transpose() << endl;
        }
        
        getSnpAnnoCntInv(data.annoMat, snpAnnoCntInv);
    }
    
    void sampleUnknowns(const unsigned iter);
    void computePiFromP(const MatrixXf &snpP, MatrixXf &snpPi);
    void initSnpPandPi(const VectorXf &pis, const unsigned numSnps, MatrixXf &snpP, MatrixXf &snpPi);
    void getSnpAnnoCntInv(const MatrixXf &annoMat, VectorXf &snpAnnoCntInv);
};



class ApproxBayesRD : public ApproxBayesRC {
public:
    
    class AnnoEffects : public vector<BayesC::SnpEffects*>, public Stat::TruncatedNormal  {
    public:
        unsigned numComp;  // number of components = number of mixture components - 1
        unsigned numAnno;  // number of annotations
        unsigned numAnnoTotal;
        
        MatrixXf values;
        MatrixXf wcorr;
        VectorXf annoDiag;
        VectorXf ssq;
        VectorXf numNonZeros;
        vector<string> colnames;
        
        VectorXf pip;
        float nnz;

        AnnoEffects(const vector<string> &header, const unsigned ndist, const MatrixXf &annoMat) {
            numComp = ndist - 1;
            colnames.resize(numComp);
            numAnno = header.size();
            unsigned numSnps = annoMat.rows();
            for (unsigned i = 0; i<numComp; ++i) {
                colnames[i] = "AnnoEffects_p" + to_string(static_cast<long long>(i + 2));
                this->push_back(new BayesC::SnpEffects(header, "Gibbs", colnames[i]));
            }
            wcorr.setZero(numSnps, numComp);
            annoDiag.setZero(numAnno);
            annoDiag[0] = numSnps;  // first annotation is intercept
            for (unsigned j=1; j<numAnno; ++j) {
                annoDiag[j] = annoMat.col(j).squaredNorm();
            }
            nnz = 0;
            ssq.setZero(numComp);
            numNonZeros.setZero(numAnno);
            numAnnoTotal = (numAnno-1)*numComp;  // leave out the intercept
            values.setZero(numAnno,numComp);
            pip.setZero(numAnno);
        }
        
        void sampleFromFC_indep(MatrixXf &z, const MatrixXf &annoMat, const VectorXf &sigmaSq, const float pi, MatrixXf &snpP);
        void sampleFromFC_joint(MatrixXf &z, const MatrixXf &annoMat, const VectorXf &sigmaSq, const float pi, MatrixXf &snpP);
        void initIntercept(const VectorXf &pis);
    };
    
    class AnnoPi : public BayesC::Pi {
    public:
        AnnoPi(): BayesC::Pi(0.1, 1, 1, "AnnoPi"){}
    };
    
    class AnnoCondProb : public ApproxBayesRC::AnnoCondProb {
    public:
        AnnoCondProb(const vector<string> &header, const unsigned numComp):
        ApproxBayesRC::AnnoCondProb(header, numComp){}
        
        void compute(const AnnoEffects &annoEffects, const vector<AnnoInfo*> &annoInfoVec);
    };
    
    class AnnoPIP : public BayesC::SnpPIP {
    public:
        AnnoPIP(const vector<string> &header, const string &lab = "AnnoPIP"): BayesC::SnpPIP(header, lab){}
    };

    AnnoEffects annoEffects;
    AnnoPi piAnno;
    AnnoCondProb annoCondProb;
    AnnoPIP annoPip;

    ApproxBayesRD(const Data &data, const bool lowRank, const float varGenotypic, const float varResidual, const VectorXf pis, const VectorXf &piPar, const VectorXf gamma, const bool estimatePi, const bool noscale, const bool hsqPercModel, const bool robustMode, const bool estimateRsqEnrich, const string &alg, const bool message = true):
    ApproxBayesRC(data, lowRank, varGenotypic, varResidual, pis, piPar, gamma, estimatePi, noscale, hsqPercModel, robustMode, estimateRsqEnrich, alg, false),
    annoEffects(data.annoNames, pis.size(), data.annoMat),
    annoCondProb(data.annoNames, annoEffects.numComp),
    annoPip(data.annoNames)
    {
        annoEffects.initIntercept(pis);
        
        paramSetVec = {&snpEffects, &snpPip, &snpHsqPep, &annoPip};
        paramSetVec.insert(paramSetVec.end(), deltaPi.begin(), deltaPi.end());
        paramSetVec.insert(paramSetVec.end(), annoEffects.begin(), annoEffects.end());
        paramSetVec.insert(paramSetVec.end(), annoCondProb.begin(), annoCondProb.end());
        paramSetVec.insert(paramSetVec.end(), annoJointProb.begin(), annoJointProb.end());
        paramSetVec.insert(paramSetVec.end(), annoGenVar.begin(), annoGenVar.end());
        paramSetVec.push_back(&annoTotalGenVar);
        paramSetVec.push_back(&annoPerSnpHsqEnrich);
        paramSetVec.push_back(&annoJointPerSnpHsqEnrich);

        paramVec = {&nnzSnp, &sigmaSq, &hsq, &vare, &piAnno};
        if (data.numFixedEffectSnps)
            paramVec.insert(paramVec.begin() + 3, &hsqFixed);
        paramVec.insert(paramVec.end(), numSnps.begin(), numSnps.end());
        paramVec.insert(paramVec.end(), Vgs.begin(), Vgs.end());
        
        paramSetToPrint.resize(0);
        paramSetToPrint.insert(paramSetToPrint.end(), annoEffects.begin(), annoEffects.end());
        paramSetToPrint.insert(paramSetToPrint.end(), annoCondProb.begin(), annoCondProb.end());
        paramSetToPrint.insert(paramSetToPrint.end(), annoJointProb.begin(), annoJointProb.end());
        paramSetToPrint.insert(paramSetToPrint.end(), annoGenVar.begin(), annoGenVar.end());
        paramSetToPrint.push_back(&annoTotalGenVar);
        paramSetToPrint.push_back(&annoPerSnpHsqEnrich);
        paramSetToPrint.push_back(&annoJointPerSnpHsqEnrich);
        paramSetToPrint.push_back(&annoPip);

        paramToPrint = {&sigmaSq, &hsq, &vare, &piAnno};
        if (data.numFixedEffectSnps)
            paramToPrint.insert(paramToPrint.begin() + 2, &hsqFixed);
        paramToPrint.insert(paramToPrint.begin(), Vgs.begin(), Vgs.end());
        paramToPrint.insert(paramToPrint.begin(), numSnps.begin(), numSnps.end());

        if (lowRankModel) {
            paramSetVec.push_back(&vargBlk);
            paramSetVec.push_back(&vareBlk);
            paramToPrint.push_back(&nBadSnps);
        }

        if (message) {
            cout << "\nSBayesRD" << endl;
            if (lowRankModel) {
                cout << "Using the low-rank model" << endl;
            }
            cout << "scale factor: " << sigmaSq.scale << endl;
            cout << "Gamma: " << gamma.transpose() << endl;
            if (noscale) {
               cout << "Fitting model assuming unscaled genotypes " << endl;
            } else {
               cout << "Fitting model assuming scaled genotypes "  << endl;
            }
            if (!hsqPercModel) cout << "The SNP effect prior is a mixture distribution with an unknown variance variable." << endl;
            if (robustMode) cout << "Using a more robust parameterisation " << endl;
        }
    }

    void sampleUnknowns(const unsigned iter);

};


// -----------------------------------------------------------------------------------------------
// ApproxBayesAPP: Bivariate Bayesian Analysis with Annotation and Pleiotropy
// -----------------------------------------------------------------------------------------------

class ApproxBayesAPP : public ApproxBayesC {
    // Bivariate summary statistics Bayesian model with:
    // - 4-state inclusion model: [1,1], [1,0], [0,1], [0,0]
    // - Annotation-based variance components
    // - LD block structure with eigen decomposition
public:
    
    // Forward declaration
    class AnnoPiBivariate;
    
    // 4-state inclusion indicator: [δ₁, δ₂] for each SNP per category
    // Structure: [d1_c1, d2_c1, ..., d_nsnp_c1, d1_c2, d2_c2, ..., d_nsnp_c2, ...] × 2 traits
    class DeltaBivariate {
    public:
        MatrixXf deltaMatrix;   // (nSNPs × nCategory) × 2: [delta1, delta2] for each SNP per category
        unsigned size;          // Number of SNPs
        vector<VectorXf> nLociCounts;  // per-category counts of SNPs in each of the 4 states

        DeltaBivariate(const vector<string> &header, const unsigned nCat):
        size(header.size()), nLociCounts(nCat, VectorXf::Zero(4)) {
            deltaMatrix.setZero(size * nCat, 2);
        }

        // Recount SNPs in each of the 4 inclusion states per annotation category.
        // State index (matching Julia): 0=[1,1], 1=[1,0], 2=[0,1], 3=[0,0]
        // Updates nLociCounts in place.
        void countLociByState(const unsigned numCategories);

        // Count SNPs with at least one non-zero delta across all categories
        unsigned countNonZero(const unsigned numCategories) const;
    };
    
    // Bivariate SNP effects: [β₁, β₂] for each SNP per category
    // Structure: [b1_c1, b2_c1, ..., b_nsnp_c1, b1_c2, b2_c2, ..., b_nsnp_c2, ...] × 2 traits
    // We use separate SnpEffects objects for each trait so MCMC can collect samples separately
    class SnpEffectsBivariate {
    public:
        vector<ApproxBayesC::SnpEffects> betaTotal;  // One per trait (size 2) - stores total SNP effect (sum over categories) for MCMC collection
        MatrixXf betaMatrix;  // (nSNPs × nCategory) × 2: [beta1, beta2] for each SNP per category - for internal computation during sampling
        MatrixXf alphaMatrix; // (nSNPs × nCategory) × 2: [alpha1, alpha2] = delta * beta
        unsigned numCategories;
        
        SnpEffectsBivariate(const vector<string> &header, const unsigned nCat): 
        numCategories(nCat),
        betaTotal(2, ApproxBayesC::SnpEffects(header)) {
            // Size is nSNPs * nCategory (flattened structure like Julia)
            unsigned totalSize = betaTotal[0].size * nCat;
            betaMatrix.setZero(totalSize, 2);
            alphaMatrix.setZero(totalSize, 2);
        }
        
        // Sample bivariate SNP effects with 4-state inclusion model
        void sampleFromFC(MatrixXf &deltaMatrix, vector<VectorXf> &wcorrBlocks, const vector<MatrixXf> &Qblocks,
                          vector<VectorXf> &whatBlocks, VectorXf &vareBlocks,
                          const vector<Matrix2f> &A_vec, const AnnoPiBivariate &Pi,
                          const VectorXf &snp2pq, const vector<Matrix2f> &Rinv_blk,
                          const vector<LDBlockInfo*> &keptLdBlockInfoVec,
                          const vector<Vector2f> &nGWASblocks, const MatrixXf &annoMat);
        
        // Update values vectors from betaMatrix/alphaMatrix after sampling
        void updateValues(void);

    };
    
    // 4-state inclusion probabilities per annotation category (matching Julia order)
    // States: 0=[1,1], 1=[1,0], 2=[0,1], 3=[0,0]
    class AnnoPiBivariate : public vector<ParamSet*> {
    public:
        vector<VectorXf> piVec;  // per category: VectorXf(4) with probabilities for [0,0], [1,0], [0,1], [1,1]
        VectorXf alphaVec;  // Dirichlet hyperparameters (size 4)
        unsigned numCategories;
        
        AnnoPiBivariate(const vector<string> &header, const unsigned nCat):
        numCategories(nCat) {
            alphaVec = VectorXf::Ones(4);
            
            for (unsigned c = 0; c < numCategories; ++c) {
                VectorXf pi_c = VectorXf::Constant(4, 0.25);  // uniform prior
                piVec.push_back(pi_c);
            }
        }
        
        // Get probability for state [δ₁, δ₂] (matching Julia order)
        // States: 0=[1,1], 1=[1,0], 2=[0,1], 3=[0,0]
        float getPi(unsigned cat, float delta1, float delta2) const {
            unsigned idx = (delta1 < 0.5f ? 0 : 1) + (delta2 < 0.5f ? 0 : 2);
            return piVec[cat][idx];
        }
        
        void sampleFromFC(const vector<VectorXf> &nLociCounts);  // counts per category (each VectorXf size 4)
    };
    
    // Marginal 4-state pi across all annotation categories, stored as vector<Parameter*>
    // so elements can be inserted into paramToPrint / paramVec directly.
    // States (matching Julia): 0=[1,1], 1=[1,0], 2=[0,1], 3=[0,0]
    class PiBivariate : public vector<Parameter*> {
    public:
        PiBivariate() {
            for (const char *name : {"Pi_00", "Pi_10", "Pi_01", "Pi_11"})
                this->push_back(new Parameter(name));
            for (auto *p : *this) p->value = 0.25f;
        }
        
        void compute(const AnnoPiBivariate &piAnno, const VectorXf &nLociAnno);
    };

    // Per-trait heritability and coheritability between traits
    class HsqBivariate : public vector<Parameter*> {
    public:
        HsqBivariate() {
            for (const char *name : {"hsq1", "hsq2", "cohsq"})
                this->push_back(new Parameter(name));
        }

        void compute(const Matrix2f &vargTotal, const vector<Matrix2f> &R_blk);
    };

    // Per-trait SNP PIPs.
    class SnpPIPBivariate : public vector<ParamSet*> {
    public:
        SnpPIPBivariate(const vector<string> &snpNames) {
            this->push_back(new BayesC::SnpPIP(snpNames, "PIP1")); // trait 1
            this->push_back(new BayesC::SnpPIP(snpNames, "PIP2")); // trait 2
        }

        BayesC::SnpPIP &pip1() { return *static_cast<BayesC::SnpPIP*>((*this)[0]); }
        BayesC::SnpPIP &pip2() { return *static_cast<BayesC::SnpPIP*>((*this)[1]); }

        void compute(const MatrixXf &deltaMatrix, const unsigned numCategories);
    };

    // Bivariate marker effect variance matrix (2×2) per annotation category
    class AnnoVarEffectsBivariate : public vector<Parameter*> {
    public:
        unsigned numCategories; // number of annotation categories

        vector<Matrix2f> A_vec;  // per category: 2×2 covariance matrix
        vector<Matrix2f> Ainv_vec;  // inverse matrices
        
        AnnoVarEffectsBivariate(const unsigned numCategories, const float varGenotypic1, 
                           const float varGenotypic2, const VectorXf &nLociAnno):
        numCategories(numCategories) {
            A_vec.resize(numCategories);
            Ainv_vec.resize(numCategories);
            
            for (unsigned c = 0; c < numCategories; ++c) {
                // Initialize with scaled heritability estimates
                float scale = nLociAnno[c] / nLociAnno.sum();
                A_vec[c] << varGenotypic1 * scale, 0.0,
                           0.0, varGenotypic2 * scale;
                Ainv_vec[c] = A_vec[c].inverse();
            }
        }
        
        void sampleFromFC(const MatrixXf &betaMatrix, const VectorXf &nLociAnno);        
    };
    
    // Per-block accumulator for genetic variance components (zeroed each iteration)
    class VarBlocks {
    public:
        vector<MatrixXf> vargBlkCat;    // per block: 2 × nCategory genetic variance per trait
        vector<VectorXf> ssqBlkCat;     // per block: (2*nCategory + nCategory) sum of alpha^2 + covariance
        vector<VectorXf> vargCovBlkCat; // per block: nCategory genetic covariance
        vector<Matrix2f> totalVargBlk;  // per block: total 2×2 genetic variance

        void resize(const unsigned nBlocks, const unsigned nCategories);

        void setZero(const unsigned nBlocks);

        void compute(const MatrixXf &alphaMatrix, const unsigned nSNPs,
                     const vector<LDBlockInfo*> &keptLdBlockInfoVec,
                     const vector<MatrixXf> &Qblocks, const unsigned numCategories);

        // Sum totalVargBlk across blocks
        Matrix2f computeTotal() const;
    };

    // Bivariate residual variance matrix (2×2) per LD block
    class ResidualVarBivariate : public ParamSet {
    public:
        vector<Matrix2f> R_blk;      // per block: 2×2 residual covariance (raw)
        vector<Matrix2f> Rinv_blk;   // per block: inverse of sample-size-normalised R (for Gibbs sampler)
        Matrix2f Rprior;
        const float df_R;
        Matrix2f scale_R;
        
        ResidualVarBivariate(const vector<string> &header, const float varPhenotypic1,
                            const float varPhenotypic2, const string &lab = "ResidualVarBivariate"):
        ParamSet(lab, header), df_R(6) {  // df = 4 + nTraits
            Rprior << varPhenotypic1, 0.0,
                     0.0, varPhenotypic2;
            scale_R = Rprior * (df_R - 3);  // df_R - nTraits - 1
            R_blk.resize(size);
            Rinv_blk.resize(size);
            for (unsigned b = 0; b < size; ++b) {
                R_blk[b] = Rprior;
                Rinv_blk[b] = Rprior.inverse();  // will be properly normalised on first sampleFromFC call
            }
        }
        
        // Normalise R_blk by sample sizes and update Rinv_blk
        void normalise(const vector<Vector2f> &nGWASblocks) {
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

        void sampleFromFC(vector<VectorXf> &wcorrBlocks, const vector<Vector2f> &nGWASblocks,
                          const VectorXf &numEigenvalBlock, const VarBlocks &varBlocks);
    };
    
    // Bivariate genetic variance matrix (2×2) per annotation category
    class AnnoGenotypicVarBivariate : public vector<ParamSet*> {
    public:
        vector<Matrix2f> G_vec;  // per category: 2×2 genetic covariance
        
        AnnoGenotypicVarBivariate(const vector<string> &header, const unsigned numCategories):
        numCategories(numCategories) {
            G_vec.resize(numCategories);
            for (unsigned c = 0; c < numCategories; ++c) {
                G_vec[c].setZero();
            }
        }
        
        void compute(const vector<VectorXf> &whatBlocks, const vector<MatrixXf> &Qblocks,
                    const vector<LDBlockInfo*> &keptLdBlockInfoVec,
                    const MatrixXf &alphaMatrix, const MatrixXf &annoMatrix);

        // Accumulate G_vec from per-block variance accumulators
        void accumulate(const VarBlocks &varBlocks);
        
        unsigned numCategories;
    };
    
    SnpEffectsBivariate snpEffects;
    DeltaBivariate delta;
    AnnoPiBivariate piAnno;  // annotation-stratified 4-state inclusion probabilities
    PiBivariate pi;          // marginal 4-state pi across annotations (for output/tracking)
    AnnoVarEffectsBivariate sigmaSq;
    ResidualVarBivariate vare;
    AnnoGenotypicVarBivariate varg;
    SnpPIPBivariate snpPip;  // per-trait SNP PIPs [PIP, PIP2]
    HsqBivariate hsq;  // per-trait heritability [hsq1, hsq2]
    
    // Precomputed constants (do not change across iterations)
    vector<Vector2f> nGWASblocks;  // per-block bivariate sample sizes
    MatrixXf annoMatIncd;          // annotation matrix for included SNPs only
    VectorXf nLociPerAnno;         // number of SNPs per annotation category (precomputed)

    // Total genetic variance (sum across categories)
    Matrix2f vargTotal;
    
    VarBlocks varBlocks;  // per-block genetic variance accumulators
    
    ApproxBayesAPP(const Data &data, const bool lowrank, const float varGenotypic1, 
              const float varGenotypic2, const float varResidual1, const float varResidual2,
              const VectorXf &nLociAnno, const bool estPi, const bool estVare,
              const bool estVara, const bool message = true)
    : ApproxBayesC(data, lowrank, varGenotypic1, varResidual1, 0.0, 0.01, 1.0, 1.0, estPi, false, false, false)
    , snpEffects(data.snpEffectNames, data.numAnnos > 0 ? data.numAnnos : 1)
    , delta(data.snpEffectNames, data.numAnnos > 0 ? data.numAnnos : 1)
    , piAnno(data.annoNames.size() > 0 ? data.annoNames : vector<string>(1, "Category1"), 
             data.numAnnos > 0 ? data.numAnnos : 1)
    , pi()
    , sigmaSq(data.numAnnos > 0 ? data.numAnnos : 1, varGenotypic1, varGenotypic2, 
              nLociAnno.size() > 0 ? nLociAnno : VectorXf::Constant(1, data.numIncdSnps))
    , vare(data.keptLdBlockInfoVec.size() > 0 ? data.ldblockNames : vector<string>(1, "ResidualVar"), 
           varResidual1, varResidual2)
    , varg(data.annoNames.size() > 0 ? data.annoNames : vector<string>(1, "Category1"), 
           data.numAnnos > 0 ? data.numAnnos : 1)
    , snpPip(data.snpEffectNames)
    , hsq()
    , estimatePi(estPi)
    , estimateVare(estVare)
    , estimateVara(estVara)
    , numCategories(data.numAnnos > 0 ? data.numAnnos : 1)
    {
        vargTotal.setZero();
        // Rename trait 2's SNP effects label so outputResults handles it separately
        snpEffects.betaTotal[1].label = "SnpEffects2";

        // Precompute per-block bivariate sample sizes (constant across iterations)
        unsigned nBlocks = data.keptLdBlockInfoVec.size();
        nGWASblocks.resize(nBlocks);
        if (data.nGWASblockBivariate.size() > 0) {
            nGWASblocks = data.nGWASblockBivariate;
        } else {
            for (unsigned b = 0; b < nBlocks; ++b) {
                float n = data.nGWASblock[b];
                nGWASblocks[b] << n, n;
            }
        }

        // Precompute annotation matrix for included SNPs (constant across iterations)
        unsigned nCat = data.numAnnos > 0 ? data.numAnnos : 1;
        annoMatIncd.resize(data.numIncdSnps, nCat);
        if (data.annoMat.cols() == 0) {
            annoMatIncd.setOnes();
        } else {
            for (unsigned i = 0; i < data.numIncdSnps; ++i) {
                unsigned originalIdx = data.incdSnpInfoVec[i]->index;
                if (originalIdx < (unsigned)data.annoMat.rows())
                    annoMatIncd.row(i) = data.annoMat.row(originalIdx);
                else
                    annoMatIncd.row(i).setZero();
            }
        }

        // Precompute normalised Rinv for first iteration
        vare.normalise(nGWASblocks);

        // Precompute number of SNPs per annotation category
        nLociPerAnno.resize(numCategories);
        for (unsigned c = 0; c < numCategories; ++c)
            nLociPerAnno[c] = annoMatIncd.col(c).sum();

        // Initialize per-block arrays
        varBlocks.resize(nBlocks, numCategories);
        
        paramSetVec = {&snpEffects.betaTotal[0], snpPip[0], &snpEffects.betaTotal[1], snpPip[1], &vare};
        paramSetVec.insert(paramSetVec.end(), piAnno.begin(), piAnno.end());
        paramSetVec.insert(paramSetVec.end(), varg.begin(), varg.end());
        paramVec = {&nnzSnp};
        paramVec.insert(paramVec.end(), hsq.begin(), hsq.end());
        paramVec.insert(paramVec.end(), pi.begin(), pi.end());    // marginal Pi_00, Pi_01, Pi_10, Pi_11
        paramVec.insert(paramVec.end(), sigmaSq.begin(), sigmaSq.end());
        paramToPrint = {&nnzSnp};
        paramToPrint.insert(paramToPrint.end(), pi.begin(), pi.end());
        paramToPrint.insert(paramToPrint.end(), hsq.begin(), hsq.end());
        
        if (lowRankModel) {
            paramSetVec.push_back(&vargBlk);
            paramSetVec.push_back(&vareBlk);
        }
        
        if (message) {
            cout << "\nSBayesAPP: Bivariate Bayesian Analysis with Annotation-stratified co-Polygenicity and Pleiotropy" << endl;
            cout << "Number of annotation categories: " << numCategories << endl;
            cout << "Estimate Pi: " << (estimatePi ? "Yes" : "No") << endl;
            cout << "Estimate residual variance: " << (estimateVare ? "Yes" : "No") << endl;
            cout << "Estimate marker variance: " << (estimateVara ? "Yes" : "No") << endl;
        }
    }
    
    void sampleUnknowns(const unsigned iter);
    
private:
    bool estimatePi;
    bool estimateVare;
    bool estimateVara;
    unsigned numCategories;
};

#endif /* model_hpp */




