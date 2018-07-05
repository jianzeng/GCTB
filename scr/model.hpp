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
    const string label;
    float value;   // sampled value
    
    Parameter(const string &label): label(label){
        value = 0.0;
    }
};

class ParamSet {
    // base class for a set of parameters of same kind, e.g. fixed effects, snp effects ...
public:
    const string label;
    const vector<string> &header;
    const unsigned size;
    VectorXf values;
        
    ParamSet(const string &label, const vector<string> &header)
    : label(label), header(header), size(int(header.size())){
        values.setZero(size);
    }
};

class Model {
public:
    unsigned numSnps;
    
    vector<ParamSet*> paramSetVec;
    vector<Parameter*> paramVec;
    vector<Parameter*> paramToPrint;
    vector<ParamSet*> paramSetToPrint;
    
    virtual void sampleUnknowns(void) = 0;
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
    
    class SnpEffects : public ParamSet, public Stat::NormalZeroMixture {
        // all snp effects has a mixture prior of a nomral distribution and a point mass at zero
    public:
        float sumSq;
        unsigned numNonZeros;
        
        enum {gibbs, hmc} algorithm;
        
        unsigned cnt;
        float mhr;

        
        SnpEffects(const vector<string> &header, const string &alg, const string &lab = "SnpEffects")
        : ParamSet(lab, header){
            sumSq = 0.0;
            numNonZeros = 0;
            if (alg=="HMC") algorithm = hmc;
            else algorithm = gibbs;
            cnt = 0;
            mhr = 0.0;
        }
        
        void sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag,
                          const float sigmaSq, const float pi, const float vare, VectorXf &ghat);
        void gibbsSampler(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag,
                          const float sigmaSq, const float pi, const float vare, VectorXf &ghat);
        void hmcSampler(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag,
                        const float sigmaSq, const float pi, const float vare, VectorXf &ghat);
        ArrayXf gradientU(const VectorXf &alpha, const MatrixXf &ZPZ, const VectorXf &ypZ,
                        const float sigmaSq, const float vare);
        float computeU(const VectorXf &alpha, const MatrixXf &ZPZ, const VectorXf &ypZ,
                       const float sigmaSq, const float vare);
        
        void sampleFromFC_omp(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag,
                              const float sigmaSq, const float pi, const float vare, VectorXf &ghat);

    };
    
    class VarEffects : public Parameter, public Stat::InvChiSq {
        // variance of snp effects has a scaled-inverse chi-square prior
    public:
        const float df;  // hyperparameter
        float scale;        // hyperparameter
        
        VarEffects(const float vg, const VectorXf &snp2pq, const float pi, const string &lab = "SigmaSq")
        : Parameter(lab), df(4)
        {
            if (myMPI::partition == "bycol") {
                int sizeFull;
                MPI_Allreduce(&myMPI::iSize, &sizeFull, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                VectorXf snp2pqFull(sizeFull);
                MPI_Allgatherv((float*)&snp2pq[0], myMPI::iSize, MPI_FLOAT, &snp2pqFull[0], &myMPI::srcounts[0], &myMPI::displs[0], MPI_FLOAT, MPI_COMM_WORLD);
                value = vg/(snp2pqFull.sum()*pi);  // derived from prior knowledge on Vg and pi
            }
            else {
                value = vg/(snp2pq.sum()*pi);  // derived from prior knowledge on Vg and pi
            }
            scale = 0.5f*value;  // due to df = 4
            
            //cout << vg << " " << snp2pq.sum() << " " << pi << endl;
        }
        
        void sampleFromFC(const float snpEffSumSq, const unsigned numSnpEff);
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
        
        Pi(const float pi, const string &lab = "Pi"): Parameter(lab), alpha(1), beta(19){  // informative prior
            value = pi;
        }
        
        void sampleFromFC(const unsigned numSnps, const unsigned numSnpEff);
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
            if (myMPI::partition == "byrow") {
                MPI_Allreduce(&n, &nobs, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
            } else {
                nobs = n;
            }
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
    
    class Heritability : public Parameter {
        // compute heritability based on sampled values of genotypic and residual variances
        // strictly speaking, this is not a model parameter
    public:
        Heritability(const string &lab = "hsq"): Parameter(lab){};
        void compute(const float genVar, const float resVar){
            value = genVar/(genVar+resVar);
        }
    };
    
    class Rounding : public Parameter {
        // re-compute ycorr to eliminate rounding errors
    public:
        unsigned count;
        
        Rounding(const string &lab = "Rounding"): Parameter(lab){
            count = 0;
        }
        void computeYcorr(const VectorXf &y, const MatrixXf &X, const MatrixXf &Z,
                          const VectorXf &fixedEffects, const VectorXf &snpEffects,
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

    
public:
    const Data &data;
    
    VectorXf ycorr;   // corrected y for mcmc sampling
    VectorXf ghat;    // predicted total genotypic values
    
    bool estimatePi;
    
    FixedEffects fixedEffects;
    SnpEffects snpEffects;
    VarEffects sigmaSq;
    ScaleVar scale;
    Pi pi;
    ResidualVar vare;
    
    GenotypicVar varg;
    Heritability hsq;
    Rounding rounding;
    NumNonZeroSnp nnzSnp;
    
    BayesC(const Data &data, const float varGenotypic, const float varResidual, const float pival, const bool estimatePi,
           const string &algorithm = "Gibbs", const bool message = true):
    data(data),
    ycorr(data.y),
    fixedEffects(data.fixedEffectNames),
    snpEffects(data.snpEffectNames, algorithm),
    sigmaSq(varGenotypic, data.snp2pq, pival),
    scale(sigmaSq.scale),
    pi(pival),
    vare(varResidual, data.numKeptInds),
    varg(varGenotypic),
    estimatePi(estimatePi)
    {
        numSnps = data.numIncdSnps;
        paramSetVec = {&snpEffects, &fixedEffects};           // for which collect mcmc samples
        paramVec = {&pi, &nnzSnp, &sigmaSq, &vare, &varg, &hsq};       // for which collect mcmc samples
        paramToPrint = {&pi, &nnzSnp, &sigmaSq, &vare, &varg, &hsq, &rounding};   // print in order
        if (message && myMPI::rank==0) {
            string alg = algorithm;
            if (alg!="HMC") alg = "Gibbs (default)";
            cout << "\nBayesC model fitted. Algorithm: " << alg << "." << endl;
            cout << "scale factor: " << sigmaSq.scale << endl;
        }
    }
    
    void sampleUnknowns(void);
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
        
        void sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag,
                          const VectorXf &sigmaSq, const float pi, const float vare, VectorXf &ghat);
    };

    class VarEffects : public ParamSet, public BayesC::VarEffects {
    public:
        VarEffects(const float vg, const VectorXf &snp2pq, const float pi):
        ParamSet("SigmaSqs", vector<string>(snp2pq.size())),
        BayesC::VarEffects(vg, snp2pq, pi){
            values.setConstant(size, value);
        }
        
        void sampleFromFC(const VectorXf &betaSq);
    };
    
    SnpEffects snpEffects;
    VarEffects sigmaSq;

    BayesB(const Data &data, const float varGenotypic, const float varResidual, const float pival,
           const bool estimatePi, const bool message = true):
    BayesC(data, varGenotypic, varResidual, pival, estimatePi, "Gibbs", false),
    snpEffects(data.snpEffectNames),
    sigmaSq(varGenotypic, data.snp2pq, pival)
    {
        paramSetVec = {&snpEffects, &fixedEffects};           // for which collect mcmc samples
        paramVec = {&pi, &nnzSnp, &vare, &varg, &hsq};       // for which collect mcmc samples
        paramToPrint = {&pi, &nnzSnp, &vare, &varg, &hsq, &rounding};   // print in order
        if (message && myMPI::rank==0) {
            cout << "\nBayesB model fitted." << endl;
            cout << "scale factor: " << sigmaSq.scale << endl;
        }
    }
    
    void sampleUnknowns(void);

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
        
        SnpEffects(const vector<string> &header, const VectorXi &windStart, const VectorXi &windSize, const unsigned snpFittedPerWindow):
        BayesC::SnpEffects(header, "Gibbs"), windStart(windStart), windSize(windSize){
            numWindows = (unsigned) windStart.size();
            windDelta.setZero(numWindows);
            localPi.setOnes(numWindows);
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

        void sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag,
                          const float sigmaSq, const float pi, const float vare, VectorXf &ghat);
    };
    
    class VarEffects : public BayesC::VarEffects {
    public:
        VarEffects(const float vg, const VectorXf &snp2pq, const float pi,
                   const VectorXf &localPi, const unsigned snpFittedPerWindow):
        BayesC::VarEffects(vg, snp2pq, pi){
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
    
    BayesN(const Data &data, const float varGenotypic, const float varResidual, const float pival,
           const bool estimatePi, const unsigned snpFittedPerWindow, const bool message = true):
    BayesC(data, varGenotypic, varResidual, pival, estimatePi, "Gibbs", false),
    snpEffects(data.snpEffectNames, data.windStart, data.windSize, snpFittedPerWindow),
    sigmaSq(varGenotypic, data.snp2pq, pival, snpEffects.localPi, snpFittedPerWindow),
    windDelta(vector<string>(snpEffects.numWindows))
    {
        paramSetVec = {&snpEffects, &fixedEffects, &windDelta};           // for which collect mcmc samples
        paramVec = {&pi, &nnzWind, &nnzSnp, &sigmaSq, &vare, &varg, &hsq};       // for which collect mcmc samples
        paramToPrint = {&pi, &nnzWind, &nnzSnp, &sigmaSq, &vare, &varg, &hsq, &rounding};   // print in order
        if (message && myMPI::rank==0) {
            cout << "\nBayesN model fitted." << endl;
            cout << "scale factor: " << sigmaSq.scale << endl;
        }
    }

    void sampleUnknowns(void);
};

// -----------------------------------------------------------------------------------------------
// Bayes R
// -----------------------------------------------------------------------------------------------

class BayesR : public BayesC {
    // Prior for snp efect pi_1 * N(0, 0) + pi_2 * N(0, sig^2_beta * gamma_2) + pi_3 * N(0, sig^2_beta * gamma_3) + pi_3 * N(0, sig^2_beta * gamma_4)
    // consider S as unknown to make inference on the relationship between MAF and effect size
public:
    
    class SnpEffects : public BayesC::SnpEffects {
    public:
      float sum2pq;
        SnpEffects(const vector<string> &header, const string &alg): BayesC::SnpEffects(header, "Gibbs"){
            sum2pq = 0.0;
        }
        
        void sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag,
                          const float sigmaSq, const VectorXf &pis,  const VectorXf &gamma,
                          const float vare, VectorXf &ghat, VectorXf &snpStore);
    };

    class ProbMixComps : public vector<Parameter*>, public Stat::Dirichlet {

        // prior probability of a snp being in any of the distributions effect has a dirichlet prior
    public:
        VectorXf alphaVec;  // hyperparameter
        VectorXf values;
        const unsigned ndist;

        ProbMixComps(const VectorXf &pis): ndist(pis.size()){  
            for (unsigned i = 0; i<ndist; ++i) {
                 //Parameter * pi = new Parameter("Pi");
                 this->push_back(new Parameter("Pi" + to_string(static_cast<long long>(i + 1))));
            }
            alphaVec.setOnes(pis.size());
            values = pis;
        }
        
        void sampleFromFC(const VectorXf snpStore, const VectorXf &pis);
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
    ProbMixComps Pis; 
    Gammas gamma;

    BayesR(const Data &data, const float varGenotypic, const float varResidual, const VectorXf pis, const VectorXf gamma, const bool estimatePi, 
           const string &algorithm, const bool message = true):
    BayesC(data, varGenotypic, varResidual, pis[0], estimatePi, "Gibbs", false),
    Pis(pis),
    gamma(gamma, vector<string>(gamma.size())),
    snpEffects(data.snpEffectNames, algorithm)
    {
        paramSetVec  = {&snpEffects, &fixedEffects};
        for (unsigned i=0; i<Pis.size(); ++i) { 
           Pis[i]->value=Pis.values[i];  
        }
        paramVec     = {&nnzSnp, &sigmaSq, &vare, &varg, &hsq};
        paramVec.insert(paramVec.begin(), Pis.begin(), Pis.end());
        paramToPrint = {&nnzSnp, &sigmaSq, &vare, &varg, &hsq, &rounding};
        paramToPrint.insert(paramToPrint.begin(), Pis.begin(), Pis.end());
        if (message && myMPI::rank==0) {
            string alg = algorithm;
            if (alg!="HMC") alg = "Gibbs (default)";
            cout << "\nBayesR model fitted. Algorithm: " << alg << "." << endl;
        }
    }   
    void sampleUnknowns(void);
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
        
        enum {random_walk, hmc} algorithm;
        
        AcceptanceRate ar;
        Parameter tuner;
        
        Sp(const unsigned m, const float var, const float start, const string &alg, const string &lab = "S"): Parameter(lab), mean(0), var(var), numSnps(m)
        , tuner(alg=="RMH" ? "varProp" : "Stepsize"){
            value = start;  // starting value
            varProp = 0.01;
            stepSize = 0.001;
            numSteps = 100;
            if (alg=="RMH") algorithm = random_walk;
            else algorithm = hmc;
        }
        
        // note that the scale factor of sigmaSq will be simultaneously updated
        void sampleFromFC(const float snpEffWtdSumSq, const unsigned numNonZeros, const float sigmaSq, const VectorXf &snpEffects,
                          const VectorXf &snp2pq, ArrayXf &snp2pqPowS, const ArrayXf &logSnp2pq,
                          const float vg, float &scale, float &sum2pqOneMinusS);
        void randomWalkMHsampler(const float snpEffWtdSumSq, const unsigned numNonZeros, const float sigmaSq, const VectorXf &snpEffects,
                                 const VectorXf &snp2pq, ArrayXf &snp2pqPowS, const ArrayXf &logSnp2pq,
                                 const float vg, float &scale, float &sum2pqOneMinusS);
        void hmcSampler(const unsigned numNonZeros, const float sigmaSq, const VectorXf &snpEffects,
                        const VectorXf &snp2pq, ArrayXf &snp2pqPowS, const ArrayXf &logSnp2pq,
                        const float vg, float &scale, float &sum2pqOneMinusS);
        float gradientU(const float S, const ArrayXf &snpEffects, const float snp2pqLogSum, const ArrayXf &snp2pq, const ArrayXf &logSnp2pq, const float sigmaSq, const float vg);
        float computeU(const float S, const ArrayXf &snpEffects, const float snp2pqLogSum, const ArrayXf &snp2pq, const ArrayXf &logSnp2pq, const float sigmaSq, const float vg, float &scale, float &U_chisq);
    };
    
    class SnpEffects : public BayesC::SnpEffects {
    public:
        float wtdSumSq;  // weighted sum of squares by 2pq^S
        float sum2pqOneMinusS;  // sum of delta_j* (2p_j q_j)^{1-S}
        
        SnpEffects(const vector<string> &header, const VectorXf &snp2pq, const float pi): BayesC::SnpEffects(header, "Gibbs") {
            wtdSumSq = 0.0;
            //sum2pqOneMinusS = 0.0;
            sum2pqOneMinusS = snp2pq.sum()*pi;  // starting value of S is 0
        }
        
        void sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag,
                          const float sigmaSq, const float pi, const float vare,
                          const ArrayXf &snp2pqPowS, const VectorXf &snp2pq,
                          const float vg, float &scale, VectorXf &ghat);
    };
    
    
public:
    float genVarPrior;
    float scalePrior;
    ArrayXf snp2pqPowS;
    const ArrayXf logSnp2pq;
    
    Sp S;
    SnpEffects snpEffects;
    
    BayesS(const Data &data, const float varGenotypic, const float varResidual, const float pival, const bool estimatePi, const float varS, const vector<float> &svalue,
           const string &algorithm, const bool message = true):
    BayesC(data, varGenotypic, varResidual, pival, estimatePi, "Gibbs", false),
    logSnp2pq(data.snp2pq.array().log()),
    S(data.numIncdSnps, varS, svalue[0], algorithm),
    snpEffects(data.snpEffectNames, data.snp2pq, pival),
    genVarPrior(varGenotypic),
    scalePrior(sigmaSq.scale)
    {
        findStartValueForS(svalue);
        snp2pqPowS = data.snp2pq.array().pow(S.value);
        sigmaSq.value = varGenotypic/((snp2pqPowS*data.snp2pq.array()).sum()*pival);
        scale.value = sigmaSq.scale = 0.5*sigmaSq.value;

        paramSetVec = {&snpEffects, &fixedEffects};
        paramVec = {&pi, &nnzSnp, &sigmaSq, &S, &vare, &varg, &hsq};
        paramToPrint = {&pi, &nnzSnp, &sigmaSq, &scale, &S, &vare, &varg, &hsq, &S.ar, &S.tuner, &rounding};
        if (message && myMPI::rank==0) {
            string alg = algorithm;
            if (alg!="RMH") alg = "HMC (default)";
            cout << "\nBayesS model fitted. Algorithm: " << alg << "." << endl;
            cout << "scale factor: " << sigmaSq.scale << endl;
        }
    }
    
    void sampleUnknowns(void);
    void findStartValueForS(const vector<float> &val);
    float computeLogLikelihood(void);
    void sampleUnknownsWarmup(void);
};


class BayesNS : public BayesS {
    // combine BayesN and BayesS primarily for speed
public:
    
    class SnpEffects : public BayesN::SnpEffects {
    public:
        float wtdSumSq;  // weighted sum of squares by 2pq^S
        float sum2pqOneMinusS;  // sum of delta_j* (2p_j q_j)^{1-S}
        
        ArrayXf varPseudoPrior;
        
        SnpEffects(const vector<string> &header, const VectorXi &windStart, const VectorXi &windSize,
                   const unsigned snpFittedPerWindow, const VectorXf &snp2pq, const float pi):
        BayesN::SnpEffects(header, windStart, windSize, snpFittedPerWindow){
            wtdSumSq = 0.0;
            sum2pqOneMinusS = 0.0;
            //sum2pqOneMinusS = snp2pq.sum()*(1.0f-pi)*(1.0f-snpFittedPerWindow/float(windSize));  // starting value of S is 0
            varPseudoPrior.setZero(size);
        }
        
        void sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag,
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
    
    BayesNS(const Data &data, const float varGenotypic, const float varResidual, const float pival,
            const bool estimatePi, const float varS, const vector<float> &svalue, const unsigned snpFittedPerWindow,
            const string &algorithm, const bool message = true):
    BayesS(data, varGenotypic, varResidual, pival, estimatePi, varS, svalue, algorithm, false),
    snpEffects(data.snpEffectNames, data.windStart, data.windSize, snpFittedPerWindow, data.snp2pq, pival),
    //S(data.numIncdSnps, "HMC"),
    sigmaSq(varGenotypic, data.snp2pq, pival, snpEffects.localPi, snpFittedPerWindow),
    scale(sigmaSq.scale),
    windDelta(vector<string>(snpEffects.numWindows))
    {
        paramSetVec = {&snpEffects, &fixedEffects, &windDelta};
        paramVec = {&pi, &nnzWind, &nnzSnp, &sigmaSq, &S, &vare, &varg, &hsq};
        paramToPrint = {&pi, &nnzWind, &nnzSnp, &sigmaSq, &scale, &S, &vare, &varg, &hsq, &S.ar, &S.tuner, &rounding};
        if (message && myMPI::rank==0) {
            string alg = algorithm;
            if (alg!="RMH") alg = "HMC (default)";
            cout << "\nBayesNS model fitted. Algorithm: " << alg << "." << endl;
            cout << "scale factor: " << sigmaSq.scale << endl;
        }
    }
    
    void sampleUnknowns(void);
};


class ApproxBayesC : public BayesC {
public:
    
    class FixedEffects : public BayesC::FixedEffects {
    public:
        FixedEffects(const vector<string> &header): BayesC::FixedEffects(header){}
        
        void sampleFromFC(const MatrixXf &XPX, const VectorXf &XPXdiag,
                          const MatrixXf &ZPX, const VectorXf &XPy,
                          const VectorXf &snpEffects, const float vare,
                          VectorXf &rcorr);
    };
    
    class SnpEffects : public BayesC::SnpEffects {
    public:
        float sum2pq;
        
        SnpEffects(const vector<string> &header): BayesC::SnpEffects(header, "Gibbs"){
            sum2pq = 0.0;
        }
        
        void sampleFromFC(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const VectorXf &snp2pq, const VectorXf &LDsamplVar,
                          const float sigmaSq, const float pi, const float vare, const float varg, const float ps, const float overdispersion);
        void sampleFromFC(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const VectorXf &snp2pq, const VectorXf &LDsamplVar,
                          const float sigmaSq, const float pi, const float vare, const float varg, const float ps, const float overdispersion);
        void hmcSampler(VectorXf &rcorr, const VectorXf &ZPy, const vector<VectorXf> &ZPZ,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const float sigmaSq, const float pi, const float vare);
        VectorXf gradientU(const VectorXf &effects, VectorXf &rcorr, const VectorXf &ZPy, const vector<VectorXf> &ZPZ,
                                                     const VectorXi &windStart, const VectorXi &windSize, const unsigned chrStart, const unsigned chrSize,
                                                     const float sigmaSq, const float vare);
        float computeU(const VectorXf &effects, const VectorXf &rcorr, const VectorXf &ZPy,                                             const float sigmaSq, const float vare);
    };
    
    class ResidualVar : public BayesC::ResidualVar {
    public:
        ResidualVar(const float vare, const unsigned nobs): BayesC::ResidualVar(vare, nobs){}
        
        //void sampleFromFC(VectorXf &rcorr, const SparseMatrix<float> &ZPZinv);
        void sampleFromFC(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr);
        void sampleFromFC(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr, const float hsq, const float phi);
        
        void sampleFromFC2(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &ghat);
        
        void randomWalkMHsampler(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr, const VectorXf &ZPZrss, const float sigmaSq, const float pi);
    };
    
    class GenotypicVar : public BayesC::GenotypicVar {
    public:
        const unsigned nobs;
        
        GenotypicVar(const float varg, const unsigned n): BayesC::GenotypicVar(varg), nobs(n){}
        void compute(const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr);
    };

    class Rounding : public BayesC::Rounding {
    public:
        void computeRcorr(const VectorXf &ZPy, const vector<SparseVector<float> > &ZPZsp,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXf &snpEffects, VectorXf &rcorr);
        void computeRcorr(const VectorXf &ZPy, const vector<VectorXf> &ZPZ,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXf &snpEffects, VectorXf &rcorr);
        void computeGhat(const MatrixXf &Z, const VectorXf &snpEffects, VectorXf &ghat);
    };
    
//    class Overdispersion : public BayesC::ResidualVar {
//    public:
//        Overdispersion(const float vare, const unsigned nobs): BayesC::ResidualVar(vare, nobs, "TauSq"){}
//        
//        void sampleFromFC(const VectorXf &y, const VectorXf &ghat);
//    };
    
    class PopulationStratification : public Parameter {
    public:
        PopulationStratification(): Parameter("PS"){}
        
        void compute(const VectorXf &rcorr, const VectorXf &ZPZdiag, const VectorXf &LDsamplVar, const float varg, const float vare);
    };
    
public:
    const Data &data;
    const float phi;   // the shrinkage parameter for heritability estimate
    const float overdispersion;
    
    VectorXf rcorr;
    VectorXf varei;   // residual variance specific to each snp
    
    bool sparse;
    bool modelPS;
    
    FixedEffects fixedEffects;
    SnpEffects snpEffects;
    BayesC::VarEffects sigmaSq;
    BayesC::Pi pi;
    ResidualVar vare;
    GenotypicVar varg;
//    BayesC::ResidualVar vare;
    Rounding rounding;
    varEffectScaled sigmaSqG;
//    Overdispersion tauSq;
    PopulationStratification ps;
    
    ApproxBayesC(const Data &data, const float varGenotypic, const float varResidual, const float pival, const bool estimatePi,
                 const float phi, const float overdispersion, const bool estimatePS, const bool message = true)
    : BayesC(data, varGenotypic, varResidual, pival, estimatePi, "Gibbs", false)
    , data(data)
    , rcorr(data.ZPy)
    , varei(data.tss.array()/data.n.array())
    , fixedEffects(data.fixedEffectNames)
    , snpEffects(data.snpEffectNames)
    , sigmaSq(varGenotypic, data.snp2pq, pival)
    , pi(pival)
    , vare(varResidual, data.numKeptInds)
    , varg(varGenotypic, data.numKeptInds)
//    , tauSq(varResidual, data.numKeptInds)
    , phi(phi)
    , overdispersion(overdispersion)
    {
        sparse = data.sparseLDM;
        modelPS = estimatePS;
        paramSetVec = {&snpEffects, &fixedEffects};
        paramVec = {&pi, &nnzSnp, &sigmaSq, &vare, &varg, &hsq, &sigmaSqG};
        paramToPrint = {&pi, &nnzSnp, &sigmaSq, &vare, &varg, &hsq, &sigmaSqG, &rounding};
        if (modelPS) {
            paramVec.push_back(&ps);
            paramToPrint.push_back(&ps);
        }
        if (message && myMPI::rank==0) {
            cout << "\nApproximate BayesC model fitted." << endl;
        }
    }
    
    void sampleUnknowns(void);
    static void ldScoreReg(const VectorXf &chisq, const VectorXf &LDscore, const VectorXf &LDsamplVar,
                           const float varg, const float vare, float &ps, float &vargj);
};



class ApproxBayesS : public BayesS {
public:
    
    class SnpEffects : public ApproxBayesC::SnpEffects {
    public:
        float wtdSumSq;  // weighted sum of squares by 2pq^S
        float sum2pqOneMinusS;  // sum of delta_j* (2p_j q_j)^{1-S}
        
        SnpEffects(const vector<string> &header, const VectorXf &snp2pq, const float pi): ApproxBayesC::SnpEffects(header) {
            wtdSumSq = 0.0;
            sum2pqOneMinusS = snp2pq.sum()*pi;  // starting value of S is 0
        }
        
        void sampleFromFC(VectorXf &rcorr,const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const float sigmaSq, const float pi, const float vare,
                          const VectorXf &snp2pqPowS, const VectorXf &snp2pq, const VectorXf &LDsamplVar,
                          const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n,
                          const float varg, const float ps, const float overdispersion);
        void sampleFromFC(VectorXf &rcorr,const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const float sigmaSq, const float pi, const float vare,
                          const VectorXf &snp2pqPowS, const VectorXf &snp2pq, const VectorXf &LDsamplVar,
                          const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n,
                          const float varg, const float ps, const float overdispersion);

        void sampleFromFC(const VectorXf &ZPy,const MatrixXf &Z, const VectorXf &ZPZdiag,
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
    };
    

public:
    VectorXf rcorr;
    VectorXf varei;   // residual variance specific to each snp
    
    VectorXf vareiMean;  ///TMP
    
    const float phi;   // the shrinkage parameter for heritability estimate
    const float overdispersion;

    bool sparse;
    bool modelPS;
    
    SnpEffects snpEffects;
    ApproxBayesC::FixedEffects fixedEffects;
    ApproxBayesC::ResidualVar vare;
    ApproxBayesC::GenotypicVar varg;
    ApproxBayesC::Rounding rounding;
    varEffectScaled sigmaSqG;
    ApproxBayesC::PopulationStratification ps;
    
//    ApproxBayesC::Overdispersion tauSq;
    
    ApproxBayesS(const Data &data, const float varGenotypic, const float varResidual, const float pival, const bool estimatePi,
                 const float phi, const float overdispersion, const bool estimatePS,
                 const float varS, const vector<float> &svalue,
                 const string &algorithm, const bool message = true)
    : BayesS(data, varGenotypic, varResidual, pival, estimatePi, varS, svalue, algorithm, false)
    , rcorr(data.ZPy)
    , varei(data.tss.array()/data.n.array())
    , snpEffects(data.snpEffectNames, data.snp2pq, pival)
    , fixedEffects(data.fixedEffectNames)
    , vare(varResidual, data.numKeptInds)
    , varg(varGenotypic, data.numKeptInds)
//    , tauSq(varResidual, data.numKeptInds)
    , phi(phi)
    , overdispersion(overdispersion)
    {
        ghat.setZero(data.Z.rows());
        sparse = data.sparseLDM;
        modelPS = estimatePS;
        paramSetVec = {&snpEffects, &fixedEffects};
        paramVec = {&pi, &nnzSnp, &sigmaSq, &S, &vare, &varg, &hsq, &sigmaSqG};
        paramToPrint = {&pi, &nnzSnp, &sigmaSq, &S, &vare, &varg, &hsq, &sigmaSqG, &S.ar, &S.tuner, &rounding};
        if (modelPS) {
            paramVec.push_back(&ps);
            paramToPrint.push_back(&ps);
        }
        if (message && myMPI::rank==0) {
            string alg = algorithm;
            if (alg!="RMH") alg = "HMC (default)";
            cout << "\nApproximate BayesS model fitted. Algorithm: " << alg << "." << endl;
        }

    }
    
    void sampleUnknowns(void);
};

// -----------------------------------------------------------------------------------------------
// Approximate Bayes R
// -----------------------------------------------------------------------------------------------

class ApproxBayesR : public ApproxBayesC {
    
public:
    
    class SnpEffects : public ApproxBayesC::SnpEffects {
    public:
        float sum2pq;
        
        SnpEffects(const vector<string> &header): ApproxBayesC::SnpEffects(header){
            sum2pq = 0.0;
            
        }
        
        void sampleFromFC(VectorXf &rcorr, const vector<SparseVector<float>> &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const VectorXf &snp2pq,
                          const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare, VectorXf &snpStore);
        void sampleFromFC(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const VectorXf &snp2pq,
                          const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare, VectorXf &snpStore);
    };
    

    class ProbMixComps : public vector<Parameter*>, public Stat::Dirichlet {

        // prior probability of a snp being in any of the distributions effect has a dirichlet prior
    public:
        VectorXf alphaVec;  // hyperparameter
        VectorXf values;
        const unsigned ndist;

        ProbMixComps(const VectorXf &pis): ndist(pis.size()){  
            for (unsigned i = 0; i<ndist; ++i) {
                 //Parameter * pi = new Parameter("Pi");
                 this->push_back(new Parameter("Pi" + to_string(static_cast<long long>(i + 1))));
            }
            alphaVec.setOnes(pis.size());
            values = pis;
        }
        
        void sampleFromFC(const VectorXf snpStore, const VectorXf &pis);
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
    ProbMixComps Pis; 
    Gammas gamma;
    ApproxBayesC::GenotypicVar varg;
    
    ApproxBayesR(const Data &data, const float varGenotypic, const float varResidual, const VectorXf pis, const VectorXf gamma, const bool estimatePi, 
                 const bool message = true):
    ApproxBayesC(data, varGenotypic, varResidual, pis[0], estimatePi, 0, 0, false, false),
    Pis(pis),
    gamma(gamma, vector<string>(gamma.size())),
    varg(varGenotypic, data.numKeptInds),
    snpEffects(data.snpEffectNames)
    {
        sparse = data.sparseLDM;
        // varg.value = varGenotypic; //// NOTE: write it into constructor!!!
        paramSetVec = {&snpEffects, &fixedEffects};
        for (unsigned i=0; i<Pis.size(); ++i) { 
           Pis[i]->value=Pis.values[i];  
        }
        paramVec     = {&nnzSnp, &sigmaSq, &vare, &varg, &hsq};
        paramVec.insert(paramVec.begin(), Pis.begin(), Pis.end());
        paramToPrint = {&nnzSnp, &sigmaSq, &vare, &varg, &hsq, &rounding};
        paramToPrint.insert(paramToPrint.begin(), Pis.begin(), Pis.end());
        if (message && myMPI::rank==0) {
            cout << "\nApproximate BayesR model fitted." << endl;
        }
    }
    
    void sampleUnknowns(void);
};

// -----------------------------------------------------------------------------------------------
// Approximate Bayes Kappa
// -----------------------------------------------------------------------------------------------

class ApproxBayesKappa : public ApproxBayesC {
    
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
    
    class Kappa : public Parameter, public Stat::Normal {
        // S parameter for genotypes or equivalently for the variance of snp effects
        
        // random-walk MH and HMC algorithms implemented
        
    public:
        const float k0;  // prior value of gamma distribution
        const float theta0;   // prior  value of gamma distribution
        // const unsigned numSnps;
        
        float varProp;     // variance of proposal normal for random walk MH
        
        
        AcceptanceRate ar;
        Parameter tuner;
        
        Kappa(const float start, const string &lab = "Kappa"): Parameter(lab), k0(2.651), theta0(0.86858896), tuner("varProp"){
            value = start;  // starting value
            varProp = 0.1;
        }

        
        void randomWalkMHsampler(const float sigmaSq, const VectorXf &snpEffects, const VectorXf &snpindist);
    };

    
    class SnpEffects : public ApproxBayesC::SnpEffects {
    public:
        float sum2pq;
        
        SnpEffects(const vector<string> &header): ApproxBayesC::SnpEffects(header){
            sum2pq = 0.0;
            
        }
        
        void sampleFromFC(VectorXf &rcorr, const vector<SparseVector<float>> &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const VectorXf &snp2pq,
                          const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare, VectorXf &snpStore, const float kappa, VectorXf &snpindist);
        void sampleFromFC(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const VectorXf &snp2pq,
                          const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare, VectorXf &snpStore, const float kappa, VectorXf &snpindist);
    };
    

    class ProbMixComps : public vector<Parameter*>, public Stat::Dirichlet {

        // prior probability of a snp being in any of the distributions effect has a dirichlet prior
    public:
        VectorXf alphaVec;  // hyperparameter
        VectorXf values;
        const unsigned ndist;

        ProbMixComps(const VectorXf &pis): ndist(pis.size()){  
            for (unsigned i = 0; i<ndist; ++i) {
                 //Parameter * pi = new Parameter("Pi");
                 this->push_back(new Parameter("Pi" + to_string(static_cast<long long>(i + 1))));
            }
            alphaVec.setOnes(pis.size());
            values = pis;
        }
        
        void sampleFromFC(const VectorXf snpStore, const VectorXf &pis);
    };

    class Gammas : public ParamSet {
        // Set of scaling factors for each of the distributions
    public:
        Gammas(const VectorXf &gamma, const vector<string> &header, const string &lab = "gamma"): ParamSet(lab, header){
            values = gamma;
        }
    };

    class SnpIndist : public ParamSet {
        // Set of scaling factors for each of the distributions
        public:
          SnpIndist(const vector<string> &header, const string &lab = "snpindist"): ParamSet(lab, header){
       }
    };

public:
    
    VectorXf snpStore;   
    SnpIndist snpindist;
    SnpEffects snpEffects;
    ProbMixComps Pis; 
    Gammas gamma;
    ApproxBayesC::GenotypicVar varg;
    Kappa kappa;
    
    ApproxBayesKappa(const Data &data, const float varGenotypic, const float varResidual, const VectorXf pis, const VectorXf gamma, const bool estimatePi, 
                     const float kappa_str, const bool message = true):
    ApproxBayesC(data, varGenotypic, varResidual, pis[0], estimatePi, 0, 0, false),
    Pis(pis),
    gamma(gamma, vector<string>(gamma.size())),
    kappa(kappa_str),
    varg(varGenotypic, data.numKeptInds),
    snpindist(data.snpEffectNames), 
    snpEffects(data.snpEffectNames)
    {
        sparse = data.sparseLDM;
        // varg.value = varGenotypic; //// NOTE: write it into constructor!!!
        paramSetVec = {&snpEffects, &fixedEffects};
        for (unsigned i=0; i<Pis.size(); ++i) { 
           Pis[i]->value=Pis.values[i];  
        }
        paramVec     = {&nnzSnp, &sigmaSq, &vare, &varg, &hsq, &kappa};
        paramVec.insert(paramVec.begin(), Pis.begin(), Pis.end());
        paramToPrint = {&nnzSnp, &sigmaSq, &vare, &varg, &hsq, &kappa, &rounding};
        paramToPrint.insert(paramToPrint.begin(), Pis.begin(), Pis.end());
        if (message && myMPI::rank==0) {
            cout << "\nApproximate Bayes Kappa model fitted." << endl;
        }
    }
    
    void sampleUnknowns(void);
};



#endif /* model_hpp */




