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
        float scale;     // hyperparameter
        
        VarEffects(const float vg, const VectorXf &snp2pq, const float pi, const bool noscale, const string &lab = "SigmaSq")
        : Parameter(lab), df(4)
        {
            // cout << "To scale or not to scale " << noscale << endl;
            // cout << "Scale value 1 " << value << endl;
            if (myMPI::partition == "bycol" && noscale == true) {
                int sizeFull;
                MPI_Allreduce(&myMPI::iSize, &sizeFull, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                VectorXf snp2pqFull(sizeFull);
                MPI_Allgatherv((float*)&snp2pq[0], myMPI::iSize, MPI_FLOAT, &snp2pqFull[0], &myMPI::srcounts[0], &myMPI::displs[0], MPI_FLOAT, MPI_COMM_WORLD);
                value = vg/(snp2pqFull.sum()*pi);  // derived from prior knowledge on Vg and pi
            } else if (myMPI::partition == "bycol" && noscale == false) {
                int sizeFull;
                MPI_Allreduce(&myMPI::iSize, &sizeFull, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                VectorXf snp2pqFull(sizeFull);
                MPI_Allgatherv((float*)&snp2pq[0], myMPI::iSize, MPI_FLOAT, &snp2pqFull[0], &myMPI::srcounts[0], &myMPI::displs[0], MPI_FLOAT, MPI_COMM_WORLD);
                value = vg/(snp2pq.size() * pi);  // derived from prior knowledge on Vg and pi
            }
            else if (noscale == true) {
                value = vg / (snp2pq.sum() * pi);  // derived from prior knowledge on Vg and pi
            } else {
                value = vg / (snp2pq.size() * pi);  // derived from prior knowledge on Vg and pi
            }
            scale = 0.5f*value;  // due to df = 4
            
//            cout << value << " " << vg << " " << snp2pq.sum() << " " << pi << endl;
        }
        
        void sampleFromFC(const float snpEffSumSq, const unsigned numSnpEff);
        void sampleFromPrior(void);
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
    
    BayesC(const Data &data, const float varGenotypic, const float varResidual, const float pival, const float piAlpha, const float piBeta, const bool estimatePi, const bool noscale,
           const string &algorithm = "Gibbs", const bool message = true):
    data(data),
    ycorr(data.y),
    fixedEffects(data.fixedEffectNames),
    snpEffects(data.snpEffectNames, algorithm),
    sigmaSq(varGenotypic, data.snp2pq, pival, noscale),
    scale(sigmaSq.scale),
    pi(pival, piAlpha, piBeta),
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
        
        void sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag,
                          const VectorXf &sigmaSq, const float pi, const float vare, VectorXf &ghat);
    };

    class VarEffects : public ParamSet, public BayesC::VarEffects {
    public:
        VarEffects(const float vg, const VectorXf &snp2pq, const float pi):
        ParamSet("SigmaSqs", vector<string>(snp2pq.size())),
        BayesC::VarEffects(vg, snp2pq, pi, false){
            values.setConstant(size, value);
        }
        
        void sampleFromFC(const VectorXf &betaSq);
    };
    
    SnpEffects snpEffects;
    VarEffects sigmaSq;

    BayesB(const Data &data, const float varGenotypic, const float varResidual, const float pival, const float piAlpha, const float piBeta,
           const bool estimatePi, const bool noscale, const bool message = true):
    BayesC(data, varGenotypic, varResidual, pival, piAlpha, piBeta, estimatePi, noscale, "Gibbs", false),
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
    
    BayesN(const Data &data, const float varGenotypic, const float varResidual, const float pival, const float piAlpha, const float piBeta,
           const bool estimatePi, const bool noscale, const unsigned snpFittedPerWindow, const bool message = true):
    BayesC(data, varGenotypic, varResidual, pival, piAlpha, piBeta, estimatePi, noscale, "Gibbs", false),
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
        
        void sampleFromFC(const VectorXf snpStore);
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

    BayesR(const Data &data, const float varGenotypic, const float varResidual, const VectorXf pis, const float piAlpha, const float piBeta, const VectorXf gamma, const bool estimatePi, const bool noscale,
           const string &algorithm, const bool message = true):
    BayesC(data, varGenotypic, varResidual, pis[0], piAlpha, piBeta, estimatePi, noscale, "Gibbs", false),
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
            cout << "scale factor: " << sigmaSq.scale << endl;
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
                          const float vg, float &scale, float &sum2pqSplusOne);
        void sampleFromPrior(void);
        void randomWalkMHsampler(const float snpEffWtdSumSq, const unsigned numNonZeros, const float sigmaSq, const VectorXf &snpEffects,
                                 const VectorXf &snp2pq, ArrayXf &snp2pqPowS, const ArrayXf &logSnp2pq,
                                 const float vg, float &scale, float &sum2pqSplusOne);
        void hmcSampler(const unsigned numNonZeros, const float sigmaSq, const VectorXf &snpEffects,
                        const VectorXf &snp2pq, ArrayXf &snp2pqPowS, const ArrayXf &logSnp2pq,
                        const float vg, float &scale, float &sum2pqSplusOne);
        float gradientU(const float S, const ArrayXf &snpEffects, const float snp2pqLogSum, const ArrayXf &snp2pq, const ArrayXf &logSnp2pq, const float sigmaSq, const float vg);
        float computeU(const float S, const ArrayXf &snpEffects, const float snp2pqLogSum, const ArrayXf &snp2pq, const ArrayXf &logSnp2pq, const float sigmaSq, const float vg, float &scale, float &U_chisq);
    };
    
    class SnpEffects : public BayesC::SnpEffects {
    public:
        float wtdSumSq;  // weighted sum of squares by 2pq^S
        float sum2pqSplusOne;  // sum of delta_j* (2p_j q_j)^{1+S}
        
        SnpEffects(const vector<string> &header, const VectorXf &snp2pq, const float pi): BayesC::SnpEffects(header, "Gibbs") {
            wtdSumSq = 0.0;
            sum2pqSplusOne = snp2pq.sum()*pi;  // starting value of S is 0
        }
        
        void sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag,
                          const float sigmaSq, const float pi, const float vare,
                          const ArrayXf &snp2pqPowS, const VectorXf &snp2pq,
                          const float vg, float &scale, VectorXf &ghat);
    };
    
    
public:
    unsigned iter;
    
    float genVarPrior;
    float scalePrior;

    ArrayXf snp2pqPowS;
    const ArrayXf logSnp2pq;
    
    Sp S;
    SnpEffects snpEffects;
    
    BayesS(const Data &data, const float varGenotypic, const float varResidual, const float pival, const float piAlpha, const float piBeta, const bool estimatePi, const float varS, const vector<float> &svalue,
           const string &algorithm, const bool message = true):
    BayesC(data, varGenotypic, varResidual, pival, piAlpha, piBeta, estimatePi, true, "Gibbs", false),
    logSnp2pq(data.snp2pq.array().log()),
    S(data.numIncdSnps, varS, svalue[0], algorithm),
    snpEffects(data.snpEffectNames, data.snp2pq, pival),
    genVarPrior(varGenotypic),
    scalePrior(sigmaSq.scale)
    {
        iter = 0;
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
        unsigned iter;
        unsigned burnin;
        
        float wtdSumSq;  // weighted sum of squares by 2pq^S
        float sum2pqSplusOne;  // sum of delta_j* (2p_j q_j)^{1+S}
        
        ArrayXf varPseudoPrior;
        
        SnpEffects(const vector<string> &header, const VectorXi &windStart, const VectorXi &windSize,
                   const unsigned snpFittedPerWindow, const VectorXf &snp2pq, const float pi):
        BayesN::SnpEffects(header, windStart, windSize, snpFittedPerWindow){
            iter - 0;
            burnin = 2000;
            wtdSumSq = 0.0;
            sum2pqSplusOne = 0.0;
            //sum2pqSplusOne = snp2pq.sum()*(1.0f-pi)*(1.0f-snpFittedPerWindow/float(windSize));  // starting value of S is 0
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
    
    BayesNS(const Data &data, const float varGenotypic, const float varResidual, const float pival, const float piAlpha, const float piBeta,
            const bool estimatePi, const float varS, const vector<float> &svalue, const unsigned snpFittedPerWindow,
            const string &algorithm, const bool message = true):
    BayesS(data, varGenotypic, varResidual, pival, piAlpha, piBeta, estimatePi, varS, svalue, algorithm, false),
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
        
        VectorXf nnzPerChr;
        VectorXi leaveout;
        
        SnpEffects(const vector<string> &header): BayesC::SnpEffects(header, "Gibbs"){
            sum2pq = 0.0;
            leaveout.setZero(size);
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
        float computeU(const VectorXf &effects, const VectorXf &rcorr, const VectorXf &ZPy, const float sigmaSq, const float vare);
    };
    
    class ResidualVar : public BayesC::ResidualVar {
    public:
        const float icrsq;
        
        ResidualVar(const float vare, const unsigned nobs, const float icrsq): BayesC::ResidualVar(vare, nobs), icrsq(icrsq) {}
        
        //void sampleFromFC(VectorXf &rcorr, const SparseMatrix<float> &ZPZinv);
//        void sampleFromFC(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr, const float varg, const float nnz);
        void sampleFromFC(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr, const float covg);
//        void sampleFromFCshrink(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr, const float hsq, const float phi);
        
//        void sampleFromFC2(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &ghat);
        
//        void randomWalkMHsampler(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr, const VectorXf &ZPZrss, const float sigmaSq, const float pi);
    };
    
    class GenotypicVar : public BayesC::GenotypicVar {
    public:
        const unsigned nobs;
        
        GenotypicVar(const float varg, const unsigned n): BayesC::GenotypicVar(varg), nobs(n){}
//        void compute(const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr);
        void compute(const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr, const float covg);
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
    
    class PopulationStratification : public Parameter, public Stat::InvChiSq {
    public:
        const float df;
        const float scale;
        
        VectorXf chrSpecific;
        
        PopulationStratification(): Parameter("PS"), df(4), scale(0.5){
            chrSpecific.setZero(22);
        }
        
        void compute(const VectorXf &rcorr, const VectorXf &ZPZdiag, const VectorXf &LDsamplVar, const float varg, const float vare, const VectorXf &chisq);
        void compute(const VectorXf &rcorr, const VectorXf &ZPZdiag, const VectorXf &LDsamplVar, const float varg, const float vare, const vector<ChromInfo*> chromInfoVec);
    };
    
    class NumResidualOutlier : public Parameter {
    public:
        ofstream out;
        unsigned iter;
        
        NumResidualOutlier(): Parameter("Nro"){
            iter = 0;
        }
        
        void compute(const VectorXf &rcorr, const VectorXf &ZPZdiag, const VectorXf &LDsamplVar, const float varg, const float vare, const vector<string> &snpName, VectorXi &leaveout, const vector<SparseVector<float> > &ZPZ, const VectorXf &ZPy);
    };
    
    class InterChrGenetCov : public Parameter {
    public:
        const float spouseCorrelation;
        const unsigned nobs;
        
        InterChrGenetCov(const float corr, const unsigned nobs): Parameter("GenCov"), spouseCorrelation(corr), nobs(nobs) {}
        
        void compute(const float ypy, const VectorXf &effects, const VectorXf &ZPy, const VectorXf &rcorr);
    };
    
    class NnzGwas : public Parameter {
    public:
        unsigned iter;
        
        NnzGwas(): Parameter("NnzGwas"){
            iter = 0;
        }
        
        void compute(const VectorXf &effects, const vector<SparseVector<float> > &ZPZ, const VectorXf &ZPZdiag);
    };
    
    class PiGwas : public Parameter {
    public:
        unsigned iter;
        PiGwas(): Parameter("PiGwas"){
            iter = 0;
        }
        
        void compute(const float nnzGwas, const unsigned numSnps);
    };
    
public:
    const Data &data;
    const float phi;   // the shrinkage parameter for heritability estimate
    const float overdispersion;
    
    VectorXf rcorr;
    VectorXf varei;   // residual variance specific to each snp
    
    bool sparse;
    bool modelPS;
    bool diagnose;
    
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
    NumResidualOutlier nro;
    InterChrGenetCov covg;
    PiGwas pigwas;
    NnzGwas nnzgwas;
    float genVarPrior;
    float scalePrior;
    bool noscale;
   
    ApproxBayesC(const Data &data, const float varGenotypic, const float varResidual, const float pival, const float piAlpha, const float piBeta, const bool estimatePi, const bool noscale,
                 const float phi, const float overdispersion, const bool estimatePS, const float icrsq, const float spouseCorrelation,
                 const bool diagnosticMode, const bool randomStart = false, const bool message = true)
    : BayesC(data, varGenotypic, varResidual, pival, piAlpha, piBeta, estimatePi, noscale, "Gibbs", false)
    , data(data)
    , rcorr(data.ZPy)
    , varei(data.tss.array()/data.n.array())
    , fixedEffects(data.fixedEffectNames)
    , snpEffects(data.snpEffectNames)
    , sigmaSq(varGenotypic, data.snp2pq, pival, noscale)
    , pi(pival, piAlpha, piBeta)
    , genVarPrior(varGenotypic)
    , noscale(noscale)
    , scalePrior(sigmaSq.scale)
    , vare(varResidual, data.numKeptInds, icrsq)
    , varg(varGenotypic, data.numKeptInds)
//    , tauSq(varResidual, data.numKeptInds)
    , phi(phi)
    , overdispersion(overdispersion)
    , covg(spouseCorrelation, data.numKeptInds)
    {
        sparse = data.sparseLDM;
        modelPS = estimatePS;
        diagnose = diagnosticMode;
        paramSetVec = {&snpEffects};
        paramVec = {&pi, &nnzSnp, &sigmaSq, &vare, &varg, &sigmaSqG, &hsq};
        paramToPrint = {&pi, &nnzSnp, &sigmaSq, &vare, &varg, &sigmaSqG, &hsq, &rounding};
//        if (sparse) {
//            paramVec.push_back(&pigwas);
//            paramVec.push_back(&nnzgwas);
//            paramToPrint.push_back(&pigwas);
//            paramToPrint.push_back(&nnzgwas);
//        }
        if (modelPS) {
            paramVec.push_back(&ps);
            paramToPrint.push_back(&ps);
        }
        if (diagnose) {
            nro.out.open((data.label+".diagnostics").c_str());
            paramVec.push_back(&nro);
            paramToPrint.push_back(&nro);
        }
        if (spouseCorrelation) {
            paramVec.push_back(&covg);
            paramToPrint.push_back(&covg);
        }
        if (message && myMPI::rank==0) {
            cout << "\nApproximate BayesC model fitted." << endl;
            cout << "scale factor: " << sigmaSq.scale << endl;
            if (noscale)
            {
               cout << "Fitting model assuming unscaled genotypes " << endl; 
            } else
            {
               cout << "Fitting model assuming scaled genotypes "  << endl;
            }
        }
        if (randomStart) sampleStartVal();
    }
    
    void sampleUnknowns(void);
    static void ldScoreReg(const VectorXf &chisq, const VectorXf &LDscore, const VectorXf &LDsamplVar,
                           const float varg, const float vare, float &ps);
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
    VectorXf varei;   // residual variance specific to each snp
    
    VectorXf vareiMean;  ///TMP
    
    const float phi;   // the shrinkage parameter for heritability estimate
    const float overdispersion;

    bool sparse;
    bool modelPS;
    bool diagnose;
    bool estimateEffectMean;

    SnpEffects snpEffects;
    ApproxBayesC::FixedEffects fixedEffects;
    ApproxBayesC::ResidualVar vare;
    ApproxBayesC::GenotypicVar varg;
    ApproxBayesC::Rounding rounding;
    varEffectScaled sigmaSqG;
    ApproxBayesC::PopulationStratification ps;
    ApproxBayesC::NumResidualOutlier nro;
    ApproxBayesC::InterChrGenetCov covg;
    ApproxBayesC::PiGwas pigwas;
    ApproxBayesC::NnzGwas nnzgwas;
    
//    ApproxBayesC::Overdispersion tauSq;
    
    MeanEffects mu;
    Smu Su;
    
    ApproxBayesS(const Data &data, const float varGenotypic, const float varResidual, const float pival, const float piAlpha, const float piBeta, const bool estimatePi,
                 const float phi, const float overdispersion, const bool estimatePS, const float icrsq, const float spouseCorrelation,
                 const float varS, const vector<float> &svalue,
                 const string &algorithm, const bool diagnosticMode, const bool randomStart = false, const bool message = true)
    : BayesS(data, varGenotypic, varResidual, pival, piAlpha, piBeta, estimatePi, varS, svalue, algorithm, false)
    , rcorr(data.ZPy)
    , varei(data.tss.array()/data.n.array())
    , snpEffects(data.snpEffectNames, data.snp2pq, pival)
    , fixedEffects(data.fixedEffectNames)
    , vare(varResidual, data.numKeptInds, icrsq)
    , varg(varGenotypic, data.numKeptInds)
//    , tauSq(varResidual, data.numKeptInds)
    , phi(phi)
    , overdispersion(overdispersion)
    , covg(spouseCorrelation, data.numKeptInds)
    , mu(data.numIncdSnps)
    {
        ghat.setZero(data.Z.rows());
        sparse = data.sparseLDM;
        modelPS = estimatePS;
        diagnose = diagnosticMode;
        
        estimateEffectMean = false;
        
        paramSetVec = {&snpEffects};
        paramVec = {&pi, &nnzSnp, &sigmaSq, &S, &vare, &varg, &sigmaSqG, &hsq};
        paramToPrint = {&pi, &nnzSnp, &sigmaSq, &S, &vare, &varg, &sigmaSqG, &hsq, &rounding};
//        if (sparse) {
//            paramVec.push_back(&pigwas);
//            paramVec.push_back(&nnzgwas);
//            paramToPrint.push_back(&pigwas);
//            paramToPrint.push_back(&nnzgwas);
//        }
        if (estimateEffectMean) {
            paramVec.push_back(&mu);
            paramVec.push_back(&Su);
            paramToPrint.push_back(&mu);
            paramToPrint.push_back(&Su);
        }
        if (modelPS) {
            paramVec.push_back(&ps);
            paramToPrint.push_back(&ps);
        }
        if (diagnose) {
            nro.out.open((data.label+".diagnostics").c_str());
            paramVec.push_back(&nro);
            paramToPrint.push_back(&nro);
        }
        if (spouseCorrelation) {
            paramVec.push_back(&covg);
            paramToPrint.push_back(&covg);
        }
        if (message && myMPI::rank==0) {
            string alg = algorithm;
            if (alg!="RMH") alg = "HMC (default)";
            cout << "\nApproximate BayesS model fitted. Algorithm: " << alg << "." << endl;
            cout << "scale factor: " << sigmaSq.scale << endl;
        }

        if (randomStart) sampleStartVal();

//        MatrixXf X(data.numIncdSnps, 2);
//        X.col(0) = VectorXf::Ones(data.numIncdSnps);
//        X.col(1) = data.LDscore;
//        VectorXf b = X.householderQr().solve(data.chisq);
//        ps.value = b[0] - 1.0;

    }
    
    void sampleUnknowns(void);
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
                          const vector<ChromInfo*> &chromInfoVec, const VectorXf &LDsamplVar, const ArrayXf &hSlT, const VectorXf &snp2pq,
                          const float sigmaSq, const float pi, const float vare, const float varg,
                          const float ps, const float overdispersion);
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
    
    ApproxBayesST(const Data &data, const float varGenotypic, const float varResidual, const float pival,
                  const float piAlpha, const float piBeta, const bool estimatePi, const float overdispersion,
                  const bool estimatePS, const float varS, const vector<float> &svalue, const bool estimateS,
                  const bool randomStart = false, const bool message = true):
    ApproxBayesS(data, varGenotypic, varResidual, pival, piAlpha, piBeta, estimatePi, 0, overdispersion, estimatePS, 0, 0, varS, svalue, "HMC", false, false),
    estimateS(estimateS),
    logLdsc(data.LDscore.array().log()),
    hSlT(snp2pqPowS),
    snpEffects(data.snpEffectNames, data.snp2pq, pival),
    S(data.numIncdSnps),
    T(data.numIncdSnps)
    {
        paramSetVec = {&snpEffects};
        paramVec = {&pi, &nnzSnp, &sigmaSq, &S, &T, &vare, &varg, &hsq};
        paramToPrint = {&pi, &nnzSnp, &sigmaSq, &S, &T, &vare, &varg, &hsq, &rounding};
        if (modelPS) {
            paramVec.push_back(&ps);
            paramToPrint.push_back(&ps);
        }
//        if (estimateS) {
//            paramToPrint.push_back(&S.ar);
//        } else {
//            paramToPrint.push_back(&T.ar);
//        }
        if (message && myMPI::rank==0) {
            cout << "\nApproximate BayesST model fitted." << endl;
        }
        
        if (randomStart) sampleStartVal();
    }

    void sampleUnknowns(void);
    void sampleStartVal(void);
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
                          const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const VectorXf &snp2pq, const VectorXf &LDsamplVar,
                          const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare, VectorXf &snpStore, 
                          const float varg, const float ps, const float overdispersion);
        void sampleFromFC(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXf &se, const VectorXf &tss, VectorXf &varei, const VectorXf &n, const VectorXf &snp2pq, const VectorXf &LDsamplVar,
                          const float sigmaSq, const VectorXf &pis, const VectorXf &gamma, const float vare, VectorXf &snpStore,
                          const float varg, const float ps, const float overdispersion);
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
        
        void sampleFromFC(const VectorXf snpStore);
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
    ApproxBayesC::FixedEffects fixedEffects;
    ApproxBayesC::ResidualVar vare;
    ApproxBayesC::GenotypicVar varg;
    ApproxBayesC::Rounding rounding;
    varEffectScaled sigmaSqG;
    ApproxBayesC::PopulationStratification ps;
    ApproxBayesC::NumResidualOutlier nro;
    ApproxBayesC::InterChrGenetCov covg;
    ApproxBayesC::PiGwas pigwas;
    ApproxBayesC::NnzGwas nnzgwas;
    ProbMixComps Pis; 
    Gammas gamma;
    float genVarPrior;
    float scalePrior;
    bool noscale;    

    const float overdispersion;

    
    ApproxBayesR(const Data &data, const float varGenotypic, const float varResidual, const VectorXf pis, const float piAlpha, const float piBeta, const VectorXf gamma, const bool estimatePi, const bool noscale, const float overdispersion, const bool estimatePS, const float spouseCorrelation, const bool diagnosticMode,
                 const bool message = true):
    ApproxBayesC(data, varGenotypic, varResidual, (1-pis[0]), piAlpha, piBeta, estimatePi, noscale, 0, overdispersion, estimatePS, 0, spouseCorrelation, false, false, false),
    Pis(pis),
    gamma(gamma, vector<string>(gamma.size())),
    vare(varResidual, data.numKeptInds, 0),
    varg(varGenotypic, data.numKeptInds),
    fixedEffects(data.fixedEffectNames),
    snpEffects(data.snpEffectNames),
    genVarPrior(varGenotypic),
    noscale(noscale),
    overdispersion(overdispersion),
    covg(spouseCorrelation, data.numKeptInds),
    scalePrior(sigmaSq.scale)
    {
        sparse = data.sparseLDM;
        // varg.value = varGenotypic; //// NOTE: write it into constructor!!!
        paramSetVec = {&snpEffects, &fixedEffects};
        // sigmaSq.value = varGenotypic/(data.snp2pq.array().sum()*(1-pis[0]));
        // scale.value = sigmaSq.scale = 0.5*sigmaSq.value;
        for (unsigned i=0; i<Pis.size(); ++i) { 
           Pis[i]->value=Pis.values[i];  
        }
        paramVec     = {&nnzSnp, &sigmaSq, &vare, &varg, &hsq};
        paramVec.insert(paramVec.begin(), Pis.begin(), Pis.end());
        paramToPrint = {&nnzSnp, &sigmaSq, &vare, &varg, &hsq, &rounding};
        paramToPrint.insert(paramToPrint.begin(), Pis.begin(), Pis.end());
        if (message && myMPI::rank==0) {
            cout << "\nApproximate BayesR model fitted." << endl;
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
    ApproxBayesR::ProbMixComps Pis;
    Gammas gamma;
    ApproxBayesC::ResidualVar vare;
    ApproxBayesC::GenotypicVar varg;
    Kappa kappa;
    float genVarPrior;
    float scalePrior;
    bool noscale;
    
    ApproxBayesKappa(const Data &data, const float varGenotypic, const float varResidual, const VectorXf pis, const float piAlpha, const float piBeta, const VectorXf gamma, const bool estimatePi, const bool noscale, const float icrsq,
                     const float kappa, const bool message = true):
    ApproxBayesC(data, varGenotypic, varResidual, (1-pis[(gamma.size()-1)]), piAlpha, piBeta, estimatePi, noscale, 0, 0, false, icrsq, 0, false, false, false),
    Pis(pis),
    gamma(gamma, vector<string>(gamma.size())),
    kappa(kappa),
    vare(varResidual, data.numKeptInds, 0),
    varg(varGenotypic, data.numKeptInds),
    snpindist(data.snpEffectNames), 
    genVarPrior(varGenotypic),
    scalePrior(sigmaSq.scale),
    noscale(noscale),
    snpEffects(data.snpEffectNames)
    {
        sparse = data.sparseLDM;
        // varg.value = varGenotypic; //// NOTE: write it into constructor!!!
        paramSetVec = {&snpEffects, &fixedEffects};
        for (unsigned i=0; i<Pis.size(); ++i) { 
           Pis[i]->value=Pis.values[i];  
        }
        //paramVec     = {&nnzSnp, &sigmaSq, &vare, &varg, &hsq, &kappa};
        paramVec.insert(paramVec.begin(), Pis.begin(), Pis.end());
        //paramToPrint = {&nnzSnp, &sigmaSq, &vare, &varg, &hsq, &kappa, &rounding};
        paramToPrint.insert(paramToPrint.begin(), Pis.begin(), Pis.end());
        if (message && myMPI::rank==0) {
            cout << "\nApproximate Bayes Kappa model fitted." << endl;
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
    
    void sampleUnknowns(void);
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
        
        void sampleFromFC(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const vector<ChromInfo*> &chromInfoVec, const VectorXf &LDsamplVar, const ArrayXf &snp2pqPowS, const VectorXf &snp2pq,
                          const Vector2f &sigmaSq, const Vector3f &pi, const float vare, const float varg,
                          const float ps, const float overdispersion, VectorXf &deltaS);
        
        void sampleFromFC(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag, const VectorXf &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec, const VectorXf &LDsamplVar, const ArrayXf &snp2pqPowS, const VectorXf &snp2pq,
                          const Vector2f &sigmaSq, const Vector3f &pi, const float vare, const float varg,
                          const float ps, const float overdispersion, VectorXf &deltaS);

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
    
    ApproxBayesSMix(const Data &data, const float varGenotypic, const float varResidual, const float pival, const float overdispersion,
                  const bool estimatePS, const float varS, const vector<float> &svalue,
                  const bool message = true):
    ApproxBayesS(data, varGenotypic, varResidual, pival, 1, 1, true, 0, overdispersion, estimatePS, 0, 0, varS, svalue, "HMC", false, false),
    snpEffects(data.snpEffectNames, data.snp2pq, 0.5*pival),
    deltaS(data.snpEffectNames),
    piMixComp(pival),
    sigmaSq(varGenotypic, pival, data.snp2pq)
    {
        paramSetVec = {&snpEffects, &deltaS};
        paramVec = {piMixComp[2], piMixComp[1], &pi, &nnzSnp, sigmaSq[1], &S, sigmaSq[0], hsqMixComp[1], hsqMixComp[0], &hsq};
        paramToPrint = {piMixComp[2], piMixComp[1], &pi, &nnzSnp, sigmaSq[1], &S, sigmaSq[0], hsqMixComp[1], hsqMixComp[0], &hsq, &rounding};
        if (modelPS) {
            paramVec.push_back(&ps);
            paramToPrint.push_back(&ps);
        }
        if (message && myMPI::rank==0) {
            cout << "\nApproximate BayesSMix model fitted." << endl;
        }
    }
    
    void sampleUnknowns(void);

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

    BayesSMix(const Data &data, const float varGenotypic, const float varResidual, const float pival, const float piAlpha, const float piBeta, const bool estimatePi, const float varS, const vector<float> &svalue, const string &algorithm, const bool message = true):
    BayesS(data, varGenotypic, varResidual, pival, piAlpha, piBeta, estimatePi, varS, svalue, "HMC", false),
    snpEffects(data.snpEffectNames, data.snp2pq, 0.5*pival),
    deltaS(data.snpEffectNames),
    piMixComp(pival),
    sigmaSq(varGenotypic, pival, data.snp2pq)
    {
        ghatMixComp.resize(2);
        paramSetVec = {&snpEffects, &deltaS, &fixedEffects};
        paramVec = {piMixComp[2], piMixComp[1], &pi, &nnzSnp, sigmaSq[1], &S, sigmaSq[0], hsqMixComp[1], hsqMixComp[0], &hsq};
        paramToPrint = {piMixComp[2], piMixComp[1], &pi, &nnzSnp, sigmaSq[1], &S, sigmaSq[0], hsqMixComp[1], hsqMixComp[0], &hsq, &rounding};
        if (message && myMPI::rank==0) {
            cout << "\nBayesSMix model fitted." << endl;
        }
    }
    
    void sampleUnknowns(void);
};



#endif /* model_hpp */




