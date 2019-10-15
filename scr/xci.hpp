//
//  xci.hpp
//  gctb
//
//  Created by Jian Zeng on 27/10/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#ifndef xci_hpp
#define xci_hpp

#include <stdio.h>
#include "gctb.hpp"

class XCI {
public:
    unsigned numKeptMales;
    unsigned numKeptFemales;
    
    XCI(){
        numKeptMales   = 0;
        numKeptFemales = 0;
    }
    
    void inputIndInfo(Data &data, const string &bedFile, const string &phenotypeFile, const string &keepIndFile,
                      const unsigned keepIndMax, const unsigned mphen, const string &covariateFile);
    void sortIndBySex(vector<IndInfo*> &indInfoVec);
    void restoreFamFileOrder(vector<IndInfo*> &indInfoVec);
    void inputSnpInfo(Data &data, const string &bedFile, const string &includeSnpFile, const string &excludeSnpFile,
                      const unsigned includeChr, const bool readGenotypes);
    void readBedFile(Data &data, const string &bedFile);

    Model* buildModel(Data &data, const string &bayesType, const float heritability, const float pi, const float piAlpha, const float piBeta, const bool estimatePi, const float piNDC, const bool estimatePiNDC, const float piGxE, const unsigned windowWidth);
    void simu(Data &data, const float pi, const float heritability, const float probNDC, const float probGxS, const bool removeQTL, const string &title, const int seed);
    void outputResults(const Data &data, const vector<McmcSamples*> &mcmcSampleVec, const string &bayesType, const string &title);
};


class BayesCXCI : public BayesC {
public:
    // y = mu + sum_j Z_j beta_j delta_j + e
    // For males,   Z_mj = X_mj
    // For females, Z_fj = X_fj with prob. p (under NDC model) or 0.5*X_fj with prob. 1-p (under FDC model)
    // p ~ U(0,1) or Beta(a,b) is the prob. of NDC model, in other word, the proportion of SNPs that escape from XCI
    // beta_j ~ N(0, sigma^2); sigma^2 ~ scaled-inverse chi-square
    // delta_j ~ Bernoulli(pi)
    
    class FixedEffects : public BayesC::FixedEffects {
    public:
        
        FixedEffects(const vector<string> &header): BayesC::FixedEffects(header){}
        
        void sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &X, const unsigned nmale, const unsigned nfemale,
                          const VectorXf &XPXdiagMale, const VectorXf &XPXdiagFemale, const float varem, const float varef);
    };
    
    class ProbNDC : public Parameter, public Stat::Beta {
    public:
        const float alpha;
        const float beta;
        
        ProbNDC(const float p): Parameter("PiNDC"), alpha(1), beta(1){  // conditional probability on nonzero SNPs, uniform prior
            value = p;
        }
        
        void sampleFromFC(const unsigned numSnps, const unsigned numNDC);
    };
        
    class Gamma : public ParamSet {
    public:
        Gamma(const vector<string> &header): ParamSet("Gamma", header){};
    };
    
    class SnpEffects : public BayesC::SnpEffects {
    public:
        SnpEffects(const vector<string> &header): BayesC::SnpEffects(header, "Gibbs"){};
        
        void sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &ZPZdiagMale,
                          const VectorXf &ZPZdiagFemale, const unsigned nmale, const unsigned nfemale, const float p,
                          const float sigmaSq, const float pi, const float varem, const float varef, VectorXf &gamma, VectorXf &ghatm, VectorXf &ghatf);
    };
    
    class ScaleVar : public BayesC::ScaleVar {
    public:
        const float sum2pq;
        
        ScaleVar(const float sum2pq, const float val): BayesC::ScaleVar(val), sum2pq(sum2pq){}
        
        void compute(const float vg, const float pi, float &scaleVar){
            value = 0.5f*vg/(sum2pq*pi);
            scaleVar = value;
        };
    };
    
    class Rounding : public BayesC::Rounding {
    public:
        Rounding(): BayesC::Rounding(){}
        void computeYcorr(const VectorXf &y, const MatrixXf &X, const MatrixXf &Z,
                          const VectorXf &gamma, const unsigned nmale, const unsigned nfemale,
                          const VectorXf &fixedEffects, const VectorXf &snpEffects,
                          VectorXf &ycorrm, VectorXf &ycorrf);
    };
    
    unsigned nmale, nfemale;
    float genVarPrior;
    float piPrior;
    bool estimatePiNDC;
    
    FixedEffects fixedEffects;
    ProbNDC piGamma;
    Parameter piNDC;
    Gamma gamma;   // indicator variable with 1: NDC, 0: FDC
    SnpEffects snpEffects;
    ScaleVar scale;
    Rounding rounding;

    VectorXf ycorrm;
    VectorXf ycorrf;
    VectorXf ghatm;
    VectorXf ghatf;

    ResidualVar varem;
    ResidualVar varef;
    GenotypicVar vargm;
    GenotypicVar vargf;
    Heritability hsqm;
    Heritability hsqf;
    
    VectorXf XPXdiagMale;
    VectorXf ZPZdiagMale;
    VectorXf XPXdiagFemale;
    VectorXf ZPZdiagFemale;
    
    BayesCXCI(const Data &data, const float varGenotypic, const float varResidual, const float pival, const float piAlpha, const float piBeta, const bool estimatePi, const float piNDCval, const bool estimatePiNDC, const unsigned nmale, const unsigned nfemale, const bool noscale, const bool message = true):
    BayesC(data, varGenotypic, varResidual, pival, piAlpha, piBeta, estimatePi, noscale, "Gibbs", false),
    ycorrm(data.y.head(nmale)),
    ycorrf(data.y.tail(nfemale)),
    fixedEffects(data.fixedEffectNames),
    piGamma(piNDCval),
    piNDC("PiNDC"),
    estimatePiNDC(estimatePiNDC),
    gamma(data.snpEffectNames),
    snpEffects(data.snpEffectNames),
    scale(data.snp2pq.sum(), sigmaSq.scale),
    genVarPrior(varGenotypic),
    piPrior(pival),
    nmale(nmale), nfemale(nfemale),
    varem(varResidual, nmale, "ResVarM"),
    varef(varResidual, nfemale, "ResVarF"),
    vargm(varGenotypic, "GenVarM"),
    vargf(varGenotypic, "GenVarF"),
    hsqm("hsqM"),
    hsqf("hsqF") {
        if (!estimatePiNDC) piGamma.value = 0.5;
        if (!estimatePiNDC && myMPI::rank==0) cout << "\n Fixed piNDC at 0.5!" << endl;
        getZPZdiag(data);
        gamma.values.setZero(data.numIncdSnps);
        paramSetVec = {&snpEffects, &gamma, &fixedEffects};
        paramVec = {&pi, &nnzSnp, &piNDC, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf};
        paramToPrint = {&pi, &nnzSnp, &piNDC, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf, &rounding};
        if (message && myMPI::rank==0)
            cout << "\nBayesCXCI model fitted." << endl;
    }
    
    void sampleUnknowns(void);
    void getZPZdiag(const Data &data);
};


class BayesBXCI : public BayesCXCI {
    // BayesB prior for the SNP effects.
public:
    
    class SnpEffects : public BayesB::SnpEffects {
    public:
        SnpEffects(const vector<string> &header): BayesB::SnpEffects(header){};
        
        void sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &ZPZdiagMale,
                          const VectorXf &ZPZdiagFemale, const unsigned nmale, const unsigned nfemale, const float p,
                          const VectorXf &sigmaSq, const float pi, const float vare, VectorXf &gamma, VectorXf &ghat);
    };

    SnpEffects snpEffects;
    BayesB::VarEffects sigmaSq;

    BayesBXCI(const Data &data, const float varGenotypic, const float varResidual, const float pival, const float piAlpha, const float piBeta, const bool estimatePi, const float piNDCval, const bool estimatePiNDC, const unsigned nmale, const unsigned nfemale, const bool noscale, const bool message = true):
    BayesCXCI(data, varGenotypic, varResidual, pival, piAlpha, piBeta, estimatePi, piNDCval, estimatePiNDC, nmale, nfemale, false),
    snpEffects(data.snpEffectNames),
    sigmaSq(varGenotypic, data.snp2pq, pival)
    {
        paramSetVec = {&snpEffects, &gamma, &fixedEffects};
        paramVec = {&pi, &nnzSnp, &piNDC, &scale, &vare, &varg, &hsq};
        paramToPrint = {&pi, &nnzSnp, &piNDC, &scale, &vare, &varg, &hsq, &rounding};
        if (message && myMPI::rank==0)
            cout << "\nBayesBXCI model fitted." << endl;
    }
    
    void sampleUnknowns(void);

};


class BayesCXCIgxs : public BayesCXCI {
public:
    // Allow for genotype-by-sex effect for each SNP. That is, the SNP effects in males and females are allowed to be different
    // y = mu + sum_j Z_j beta_j delta_j + e
    // For males,   Z_mj = X_mj
    // For females, Z_fj = X_fj with prob. p (under NDC model) or 0.5*X_fj with prob. 1-p (under FDC model)
    // p ~ U(0,1) or Beta(a,b) is the prob. of NDC model, in other word, the proportion of SNPs that escape from XCI
    // beta_j ~ N(0, sigma^2)*pi1 + N(c(0,0), I sigma^2)*pi2 + 0*(1-pi1-pi2); sigma^2 ~ scaled-inverse chi-square
    // pi ~ Dirichlet(1)
    
    
    class SnpEffects : public BayesCXCI::SnpEffects {
    public:
        MatrixXf values;   // 1st column: male effects; 2nd column: female effects
        Vector3f numSnpMixComp;
        
        SnpEffects(const vector<string> &header): BayesCXCI::SnpEffects(header){
            values.setZero(header.size(), 2);
        };
        
        void sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &ZPZdiagMale,
                          const VectorXf &ZPZdiagFemale, const unsigned nmale, const unsigned nfemale, const float piNDC,
                          const float sigmaSq, const Vector3f &pis, const float varem, const float varef,
                          VectorXf &gamma, VectorXf &deltaGxS, VectorXf &ghatm, VectorXf &ghatf);
    };
    
    class ProbMixComps : public BayesR::ProbMixComps {
    public:
        ProbMixComps(const VectorXf &pis): BayesR::ProbMixComps(pis){}
    };
    
    class ProbGxS : public Parameter {
    public:
        ProbGxS(const float p): Parameter("PiGxS"){value = p;}
    };
    
    class SnpEffectsMale : public ParamSet {
    public:
        SnpEffectsMale(const vector<string> &header): ParamSet("SnpEffectsMale", header){}
    };
    
    class SnpEffectsFemale : public ParamSet {
    public:
        SnpEffectsFemale(const vector<string> &header): ParamSet("SnpEffectsFemale", header){}
    };
    
    class DeltaGxS : public ParamSet {
    public:
        DeltaGxS(const vector<string> &header): ParamSet("DeltaGxS", header){};
    };
    
    class Rounding : public BayesCXCI::Rounding {
    public:
        Rounding(): BayesCXCI::Rounding(){}
        void computeYcorr(const VectorXf &y, const MatrixXf &X, const MatrixXf &Z,
                          const VectorXf &gamma, const unsigned nmale, const unsigned nfemale,
                          const VectorXf &fixedEffects, const MatrixXf &snpEffects,
                          VectorXf &ycorrm, VectorXf &ycorrf);
    };
    
    SnpEffects snpEffects;
    SnpEffectsMale snpEffectsMale;
    SnpEffectsFemale snpEffectsFemale;
    ProbMixComps pis;
    ProbGxS piGxS;
    DeltaGxS deltaGxS;
    Rounding rounding;
    
    BayesCXCIgxs(const Data &data, const float varGenotypic, const float varResidual, const Vector3f &pival, const bool estimatePi, const float piNDCval, const bool estimatePiNDC, const unsigned nmale, const unsigned nfemale, const bool noscale, const bool message = true):
    BayesCXCI(data, varGenotypic, varResidual, pival[1]+pival[2], 1, 1, estimatePi, piNDCval, estimatePiNDC, nmale, nfemale, noscale, false),
    snpEffects(data.snpEffectNames),
    snpEffectsMale(data.snpEffectNames),
    snpEffectsFemale(data.snpEffectNames),
    pis(pival),
    piGxS(pival[2]),
    deltaGxS(data.snpEffectNames) {
        paramSetVec = {&snpEffectsMale, &snpEffectsFemale, &gamma, &deltaGxS, &fixedEffects};
        paramVec = {&pi, &nnzSnp, &piNDC, &piGxS, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf};
        paramToPrint = {&pi, &nnzSnp, &piNDC, &piGxS, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf, &rounding};
        if (message && myMPI::rank==0) {
            cout << "\nBayesCXCIgxs model fitted." << endl;
            cout << "sigmaSq: " << sigmaSq.value << endl;
        }
    }
    
    void sampleUnknowns(void);
};

class SBayesCXCIgxs : public BayesCXCIgxs {
public:
    // Same as BayesCXCIgxs but make use of LD matrix for RHS updating scheme
    class SnpEffects : public BayesCXCI::SnpEffects {
    public:
        MatrixXf values;   // 1st column: male effects; 2nd column: female effects
        Vector3f numSnpMixComp;
        
        SnpEffects(const vector<string> &header): BayesCXCI::SnpEffects(header){
            values.setZero(header.size(), 2);
        };
        
        void sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &ZPZdiagMale,
                          const VectorXf &ZPZdiagFemale, const unsigned nmale, const unsigned nfemale, const float piNDC,
                          const float sigmaSq, const Vector3f &pis, const float varem, const float varef,
                          VectorXf &gamma, VectorXf &deltaGxS, VectorXf &ghatm, VectorXf &ghatf);
    };

};

class BayesNXCIgxs : public BayesCXCIgxs {
public:
    
    class SnpEffects : public BayesCXCIgxs::SnpEffects {
    public:
        unsigned numWindows;
        unsigned numNonZeroWind;
        
        const VectorXi &windStart;
        const VectorXi &windSize;
        
        VectorXf windDeltaEff;
        VectorXf windDeltaNDC;
        VectorXf windDeltaGxS;
        
        
        SnpEffects(const vector<string> &header, const VectorXi &windStart, const VectorXi &windSize):
        BayesCXCIgxs::SnpEffects(header), windStart(windStart), windSize(windSize){
            numWindows = (unsigned) windStart.size();
            windDeltaEff.setZero(numWindows);
            windDeltaNDC.setZero(numWindows);
            windDeltaGxS.setZero(numWindows);
        };
        
        void sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z, const vector<MatrixXf> &ZPZdiagMale,
                          const vector<MatrixXf> &ZPZdiagFemale, const unsigned nmale, const unsigned nfemale, const float piNDC,
                          const float sigmaSq, const Vector3f &pis, const float varem, const float varef,
                          VectorXf &gamma, VectorXf &deltaGxS, VectorXf &ghatm, VectorXf &ghatf);
    };
    
    SnpEffects snpEffects;
    BayesN::NumNonZeroWind nnzWind;
    
    BayesNXCIgxs(const Data &data, const float varGenotypic, const float varResidual, const Vector3f &pival, const bool estimatePi, const float piNDCval, const bool estimatePiNDC, const unsigned nmale, const unsigned nfemale, const bool noscale, const bool message = true):
    BayesCXCIgxs(data, varGenotypic, varResidual, pival, estimatePi, piNDCval, estimatePiNDC, nmale, nfemale, noscale, false),
    snpEffects(data.snpEffectNames, data.windStart, data.windSize)
    {
        getZPZblockDiag(data);
        paramSetVec = {&snpEffectsMale, &snpEffectsFemale, &gamma, &deltaGxS, &fixedEffects};
        paramVec = {&pi, &nnzWind, &nnzSnp, &piNDC, &piGxS, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf};
        paramToPrint = {&pi, &nnzWind, &nnzSnp, &piNDC, &piGxS, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf, &rounding};
        if (message && myMPI::rank==0)
            cout << "\nBayesNXCIgxs model fitted." << endl;
    }
    
    vector<MatrixXf> ZPZblockDiagMale;
    vector<MatrixXf> ZPZblockDiagFemale;

    void sampleUnknowns(void);
    void getZPZblockDiag(const Data &data);
};


class BayesXgxs : public BayesCXCIgxs {
public:
    
    class SnpEffects : public BayesCXCIgxs::SnpEffects {
    public:
        Vector3f numSnpBetaMix;
        Vector3f numSnpDosageMix;
        
        SnpEffects(const vector<string> &header): BayesCXCIgxs::SnpEffects(header){};
        
        void sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z, const VectorXf &ZPZdiagMale,
                          const VectorXf &ZPZdiagFemale, const unsigned nmale, const unsigned nfemale,
                          const float sigmaSq, const Vector3f &piDosage, const Vector3f &piBeta, const float varem, const float varef,
                          VectorXf &gamma, VectorXf &deltaNDC, VectorXf &deltaFDC, VectorXf &deltaGxS, VectorXf &ghatm, VectorXf &ghatf);
    };

    class Prob : public Parameter {
    public:
        Prob(const string &label, const float p): Parameter(label){value = p;}
    };
    
    class Delta : public ParamSet {
    public:
        Delta(const string &label, const vector<string> &header): ParamSet(label, header){};
    };

    SnpEffects snpEffects;
    ProbMixComps piDosage;
    ProbMixComps piBeta;
    Prob piNDC;
    Prob piFDC;
    Prob piGxS;
    Delta deltaNDC;
    Delta deltaFDC;
    Delta deltaGxS;
    
    BayesXgxs(const Data &data, const float varGenotypic, const float varResidual, const Vector3f &piBetaVal, const Vector3f &piDosageVal, const bool estimatePi, const unsigned nmale, const unsigned nfemale, const bool noscale, const bool message = true):
    BayesCXCIgxs(data, varGenotypic, varResidual, piBetaVal, estimatePi, piDosageVal[1], nmale, nfemale, noscale, false),
    snpEffects(data.snpEffectNames),
    piDosage(piDosageVal),
    piBeta(piBetaVal),
    piNDC("PiNDC",piDosageVal[1]),
    piFDC("PiFDC",piDosageVal[2]),
    piGxS("PiGxS",piBetaVal[2]),
    deltaNDC("DeltaNDC", data.snpEffectNames),
    deltaFDC("DeltaFDC", data.snpEffectNames),
    deltaGxS("DeltaGxS", data.snpEffectNames) {
        paramSetVec = {&snpEffectsMale, &snpEffectsFemale, &deltaNDC, &deltaFDC, &deltaGxS, &fixedEffects};
        paramVec = {&pi, &nnzSnp, &piNDC, &piFDC, &piGxS, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf};
        paramToPrint = {&pi, &nnzSnp, &piNDC, &piFDC, &piGxS, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf, &rounding};
        if (message && myMPI::rank==0) {
            cout << "\nBayesXgxs model fitted." << endl;
            cout << "sigmaSq: " << sigmaSq.value << endl;
        }
    }
    
    void sampleUnknowns(void);
};

#endif /* xci_hpp */
