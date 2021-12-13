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
                      const unsigned keepIndMax, const unsigned mphen, const string &covariateFile, const bool femaleOnly = false);
    void sortIndBySex(vector<IndInfo*> &indInfoVec);
    void restoreFamFileOrder(vector<IndInfo*> &indInfoVec);
    void inputSnpInfo(Data &data, const string &bedFile, const string &includeSnpFile, const string &excludeSnpFile,
                      const unsigned includeChr, const string &annotationFile, const string &windowFile, const bool readGenotypes);
    void readBedFile(Data &data, const string &bedFile);

    Model* buildModel(Data &data, const string &bayesType, const float heritability, const float pi, const VectorXf &piPar, const bool estimatePi, const float piNDC, const Vector2f &piNDCpar, const bool estimatePiNDC, const float piGxE, const bool estimatePiGxE, const unsigned windowWidth);
    
    Model* buildModelStageOne(Data &data, const string &bayesType, const float heritability, const float pi, const VectorXf &piPar, const bool estimatePi, const float piNDC, const Vector2f &piNDCpar, const bool estimatePiNDC);
    Model* buildModelStageTwo(Data &data, const string &bayesType, const float heritability, const float pi, const VectorXf &piPar, const bool estimatePi, const float piNDC, const Vector2f &piNDCpar, const bool estimatePiNDC, const string &snpResFile, const float piGxE, const bool estimatePiGxE);

    vector<McmcSamples*> multi_chain_mcmc(Data &data, const string &bayesType, const float heritability, const float pi, const VectorXf &piPar, const bool estimatePi, const float piNDC, const Vector2f &piNDCpar, const bool estimatePiNDC, const float piGxE, const bool estimatePiGxE, const unsigned numChains, const unsigned chainLength, const unsigned burnin, const unsigned thin, const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior);

    void simu(Data &data, const float pi, const float heritability, const float probNDC, const float probGxS, const bool removeQTL, const string &title, const int seed);
    void outputResults(const Data &data, const vector<McmcSamples*> &mcmcSampleVec, const string &bayesType, const string &title);
    void readSnpPiNDC(VectorXf &snpPiNDC, const Data &data, const string &snpResFile);
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
        
        ProbNDC(const float p, const Vector2f &par): Parameter("PiNDC"), alpha(par[0]), beta(par[1]){  // conditional probability on nonzero SNPs, uniform prior
            value = p;
        }
        
        void sampleFromFC(const unsigned numSnps, const unsigned numNDC);
        void sampleFromPrior(void);
    };
        
    class DeltaNDC : public ParamSet {
    public:
        DeltaNDC(const vector<string> &header): ParamSet("DeltaNDC", header){};
    };
    
    class SnpEffects : public BayesC::SnpEffects {
    public:
        SnpEffects(const vector<string> &header): BayesC::SnpEffects(header, "Gibbs"){};
        
        void sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &ZPZdiagMale, const VectorXf &ZPZdiagFemale,
                          const VectorXf &ZPZdiagMaleRank, const VectorXf &ZPZdiagFemaleRank, const unsigned nmale, const unsigned nfemale, const float p,
                          const float sigmaSq, const float pi, const float varem, const float varef, VectorXf &deltaNDC, VectorXf &ghatm, VectorXf &ghatf);
    };
    
    class VarEffects : public BayesC::VarEffects {
    public:
        // a priori assuming full dosage compensation model and the genetic variance computed from females
        // so the genetic variance is 0.25*m*pi*2pq*sigmaSq
        
        VarEffects(const float vg, const VectorXf &snp2pq, const float pi, const float piNDC, const bool noscale, const string &lab = "SigmaSq"):
        BayesC::VarEffects(vg, snp2pq, pi, noscale, lab) {
            if (noscale == true) {
                //value = vg / ((0.25*(1.0-piNDC) + piNDC) * snp2pq.sum() * pi);  // derived from prior knowledge on Vg and pi
                value = vg / (0.25 * snp2pq.sum() * pi);  // derived from prior knowledge on Vg and pi
            } else {
                //value = vg / ((0.25*(1.0-piNDC) + piNDC) * snp2pq.size() * pi);  // derived from prior knowledge on Vg and pi
                value = vg / (0.25 * snp2pq.size() * pi);  // derived from prior knowledge on Vg and pi
            }
            
            scale = 0.5f*value;  // due to df = 4

            cout << "sigmaSq " << value << " " << vg << " " << snp2pq.sum() << " " << pi << " " << noscale << endl;
        }
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
                          const VectorXf &deltaNDC, const unsigned nmale, const unsigned nfemale,
                          const VectorXf &fixedEffects, const VectorXf &snpEffects,
                          VectorXf &ycorrm, VectorXf &ycorrf);
    };
    
    unsigned nmale, nfemale;
    float genVarPrior;
    float piPrior;
    bool estimatePiNDC;
    
    FixedEffects fixedEffects;
    ProbNDC piDeltaNDC;
    Parameter piNDC;
    DeltaNDC deltaNDC;   // indicator variable with 1: NDC, 0: FDC
    SnpEffects snpEffects;
    VarEffects sigmaSq;
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
    VectorXf ZPZdiagMaleRank;   // for MPI
    VectorXf XPXdiagFemale;
    VectorXf ZPZdiagFemale;
    VectorXf ZPZdiagFemaleRank; // for MPI
    
    BayesCXCI(const Data &data, const float varGenotypic, const float varResidual, const float pival, const float piAlpha, const float piBeta, const bool estimatePi, const float piNDCval, const Vector2f &piNDCpar, const bool estimatePiNDC, const unsigned nmale, const unsigned nfemale, const bool noscale, const bool message = true, const bool randomStart = false):
    BayesC(data, varGenotypic, varResidual, 0.0, pival, piAlpha, piBeta, estimatePi, noscale, "Gibbs", false),
    ycorrm(data.y.head(nmale)),
    ycorrf(data.y.tail(nfemale)),
    fixedEffects(data.fixedEffectNames),
    piDeltaNDC(piNDCval, piNDCpar),
    piNDC("PiNDC"),
    estimatePiNDC(estimatePiNDC),
    deltaNDC(data.snpEffectNames),
    snpEffects(data.snpEffectNames),
    sigmaSq(varGenotypic, data.snp2pq, pival, piNDCval, noscale),
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
        getZPZdiag(data);
        deltaNDC.values.setZero(data.numIncdSnps);
        paramSetVec = {&snpEffects, &deltaNDC, &fixedEffects};
        paramVec = {&pi, &nnzSnp, &piNDC, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf};
        paramToPrint = {&pi, &nnzSnp, &piNDC, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf, &rounding};
        if (message)
            cout << "\nBayesCXCI model fitted." << endl;
        if (randomStart) sampleStartVal();
    }
    
    void sampleUnknowns(void);
    void sampleStartVal(void);
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
                          const VectorXf &sigmaSq, const float pi, const float vare, VectorXf &deltaNDC, VectorXf &ghat);
    };

    SnpEffects snpEffects;
    BayesB::VarEffects sigmaSq;

    BayesBXCI(const Data &data, const float varGenotypic, const float varResidual, const float pival, const float piAlpha, const float piBeta, const bool estimatePi, const float piNDCval, const Vector2f &piNDCpar, const bool estimatePiNDC, const unsigned nmale, const unsigned nfemale, const bool noscale, const bool message = true):
    BayesCXCI(data, varGenotypic, varResidual, pival, piAlpha, piBeta, estimatePi, piNDCval, piNDCpar, estimatePiNDC, nmale, nfemale, false),
    snpEffects(data.snpEffectNames),
    sigmaSq(varGenotypic, data.snp2pq, pival, false)
    {
        paramSetVec = {&snpEffects, &deltaNDC, &fixedEffects};
        paramVec = {&pi, &nnzSnp, &piNDC, &scale, &vare, &varg, &hsq};
        paramToPrint = {&pi, &nnzSnp, &piNDC, &scale, &vare, &varg, &hsq, &rounding};
        if (message)
            cout << "\nBayesBXCI model fitted." << endl;
    }
    
    void sampleUnknowns(void);

};


class SBayesCXCI : public BayesCXCI {
    // efficient BayesCXCI using right-hand-side updating strategy
public:
    class SnpEffects : public BayesCXCI::SnpEffects {
    public:
        SnpEffects(const vector<string> &header): BayesCXCI::SnpEffects(header){};
        
        void sampleFromFC(VectorXf &rcorrm, VectorXf &rcorrf, const MatrixXf &ZPZ, const MatrixXf &ZPZmale, const MatrixXf &ZPZfemale, const float piNDC,
                          const float sigmaSq, const float pi, const float varem, const float varef, VectorXf &deltaNDC, VectorXf &ghatm, VectorXf &ghatf);
    };
    
    void sampleUnknowns(void);

    VectorXf rcorrm;
    VectorXf rcorrf;
    
    SnpEffects snpEffects;

    SBayesCXCI(const Data &data, const float varGenotypic, const float varResidual, const float pival, const VectorXf &piPar, const bool estimatePi, const float piNDCval, const Vector2f &piNDCpar, const bool estimatePiNDC, const unsigned nmale, const unsigned nfemale, const bool noscale, const bool message = true):
    BayesCXCI(data, varGenotypic, varResidual, pival, piPar[0], piPar[1], estimatePi, piNDCval, piNDCpar, estimatePiNDC, nmale, nfemale, noscale, false),
    snpEffects(data.snpEffectNames) {
        paramSetVec = {&snpEffects, &deltaNDC, &fixedEffects};
        paramVec = {&pi, &nnzSnp, &piNDC, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf};
        paramToPrint = {&pi, &nnzSnp, &piNDC, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf, &rounding};
        if (message)
            cout << "\nSBayesCXCI model fitted." << endl;
    }
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
        VectorXf delta;
        float numNDCandNoGXS;

        SnpEffects(const vector<string> &header): BayesCXCI::SnpEffects(header){
            values.setZero(header.size(), 2);
            delta.setZero(size);
        };
        
        void sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &ZPZdiagMale, const VectorXf &ZPZdiagFemale,
                          const VectorXf &ZPZdiagMaleRank, const VectorXf &ZPZdiagFemaleRank, const unsigned nmale, const unsigned nfemale, const float piNDC, VectorXf &snpPiNDC,
                          const float sigmaSq, const Vector3f &pis, const float varem, const float varef,
                          VectorXf &deltaNDC, VectorXf &deltaGxS, VectorXf &ghatm, VectorXf &ghatf);
        
        void sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z, const VectorXf &ZPZdiag, const VectorXf &ZPZdiagMale, const VectorXf &ZPZdiagFemale,
                          const VectorXf &ZPZdiagMaleRank, const VectorXf &ZPZdiagFemaleRank, const unsigned nmale, const unsigned nfemale,
                          VectorXf &snpPiNDC, const VectorXf &logSnpPiNDC, const VectorXf &logSnpPiNDCcomp,
                          const float sigmaSq, const Vector3f &pis, const float varem, const float varef,
                          VectorXf &deltaNDC, VectorXf &deltaGxS, VectorXf &ghatm, VectorXf &ghatf);
    };
    
    class ProbMixComps : public BayesR::ProbMixComps {
    public:
        ProbMixComps(const VectorXf &pis, const VectorXf &piPar): BayesR::ProbMixComps(pis, piPar){
            cout << "alphaVec " << alphaVec.transpose() << endl;
        }
        
        void getValues(VectorXf &pis);
        void sampleFromPrior(void);
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
                          const VectorXf &deltaNDC, const unsigned nmale, const unsigned nfemale,
                          const VectorXf &fixedEffects, const MatrixXf &snpEffects,
                          VectorXf &ycorrm, VectorXf &ycorrf);
    };
    
    class AnnoEffects : public vector<Parameter*>, public Stat::Normal {
    public:
        unsigned numAnnos;
        VectorXf values;
        float varProp;  // proposal variance
        float sigmaSq;  // prior variance

        BayesS::AcceptanceRate ar;

        AnnoEffects(const vector<string> &annoNames, const float piNDC){
            if (annoNames.size()) {
                numAnnos = annoNames.size();
                for (unsigned i = 0; i<numAnnos; ++i) {
                    this->push_back(new Parameter(annoNames[i]+"Effect"));
                }
                values.setZero(numAnnos);
                values[0] = log(piNDC/(1.0f-piNDC));  // initial value of intercept
                varProp = 0.1;
                sigmaSq = 5;
            }
        }
        
        void sampleFromFC(VectorXf &snpPiNDC, VectorXf &logSnpPiNDC, VectorXf &logSnpPiNDCcomp, const MatrixXf &annoMat, const MatrixXf &APA, VectorXf &deltaNDC, VectorXf &deltaBeta);
    };
        
    class DeltaNDC2 : public ParamSet {
    public:
        DeltaNDC2(const vector<string> &header): ParamSet("DeltaNDC2", header){};
        
        void compute(const VectorXf &snpPiNDC, const float &piNDC);
    };
    
    class DeltaWindow : public ParamSet {
    public:
        DeltaWindow(const string &label, const vector<string> &header): ParamSet(label, header){};
    };

    SnpEffects snpEffects;
    SnpEffectsMale snpEffectsMale;
    SnpEffectsFemale snpEffectsFemale;
    ProbMixComps pis;
    ProbGxS piGxS;
    DeltaGxS deltaGxS;
    AnnoEffects annoEffects;
    Rounding rounding;
    
    DeltaNDC2 deltaNDC2;
    DeltaWindow deltaPIPwind;  // at least one of the SNPs in the window has an effect
    DeltaWindow deltaNDCwind;  // at least one of the SNPs in the window has an effect and is under NDC
    DeltaWindow deltaGxSwind;  // at least one of the SNPs in the window has a GxS effect

    bool estimatePiGxS;
    float piGxSgiven;
    float snpPiNDCgiven;
    
    VectorXf snpPiNDC;
    VectorXf logSnpPiNDC;
    VectorXf logSnpPiNDCcomp;

    BayesCXCIgxs(const Data &data, const float varGenotypic, const float varResidual, const Vector3f &pival, const VectorXf &piPar, const bool estimatePi, const float piNDCval, const Vector2f &piNDCpar, const bool estimatePiNDC, const bool estimatePiGxE, const unsigned nmale, const unsigned nfemale, const bool noscale, const bool message = true, const bool randomStart = false):
    BayesCXCI(data, varGenotypic, varResidual, pival[1]+pival[2], 1, 1, estimatePi, piNDCval, piNDCpar, estimatePiNDC, nmale, nfemale, noscale, false),
    snpEffects(data.snpEffectNames),
    snpEffectsMale(data.snpEffectNames),
    snpEffectsFemale(data.snpEffectNames),
    pis(pival, piPar.head(3)),
    piGxS(pival[2]),
    estimatePiGxS(estimatePiGxE),
    piGxSgiven(pival[2]/(1.0-pival[0])),
    deltaGxS(data.snpEffectNames),
    annoEffects(data.annoNames, piNDCval),
    deltaNDC2(data.snpEffectNames),
    deltaPIPwind("DeltaPIPwind", data.snpEffectNames),
    deltaNDCwind("DeltaNDCwind", data.snpEffectNames),
    deltaGxSwind("DeltaGxSwind", data.snpEffectNames){
        snpPiNDCgiven = false;
        paramSetVec = {&snpEffectsMale, &snpEffectsFemale, &deltaNDC, &deltaGxS, &fixedEffects, &deltaNDC2};
        paramVec = {&pi, &nnzSnp, &piNDC, &piGxS, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf};
        paramToPrint = {&pi, &nnzSnp, &piNDC, &piGxS, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf, &rounding};
        snpPiNDC.setConstant(data.numIncdSnps, piNDCval);
        if (data.numAnnos) {
            logSnpPiNDC = snpPiNDC.array().log();
            logSnpPiNDCcomp = (1.0f - snpPiNDC.array()).log();
            paramVec.insert(paramVec.end(), annoEffects.begin(), annoEffects.end());
            paramToPrint.insert(paramToPrint.end(), annoEffects.begin(), annoEffects.end());
        }
        if (data.numWindows) {
            paramSetVec.push_back(&deltaPIPwind);
            paramSetVec.push_back(&deltaNDCwind);
            paramSetVec.push_back(&deltaGxSwind);
        }
        if (message) {
            cout << "\nBayesCXCIgxs model fitted." << endl;
            cout << "sigmaSq: " << sigmaSq.value << endl;
        }
        
        if (randomStart) sampleStartVal();

    }
    
    BayesCXCIgxs(const Data &data, const float varGenotypic, const float varResidual, const Vector3f &pival, const VectorXf &piPar, const bool estimatePi, const float piNDCval, const Vector2f &piNDCpar, const bool estimatePiNDC, const VectorXf &snpPiNDC, const bool estimatePiGxE, const unsigned nmale, const unsigned nfemale, const bool noscale, const bool message = true):
    BayesCXCI(data, varGenotypic, varResidual, pival[1]+pival[2], 1, 1, estimatePi, piNDCval, piNDCpar, estimatePiNDC, nmale, nfemale, noscale, false),
    snpEffects(data.snpEffectNames),
    snpEffectsMale(data.snpEffectNames),
    snpEffectsFemale(data.snpEffectNames),
    pis(pival, piPar.head(3)),
    piGxS(pival[2]),
    estimatePiGxS(estimatePiGxE),
    piGxSgiven(pival[2]/(1.0-pival[0])),
    deltaGxS(data.snpEffectNames),
    snpPiNDC(snpPiNDC),
    annoEffects(data.annoNames, piNDCval),
    deltaNDC2(data.snpEffectNames),
    deltaPIPwind("DeltaPIPwind", data.snpEffectNames),
    deltaNDCwind("DeltaNDCwind", data.snpEffectNames),
    deltaGxSwind("DeltaGxSwind", data.snpEffectNames){
        snpPiNDCgiven = true;
        logSnpPiNDC = snpPiNDC.array().log();
        logSnpPiNDCcomp = (1.0f - snpPiNDC.array()).log();
        paramSetVec = {&snpEffectsMale, &snpEffectsFemale, &deltaNDC, &deltaGxS, &fixedEffects, &deltaNDC2};
        paramVec = {&pi, &nnzSnp, &piNDC, &piGxS, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf};
        paramToPrint = {&pi, &nnzSnp, &piNDC, &piGxS, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf, &rounding};
        if (data.numWindows) {
            paramSetVec.push_back(&deltaPIPwind);
            paramSetVec.push_back(&deltaNDCwind);
            paramSetVec.push_back(&deltaGxSwind);
        }
        if (message) {
            cout << "\nBayesCXCIgxs model fitted." << endl;
            cout << "sigmaSq: " << sigmaSq.value << endl;
        }
    }

    
    void sampleUnknowns(void);
    void sampleStartVal(void);
    void computeWindowDelta(const Data &data,
                            const VectorXf &deltaPIPsnp, const VectorXf &deltaNDCsnp, const VectorXf &deltaGxSsnp,
                            VectorXf &deltaPIPwind, VectorXf &deltaNDCwind, VectorXf &deltaGxSwind);
};

class BayesCXCIgxs2 : public BayesCXCI {
public:
    // y = mu + sum_j Z_j beta_j + e
    // For males,   Z_mj = X_mj
    // For females, Z_fj = X_fj with prob. piNDC (under NDC model) or 0.5*X_fj with prob. 1-piNDC (under FDC model)
    // beta_j ~ mixture*piSNP + 0*(1-piSNP); mixture = N(0, I sigma^2)*piGXS + N(0, 1 sigma^2)*(1-piGXS)
    // piSNP, piGXS, piNDC ~ i.i.d. Beta(0.1, 1)
    // sigma^2 ~ scaled-inverse chi-square
    
    class SnpEffects : public BayesCXCI::SnpEffects {
    public:
        MatrixXf values;   // 1st column: male effects; 2nd column: female effects
        unsigned numGXS;
        unsigned numNDC;
        VectorXf deltaSNP;
        
        SnpEffects(const vector<string> &header): BayesCXCI::SnpEffects(header){
            values.setZero(header.size(), 2);
        };

        void sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z,
                          const VectorXf &ZPZdiagMaleRank, const VectorXf &ZPZdiagFemaleRank,
                          const VectorXf &ZPZdiagMale, const VectorXf &ZPZdiagFemale,
                          const unsigned nmale, const unsigned nfemale,
                          const float piSNP, const float piGXS,
                          VectorXf &snpPiNDC, const VectorXf &logSnpPiNDC, const VectorXf &logSnpPiNDCcomp,
                          const float sigmaSq, const float varem, const float varef,
                          VectorXf &deltaNDC, VectorXf &deltaGxS, VectorXf &ghatm, VectorXf &ghatf);
    };
    
    SnpEffects snpEffects;
    BayesCXCIgxs::SnpEffectsMale snpEffectsMale;
    BayesCXCIgxs::SnpEffectsFemale snpEffectsFemale;
    BayesCXCIgxs::DeltaGxS deltaGXS;
    BayesCXCIgxs::Rounding rounding;
    BayesCXCIgxs::AnnoEffects annoEffects;

    BayesC::Pi piSNP;
    BayesC::Pi piGXS;
    BayesC::Pi piNDC;
    
    bool estimatePiSNP;
    bool estimatePiNDC;
    bool estimatePiGXS;
    
    VectorXf snpPiNDC;
    VectorXf logSnpPiNDC;
    VectorXf logSnpPiNDCcomp;

    BayesCXCIgxs2(const Data &data, const float varGenotypic, const float varResidual, const float piSNPval, const bool estimatePiSNP, const float piGXSval, const bool estimatePiGXS, const float piNDCval, const bool estimatePiNDC, const Vector2f &piNDCpar, const unsigned nmale, const unsigned nfemale, const bool noscale, const bool message = true, const bool randomStart = false):
    BayesCXCI(data, varGenotypic, varResidual, piSNPval, 1, 1, estimatePiSNP, piNDCval, piNDCpar, estimatePiNDC, nmale, nfemale, noscale, false),
    snpEffects(data.snpEffectNames),
    snpEffectsMale(data.snpEffectNames),
    snpEffectsFemale(data.snpEffectNames),
    piSNP(piSNPval, 0.1, 1, "PiSnp"),
    piGXS(piGXSval, 0.1, 1, "PiGxS"),
    piNDC(piNDCval, 0.1, 1, "PiNDC"),
    estimatePiSNP(estimatePiSNP),
    estimatePiGXS(estimatePiGXS),
    estimatePiNDC(estimatePiNDC),
    deltaGXS(data.snpEffectNames),
    annoEffects(data.annoNames, piNDCval){
        snpPiNDC.setConstant(data.numIncdSnps, piNDCval);
        logSnpPiNDC = snpPiNDC.array().log();
        logSnpPiNDCcomp = (1.0f - snpPiNDC.array()).log();
        paramSetVec = {&snpEffectsMale, &snpEffectsFemale, &deltaNDC, &deltaGXS, &fixedEffects};
        paramVec = {&nnzSnp, &piSNP, &piNDC, &piGXS, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf};
        paramToPrint = {&nnzSnp, &piSNP, &piNDC, &piGXS, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf, &rounding};
        if (message) {
            cout << "\nBayesCXCIgxs2 model fitted." << endl;
            cout << "sigmaSq: " << sigmaSq.value << endl;
        }
        
        if (randomStart) sampleStartVal();
    }
    
    void sampleUnknowns(void);
    void sampleStartVal(void);
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
                          VectorXf &deltaNDC, VectorXf &deltaGxS, VectorXf &ghatm, VectorXf &ghatf);
    };

};

class BayesNXCIgxs : public BayesCXCIgxs {
public:
    
    class SnpEffects : public BayesCXCIgxs::SnpEffects {
    public:
        unsigned numWindows;
        unsigned numNonZeroWind;
        float numWindNDCandNoGXS;
        
        const VectorXi &windStart;
        const VectorXi &windSize;
        
        Vector3f numWindMixComp;
        VectorXf windDeltaNDC;
        VectorXf windDeltaGxS;
        VectorXf windDelta;
                
        SnpEffects(const vector<string> &header, const VectorXi &windStart, const VectorXi &windSize):
        BayesCXCIgxs::SnpEffects(header), windStart(windStart), windSize(windSize){
            numWindows = (unsigned) windSize.size();
            windDeltaNDC.setZero(numWindows);
            windDeltaGxS.setZero(numWindows);
            windDelta.setZero(numWindows);
        };
        
        void sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z,
                          const vector<MatrixXf> &ZPZdiagMale, const vector<MatrixXf> &ZPZdiagFemale,
                          const vector<vector<unsigned> > &windowSnpIdxVec,
                          const unsigned nmale, const unsigned nfemale, const float piNDC,
                          const float sigmaSq, const Vector3f &pis, const float varem, const float varef,
                          VectorXf &deltaNDC, VectorXf &deltaGxS, VectorXf &ghatm, VectorXf &ghatf);

//        void sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z,
//                          const vector<MatrixXf> &ZPZdiagMale, const vector<MatrixXf> &ZPZdiagFemale,
//                          const unsigned nmale, const unsigned nfemale, const float piNDC,
//                          const float sigmaSq, const Vector3f &pis, const float varem, const float varef,
//                          VectorXf &deltaNDC, VectorXf &deltaGxS, VectorXf &ghatm, VectorXf &ghatf);
    };
    
    SnpEffects snpEffects;
    BayesN::NumNonZeroWind nnzWind;
    
    BayesNXCIgxs(const Data &data, const float varGenotypic, const float varResidual, const Vector3f &pival, const VectorXf &piPar, const bool estimatePi, const float piNDCval, const Vector2f &piNDCpar, const bool estimatePiNDC, const bool estimatePiGxE, const unsigned nmale, const unsigned nfemale, const bool noscale, const bool message = true):
    BayesCXCIgxs(data, varGenotypic, varResidual, pival, piPar, estimatePi, piNDCval, piNDCpar, estimatePiNDC, estimatePiGxE, nmale, nfemale, noscale, false),
    snpEffects(data.snpEffectNames, data.windStart, data.windSize)
    {
        getZPZblockDiag(data);
        paramSetVec = {&snpEffectsMale, &snpEffectsFemale, &deltaNDC, &deltaGxS, &fixedEffects, &deltaNDC2};
        paramVec = {&pi, &nnzWind, &nnzSnp, &piNDC, &piGxS, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf};
        paramToPrint = {&pi, &nnzWind, &nnzSnp, &piNDC, &piGxS, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf, &rounding};
        if (message)
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
                          VectorXf &deltaNDC, VectorXf &deltaFDC, VectorXf &deltaGxS, VectorXf &ghatm, VectorXf &ghatf);
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
    
    BayesXgxs(const Data &data, const float varGenotypic, const float varResidual, const Vector3f &piBetaVal, const VectorXf &piPar, const Vector3f &piDosageVal, const Vector2f &piNDCpar, const bool estimatePi, const bool estimatePiGxE, const unsigned nmale, const unsigned nfemale, const bool noscale, const bool message = true):
    BayesCXCIgxs(data, varGenotypic, varResidual, piBetaVal, piPar, estimatePi, piDosageVal[1], piNDCpar, estimatePiGxE, nmale, nfemale, noscale, false),
    snpEffects(data.snpEffectNames),
    piDosage(piDosageVal, piPar),
    piBeta(piBetaVal, piPar),
    piNDC("PiNDC",piDosageVal[1]),
    piFDC("PiFDC",piDosageVal[2]),
    piGxS("PiGxS",piBetaVal[2]),
    deltaNDC("DeltaNDC", data.snpEffectNames),
    deltaFDC("DeltaFDC", data.snpEffectNames),
    deltaGxS("DeltaGxS", data.snpEffectNames) {
        paramSetVec = {&snpEffectsMale, &snpEffectsFemale, &deltaNDC, &deltaFDC, &deltaGxS, &fixedEffects};
        paramVec = {&pi, &nnzSnp, &piNDC, &piFDC, &piGxS, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf};
        paramToPrint = {&pi, &nnzSnp, &piNDC, &piFDC, &piGxS, &sigmaSq, &vargm, &vargf, &varem, &varef, &hsqm, &hsqf, &rounding};
        if (message) {
            cout << "\nBayesXgxs model fitted." << endl;
            cout << "sigmaSq: " << sigmaSq.value << endl;
        }
    }
    
    void sampleUnknowns(void);
};

#endif /* xci_hpp */
