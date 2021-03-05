//
//  stratify.hpp
//  gctb
//
//  Created by Jian Zeng on 25/06/2018.
//  Copyright © 2018 Jian Zeng. All rights reserved.
//

#ifndef stratify_hpp
#define stratify_hpp

#include <stdio.h>
#include "data.hpp"
#include "mcmc.hpp"


class StratApproxBayesS : public ApproxBayesS {  // annotation stratified analysis
public:
    
    class VarEffectStratified : public ParamSet, public Stat::InvChiSq {
    public:
        const float df;
        VectorXf scales;
        
        VarEffectStratified(const vector<string> &header, const vector<AnnoInfo*> annoVec,
                            const float vg, const float pi, const string &lab = "SigmaSq_Stratified"):
        ParamSet(lab, header), df(4) {
            scales.resize(size);
            for (unsigned i=0; i<size; ++i) {
                AnnoInfo *anno = annoVec[i];
                values[i] = vg*anno->fraction/(anno->snp2pq.sum()*pi);
                scales[i] = 0.5*values[i];                
//                cout << i << " " << values[i] << " " << anno->fraction << " " << anno->snp2pq.sum() << " " << pi << endl;
            }
        }
        
        void sampleFromFC(const VectorXf &snpEffSumSq, const VectorXf &numSnpEff);
        void sampleFromPrior(void);
    };
    
    class VarEffectEnrichment : public ParamSet {
    public:
        VarEffectEnrichment(const vector<string> &header, const string &lab = "SigmaSq_Enrichment"): ParamSet(lab, header){}
        
        void compute(const VectorXf &sigmaSqStrat, const float sigmaSq);
    };
    
    class ScaleVarStratified : public ParamSet {
    public:
        ScaleVarStratified(const vector<string> &header, const string &lab = "ScaleVar_Stratified"): ParamSet(lab, header){}
    };

    class PiStratified : public ParamSet, public Stat::Beta {
    public:
        const float alpha;
        const float beta;
        
        PiStratified(const vector<string> &header, const float pi, const float alpha, const float beta, const string &lab = "Pi_Stratified"):
        ParamSet(lab, header), alpha(alpha), beta(beta) {
            values.setConstant(size, pi);
        }
        
        void sampleFromFC(const vector<unsigned> &numSnps, const VectorXf &numSnpEff);
        void sampleFromPrior(void);
    };
    
    class PiEnrichment : public ParamSet {
    public:
        VectorXf expectation;
        
        PiEnrichment(const vector<string> &header, const vector<AnnoInfo*> &annoVec, const string &lab = "Pi_Enrichment"): ParamSet(lab, header) {
            expectation.resize(size);
            for (unsigned i=0; i<size; ++i) {
                AnnoInfo *anno = annoVec[i];
                expectation[i] = anno->fraction;
            }
        }
        
        void compute(const VectorXf &nnzStrat, const float nnzTotal);
    };
    
    class NnzStratified : public ParamSet {  // number of non-zero SNP effects
    public:
        NnzStratified(const vector<string> &header, const string &lab = "Nnz_Stratified"): ParamSet(lab, header) {};
        void getValues(const VectorXf &nnz) {values = nnz;};
    };
    
    class PropNnzStratified : public ParamSet {  // propotion of nnz in each annotation among nnz in the whole genome
    public:
        PropNnzStratified(const vector<string> &header, const string &lab = "PropNnz_Stratified"): ParamSet(lab, header) {};
        void compute(const VectorXf &nnzStrat, const float nnzGw) {values = nnzStrat/nnzGw;};
    };
    
    class HeritabilityStratified : public ParamSet {
    public:
        const unsigned sampleSize;
        
        HeritabilityStratified(const vector<string> &header, const unsigned n, const string &lab = "hsq_Stratified"):
        ParamSet(lab, header), sampleSize(n) {}
        
        void compute(const VectorXf &sigmaSq, const VectorXf &sum2pqSplusOne, const float genVar, const float resVar);
        
        void compute(const vector<VectorXf> &snpEffectPerAnno, const vector<SpMat> &annowiseZPZsp, const vector<VectorXf> &annowiseZPZdiag, const float genVar, const float resVar);
    };
    
    class PropHeritabilityStratified : public ParamSet {
    public:
        PropHeritabilityStratified(const vector<string> &header, const string &lab = "PropHsq_Stratified"): ParamSet(lab, header) {};
        void compute(const VectorXf &hsqStrat, const float hsqGw) {values = hsqStrat/hsqGw;};
    };
    
    class PerSnpHeritabilityEnrichment : public ParamSet {
    public:
        PerSnpHeritabilityEnrichment(const vector<string> &header, const string &lab = "PerSnpHsq_Enrichment"): ParamSet(lab, header) {}
        
        void compute(const VectorXf &hsqStrat, const VectorXf &expFrac, const float hsqTotal);
    };
    
    class PerNzHeritabilityEnrichment : public ParamSet {  // per non-zero effect heritability enrichment
    public:
        PerNzHeritabilityEnrichment(const vector<string> &header, const string &lab = "PerNzHsq_Enrichment"): ParamSet(lab, header) {}
        
        void compute(const VectorXf &hsqStrat, const VectorXf &nnzStrat, const float hsqTotal, const float nnzTotal);
    };
    
    class SpStratified : public ParamSet {
    public:
        const float var;  // prior
        VectorXf stepSize;
        unsigned numSteps;
        
        vector<VectorXf> snp2pqLog;
        vector<BayesS::AcceptanceRate*> ars;
        
        SpStratified(const vector<string> &header, const vector<AnnoInfo*> &annoInfoVec, const float var, const string &lab = "S_Stratified"):
        ParamSet(lab, header), var(var) {
            stepSize.setConstant(size, 0.001);
            numSteps = 100;
            ars.resize(size);
            snp2pqLog.resize(size);
            for (unsigned i=0; i<size; ++i) {
                snp2pqLog[i] = annoInfoVec[i]->snp2pq.array().log();
                ars[i] = new BayesS::AcceptanceRate();
            }
        }
        
        void sampleFromFC(const vector<VectorXf> &snpEffects, const VectorXf &numNonZeros,
                          const VectorXf &sigmaSq, const VectorXf &hsq, const float genVar, const float resVar,
                          const vector<AnnoInfo*> &annoInfoVec, VectorXf &scales, VectorXf &sum2pqSplusOneVec);
        void hmcSampler(const unsigned annoIdx, const VectorXf &snpEffects, const VectorXf &snp2pq, const VectorXf &snp2pqLog,
                        const float sigmaSq, const float varg, float &scale, float &sum2pqSplusOne, float &value);
        float gradientU(const float S, const VectorXf &snpEffects, const float snp2pqLogSum, const VectorXf &snp2pq,
                        const VectorXf &Snp2pqLog, const float sigmaSq);
        float computeU(const float S, const VectorXf &snpEffects, const float snp2pqLogSum, const VectorXf &snp2pq, const float sigmaSq);
        void sampleFromPrior(void);
    };
    
    class SpEnrichment : public ParamSet {
    public:
        SpEnrichment(const vector<string> &header, const string &lab = "S_Enrichment"): ParamSet(lab, header){}
        
        void compute(const VectorXf &Sstrat, const float S);
    };
    
    class SnpAnnoMembership : public ParamSet {
    public:
        VectorXf numAnnoPerSnpVec;
        
        SnpAnnoMembership(const vector<string> &header, const VectorXf &numAnnoVec, const string &lab = "SnpAnnoMembershipDelta"): ParamSet(lab, header), numAnnoPerSnpVec(numAnnoVec){
        }
        
        void getValues(const VectorXf &snpAnnoVec);
    };
    
    class PerSnpPi : public ParamSet, public Stat::Normal {
    public:
        VectorXf beta;  // intercept + annotations
        
        float varProp;  // proposal variance
        float sigmaSq;  // prior variance of betas

        AcceptanceRate ar;

        PerSnpPi(const vector<string> &header, const unsigned numAnnos, const float pival, const string &lab = "PerSnpPi"): ParamSet(lab, header){
            varProp = 0.01;
            sigmaSq = 10;
//            beta.setZero(numAnnos+1);
            beta.setZero(numAnnos);
            beta[0] = log(pival/(1.0f-pival)); // intercept
            values.setConstant(size, pival);
        }
        
        void sampleFromFC(const MatrixXf &annoMat, const MatrixXf &APA, const VectorXf &snpEffects);
        void computeFromAnnoPi(const MatrixXf &annoMat, const VectorXf &piStrat);
    };
    
    class SnpEffects : public ApproxBayesS::SnpEffects {
    public:
        vector<VectorXf> valuesPerAnno;
        
        float sum2pqBetaSq;
        
        VectorXf wtdSumSqPerAnno;
        VectorXf numNonZeroPerAnno;
        VectorXf sum2pqSplusOnePerAnno;
        VectorXf snpAnnoVec;
        VectorXf sum2pqBetaSqAnno;

        SnpEffects(const vector<string> &header, const VectorXf &snp2pq, const float pi, const vector<AnnoInfo*> &annoVec):
        ApproxBayesS::SnpEffects(header, snp2pq, pi) {
            sum2pqBetaSq = 0.0;
            long numAnnos = annoVec.size();
            valuesPerAnno.resize(numAnnos);
            wtdSumSqPerAnno.setZero(numAnnos);
            numNonZeroPerAnno.setZero(numAnnos);
            sum2pqSplusOnePerAnno.setZero(numAnnos);
            sum2pqBetaSqAnno.setZero(numAnnos);
            for (unsigned i=0; i<numAnnos; ++i) {
                valuesPerAnno[i].setZero(annoVec[i]->size);
                sum2pqSplusOnePerAnno[i] = annoVec[i]->snp2pq.sum()*pi;
            }
            snpAnnoVec.setZero(size);
        }
        
        // assuming a mixture distribution in light of overlapping annotations
        void sampleFromFCMixture(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag,
                          const vector<ChromInfo*> &chromInfoVec, const vector<SnpInfo*> &incdSnpInfoVec,
                          const VectorXf &snp2pq, const VectorXf &LDsamplVar, const unsigned numAnnos,
                          const VectorXf &sigmaSq, const VectorXf &pi, const VectorXf &S, const float Sgw,
                          const float varg, const float vare, const float ps, const float overdispersion);
        void sampleFromFCMixture(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag,
                          const VectorXi &windStart, const VectorXi &windSize,
                          const vector<ChromInfo*> &chromInfoVec, const vector<SnpInfo*> &incdSnpInfoVec,
                          const VectorXf &snp2pq, const VectorXf &LDsamplVar, const unsigned numAnnos,
                          const VectorXf &sigmaSq, const VectorXf &pi, const VectorXf &S, const float Sgw,
                          const float varg, const float vare, const float ps, const float overdispersion);

        // assuming a linear model for overlapping annotations
        void sampleFromFCLinear(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag,
                          const vector<ChromInfo*> &chromInfoVec, const vector<SnpInfo*> &incdSnpInfoVec,
                          const VectorXf &snp2pq, const VectorXf &LDsamplVar, const unsigned numAnnos,
                          const VectorXf &sigmaSq, const VectorXf &pi, const VectorXf &S, const float Sgw,
                          const float varg, const float vare, const float ps, const float overdispersion);
        void sampleFromFCLinear(VectorXf &rcorr, const vector<VectorXf> &ZPZ, const VectorXf &ZPZdiag,
                          const VectorXi &windStart, const VectorXi &windSize,
                          const vector<ChromInfo*> &chromInfoVec, const vector<SnpInfo*> &incdSnpInfoVec,
                          const VectorXf &snp2pq, const VectorXf &LDsamplVar, const unsigned numAnnos,
                          const VectorXf &sigmaSq, const VectorXf &pi, const VectorXf &S, const float Sgw,
                          const float varg, const float vare, const float ps, const float overdispersion);
    };
    
    SnpEffects snpEffects;
    SnpAnnoMembership snpAnnoMembership;
    VarEffectStratified sigmaSqStrat;
    VarEffectEnrichment sigmaSqEnrich;
    PiStratified piStrat;
    PiEnrichment piEnrich;
    NnzStratified nnzStrat;
    PropNnzStratified propNnzStrat;
    HeritabilityStratified hsqStrat;
    PropHeritabilityStratified propHsqStrat;
    PerSnpHeritabilityEnrichment perSnpHsqEnrich;
    PerNzHeritabilityEnrichment perNzHsqEnrich;
    SpStratified Sstrat;
    SpEnrichment Senrich;
    
    ScaleVarStratified scaleStrat;
    PerSnpPi perSnpPi;
    
    enum {linear, mixture} model;
    
    StratApproxBayesS(const Data &data, const float varGenotypic, const float varResidual, const float pival, const float piAlpha, const float piBeta, const bool estimatePi,
                      const float phi, const float overdispersion, const bool estimatePS, const float icrsq, const float spouseCorrelation,
                      const float varS, const vector<float> &svalue,
                      const string &algorithm, const bool robustMode, const bool randomStart = false, const bool message = true):
    ApproxBayesS(data, varGenotypic, varResidual, pival, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, svalue, "HMC", false, robustMode, randomStart, false),
    snpEffects(data.snpEffectNames, data.snp2pq, pival, data.annoInfoVec),
    snpAnnoMembership(data.snpAnnoPairNames, data.numAnnoPerSnpVec),
    sigmaSqStrat(data.annoNames, data.annoInfoVec, varGenotypic, pival),
    sigmaSqEnrich(data.annoNames),
    piStrat(data.annoNames, pival, piAlpha, piBeta),
    piEnrich(data.annoNames, data.annoInfoVec),
    nnzStrat(data.annoNames),
    propNnzStrat(data.annoNames),
    hsqStrat(data.annoNames, data.numKeptInds),
    propHsqStrat(data.annoNames),
    perSnpHsqEnrich(data.annoNames),
    perNzHsqEnrich(data.annoNames),
    Sstrat(data.annoNames, data.annoInfoVec, varS),
    Senrich(data.annoNames),
    scaleStrat(data.annoNames),
    perSnpPi(data.snpEffectNames, data.numAnnos, pival)
    {
        if (algorithm == "linear") model = linear;
        else model = mixture;
        paramSetVec = {&snpEffects, &piStrat, &piEnrich, &sigmaSqStrat, &propNnzStrat, &propHsqStrat, &perSnpHsqEnrich, &perNzHsqEnrich, &Sstrat, &Senrich, &snpAnnoMembership};
        paramVec = {&pi, &nnzSnp, &sigmaSq, &S, &vare, &varg, &hsq};
        paramSetToPrint = {&piStrat, &piEnrich, &sigmaSqStrat, &propNnzStrat, &propHsqStrat, &perSnpHsqEnrich, &perNzHsqEnrich, &Sstrat, &Senrich, &snpAnnoMembership};
        paramToPrint = {&pi, &nnzSnp, &sigmaSq, &S, &vare, &varg, &hsq, &rounding};
        if (modelPS) {
            paramVec.push_back(&ps);
            paramToPrint.push_back(&ps);
        }
        if (spouseCorrelation) {
            paramVec.push_back(&covg);
            paramToPrint.push_back(&covg);
        }
        if (message) {
//            string alg = algorithm;
//            if (alg!="RWMH" && alg!="Reg") alg = "HMC";
            cout << "\nAnnotation-stratified summary-data-based BayesS model fitted." << endl;
            if (model == linear) cout << "  Linear model" << endl;
            if (model == mixture) cout << "  Mixture model" << endl;
        }
        if (randomStart) sampleStartVal();
    }
    
    void sampleUnknowns(void);
    void sampleStartVal(void);
};



///// post hoc stratified analysis based on MCMC samples of SNP effects

class PostHocStratifyS : public StratApproxBayesS {
public:
    
    class SnpEffects : public StratApproxBayesS::SnpEffects {
    public:
        
        SnpEffects(const vector<string> &header, const VectorXf &snp2pq, const vector<AnnoInfo*> &annoVec):
        StratApproxBayesS::SnpEffects(header, snp2pq, 0.01, annoVec){
            values.resize(annoVec.size());
        }
        
        void getValues(const SparseVector<float> &snpEffects, const vector<SnpInfo*> &snpInfoVec, const vector<AnnoInfo*> &annoInfoVec, const VectorXf &snp2pq, const VectorXf &S, const float Sgw);
    };
    
    class PiStratified : public StratApproxBayesS::PiStratified {
    public:
        
        PiStratified(const vector<string> &header, const string &lab = "Pi_Stratified"): StratApproxBayesS::PiStratified(header, 0.01, 1, 1, lab){}
        
        void compute(const vector<unsigned> &numSnps, const VectorXf &numSnpEff);
    };

    class VarEffectStratified : public StratApproxBayesS::VarEffectStratified {
    public:
        
        VarEffectStratified(const vector<string> &header, const vector<AnnoInfo*> annoVec, const float vg, const string &lab = "SigmaSq_Stratified"):
        StratApproxBayesS::VarEffectStratified(header, annoVec, vg, 0.01, lab){}
        
        void compute(const VectorXf &snpEffSumSq, const VectorXf &numSnpEff);

    };
    
    SnpEffects snpEffects;
    PiStratified piStrat;
    VarEffectStratified sigmaSqStrat;
    
    const McmcSamples &snpEffectsMcmc;
    const McmcSamples &hsqMcmc;
    
    const unsigned thin;
    
    unsigned iter;
    
    PostHocStratifyS(const Data &data, const McmcSamples &snpEffectsMcmc, const McmcSamples &hsqMcmc, const unsigned thin, const float hsqhat, const bool message = true):
    StratApproxBayesS(data, hsqhat, 1.0-hsqhat, 0.01, 1, 1, true, 0, 0, 0, 0, 0, 1, vector<float>(1,0), "HMC", false),
    snpEffects(data.snpEffectNames, data.snp2pq, data.annoInfoVec),
    piStrat(data.annoNames),
    sigmaSqStrat(data.annoNames, data.annoInfoVec, hsqhat),
    snpEffectsMcmc(snpEffectsMcmc),
    hsqMcmc(hsqMcmc),
    thin(thin)
    {
        iter = 0;
        paramVec = {&pi, &nnzSnp, &sigmaSq, &S, &hsq};
        paramToPrint = {&pi, &nnzSnp, &sigmaSq, &S, &hsq};
        paramSetVec = {&piStrat, &piEnrich, &propNnzStrat, &propHsqStrat, &perSnpHsqEnrich, &perNzHsqEnrich, &Sstrat, &Senrich};
        paramSetToPrint = {&piStrat, &piEnrich, &propNnzStrat, &propHsqStrat, &perSnpHsqEnrich, &perNzHsqEnrich, &Sstrat, &Senrich};
        if (message) {
            cout << "\nPost hoc Annotation-stratified summary-data-based BayesS analysis: " << endl;
        }
    }
    
    void sampleUnknowns(void);
};


class PostHocStratifySMix : public PostHocStratifyS {
public:
    
    class DeltaS : public ApproxBayesSMix::DeltaS {
    public:
        
        float sum;
        VectorXf sumPerAnno;
        
        DeltaS(const vector<string> &header): ApproxBayesSMix::DeltaS(header){}
        
        void getValues(const SparseVector<float> &deltaS, const vector<AnnoInfo*> &annoInfoVec);
    };
    
    class PiS : public BayesC::Pi {
    public:
        
        PiS(const string &lab = "PiS"): BayesC::Pi(0.01, 1, 1, lab){}
    };
    
    class PiSstratified : public PostHocStratifyS::PiStratified {
    public:
        
        PiSstratified(const vector<string> &header, const string &lab = "PiS_Stratified"): PostHocStratifyS::PiStratified(header, lab){}
    };
    
    class PiSenrichment : public StratApproxBayesS::PiEnrichment {
    public:
        
        PiSenrichment(const vector<string> &header, const vector<AnnoInfo*> &annoVec, const string &lab = "PiS_Enrichment"): StratApproxBayesS::PiEnrichment(header, annoVec, lab){}
    };
    
    DeltaS deltaS;
    PiS piS;
    PiSstratified piSstrat;
    PiSenrichment piSenrich;
    
    const McmcSamples &deltaSmcmc;

    PostHocStratifySMix(const Data &data, const McmcSamples &snpEffectsMcmc, const McmcSamples &hsqMcmc, const McmcSamples &deltaSmcmc, const unsigned thin, const float hsqhat, const bool message = true):
    PostHocStratifyS(data, snpEffectsMcmc, hsqMcmc, thin, hsqhat, false),
    deltaSmcmc(deltaSmcmc),
    deltaS(data.snpEffectNames),
    piSstrat(data.annoNames),
    piSenrich(data.annoNames, data.annoInfoVec)
    {
        paramVec = {&pi, &piS, &hsq};
        paramToPrint = {&pi, &piS, &hsq};
        paramSetVec = {&piStrat, &piEnrich, &propNnzStrat, &propHsqStrat, &perSnpHsqEnrich, &perNzHsqEnrich, &piSstrat, &piSenrich};
        paramSetToPrint = {&piStrat, &piEnrich, &propNnzStrat, &propHsqStrat, &perSnpHsqEnrich, &perNzHsqEnrich, &piSstrat, &piSenrich};
        if (message) {
            cout << "\nPost hoc Annotation-stratified summary-data-based BayesSMix analysis: " << endl;
        }

    }

    void sampleUnknowns(void);
};



#endif /* stratify_hpp */
