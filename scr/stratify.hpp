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
    
    class HeritabilityStratified : public ParamSet {
    public:
        const unsigned sampleSize;
        
        HeritabilityStratified(const vector<string> &header, const unsigned n, const string &lab = "hsq_Stratified"):
        ParamSet(lab, header), sampleSize(n) {}
        
        void compute(const VectorXf &sigmaSq, const VectorXf &sum2pqSplusOne, const float genVar, const float resVar);
        
        void compute(const float genVar, const float resVar, const VectorXf &snpEffects,
                     const vector<SparseVector<float> > &ZPZsp, const vector<AnnoInfo*> &annoInfoVec);
        
        void compute(const VectorXf &snpEffects, const vector<SparseMatrix<float> > &annowiseZPZsp, const vector<VectorXf> &annowiseZPZdiag, const vector<AnnoInfo*> &annoInfoVec, const float genVar, const float resVar);
    };
    
    class TotalHeritabilityEnrichment : public ParamSet {
    public:
        TotalHeritabilityEnrichment(const vector<string> &header, const string &lab = "TotalHsq_Enrichment"): ParamSet(lab, header) {}
        
        void compute(const VectorXf &hsqStrat, const VectorXf &expFrac, const float hsqTotal);
    };
    
    class PerSnpHeritabilityEnrichment : public ParamSet {
    public:
        PerSnpHeritabilityEnrichment(const vector<string> &header, const string &lab = "PerSnpHsq_Enrichment"): ParamSet(lab, header) {}
        
        void compute(const VectorXf &hsqStrat, const VectorXf &nnzStrat, const float hsqTotal, const float nnzTotal);
    };
    
    class SpStratified : public ParamSet {
    public:
        const float var;  // prior
        VectorXf stepSize;
        unsigned numSteps;
        
        const VectorXf snp2pqLog;
        vector<BayesS::AcceptanceRate*> ars;
        
        SpStratified(const vector<string> &header, const VectorXf &snp2pq, const unsigned numAnnos, const float var, const string &lab = "S_Stratified"):
        ParamSet(lab, header), snp2pqLog(snp2pq.array().log()), var(var) {
            stepSize.setConstant(numAnnos, 0.001);
            numSteps = 100;
            ars.resize(numAnnos);
            for (unsigned i=0; i<numAnnos; ++i) ars[i] = new BayesS::AcceptanceRate();
        }
        
        void sampleFromFC(const VectorXf &snpEffects, const VectorXf &numNonZeros, const VectorXf &snp2pq,
                          const VectorXf &sigmaSq, const VectorXf &hsq, const float genVar, const float resVar,
                          const vector<AnnoInfo*> &annoInfoVec, VectorXf &scales, VectorXf &sum2pqSplusOneVec);
        void hmcSampler(const unsigned annoIdx, const VectorXf &snpEffects, const VectorXf &snp2pq, const VectorXf &snp2pqLog,
                        const float sigmaSq, const float varg, float &scale, float &sum2pqSplusOne, float &value);
        float gradientU(const float S, const VectorXf &snpEffects, const float snp2pqLogSum, const VectorXf &snp2pq,
                        const VectorXf &Snp2pqLog, const float sigmaSq);
        float computeU(const float S, const VectorXf &snpEffects, const float snp2pqLogSum, const VectorXf &snp2pq, const float sigmaSq);
    };
    
    class SpEnrichment : public ParamSet {
    public:
        SpEnrichment(const vector<string> &header, const string &lab = "S_Enrichment"): ParamSet(lab, header){}
        
        void compute(const VectorXf &Sstrat, const float S);
    };
    
    class SnpEffects : public ApproxBayesS::SnpEffects {
    public:
        VectorXf wtdSumSq;
        VectorXf numNonZeros;
        VectorXf sum2pqSplusOneVec;
        
        SnpEffects(const vector<string> &header, const VectorXf &snp2pq, const float pi, const vector<AnnoInfo*> &annoVec):
        ApproxBayesS::SnpEffects(header, snp2pq, pi) {
            long numAnnos = annoVec.size();
            wtdSumSq.setZero(numAnnos);
            numNonZeros.setZero(numAnnos);
            sum2pqSplusOneVec.setZero(numAnnos);
            for (unsigned i=0; i<numAnnos; ++i) {
                sum2pqSplusOneVec[i] = annoVec[i]->snp2pq.sum()*pi;
            }
        }
        
        void sampleFromFC(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag,
                          const vector<ChromInfo*> &chromInfoVec, const vector<SnpInfo*> &incdSnpInfoVec,
                          const VectorXf &snp2pq, const VectorXf &LDsamplVar, const unsigned numAnnos,
                          const VectorXf &sigmaSq, const VectorXf &pi, const VectorXf &S,
                          const float varg, const float vare, const float ps, const float overdispersion);
    };
    
    vector<SparseMatrix<float> > annowiseZPZsp;
    vector<VectorXf> annowiseZPZdiag;
    
    SnpEffects snpEffects;
    VarEffectStratified sigmaSqStrat;
    VarEffectEnrichment sigmaSqEnrich;
    PiStratified piStrat;
    PiEnrichment piEnrich;
    NnzStratified nnzStrat;
    HeritabilityStratified hsqStrat;
    TotalHeritabilityEnrichment totalHsqEnrich;
    PerSnpHeritabilityEnrichment perSnpHsqEnrich;
    SpStratified Sstrat;
    SpEnrichment Senrich;
    
    ScaleVarStratified scaleStrat;
    
    StratApproxBayesS(const Data &data, const float varGenotypic, const float varResidual, const float pival, const float piAlpha, const float piBeta, const bool estimatePi,
                      const float phi, const float overdispersion, const bool estimatePS, const float icrsq, const float spouseCorrelation,
                      const float varS, const vector<float> &svalue,
                      const string &algorithm, const bool message = true):
    ApproxBayesS(data, varGenotypic, varResidual, pival, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, svalue, algorithm, false, false),
    snpEffects(data.snpEffectNames, data.snp2pq, pival, data.annoInfoVec),
    sigmaSqStrat(data.annoNames, data.annoInfoVec, varGenotypic, pival),
    sigmaSqEnrich(data.annoNames),
    piStrat(data.annoNames, pival, piAlpha, piBeta),
    piEnrich(data.annoNames, data.annoInfoVec),
    nnzStrat(data.annoNames),
    hsqStrat(data.annoNames, data.numKeptInds),
    totalHsqEnrich(data.annoNames),
    perSnpHsqEnrich(data.annoNames),
    Sstrat(data.annoNames, data.snp2pq, data.numAnnos, varS),
    Senrich(data.annoNames),
    scaleStrat(data.annoNames)
    {
        paramSetVec = {&snpEffects, &piStrat, &piEnrich, &nnzStrat, &sigmaSqStrat, &sigmaSqEnrich, &hsqStrat, &totalHsqEnrich, &perSnpHsqEnrich, &Sstrat, &Senrich};
        paramVec = {&pi, &nnzSnp, &sigmaSq, &S, &vare, &varg, &hsq};
        paramSetToPrint = {&piStrat, &piEnrich, &nnzStrat, &sigmaSqStrat, &sigmaSqEnrich, &hsqStrat, &totalHsqEnrich, &perSnpHsqEnrich, &Sstrat, &Senrich};
        paramToPrint = {&pi, &nnzSnp, &sigmaSq, &S, &vare, &varg, &hsq, &rounding};
        if (modelPS) {
            paramVec.push_back(&ps);
            paramToPrint.push_back(&ps);
        }
        if (spouseCorrelation) {
            paramVec.push_back(&covg);
            paramToPrint.push_back(&covg);
        }
        if (message && myMPI::rank==0) {
//            string alg = algorithm;
//            if (alg!="RMH") alg = "HMC (default)";
            cout << "\nAnnotation-stratified summary-data-based BayesS model fitted." << endl;
        }
//        makeAnnowiseSparseLDM(data.ZPZsp, data.annoInfoVec, data.incdSnpInfoVec);
    }
    
    void sampleUnknowns(void);
    void makeAnnowiseSparseLDM(const vector<SparseVector<float> > &ZPZsp, const vector<AnnoInfo*> &annoInfoVec, const vector<SnpInfo*> &snpInfoVec);
};



///// post hoc stratified analysis based on MCMC samples of SNP effects

class PostHocStratify : public StratApproxBayesS {
public:
    
    class SnpEffects : public StratApproxBayesS::SnpEffects {
    public:
        vector<VectorXf> values;
        
        SnpEffects(const vector<string> &header, const VectorXf &snp2pq, const float pi, const vector<AnnoInfo*> &annoVec):
        StratApproxBayesS::SnpEffects(header, snp2pq, pi, annoVec){
            values.resize(annoVec.size());
        }
        
        void getValues(const SparseVector<float> &snpEffects, const vector<AnnoInfo*> &annoInfoVec, const VectorXf &snp2pq, const VectorXf &S);
    };

    class Heritability : public StratApproxBayesS::HeritabilityStratified {
    public:
        
        Heritability(const vector<string> &header, const unsigned n):HeritabilityStratified(header, n){}
        
        void compute(const vector<VectorXf> &snpEffects, const vector<SparseMatrix<float> > &annowiseZPZsp,
                                 const vector<VectorXf> &annowiseZPZdiag, const float genVar, const float resVar);
    };
    
    class Pi : public StratApproxBayesS::PiStratified {
    public:
        
        Pi(const vector<string> &header, const float pi, const float alpha, const float beta): StratApproxBayesS::PiStratified(header, pi, alpha, beta){}

        void compute(const vector<unsigned> &numSnps, const VectorXf &numSnpEff);
    };
    
    class VarEffects : public StratApproxBayesS::VarEffectStratified {
    public:
        
        VarEffects(const vector<string> &header, const vector<AnnoInfo*> annoVec, const float vg, const float pi):
        StratApproxBayesS::VarEffectStratified(header, annoVec, vg, pi){}
        
        void compute(const VectorXf &snpEffSumSq, const VectorXf &numSnpEff);

    };
    
    class Sp : public StratApproxBayesS::SpStratified {
    public:
        
        vector<VectorXf> snp2pqLog;
        
        Sp(const vector<string> &header, const VectorXf &snp2pq, const vector<AnnoInfo*> &annoInfoVec, const float var): StratApproxBayesS::SpStratified(header, snp2pq, annoInfoVec.size(), var){
            snp2pqLog.resize(size);
            for (unsigned i=0; i<size; ++i) {
                snp2pqLog[i] = annoInfoVec[i]->snp2pq.array().log();
            }
        }
        
        void sampleFromFC(const vector<VectorXf> &snpEffects,  const VectorXf &numNonZeros, const VectorXf &sigmaSq, const VectorXf &hsq, const float genVar, const float resVar, const vector<AnnoInfo*> &annoInfoVec, VectorXf &scales, VectorXf &sum2pqSplusOneVec);
    };
    
    SnpEffects snpEffects;
    Heritability hsqStrat;
    Pi piStrat;
    VarEffects sigmaSqStrat;
    Sp Sstrat;
    
    const McmcSamples &snpEffectsMcmc;
    const McmcSamples &vargMcmc;
    const McmcSamples &vareMcmc;
    const McmcSamples &sigmaSqMcmc;
    const McmcSamples &piMcmc;
    
    const unsigned thin;
    
    unsigned iter;
    
    PostHocStratify(const Data &data, const McmcSamples &snpEffectsMcmc, const McmcSamples &vargMcmc, const McmcSamples &vareMcmc, const McmcSamples &sigmaSqMcmc, const McmcSamples &piMcmc, const unsigned thin, const float varGenotypic, const float varResidual, const float pival, const float piAlpha, const float piBeta, const float varS, const vector<float> &svalue, const bool message = true):
    StratApproxBayesS(data, varGenotypic, varResidual, pival, piAlpha, piBeta, true, 0, 0, 0, 0, 0, varS, svalue, "HMC", false),
    snpEffects(data.snpEffectNames, data.snp2pq, pival, data.annoInfoVec),
    hsqStrat(data.annoNames, data.numKeptInds),
    piStrat(data.annoNames, pival, piAlpha, piBeta),
    sigmaSqStrat(data.annoNames, data.annoInfoVec, varGenotypic, pival),
    Sstrat(data.annoNames, data.snp2pq, data.annoInfoVec, varS),
    snpEffectsMcmc(snpEffectsMcmc),
    vargMcmc(vargMcmc),
    vareMcmc(vareMcmc),
    sigmaSqMcmc(sigmaSqMcmc),
    piMcmc(piMcmc),
    thin(thin)
    {
        iter = 0;
        paramVec.clear();
        paramToPrint.clear();
        paramSetVec = {&piStrat, &piEnrich, &nnzStrat, &sigmaSqStrat, &sigmaSqEnrich, &hsqStrat, &totalHsqEnrich, &perSnpHsqEnrich, &Sstrat};
        paramSetToPrint = {&piStrat, &piEnrich, &nnzStrat, &sigmaSqStrat, &sigmaSqEnrich, &hsqStrat, &totalHsqEnrich, &perSnpHsqEnrich, &Sstrat};
        if (message && myMPI::rank==0) {
            cout << "\nPost hoc Annotation-stratified summary-data-based BayesS analysis: " << endl;
        }
    }
    
    void sampleUnknowns(void);
};






#endif /* stratify_hpp */
