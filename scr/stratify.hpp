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
        
        PiStratified(const vector<string> &header, const float pi, const string &lab = "Pi_Stratified"):
        ParamSet(lab, header), alpha(1), beta(19) {
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
        
        void compute(const float genVar, const float resVar, const VectorXf &snpEffects,
                     const vector<SparseVector<float> > &ZPZsp, const vector<AnnoInfo*> &annoInfoVec);
        
        void compute(const VectorXf &snpEffects, const vector<SparseMatrix<float> > &annowiseZPZsp, const vector<VectorXf> &annowiseZPZdiag, const vector<AnnoInfo*> &annoInfoVec, const float genVar, const float resVar);
    };
    
    class TotalHeritabilityEnrichment : public ParamSet {
    public:
        TotalHeritabilityEnrichment(const vector<string> &header, const string &lab = "TotalHsq_enrichment"): ParamSet(lab, header) {}
        
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
                          const vector<AnnoInfo*> &annoInfoVec, VectorXf &scales);
        void hmcSampler(const unsigned annoIdx, const VectorXf &snpEffects, const VectorXf &snp2pq, const VectorXf &snp2pqLog,
                        const float sigmaSq, const float varg, float &scale, float &value);
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
        
        SnpEffects(const vector<string> &header, const VectorXf &snp2pq, const float pi, const unsigned numAnnos):
        ApproxBayesS::SnpEffects(header, snp2pq, pi) {
            wtdSumSq.setZero(numAnnos);
            numNonZeros.setZero(numAnnos);
        }
        
        void sampleFromFC(VectorXf &rcorr, const vector<SparseVector<float> > &ZPZsp, const VectorXf &ZPZdiag,
                          const vector<ChromInfo*> &chromInfoVec, const vector<SnpInfo*> &incdSnpInfoVec,
                          const VectorXf &snp2pq, const VectorXf &LDsamplVar, const unsigned numAnnos,
                          const VectorXf &sigmaSq, const VectorXf &pi, const VectorXf &S,
                          const float varg, const float vare, const float ps, const float overdispersion);
    };
    
public:    
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

    StratApproxBayesS(const Data &data, const float varGenotypic, const float varResidual, const float pival, const bool estimatePi,
                      const float phi, const float overdispersion, const bool estimatePS, const float icrsq,
                      const float varS, const vector<float> &svalue,
                      const string &algorithm, const bool message = true):
    ApproxBayesS(data, varGenotypic, varResidual, pival, estimatePi, phi, overdispersion, estimatePS, icrsq, varS, svalue, algorithm, false),
    snpEffects(data.snpEffectNames, data.snp2pq, pival, data.numAnnos),
    sigmaSqStrat(data.annoNames, data.annoInfoVec, varGenotypic, pival),
    sigmaSqEnrich(data.annoNames),
    piStrat(data.annoNames, pival),
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
        if (message && myMPI::rank==0) {
//            string alg = algorithm;
//            if (alg!="RMH") alg = "HMC (default)";
            cout << "\nAnnotation-stratified summary-data-based BayesS model fitted." << endl;
        }
        makeAnnowiseSparseLDM(data.ZPZsp, data.annoInfoVec, data.incdSnpInfoVec);
    }
    
    void sampleUnknowns(void);
    void makeAnnowiseSparseLDM(const vector<SparseVector<float> > &ZPZsp, const vector<AnnoInfo*> &annoInfoVec, const vector<SnpInfo*> &snpInfoVec);
};

#endif /* stratify_hpp */
