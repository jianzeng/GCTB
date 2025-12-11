//
//  multichain.hpp
//  gctb
//
//  Created by Jian Zeng on 24/11/2024.
//

#include "model.hpp"
#include "options.hpp"

using namespace std;


class MultiChainParameter : public Parameter {
public:
    vector<Parameter*> chainVec;
            
    MultiChainParameter(const string &label, const unsigned numChains): Parameter(label){
        Parameter::numChains = numChains;
        perChainValue.setZero(numChains);
    }
    
    void getValues(void);
};

class MultiChainParamSet : public ParamSet {
public:
    vector<ParamSet*> chainVec;
    
    MultiChainParamSet(const string &label, const vector<string> &header, const unsigned numChains): ParamSet(label, header){
        ParamSet::numChains = numChains;
        perChainValues.setZero(size, numChains);
    }
    
    void getValues(void);
};

class MultiChainParamVec : public vector<MultiChainParameter*> {
public:
    unsigned numParams;
    
    MultiChainParamVec(const string &label, const unsigned numParams, const unsigned numChains): numParams(numParams){
        for (unsigned i=0; i<numParams; ++i) {
            this->push_back(new MultiChainParameter(label + to_string(static_cast<long long>(i + 1)), numChains));
        }
    }

    void getValues(void);
};

class MultiChainParamSetVec : public vector<MultiChainParamSet*> {
public:
    unsigned numParams;
    
    MultiChainParamSetVec(const string &label, const vector<string> &header, const unsigned numParams, const unsigned numChains): numParams(numParams){
        for (unsigned i=0; i<numParams; ++i) {
            this->push_back(new MultiChainParamSet(label + to_string(static_cast<long long>(i + 1)), header, numChains));
        }
    }

    void getValues(void);
};

class MultiChainSBayesR : public ApproxBayesR {
public:
    
    class ChainVecSBayesR : public vector<ApproxBayesR*> {
    public:
        ChainVecSBayesR(const Data &data, const Options &opt){
            for (unsigned i=0; i<opt.numChains; ++i) {
                this->push_back(new ApproxBayesR(data, data.lowRankModel, data.varGenotypic, data.varResidual, opt.pis, opt.piPar, opt.gamma, opt.estimatePi, opt.noscale, opt.hsqPercModel, opt.robustMode, opt.algorithm, false));
            }
        }
    };
    
    class Heritability : public MultiChainParameter {
    public:
        Heritability(const ChainVecSBayesR &chains): MultiChainParameter("hsq", chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->hsq);
            }
        }
    };
    
    class SnpPIP : public MultiChainParamSet {
    public:
        SnpPIP(const vector<string> &header, const ChainVecSBayesR &chains): MultiChainParamSet("PIP", header, chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->snpPip);
            }
        }
    };
    
    class SnpEffects : public MultiChainParamSet {
    public:
        SnpEffects(const vector<string> &header, const ChainVecSBayesR &chains): MultiChainParamSet("SnpEffects", header, chains.size()){
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->snpEffects);
            }
        }
    };
    
    class DeltaPi : public MultiChainParamSetVec {
    public:
        DeltaPi(const vector<string> &header, const unsigned numDist, const ChainVecSBayesR &chains):
        MultiChainParamSetVec("DeltaPi", header, numDist, chains.size()) {
            for (unsigned i=0; i<numDist; ++i) {
                for (unsigned j=0; j<chains.size(); ++j) {
                    (*this)[i]->chainVec.push_back(chains[j]->deltaPi[i]);
                }
            }
        }
    };
    
    class NumSnpMixComps : public MultiChainParamVec {
    public:
        NumSnpMixComps(const unsigned numDist, const ChainVecSBayesR &chains):
        MultiChainParamVec("NumSnp", numDist, chains.size()){
            for (unsigned i=0; i<numDist; ++i) {
                for (unsigned j=0; j<chains.size(); ++j) {
                    (*this)[i]->chainVec.push_back(chains[j]->numSnps[i]);
                }
            }
        }
    };
    
    class VgMixComps : public MultiChainParamVec {
    public:
        VgMixComps(const unsigned numDist, const ChainVecSBayesR &chains):
        MultiChainParamVec("Vg", numDist, chains.size()){
            for (unsigned i = 0; i<numDist; ++i) {
                for (unsigned j=0; j<chains.size(); ++j) {
                    (*this)[i]->chainVec.push_back(chains[j]->Vgs[i]);
                }
            }
        }
    };
    
    class NumBadSnps : public MultiChainParameter {
    public:
        vector<ApproxBayesC::NumBadSnps*> nBadSnpVec;
        set<string> badSnpSet;
        set<unsigned> badSnpIdxSet;

        ofstream out;

        NumBadSnps(const string &title, const ChainVecSBayesR &chains): MultiChainParameter("NumSkeptSnp", chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
               chainVec.push_back(&chains[i]->nBadSnps);
               nBadSnpVec.push_back(&chains[i]->nBadSnps);
               chains[i]->nBadSnps.writeTxt = false;
               chains[i]->nBadSnps.out.close();
            }
            string filename = title + ".skepticalSNPs";
            out.open(filename.c_str());
        }
        void output(void);
    };
    
    class NumHighPIPs : public Parameter {
    public:
        float threshold;
        
        NumHighPIPs(const string &lab = "NumHighPIP"): Parameter(lab){
            threshold = 0.9;
        }
        
        void getValue(const VectorXf &PIP);
    };
    
    class SnpHsqPEP : public MultiChainParamSet {
    public:
        SnpHsqPEP(const vector<string> &header, const ChainVecSBayesR &chains): MultiChainParamSet("PEP", header, chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->snpHsqPep);
            }
        }
    };

    unsigned numChains;
    
    // for nested OMP
    unsigned numThreadTotal;
    unsigned numThreadLevel1;
    unsigned numThreadLevel2;

    ChainVecSBayesR chainVec;
    Heritability hsq;
    SnpPIP pip;
    SnpEffects snpEffects;
    DeltaPi deltaPi;
    NumSnpMixComps numSnpMix;
    VgMixComps vgMix;
    NumBadSnps nBadSnps;
    NumHighPIPs nHighPips;
    SnpHsqPEP snpHsqPep;
    
    MultiChainSBayesR(const Data &data, const Options &opt, const bool message = true):
    ApproxBayesR(data, data.lowRankModel, data.varGenotypic, data.varResidual, opt.pis, opt.piPar, opt.gamma, opt.estimatePi, opt.noscale, opt.hsqPercModel, opt.robustMode, opt.algorithm, false),
    numChains(opt.numChains),
    chainVec(data, opt),
    hsq(chainVec),
    pip(data.snpEffectNames, chainVec),
    snpEffects(data.snpEffectNames, chainVec),
    deltaPi(data.snpEffectNames, opt.gamma.size(), chainVec),
    numSnpMix(opt.gamma.size(), chainVec),
    vgMix(opt.gamma.size(), chainVec),
    nBadSnps(opt.title, chainVec),
    nHighPips(),
    snpHsqPep(data.snpEffectNames, chainVec)
    {
        
        // for nested OMP
        omp_set_max_active_levels(1);  // reset previous nested parallelism if any
        omp_set_max_active_levels(2);  // Enable nested parallelism
        numThreadTotal = omp_get_max_threads();
        numThreadLevel1 = std::min(numChains, numThreadTotal);
        numThreadLevel2 = std::floor(numThreadTotal/numThreadLevel1);
        
        //cout << "numThreadTotal " << numThreadTotal << " numThreadLevel1 " << numThreadLevel1 << " numThreadLevel2 " << numThreadLevel2 << endl;

        paramVec    = {&hsq};
        paramVec.insert(paramVec.end(), numSnpMix.begin(), numSnpMix.end());
        paramVec.insert(paramVec.end(), vgMix.begin(), vgMix.end());
        
        paramSetVec = {&snpEffects, &pip, &snpHsqPep};
        paramSetVec.insert(paramSetVec.end(), deltaPi.begin(), deltaPi.end());
        
        paramToPrint = {&hsq, &nHighPips, &nBadSnps};
        paramToPrint.insert(paramToPrint.begin(), vgMix.begin(), vgMix.end());
        paramToPrint.insert(paramToPrint.begin(), numSnpMix.begin(), numSnpMix.end());

        if (message) {
            cout << "\nMulti-chain SBayesR (" << numChains << " chains)" << endl;
            if (lowRankModel) {
                cout << "Using the low-rank model" << endl;
            }
            cout << "Gamma: " << gamma.values.transpose() << endl;
            if (!hsqPercModel) cout << "The SNP effect prior is a mixture distribution with an unknown variance variable." << endl;
            if (numThreadTotal == 1) {
                cout << "\nSUGGESTION: Enabling multi-threading is recommended when using multiple chains. You can set this by --thread [any value that is the multiple of the number of chains].\n" << endl;
            } else {
                cout << "Using nested multi-threading (" << numThreadTotal << " threads in total):\n  Level 1: " << numThreadLevel1 << " threads\n    Level 2: " << numThreadLevel2 << " threads" << endl;
            }
            cout << endl;
        }
    }
    
    void sampleUnknowns(const unsigned iter);
};


class MultiModelSBayesR : public MultiChainSBayesR {
public:
    
    class ModelVecSBayesR : public vector<ApproxBayesR*> {
    public:
        ModelVecSBayesR(const Data &data, const Options &opt){
            unsigned numModels = opt.numDist - 1;
            VectorXf gamma = opt.gamma;
            VectorXf pis = opt.pis;
            VectorXf piPar = opt.piPar;
            for (unsigned i=0; i<numModels; ++i) {
                this->push_back(new ApproxBayesR(data, data.lowRankModel, data.varGenotypic, data.varResidual, pis, piPar, gamma, opt.estimatePi, opt.noscale, opt.hsqPercModel, opt.robustMode, opt.algorithm, false));
                
                Gadget::removeSecondElement(gamma);
                Gadget::removeSecondElement(pis);
                Gadget::removeSecondElement(piPar);
            }
        }
    };
    
    class Heritability : public MultiChainParameter {
    public:
        Heritability(const ModelVecSBayesR &models): MultiChainParameter("hsq", models.size()){
           for (unsigned i=0; i<numChains; ++i) {
               models[i]->hsq.label += "_M" + to_string(i+1);
                chainVec.push_back(&models[i]->hsq);
            }
        }
    };

    class NumNonZeroSnp : public MultiChainParameter {
    public:
        NumNonZeroSnp(const ModelVecSBayesR &models): MultiChainParameter("NnzSnp", models.size()){
           for (unsigned i=0; i<numChains; ++i) {
               models[i]->nnzSnp.label += "_M" + to_string(i+1);
                chainVec.push_back(&models[i]->nnzSnp);
            }
        }
    };


    unsigned numModels;

    ModelVecSBayesR modelVec;
    Heritability hsq;
    NumNonZeroSnp nnzSnp;

    MultiModelSBayesR(const Data &data, const Options &opt, const bool message = true):
    MultiChainSBayesR(data, opt, false),
    numModels(opt.numDist - 1),
    modelVec(data, opt),
    hsq(modelVec),
    nnzSnp(modelVec)
    {
        
        // for nested OMP
        omp_set_max_active_levels(1);  // reset previous nested parallelism if any
        omp_set_max_active_levels(2);  // Enable nested parallelism
        numThreadTotal = omp_get_max_threads();
        numThreadLevel1 = std::min(numModels, numThreadTotal);
        numThreadLevel2 = std::floor(numThreadTotal/numThreadLevel1);

        paramVec.resize(0);
        paramToPrint.resize(0);
        
        paramVec.insert(paramVec.end(), hsq.chainVec.begin(), hsq.chainVec.end());
        paramVec.insert(paramVec.end(), nnzSnp.chainVec.begin(), nnzSnp.chainVec.end());
        
        paramToPrint.insert(paramToPrint.end(), hsq.chainVec.begin(), hsq.chainVec.end());
        paramToPrint.insert(paramToPrint.end(), nnzSnp.chainVec.begin(), nnzSnp.chainVec.end());

        if (message) {
            cout << "Running " << numModels << " SBayesR models" << endl;
            if (opt.algorithm == "TGS_thin") cout << "Using tempered Gibbs sampling (TGS)" << endl;
            for (unsigned i=0; i<numModels; ++i) {
                unsigned numComp = modelVec[i]->gamma.values.size();
                cout << "  Model " << i+1 << " (M" << i+1 << "): " << numComp << " components with gamma = [";
                for (unsigned k=0; k<numComp; ++k) {
                    if (k==0) cout << modelVec[i]->gamma.values[k];
                    else cout << ", " << modelVec[i]->gamma.values[k];
                }
                cout << "]" << endl;
            }
            if (!hsqPercModel) cout << "The SNP effect prior is a mixture distribution with an unknown variance variable." << endl;
            //cout << endl;
        }

    }

    void sampleUnknowns(const unsigned iter);

};



class MultiChainSBayesRC : public MultiChainSBayesR {
public:
    
    class ChainVecSBayesRC : public vector<ApproxBayesRC*> {
    public:
        ChainVecSBayesRC(const Data &data, const Options &opt){
            for (unsigned i=0; i<opt.numChains; ++i) {
                this->push_back(new ApproxBayesRC(data, data.lowRankModel, data.varGenotypic, data.varResidual, opt.pis, opt.piPar, opt.gamma, opt.estimatePi, opt.noscale, opt.hsqPercModel, opt.robustMode, opt.estimateRsqEnrich, opt.algorithm, false));
            }
        }
    };
    
    class Heritability : public MultiChainParameter {
    public:
        Heritability(const ChainVecSBayesRC &chains): MultiChainParameter("hsq", chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->hsq);
            }
        }
    };
    
    class SnpPIP : public MultiChainParamSet {
    public:
        
        MatrixXf perChainMean;
        MatrixXf perChainSqrMean;
        VectorXf GelmanRubinStat;
        VectorXf posteriorMean;
        unsigned cntSample;
        
        vector<vector<int> > selectedSnpForTGS;   // select SNPs with high GelmanRubin R statistics and non-trivial PIP

        SnpPIP(const vector<string> &header, const ChainVecSBayesRC &chains): MultiChainParamSet("PIP", header, chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->snpPip);
            }
            perChainMean.setZero(size, numChains);
            perChainSqrMean.setZero(size, numChains);
            GelmanRubinStat.setZero(size);
            posteriorMean.setZero(size);
            cntSample = 0;
        }
        
        void computeGelmanRubinStat(const unsigned iter);
        void selectUnconvergedSnps(const vector<LDBlockInfo*> &keptLdBlockInfoVec);
    };
    
    class SnpEffects : public MultiChainParamSet {
    public:
        SnpEffects(const vector<string> &header, const ChainVecSBayesRC &chains): MultiChainParamSet("SnpEffects", header, chains.size()){
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->snpEffects);
            }
        }
    };
    
    class DeltaPi : public MultiChainParamSetVec {
    public:
        DeltaPi(const vector<string> &header, const unsigned numDist, const ChainVecSBayesRC &chains):
        MultiChainParamSetVec("DeltaPi", header, numDist, chains.size()) {
            for (unsigned i=0; i<numDist; ++i) {
                for (unsigned j=0; j<chains.size(); ++j) {
                    (*this)[i]->chainVec.push_back(chains[j]->deltaPi[i]);
                }
            }
        }
    };
    
    class NumSnpMixComps : public MultiChainParamVec {
    public:
        NumSnpMixComps(const unsigned numDist, const ChainVecSBayesRC &chains):
        MultiChainParamVec("NumSnp", numDist, chains.size()){
            for (unsigned i=0; i<numDist; ++i) {
                for (unsigned j=0; j<chains.size(); ++j) {
                    (*this)[i]->chainVec.push_back(chains[j]->numSnps[i]);
                }
            }
        }
    };
    
    class VgMixComps : public MultiChainParamVec {
    public:
        VgMixComps(const unsigned numDist, const ChainVecSBayesRC &chains):
        MultiChainParamVec("Vg", numDist, chains.size()){
            for (unsigned i = 0; i<numDist; ++i) {
                for (unsigned j=0; j<chains.size(); ++j) {
                    (*this)[i]->chainVec.push_back(chains[j]->Vgs[i]);
                }
            }
        }
    };
    
    class NumBadSnps : public MultiChainParameter {
    public:
        vector<ApproxBayesC::NumBadSnps*> nBadSnpVec;
        set<string> badSnpSet;
        set<unsigned> badSnpIdxSet;

        ofstream out;

        NumBadSnps(const string &title, const ChainVecSBayesRC &chains): MultiChainParameter("NumSkeptSnp", chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
               chainVec.push_back(&chains[i]->nBadSnps);
               nBadSnpVec.push_back(&chains[i]->nBadSnps);
               chains[i]->nBadSnps.writeTxt = false;
               chains[i]->nBadSnps.out.close();
            }
            string filename = title + ".skepticalSNPs";
            out.open(filename.c_str());
        }
        void output(void);
    };
    
    class AnnoEffects : public MultiChainParamSetVec {
    public:
        AnnoEffects(const vector<string> &header, const unsigned numDist, const ChainVecSBayesRC &chains):
        MultiChainParamSetVec("AnnoEffects", header, numDist, chains.size()) {
            for (unsigned i=0; i<numDist; ++i) {
                for (unsigned j=0; j<chains.size(); ++j) {
                    (*this)[i]->chainVec.push_back(chains[j]->annoEffects[i]);
                }
            }
        }
    };
        
    class AnnoJointProb : public MultiChainParamSetVec {
    public:
        AnnoJointProb(const vector<string> &header, const unsigned numDist, const ChainVecSBayesRC &chains):
        MultiChainParamSetVec("AnnoJointProb", header, numDist, chains.size()) {
            for (unsigned i=0; i<numDist; ++i) {
                for (unsigned j=0; j<chains.size(); ++j) {
                    (*this)[i]->chainVec.push_back(chains[j]->annoJointProb[i]);
                }
            }
        }
    };

    class AnnoTotalGenVar : public MultiChainParamSet {
    public:
        AnnoTotalGenVar(const vector<string> &header, const ChainVecSBayesRC &chains): MultiChainParamSet("AnnoTotalGenVar", header, chains.size()){
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->annoTotalGenVar);
            }
        }
    };
    
    class AnnoPerSnpHsqEnrichment : public MultiChainParamSet {
    public:
        AnnoPerSnpHsqEnrichment(const vector<string> &header, const ChainVecSBayesRC &chains): MultiChainParamSet("Marginal_Heritability_Enrichment", header, chains.size()){
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->annoPerSnpHsqEnrich);
            }
        }
    };

    class AnnoPerSnpRsqEnrichment : public MultiChainParamSet {
    public:
        AnnoPerSnpRsqEnrichment(const vector<string> &header, const ChainVecSBayesRC &chains): MultiChainParamSet("Marginal_Predictability_Enrichment", header, chains.size()){
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->annoPerSnpRsqEnrich);
            }
        }
    };

    class AnnoJointPerSnpHsqEnrichment : public MultiChainParamSet {
    public:
        AnnoJointPerSnpHsqEnrichment(const vector<string> &header, const ChainVecSBayesRC &chains): MultiChainParamSet("Joint_Heritability_Enrichment", header, chains.size()){
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->annoJointPerSnpHsqEnrich);
            }
        }
    };
    
    class AnnoJointPerSnpRsqEnrichment : public MultiChainParamSet {
    public:
        AnnoJointPerSnpRsqEnrichment(const vector<string> &header, const ChainVecSBayesRC &chains): MultiChainParamSet("Joint_Predictability_Enrichment", header, chains.size()){
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->annoJointPerSnpRsqEnrich);
            }
        }
    };

    class SnpHsqPEP : public MultiChainParamSet {
    public:
        SnpHsqPEP(const vector<string> &header, const ChainVecSBayesRC &chains): MultiChainParamSet("PEP", header, chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->snpHsqPep);
            }
        }
    };

    unsigned numChains;
    
    ChainVecSBayesRC chainVec;
    Heritability hsq;
    SnpPIP pip;
    SnpEffects snpEffects;
    DeltaPi deltaPi;
    NumSnpMixComps numSnpMix;
    VgMixComps vgMix;
    NumBadSnps nBadSnps;
    AnnoEffects annoEffects;
    AnnoJointProb annoJointProb;
    AnnoTotalGenVar annoTotalGenVar;
    AnnoPerSnpHsqEnrichment annoPerSnpHsqEnrich;
    AnnoPerSnpRsqEnrichment annoPerSnpRsqEnrich;
    AnnoJointPerSnpHsqEnrichment annoJointPerSnpHsqEnrich;
    AnnoJointPerSnpRsqEnrichment annoJointPerSnpRsqEnrich;
    MultiChainSBayesR::NumHighPIPs nHighPips;
    SnpHsqPEP snpHsqPep;
    
    const vector<LDBlockInfo*> &keptLdBlockInfoVec;
    
    bool estimateRsqEnrich;
    
    MultiChainSBayesRC(const Data &data, const Options &opt, const bool message = true):
    MultiChainSBayesR(data, opt, false),
    numChains(opt.numChains),
    chainVec(data, opt),
    hsq(chainVec),
    pip(data.snpEffectNames, chainVec),
    snpEffects(data.snpEffectNames, chainVec),
    deltaPi(data.snpEffectNames, opt.gamma.size(), chainVec),
    numSnpMix(opt.gamma.size(), chainVec),
    vgMix(opt.gamma.size(), chainVec),
    nBadSnps(opt.title, chainVec),
    nHighPips(),
    annoEffects(data.annoNames, opt.gamma.size()-1, chainVec),
    annoJointProb(data.annoNames, opt.gamma.size(), chainVec),
    annoTotalGenVar(data.annoNames, chainVec),
    annoPerSnpHsqEnrich(data.annoNames, chainVec),
    annoPerSnpRsqEnrich(data.annoNames, chainVec),
    annoJointPerSnpHsqEnrich(data.annoNames, chainVec),
    annoJointPerSnpRsqEnrich(data.annoNames, chainVec),
    snpHsqPep(data.snpEffectNames, chainVec),
    keptLdBlockInfoVec(data.keptLdBlockInfoVec),
    estimateRsqEnrich(opt.estimateRsqEnrich)
    {

        paramVec    = {&hsq};
        paramVec.insert(paramVec.end(), numSnpMix.begin(), numSnpMix.end());
        paramVec.insert(paramVec.end(), vgMix.begin(), vgMix.end());

        paramSetVec = {&snpEffects, &pip, &snpHsqPep, &annoTotalGenVar, &annoPerSnpHsqEnrich};
        if (estimateRsqEnrich) {
            paramSetVec.push_back(&annoJointPerSnpHsqEnrich);
            paramSetVec.push_back(&annoPerSnpRsqEnrich);
            paramSetVec.push_back(&annoJointPerSnpRsqEnrich);
        }
        paramSetVec.insert(paramSetVec.end(), deltaPi.begin(), deltaPi.end());
        paramSetVec.insert(paramSetVec.end(), annoEffects.begin(), annoEffects.end());
        paramSetVec.insert(paramSetVec.end(), annoJointProb.begin(), annoJointProb.end());

        paramToPrint = {&hsq, &nHighPips, &nBadSnps};
        paramToPrint.insert(paramToPrint.begin(), vgMix.begin(), vgMix.end());
        paramToPrint.insert(paramToPrint.begin(), numSnpMix.begin(), numSnpMix.end());
        
        paramSetToPrint.resize(0);
        paramSetToPrint.insert(paramSetToPrint.end(), annoEffects.begin(), annoEffects.end());
        paramSetToPrint.insert(paramSetToPrint.end(), annoJointProb.begin(), annoJointProb.end());
        paramSetToPrint.push_back(&annoTotalGenVar);
        paramSetToPrint.push_back(&annoPerSnpHsqEnrich);
        if (estimateRsqEnrich) {
            paramSetToPrint.push_back(&annoJointPerSnpHsqEnrich);
            paramSetToPrint.push_back(&annoPerSnpRsqEnrich);
            paramSetToPrint.push_back(&annoJointPerSnpRsqEnrich);
        }

        if (message) {
            cout << "\nMulti-chain SBayesRC (" << numChains << " chains)" << endl;
            if (lowRankModel) {
                cout << "Using the low-rank model" << endl;
            }
            if (opt.algorithm == "TGS_thin") cout << "Using tempered Gibbs sampling (TGS)" << endl;
            cout << "Gamma: " << gamma.values.transpose() << endl;
            if (!hsqPercModel) cout << "The SNP effect prior is a mixture distribution with an unknown variance variable." << endl;
            if (numThreadTotal == 1) {
                cout << "\nSUGGESTION: Enabling multi-threading is recommended when using multiple chains. You can set this by --thread [any value that is the multiple of the number of chains].\n" << endl;
            } else {
                cout << "Using nested multi-threading (" << numThreadTotal << " threads in total):\n  Level 1: " << numThreadLevel1 << " threads\n    Level 2: " << numThreadLevel2 << " threads" << endl;
            }
            cout << endl;
        }
    }
    
    void sampleUnknowns(const unsigned iter);
};




class MultiChainSBayesRD : public MultiChainSBayesRC {
public:
    
    class ChainVecSBayesRD : public vector<ApproxBayesRD*> {
    public:
        ChainVecSBayesRD(const Data &data, const Options &opt){
            for (unsigned i=0; i<opt.numChains; ++i) {
                this->push_back(new ApproxBayesRD(data, data.lowRankModel, data.varGenotypic, data.varResidual, opt.pis, opt.piPar, opt.gamma, opt.estimatePi, opt.noscale, opt.hsqPercModel, opt.robustMode, opt.estimateRsqEnrich, opt.algorithm, false));
            }
        }
    };
    
    class Heritability : public MultiChainParameter {
    public:
        Heritability(const ChainVecSBayesRD &chains): MultiChainParameter("hsq", chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->hsq);
            }
        }
    };
    
    class SnpPIP : public MultiChainParamSet {
    public:
        SnpPIP(const vector<string> &header, const ChainVecSBayesRD &chains): MultiChainParamSet("PIP", header, chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->snpPip);
            }
        }
    };
    
    class SnpEffects : public MultiChainParamSet {
    public:
        SnpEffects(const vector<string> &header, const ChainVecSBayesRD &chains): MultiChainParamSet("SnpEffects", header, chains.size()){
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->snpEffects);
            }
        }
    };
    
    class DeltaPi : public MultiChainParamSetVec {
    public:
        DeltaPi(const vector<string> &header, const unsigned numDist, const ChainVecSBayesRD &chains):
        MultiChainParamSetVec("DeltaPi", header, numDist, chains.size()) {
            for (unsigned i=0; i<numDist; ++i) {
                for (unsigned j=0; j<chains.size(); ++j) {
                    (*this)[i]->chainVec.push_back(chains[j]->deltaPi[i]);
                }
            }
        }
    };
    
    class NumSnpMixComps : public MultiChainParamVec {
    public:
        NumSnpMixComps(const unsigned numDist, const ChainVecSBayesRD &chains):
        MultiChainParamVec("NumSnp", numDist, chains.size()){
            for (unsigned i=0; i<numDist; ++i) {
                for (unsigned j=0; j<chains.size(); ++j) {
                    (*this)[i]->chainVec.push_back(chains[j]->numSnps[i]);
                }
            }
        }
    };
    
    class VgMixComps : public MultiChainParamVec {
    public:
        VgMixComps(const unsigned numDist, const ChainVecSBayesRD &chains):
        MultiChainParamVec("Vg", numDist, chains.size()){
            for (unsigned i = 0; i<numDist; ++i) {
                for (unsigned j=0; j<chains.size(); ++j) {
                    (*this)[i]->chainVec.push_back(chains[j]->Vgs[i]);
                }
            }
        }
    };
    
    class NumBadSnps : public MultiChainParameter {
    public:
        vector<ApproxBayesC::NumBadSnps*> nBadSnpVec;
        set<string> badSnpSet;
        set<unsigned> badSnpIdxSet;
        
        ofstream out;

        NumBadSnps(const string &title, const ChainVecSBayesRD &chains): MultiChainParameter("NumSkeptSnp", chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
               chainVec.push_back(&chains[i]->nBadSnps);
               nBadSnpVec.push_back(&chains[i]->nBadSnps);
               chains[i]->nBadSnps.writeTxt = false;
               chains[i]->nBadSnps.out.close();
            }
            string filename = title + ".skepticalSNPs";
            out.open(filename.c_str());
        }
        void output(void);
    };
    
    class AnnoEffects : public MultiChainParamSetVec {
    public:
        AnnoEffects(const vector<string> &header, const unsigned numDist, const ChainVecSBayesRD &chains):
        MultiChainParamSetVec("AnnoEffects", header, numDist, chains.size()) {
            for (unsigned i=0; i<numDist; ++i) {
                for (unsigned j=0; j<chains.size(); ++j) {
                    (*this)[i]->chainVec.push_back(chains[j]->annoEffects[i]);
                }
            }
        }
    };
        
    class AnnoJointProb : public MultiChainParamSetVec {
    public:
        AnnoJointProb(const vector<string> &header, const unsigned numDist, const ChainVecSBayesRD &chains):
        MultiChainParamSetVec("AnnoJointProb", header, numDist, chains.size()) {
            for (unsigned i=0; i<numDist; ++i) {
                for (unsigned j=0; j<chains.size(); ++j) {
                    (*this)[i]->chainVec.push_back(chains[j]->annoJointProb[i]);
                }
            }
        }
    };

    class AnnoTotalGenVar : public MultiChainParamSet {
    public:
        AnnoTotalGenVar(const vector<string> &header, const ChainVecSBayesRD &chains): MultiChainParamSet("AnnoTotalGenVar", header, chains.size()){
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->annoTotalGenVar);
            }
        }
    };
    
    class AnnoPerSnpHsqEnrichment : public MultiChainParamSet {
    public:
        AnnoPerSnpHsqEnrichment(const vector<string> &header, const ChainVecSBayesRD &chains): MultiChainParamSet("Marginal_Heritability_Enrichment", header, chains.size()){
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->annoPerSnpHsqEnrich);
            }
        }
    };

    class AnnoJointPerSnpHsqEnrichment : public MultiChainParamSet {
    public:
        AnnoJointPerSnpHsqEnrichment(const vector<string> &header, const ChainVecSBayesRD &chains): MultiChainParamSet("Joint_Heritability_Enrichment", header, chains.size()){
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->annoJointPerSnpHsqEnrich);
            }
        }
    };
        
    class SnpHsqPEP : public MultiChainParamSet {
    public:
        SnpHsqPEP(const vector<string> &header, const ChainVecSBayesRD &chains): MultiChainParamSet("PEP", header, chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->snpHsqPep);
            }
        }
    };
    
    class AnnoPi : public MultiChainParameter {
    public:
        AnnoPi(const ChainVecSBayesRD &chains): MultiChainParameter("AnnoPi", chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->piAnno);
            }
        }
    };
    
    class AnnoPIP : public MultiChainParamSet {
    public:
        AnnoPIP(const vector<string> &header, const ChainVecSBayesRD &chains): MultiChainParamSet("AnnoPIP", header, chains.size()){
            for (unsigned i=0; i<numChains; ++i) {
                 chainVec.push_back(&chains[i]->annoPip);
             }
         }
    };
    
    
    ChainVecSBayesRD chainVec;
    Heritability hsq;
    SnpPIP pip;
    SnpEffects snpEffects;
    DeltaPi deltaPi;
    NumSnpMixComps numSnpMix;
    VgMixComps vgMix;
    NumBadSnps nBadSnps;
    AnnoEffects annoEffects;
    AnnoJointProb annoJointProb;
    AnnoTotalGenVar annoTotalGenVar;
    AnnoPerSnpHsqEnrichment annoPerSnpHsqEnrich;
    AnnoJointPerSnpHsqEnrichment annoJointPerSnpHsqEnrich;
    SnpHsqPEP snpHsqPep;
    AnnoPi piAnno;
    AnnoPIP annoPip;
    
    bool estimateRsqEnrich;
    
    MultiChainSBayesRD(const Data &data, const Options &opt, const bool message = true):
    MultiChainSBayesRC(data, opt, false),
    chainVec(data, opt),
    hsq(chainVec),
    pip(data.snpEffectNames, chainVec),
    snpEffects(data.snpEffectNames, chainVec),
    deltaPi(data.snpEffectNames, opt.gamma.size(), chainVec),
    numSnpMix(opt.gamma.size(), chainVec),
    vgMix(opt.gamma.size(), chainVec),
    nBadSnps(opt.title, chainVec),
    annoEffects(data.annoNames, opt.gamma.size()-1, chainVec),
    annoJointProb(data.annoNames, opt.gamma.size(), chainVec),
    annoTotalGenVar(data.annoNames, chainVec),
    annoPerSnpHsqEnrich(data.annoNames, chainVec),
    annoJointPerSnpHsqEnrich(data.annoNames, chainVec),
    snpHsqPep(data.snpEffectNames, chainVec),
    piAnno(chainVec),
    annoPip(data.annoNames, chainVec),
    estimateRsqEnrich(opt.estimateRsqEnrich)
    {
        
        paramVec    = {&hsq, &piAnno};
        paramVec.insert(paramVec.end(), numSnpMix.begin(), numSnpMix.end());
        paramVec.insert(paramVec.end(), vgMix.begin(), vgMix.end());

        paramSetVec = {&snpEffects, &pip, &snpHsqPep, &annoTotalGenVar, &annoPerSnpHsqEnrich, &annoPip};
        paramSetVec.push_back(&annoJointPerSnpHsqEnrich);
        paramSetVec.insert(paramSetVec.end(), deltaPi.begin(), deltaPi.end());
        paramSetVec.insert(paramSetVec.end(), annoEffects.begin(), annoEffects.end());
        paramSetVec.insert(paramSetVec.end(), annoJointProb.begin(), annoJointProb.end());

        paramToPrint = {&hsq, &nHighPips, &nBadSnps, &piAnno};
        paramToPrint.insert(paramToPrint.begin(), vgMix.begin(), vgMix.end());
        paramToPrint.insert(paramToPrint.begin(), numSnpMix.begin(), numSnpMix.end());
        
        paramSetToPrint.resize(0);
        paramSetToPrint.insert(paramSetToPrint.end(), annoEffects.begin(), annoEffects.end());
        paramSetToPrint.insert(paramSetToPrint.end(), annoJointProb.begin(), annoJointProb.end());
        paramSetToPrint.push_back(&annoTotalGenVar);
        paramSetToPrint.push_back(&annoPerSnpHsqEnrich);
        paramSetToPrint.push_back(&annoJointPerSnpHsqEnrich);
        paramSetToPrint.push_back(&annoPip);

        if (message) {
            cout << "\nMulti-chain SBayesRD (" << numChains << " chains)" << endl;
            if (lowRankModel) {
                cout << "Using the low-rank model" << endl;
            }
            cout << "Gamma: " << gamma.values.transpose() << endl;
            if (!hsqPercModel) cout << "The SNP effect prior is a mixture distribution with an unknown variance variable." << endl;
            if (numThreadTotal == 1) {
                cout << "\nSUGGESTION: Enabling multi-threading is recommended when using multiple chains. You can set this by --thread [any value that is the multiple of the number of chains].\n" << endl;
            } else {
                cout << "Using nested multi-threading (" << numThreadTotal << " threads in total):\n  Level 1: " << numThreadLevel1 << " threads\n    Level 2: " << numThreadLevel2 << " threads" << endl;
            }
            cout << endl;
        }
    }
    
    void sampleUnknowns(const unsigned iter);
};

class MultiChainSBayesS : public ApproxBayesS {
public:
    
    class ChainVecSBayesS : public vector<ApproxBayesS*> {
    public:
        ChainVecSBayesS(const Data &data, const Options &opt){
            for (unsigned i=0; i<opt.numChains; ++i) {
                this->push_back(new ApproxBayesS(data, data.lowRankModel, data.varGenotypic, data.varResidual, opt.pi, opt.piAlpha, opt.piBeta, opt.estimatePi, opt.varS, opt.S, opt.algorithm, opt.noscale, false));
            }
        }
    };
    
    class Heritability : public MultiChainParameter {
    public:
        Heritability(const ChainVecSBayesS &chains): MultiChainParameter("hsq", chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->hsq);
            }
        }
    };
    
    class SnpPIP : public MultiChainParamSet {
    public:
        SnpPIP(const vector<string> &header, const ChainVecSBayesS &chains): MultiChainParamSet("PIP", header, chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->snpPip);
            }
        }
    };
    
    class SnpEffects : public MultiChainParamSet {
    public:
        SnpEffects(const vector<string> &header, const ChainVecSBayesS &chains): MultiChainParamSet("SnpEffects", header, chains.size()){
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->snpEffects);
            }
        }
    };
    
    class SParameter : public MultiChainParameter {
    public:
        SParameter(const ChainVecSBayesS &chains): MultiChainParameter("S", chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->S);
            }
        }
    };
    
    class Pi : public MultiChainParameter {
    public:
        Pi(const ChainVecSBayesS &chains): MultiChainParameter("Pi", chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->pi);
            }
        }
    };
    
    class NnzSnp : public MultiChainParameter {
    public:
        NnzSnp(const ChainVecSBayesS &chains): MultiChainParameter("NnzSnp", chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->nnzSnp);
            }
        }
    };
    
    class SigmaSq : public MultiChainParameter {
    public:
        SigmaSq(const ChainVecSBayesS &chains): MultiChainParameter("SigmaSq", chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->sigmaSq);
            }
        }
    };
    
    class GenVar : public MultiChainParameter {
    public:
        GenVar(const ChainVecSBayesS &chains): MultiChainParameter("GenVar", chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->varg);
            }
        }
    };
    
    class ResVar : public MultiChainParameter {
    public:
        ResVar(const ChainVecSBayesS &chains): MultiChainParameter("ResVar", chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->vare);
            }
        }
    };
    
    class NumBadSnps : public MultiChainParameter {
    public:
        vector<ApproxBayesC::NumBadSnps*> nBadSnpVec;
        set<string> badSnpSet;
        set<unsigned> badSnpIdxSet;

        ofstream out;

        NumBadSnps(const string &title, const ChainVecSBayesS &chains): MultiChainParameter("NumSkeptSnp", chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
               chainVec.push_back(&chains[i]->nBadSnps);
               nBadSnpVec.push_back(&chains[i]->nBadSnps);
               chains[i]->nBadSnps.writeTxt = false;
               chains[i]->nBadSnps.out.close();
            }
            string filename = title + ".skepticalSNPs";
            out.open(filename.c_str());
        }
        void output(void);
    };
    
    class NumHighPIPs : public Parameter {
    public:
        float threshold;
        
        NumHighPIPs(const string &lab = "NumHighPIP"): Parameter(lab){
            threshold = 0.9;
        }
        
        void getValue(const VectorXf &PIP);
    };

    unsigned numChains;
    
    // for nested OMP
    unsigned numThreadTotal;
    unsigned numThreadLevel1;
    unsigned numThreadLevel2;

    ChainVecSBayesS chainVec;
    Heritability hsq;
    SnpPIP pip;
    SnpEffects snpEffects;
    SParameter S;
    Pi pi;
    NnzSnp nnzSnp;
    SigmaSq sigmaSq;
    GenVar varg;
    ResVar vare;
    NumBadSnps nBadSnps;
    NumHighPIPs nHighPips;
    
    MultiChainSBayesS(const Data &data, const Options &opt, const bool message = true):
    ApproxBayesS(data, data.lowRankModel, data.varGenotypic, data.varResidual, opt.pi, opt.piAlpha, opt.piBeta, opt.estimatePi, opt.varS, opt.S, opt.algorithm, opt.robustMode, false),
    numChains(opt.numChains),
    chainVec(data, opt),
    hsq(chainVec),
    pip(data.snpEffectNames, chainVec),
    snpEffects(data.snpEffectNames, chainVec),
    S(chainVec),
    pi(chainVec),
    nnzSnp(chainVec),
    sigmaSq(chainVec),
    varg(chainVec),
    vare(chainVec),
    nBadSnps(opt.title, chainVec),
    nHighPips()
    {
        
        // for nested OMP
        omp_set_max_active_levels(1);  // reset previous nested parallelism if any
        omp_set_max_active_levels(2);  // Enable nested parallelism
        numThreadTotal = omp_get_max_threads();
        numThreadLevel1 = std::min(numChains, numThreadTotal);
        numThreadLevel2 = std::floor(numThreadTotal/numThreadLevel1);
        
        paramVec    = {&pi, &nnzSnp, &sigmaSq, &S, &varg, &vare, &hsq};
        paramSetVec = {&snpEffects, &pip};
        
        paramToPrint = {&pi, &nnzSnp, &sigmaSq, &S, &varg, &vare, &hsq, &nHighPips, &nBadSnps};

        if (message) {
            cout << "\nMulti-chain SBayesS (" << numChains << " chains)" << endl;
            if (lowRankModel) {
                cout << "Using the low-rank model" << endl;
            }
            if (numThreadTotal == 1) {
                cout << "\nSUGGESTION: Enabling multi-threading is recommended when using multiple chains. You can set this by --thread [any value that is the multiple of the number of chains].\n" << endl;
            } else {
                cout << "Using nested multi-threading (" << numThreadTotal << " threads in total):\n  Level 1: " << numThreadLevel1 << " threads\n    Level 2: " << numThreadLevel2 << " threads" << endl;
            }
            cout << endl;
        }
    }
    
    void sampleUnknowns(const unsigned iter);
};
