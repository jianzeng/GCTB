//
//  multichains.hpp
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
    
    class ChainVec : public vector<ApproxBayesR*> {
    public:
        ChainVec(const Data &data, const Options &opt){
            for (unsigned i=0; i<opt.numChains; ++i) {
                this->push_back(new ApproxBayesR(data, data.lowRankModel, data.varGenotypic, data.varResidual, opt.pis, opt.piPar, opt.gamma, opt.estimatePi, opt.estimateSigmaSq, opt.noscale, opt.hsqPercModel, opt.overdispersion, opt.estimatePS, opt.spouseCorrelation, opt.diagnosticMode, opt.robustMode, opt.algorithm, opt.nDistAuto, opt.nDistAutoThreshold, false));
            }
        }
    };
    
    class Heritability : public MultiChainParameter {
    public:
        Heritability(const ChainVec &chains): MultiChainParameter("hsq", chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->hsq);
            }
        }
    };
    
    class SnpPIP : public MultiChainParamSet {
    public:
        SnpPIP(const vector<string> &header, const ChainVec &chains): MultiChainParamSet("PIP", header, chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->snpPip);
            }
        }
    };
    
    class SnpEffects : public MultiChainParamSet {
    public:
        SnpEffects(const vector<string> &header, const ChainVec &chains): MultiChainParamSet("SnpEffects", header, chains.size()){
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->snpEffects);
            }
        }
    };
    
    class DeltaPi : public MultiChainParamSetVec {
    public:
        DeltaPi(const vector<string> &header, const unsigned numDist, const ChainVec &chains):
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
        NumSnpMixComps(const unsigned numDist, const ChainVec &chains):
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
        VgMixComps(const unsigned numDist, const ChainVec &chains):
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
        
        ofstream out;

        NumBadSnps(const string &title, const ChainVec &chains): MultiChainParameter("NumBadSnps", chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
               chainVec.push_back(&chains[i]->nBadSnps);
               nBadSnpVec.push_back(&chains[i]->nBadSnps);
               chains[i]->nBadSnps.writeTxt = false;
               chains[i]->nBadSnps.out.close();
            }
            string filename = title + ".badSNPlist";
            out.open(filename.c_str());
        }
        void output(void);
    };
    

    unsigned numChains;
    
    ChainVec chainVec;
    Heritability hsq;
    SnpPIP pip;
    SnpEffects snpEffects;
    DeltaPi deltaPi;
    NumSnpMixComps numSnpMix;
    VgMixComps vgMix;
    NumBadSnps nBadSnps;
    
    MultiChainSBayesR(const Data &data, const Options &opt, const bool message = true):
    ApproxBayesR(data, data.lowRankModel, data.varGenotypic, data.varResidual, opt.pis, opt.piPar, opt.gamma, opt.estimatePi, opt.estimateSigmaSq, opt.noscale, opt.hsqPercModel, opt.overdispersion, opt.estimatePS, opt.spouseCorrelation, opt.diagnosticMode, opt.robustMode, opt.algorithm, opt.nDistAuto, opt.nDistAutoThreshold, false),
    numChains(opt.numChains),
    chainVec(data, opt),
    hsq(chainVec),
    pip(data.snpEffectNames, chainVec),
    snpEffects(data.snpEffectNames, chainVec),
    deltaPi(data.snpEffectNames, opt.gamma.size(), chainVec),
    numSnpMix(opt.gamma.size(), chainVec),
    vgMix(opt.gamma.size(), chainVec),
    nBadSnps(opt.title, chainVec)
    {
        
        paramVec    = {&hsq};
        paramVec.insert(paramVec.end(), numSnpMix.begin(), numSnpMix.end());
        paramVec.insert(paramVec.end(), vgMix.begin(), vgMix.end());
        
        paramSetVec = {&snpEffects, &pip};
        paramSetVec.insert(paramSetVec.end(), deltaPi.begin(), deltaPi.end());
        
        paramToPrint = {&hsq, &nBadSnps};
        paramToPrint.insert(paramToPrint.begin(), vgMix.begin(), vgMix.end());
        paramToPrint.insert(paramToPrint.begin(), numSnpMix.begin(), numSnpMix.end());

        if (message) {
            cout << "\nMulti-chain SBayesR (" << numChains << " chains)" << endl;
            if (lowRankModel) {
                cout << "Using the low-rank model" << endl;
            }
            cout << "Gamma: " << gamma.values.transpose() << endl;
            if (nDistAuto) cout << "The number of mixture components will be automatically assessed at iteration 500." << endl;
            if (!hsqPercModel) cout << "The SNP effect prior is a mixture distribution with an unknown variance variable." << endl;
        }
    }
    
    void sampleUnknowns(void);
};




class MultiChainSBayesRC : public ApproxBayesRC {
public:
    
    class ChainVec : public vector<ApproxBayesRC*> {
    public:
        ChainVec(const Data &data, const Options &opt){
            for (unsigned i=0; i<opt.numChains; ++i) {
                this->push_back(new ApproxBayesRC(data, data.lowRankModel, data.varGenotypic, data.varResidual, opt.pis, opt.piPar, opt.gamma, opt.estimatePi, opt.estimateSigmaSq, opt.noscale, opt.hsqPercModel, opt.perSnpGV, opt.overdispersion, opt.estimatePS, opt.spouseCorrelation, opt.diagnosticMode, opt.robustMode, opt.algorithm, opt.nDistAuto, opt.nDistAutoThreshold, false));
            }
        }
    };
    
    class Heritability : public MultiChainParameter {
    public:
        Heritability(const ChainVec &chains): MultiChainParameter("hsq", chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->hsq);
            }
        }
    };
    
    class SnpPIP : public MultiChainParamSet {
    public:
        SnpPIP(const vector<string> &header, const ChainVec &chains): MultiChainParamSet("PIP", header, chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->snpPip);
            }
        }
    };
    
    class SnpEffects : public MultiChainParamSet {
    public:
        SnpEffects(const vector<string> &header, const ChainVec &chains): MultiChainParamSet("SnpEffects", header, chains.size()){
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->snpEffects);
            }
        }
    };
    
    class DeltaPi : public MultiChainParamSetVec {
    public:
        DeltaPi(const vector<string> &header, const unsigned numDist, const ChainVec &chains):
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
        NumSnpMixComps(const unsigned numDist, const ChainVec &chains):
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
        VgMixComps(const unsigned numDist, const ChainVec &chains):
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
        
        ofstream out;

        NumBadSnps(const string &title, const ChainVec &chains): MultiChainParameter("NumBadSnps", chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
               chainVec.push_back(&chains[i]->nBadSnps);
               nBadSnpVec.push_back(&chains[i]->nBadSnps);
               chains[i]->nBadSnps.writeTxt = false;
               chains[i]->nBadSnps.out.close();
            }
            string filename = title + ".badSNPlist";
            out.open(filename.c_str());
        }
        void output(void);
    };
    
    class AnnoEffects : public MultiChainParamSetVec {
    public:
        AnnoEffects(const vector<string> &header, const unsigned numDist, const ChainVec &chains):
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
        AnnoJointProb(const vector<string> &header, const unsigned numDist, const ChainVec &chains):
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
        AnnoTotalGenVar(const vector<string> &header, const ChainVec &chains): MultiChainParamSet("AnnoTotalGenVar", header, chains.size()){
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->annoTotalGenVar);
            }
        }
    };
    
    class AnnoPerSnpHsqEnrichment : public MultiChainParamSet {
    public:
        AnnoPerSnpHsqEnrichment(const vector<string> &header, const ChainVec &chains): MultiChainParamSet("AnnoPerSnpHsqEnrichment", header, chains.size()){
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->annoPerSnpHsqEnrich);
            }
        }
    };

    unsigned numChains;
    
    ChainVec chainVec;
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
    
    MultiChainSBayesRC(const Data &data, const Options &opt, const bool message = true):
    ApproxBayesRC(data, data.lowRankModel, data.varGenotypic, data.varResidual, opt.pis, opt.piPar, opt.gamma, opt.estimatePi, opt.estimateSigmaSq, opt.noscale, opt.hsqPercModel, opt.perSnpGV, opt.overdispersion, opt.estimatePS, opt.spouseCorrelation, opt.diagnosticMode, opt.robustMode, opt.algorithm, opt.nDistAuto, opt.nDistAutoThreshold, false),
    numChains(opt.numChains),
    chainVec(data, opt),
    hsq(chainVec),
    pip(data.snpEffectNames, chainVec),
    snpEffects(data.snpEffectNames, chainVec),
    deltaPi(data.snpEffectNames, opt.gamma.size(), chainVec),
    numSnpMix(opt.gamma.size(), chainVec),
    vgMix(opt.gamma.size(), chainVec),
    nBadSnps(opt.title, chainVec),
    annoEffects(data.annoNames, opt.gamma.size()-1, chainVec),
    annoJointProb(data.annoNames, opt.gamma.size()-1, chainVec),
    annoTotalGenVar(data.annoNames, chainVec),
    annoPerSnpHsqEnrich(data.annoNames, chainVec)
    {
        
        paramVec    = {&hsq};
        paramVec.insert(paramVec.end(), numSnpMix.begin(), numSnpMix.end());
        paramVec.insert(paramVec.end(), vgMix.begin(), vgMix.end());

        paramSetVec = {&snpEffects, &pip, &annoTotalGenVar, &annoPerSnpHsqEnrich};
        paramSetVec.insert(paramSetVec.end(), deltaPi.begin(), deltaPi.end());
        paramSetVec.insert(paramSetVec.end(), annoEffects.begin(), annoEffects.end());
        paramSetVec.insert(paramSetVec.end(), annoJointProb.begin(), annoJointProb.end());

        paramToPrint = {&hsq, &nBadSnps};
        paramToPrint.insert(paramToPrint.begin(), vgMix.begin(), vgMix.end());
        paramToPrint.insert(paramToPrint.begin(), numSnpMix.begin(), numSnpMix.end());
        
        paramSetToPrint.resize(0);
        paramSetToPrint.insert(paramSetToPrint.begin(), annoEffects.begin(), annoEffects.end());
        paramSetToPrint.insert(paramSetToPrint.begin(), annoJointProb.begin(), annoJointProb.end());
        paramSetToPrint.push_back(&annoTotalGenVar);
        paramSetToPrint.push_back(&annoPerSnpHsqEnrich);

        if (message) {
            cout << "\nMulti-chain SBayesRC (" << numChains << " chains)" << endl;
            if (lowRankModel) {
                cout << "Using the low-rank model" << endl;
            }
            cout << "Gamma: " << gamma.values.transpose() << endl;
            if (nDistAuto) cout << "The number of mixture components will be automatically assessed at iteration 500." << endl;
            if (!hsqPercModel) cout << "The SNP effect prior is a mixture distribution with an unknown variance variable." << endl;
        }
    }
    
    void sampleUnknowns(void);
};
