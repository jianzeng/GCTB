//
//  multichain.hpp
//  gctb
//
//  Created by Jian Zeng on 24/11/2024.
//

#include "model.hpp"
#include "options.hpp"
#include <memory>

using namespace std;

inline const vector<string> &emptyMultiChainHeader(){
    static const vector<string> header;
    return header;
}


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

class MultiChainSBayesR : public Model {
public:
    
    class ChainVecSBayesR : public vector<ApproxBayesR*> {
    public:
        ChainVecSBayesR(const Data &data, const Options &opt, const bool constructChains = true){
            if (!constructChains) return;
            for (unsigned i=0; i<opt.numChains; ++i) {
                this->push_back(new ApproxBayesR(data, data.lowRankModel, data.varGenotypic, data.varResidual, opt.pis, opt.piPar, opt.gamma, opt.estimatePi, opt.noscale, opt.hsqPercModel, opt.robustMode, opt.algorithm, false));
            }
        }

        ~ChainVecSBayesR(){
            for (unsigned i=0; i<this->size(); ++i) delete (*this)[i];
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

    class HeritabilityFixed : public MultiChainParameter {
    public:
        template<typename ChainVec>
        HeritabilityFixed(const ChainVec &chains): MultiChainParameter("hsqFixed", (unsigned)chains.size()){
            for (unsigned i = 0; i < numChains; ++i)
                chainVec.push_back(&chains[i]->hsqFixed);
        }
    };

    class HeritabilityRandom : public MultiChainParameter {
    public:
        template<typename ChainVec>
        HeritabilityRandom(const ChainVec &chains): MultiChainParameter("hsqRandom", (unsigned)chains.size()){
            for (unsigned i = 0; i < numChains; ++i)
                chainVec.push_back(&chains[i]->hsqRandom);
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
        bool enabled;

        NumBadSnps(const string &title, const ChainVecSBayesR &chains, const bool enabled):
        MultiChainParameter("NumSkeptSnp", chains.size()), enabled(enabled){
           for (unsigned i=0; i<numChains; ++i) {
               chainVec.push_back(&chains[i]->nBadSnps);
               nBadSnpVec.push_back(&chains[i]->nBadSnps);
               chains[i]->nBadSnps.writeTxt = false;
               chains[i]->nBadSnps.out.close();
            }
            if (enabled && numChains) {
                string filename = title + ".skepticalSNPs";
                out.open(filename.c_str());
            }
        }
        void output(void);
    };
    
    class NumHighPIPs : public Parameter {
    public:
        float threshold;
        unsigned burnin;
        unsigned thin;
        unsigned previewThin;
        unsigned cntSample;
        bool posteriorPhase;
        bool excludeSkepticalSnps;
        VectorXf runningPip;

        NumHighPIPs(const unsigned burnin, const unsigned thin, const unsigned outputFreq,
                    const float threshold, const bool excludeSkepticalSnps,
                    const string &lab = "NumHighPIP"):
        Parameter(lab), threshold(threshold), burnin(burnin), thin(std::max(1u, thin)),
        previewThin(std::max(1u, std::min(thin, outputFreq))), cntSample(0),
        posteriorPhase(false), excludeSkepticalSnps(excludeSkepticalSnps) {}

        void getValue(const VectorXf &PIP, const unsigned iter, const set<unsigned> &skepticalSnpIdx);
    };
    
    class SnpHsqPEP : public MultiChainParamSet {
    public:
        SnpHsqPEP(const vector<string> &header, const ChainVecSBayesR &chains): MultiChainParamSet("PEP", header, chains.size()){
           for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->snpHsqPep);
            }
        }
    };

    std::unique_ptr<HeritabilityFixed> mcHsqFixed;
    std::unique_ptr<HeritabilityRandom> mcHsqRandom;

    bool lowRankModel;
    bool hsqPercModel;
    VectorXf gammaValues;
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
    
    MultiChainSBayesR(const Data &data, const Options &opt, const bool message = true, const bool constructChains = true):
    lowRankModel(data.lowRankModel),
    hsqPercModel(opt.hsqPercModel),
    gammaValues(opt.gamma),
    numChains(opt.numChains),
    chainVec(data, opt, constructChains),
    hsq(chainVec),
    pip(constructChains ? data.snpEffectNames : emptyMultiChainHeader(), chainVec),
    snpEffects(constructChains ? data.snpEffectNames : emptyMultiChainHeader(), chainVec),
    deltaPi(constructChains ? data.snpEffectNames : emptyMultiChainHeader(), opt.gamma.size(), chainVec),
    numSnpMix(opt.gamma.size(), chainVec),
    vgMix(opt.gamma.size(), chainVec),
    nBadSnps(opt.title, chainVec, opt.removeSkepticalSnps),
    nHighPips(opt.burnin, opt.thin, opt.outputFreq, opt.pipThreshold, opt.removeSkepticalSnps),
    snpHsqPep(constructChains ? data.snpEffectNames : emptyMultiChainHeader(), chainVec)
    {
        bayesType = "R";
        
        // for nested OMP
        omp_set_max_active_levels(1);  // reset previous nested parallelism if any
        omp_set_max_active_levels(2);  // Enable nested parallelism
        numThreadTotal = omp_get_max_threads();
        numThreadLevel1 = std::min(numChains, numThreadTotal);
        numThreadLevel2 = std::floor(numThreadTotal/numThreadLevel1);
        
        //cout << "numThreadTotal " << numThreadTotal << " numThreadLevel1 " << numThreadLevel1 << " numThreadLevel2 " << numThreadLevel2 << endl;

        if (data.numFixedEffectSnps && constructChains) {
            mcHsqFixed.reset(new HeritabilityFixed(chainVec));
            mcHsqRandom.reset(new HeritabilityRandom(chainVec));
        }

        paramVec.clear();
        if (mcHsqFixed) {
            paramVec.push_back(mcHsqFixed.get());
            paramVec.push_back(mcHsqRandom.get());
        }
        paramVec.push_back(&hsq);
        paramVec.insert(paramVec.end(), numSnpMix.begin(), numSnpMix.end());
        paramVec.insert(paramVec.end(), vgMix.begin(), vgMix.end());
        
        paramSetVec = {&snpEffects, &pip, &snpHsqPep};
        paramSetVec.insert(paramSetVec.end(), deltaPi.begin(), deltaPi.end());
        
        paramToPrint.clear();
        paramToPrint.insert(paramToPrint.end(), numSnpMix.begin(), numSnpMix.end());
        paramToPrint.insert(paramToPrint.end(), vgMix.begin(), vgMix.end());
        if (mcHsqFixed) {
            paramToPrint.push_back(mcHsqFixed.get());
            paramToPrint.push_back(mcHsqRandom.get());
        }
        paramToPrint.push_back(&hsq);
        paramToPrint.push_back(&nHighPips);
        paramToPrint.push_back(&nBadSnps);

        if (message) {
            cout << "\nMulti-chain SBayesR (" << numChains << " chains)" << endl;
            if (lowRankModel) {
                cout << "Using the low-rank model" << endl;
            }
            cout << "Gamma: " << gammaValues.transpose() << endl;
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
    void sampleStartVal(void) {
        for (unsigned i=0; i<chainVec.size(); ++i) chainVec[i]->sampleStartVal();
    }
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

        ~ModelVecSBayesR(){
            for (unsigned i=0; i<this->size(); ++i) delete (*this)[i];
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

    MultiModelSBayesR(const Data &data, const Options &opt, const bool message = true, const unsigned shortChainLength = 0):
    MultiChainSBayesR(data, opt, false, false),
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

        // Parent MultiChainSBayesR paramSetVec points at chainVec; we only sample modelVec.
        // MCMC::collectSamples would read stale/wrong-sized ParamSets → crash. Keep model scalars
        // plus model-specific SNP effects for held-out marginal-effect validation.
        paramSetVec.clear();
        paramSetToPrint.clear();
        for (unsigned i=0; i<numModels; ++i) {
            modelVec[i]->snpEffects.label += "_M" + to_string(static_cast<long long>(i + 1));
            paramSetVec.push_back(&modelVec[i]->snpEffects);
        }

        paramVec.resize(0);
        paramToPrint.resize(0);
        
        paramVec.insert(paramVec.end(), hsq.chainVec.begin(), hsq.chainVec.end());
        paramVec.insert(paramVec.end(), nnzSnp.chainVec.begin(), nnzSnp.chainVec.end());
        
        paramToPrint.insert(paramToPrint.end(), hsq.chainVec.begin(), hsq.chainVec.end());
        paramToPrint.insert(paramToPrint.end(), nnzSnp.chainVec.begin(), nnzSnp.chainVec.end());

        if (message) {
            cout << "Running " << numModels << " SBayesR models";
            if (shortChainLength) cout << " with short chains (" << shortChainLength << " iterations)";
            cout << endl;
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
        ChainVecSBayesRC(const Data &data, const Options &opt, const bool constructChains = true){
            if (!constructChains) return;
            for (unsigned i=0; i<opt.numChains; ++i) {
                this->push_back(new ApproxBayesRC(data, data.lowRankModel, data.varGenotypic, data.varResidual,
                                                  opt.pis, opt.piPar, opt.gamma, opt.estimatePi, opt.noscale,
                                                  opt.hsqPercModel, opt.robustMode, false, opt.algorithm, false,
                                                  false, 100, 2, 0.5f, 0.2f, 0.1f,
                                                  opt.estimateContextVar));
            }
        }

        ~ChainVecSBayesRC(){
            for (unsigned i=0; i<this->size(); ++i) delete (*this)[i];
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
        bool enabled;

        NumBadSnps(const string &title, const ChainVecSBayesRC &chains, const bool enabled):
        MultiChainParameter("NumSkeptSnp", chains.size()), enabled(enabled){
           for (unsigned i=0; i<numChains; ++i) {
               chainVec.push_back(&chains[i]->nBadSnps);
               nBadSnpVec.push_back(&chains[i]->nBadSnps);
               chains[i]->nBadSnps.writeTxt = false;
               chains[i]->nBadSnps.out.close();
            }
            if (enabled && numChains) {
                string filename = title + ".skepticalSNPs";
                out.open(filename.c_str());
            }
        }
        void output(void);
    };
    
    class AnnoEffects : public vector<MultiChainParamSet*> {
    public:
        AnnoEffects(const vector<string> &header, const unsigned numDist, const ChainVecSBayesRC &chains):
        vector<MultiChainParamSet*>() {
            for (unsigned i=0; i<numDist; ++i) {
                this->push_back(new MultiChainParamSet("AnnoEffects_p" + to_string(static_cast<long long>(i + 2)), header, chains.size()));
                for (unsigned j=0; j<chains.size(); ++j) {
                    (*this)[i]->chainVec.push_back(chains[j]->annoEffects[i]);
                }
            }
        }

        void getValues(void) {
            for (unsigned i=0; i<this->size(); ++i) (*this)[i]->getValues();
        }
    };
        
    class AnnoJointProb : public vector<MultiChainParamSet*> {
    public:
        AnnoJointProb(const vector<string> &header, const unsigned numDist, const ChainVecSBayesRC &chains):
        vector<MultiChainParamSet*>() {
            for (unsigned i=0; i<numDist; ++i) {
                this->push_back(new MultiChainParamSet("AnnoJointProb_pi" + to_string(static_cast<long long>(i + 1)), header, chains.size()));
                for (unsigned j=0; j<chains.size(); ++j) {
                    (*this)[i]->chainVec.push_back(chains[j]->annoJointProb[i]);
                }
            }
        }

        void getValues(void) {
            for (unsigned i=0; i<this->size(); ++i) (*this)[i]->getValues();
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
        AnnoPerSnpRsqEnrichment(const vector<string> &header, const ChainVecSBayesRC &chains):
        MultiChainParamSet("Marginal_Predictability_Enrichment", header, chains.size()) {}
    };

    class AnnoConditionalHsqEnrichment : public MultiChainParamSet {
    public:
        AnnoConditionalHsqEnrichment(const vector<string> &header, const ChainVecSBayesRC &chains): MultiChainParamSet("Conditional_Heritability_Enrichment", header, chains.size()){
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->annoConditionalHsqEnrich);
            }
        }
    };

    class AnnoConditionalPolygenicityEnrichment : public MultiChainParamSet {
    public:
        AnnoConditionalPolygenicityEnrichment(const vector<string> &header, const ChainVecSBayesRC &chains):
        MultiChainParamSet("Conditional_Polygenicity_Enrichment", header, chains.size()) {
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->annoConditionalPolygenicityEnrich);
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
    AnnoConditionalHsqEnrichment annoConditionalHsqEnrich;
    AnnoConditionalPolygenicityEnrichment annoConditionalPolygenicityEnrich;
    MultiChainSBayesR::NumHighPIPs nHighPips;
    SnpHsqPEP snpHsqPep;

    std::unique_ptr<MultiChainSBayesR::HeritabilityFixed> rcMcHsqFixed;
    std::unique_ptr<MultiChainSBayesR::HeritabilityRandom> rcMcHsqRandom;
    std::unique_ptr<MultiChainParameter> contextVarianceRatioE8;
    std::unique_ptr<MultiChainParameter> contextVarianceE8;
    
    const vector<LDBlockInfo*> &keptLdBlockInfoVec;
    const vector<MatrixXf> &Qblocks;
    const vector<QuantizedEigenBlock> &quantizedEigenBlocks;
    const vector<VectorXf> &contextOriginalNormWeights;
    const MatrixXf &annoMat;
    const vector<AnnoInfo*> &annoInfoVec;

    bool estimateRsqEnrich;
    bool quantizedRsqEnrich;
    unsigned rsqEnrichBurnin;
    unsigned rsqEnrichFreq;
    
    MultiChainSBayesRC(const Data &data, const Options &opt, const bool message = true, const bool constructChains = true):
    MultiChainSBayesR(data, opt, false, false),
    numChains(opt.numChains),
    chainVec(data, opt, constructChains),
    hsq(chainVec),
    pip(constructChains ? data.snpEffectNames : emptyMultiChainHeader(), chainVec),
    snpEffects(constructChains ? data.snpEffectNames : emptyMultiChainHeader(), chainVec),
    deltaPi(constructChains ? data.snpEffectNames : emptyMultiChainHeader(), opt.gamma.size(), chainVec),
    numSnpMix(opt.gamma.size(), chainVec),
    vgMix(opt.gamma.size(), chainVec),
    nBadSnps(opt.title, chainVec, opt.removeSkepticalSnps),
    nHighPips(opt.burnin, opt.thin, opt.outputFreq, opt.pipThreshold, opt.removeSkepticalSnps),
    annoEffects(data.annoNames, opt.gamma.size()-1, chainVec),
    annoJointProb(data.annoNames, opt.gamma.size(), chainVec),
    annoTotalGenVar(data.annoNames, chainVec),
    annoPerSnpHsqEnrich(data.annoNames, chainVec),
    annoPerSnpRsqEnrich(data.annoNames, chainVec),
    annoConditionalHsqEnrich(data.annoNames, chainVec),
    annoConditionalPolygenicityEnrich(data.annoNames, chainVec),
    snpHsqPep(constructChains ? data.snpEffectNames : emptyMultiChainHeader(), chainVec),
    keptLdBlockInfoVec(data.keptLdBlockInfoVec),
    Qblocks(data.Qblocks),
    quantizedEigenBlocks(data.quantizedEigenBlocks),
    contextOriginalNormWeights(data.marginalizedContextOriginalNormWeights),
    annoMat(data.annoMat),
    annoInfoVec(data.annoInfoVec),
    estimateRsqEnrich(opt.estimateRsqEnrich),
    quantizedRsqEnrich(chainVec.size() && chainVec[0]->quantizedLowRankModel),
    rsqEnrichBurnin(opt.burnin),
    rsqEnrichFreq(opt.thin ? opt.thin : 1)
    {
        bayesType = "RC";

        if (estimateRsqEnrich && !data.lowRankModel) {
            throw("Error: --estimate-rsq-enrich requires eigen-decomposed LD data.");
        }
        if (estimateRsqEnrich && opt.numChains < 2) {
            throw("Error: --estimate-rsq-enrich requires at least two MCMC chains.");
        }
        if (data.numFixedEffectSnps && constructChains) {
            rcMcHsqFixed.reset(new MultiChainSBayesR::HeritabilityFixed(chainVec));
            rcMcHsqRandom.reset(new MultiChainSBayesR::HeritabilityRandom(chainVec));
        }
        if (opt.estimateContextVar && constructChains) {
            contextVarianceRatioE8.reset(new MultiChainParameter("CtxRatioE8", chainVec.size()));
            contextVarianceE8.reset(new MultiChainParameter("CtxVarE8", chainVec.size()));
            for (unsigned i=0; i<chainVec.size(); ++i) {
                contextVarianceRatioE8->chainVec.push_back(&chainVec[i]->contextVarianceRatioE8);
                contextVarianceE8->chainVec.push_back(&chainVec[i]->contextVarianceE8);
            }
        }

        paramVec.clear();
        if (rcMcHsqFixed) {
            paramVec.push_back(rcMcHsqFixed.get());
            paramVec.push_back(rcMcHsqRandom.get());
        }
        paramVec.push_back(&hsq);
        if (contextVarianceRatioE8) {
            paramVec.push_back(contextVarianceRatioE8.get());
            paramVec.push_back(contextVarianceE8.get());
        }
        paramVec.insert(paramVec.end(), numSnpMix.begin(), numSnpMix.end());
        paramVec.insert(paramVec.end(), vgMix.begin(), vgMix.end());

        paramSetVec = {&snpEffects, &pip, &snpHsqPep, &annoTotalGenVar, &annoPerSnpHsqEnrich,
                       &annoConditionalHsqEnrich, &annoConditionalPolygenicityEnrich};
        if (estimateRsqEnrich) {
            paramSetVec.push_back(&annoPerSnpRsqEnrich);
        }
        paramSetVec.insert(paramSetVec.end(), deltaPi.begin(), deltaPi.end());
        paramSetVec.insert(paramSetVec.end(), annoEffects.begin(), annoEffects.end());
        paramSetVec.insert(paramSetVec.end(), annoJointProb.begin(), annoJointProb.end());

        paramToPrint.clear();
        paramToPrint.insert(paramToPrint.end(), numSnpMix.begin(), numSnpMix.end());
        paramToPrint.insert(paramToPrint.end(), vgMix.begin(), vgMix.end());
        if (rcMcHsqFixed) {
            paramToPrint.push_back(rcMcHsqFixed.get());
            paramToPrint.push_back(rcMcHsqRandom.get());
        }
        paramToPrint.push_back(&hsq);
        if (contextVarianceRatioE8) {
            paramToPrint.push_back(contextVarianceRatioE8.get());
            paramToPrint.push_back(contextVarianceE8.get());
        }
        paramToPrint.push_back(&nHighPips);
        paramToPrint.push_back(&nBadSnps);
        
        paramSetToPrint.resize(0);
        paramSetToPrint.insert(paramSetToPrint.end(), annoEffects.begin(), annoEffects.end());
        paramSetToPrint.insert(paramSetToPrint.end(), annoJointProb.begin(), annoJointProb.end());
        paramSetToPrint.push_back(&annoTotalGenVar);
        paramSetToPrint.push_back(&annoPerSnpHsqEnrich);
        paramSetToPrint.push_back(&annoConditionalHsqEnrich);
        paramSetToPrint.push_back(&annoConditionalPolygenicityEnrich);
        if (estimateRsqEnrich) {
            paramSetToPrint.push_back(&annoPerSnpRsqEnrich);
        }

        if (message) {
            cout << "\nMulti-chain SBayesRC (" << numChains << " chains)" << endl;
            if (lowRankModel) {
                cout << "Using the low-rank model" << endl;
            }
            if (opt.algorithm == "TGS_thin") cout << "Using tempered Gibbs sampling (TGS)" << endl;
            if (estimateRsqEnrich)
                cout << "Prediction-accuracy enrichment will be evaluated at every retained MCMC sample (every "
                     << rsqEnrichFreq << " iterations after burn-in)." << endl;
            if (data.integrateOutContext) {
                cout << "Using marginalized-context tiled LD with "
                     << (opt.estimateContextVar ? "chain-specific estimated" : "fixed")
                     << " context variance." << endl;
            }
            cout << "Gamma: " << gammaValues.transpose() << endl;
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
    void sampleStartVal(void) {
        for (unsigned i=0; i<chainVec.size(); ++i) chainVec[i]->sampleStartVal();
    }
    void computePredictionAccuracyEnrichment(void);
};




class MultiChainSBayesRD : public MultiChainSBayesRC {
public:
    
    class ChainVecSBayesRD : public vector<ApproxBayesRD*> {
    public:
        ChainVecSBayesRD(const Data &data, const Options &opt){
            bool sampleAnnoEffectsIndep = opt.bayesType == "RD_indep";
            for (unsigned i=0; i<opt.numChains; ++i) {
                this->push_back(new ApproxBayesRD(data, data.lowRankModel, data.varGenotypic, data.varResidual, opt.pis, opt.piPar, opt.gamma, opt.estimatePi, opt.noscale, opt.hsqPercModel, opt.robustMode, opt.estimateRsqEnrich, opt.algorithm, sampleAnnoEffectsIndep, false));
            }
        }

        ~ChainVecSBayesRD(){
            for (unsigned i=0; i<this->size(); ++i) delete (*this)[i];
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
        bool enabled;

        NumBadSnps(const string &title, const ChainVecSBayesRD &chains, const bool enabled):
        MultiChainParameter("NumSkeptSnp", chains.size()), enabled(enabled){
           for (unsigned i=0; i<numChains; ++i) {
               chainVec.push_back(&chains[i]->nBadSnps);
               nBadSnpVec.push_back(&chains[i]->nBadSnps);
               chains[i]->nBadSnps.writeTxt = false;
               chains[i]->nBadSnps.out.close();
            }
            if (enabled && numChains) {
                string filename = title + ".skepticalSNPs";
                out.open(filename.c_str());
            }
        }
        void output(void);
    };
    
    class AnnoEffects : public vector<MultiChainParamSet*> {
    public:
        AnnoEffects(const vector<string> &header, const unsigned numDist, const ChainVecSBayesRD &chains):
        vector<MultiChainParamSet*>() {
            for (unsigned i=0; i<numDist; ++i) {
                this->push_back(new MultiChainParamSet("AnnoEffects_p" + to_string(static_cast<long long>(i + 2)), header, chains.size()));
                for (unsigned j=0; j<chains.size(); ++j) {
                    (*this)[i]->chainVec.push_back(chains[j]->annoEffects[i]);
                }
            }
        }

        void getValues(void) {
            for (unsigned i=0; i<this->size(); ++i) (*this)[i]->getValues();
        }
    };
        
    class AnnoJointProb : public vector<MultiChainParamSet*> {
    public:
        AnnoJointProb(const vector<string> &header, const unsigned numDist, const ChainVecSBayesRD &chains):
        vector<MultiChainParamSet*>() {
            for (unsigned i=0; i<numDist; ++i) {
                this->push_back(new MultiChainParamSet("AnnoJointProb_pi" + to_string(static_cast<long long>(i + 1)), header, chains.size()));
                for (unsigned j=0; j<chains.size(); ++j) {
                    (*this)[i]->chainVec.push_back(chains[j]->annoJointProb[i]);
                }
            }
        }

        void getValues(void) {
            for (unsigned i=0; i<this->size(); ++i) (*this)[i]->getValues();
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

    class AnnoConditionalHsqEnrichment : public MultiChainParamSet {
    public:
        AnnoConditionalHsqEnrichment(const vector<string> &header, const ChainVecSBayesRD &chains): MultiChainParamSet("Conditional_Heritability_Enrichment", header, chains.size()){
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->annoConditionalHsqEnrich);
            }
        }
    };

    class AnnoConditionalPolygenicityEnrichment : public MultiChainParamSet {
    public:
        AnnoConditionalPolygenicityEnrichment(const vector<string> &header, const ChainVecSBayesRD &chains):
        MultiChainParamSet("Conditional_Polygenicity_Enrichment", header, chains.size()) {
            for (unsigned i=0; i<numChains; ++i) {
                chainVec.push_back(&chains[i]->annoConditionalPolygenicityEnrich);
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

    class AnnoPis : public vector<MultiChainParameter*> {
    public:
        AnnoPis(const unsigned numComp, const ChainVecSBayesRD &chains):
        vector<MultiChainParameter*>() {
            for (unsigned i=0; i<numComp; ++i) {
                this->push_back(new MultiChainParameter("AnnoPi_p" + to_string(static_cast<long long>(i + 2)), chains.size()));
                for (unsigned j=0; j<chains.size(); ++j) {
                    (*this)[i]->chainVec.push_back(chains[j]->piAnno[i]);
                }
            }
        }

        void getValues(void) {
            for (unsigned i=0; i<this->size(); ++i) (*this)[i]->getValues();
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

    class AnnoPIPs : public vector<MultiChainParamSet*> {
    public:
        AnnoPIPs(const vector<string> &header, const unsigned numComp, const ChainVecSBayesRD &chains):
        vector<MultiChainParamSet*>() {
            for (unsigned i=0; i<numComp; ++i) {
                this->push_back(new MultiChainParamSet("AnnoPIP_p" + to_string(static_cast<long long>(i + 2)), header, chains.size()));
                for (unsigned j=0; j<chains.size(); ++j) {
                    (*this)[i]->chainVec.push_back(chains[j]->annoPips[i]);
                }
            }
        }

        void getValues(void) {
            for (unsigned i=0; i<this->size(); ++i) (*this)[i]->getValues();
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
    AnnoConditionalHsqEnrichment annoConditionalHsqEnrich;
    AnnoConditionalPolygenicityEnrichment annoConditionalPolygenicityEnrich;
    SnpHsqPEP snpHsqPep;
    AnnoPi piAnno;
    AnnoPis piAnnos;
    AnnoPIP annoPip;
    AnnoPIPs annoPips;
    
    bool estimateRsqEnrich;
    bool sampleAnnoEffectsIndep;
    
    MultiChainSBayesRD(const Data &data, const Options &opt, const bool message = true):
    MultiChainSBayesRC(data, opt, false, false),
    chainVec(data, opt),
    hsq(chainVec),
    pip(data.snpEffectNames, chainVec),
    snpEffects(data.snpEffectNames, chainVec),
    deltaPi(data.snpEffectNames, opt.gamma.size(), chainVec),
    numSnpMix(opt.gamma.size(), chainVec),
    vgMix(opt.gamma.size(), chainVec),
    nBadSnps(opt.title, chainVec, opt.removeSkepticalSnps),
    annoEffects(data.annoNames, opt.gamma.size()-1, chainVec),
    annoJointProb(data.annoNames, opt.gamma.size(), chainVec),
    annoTotalGenVar(data.annoNames, chainVec),
    annoPerSnpHsqEnrich(data.annoNames, chainVec),
    annoConditionalHsqEnrich(data.annoNames, chainVec),
    annoConditionalPolygenicityEnrich(data.annoNames, chainVec),
    snpHsqPep(data.snpEffectNames, chainVec),
    piAnno(chainVec),
    piAnnos(opt.gamma.size()-1, chainVec),
    annoPip(data.annoNames, chainVec),
    annoPips(data.annoNames, opt.gamma.size()-1, chainVec),
    estimateRsqEnrich(opt.estimateRsqEnrich),
    sampleAnnoEffectsIndep(opt.bayesType == "RD_indep")
    {
        bayesType = opt.bayesType;
        
        paramVec    = {&hsq};
        if (sampleAnnoEffectsIndep)
            paramVec.insert(paramVec.end(), piAnnos.begin(), piAnnos.end());
        else
            paramVec.push_back(&piAnno);
        paramVec.insert(paramVec.end(), numSnpMix.begin(), numSnpMix.end());
        paramVec.insert(paramVec.end(), vgMix.begin(), vgMix.end());

        paramSetVec = {&snpEffects, &pip, &snpHsqPep, &annoTotalGenVar, &annoPerSnpHsqEnrich,
                       &annoConditionalHsqEnrich, &annoConditionalPolygenicityEnrich, &annoPip};
        paramSetVec.insert(paramSetVec.end(), annoPips.begin(), annoPips.end());
        paramSetVec.insert(paramSetVec.end(), deltaPi.begin(), deltaPi.end());
        paramSetVec.insert(paramSetVec.end(), annoEffects.begin(), annoEffects.end());
        paramSetVec.insert(paramSetVec.end(), annoJointProb.begin(), annoJointProb.end());

        paramToPrint = {&hsq, &nHighPips};
        paramToPrint.push_back(&nBadSnps);
        if (sampleAnnoEffectsIndep)
            paramToPrint.insert(paramToPrint.end(), piAnnos.begin(), piAnnos.end());
        else
            paramToPrint.push_back(&piAnno);
        paramToPrint.insert(paramToPrint.begin(), vgMix.begin(), vgMix.end());
        paramToPrint.insert(paramToPrint.begin(), numSnpMix.begin(), numSnpMix.end());
        
        paramSetToPrint.resize(0);
        paramSetToPrint.insert(paramSetToPrint.end(), annoEffects.begin(), annoEffects.end());
        paramSetToPrint.insert(paramSetToPrint.end(), annoJointProb.begin(), annoJointProb.end());
        paramSetToPrint.push_back(&annoTotalGenVar);
        paramSetToPrint.push_back(&annoPerSnpHsqEnrich);
        paramSetToPrint.push_back(&annoConditionalHsqEnrich);
        paramSetToPrint.push_back(&annoConditionalPolygenicityEnrich);
        if (sampleAnnoEffectsIndep) {
            paramSetToPrint.insert(paramSetToPrint.end(), annoPips.begin(), annoPips.end());
        } else {
            paramSetToPrint.push_back(&annoPip);
        }

        if (message) {
            cout << "\nMulti-chain SBayesRD (" << numChains << " chains)" << endl;
            if (lowRankModel) {
                cout << "Using the low-rank model" << endl;
            }
            cout << "Gamma: " << gammaValues.transpose() << endl;
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
    void sampleStartVal(void) {
        for (unsigned i=0; i<chainVec.size(); ++i) chainVec[i]->sampleStartVal();
    }
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
        bool enabled;

        NumBadSnps(const string &title, const ChainVecSBayesS &chains, const bool enabled):
        MultiChainParameter("NumSkeptSnp", chains.size()), enabled(enabled){
           for (unsigned i=0; i<numChains; ++i) {
               chainVec.push_back(&chains[i]->nBadSnps);
               nBadSnpVec.push_back(&chains[i]->nBadSnps);
               chains[i]->nBadSnps.writeTxt = false;
               chains[i]->nBadSnps.out.close();
            }
            if (enabled && numChains) {
                string filename = title + ".skepticalSNPs";
                out.open(filename.c_str());
            }
        }
        void output(void);
    };
    
    class NumHighPIPs : public Parameter {
    public:
        float threshold;
        unsigned burnin;
        unsigned thin;
        unsigned previewThin;
        unsigned cntSample;
        bool posteriorPhase;
        bool excludeSkepticalSnps;
        VectorXf runningPip;

        NumHighPIPs(const unsigned burnin, const unsigned thin, const unsigned outputFreq,
                    const float threshold, const bool excludeSkepticalSnps,
                    const string &lab = "NumHighPIP"):
        Parameter(lab), threshold(threshold), burnin(burnin), thin(std::max(1u, thin)),
        previewThin(std::max(1u, std::min(thin, outputFreq))), cntSample(0),
        posteriorPhase(false), excludeSkepticalSnps(excludeSkepticalSnps) {}

        void getValue(const VectorXf &PIP, const unsigned iter, const set<unsigned> &skepticalSnpIdx);
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
    nBadSnps(opt.title, chainVec, opt.removeSkepticalSnps),
    nHighPips(opt.burnin, opt.thin, opt.outputFreq, opt.pipThreshold, opt.removeSkepticalSnps)
    {
        
        // for nested OMP
        omp_set_max_active_levels(1);  // reset previous nested parallelism if any
        omp_set_max_active_levels(2);  // Enable nested parallelism
        numThreadTotal = omp_get_max_threads();
        numThreadLevel1 = std::min(numChains, numThreadTotal);
        numThreadLevel2 = std::floor(numThreadTotal/numThreadLevel1);
        
        paramVec    = {&pi, &nnzSnp, &sigmaSq, &S, &varg, &vare, &hsq};
        paramSetVec = {&snpEffects, &pip};
        
        paramToPrint = {&pi, &nnzSnp, &sigmaSq, &S, &varg, &vare, &hsq, &nHighPips};
        paramToPrint.push_back(&nBadSnps);

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
