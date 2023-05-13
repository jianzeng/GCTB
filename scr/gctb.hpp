//
//  gctb.hpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#ifndef amber_hpp
#define amber_hpp

#include <stdio.h>
#include <omp.h>
#include "options.hpp"
#include "data.hpp"
#include "model.hpp"
#include "mcmc.hpp"
#include "hsq.hpp"
#include "predict.hpp"
#include "stratify.hpp"

class GCTB {
public:
    Options &opt;

    GCTB(Options &options): opt(options){};
    
    void inputIndInfo(Data &data, const string &bedFile, const string &phenotypeFile, const string &keepIndFile,
                      const unsigned keepIndMax, const unsigned mphen, const string &covariateFile, const string &randomCovariateFile, const string &residualDiagFile);
    void inputSnpInfo(Data &data, const string &bedFile, const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile,
                      const unsigned includeChr, const bool excludeAmbiguousSNP, const string &skeletonSnpFile, const string &geneticMapFile, const string &ldBlockInfoFile, const unsigned includeBlock,
                      const string &annotationFile, const bool transpose, const string &continuousAnnoFile, const unsigned flank, const string &eQTLFile,
                      const float mafmin, const float mafmax, const bool noscale, const bool readGenotypes);
    void inputSnpInfo(Data &data, const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile,
                      const string &gwasSummaryFile, const string &ldmatrixFile, const unsigned includeChr, const bool excludeAmbiguousSNP,
                      const string &skeletonSnpFile, const string &geneticMapFile, const float genMapN, const string &annotationFile, const bool transpose, const string &continuousAnnoFile, const unsigned flank, const string &eQTLFile, const string &ldscoreFile, const string &windowFile,
                      const bool multiLDmatrix, const bool excludeMHC, const float afDiff, const float mafmin, const float mafmax, const float pValueThreshold, const float rsqThreshold, const bool sampleOverlap, const bool imputeN, const bool noscale, const bool binSnp, const bool readLDMfromTxtFile);
    // this function read eigen matrices
    void inputSnpInfo(Data &data, const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile,
                            const string &gwasSummaryFile, const string &eigenMatrixFile, const string &ldBlockInfoFile,
                            const unsigned includeChr, const bool excludeAmbiguousSNP,
                            const string &annotationFile, const bool transpose,
                            const string &continuousAnnoFile, const unsigned flank, const string &eQTLFile, const string &ldscoreFile,
                            const float eigenCutoff, const bool excludeMHC,
                            const float afDiff, const float mafmin, const float mafmax, const float pValueThreshold, const float rsqThreshold,
                      const bool sampleOverlap, const bool imputeN, const bool noscale, const bool readLDMfromTxtFile, const bool imputeSummary, const unsigned includeBlock);
    
    void truncBlockEigen(float cutShresh=1e-6);
    
    void inputSnpInfo(Data &data, const string &bedFile, const string &gwasSummaryFile, const float afDiff, const float mafmin, const float mafmax, const float pValueThreshold, const bool sampleOverlap, const bool imputeN, const bool noscale);

    Model* buildModel(Data &data, const string &bedFile, const string &gwasFile, const string &bayesType, const unsigned windowWidth,
                      const float heritability, const float propVarRandom, const float pi, const float piAlpha, const float piBeta, const bool estimatePi, const bool noscale,
                      const VectorXf &pis, const VectorXf &piPar, const VectorXf &gamma, const bool estimateSigmaSq,
                      const float phi, const float kappa, const string &algorithm, const unsigned snpFittedPerWindow,
                      const float varS, const vector<float> &S, const float overdispersion, const bool estimatePS,
                      const float icrsq, const float spouseCorrelation, const bool diagnosticMode, const bool originalModel, const bool perSnpGV, const bool robustMode);
    vector<McmcSamples*> runMcmc(Model &model, const unsigned chainLength, const unsigned burnin, const unsigned thin, const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior);
    void saveMcmcSamples(const vector<McmcSamples*> &mcmcSampleVec, const string &filename);
    void outputResults(const Data &data, const vector<McmcSamples*> &mcmcSampleVec, const string &bayesType, const bool noscale, const string &filename);

    McmcSamples* inputMcmcSamples(const string &mcmcSampleFile, const string &label, const string &fileformat);
    void estimateHsq(const Data &data, const McmcSamples &snpEffects, const McmcSamples &resVar, const string &filename, const unsigned outputFreq);
    void estimatePi(const Data &data, const McmcSamples &snpEffects, const McmcSamples &genVar, const string &filename, const unsigned outputFreq);
    void predict(const Data &data, const string &filename);

    void clearGenotypes(Data &data);
    void stratify(Data &data, const string &ldmatrixFile, const bool multiLDmat, const string &geneticMapFile, const float genMapN, const string &snpResFile, const string &mcmcSampleFile, const string &annotationFile, const bool transpose, const string &continuousAnnoFile, const unsigned flank, const string &eQTLFile, const string &gwasSummaryFile, const float pValueThreshold, const bool imputeN, const string &filename, const string &bayesType, unsigned chainLength, unsigned burnin, const unsigned thin, const unsigned outputFreq);
    
    vector<McmcSamples*> multi_chain_mcmc(Data &data, const string &bayesType, const unsigned windowWidth, const float heritability, const float propVarRandom, const float pi, const float piAlpha, const float piBeta, const bool estimatePi, const VectorXf &pis, const VectorXf &gamma, const float phi, const float kappa, const string &algorithm, const unsigned snpFittedPerWindow, const float varS, const vector<float> &S, const float overdispersion, const bool estimatePS, const float icrsq, const float spouseCorrelation, const bool diagnosticMode, const bool robustMode, const unsigned numChains, const unsigned chainLength, const unsigned burnin, const unsigned thin, const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior);
    
    void solveSnpEffectsByConjugateGradientMethod(Data &data, const float lambda, const string &filename) const;
    
    void pip2p(const Data &data, const VectorXf &pip, const float propNull, VectorXf &pval);
    
    float tuneEigenCutoff(Data &data, const Options &opt);
};

#endif /* amber_hpp */
