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
#include <mpi.h>
#include <omp.h>
#include "options.hpp"
#include "data.hpp"
#include "model.hpp"
#include "mcmc.hpp"
#include "hsq.hpp"
#include "predict.hpp"
#include "mympi.hpp"
#include "stratify.hpp"

class GCTB {
public:
    Options &opt;

    GCTB(Options &options): opt(options){};
    
    void inputIndInfo(Data &data, const string &bedFile, const string &phenotypeFile, const string &keepIndFile,
                      const unsigned keepIndMax, const unsigned mphen, const string &covariateFile);
    void inputSnpInfo(Data &data, const string &bedFile, const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile,
                      const unsigned includeChr, const bool excludeAmbiguousSNP, const string &skeletonSnpFile, const string &geneticMapFile, const float mafmin, const float mafmax, const bool noscale, const bool readGenotypes);
    void inputSnpInfo(Data &data, const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile,
                      const string &gwasSummaryFile, const string &ldmatrixFile, const unsigned includeChr, const bool excludeAmbiguousSNP,
                      const string &skeletonSnpFile, const string &geneticMapFile, const float genMapN, const string &annotationFile, const bool transpose, const string &continuousAnnoFile, const unsigned flank, const string &eQTLFile, const string &ldscoreFile,
                      const bool multiLDmatrix, const bool excludeMHC, const float afDiff, const float mafmin, const float mafmax, const bool sampleOverlap, const bool imputeN, const string &bayesType, const bool noscale);
    
    void inputSnpInfo(Data &data, const string &bedFile, const string &gwasSummaryFile, const float afDiff, const float mafmin, const float mafmax, const bool sampleOverlap, const bool imputeN, const string &bayesType, const bool noscale);

    Model* buildModel(Data &data, const string &bedFile, const string &gwasFile, const string &bayesType, const unsigned windowWidth,
                      const float heritability, const float pi, const float piAlpha, const float piBeta, const bool estimatePi, const bool noscale, const VectorXf &pis, const VectorXf &gamma,
                      const float phi, const float kappa, const string &algorithm, const unsigned snpFittedPerWindow,
                      const float varS, const vector<float> &S, const float overdispersion, const bool estimatePS,
                      const float icrsq, const float spouseCorrelation, const bool diagnosticMode, const bool originalModel);
    vector<McmcSamples*> runMcmc(Model &model, const unsigned chainLength, const unsigned burnin, const unsigned thin, const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior);
    void saveMcmcSamples(const vector<McmcSamples*> &mcmcSampleVec, const string &filename);
    void outputResults(const Data &data, const vector<McmcSamples*> &mcmcSampleVec, const string &bayesType, const bool noscale, const string &filename);

    McmcSamples* inputMcmcSamples(const string &mcmcSampleFile, const string &label, const string &fileformat);
    void estimateHsq(const Data &data, const McmcSamples &snpEffects, const string &filename);
    void predict(const Data &data, const string &filename);

    void clearGenotypes(Data &data);
    void stratify(Data &data, const string &ldmatrixFile, const bool multiLDmat, const string &geneticMapFile, const float genMapN, const string &snpResFile, const string &mcmcSampleFile, const string &annotationFile, const bool transpose, const string &continuousAnnoFile, const unsigned flank, const string &eQTLFile, const string &gwasSummaryFile, const bool imputeN, const string &filename, const string &bayesType, unsigned chainLength, unsigned burnin, const unsigned thin, const unsigned outputFreq);
    
    vector<McmcSamples*> multi_chain_mcmc(Data &data, const string &bayesType, const unsigned windowWidth, const float heritability, const float pi, const float piAlpha, const float piBeta, const bool estimatePi, const VectorXf &pis, const VectorXf &gamma, const float phi, const float kappa, const string &algorithm, const unsigned snpFittedPerWindow, const float varS, const vector<float> &S, const float overdispersion, const bool estimatePS, const float icrsq, const float spouseCorrelation, const bool diagnosticMode, const unsigned numChains, const unsigned chainLength, const unsigned burnin, const unsigned thin, const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior);
};

#endif /* amber_hpp */
