//
//  options.hpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#ifndef options_hpp
#define options_hpp
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

#include <iostream>
#include <sstream>
#include <fstream>
#include <set>
#include <stdio.h>
#include <cstring>
#include <string>
#include <limits.h>
#include <omp.h>
#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <boost/format.hpp>
#include "gadgets.hpp"

using namespace std;
using namespace boost;
using namespace Eigen;

const unsigned Megabase = 1e6;

class Options {
public:
    unsigned numChains;
    unsigned chainLength;
    unsigned burnin;
    unsigned outputFreq;
    unsigned seed;
    unsigned numThread;
    unsigned mphen; // triat order id in the phenotype file for analysis
    unsigned windowWidth; // in mega-base unit
    unsigned keepIndMax;  // the maximum number of individuals kept for analysis
    unsigned snpFittedPerWindow;    // for BayesN
    unsigned thin;  // save every this th sampled value in MCMC
    unsigned includeChr;  // chromosome to include
    unsigned ndists; // Number of distributions for base Bayes R
    unsigned flank;
    unsigned includeBlock;  // block to include
    
    float pi;
    float piAlpha;
    float piBeta;
    float heritability;
//    float varGenotypic;
//    float varResidual;
    float propVarRandom;  // proportion of variance explained by random covariate effects
    float varS; // prior variance of S in BayesS and BayesNS
    vector<float> S;    // starting value of S in BayesS and BayesNS
    float LDthreshold;  // used to define the two ends of per-SNP LD window in the banded LD matrix
    float chisqThreshold;  // significance threshold for nonzero LD chi-square test
    float piNDC;  // proportion of X-lined SNPs under no dosage compensation model (escape from X-chromosome inactivation)
    float piGxE;  // pi for genotype-by-env effects
    float phi;   // a shrinkage parameter for the heritability estimate in sbayes
    float overdispersion;
    float kappa;     // for Luke's kappa model
    float effpopNE;  // for shrunk LDM
    float genMapN;   // for shrunk LDM
    float cutOff;    // for shrunk LDM
    float icrsq;  // average inter-chromosome r^2 across SNPs
    float spouseCorrelation;
    float afDiff; // filtering SNPs by the allele frequency difference in LD and GWAS samples
    float mafmin;  // lower bound of maf
    float mafmax;  // upper bound of maf
    float lambda;  // for conjugate gradient
    float rsqThreshold;
    float pValueThreshold;
    
    bool estimatePi;
    bool estimateSigmaSq; // variance of SNP effects
    bool estimatePiNDC;  // for XCI
    bool estimatePiGxE;  // for XCI
    bool estimateScale;
    bool writeBinPosterior;
    bool writeTxtPosterior;
    bool outputResults;
    bool multiLDmat;
    bool multiThreadEigen;
    bool writeLdmTxt;      // write ldm to txt file
    bool readLdmTxt;      // read ldm from a txt file
    bool excludeMHC;  // exclude SNPs in the MHC region
    bool directPrune; // direct prune ldm
    bool estimatePS;  // estimate population stratification in sbayes
    bool diagnosticMode; // for sbayes
    bool jackknife;   // jackknife estimate for LD sampling variance
    bool excludeAmbiguousSNP;  // exlcude ambiguous SNPs with A/T or G/C alleles
    bool transpose;   // transpose the annotation file
    bool sampleOverlap;  // whether LD ref is the same as GWAS sample
    bool imputeN;  // impute per-SNP sample size
    bool noscale;
    bool simuMode; // simulation mode
    bool originalModel; // original BayesR model
    bool twoStageModel;  // two-step approach for estimating X-chr dosage model and G by sex
    bool binSnp;  // bin SNPs
    bool robustMode;  // use the robust parameterisation in SBayes models
    bool perSnpGV;
    bool mergeLdm;
    bool imputeSummary;
    
    string eigCutMethod = "value";
    float eigThreshold = 0.001;

    // Bayes R defauls
    VectorXf gamma;  // Default scaling parameters for Bayes R
    VectorXf pis;    // Default pis for Bayes R
    
    // hyperparameters for the prior distributions
    VectorXf piPar;
    Vector2f piNDCpar;
    
    // for low-rank model
    VectorXf eigenCutoff;
    
    string title;
    string analysisType;
    string bayesType;
    string algorithm;
    string optionFile;
    string phenotypeFile;
    string covariateFile;
    string randomCovariateFile;
    string bedFile;
    string alleleFreqFile;
    string includeSnpFile;
    string excludeSnpFile;
    string excludeRegionFile;
    string geneticMapFile;
    string keepIndFile;
    string snpResFile;
    string mcmcSampleFile;
    string gwasSummaryFile;
    string ldmatrixFile;
    string skeletonSnpFile;
    string annotationFile;
    string continuousAnnoFile;
    string ldscoreFile;
    string eQTLFile;
    string snpRange;
    string partParam;
    string outLDmatType;
    string windowFile;
    string residualDiagFile;
    string eigenMatrixFile;
    string ldBlockInfoFile;
    
    Options(){
        numChains               = 1;
        chainLength             = 3000;
        burnin                  = 1000;
        outputFreq              = 100;
        seed                    = 0;
        numThread               = 1;
        mphen                   = 1;
        keepIndMax              = UINT_MAX;
        snpFittedPerWindow      = 2;
        thin                    = 10;
        includeChr              = 0;
        includeBlock            = 0;
                
        windowWidth             = 0*Megabase;
        pi                      = 0.05;
        piAlpha                 = 1;
        piBeta                  = 1;
        heritability            = 0.1;
//        varGenotypic            = 1.0;
//        varResidual             = 1.0;
        propVarRandom           = 0.05;
        varS                    = 1.0;
        S.resize(1);
        S[0]                    = 0.0;
        LDthreshold             = 0.0;
        chisqThreshold          = 10;
        piNDC                   = 0.15;
        piGxE                   = 0.05;
        phi                     = 0;
        overdispersion          = 0;
        // Shrunk matrix defaults
        effpopNE                = 11490.672741;
        cutOff                  = 1e-5;
        icrsq                   = 0;
        spouseCorrelation       = 0;
        afDiff                  = 999;
        mafmin                  = 0;
        mafmax                  = 0;
        flank                   = 0;
        genMapN                 = 183; // Sample size of CEU population
        lambda                  = 1e6;
        rsqThreshold            = 1.0;
        pValueThreshold         = 1.0;

        // Bayes R defaults
        ndists                  = 5;
        gamma.resize(ndists);
        gamma                   << 0.0, 0.001, 0.01, 0.1, 1;
        pis.resize(ndists);                      
        pis                     << 0.95, 0.02, 0.01, 0.01, 0.01;
        // Kappa defaults
        kappa                   = 10;
        
        piPar.setOnes(ndists);
        piNDCpar.setOnes(2);
        
        eigenCutoff.resize(4);
        eigenCutoff             << 0.995, 0.99, 0.95, 0.9;

        estimatePi              = true;
        estimateSigmaSq         = true;
        estimatePiNDC           = true;
        estimatePiGxE           = true;
        estimateScale           = false;
        writeBinPosterior       = true;
        writeTxtPosterior       = true;
        outputResults           = true;
        multiLDmat              = false;
        multiThreadEigen        = false;
        writeLdmTxt             = false;
        readLdmTxt              = false;
        excludeMHC              = false;
        directPrune             = false;
        estimatePS              = false;
        diagnosticMode          = false;
        jackknife               = false;
        excludeAmbiguousSNP     = false;
        transpose               = false;
        sampleOverlap           = false;
        imputeN                 = false;
        noscale                 = false; // Scale the genotypes or not. Default is scaling 0
        simuMode                = false;
        originalModel           = false;
        twoStageModel           = false;
        binSnp                  = false;
        robustMode              = false;
        perSnpGV                = false;
        mergeLdm                = false;
        imputeSummary           = false;
        
        title                   = "gctb";
        analysisType            = "Bayes";
        bayesType               = "C";
        algorithm               = "";
        optionFile              = "";
        phenotypeFile           = "";
        covariateFile           = "";
        randomCovariateFile     = "";
        bedFile                 = "";
        alleleFreqFile          = "";
        includeSnpFile          = "";
        excludeSnpFile          = "";
        excludeRegionFile       = "";
        geneticMapFile          = "";
        keepIndFile             = "";
        snpResFile              = "";
        mcmcSampleFile          = "";
        gwasSummaryFile         = "";
        ldmatrixFile            = "";
        skeletonSnpFile         = "";
        annotationFile          = "";
        continuousAnnoFile      = "";
        ldscoreFile             = "";
        eQTLFile                = "";
        snpRange                = "";
        partParam               = "";
        windowFile              = "";
        residualDiagFile        = "";
        eigenMatrixFile         = "";
        ldBlockInfoFile         = "";
        outLDmatType            = "sparse";
    }
    
    void inputOptions(const int argc, const char* argv[]);
    
private:
    void readFile(const string &file);
    void makeTitle(void);
    void seedEngine(void);
    void setThread(void);
};

#endif /* options_hpp */
