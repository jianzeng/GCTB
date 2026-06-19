//
//  options.hpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#ifndef options_hpp
#define options_hpp

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
    unsigned numDist; // Number of distributions for base Bayes R
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
    float pipThreshold;
    float pepThreshold;
    float nDistAutoThreshold;
    
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
    bool hsqPercModel; // heritability percentage model: beta ~ sum pi_k N(0, gamma_k % hsq)
    bool twoStageModel;  // two-step approach for estimating X-chr dosage model and G by sex
    bool binSnp;  // bin SNPs
    bool robustMode;  // use the robust parameterisation in SBayes models
    bool perSnpGV;
    bool mergeLdm;
    bool imputeSummary;
    /// Optional SNP ID list (one ID per line; # comments): do not alter b/se/N for these SNPs after GWAS read (outlier-N filter, N imputation, LD imputation).
    string keepSumstatIntactFile;
    bool nDistAuto;  // automatically determine the number of mixture distributions
    bool nDistAutoByPred;  // if true with nDistAuto, select K by pseudo CV prediction accuracy (else heritability)
    bool writeWandQ; // output w and Q in text format
    bool estimateRsqEnrich;  // estimate prediction R2 enrichment
    GwasScalarMode gwasScalarMode;  // --gwas-scalar: effect scaling for summary-data GWAS (see Data::gwasScalarMode)
    int eigenMatrixQuantBits;    // 0 for float .eigen.bin; 4/8/16 for U-column quantized .eigen.q*.bin
    bool eigenMatrixQ8Entropy;   // true for .eigen.q8e.bin
    int quantEigenBits;
    bool quantEigenEntropy;
    string quantEigenInputDir;
    string quantEigenOutputDir;
    
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
    string optionsSummary;
    string analysisType;
    string bayesType;
    string algorithm;
    string optionFile;
    string phenotypeFile;
    string covariateFile;
    string randomCovariateFile;
    string bedFile;
    string genoTxtFile;
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
    string plinkLDtxtfile;
    string plinkLDbinfile;
    string plinkAFfile;
    string pairwiseLDfile;
    string ldfriendFile;
    string geneMapFile;
    string genomeBuild;
    string skipSnpFile;
    /// If true with --skip, set GWAS marginal effects (b, b2) to zero for those SNPs after read/impute
    bool setZeroGwasForSkip;
    /// If true with --skip, re-eigen-decompose the LD submatrix of non-skipped SNPs (needs block*.ldm.bin or block*.eigen.bin per block)
    bool recomputeEigen;
    /// Block LD folder with block*.ldm.bin + ldm/snp info: full dense LD for SBayes (mutually exclusive with --ldm-eigen for the same role)
    string ldmBlockDir;
    /// Set via --fixed-effect; one SNP ID per line (first column if tab/space separated); lines starting with # ignored
    string fitSnpsAsFixedEffectsFile;
    string label;
    
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
        pi                      = 0.01;
        piAlpha                 = 1;
        piBeta                  = 1;
        heritability            = 0.5;
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
        pipThreshold            = 0.9;
        pepThreshold            = 0.7;
        nDistAutoThreshold      = 0.5;

        // Bayes R defaults
        numDist                  = 5;
        gamma.resize(numDist);
        gamma                   << 0.0, 0.001, 0.01, 0.1, 1;
        pis.resize(numDist);
        pis                     << 0.99, 0.005, 0.003, 0.001, 0.001;
        // Kappa defaults
        kappa                   = 10;
        
        piPar.setOnes(numDist);
        piNDCpar.setOnes(2);
        
        eigenCutoff.resize(4);
        eigenCutoff             << 0.995, 0.99, 0.95, 0.9;

        estimatePi              = true;
        estimateSigmaSq         = true;
        estimatePiNDC           = true;
        estimatePiGxE           = true;
        estimateScale           = false;
        writeBinPosterior       = false;
        writeTxtPosterior       = false;
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
        hsqPercModel            = true;
        twoStageModel           = false;
        binSnp                  = false;
        robustMode              = false;
        perSnpGV                = false;
        mergeLdm                = false;
        imputeSummary           = false;
        keepSumstatIntactFile   = "";
        nDistAuto               = false;
        nDistAutoByPred         = false;
        writeWandQ              = false;
        estimateRsqEnrich       = false;
        gwasScalarMode          = GwasScalarMode::LeastSquares;
        eigenMatrixQuantBits    = 0;
        eigenMatrixQ8Entropy    = false;
        quantEigenBits          = 8;
        quantEigenEntropy       = false;
        quantEigenInputDir      = "";
        quantEigenOutputDir     = "";
        
        title                   = "gctb";
        optionsSummary          = "";
        analysisType            = "Bayes";
        bayesType               = "C";
        algorithm               = "";
        optionFile              = "";
        phenotypeFile           = "";
        covariateFile           = "";
        randomCovariateFile     = "";
        bedFile                 = "";
        genoTxtFile             = "";
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
        ldmBlockDir             = "";
        ldBlockInfoFile         = "";
        plinkLDtxtfile          = "";
        plinkLDbinfile          = "";
        plinkAFfile             = "";
        pairwiseLDfile          = "";
        ldfriendFile            = "";
        geneMapFile             = "";
        outLDmatType            = "sparse";
        genomeBuild             = "hg19";
        skipSnpFile             = "";
        setZeroGwasForSkip      = false;
        recomputeEigen          = false;
        fitSnpsAsFixedEffectsFile = "";
        label                   = "";
    }
    
    void inputOptions(const int argc, const char* argv[]);
    const string& getOptionsSummary(void) const { return optionsSummary; }
    void setThread(void);
    
private:
    void readFile(const string &file);
    void makeTitle(void);
    void seedEngine(void);
};

#endif /* options_hpp */
