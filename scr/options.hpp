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
#include <mpi.h>
#include <omp.h>
#include <boost/format.hpp>
#include <Eigen/Core>
#include <Eigen/Eigen>
#include "mympi.hpp"
#include "gadgets.hpp"

using namespace std;
using namespace boost;
using namespace Eigen;

const unsigned Megabase = 1e6;

class Options {
public:
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
    
    float pi;
    float piAlpha;
    float piBeta;
    float heritability;
//    float varGenotypic;
//    float varResidual;
    float varS; // prior variance of S in BayesS and BayesNS
    vector<float> S;    // starting value of S in BayesS and BayesNS
    float LDthreshold;  // used to define the two ends of per-SNP LD window in the banded LD matrix
    float chisqThreshold;  // significance threshold for nonzero LD chi-square test
    float piNDC;  // proportion of X-lined SNPs under no dosage compensation model (escape from X-chromosome inactivation)
    float phi;   // a shrinkage parameter for the heritability estimate in sbayes
    float overdispersion;
    float kappa;     // for Luke's kappa model
    float effpopNE;  // for shrunk LDM
    float cutOff;    // for shrunk LDM
    float icrsq;  // average inter-chromosome r^2 across SNPs
    float spouseCorrelation;
    float afDiff; // filtering SNPs by the allele frequency difference in LD and GWAS samples
    float mafmin;  // lower bound of maf
    float mafmax;  // upper bound of maf
    
    bool estimatePi;
    bool estimateScale;
    bool writeBinPosterior;
    bool outputResults;
    bool multiLDmat;
    bool multiThreadEigen;
    bool writeLdmTxt;      // write ldm to txt file
    bool excludeMHC;  // exclude SNPs in the MHC region
    bool directPrune; // direct prune ldm
    bool estimatePS;  // estimate population stratification in sbayes
    bool diagnosticMode; // for sbayes
    bool jackknife;   // jackknife estimate for LD sampling variance
    bool excludeAmbiguousSNP;  // exlcude ambiguous SNPs with A/T or G/C alleles

    // Bayes R defauls
    VectorXf gamma;  // Default scaling parameters for Bayes R
    VectorXf pis;    // Default pis for Bayes R
    
    string title;
    string analysisType;
    string bayesType;
    string algorithm;
    string optionFile;
    string phenotypeFile;
    string covariateFile;
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
    string snpRange;
    string outLDmatType;
    
    Options(){
        chainLength             = 50000;
        burnin                  = 20000;
        outputFreq              = 100;
        seed                    = 0;
        numThread               = 1;
        mphen                   = 1;
        keepIndMax              = UINT_MAX;
        snpFittedPerWindow      = 2;
        thin                    = 10;
        includeChr              = 0;
                
        windowWidth             = 0*Megabase;
        pi                      = 0.05;
        piAlpha                 = 1;
        piBeta                  = 19;
        heritability            = 0.1;
//        varGenotypic            = 1.0;
//        varResidual             = 1.0;
        varS                    = 1.0;
        S.resize(1);
        S[0]                    = 0.0;
        LDthreshold             = 0.0;
        chisqThreshold          = 10;
        piNDC                   = 0.15;
        phi                     = 0;
        overdispersion          = 0;
        // Shrunk matrix defaults
        effpopNE                = 11490.672741;
        cutOff                  = 1e-5;
        icrsq                   = 0;
        spouseCorrelation       = 0;
        afDiff                  = 0.05;
        mafmin                  = 0;
        mafmax                  = 0;

        // Bayes R defaults
        ndists                  = 4;
        gamma.resize(ndists);
        gamma                   << 0.0, 0.01, 0.1, 1;  
        pis.resize(ndists);                      
        pis                     << 0.95, 0.03, 0.01, 0.01;
        // Kappa defaults
        kappa                   = 10;

        estimatePi              = true;
        estimateScale           = false;
        writeBinPosterior       = true;
        outputResults           = true;
        multiLDmat              = false;
        multiThreadEigen        = false;
        writeLdmTxt             = false;
        excludeMHC              = false;
        directPrune             = false;
        estimatePS              = false;
        diagnosticMode          = false;
        jackknife               = false;
        excludeAmbiguousSNP     = false;
        
        title                   = "gctb";
        analysisType            = "Bayes";
        bayesType               = "C";
        algorithm               = "";
        optionFile              = "";
        phenotypeFile           = "";
        covariateFile           = "";
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
        snpRange                = "";
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
