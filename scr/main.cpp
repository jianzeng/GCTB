//
//  main.cpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include <iostream>
#include "gctb.hpp"
#include "xci.hpp"
#include "vgmaf.hpp"

using namespace std;


int main(int argc, const char * argv[]) {
    
    cout << "******************************************************************\n";
    cout << "* GCTB 2.05beta                                                  *\n";
    cout << "* Genome-wide Complex Trait Bayesian analysis                    *\n";
    cout << "* Authors: Jian Zeng, Luke Lloyd-Jones, Zhili Zheng, Shouye Liu  *\n";
    cout << "* MIT License                                                    *\n";
    cout << "******************************************************************\n";
    
    Gadget::Timer timer;
    timer.setTime();
    cout << "\nAnalysis started: " << timer.getDate();
    
    if (argc < 2){
        cerr << " \nDid you forget to give the input parameters?\n" << endl;
        exit(1);
    }
    
    try {
        
        Options opt;
        opt.inputOptions(argc, argv);
        
        if (opt.seed) Stat::seedEngine(opt.seed);
        else          Stat::seedEngine(011415);  // fix the random seed if not given due to the use of MPI
        
//        cout << "==========" << opt.seed << " " << Stat::ranf() << " " << Stat::snorm() << endl;
        
        Data data;
        data.title = opt.title;
        bool readGenotypes;
        
        GCTB gctb(opt);


        if (opt.analysisType == "Bayes") {
            if (opt.numChains > 1) {
                throw(" Error: multi-chain MCMC is not yet available for individual-level-data analysis.");
            }
            readGenotypes = false;
            gctb.inputIndInfo(data, opt.bedFile, opt.phenotypeFile, opt.keepIndFile, opt.keepIndMax,
                               opt.mphen, opt.covariateFile, opt.randomCovariateFile, opt.residualDiagFile);
            gctb.inputSnpInfo(data, opt.bedFile, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.ldBlockInfoFile, opt.includeBlock, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.mafmin, opt.mafmax, opt.noscale, readGenotypes);
            
            Model *model = gctb.buildModel(data, opt.bedFile, "", opt.bayesType, opt.windowWidth,
                                            opt.heritability, opt.propVarRandom, opt.pi, opt.piAlpha, opt.piBeta, opt.estimatePi, opt.noscale, opt.pis, opt.piPar, opt.gamma, opt.estimateSigmaSq, opt.phi, opt.kappa,
                                            opt.algorithm, opt.snpFittedPerWindow, opt.varS, opt.S, opt.overdispersion, opt.estimatePS, opt.icrsq, opt.spouseCorrelation, opt.diagnosticMode, opt.originalModel, opt.perSnpGV, opt.robustMode);
            vector<McmcSamples*> mcmcSampleVec = gctb.runMcmc(*model, opt.chainLength, opt.burnin, opt.thin,
                                                               opt.outputFreq, opt.title, opt.writeBinPosterior, opt.writeTxtPosterior);
            //gctb.saveMcmcSamples(mcmcSampleVec, opt.title);
            gctb.clearGenotypes(data);
            if (opt.outputResults) gctb.outputResults(data, mcmcSampleVec, opt.bayesType, opt.noscale, opt.title);
        }
        else if (opt.analysisType == "LDmatrix") {
            readGenotypes = false;
            if (opt.ldmatrixFile.empty()) { // make LD matrix from genotypes
                gctb.inputIndInfo(data, opt.bedFile, opt.bedFile + ".fam", opt.keepIndFile, opt.keepIndMax,
                                  opt.mphen, opt.covariateFile, opt.randomCovariateFile, opt.residualDiagFile);
                gctb.inputSnpInfo(data, opt.bedFile, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.ldBlockInfoFile, opt.includeBlock, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.mafmin, opt.mafmax, opt.noscale, readGenotypes);
                if (opt.outLDmatType == "shrunk") {
                    data.makeshrunkLDmatrix(opt.bedFile + ".bed", opt.outLDmatType, opt.snpRange, opt.title, opt.writeLdmTxt, opt.effpopNE, opt.cutOff, opt.genMapN);
                } else if (opt.outLDmatType == "block") {
                    data.makeBlockLDmatrix(opt.bedFile + ".bed", opt.outLDmatType, opt.includeBlock, opt.title, opt.writeLdmTxt);
                }
                else {
                    string snpRange = opt.snpRange;
                    if(!opt.partParam.empty()){
                        snpRange = data.partLDMatrix(opt.partParam, opt.title, opt.outLDmatType);
                    }
                    data.makeLDmatrix(opt.bedFile + ".bed", opt.outLDmatType, opt.chisqThreshold, opt.LDthreshold, opt.windowWidth, snpRange, opt.title, opt.writeLdmTxt);
                }
            }
//            else if (opt.ldmatrixFile.empty() != 1 && opt.outLDmatType == "shrunk" || opt.outLDmatType == "sparseshrunk") { // make shrunk LD matrix from other LDM
//                gctb.inputSnpInfo(data, opt.includeSnpFile, opt.excludeSnpFile, "", opt.ldmatrixFile, opt.includeChr, opt.multiLDmat, opt.geneticMapFile);
//                data.resizeLDmatrix(opt.outLDmatType, opt.chisqThreshold, opt.windowWidth, opt.LDthreshold, opt.effpopNE, opt.cutOff, opt.genMapN);
//                data.outputLDmatrix(opt.outLDmatType, opt.title);
//            }
//            else if (opt.ldmatrixFile.empty() == 1 && opt.outLDmatType != "shrunk") { // make LD matrix from genotypes
//                gctb.inputIndInfo(data, opt.bedFile, opt.bedFile + ".fam", opt.keepIndFile, opt.keepIndMax,
//                                  opt.mphen, opt.covariateFile);
//                gctb.inputSnpInfo(data, opt.bedFile, opt.includeSnpFile, opt.excludeSnpFile, opt.includeChr, readGenotypes);
//                data.makeLDmatrix(opt.bedFile + ".bed", opt.outLDmatType, opt.chisqThreshold, opt.LDthreshold, opt.windowWidth, opt.snpRange, opt.title);
//            }
            else { // manipulate an existing LD matrix or merge existing LD matrices
                if (opt.mergeLdm) {
                    data.mergeLdmInfo(opt.outLDmatType, opt.ldmatrixFile);
                }
                else if (opt.directPrune) {
                    data.directPruneLDmatrix(opt.ldmatrixFile, opt.outLDmatType, opt.chisqThreshold, opt.title, opt.writeLdmTxt);
                }
                else if (opt.jackknife) {
                    readGenotypes = true;
                    gctb.inputIndInfo(data, opt.bedFile, opt.bedFile + ".fam", opt.keepIndFile, opt.keepIndMax, opt.mphen, opt.covariateFile, opt.randomCovariateFile, opt.residualDiagFile);
                    gctb.inputSnpInfo(data, opt.bedFile, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.ldBlockInfoFile, opt.includeBlock, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.mafmin, opt.mafmax, opt.noscale, readGenotypes);
                    data.jackknifeLDmatrix(opt.ldmatrixFile, opt.outLDmatType, opt.title, opt.writeLdmTxt);
                }
                else if (opt.binSnp) {
                    gctb.inputSnpInfo(data, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, "", opt.ldmatrixFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.genMapN, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.ldscoreFile, opt.windowFile, opt.multiLDmat, opt.excludeMHC, opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.rsqThreshold, opt.sampleOverlap, opt.imputeN, opt.noscale, opt.binSnp, opt.readLdmTxt);
                    data.binSnpByLDrsq(opt.rsqThreshold, opt.title);
                }
                else {
                    gctb.inputSnpInfo(data, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, "", opt.ldmatrixFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.genMapN, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.ldscoreFile, opt.windowFile, opt.multiLDmat, opt.excludeMHC, opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.rsqThreshold, opt.sampleOverlap, opt.imputeN, opt.noscale, opt.binSnp, opt.readLdmTxt);
                    data.resizeLDmatrix(opt.outLDmatType, opt.chisqThreshold, opt.windowWidth, opt.LDthreshold, opt.effpopNE, opt.cutOff, opt.genMapN);
                    data.outputLDmatrix(opt.outLDmatType, opt.title, opt.writeLdmTxt);
                }
            }
        }
        else if (opt.analysisType == "LDmatrixEigen") {
            readGenotypes = false;
            if (opt.eigenMatrixFile.empty()) { // perform eigen decomposition for the blocked LD matrices
                //gctb.inputSnpInfo(data, opt.bedFile, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.mafmin, opt.mafmax, opt.noscale, readGenotypes);
                //data.getEigenDataFromFullLDM(opt.title, opt.eigenCutoff);
                data.readBlockLDmatrixAndDoEigenDecomposition(opt.ldmatrixFile, opt.includeBlock, opt.eigenCutoff.maxCoeff(), opt.writeLdmTxt);
            }
            else { // merge existing eigen matrices
                //gctb.inputSnpInfo(data, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, "", opt.eigenMatrixFile, opt.ldBlockInfoFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.ldscoreFile, opt.eigenCutoff, opt.excludeMHC, opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.rsqThreshold, opt.sampleOverlap, opt.imputeN, opt.noscale, opt.readLdmTxt);
            }
        }
        else if (opt.analysisType == "ImputeSumStats") {
            readGenotypes = false;
            if (opt.eigenMatrixFile.empty()) {
                throw("Error: --impute-summary requires the results of eigen-decomposition of LD matrices as input.");
            }
            gctb.inputSnpInfo(data, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile,
                              opt.gwasSummaryFile, opt.eigenMatrixFile, opt.ldBlockInfoFile,
                              opt.includeChr, opt.excludeAmbiguousSNP,
                              opt.annotationFile, opt.transpose,
                              opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.ldscoreFile,
                              opt.eigenCutoff.maxCoeff(), opt.excludeMHC,
                              opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.rsqThreshold,
                              opt.sampleOverlap, opt.imputeN, opt.noscale, opt.readLdmTxt, opt.imputeSummary, opt.includeBlock);
        }
        else if (opt.analysisType == "MergeGwasSummary") {
            if (opt.outLDmatType == "block") {
                data.mergeBlockGwasSummary(opt.gwasSummaryFile, opt.title);
            }
        }
        else if (opt.analysisType == "SBayes") {
            if (!opt.ldmatrixFile.empty()) {
                gctb.inputSnpInfo(data, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, opt.gwasSummaryFile, opt.ldmatrixFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.genMapN, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.ldscoreFile, opt.windowFile, opt.multiLDmat, opt.excludeMHC, opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.rsqThreshold, opt.sampleOverlap, opt.imputeN, opt.noscale, opt.binSnp, opt.readLdmTxt);
            } else if (!opt.eigenMatrixFile.empty()) {  // low-rank model
                gctb.inputSnpInfo(data, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile,
                                  opt.gwasSummaryFile, opt.eigenMatrixFile, opt.ldBlockInfoFile,
                                  opt.includeChr, opt.excludeAmbiguousSNP,
                                  opt.annotationFile, opt.transpose,
                                  opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.ldscoreFile,
                                  opt.eigenCutoff.maxCoeff(), opt.excludeMHC,
                                  opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.rsqThreshold,
                                  opt.sampleOverlap, opt.imputeN, opt.noscale, opt.readLdmTxt, opt.imputeSummary, opt.includeBlock);
                float bestEigenCutoff = gctb.tuneEigenCutoff(data, opt);
                data.readEigenMatrixBinaryFileAndMakeWandQ(opt.eigenMatrixFile, bestEigenCutoff, data.gwasEffectInBlock, data.numKeptInds, true);
                //data.readEigenMatrixBinaryFile(opt.eigenMatrixFile, bestEigenCutoff);
                //data.constructWandQ(data.gwasEffectInBlock, data.numKeptInds);
            } else {
                gctb.inputSnpInfo(data, opt.bedFile, opt.gwasSummaryFile, opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.sampleOverlap, opt.imputeN, opt.noscale);
            }
            
            data.label = opt.title;
            if (opt.numChains > 1) {
                vector<McmcSamples*> mcmcSampleVec = gctb.multi_chain_mcmc(data, opt.bayesType, opt.windowWidth, opt.heritability, opt.propVarRandom, opt.pi, opt.piAlpha, opt.piBeta, opt.estimatePi, opt.pis, opt.gamma, opt.phi, opt.kappa, opt.algorithm, opt.snpFittedPerWindow, opt.varS, opt.S, opt.overdispersion, opt.estimatePS, opt.icrsq, opt.spouseCorrelation, opt.diagnosticMode, opt.robustMode, opt.numChains, opt.chainLength, opt.burnin, opt.thin, opt.outputFreq, opt.title, opt.writeBinPosterior, opt.writeTxtPosterior);
                if (opt.outputResults) gctb.outputResults(data, mcmcSampleVec, opt.bayesType, opt.noscale, opt.title);
            } else {
                Model *model = gctb.buildModel(data, opt.bedFile, opt.gwasSummaryFile, opt.bayesType, opt.windowWidth,
                                               opt.heritability, opt.propVarRandom, opt.pi, opt.piAlpha, opt.piBeta, opt.estimatePi, opt.noscale, opt.pis, opt.piPar, opt.gamma, opt.estimateSigmaSq, opt.phi, opt.kappa,
                                               opt.algorithm, opt.snpFittedPerWindow, opt.varS, opt.S, opt.overdispersion, opt.estimatePS, opt.icrsq, opt.spouseCorrelation, opt.diagnosticMode, opt.originalModel, opt.perSnpGV, opt.robustMode);
                vector<McmcSamples*> mcmcSampleVec;
                
                try{
                    mcmcSampleVec = gctb.runMcmc(*model, opt.chainLength, opt.burnin, opt.thin,
                                                                  opt.outputFreq, opt.title, opt.writeBinPosterior, opt.writeTxtPosterior);
                }
                catch(const string &err_msg) {
                    cout << err_msg << endl;
                    if (opt.robustMode) {
                        cout << "Please refer to our website (https://cnsgenomics.com/software/gctb) for further information on this problem." << endl;
                        exit(1);
                    } else {
                        cout << "\nRestarting MCMC with a more robust parameterisation for SBayes" << opt.bayesType << " ..." << endl;
                        cout << "Please refer to our website (https://cnsgenomics.com/software/gctb) for more information." << endl;
                       opt.robustMode = true;
                        Model *model = gctb.buildModel(data, opt.bedFile, opt.gwasSummaryFile, opt.bayesType, opt.windowWidth,
                                                       opt.heritability, opt.propVarRandom, opt.pi, opt.piAlpha, opt.piBeta, opt.estimatePi, opt.noscale, opt.pis, opt.piPar, opt.gamma, opt.estimateSigmaSq, opt.phi, opt.kappa,
                                                       opt.algorithm, opt.snpFittedPerWindow, opt.varS, opt.S, opt.overdispersion, opt.estimatePS, opt.icrsq, opt.spouseCorrelation, opt.diagnosticMode, opt.originalModel, opt.perSnpGV, opt.robustMode);

                        mcmcSampleVec = gctb.runMcmc(*model, opt.chainLength, opt.burnin, opt.thin,
                                                     opt.outputFreq, opt.title, opt.writeBinPosterior, opt.writeTxtPosterior);
                    }
                }
                catch (const char *err_msg) {
                    cout << err_msg << endl;
                    if (opt.robustMode) {
                        cout << "Please refer to our website (https://cnsgenomics.com/software/gctb) for further information on this problem." << endl;
                        exit(1);
                    } else {
                        cout << "\nRestarting MCMC with a more robust parameterisation for SBayes" << opt.bayesType << " ..." << endl;
                        cout << "Please refer to our website (https://cnsgenomics.com/software/gctb) for more information." << endl;
                        opt.robustMode = true;
                        Model *model = gctb.buildModel(data, opt.bedFile, opt.gwasSummaryFile, opt.bayesType, opt.windowWidth,
                                                       opt.heritability, opt.propVarRandom, opt.pi, opt.piAlpha, opt.piBeta, opt.estimatePi, opt.noscale, opt.pis, opt.piPar, opt.gamma, opt.estimateSigmaSq, opt.phi, opt.kappa,
                                                       opt.algorithm, opt.snpFittedPerWindow, opt.varS, opt.S, opt.overdispersion, opt.estimatePS, opt.icrsq, opt.spouseCorrelation, opt.diagnosticMode, opt.originalModel, opt.perSnpGV, opt.robustMode);

                        mcmcSampleVec = gctb.runMcmc(*model, opt.chainLength, opt.burnin, opt.thin,
                                                     opt.outputFreq, opt.title, opt.writeBinPosterior, opt.writeTxtPosterior);
                    }
                }

                if (opt.outputResults) gctb.outputResults(data, mcmcSampleVec, opt.bayesType, opt.noscale, opt.title);
            }
        }
        else if (opt.analysisType == "ConjugateGradient") {
            gctb.inputSnpInfo(data, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, opt.gwasSummaryFile, opt.ldmatrixFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.genMapN, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.ldscoreFile, opt.windowFile, opt.multiLDmat, opt.excludeMHC, opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.rsqThreshold, opt.sampleOverlap, opt.imputeN, opt.noscale, opt.binSnp, opt.readLdmTxt);
            gctb.solveSnpEffectsByConjugateGradientMethod(data, opt.lambda, opt.title + ".snpRes");
        }
        else if (opt.analysisType == "Stratify") { // post hoc stratified analysis
            gctb.stratify(data, opt.ldmatrixFile, opt.multiLDmat, opt.geneticMapFile, opt.genMapN, opt.snpResFile, opt.mcmcSampleFile, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.gwasSummaryFile, opt.pValueThreshold, opt.imputeN, opt.title, opt.bayesType, opt.chainLength, opt.burnin, opt.thin, opt.outputFreq);
        }
        else if (opt.analysisType == "hsq") {
            if (opt.ldmatrixFile.empty()) {
                readGenotypes = true;
                gctb.inputIndInfo(data, opt.bedFile, opt.phenotypeFile, opt.keepIndFile, opt.keepIndMax,
                               opt.mphen, opt.covariateFile, opt.randomCovariateFile, opt.residualDiagFile);
                gctb.inputSnpInfo(data, opt.bedFile, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.ldBlockInfoFile, opt.includeBlock, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.mafmin, opt.mafmax, opt.noscale, readGenotypes);
            } else {
                gctb.inputSnpInfo(data, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, opt.gwasSummaryFile, opt.ldmatrixFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.genMapN, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.ldscoreFile, opt.windowFile, opt.multiLDmat, opt.excludeMHC, opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.rsqThreshold, opt.sampleOverlap, opt.imputeN, opt.noscale, opt.binSnp, opt.readLdmTxt);
                if (data.sparseLDM) data.getZPZspmat();
                else data.getZPZmat();
            }
            McmcSamples *snpEffects = gctb.inputMcmcSamples(opt.mcmcSampleFile, "SnpEffects", "bin");
            McmcSamples *resVar = gctb.inputMcmcSamples(opt.mcmcSampleFile, "ResVar", "txt");
            gctb.estimateHsq(data, *snpEffects, *resVar, opt.title, opt.outputFreq);
        }
        else if (opt.analysisType == "Pi") {
            if (opt.ldmatrixFile.empty()) {
                readGenotypes = true;  // need this to calculate allele frequencies
                gctb.inputIndInfo(data, opt.bedFile, opt.phenotypeFile, opt.keepIndFile, opt.keepIndMax,
                               opt.mphen, opt.covariateFile, opt.randomCovariateFile, opt.residualDiagFile);
                gctb.inputSnpInfo(data, opt.bedFile, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.ldBlockInfoFile, opt.includeBlock, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.mafmin, opt.mafmax, opt.noscale, readGenotypes);
            } else {
                gctb.inputSnpInfo(data, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, opt.gwasSummaryFile, opt.ldmatrixFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.genMapN, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.ldscoreFile, opt.windowFile, opt.multiLDmat, opt.excludeMHC, opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.rsqThreshold, opt.sampleOverlap, opt.imputeN, opt.noscale, opt.binSnp, opt.readLdmTxt);
                //if (data.sparseLDM) data.getZPZspmat();
                //else data.getZPZmat();
            }
            McmcSamples *snpEffects = gctb.inputMcmcSamples(opt.mcmcSampleFile, "SnpEffects", "bin");
            McmcSamples *genVar = gctb.inputMcmcSamples(opt.mcmcSampleFile, "GenVar", "txt");
            gctb.estimatePi(data, *snpEffects, *genVar, opt.title, opt.outputFreq);
        }
        else if (opt.analysisType == "Predict") {
            readGenotypes = true;
            gctb.inputIndInfo(data, opt.bedFile, opt.phenotypeFile, opt.keepIndFile, opt.keepIndMax,
                               opt.mphen, opt.covariateFile, opt.randomCovariateFile, opt.residualDiagFile);
            gctb.inputSnpInfo(data, opt.bedFile, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.ldBlockInfoFile, opt.includeBlock, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.mafmin, opt.mafmax, opt.noscale, readGenotypes);
            
            data.inputSnpResults(opt.snpResFile);
            gctb.predict(data, opt.title);
        }
        else if (opt.analysisType == "Summarize") {  // ad hoc method for producing summary from binary MCMC samples of SNP effects
            readGenotypes = true;
            gctb.inputIndInfo(data, opt.bedFile, opt.phenotypeFile, opt.keepIndFile, opt.keepIndMax,
                               opt.mphen, opt.covariateFile, opt.randomCovariateFile, opt.residualDiagFile);
            gctb.inputSnpInfo(data, opt.bedFile, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.ldBlockInfoFile, opt.includeBlock, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.mafmin, opt.mafmax, opt.noscale, readGenotypes);
            gctb.clearGenotypes(data);
            McmcSamples *snpEffects = gctb.inputMcmcSamples(opt.mcmcSampleFile, "SnpEffects", "bin");
            data.summarizeSnpResults(snpEffects->datMatSp, opt.title + ".snpRes");
        }
        else if (opt.analysisType == "XCI") {  // ad hoc method for X chromosome inactivation project
            XCI xci;
            readGenotypes = true;
            xci.inputIndInfo(data, opt.bedFile, opt.phenotypeFile, opt.keepIndFile, opt.keepIndMax,
                             opt.mphen, opt.covariateFile);
            xci.inputSnpInfo(data, opt.bedFile, opt.includeSnpFile, opt.excludeSnpFile, opt.includeChr, opt.annotationFile, opt.windowFile, readGenotypes);
            if (opt.simuMode) {
                xci.simu(data, opt.pi, opt.heritability, opt.piNDC, opt.piGxE, false, opt.title, opt.seed);  // ad hoc simulation to test BayesXCI method
            }
            else {
                if (opt.twoStageModel) {
                    // Stage 1: estimate NDC using female data only
                    Model *model = xci.buildModelStageOne(data, "C", opt.heritability, opt.pi, opt.piPar, opt.estimatePi, opt.piNDC, opt.piNDCpar, opt.estimatePiNDC);
                    vector<McmcSamples*> mcmcSampleVec = gctb.runMcmc(*model, opt.chainLength, opt.burnin, opt.thin,
                                                                      opt.outputFreq, opt.title + ".stage1", opt.writeBinPosterior, opt.writeTxtPosterior);
                    gctb.saveMcmcSamples(mcmcSampleVec, opt.title + ".stage1");
                    gctb.outputResults(data, mcmcSampleVec, "C", true, opt.title + ".stage1");
                    xci.outputResults(data, mcmcSampleVec, "C", opt.title + ".stage1");
                    delete model;
                    // Stage 2: estimate GxS using both male and female data
                    model = xci.buildModelStageTwo(data, "Cgxs", opt.heritability, opt.pi, opt.piPar, opt.estimatePi, opt.piNDC, opt.piNDCpar, opt.estimatePiNDC, opt.title + ".stage1.snpRes", opt.piGxE, opt.estimatePiGxE);
                    mcmcSampleVec = gctb.runMcmc(*model, opt.chainLength, opt.burnin, opt.thin,
                                                                      opt.outputFreq, opt.title + ".stage2", opt.writeBinPosterior, opt.writeTxtPosterior);
                    gctb.saveMcmcSamples(mcmcSampleVec, opt.title + ".stage2");
                    gctb.clearGenotypes(data);
                    gctb.outputResults(data, mcmcSampleVec, "Cgxs", true, opt.title + ".stage2");
                    xci.outputResults(data, mcmcSampleVec, "Cgxs", opt.title + ".stage2");
                    
                }
                else {
                    if (opt.numChains > 1) {  // multi chains
                        vector<McmcSamples*> mcmcSampleVec = xci.multi_chain_mcmc(data, opt.bayesType, opt.heritability, opt.pi, opt.piPar, opt.estimatePi, opt.piNDC, opt.piNDCpar, opt.estimatePiNDC, opt.piGxE, opt.estimatePiGxE, opt.numChains, opt.chainLength, opt.burnin, opt.thin, opt.outputFreq, opt.title, opt.writeBinPosterior, opt.writeTxtPosterior);
                        gctb.saveMcmcSamples(mcmcSampleVec, opt.title);
                        gctb.clearGenotypes(data);
                        gctb.outputResults(data, mcmcSampleVec, opt.bayesType, true, opt.title);
                        xci.outputResults(data, mcmcSampleVec, opt.bayesType, opt.title);
                    } else {
                        Model *model = xci.buildModel(data, opt.bayesType, opt.heritability, opt.pi, opt.piPar, opt.estimatePi, opt.piNDC, opt.piNDCpar, opt.estimatePiNDC, opt.piGxE, opt.estimatePiGxE, opt.windowWidth);
                        vector<McmcSamples*> mcmcSampleVec = gctb.runMcmc(*model, opt.chainLength, opt.burnin, opt.thin,
                                                                          opt.outputFreq, opt.title, opt.writeBinPosterior, opt.writeTxtPosterior);
                        gctb.saveMcmcSamples(mcmcSampleVec, opt.title);
                        gctb.clearGenotypes(data);
                        gctb.outputResults(data, mcmcSampleVec, opt.bayesType, true, opt.title);
                        xci.outputResults(data, mcmcSampleVec, opt.bayesType, opt.title);
                    }
                }
            }
        }
        else if (opt.analysisType == "VGMAF") {  // ad hoc method for cumulative Vg against MAF to detect selection
            readGenotypes = true;
            gctb.inputIndInfo(data, opt.bedFile, opt.phenotypeFile, opt.keepIndFile, opt.keepIndMax,
                               opt.mphen, opt.covariateFile, opt.randomCovariateFile, opt.residualDiagFile);
            gctb.inputSnpInfo(data, opt.bedFile, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.ldBlockInfoFile, opt.includeBlock, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.mafmin, opt.mafmax, opt.noscale, readGenotypes);
            VGMAF vgmaf;
            if (opt.bayesType == "Simu") {
                vgmaf.simulate(data, opt.title);
            } else {
                McmcSamples *snpEffects = gctb.inputMcmcSamples(opt.mcmcSampleFile, "SnpEffects", "bin");
                vgmaf.compute(data, *snpEffects, opt.burnin, opt.thin, opt.title);
            }
        }
        
        else if (opt.analysisType == "OutputEffectSamples") { // for now an ad hoc method to output the MCMC SNP effect samples in text file
            McmcSamples *snpEffects = gctb.inputMcmcSamples(opt.mcmcSampleFile, "SnpEffects", "bin");
            data.outputSnpEffectSamples(snpEffects->datMatSp, opt.burnin, opt.outputFreq, opt.snpResFile, opt.title + ".snpEffectSamples");
        }
        
        else {
            throw(" Error: Wrong analysis type: " + opt.analysisType);
        }
    }
    catch (const string &err_msg) {
        cerr << "\n" << err_msg << endl;
    }
    catch (const char *err_msg) {
        cerr << "\n" << err_msg << endl;
    }
    
    timer.getTime();
    
    cout << "\nAnalysis finished: " << timer.getDate();
    cout << "Computational time: "  << timer.format(timer.getElapse()) << endl;

    return 0;
}
