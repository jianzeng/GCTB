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
    
    cout << "*********************************************************\n";
    cout << "* GCTB 2.5.5.5                                          *\n";
    cout << "* Genome-wide Complex Trait Bayesian analysis           *\n";
    cout << "* For inquiries, contact: Jian Zeng <j.zeng@uq.edu.au>  *\n";
    cout << "* Last updated: 12 Dec, 2025                            *\n";
    cout << "* MIT License                                           *\n";
    cout << "*********************************************************\n";
    
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
        data.gwasScalarMode = opt.gwasScalarMode;
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
            
            Model *model = gctb.buildModel(data, opt, opt.bedFile, "", opt.bayesType, opt.windowWidth,
                                            opt.heritability, opt.propVarRandom, opt.pi, opt.piAlpha, opt.piBeta, opt.estimatePi, opt.noscale, opt.pis, opt.piPar, opt.gamma, opt.estimateSigmaSq, opt.phi, opt.kappa,
                                            opt.algorithm, opt.snpFittedPerWindow, opt.varS, opt.S, opt.overdispersion, opt.estimatePS, opt.icrsq, opt.spouseCorrelation, opt.diagnosticMode, opt.hsqPercModel, opt.perSnpGV, opt.robustMode, opt.nDistAuto, opt.estimateRsqEnrich);
            vector<McmcSamples*> mcmcSampleVec = gctb.runMcmc(*model, 1, opt.chainLength, opt.burnin, opt.thin,
                                                               opt.outputFreq, opt.title, opt.writeBinPosterior, opt.writeTxtPosterior);
            //gctb.saveMcmcSamples(mcmcSampleVec, opt.title);
            gctb.clearGenotypes(data);
            if (opt.outputResults) gctb.outputResults(data, mcmcSampleVec, opt.bayesType, opt.noscale, opt.title);
        }
        else if (opt.analysisType == "LDmatrix") {
            readGenotypes = false;
            if (!opt.plinkLDtxtfile.empty() || !opt.plinkLDbinfile.empty()) {  // read from PLINK output file
                data.readPlinkAFfile(opt.plinkAFfile);
                if (!opt.plinkLDtxtfile.empty()) data.readPlinkLDtxtfile(opt.plinkLDtxtfile);
                if (!opt.plinkLDbinfile.empty()) data.readPlinkLDbinfile(opt.plinkLDbinfile);
                data.outputLDmatrix("full", opt.title, opt.writeLdmTxt);
            }
            else if (opt.ldmatrixFile.empty()) { // make LD matrix from genotypes (or from eigen data for block LDM)
                if (opt.outLDmatType == "block" && !opt.eigenMatrixFile.empty()) {
                    data.makeBlockLDmatrixFromEigen(opt.eigenMatrixFile, opt.eigenCutoff.maxCoeff(), opt.outLDmatType, opt.includeBlock, opt.title, opt.writeLdmTxt);
                } else {
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
                    data.mergeLdmInfo(opt.outLDmatType, opt.ldmatrixFile, true);
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
                else if (opt.outLDmatType == "block") { // resize block matrix
                    if (!opt.includeSnpFile.empty()) {
                        data.resizeBlockLDmatrix(opt.ldmatrixFile, opt.outLDmatType, opt.includeSnpFile, opt.title, opt.writeLdmTxt);
                    }
                    else if (opt.writeLdmTxt) {
                        data.outputBlockLDmatrixTxt(opt.ldmatrixFile, opt.includeBlock);
                    }
                }
                else if (opt.outLDmatType == "blockTri") { // only save the lower triangular matrix
                    data.convertToBlockTriangularMatrix(opt.ldmatrixFile, opt.writeLdmTxt, opt.title);
                }
                else if (opt.outLDmatType == "blockSparse") {
                    if (opt.includeSnpFile.empty()) {  // make block sparse LD matrix from block full LD matrix
                        data.readBlockLDmatrixAndMakeItSparse(opt.ldmatrixFile, opt.includeBlock, opt.chisqThreshold, opt.writeLdmTxt);
                    } else {  // resize block sparse matrix
                        data.resizeBlockLDmatrix(opt.ldmatrixFile, opt.outLDmatType, opt.includeSnpFile, opt.title, opt.writeLdmTxt);
                    }
                }
                else {  // resize existing non-block LD matrix
                    gctb.inputSnpInfo(data, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, "", opt.ldmatrixFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.genMapN, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.ldscoreFile, opt.windowFile, opt.multiLDmat, opt.excludeMHC, opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.rsqThreshold, opt.sampleOverlap, opt.imputeN, opt.noscale, opt.binSnp, opt.readLdmTxt);
                    data.resizeLDmatrix(opt.outLDmatType, opt.chisqThreshold, opt.windowWidth, opt.LDthreshold, opt.effpopNE, opt.cutOff, opt.genMapN);
                    data.outputLDmatrix(opt.outLDmatType, opt.title, opt.writeLdmTxt);
                }
            }
        }
        else if (opt.analysisType == "LDmatrixEigen") {
            readGenotypes = false;
            if (opt.eigenMatrixFile.empty()) { // perform eigen decomposition for the blocked LD matrices
                if (!opt.gwasSummaryFile.empty()) { // match LD ref SNPs with GWAS SNPs
                    gctb.inputSnpInfo(data, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile,
                                      opt.gwasSummaryFile, opt.ldmatrixFile, opt.ldBlockInfoFile,
                                      opt.includeChr, opt.excludeAmbiguousSNP,
                                      opt.annotationFile, opt.transpose,
                                      opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.ldscoreFile,
                                      opt.eigenCutoff.maxCoeff(), opt.excludeMHC,
                                      opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.rsqThreshold,
                                      opt.sampleOverlap, opt.imputeN, opt.noscale, opt.readLdmTxt, opt.imputeSummary, opt.keepSumstatIntactFile, opt.includeBlock, opt.skipSnpFile, opt.setZeroGwasForSkip, false);
                    data.resizeBlockLDmatrixAndDoEigenDecomposition(opt.ldmatrixFile, opt.eigenCutoff.maxCoeff(), 0.5, opt.title, opt.writeLdmTxt);
                } else {
                    //gctb.inputSnpInfo(data, opt.bedFile, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.mafmin, opt.mafmax, opt.noscale, readGenotypes);
                    //data.getEigenDataFromFullLDM(opt.title, opt.eigenCutoff);
                    data.readBlockLDmatrixAndDoEigenDecomposition(opt.ldmatrixFile, opt.includeBlock, opt.eigenCutoff.maxCoeff(), opt.writeLdmTxt);
                }
            }
            else { // merge existing eigen matrices
                //gctb.inputSnpInfo(data, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, "", opt.eigenMatrixFile, opt.ldBlockInfoFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.ldscoreFile, opt.eigenCutoff, opt.excludeMHC, opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.rsqThreshold, opt.sampleOverlap, opt.imputeN, opt.noscale, opt.readLdmTxt);
            }
        }
        else if (opt.analysisType == "sparseLDmatrixEigen") {
            readGenotypes = false;
            data.readSparseBlockLDmatrixAndDoEigenDecomposition(opt.ldmatrixFile, opt.includeBlock, opt.eigenCutoff.maxCoeff(), opt.writeLdmTxt);
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
                              opt.sampleOverlap, opt.imputeN, opt.noscale, opt.readLdmTxt, opt.imputeSummary, opt.keepSumstatIntactFile, opt.includeBlock, opt.skipSnpFile, opt.setZeroGwasForSkip);
        }
        else if (opt.analysisType == "MergeGwasSummary") {
            if (opt.outLDmatType == "block") {
                data.mergeBlockGwasSummary(opt.gwasSummaryFile, opt.title);
            }
        }
        else if (opt.analysisType == "Convert") {
            if (!opt.eigenMatrixFile.empty()) {
                data.readEigenMatrix(opt.eigenMatrixFile, opt.eigenCutoff.maxCoeff());
                data.inputMatchedSnpResults(opt.snpResFile);
                data.convert(opt.eigenMatrixFile, opt.includeSnpFile, opt.title);
            }
        }
        else if (opt.analysisType == "GetLD") {
            if (!opt.eigenMatrixFile.empty()) {
                data.readEigenMatrix(opt.eigenMatrixFile, opt.eigenCutoff.maxCoeff());
                data.getLDfromEigenMatrix(opt.eigenMatrixFile, opt.rsqThreshold, opt.title);
            }
        }
        else if (opt.analysisType == "GetLDfriends") {
            if (opt.rsqThreshold == 1.0) {
                opt.rsqThreshold = 0.5;  // default value
            }
            data.getLDfriends(opt.pairwiseLDfile, opt.rsqThreshold, opt.title);
        }
        else if (opt.analysisType == "SBayes" || opt.analysisType == "GWFM") {
            bool isBivariate = opt.gwasSummaryFile.find(',') != string::npos;
            float bestEigenCutoff = opt.eigenCutoff.size() ? opt.eigenCutoff[0] : 0.995f;
            if (!opt.ldmatrixFile.empty()) {
                gctb.inputSnpInfo(data, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, opt.gwasSummaryFile, opt.ldmatrixFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.genMapN, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.ldscoreFile, opt.windowFile, opt.multiLDmat, opt.excludeMHC, opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.rsqThreshold, opt.sampleOverlap, opt.imputeN, opt.noscale, opt.binSnp, opt.readLdmTxt);
            } else if (!opt.ldmBlockDir.empty()) {
                if (!opt.eigenMatrixFile.empty()) {
                    throw("Error: use only one of --ldm-eigen and --ldm-block (both specify block LD reference folders).");
                }
                data.mergeLdmInfo("block", opt.ldmBlockDir, false);
                gctb.inputSnpInfo(data, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile,
                                  opt.gwasSummaryFile, opt.ldmBlockDir, opt.ldBlockInfoFile,
                                  opt.includeChr, opt.excludeAmbiguousSNP,
                                  opt.annotationFile, opt.transpose,
                                  opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.ldscoreFile,
                                  opt.eigenCutoff.maxCoeff(), opt.excludeMHC,
                                  opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.rsqThreshold,
                                  opt.sampleOverlap, opt.imputeN, opt.noscale, opt.readLdmTxt, opt.imputeSummary, opt.keepSumstatIntactFile, opt.includeBlock, opt.skipSnpFile, opt.setZeroGwasForSkip, true);
                if (opt.analysisType == "GWFM") {
                    data.inputPairwiseLD(opt.ldmBlockDir + "/" + opt.pairwiseLDfile, 0.95);
                }
                cout << "Using full block LD matrices (block*.ldm.bin) from [" << opt.ldmBlockDir << "]; eigen low-rank W/Q not built." << endl;
            } else if (!opt.eigenMatrixFile.empty()) {  // low-rank model
                data.recomputeEigen = opt.recomputeEigen;
                data.mergeLdmInfo("block", opt.eigenMatrixFile, false); // if each block has its own .info file, then merge them
                gctb.inputSnpInfo(data, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile,
                                  opt.gwasSummaryFile, opt.eigenMatrixFile, opt.ldBlockInfoFile,
                                  opt.includeChr, opt.excludeAmbiguousSNP,
                                  opt.annotationFile, opt.transpose,
                                  opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.ldscoreFile,
                                  opt.eigenCutoff.maxCoeff(), opt.excludeMHC,
                                  opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.rsqThreshold,
                                  opt.sampleOverlap, opt.imputeN, opt.noscale, opt.readLdmTxt, opt.imputeSummary, opt.keepSumstatIntactFile, opt.includeBlock, opt.skipSnpFile, opt.setZeroGwasForSkip, false);
                if (opt.analysisType == "GWFM") {
                    data.inputPairwiseLD(opt.eigenMatrixFile+"/"+opt.pairwiseLDfile, 0.95);  // for TGS sampling
                }
                // Check if bivariate (comma-separated GWAS files) - skip tuning for bivariate
                if (isBivariate) {
                    // Skip tuning for bivariate - use first cutoff value
                    bestEigenCutoff = opt.eigenCutoff[0];
                    cout << "Using eigen cutoff " << bestEigenCutoff << " for bivariate analysis (tuning not yet supported for bivariate)." << endl;
                } else {
                    bestEigenCutoff = opt.eigenCutoff.size() > 1 ? gctb.tuneEigenCutoff(data, opt) : opt.eigenCutoff[0];
                    // Single cutoff: buildMMEeigen (via inputSnpInfo) already called readEigenMatrixBinaryFileAndMakeWandQ; repeating it would redo all work with no console progress (header messages are suppressed when wcorrBlocks is non-empty).
                    if (opt.eigenCutoff.size() > 1) {
                        data.readEigenMatrixBinaryFileAndMakeWandQ(opt.eigenMatrixFile, bestEigenCutoff, data.gwasEffectInBlock, data.nGWASblock, opt.noscale, false);
                    }
                }
                if (opt.writeWandQ) data.outputWandQ("w_and_Q");
                //data.readEigenMatrixBinaryFile(opt.eigenMatrixFile, bestEigenCutoff);
                //data.constructWandQ(data.gwasEffectInBlock, data.numKeptInds);
            } else {
                gctb.inputSnpInfo(data, opt.bedFile, opt.gwasSummaryFile, opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.sampleOverlap, opt.imputeN, opt.noscale);
            }
            
            if (!opt.fitSnpsAsFixedEffectsFile.empty())
                data.readFixedEffectSnpFile(opt.fitSnpsAsFixedEffectsFile);
            
            data.label = opt.title;
//            if (opt.numChains > 1) {
//                vector<McmcSamples*> mcmcSampleVec = gctb.multi_chain_mcmc(data, opt.bayesType, opt.windowWidth, opt.heritability, opt.propVarRandom, opt.pi, opt.piAlpha, opt.piBeta, opt.estimatePi, opt.pis, opt.gamma, opt.phi, opt.kappa, opt.algorithm, opt.snpFittedPerWindow, opt.varS, opt.S, opt.overdispersion, opt.estimatePS, opt.icrsq, opt.spouseCorrelation, opt.diagnosticMode, opt.robustMode, opt.numChains, opt.chainLength, opt.burnin, opt.thin, opt.outputFreq, opt.title, opt.writeBinPosterior, opt.writeTxtPosterior);
//                if (opt.outputResults) gctb.outputResults(data, mcmcSampleVec, opt.bayesType, opt.noscale, opt.title);
//            } else {
                
            
            if (opt.nDistAuto) {
                if (opt.nDistAutoByPred && !opt.eigenMatrixFile.empty() && !isBivariate && opt.ldmBlockDir.empty())
                    gctb.findBestFitModelByPredictionAccuracy(data, opt, bestEigenCutoff);
                else {
                    if (opt.nDistAutoByPred) {
                        if (opt.eigenMatrixFile.empty())
                            cout << "Note: --n-dist-auto-pred requires --ldm-eigen; using heritability-based --n-dist-auto instead." << endl;
                        else if (isBivariate)
                            cout << "Note: --n-dist-auto-pred is not supported for bivariate analysis; using heritability-based --n-dist-auto instead." << endl;
                    }
                    gctb.findBestFitModel(data, opt);
                }
            }
            
            Model *model = gctb.buildModel(data, opt, opt.bedFile, opt.gwasSummaryFile, opt.bayesType, opt.windowWidth,
                                           opt.heritability, opt.propVarRandom, opt.pi, opt.piAlpha, opt.piBeta, opt.estimatePi, opt.noscale, opt.pis, opt.piPar, opt.gamma, opt.estimateSigmaSq, opt.phi, opt.kappa,
                                           opt.algorithm, opt.snpFittedPerWindow, opt.varS, opt.S, opt.overdispersion, opt.estimatePS, opt.icrsq, opt.spouseCorrelation, opt.diagnosticMode, opt.hsqPercModel, opt.perSnpGV, opt.robustMode, opt.nDistAuto, opt.estimateRsqEnrich);
            
//            vector<McmcSamples*> mcmcSampleVec = gctb.runMcmc(*model, opt.numChains, opt.chainLength, opt.burnin, opt.thin,
//                                                              opt.outputFreq, opt.title, opt.writeBinPosterior, opt.writeTxtPosterior);

            MCMC *mcmc = new MCMC();
            vector<McmcSamples*> mcmcSampleVec = mcmc->run(*model, opt.numChains, opt.chainLength, opt.burnin, opt.thin, true,
                                                           opt.outputFreq, opt.title, opt.writeBinPosterior, opt.writeTxtPosterior);

            if (mcmc->action == mcmc->restart_and_use_robust_model) {
                cout << "\nRestarting MCMC with a more robust parameterisation for SBayes" << opt.bayesType << " ..." << endl;
                cout << "Please refer to our website (https://cnsgenomics.com/software/gctb) for more information." << endl;
                opt.robustMode = true;
                model = gctb.buildModel(data, opt, opt.bedFile, opt.gwasSummaryFile, opt.bayesType, opt.windowWidth,
                                        opt.heritability, opt.propVarRandom, opt.pi, opt.piAlpha, opt.piBeta, opt.estimatePi, opt.noscale, opt.pis, opt.piPar, opt.gamma, opt.estimateSigmaSq, opt.phi, opt.kappa,
                                        opt.algorithm, opt.snpFittedPerWindow, opt.varS, opt.S, opt.overdispersion, opt.estimatePS, opt.icrsq, opt.spouseCorrelation, opt.diagnosticMode, opt.hsqPercModel, opt.perSnpGV, opt.robustMode, opt.nDistAuto, opt.estimateRsqEnrich);
                mcmc = new MCMC();
                vector<McmcSamples*> mcmcSampleVec = mcmc->run(*model, opt.numChains, opt.chainLength, opt.burnin, opt.thin, true,
                                                               opt.outputFreq, opt.title, opt.writeBinPosterior, opt.writeTxtPosterior);
            }
            
            if (opt.outputResults) gctb.outputResults(data, mcmcSampleVec, opt.bayesType, opt.noscale, opt.title);
            
            if (opt.analysisType == "GWFM") {
                McmcSamples *snpEffects = gctb.inputMcmcSamples(opt.title, "SnpEffects", "bin");
                data.inputNewSnpResults(opt.title + ".snpRes");
                data.inputPairwiseLD(opt.eigenMatrixFile+"/"+opt.pairwiseLDfile, opt.rsqThreshold);
                if (!opt.geneMapFile.empty()) data.readGeneMapFile(opt.geneMapFile, opt.flank, opt.genomeBuild);
                gctb.calcCredibleSets(data, *snpEffects, opt.pipThreshold, opt.pepThreshold, opt.title);
            }
        }
        else if (opt.analysisType == "CalcRsqEnrichment") {
            data.inputNewSnpResults(opt.snpResFile + ".snpRes");
            data.readAnnotationFile(opt.annotationFile, opt.transpose, true);
            data.setAnnoInfoVec();
            data.calcMarignalEnrichmentJackknife("rsq", opt.title);
            data.calcJointEnrichmentJackknifeLM("rsq", opt.title);
        }
        else if (opt.analysisType == "CalcHsqEnrichment") {
            data.inputNewSnpResults(opt.snpResFile + ".snpRes");
            data.readAnnotationFile(opt.annotationFile, opt.transpose, true);
            data.setAnnoInfoVec();
            data.calcMarignalEnrichmentJackknife("hsq", opt.title);
            data.calcJointEnrichmentJackknifeLM("hsq", opt.title);
        }
        else if (opt.analysisType == "CalcPipEnrichment") {
            data.inputNewSnpResults(opt.snpResFile + ".snpRes");
            data.readAnnotationFile(opt.annotationFile, opt.transpose, true);
            data.setAnnoInfoVec();
            data.calcMarignalEnrichmentJackknife("pip", opt.title);
            data.calcJointEnrichmentJackknifeLM("pip", opt.title);
        }
        else if (opt.analysisType == "CalcRsqEnrichmentPart") {
            data.inputNewSnpResults(opt.snpResFile + ".snpRes");
            data.readAnnotationFile(opt.annotationFile, opt.transpose, true);
            data.setAnnoInfoVec();
            data.calcMarignalEnrichmentJackknife("rsq", opt.title);
            data.calcJointEnrichmentJackknife("rsq", opt.title);
        }
        else if (opt.analysisType == "CalcHsqEnrichmentPart") {
            data.inputNewSnpResults(opt.snpResFile + ".snpRes");
            data.readAnnotationFile(opt.annotationFile, opt.transpose, true);
            data.setAnnoInfoVec();
            data.calcMarignalEnrichmentJackknife("hsq", opt.title);
            data.calcJointEnrichmentJackknife("hsq", opt.title);
        }
        else if (opt.analysisType == "CalcPipEnrichmentPart") {
            data.inputNewSnpResults(opt.snpResFile + ".snpRes");
            data.readAnnotationFile(opt.annotationFile, opt.transpose, true);
            data.setAnnoInfoVec();
            data.calcMarignalEnrichmentJackknife("pip", opt.title);
            data.calcJointEnrichmentJackknife("pip", opt.title);
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
        else if (opt.analysisType == "WindowPIP") {
            McmcSamples *snpEffects = gctb.inputMcmcSamples(opt.mcmcSampleFile, "SnpEffects", "bin");
            gctb.getWindowPIP(data, *snpEffects, opt.mcmcSampleFile + ".snpRes", opt.windowWidth, 0.5*opt.windowWidth, opt.title);
        }
        else if (opt.analysisType == "CS") {
//            if (!opt.eigenMatrixFile.empty()) {
//                McmcSamples *snpEffects = gctb.inputMcmcSamples(opt.mcmcSampleFile, "SnpEffects", "bin");
//                gctb.calcCredibleSets(data, opt.mcmcSampleFile + ".snpRes", *snpEffects, opt.eigenMatrixFile, opt.eigenCutoff.maxCoeff(), opt.pipThreshold, opt.title);
//            } 
            if (!opt.pairwiseLDfile.empty()) {
                McmcSamples *snpEffects = gctb.inputMcmcSamples(opt.mcmcSampleFile, "SnpEffects", "bin");
                data.inputNewSnpResults(opt.mcmcSampleFile + ".snpRes");
                data.inputPairwiseLD(opt.eigenMatrixFile+"/"+opt.pairwiseLDfile, opt.rsqThreshold);
                if (!opt.geneMapFile.empty()) data.readGeneMapFile(opt.geneMapFile, opt.flank, opt.genomeBuild);
                gctb.calcCredibleSets(data, *snpEffects, opt.pipThreshold, opt.pepThreshold, opt.title);
            }
            else if (!opt.ldfriendFile.empty()) {
                McmcSamples *snpEffects = gctb.inputMcmcSamples(opt.mcmcSampleFile, "SnpEffects", "bin");
                data.inputNewSnpResults(opt.mcmcSampleFile + ".snpRes");
                data.inputLDfriends(opt.eigenMatrixFile+"/"+opt.ldfriendFile);
                if (!opt.geneMapFile.empty()) data.readGeneMapFile(opt.geneMapFile, opt.flank, opt.genomeBuild);
                gctb.calcCredibleSets(data, *snpEffects, opt.pipThreshold, opt.pepThreshold, opt.title);
            } else {
                int windowWidth = 100000;
                if (opt.windowWidth) windowWidth = opt.windowWidth;
                McmcSamples *snpEffects = gctb.inputMcmcSamples(opt.mcmcSampleFile, "SnpEffects", "bin");
                gctb.calcCredibleSets(data, opt.mcmcSampleFile + ".snpRes", *snpEffects, opt.pipThreshold, opt.pepThreshold, windowWidth, opt.title);
            }
        }
        else if (opt.analysisType == "Bin2Txt") {
            McmcSamples *mcmcSamples = gctb.inputMcmcSamples(opt.mcmcSampleFile, opt.label, "bin");
            mcmcSamples->writeMatSpTxt(opt.mcmcSampleFile);
        }
        else if (opt.analysisType == "Print") {
            if (!opt.eigenMatrixFile.empty()) {
                data.readEigenMatrix(opt.eigenMatrixFile, opt.eigenCutoff.maxCoeff(), true, true, opt.title);
            }
        }
        else if (opt.analysisType == "Predict") {
            readGenotypes = true;
            gctb.inputIndInfo(data, opt.bedFile, opt.phenotypeFile, opt.keepIndFile, opt.keepIndMax,
                               opt.mphen, opt.covariateFile, opt.randomCovariateFile, opt.residualDiagFile);
            gctb.inputSnpInfo(data, opt.bedFile, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.ldBlockInfoFile, opt.includeBlock, opt.annotationFile, opt.transpose, opt.continuousAnnoFile, opt.flank, opt.eQTLFile, opt.mafmin, opt.mafmax, opt.noscale, readGenotypes);
            
            data.inputMatchedSnpResults(opt.snpResFile);
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
                    vector<McmcSamples*> mcmcSampleVec = gctb.runMcmc(*model, 1, opt.chainLength, opt.burnin, opt.thin,
                                                                      opt.outputFreq, opt.title + ".stage1", opt.writeBinPosterior, opt.writeTxtPosterior);
                    gctb.saveMcmcSamples(mcmcSampleVec, opt.title + ".stage1");
                    gctb.outputResults(data, mcmcSampleVec, "C", true, opt.title + ".stage1");
                    xci.outputResults(data, mcmcSampleVec, "C", opt.title + ".stage1");
                    delete model;
                    // Stage 2: estimate GxS using both male and female data
                    model = xci.buildModelStageTwo(data, "Cgxs", opt.heritability, opt.pi, opt.piPar, opt.estimatePi, opt.piNDC, opt.piNDCpar, opt.estimatePiNDC, opt.title + ".stage1.snpRes", opt.piGxE, opt.estimatePiGxE);
                    mcmcSampleVec = gctb.runMcmc(*model, 1, opt.chainLength, opt.burnin, opt.thin,
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
                        vector<McmcSamples*> mcmcSampleVec = gctb.runMcmc(*model, 1, opt.chainLength, opt.burnin, opt.thin,
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
    catch (const std::string &err_msg) {
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
