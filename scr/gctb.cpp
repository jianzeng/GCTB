//
//  gctb.cpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "gctb.hpp"
#include "stat.hpp"
#include <Eigen/Core>
#include <cmath>

void GCTB::inputIndInfo(Data &data, const string &bedFile, const string &phenotypeFile, const string &keepIndFile, const unsigned keepIndMax, const unsigned mphen, const string &covariateFile, const string &randomCovariateFile, const string &residualDiagFile){
    data.readFamFile(bedFile + ".fam");
    data.readPhenotypeFile(phenotypeFile, mphen);
    data.readCovariateFile(covariateFile);
    data.readRandomCovariateFile(randomCovariateFile);
    data.readResidualDiagFile(residualDiagFile);
    data.keepMatchedInd(keepIndFile, keepIndMax);
}

void GCTB::inputSnpInfo(Data &data, const string &bedFile, const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile, const unsigned includeChr, const bool excludeAmbiguousSNP, const string &skeletonSnpFile, const string &geneticMapFile,  const string &ldBlockInfoFile, const unsigned includeBlock, const string &annotationFile, const bool transpose, const string &continuousAnnoFile, const unsigned flank, const string &eQTLFile, const float mafmin, const float mafmax, const bool noscale, const bool readGenotypes){
    data.readBimFile(bedFile + ".bim");
    if (!includeSnpFile.empty()) data.includeSnp(includeSnpFile);
    if (!excludeSnpFile.empty()) data.excludeSnp(excludeSnpFile);
    if (includeChr) data.includeChr(includeChr);
    if (excludeAmbiguousSNP) data.excludeAmbiguousSNP();
//    if (mafmin || mafmax) data.excludeSNPwithMaf(mafmin, mafmax);  // need to read in genotype data first
    if (!excludeRegionFile.empty()) data.excludeRegion(excludeRegionFile);
    if (!skeletonSnpFile.empty()) data.includeSkeletonSnp(skeletonSnpFile);
    if (!geneticMapFile.empty()) data.readGeneticMapFile(geneticMapFile);
    if (!annotationFile.empty())
        data.readAnnotationFile(annotationFile, transpose, true);
    else if (!continuousAnnoFile.empty())
        data.readAnnotationFileFormat2(continuousAnnoFile, flank*1000, eQTLFile);
    if (!ldBlockInfoFile.empty()) data.readLDBlockInfoFile(ldBlockInfoFile);
    if (includeBlock) data.includeBlock(includeBlock);
    data.includeMatchedSnp();
    if (data.numAnnos) data.setAnnoInfoVec();
//    data.makeWindowAnno(annotationFile, 5e5);
    if (readGenotypes) data.readBedFile(noscale, bedFile + ".bed");
}

void GCTB::inputSnpInfo(Data &data, const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile, const string &gwasSummaryFile, const string &ldmatrixFile, const unsigned includeChr, const bool excludeAmbiguousSNP, const string &skeletonSnpFile, const string &geneticMapFile, const float genMapN, const string &annotationFile, const bool transpose, const string &continuousAnnoFile, const unsigned flank, const string &eQTLFile, const string &ldscoreFile, const string &windowFile, const bool multiLDmat, const bool excludeMHC, const float afDiff, const float mafmin, const float mafmax, const float pValueThreshold, const float rsqThreshold, const bool sampleOverlap, const bool imputeN, const bool noscale, const bool binSnp, const bool readLDMfromTxtFile){
    if (multiLDmat)
        data.readMultiLDmatInfoFile(ldmatrixFile);
    else
        data.readLDmatrixInfoFile(ldmatrixFile + ".info");
    if (!includeSnpFile.empty()) data.includeSnp(includeSnpFile);
    if (!excludeSnpFile.empty()) data.excludeSnp(excludeSnpFile);
    if (includeChr) data.includeChr(includeChr);
    if (excludeAmbiguousSNP) data.excludeAmbiguousSNP();
    if (!excludeRegionFile.empty()) data.excludeRegion(excludeRegionFile);
    if (excludeMHC) data.excludeMHC();
    if (!skeletonSnpFile.empty()) data.includeSkeletonSnp(skeletonSnpFile);
    if (!geneticMapFile.empty()) data.readGeneticMapFile(geneticMapFile);
    if (!annotationFile.empty())
        data.readAnnotationFile(annotationFile, transpose, true);
    else if (!continuousAnnoFile.empty())
        data.readAnnotationFileFormat2(continuousAnnoFile, flank*1000, eQTLFile);
    if (!ldscoreFile.empty()) data.readLDscoreFile(ldscoreFile);
    if (!windowFile.empty()) data.readWindowFile(windowFile);
    if (!gwasSummaryFile.empty()) data.readGwasSummaryFile(gwasSummaryFile, afDiff, mafmin, mafmax, pValueThreshold, imputeN, true);
    data.includeMatchedSnp();
    if (readLDMfromTxtFile) {
        data.readLDmatrixTxtFile(ldmatrixFile + ".txt");
    } else {
//        if (geneticMapFile.empty()) {
            if (multiLDmat)
                data.readMultiLDmatBinFile(ldmatrixFile);
            else
                data.readLDmatrixBinFile(ldmatrixFile + ".bin");
//        } else {
//            if (multiLDmat)
//                data.readMultiLDmatBinFileAndShrink(ldmatrixFile, genMapN);
//            else
//                data.readLDmatrixBinFileAndShrink(ldmatrixFile + ".bin");
//        }
    }
    
    if (rsqThreshold < 1.0 && !binSnp) {
        data.filterSnpByLDrsq(rsqThreshold);
        data.includeMatchedSnp();
        if (geneticMapFile.empty()) {  // need to read LD data again after LD filtering
            if (multiLDmat)
                data.readMultiLDmatBinFile(ldmatrixFile);
            else
                data.readLDmatrixBinFile(ldmatrixFile + ".bin");
        } else {
            if (multiLDmat)
                data.readMultiLDmatBinFileAndShrink(ldmatrixFile, genMapN);
            else
                data.readLDmatrixBinFileAndShrink(ldmatrixFile + ".bin");
        }
    }
    if (!gwasSummaryFile.empty()) data.buildSparseMME(sampleOverlap, noscale);
    if (!windowFile.empty()) data.binSnpByWindowID();
}

// this function read eigen matrices
void GCTB::inputSnpInfo(Data &data, const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile,
                        const string &gwasSummaryFile, const string &eigenMatrixFile, const string &ldBlockInfoFile,
                        const unsigned includeChr, const bool excludeAmbiguousSNP,
                        const string &annotationFile, const bool transpose,
                        const string &continuousAnnoFile, const unsigned flank, const string &eQTLFile, const string &ldscoreFile,
                        const float eigenCutoff, const bool excludeMHC,
                        const float afDiff, const float mafmin, const float mafmax, const float pValueThreshold, const float rsqThreshold,
                        const bool sampleOverlap, const bool imputeN, const bool noscale, const bool readLDMfromTxtFile, const bool imputeSummary, const unsigned includeBlock, const string &skipSnpFile, const bool setZeroGwasForSkip, const bool buildMME){
    data.readEigenMatrix(eigenMatrixFile, eigenCutoff);
    if (!includeSnpFile.empty()) data.includeSnp(includeSnpFile);
    if (!excludeSnpFile.empty()) data.excludeSnp(excludeSnpFile);
    if (includeChr) data.includeChr(includeChr);
    if (includeBlock) data.includeBlock(includeBlock);
    if (excludeAmbiguousSNP) data.excludeAmbiguousSNP();
    if (!excludeRegionFile.empty()) data.excludeRegion(excludeRegionFile);
    if (excludeMHC) data.excludeMHC();
    if (!skipSnpFile.empty()) data.skipSnp(skipSnpFile);
    if (setZeroGwasForSkip && skipSnpFile.empty())
        cout << "Warning: --set-zero-gwas-for-skip has no effect without --skip <file>." << endl;
    if (!annotationFile.empty())
        data.readAnnotationFile(annotationFile, transpose, true);
    else if (!continuousAnnoFile.empty())
        data.readAnnotationFileFormat2(continuousAnnoFile, flank*1000, eQTLFile);
    if (!ldscoreFile.empty()) data.readLDscoreFile(ldscoreFile);
    if (!gwasSummaryFile.empty()) {
        bool removeOutlierN = imputeSummary;
        vector<string> gwasFiles;
        Gadget::Tokenizer token;
        token.getTokens(gwasSummaryFile, ",");
        for (unsigned i=0; i<token.size(); ++i) gwasFiles.push_back(token[i]);
        
        if (gwasFiles.size() == 1) {
            data.readGwasSummaryFile(gwasFiles[0], afDiff, mafmin, mafmax, pValueThreshold, imputeN, removeOutlierN);
        } else if (gwasFiles.size() == 2) {
            data.readBivariateGwasSummaryFile(gwasFiles[0], gwasFiles[1], afDiff, mafmin, mafmax, pValueThreshold, imputeN, removeOutlierN);
        } else {
            throw("Error: only support bivariate analysis for now.");
        }

        if (imputeSummary) {
            //data.includeMatchedBlocks();
            //data.scaleGwasEffects();
            data.readEigenMatrixBinaryFile(eigenMatrixFile, eigenCutoff);
            data.impG(includeBlock);
            if (setZeroGwasForSkip && !skipSnpFile.empty()) data.applySetZeroGwasForSkip();
            return;
        }
        data.includeMatchedSnp();
        if (setZeroGwasForSkip && !skipSnpFile.empty()) data.applySetZeroGwasForSkip();
    }
    

    /// partition ld into blocks
//    if(!ldBlockInfoFile.empty()) data.readLDBlockInfoFile(ldBlockInfoFile);
        
    if(!gwasSummaryFile.empty() && buildMME) {
        Gadget::Tokenizer token;
        token.getTokens(gwasSummaryFile, ",");
        if (token.size() == 1)
            data.buildMMEeigen(eigenMatrixFile, sampleOverlap, eigenCutoff, noscale);
        else if (token.size() == 2)
            data.buildMMEeigenBivariate(eigenMatrixFile, sampleOverlap, eigenCutoff, noscale);
    }
}


void GCTB::inputSnpInfo(Data &data, const string &bedFile, const string &gwasSummaryFile, const float afDiff, const float mafmin, const float mafmax, const float pValueThreshold, const bool sampleOverlap, const bool imputeN, const bool noscale){
    data.readFamFile(bedFile + ".fam");
    data.readBimFile(bedFile + ".bim");

    data.keptIndInfoVec = data.makeKeptIndInfoVec(data.indInfoVec);
    data.numKeptInds =  (unsigned) data.keptIndInfoVec.size();
    
    data.readGwasSummaryFile(gwasSummaryFile, afDiff, mafmin, mafmax, pValueThreshold, imputeN, true);
    data.includeMatchedSnp();
    data.readBedFile(noscale, bedFile + ".bed");
    data.buildSparseMME(sampleOverlap, noscale);
}

Model* GCTB::buildModel(Data &data, const Options &opt, const string &bedFile, const string &gwasFile, const string &bayesType, const unsigned windowWidth,
                        const float heritability, const float propVarRandom, const float pi, const float piAlpha, const float piBeta, const bool estimatePi, const bool noscale,
                        const VectorXf &pis, const VectorXf &piPar, const VectorXf &gamma, const bool estimateSigmaSq,
                        const float phi, const float kappa, const string &algorithm, const unsigned snpFittedPerWindow,
                        const float varS, const vector<float> &S, const float overdispersion, const bool estimatePS,
                        const float icrsq, const float spouseCorrelation, const bool diagnosticMode, const bool hsqPercModel, const bool perSnpGV, const bool robustMode, const bool nDistAuto, const bool estimateRsqEnrich){
    
    data.initVariances(heritability, propVarRandom);

    if (opt.numChains > 1) {
        if (bayesType == "R")
            return new MultiChainSBayesR(data, opt);
        else if (bayesType == "RC")
            return new MultiChainSBayesRC(data, opt);
        else if (bayesType == "RD")
            return new MultiChainSBayesRD(data, opt);
        else if (bayesType == "S")
            return new MultiChainSBayesS(data, opt);
        else
            throw(" Error: " + bayesType + " is not available for multi-chain analysis.");
    }
    if (!gwasFile.empty()) {
        if (data.numAnnos) {
            if (bayesType == "S")
                return new StratApproxBayesS(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, varS, S, algorithm, robustMode, noscale);
            else if (bayesType == "RC")
                return new ApproxBayesRC(data, data.lowRankModel, data.varGenotypic, data.varResidual, pis, piPar, gamma, estimatePi, noscale, hsqPercModel, robustMode, estimateRsqEnrich, algorithm);
            else if (bayesType == "RD")
                return new ApproxBayesRD(data, data.lowRankModel, data.varGenotypic, data.varResidual, pis, piPar, gamma, estimatePi, noscale, hsqPercModel, robustMode, estimateRsqEnrich, algorithm);
            else if (bayesType == "APP") {
                VectorXf nLociAnno(data.numAnnos);
                for(unsigned i=0; i<data.numAnnos; ++i) nLociAnno[i] = data.numSnpAnnoVec[i];
                return new ApproxBayesAPP(data, data.lowRankModel, data.varGenotypic, data.varGenotypic, data.varResidual, data.varResidual, nLociAnno, estimatePi, true, true);
            }
            else
                throw(" Error: Wrong bayes type: " + bayesType + " in the annotation-stratified summary-data-based Bayesian analysis.");
        }
        else {
            if (bayesType == "C")
                return new ApproxBayesC(data, data.lowRankModel, data.varGenotypic, data.varResidual, data.varRandom, pi, piAlpha, piBeta, estimatePi, noscale, robustMode);
            else if (bayesType == "B")
            return new ApproxBayesB(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, noscale);
            else if (bayesType == "S")
                return new ApproxBayesS(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, varS, S, algorithm, noscale);
            else if (bayesType == "ST")
                return new ApproxBayesST(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, varS, S, true, noscale);
            else if (bayesType == "T")
                return new ApproxBayesST(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, varS, S, false, noscale);
            else if (bayesType == "SMix")
                return new ApproxBayesSMix(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, varS, S, noscale);
            else if (bayesType == "R")
                return new ApproxBayesR(data, data.lowRankModel, data.varGenotypic, data.varResidual, pis, piPar, gamma, estimatePi, noscale, hsqPercModel, robustMode, algorithm);
            else if (bayesType == "RS")
                return new ApproxBayesRS(data, data.lowRankModel, data.varGenotypic, data.varResidual, pis, piPar, gamma, estimatePi, varS, S, noscale, hsqPercModel, robustMode, algorithm);
            else if (bayesType == "APP")
                return new ApproxBayesAPP(data, data.lowRankModel, data.varGenotypic, data.varGenotypic, data.varResidual, data.varResidual, VectorXf(), estimatePi, true, true);
            else
                throw(" Error: Wrong bayes type: " + bayesType + " in the summary-data-based Bayesian analysis.");
        }
    }
    if (data.numAnnos) {
        if (bayesType == "RC") {
            data.readBedFile(noscale, bedFile + ".bed");
            return new BayesRC(data, data.varGenotypic, data.varResidual, data.varRandom, pis, piPar, gamma, estimatePi, noscale, hsqPercModel, estimateRsqEnrich, "Gibbs");
        }
        else
            throw(" Error: Wrong bayes type: " + bayesType + " in the annotation-stratified Bayesian analysis.");
    }
    if (bayesType == "B") {
        data.readBedFile(noscale, bedFile + ".bed");
        return new BayesB(data, data.varGenotypic, data.varResidual, data.varRandom, pi, piAlpha, piBeta, estimatePi, noscale);
    }
    if (bayesType == "C") {
        data.readBedFile(noscale, bedFile + ".bed");
        return new BayesC(data, data.varGenotypic, data.varResidual, data.varRandom, pi, piAlpha, piBeta, estimatePi, noscale, algorithm);
    } 
    if (bayesType == "R") {
        data.readBedFile(noscale, bedFile + ".bed");
        return new BayesR(data, data.varGenotypic, data.varResidual, data.varRandom, pis, piPar, gamma, estimatePi, noscale, hsqPercModel, algorithm);
    }
    else if (bayesType == "S") {
        data.readBedFile(noscale, bedFile + ".bed");
        return new BayesS(data, data.varGenotypic, data.varResidual, data.varRandom, pi, piAlpha, piBeta, estimatePi, varS, S, algorithm, noscale);
    }
    else if (bayesType == "SMix") {
        data.readBedFile(noscale, bedFile + ".bed");
        return new BayesSMix(data, data.varGenotypic, data.varResidual, data.varRandom, pi, piAlpha, piBeta, estimatePi, varS, S, algorithm);
    }
    else if (bayesType == "N") {
        data.readBedFile(noscale, bedFile + ".bed");
        data.getNonoverlapWindowInfo(windowWidth);
        return new BayesN(data, data.varGenotypic, data.varResidual, data.varRandom, pi, piAlpha, piBeta, estimatePi, noscale, snpFittedPerWindow);
    }
    else if (bayesType == "NS") {
        data.readBedFile(noscale, bedFile + ".bed");
        data.getNonoverlapWindowInfo(windowWidth);
        return new BayesNS(data, data.varGenotypic, data.varResidual, data.varRandom, pi, piAlpha, piBeta, estimatePi, varS, S, snpFittedPerWindow, algorithm);
    }
    else if (bayesType == "RS") {
        data.readBedFile(noscale, bedFile + ".bed");
        return new BayesRS(data, data.varGenotypic, data.varResidual, data.varRandom, pis, piPar, gamma, estimatePi, varS, S, noscale, hsqPercModel, algorithm);
    }
    else if (bayesType == "Cap") {
        //data.readBedFile(bedFile + ".bed");
        data.buildSparseMME(bedFile + ".bed", windowWidth);
        return new ApproxBayesC(data, data.lowRankModel, data.varGenotypic, data.varResidual, data.varRandom, pi, piAlpha, piBeta, estimatePi, noscale, robustMode);
    }
    else if (bayesType == "Sap") {
        data.buildSparseMME(bedFile + ".bed", windowWidth);
        return new ApproxBayesS(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, varS, S, algorithm, noscale);
    }
    else {
        throw(" Error: Wrong bayes type: " + bayesType);
    }
}

vector<McmcSamples*> GCTB::runMcmc(Model &model, const unsigned numChains, const unsigned chainLength, const unsigned burnin, const unsigned thin, const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior){
    MCMC mcmc;
    return mcmc.run(model, numChains, chainLength, burnin, thin, true, outputFreq, title, writeBinPosterior, writeTxtPosterior);
}

void GCTB::findBestFitModel(Data &data, Options &opt){
    cout << "\nComparing models with different number of components ..." << endl;
    if (opt.bayesType != "R" && opt.bayesType != "RC") {
        throw(" Error: --n-dist-auto is available only in R or RC model. Current model is " + opt.bayesType + ".");
    }
    data.initVariances(opt.heritability, opt.propVarRandom);
    
    MultiModelSBayesR model(data, opt);
    
    unsigned numChains = 1;
    unsigned chainLength = 500;
    unsigned burnin = 100;
    bool print = true;
    bool writeBinPosterior = false;
    bool writeTxtPosterior = false;
    
    MCMC mcmc;
    vector<McmcSamples*> mcmcSampleVec = mcmc.run(model, numChains, chainLength, burnin, opt.thin, print, opt.outputFreq, opt.title, writeBinPosterior, writeTxtPosterior);

    vector<float> hsqMeanVec;
    vector<float> hsqSDVec;
    map<float, int, std::greater<float> > hsqMap;
    unsigned idx = 0;
    for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
        McmcSamples *mcmcSamples = mcmcSampleVec[i];
        Gadget::Tokenizer token;
        token.getTokens(mcmcSamples->label, "_");
        if (token.front() == "hsq") {
            hsqMeanVec.push_back(mcmcSamples->mean()[0]);
            hsqSDVec.push_back(mcmcSamples->sd()[0]);
            hsqMap[hsqMeanVec[idx]] = idx;
            ++idx;
        }
    }
    map<float, int, std::greater<float> >::iterator it, it2;
    while (hsqMap.size() > 1) {
        it = hsqMap.begin();
        it2 = it;
        ++it2;
        if (it->first - hsqSDVec[it->second] > it2->first) break;
        else {
            if (it->second < it2->second) hsqMap.erase(it);
            else hsqMap.erase(it2);
        }
    }
    it = hsqMap.begin();
    unsigned selectedModelIdx = it->second;
    
    opt.numDist = model.modelVec[selectedModelIdx]->gamma.values.size();
    opt.gamma = model.modelVec[selectedModelIdx]->gamma.values;
    opt.pis = model.modelVec[selectedModelIdx]->Pis.values;

    cout << "\nModel " << selectedModelIdx+1 << " (" << opt.numDist << "-component model) is selected because a more complex model did not explain a significantly higher SNP-based heritability." << endl;
    
}

void GCTB::findBestFitModelByPredictionAccuracy(Data &data, Options &opt, const float eigenCutoffUsed){
    const float predRelMin = 1.05f;
    
    cout << "\nComparing models with different numbers of mixture components ..." << endl;
    cout << "Criterion (combined): keep K over K-1 only if (1) hsq is significantly higher than for the simpler model (mean[K] - SD[K] > mean[K-1], as in --n-dist-auto), and (2) pseudo CV r[K] > " << predRelMin << " * r[K-1]." << endl;
    if (opt.bayesType != "R" && opt.bayesType != "RC") {
        throw(" Error: --n-dist-auto-pred is available only in R or RC model. Current model is " + opt.bayesType + ".");
    }
    if (!data.lowRankModel || opt.eigenMatrixFile.empty()) {
        throw(" Error: --n-dist-auto-pred requires a low-rank eigen LD matrix (--ldm-eigen).");
    }
    if (!data.pseudoGwasNtrnBlock.size()) {
        throw(" Error: --n-dist-auto-pred requires pseudo summary statistics from the eigen data build (buildMMEeigen).");
    }
    
    unsigned numModels = opt.numDist - 1;
    if (!numModels) {
        throw(" Error: numDist must be at least 2 for automatic model selection.");
    }
    
    // Stat::engine is thread_local; seedEngine() only seeds the current thread. Parallel
    // MultiModelSBayesR would use unseeded per-thread RNG → different hsq/r for R vs RC runs.
    // Nested OpenMP + Eigen threads can also interleave cout from workers vs main (misplaced lines).
    // Force single-thread OpenMP/Eigen, disable nested parallel regions, re-seed RNG.
    struct OmpEigenRestore {
        int ompThreads;
        int dyn;
        int maxLevels;
        int eigenThreads;
        OmpEigenRestore() {
            ompThreads = omp_get_max_threads();
            dyn = omp_get_dynamic();
            maxLevels = omp_get_max_active_levels();
            eigenThreads = Eigen::nbThreads();
            omp_set_dynamic(0);
            omp_set_max_active_levels(1);
            omp_set_num_threads(1);
            Eigen::setNbThreads(1);
        }
        ~OmpEigenRestore() {
            Eigen::setNbThreads(eigenThreads);
            omp_set_num_threads(ompThreads);
            omp_set_max_active_levels(maxLevels);
            omp_set_dynamic(dyn);
        }
    };
    OmpEigenRestore ompEigenGuard;
    unsigned rngBase = (unsigned)(opt.seed ? opt.seed : 011415);
    rngBase ^= 0x6E446973u;
    if (!rngBase) rngBase = 4877u;
    Stat::seedEngine((int)(rngBase & 0x7FFFFFFFu));
    
    Gadget::Timer timer;
    timer.setTime();
    
    // --- (1) Joint multi-model MCMC on real GWAS for hsq (same as findBestFitModel) ---
    data.initVariances(opt.heritability, opt.propVarRandom);
    MultiModelSBayesR jointModel(data, opt);
    // Constructor sets omp_set_max_active_levels(2) for nested model parallelism — re-lock for pseudo CV.
    omp_set_dynamic(0);
    omp_set_max_active_levels(1);
    omp_set_num_threads(1);
    Eigen::setNbThreads(1);
    unsigned chainLengthJoint = 500;
    unsigned burninJoint = 100;
    bool printJoint = true;
    MCMC mcmcJoint;
    vector<McmcSamples*> mcmcSampleVec = mcmcJoint.run(jointModel, 1, chainLengthJoint, burninJoint, opt.thin, printJoint, opt.outputFreq, opt.title, false, false);
    
    vector<float> hsqMeanVec(numModels, NAN);
    vector<float> hsqSDVec(numModels, NAN);
    for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
        McmcSamples *mcmcSamples = mcmcSampleVec[i];
        Gadget::Tokenizer token;
        token.getTokens(mcmcSamples->label, "_");
        if (token.size() >= 2 && token.front() == "hsq" && token[1].size() > 1 && token[1][0] == 'M') {
            unsigned mid = (unsigned) stoi(token[1].substr(1)) - 1;
            if (mid < numModels) {
                hsqMeanVec[mid] = mcmcSamples->mean()[0];
                hsqSDVec[mid] = mcmcSamples->sd()[0];
            }
        }
    }
    for (unsigned i=0; i<numModels; ++i) {
        if (std::isnan(hsqMeanVec[i]) || std::isnan(hsqSDVec[i])) {
            throw(" Error: could not read hsq posterior for model M" + to_string(static_cast<long long>(i + 1)) + " in --n-dist-auto-pred.");
        }
    }
    
    // --- (2) Pseudo CV prediction r per model (short MCMC on pseudo training summaries) ---
    VectorXf nGWASblockSaved = data.nGWASblock;
    data.nGWASblock = data.pseudoGwasNtrnBlock;
    
    VectorXf cor(numModels);
    
#pragma omp critical(gctb_ndist_cout)
    {
        cout << "\nPseudo cross-validation uses SBayesR (R) for each candidate K, including when the final analysis is RC." << endl;
        cout << boost::format("%12s %6s %12s %12s %25s\n") % "Model" % "K" % "hsq mean" % "hsq SD" % "Pseudo CV r";
        cout.flush();
    }
    
    for (unsigned i=0; i<numModels; ++i) {
        VectorXf gamma = opt.gamma;
        VectorXf pis = opt.pis;
        VectorXf piPar = opt.piPar;
        for (unsigned j=0; j<i; ++j) {
            Gadget::removeSecondElement(gamma);
            Gadget::removeSecondElement(pis);
            Gadget::removeSecondElement(piPar);
        }
        
#pragma omp critical(gctb_ndist_cout)
        {
            //cout << "  Pseudo CV short MCMC for model M" << (i+1) << " (K=" << gamma.size() << ") ... " << endl;
        }
        
        data.readEigenMatrixBinaryFileAndMakeWandQ(opt.eigenMatrixFile, eigenCutoffUsed, data.pseudoGwasEffectTrn, data.pseudoGwasNtrnBlock, opt.noscale, false);
        data.initVariances(opt.heritability, opt.propVarRandom);
        
        bool print = false;
        // Always ApproxBayesR here: matches joint MultiModelSBayesR and avoids RC cost / threading issues in this inner loop.
        Model *modeli = new ApproxBayesR(data, data.lowRankModel, data.varGenotypic, data.varResidual, pis, piPar, gamma, opt.estimatePi, opt.noscale, opt.hsqPercModel, opt.robustMode, opt.algorithm, print);
        
        MCMC mcmc;
        unsigned chainLength = 150;
        unsigned burnin = 100;
        unsigned thin = 1;
        vector<McmcSamples*> mcmcSampleVeci = mcmc.run(*modeli, 1, chainLength, burnin, thin, print, opt.outputFreq, opt.title, print, print);
        delete modeli;
        
        VectorXf betaMean;
        for (unsigned k=0; k<mcmcSampleVeci.size(); ++k) {
            if (mcmcSampleVeci[k]->label == "SnpEffects") {
                betaMean = mcmcSampleVeci[k]->posteriorMean;
                break;
            }
        }
        if (!betaMean.size()) {
            throw(" Error: SnpEffects posterior not found in MCMC output for --n-dist-auto-pred.");
        }
        
        cor[i] = betaMean.dot(data.b_val) / sqrt(betaMean.squaredNorm() * data.varPhenotypic);
        
        string modelTag = "M" + to_string(i+1);
#pragma omp critical(gctb_ndist_cout)
        {
            cout << boost::format("%12s %6s %12s %12s %25s\n") % modelTag % gamma.size() % hsqMeanVec[i] % hsqSDVec[i] % cor[i];
            cout.flush();
        }
    }
    
    // --- (3) Sequential simplification: start at largest K (model index 0), drop to K-1 unless both hsq and r justify keeping K ---
    unsigned selectedIdx = 0;
    while (selectedIdx < numModels - 1) {
        unsigned simpler = selectedIdx + 1;
        bool hsqClearlyHigher = (hsqMeanVec[selectedIdx] - hsqSDVec[selectedIdx] > hsqMeanVec[simpler]);
        bool predClearlyHigher = (cor[selectedIdx] > cor[simpler] * predRelMin);
        if (hsqClearlyHigher && predClearlyHigher)
            break;
#pragma omp critical(gctb_ndist_cout)
        {
            cout << "  -> K=" << (opt.numDist - selectedIdx) << " does not beat K=" << (opt.numDist - simpler)
                 << " on both hsq (significant) and r (>" << predRelMin << " x simpler); trying simpler model." << endl;
        }
        ++selectedIdx;
    }
    
    VectorXf gammaSel = opt.gamma;
    VectorXf pisSel = opt.pis;
    VectorXf piParSel = opt.piPar;
    for (unsigned j=0; j<selectedIdx; ++j) {
        Gadget::removeSecondElement(gammaSel);
        Gadget::removeSecondElement(pisSel);
        Gadget::removeSecondElement(piParSel);
    }
    
    opt.numDist = gammaSel.size();
    opt.gamma = gammaSel;
    opt.pis = pisSel;
    opt.piPar = piParSel;
    
    data.nGWASblock = nGWASblockSaved;
    data.readEigenMatrixBinaryFileAndMakeWandQ(opt.eigenMatrixFile, eigenCutoffUsed, data.gwasEffectInBlock, data.nGWASblock, opt.noscale, false);
    
    timer.getTime();
    cout << "\nModel " << selectedIdx+1 << " (" << opt.numDist << "-component model) is selected (hsq + pseudo CV r vs simpler model; time used: " << timer.format(timer.getElapse()) << ")." << endl;
}

vector<McmcSamples*> GCTB::multi_chain_mcmc(Data &data, const string &bayesType, const unsigned windowWidth, const float heritability, const float propVarRandom, const float pi, const float piAlpha, const float piBeta, const bool estimatePi, const VectorXf &pis, const VectorXf &gamma, const float phi, const float kappa, const string &algorithm, const unsigned snpFittedPerWindow, const float varS, const vector<float> &S, const float overdispersion, const bool estimatePS, const float icrsq, const float spouseCorrelation, const bool diagnosticMode, const bool robustMode, const unsigned numChains, const unsigned chainLength, const unsigned burnin, const unsigned thin, const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior){
    
    data.initVariances(heritability, propVarRandom);

    vector<Model*> modelVec(numChains);
    
    for (unsigned i=0; i<numChains; ++i) {
        if (data.numAnnos) {
            if (bayesType == "S")
                modelVec[i] = new StratApproxBayesS(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, varS, S, algorithm, robustMode, false, !i);
            else
                throw(" Error: " + bayesType + " is not available in the multi-chain annotation-stratified Bayesian analysis.");
        }
        else {
            if (bayesType == "C")
                modelVec[i] = new ApproxBayesC(data, data.lowRankModel, data.varGenotypic, data.varResidual, data.varRandom, pi, piAlpha, piBeta, estimatePi, false, robustMode, !i);
            else if (bayesType == "S")
                modelVec[i] = new ApproxBayesS(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, varS, S, algorithm, false, !i);
            else if (bayesType == "ST")
                modelVec[i] = new ApproxBayesST(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, varS, S, true, false, !i);
            else if (bayesType == "T")
                modelVec[i] = new ApproxBayesST(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, varS, S, false, false, !i);
            else
                throw(" Error: " + bayesType + " is not currently available in the multi-chain Bayesian analysis.");
        }
    }
    
    vector<vector<McmcSamples*> > mcmcSampleVecChain;
    mcmcSampleVecChain.resize(numChains);
    
    cout << numChains << "-chain ";

//#pragma omp parallel for
    for (unsigned i=0; i<numChains; ++i) {
        MCMC mcmc;
        bool print = true;
        mcmcSampleVecChain[i] = mcmc.run(*modelVec[i], 1, chainLength, burnin, thin, print, outputFreq, title, (writeBinPosterior && print), (writeTxtPosterior && print));
    }
    
    if (numChains) {
        MCMC mcmc;
        mcmc.convergeDiagGelmanRubin(*modelVec[0], mcmcSampleVecChain, title);
    }
    
    return mcmcSampleVecChain[0];
}

void GCTB::saveMcmcSamples(const vector<McmcSamples*> &mcmcSampleVec, const string &filename){
    for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
        McmcSamples *mcmcSamples = mcmcSampleVec[i];
        if (mcmcSamples->label == "SnpEffects" )  continue;
        if (mcmcSamples->label == "WindowDelta") continue;
        mcmcSamples->writeDataTxt(filename);
    }
}

void GCTB::outputResults(Data &data, const vector<McmcSamples*> &mcmcSampleVec, const string &bayesType, const bool noscale, const string &filename){
    
    string unconvergedSnpFile = filename + ".skepticalSNPs";
    ifstream in(unconvergedSnpFile.c_str());
    if (in) data.readUnconvergedSnplist(unconvergedSnpFile);
        
    // For APP (bivariate), handle output separately
    if (bayesType == "APP") {
        McmcSamples *snpEff1 = NULL, *snpEff2 = NULL, *pip1 = NULL, *pip2 = NULL;
        vector<McmcSamples*> appPar;
        for (unsigned i = 0; i < mcmcSampleVec.size(); ++i) {
            if (mcmcSampleVec[i]->label == "SnpEffects")  snpEff1 = mcmcSampleVec[i];
            else if (mcmcSampleVec[i]->label == "SnpEffects2") snpEff2 = mcmcSampleVec[i];
            else if (mcmcSampleVec[i]->label == "PIP")    pip1 = mcmcSampleVec[i];
            else if (mcmcSampleVec[i]->label == "PIP2")   pip2 = mcmcSampleVec[i];
            else appPar.push_back(mcmcSampleVec[i]);
        }
        if (snpEff1 && pip1 && snpEff2 && pip2)
            data.outputBivariateSnpResults(snpEff1->posteriorMean, snpEff1->posteriorSqrMean, pip1->posteriorMean,
                                           snpEff2->posteriorMean, snpEff2->posteriorSqrMean, pip2->posteriorMean,
                                           noscale, filename + ".snpRes");
        else if (snpEff1 && pip1)
            data.outputSnpResults(snpEff1->posteriorMean, snpEff1->posteriorSqrMean, pip1->posteriorMean, noscale, filename + ".snpRes");
        for (McmcSamples *s : appPar) s->writeDataTxt(filename);
        return;
    }

    vector<McmcSamples*> mcmcSamplesPar;
    for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
        McmcSamples *mcmcSamples = mcmcSampleVec[i];
        if (mcmcSamples->label == "SnpEffects") {
            McmcSamples *pip = NULL;
            for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
                if (mcmcSampleVec[i]->label == "PIP") {
                    pip = mcmcSampleVec[i];
                    break;
                }
            }
            if (pip == NULL) {
                throw("Error: couldn't find MCMC samples for PIP!\n");
            }
            data.outputSnpResults(mcmcSamples->posteriorMean, mcmcSamples->posteriorSqrMean, pip->posteriorMean, noscale, filename + ".snpRes");
        }
        else if (mcmcSamples->label == "CovEffects") {
//            if (mcmcSamples->datMat.size()) data.outputFixedEffects(mcmcSamples->datMat, filename + ".covRes");
            data.outputFixedEffects(mcmcSamples->mean(), mcmcSamples->sd(), filename + ".covRes");
        }
        else if (mcmcSamples->label == "RandCovEffects") {
//            data.outputRandomEffects(mcmcSamples->datMat, filename + ".randCovRes");
            data.outputRandomEffects(mcmcSamples->mean(), mcmcSamples->sd(), filename + ".randCovRes");
//        }
//        else if (mcmcSamples->label == "WindowDelta") {
//            //mcmcSamples->readDataBin(mcmcSamples->filename);
//            data.outputWindowResults(mcmcSamples->posteriorMean, filename + ".window");
        } else {
            mcmcSamplesPar.push_back(mcmcSamples);
        }
    }
    
    if (bayesType == "SMix") {
        McmcSamples *snpEffects = NULL;
        McmcSamples *delta = NULL;
        McmcSamples *pip = NULL;
        for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
            if (mcmcSampleVec[i]->label == "SnpEffects") snpEffects = mcmcSampleVec[i];
            if (mcmcSampleVec[i]->label == "DeltaS") delta = mcmcSampleVec[i];
            if (mcmcSampleVec[i]->label == "PIP") pip = mcmcSampleVec[i];
        }
        string newfilename = filename + ".snpRes";
        ofstream out(newfilename.c_str());
        out << boost::format("%6s %20s %6s %12s %8s %12s %12s %12s %8s %8s\n")
        % "Index"
        % "Name"
        % "Chrom"
        % "Position"
        % "GeneFrq"
        % "Effect"
        % "SE"
        % "VarExplained"
        % "PIP"
        % "PiS";
        for (unsigned i=0; i<data.numIncdSnps; ++i) {
            SnpInfo *snp = data.incdSnpInfoVec[i];
            out << boost::format("%6s %20s %6s %12s %8.3f %12.6f %12.6f %12.6e %8.3f %8.3f\n")
            % (i+1)
            % snp->ID
            % snp->chrom
            % snp->physPos
            % snp->af
            % snpEffects->posteriorMean[i]
            % sqrt(snpEffects->posteriorSqrMean[i]-snpEffects->posteriorMean[i]*snpEffects->posteriorMean[i])
            % (2.0*snp->af*(1.0-snp->af)*snpEffects->posteriorSqrMean[i])
            % pip->posteriorMean[i]
            % delta->posteriorMean[i];
        }
        out.close();
    }
    else if (bayesType == "RC" || bayesType == "R" || bayesType == "RD") {
        McmcSamples *snpEffects = NULL;
        McmcSamples *pip = NULL;
        McmcSamples *pep = NULL;
        vector<McmcSamples*> deltaPiVec;
        for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
            if (mcmcSampleVec[i]->label == "SnpEffects") snpEffects = mcmcSampleVec[i];
            if (mcmcSampleVec[i]->label == "PIP") pip = mcmcSampleVec[i];
            if (mcmcSampleVec[i]->label == "PEP") pep = mcmcSampleVec[i];
            if (mcmcSampleVec[i]->label.substr(0, 7) == "DeltaPi") deltaPiVec.push_back(mcmcSampleVec[i]);
        }
        string newfilename = filename + ".snpRes";
        ofstream out(newfilename.c_str());
        out << boost::format("%6s %20s %6s %12s %6s %6s %12s %12s %12s %12s")
        % "Index"
        % "Name"
        % "Chrom"
        % "Position"
        % "A1"
        % "A2"
        % "A1Frq"
        % "A1Effect"
        % "SE"
        % "VarExplained";
        if (pep) out << boost::format(" %12s") % "PEP";
        for (unsigned i=0; i<deltaPiVec.size(); ++i) {
            out << boost::format(" %12s") % deltaPiVec[i]->label.substr(5);
        }
//        out << boost::format(" %14s %14s") % "PIP" % "Pvalue";
        out << boost::format(" %14s") % "PIP";
        if (pip) {
            if(pip->numChains > 1) out << boost::format(" %12s") % "GelmanRubin_R";
        }
        out << endl;
        
        // estimate P value from PIP
        VectorXf pip_vec = 1.0 - deltaPiVec[0]->posteriorMean.array();
        McmcSamples *numSnp1 = NULL;
        for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
            McmcSamples *mcmcSamples = mcmcSampleVec[i];
            if (mcmcSamples->label == "NumSnp1") numSnp1 = mcmcSampleVec[i];
        }
        float propNull = numSnp1->posteriorMean[0]/(float)data.numIncdSnps;
        VectorXf pval(data.numIncdSnps);
        pip2p(data, pip_vec, propNull, pval);
        // END
        
        for (unsigned i=0, idx=0; i<data.numSnps; ++i) {
            SnpInfo *snp = data.snpInfoVec[i];
            if(!data.fullSnpFlag[i]) continue;
            float sqrt2pq = sqrt(2.0*snp->af*(1.0-snp->af));
            if (snp->gwas_scalar) sqrt2pq = snp->gwas_scalar;
            float effect = (snp->flipped ? - snpEffects->posteriorMean[idx] : snpEffects->posteriorMean[idx]);
            float varExp = snpEffects->posteriorSqrMean[idx];
            float se = sqrt(snpEffects->posteriorSqrMean[idx]-snpEffects->posteriorMean[idx]*snpEffects->posteriorMean[idx]);
            if (snp->unconverged) {
                effect = 0.0;
                varExp = 0.0;
                se = 0.0;
            }
            out << boost::format("%6s %20s %6s %12s %6s %6s %12.6f %12.6f %12.6f %12.6e")
            % (i+1)
            % snp->ID
            % snp->chrom
            % snp->physPos
            % (snp->flipped ? snp->a2 : snp->a1)
            % (snp->flipped ? snp->a1 : snp->a2)
            % (snp->flipped ? 1.0-snp->af : snp->af)
            % (noscale ? effect : effect/sqrt2pq)
            % (noscale ? se : se/sqrt2pq)
            % (noscale ? sqrt2pq*sqrt2pq*varExp : varExp);
            if (snp->unconverged) {
                if (pep) {
                    out << boost::format(" %12.6f") % 0;
                }
                for (unsigned j = 0; j < deltaPiVec.size(); ++j) {
                    if (j==0) out << boost::format(" %12.0f") % 1;
                    else out << boost::format(" %12.0f") % 0;
                }
                out << " " << setw(14) << 0;
                if (pip) {
                    if (pip->numChains > 1) out << " " << setw(12) << 0;
                }
            } else {
                if (pep) {
                    out << boost::format(" %12.6f") % pep->posteriorMean[idx];
                }
                for (unsigned j = 0; j < deltaPiVec.size(); ++j) {
                    out << boost::format(" %12.6f") % deltaPiVec[j]->posteriorMean[idx];
                }
                //out << " " << setw(14) << (pip_vec[idx]) << " " << setw(14) << pval[idx];
                out << " " << setw(14) << (pip_vec[idx]);
                if (pip) {
                    if (pip->numChains > 1) out << " " << setw(12) << pip->GelmanRubinStat[idx];
                }
            }
            out << endl;
            ++idx;
        }
        out.close();
    }
}

McmcSamples* GCTB::inputMcmcSamples(const string &mcmcSampleFile, const string &label, const string &fileformat){
    cout << "Reading MCMC samples for " << label << endl;
    McmcSamples *mcmcSamples = new McmcSamples(label);
//    if (fileformat == "bin") mcmcSamples->readDataBin(mcmcSampleFile + "." + label);
    if (fileformat == "bin") mcmcSamples->readDataBin(mcmcSampleFile);
//    if (fileformat == "txt") mcmcSamples->readDataTxt(mcmcSampleFile + "." + label);
    if (fileformat == "txt") mcmcSamples->readDataTxt(mcmcSampleFile + ".Par", label);
    return mcmcSamples;
}

void GCTB::estimateHsq(const Data &data, const McmcSamples &snpEffects, const McmcSamples &resVar, const string &filename, const unsigned outputFreq){
    Heritability hsq(snpEffects.nrow);
    //float phenVar = Gadget::calcVariance(data.y);
    hsq.getEstimate(data, snpEffects, resVar, outputFreq);
    hsq.writeRes(filename);
    hsq.writeMcmcSamples(filename);
}

void GCTB::estimatePi(const Data &data, const McmcSamples &snpEffects, const McmcSamples &genVar, const string &filename, const unsigned outputFreq){
    Polygenicity pi(snpEffects.nrow);
    pi.getEstimate(data, snpEffects, genVar, outputFreq);
    pi.writeRes(filename);
    //pi.writeMcmcSamples(filename);
}

void GCTB::predict(const Data &data, const string &filename){
    Predict pred;
    pred.getAccuracy(data, filename + ".predRes");
    pred.writeRes(data, filename + ".ghat");
}

void GCTB::getWindowPIP(Data &data, McmcSamples &snpEffects, const string &snpResFile, const int windowWidth, const int stepSize, const string &title){
    string filename1 = title + ".snpEffectSamples.txt";
    string filename2 = title + ".windowPIP";
    string filename3 = title + ".credibleSet";
    ofstream out1(filename1.c_str());
    ofstream out2(filename2.c_str());
    ofstream out3(filename3.c_str());

    data.inputNewSnpResults(snpResFile);
    data.getOverlapWindows(windowWidth, stepSize);
    
    MatrixXf windowDelta;
    windowDelta.setZero(snpEffects.nrow, data.numWindows);
        
    //cout << snpEffects.datMatSp.nonZeros() << endl;
    
    VectorXf snpPip;
    snpPip.setZero(data.numSnps);
    
    out1 << boost::format("%8s %8s %12s\n") % "Iter" % "SnpIdx" % "Effect";
    
    for (int k=0; k<snpEffects.datMatSp.outerSize(); ++k) {
        for (SpMat::InnerIterator it(snpEffects.datMatSp,k); it; ++it) {
            //it.value();
            unsigned iter = it.row();   // row index
            unsigned snpIdx = it.col();   // col index (here it is equal to k)
            unsigned winIdx = data.snpInfoVec[snpIdx]->window;
            //cout << "iter " << iter << " snpIdx " << snpIdx << " winIdx " << winIdx << " value " << it.value() << endl;
            windowDelta(iter, winIdx) = 1;
            snpPip[snpIdx]++;
            out1 << boost::format("%8s %8s %12s\n") % iter % snpIdx % it.value();
        }
    }
    
    snpPip /= snpEffects.nrow;
        
    VectorXf windowPip = windowDelta.colwise().mean();
    
    //cout << "windowPip " << endl << windowPip.block(0,10,0,10) << endl;
    //cout << windowPip.head(10).transpose() << endl << endl;
    //cout << snpPip.head(10).transpose() << endl;

    out2 << boost::format("%6s %12s %12s %8s %8s\n")
    % "Index"
    % "Start"
    % "End"
    % "Size"
    % "PIP";
    for (unsigned i=0; i<data.numWindows; ++i) {
        out2 << boost::format("%6s %12s %12s %8s %8.6f\n")
        % (i+1)
        % data.snpInfoVec[data.windStart[i]]->ID
        % data.snpInfoVec[data.windStart[i] + data.windSize[i] - 1]->ID
        % data.windSize[i]
        % windowPip[i];
    }
    out2.close();

    out3 << boost::format("%6s %12s %12s\n") % "Window" % "90_Credible_Set" % "SNP_PIP";

    map<int, vector<SnpInfo*> > credibleSet;
    for (unsigned i=0; i<data.numWindows; ++i) {
        if (windowPip[i] > 0.9) {
            VectorXf snpPipWin = snpPip.segment(data.windStart[i], data.windSize[i]);
            std::sort(snpPipWin.data(), snpPipWin.data() + snpPipWin.size(), std::greater<float>());
            float cumPip = 0.0;
            for (unsigned j=0; j<data.windSize[i]; ++j){
                cumPip += snpPipWin[j];
                unsigned snpIdx = data.windStart[i] + j;
                credibleSet[i].push_back(data.snpInfoVec[snpIdx]);
                out3 << boost::format("%6s %12s %12s\n")
                % (i+1)
                % data.snpInfoVec[snpIdx]->ID
                % snpPipWin[j];
                //cout << i << " " << snpIdx << " " << snpPipWin[j] << " " << data.windSize[i] << endl;
                if (cumPip > 0.9) {
                    break;
                }
            }
        }
    }
    out3.close();
}


void GCTB::calcCredibleSets(Data &data, const string &snpResFile, McmcSamples &snpEffects, const float pipThreshold, const float pepThreshold, const int windowWidth, const string &title){
    string alphaStr = to_string(int(pipThreshold*100));
    string windowWidthStr = to_string(int(windowWidth/1000));
    
    string filename1 = title + "." + windowWidthStr + "kb_" + alphaStr + "_CS.txt";
    string filename2 = title + "." + windowWidthStr + "kb_" + alphaStr + "_CS_summary.txt";
    string filename3 = title + ".genomewide_" + alphaStr + "_CS.txt";
    string filename4 = title + ".genomewide_CS_summary.txt";
    string filename5 = title + "." + windowWidthStr + "kb_WPEP.txt";
    ofstream out1(filename1.c_str());
    ofstream out2(filename2.c_str());
    ofstream out3(filename3.c_str());
    ofstream out4(filename4.c_str());
    ofstream out5(filename5.c_str());

    data.inputNewSnpResults(snpResFile);
    data.getNonoverlapWindowInfo(windowWidth);
    
    string unconvergedSnpFile = title + ".skepticalSNPs";
    ifstream in(unconvergedSnpFile.c_str());
    if (in) data.readUnconvergedSnplist(unconvergedSnpFile);
    for (unsigned j=0; j<data.numSnps; ++j) {
        SnpInfo *snpj = data.snpInfoVec[j];
        if (snpj->unconverged) snpEffects.datMatSp.col(j) *= 0;
    }
    
    // get total genetic variance over MCMC iterations
    VectorXf totalVar;
    totalVar.setZero(snpEffects.nrow);

    for (int k=0; k<snpEffects.datMatSp.outerSize(); ++k) {
        for (SpMat::InnerIterator it(snpEffects.datMatSp,k); it; ++it) {
            totalVar(it.row()) += it.value() * it.value();
        }
    }
    
    // calculate variance explained by each SNP
    for (unsigned j=0; j<data.numSnps; ++j) {
        SnpInfo *snpj = data.snpInfoVec[j];
        VectorXf betaj = snpEffects.datMatSp.col(j);
        snpj->varExplained = (betaj.array().square()/totalVar.array()).mean();  // per-SNP variance explained is the mean of MCMC samples of variance explained
    }
    
    // a parsimonous way of calculating per-SNP variance explained using posterior mean of SNP effects.
    // But this is wrong! The calculation is only for comparison purpose.
    double vg = 0;
    for (unsigned j=0; j<data.numSnps; ++j){
        SnpInfo *snp = data.snpInfoVec[j];
        vg += 2.0*snp->af*(1.0-snp->af)*snp->effect*snp->effect;
    }
    
    out5 << boost::format("%12s %12s %12s %12s %12s\n")
    % "Window"
    % "Size"
    % "PVE"
    % "PVE_Enrich"
    % "PVE_Enrich_PP";

    // calculate per-window variance explained and the probability of per-window heritability enrichment
    VectorXf windowVar;
    VectorXf windowVarImproper;
    windowVar.setZero(data.numWindows);
    windowVarImproper.setZero(data.numWindows);
    vector<WindowInfo*> windowInfoVec;
    windowInfoVec.resize(data.numWindows);
    for (unsigned i=0; i<data.numWindows; ++i) {
        vector<SnpInfo*> snpveci(data.windSize[i]);
        VectorXf windowVarMcmc;
        windowVarMcmc.setZero(snpEffects.nrow);
        for (unsigned j=0; j<data.windSize[i]; ++j){
            SnpInfo *snp = data.snpInfoVec[data.windStart[i]+j];
            //if (i==952) cout << "SNP " << j << " " << snp->index << " " << snp->ID << " " << data.windSize[i] << " " << data.numSnps << endl;
            snpveci[j] = snp;
            windowVar[i] += snp->varExplained;
            windowVarImproper[i] += 2.0*snp->af*(1.0-snp->af)*snp->effect*snp->effect;

            VectorXf betaj = snpEffects.datMatSp.col(snp->index-1);
            VectorXf varj = betaj.array().square();
            windowVarMcmc += varj;
        }
        windowVarMcmc = windowVarMcmc.array()/totalVar.array();
        WindowInfo *window = new WindowInfo(i+1, snpveci);
        window->propGenVarMcmc = windowVarMcmc;
        window->calcVarEnrichPP(float(data.numWindows));
        windowInfoVec[i] = window;
        
        out5 << boost::format("%12s %12s %12.6f %12.6f %12.6f\n")
        % window->index
        % window->size
        % window->propGenVar
        % window->genVarEnrich
        % window->genVarEnrichPP;
    }
    windowVarImproper /= vg;
    
    out5.close();
    
    //cout << "windowVarEnrichPP\n" << windowVarEnrichPP << endl;
        
    // Calculate credible sets per window.
    // first select individual SNPs with PIP > pipThreshold (1-SNP CS). Can be as many as 1-SNP CS per window.
    // then find secondary CS conditional on the 1-SNP CS. Max 1 secondary CS per window, and max 10 SNPs per CS.
    map<unsigned, vector<CredibleSetInfo*> > winCSmap;
    vector<CredibleSetInfo*> CSvec;
    unsigned numSingleSnpCS = 0;
    unsigned numCS = 0;
    unsigned sumCSsize = 0;
    for (unsigned i=0; i<data.numWindows; ++i) {
        //cout << "i " << i << " " << data.windSize[i] << " " << data.windStart[i] << endl;
        vector<SnpInfo*> snpveci(data.windSize[i]);
        for (unsigned j=0; j<data.windSize[i]; ++j) {
            snpveci[j] = data.snpInfoVec[data.windStart[i]+j];
        }
        //cout << "snpveci.size " << snpveci.size() << endl;
        std::sort(snpveci.begin(), snpveci.end(), &GCTB::comparePIP);
        // find out individual SNPs with PIP > pipThreshold
        WindowInfo *window = windowInfoVec[i];
        unsigned sumCSsnpWindowi = 0;
        double cumPip = 0.0;
        double propVar = 0.0;
        vector<SnpInfo*> topSnps;
        unsigned topSnpSize = 0;
        for (unsigned j=0; j<data.windSize[i]; ++j){
            SnpInfo *snp = snpveci[j];
            
//            if (snp->ID == "rs6952746") cout << "rs6952746: window_index " << window->index << " window_size " << window->size << " window_var " << window->propGenVar << " window_var_enrich " << window->genVarEnrich << " window_var_enrich_PP " << window->genVarEnrichPP << endl;
//            if (snp->ID == "rs7374952") cout << "rs7374952: window_index " << window->index << " window_size " << window->size << " window_var " << window->propGenVar << " window_var_enrich " << window->genVarEnrich << " window_var_enrich_PP " << window->genVarEnrichPP << endl;
//            if (snp->ID == "rs38304") cout << "rs38304: window_index " << window->index << " window_size " << window->size << " window_var " << window->propGenVar << " window_var_enrich " << window->genVarEnrich << " window_var_enrich_PP " << window->genVarEnrichPP << endl;
            //if (window->index == 11411) cout << snp->ID << endl;
            
            VectorXf betaj = snpEffects.datMatSp.col(snp->index-1);
            VectorXf varj = betaj.array().square()/totalVar.array();
            window->propGenVarMcmc -= varj;
                        
            if (snp->pip > pipThreshold) {  // single-SNP CS
                vector<SnpInfo*> singleSnp;
                singleSnp.push_back(snp);
                CredibleSetInfo *cs = new CredibleSetInfo(++numCS, pipThreshold, snp->pip, snp->varExplained, singleSnp);
                winCSmap[i].push_back(cs);
                CSvec.push_back(cs);
                ++numSingleSnpCS;
                ++sumCSsize;
                
                cs->windSize = data.windSize[i] - sumCSsnpWindowi;
                cs->windPropGenVar = window->propGenVar;
                cs->windGenVarEnrich = window->genVarEnrich;
                cs->windGenVarEnrichPP = window->genVarEnrichPP;
                window->calcVarEnrichPP(float(data.numWindows));
                ++sumCSsnpWindowi;

            } else { // find multi-SNP credible sets using an iterative approach that tests if the remaining window is still enriched in variance explained
                if (window->genVarEnrichPP > pepThreshold) {
                    cumPip += snpveci[j]->pip;
                    propVar += snpveci[j]->varExplained;
                    topSnps.push_back(snpveci[j]);
                    ++topSnpSize;
                    
                    if (cumPip > pipThreshold) { // && topSnpSize <= 5) {  // secondary multi-SNP CS
                        CredibleSetInfo *cs = new CredibleSetInfo(++numCS, pipThreshold, cumPip, propVar, topSnps);
                        winCSmap[i].push_back(cs);
                        CSvec.push_back(cs);
                        sumCSsize += cs->size;
                        
                        cs->windSize = data.windSize[i] - sumCSsnpWindowi;
                        cs->windPropGenVar = window->propGenVar;
                        cs->windGenVarEnrich = window->genVarEnrich;
                        cs->windGenVarEnrichPP = window->genVarEnrichPP;
                        window->calcVarEnrichPP(float(data.numWindows));
                        sumCSsnpWindowi += cs->size;
                        
                        cumPip = 0.0;
                        propVar = 0.0;
                        topSnps.resize(0);
                        topSnpSize = 0;
                        
                        break;
                    }
                    //else if (topSnpSize > 5) break;
                } else break;
            }
        }
    }
    
    out1 << boost::format("%12s %12s %12s %12s %12s %16s %12s %12s %12s %12s %12s\n")
    % "Window"
    % "Size"
    % "PVE"
    % "PVE_Faulty"
    % "PVE_Enrich"
    % "PVE_Enrich_PP"
    % "CS_Index"
    % "CS_Size"
    % "CS_PIP"
    % "CS_PVE"
    % "CS_SNPs";
    
    map<unsigned, vector<CredibleSetInfo*> >::iterator it, end = winCSmap.end();
    for (it=winCSmap.begin(); it!=end; ++it) {
        unsigned winIdx = it->first;
        for (unsigned i=0; i<it->second.size(); ++i) {
            CredibleSetInfo *cs = it->second[i];
            out1 << boost::format("%12s %12s %12.6f %12.6f %12.4f %16.4f %12s %12s %12.6f %12.6f ")
            % (winIdx + 1)
            % cs->windSize
            % cs->windPropGenVar
            % windowVarImproper[winIdx]
            % cs->windGenVarEnrich //(windowVar[winIdx]*(data.numIncdSnps/float(data.windSize[winIdx])))
            % cs->windGenVarEnrichPP
            % cs->index
            % cs->size
            % cs->sumPIP
            % cs->propVar;
            for (unsigned k=0; k<cs->size; ++k) {
                SnpInfo *snpk = cs->snpVec[k];
                if (k==0) out1 << "\t" << snpk->ID;
                else out1 << "," << snpk->ID;
            }
            out1 << endl;
        }
    }
    out1.close();
        
    
    // summarise the window CS results
    // calculate the power, CS size, and prop hsq for all local credible sets
    
    VectorXd snpPip(data.numSnps);
    for (unsigned i=0; i<data.numSnps; ++i) {
        snpPip[i] = data.snpInfoVec[i]->pip;
    }

    double nnz = data.numSnps * snpPip.mean();
    

    std::sort(CSvec.begin(), CSvec.end(), &GCTB::compareCS);
    
    double cumsumPip = 0.0;
    double cumsumPropVar = 0.0;
    for (unsigned i=0; i<CSvec.size(); ++i) {
        CredibleSetInfo *cs = CSvec[i];
        cumsumPip += cs->sumPIP;
        cumsumPropVar += cs->propVar;
    }

    out2 << boost::format("%50s %12s\n") % "PIP threshold: " % pipThreshold;
    out2 << boost::format("%50s %12s\n") % "PEP threshold: " % pepThreshold;
    out2 << boost::format("%50s %12s\n") % "Number of 1-SNP credible sets: " % numSingleSnpCS;
    out2 << boost::format("%50s %12s\n") % "Number of multi-SNP credible sets: " % (numCS-numSingleSnpCS);
    out2 << boost::format("%50s %12s\n") % "Total number of SNPs in credible sets: " % sumCSsize;
    out2 << boost::format("%50s %12.1f\n") % "Average credible set size: " % (sumCSsize/float(numCS));
    out2 << boost::format("%50s %12.1f\n") % "Estimated total number of causal variants: " %nnz;
    out2 << boost::format("%50s %12.4f\n") % "Estimated power: " % (cumsumPip/nnz);
    out2 << boost::format("%50s %12.1f\n") % "Estimated number of identified causal variants: " % cumsumPip;
    out2 << boost::format("%50s %12.4f\n") % "Estimated proportion of variance explained: " % cumsumPropVar;
        
    out2.close();

    
    out3 << boost::format("%12s %12s %12s\n") % "SNP" % "PIP" % "PropVar";

    std::sort(data.snpInfoVec.begin(), data.snpInfoVec.end(), &GCTB::comparePIP);

    double cumPip = 0.0;
    for (unsigned j=0; j<data.numSnps; ++j){
        SnpInfo *snp = data.snpInfoVec[j];
        cumPip += snp->pip;
        out3 << boost::format("%12s %12s %12s\n") % snp->ID % snp->pip % snp->varExplained;
        if (cumPip > pipThreshold*nnz) break;
    }
    out3.close();
    
    
    out4 << boost::format("%8s %12s %12s\n") % "Threshold" % "CS_size" % "Prop_hsq";
    
    VectorXf threshold_vec(12);
    threshold_vec << 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1;
    for (unsigned k=0; k<threshold_vec.size(); ++k) {
        double cumPip = 0.0;
        double propVar = 0.0;
        int cs_size = 0;
        for (unsigned j=0; j<data.numSnps; ++j){
            SnpInfo *snp = data.snpInfoVec[j];
            cumPip += snp->pip;
            propVar += snp->varExplained;
            ++cs_size;
            if (cumPip > threshold_vec[k]*nnz) {
                //cout << threshold_vec[k] << " " << cs_size << " " << propVar << " " << cumPip << " " << nnz << " " << threshold_vec[k]*nnz << " " << snpPip.sum() << " " << snpPip[0] << endl;
                out4 << boost::format("%8s %12s %12.6f\n")
                % threshold_vec[k]
                % cs_size
                % propVar;
                break;
            }
            if (j == (data.numSnps-1)) {
                out4 << boost::format("%8s %12s %12.6f\n")
                % threshold_vec[k]
                % data.numSnps
                % 1.0;
            }
        }
    }
    out4.close();
    
    
    cout << "The estimated total number of causal variants is " << nnz << "." << endl;
    cout << "Identified " << numCS << " credible sets in " << data.numWindows << " " << windowWidthStr << "kb windows (including " << numSingleSnpCS << " single-SNP credible sets)." << endl;
        
    cout << "Output " << windowWidthStr << "kb window WPEP results into [" + filename5 + "]." << endl;
    cout << "Output " << windowWidthStr << "kb window credible set results into [" + filename1 + "]." << endl;
    cout << "Output " << windowWidthStr << "kb window credible set result summary into [" + filename2 + "]." << endl;
    cout << "Output genome-wide credible set result into [" + filename3 + "]." << endl;
    cout << "Output genome-wide credible set result summary into [" + filename4 + "]." << endl;
    
    cout << endl << "Summary:" << endl;
    cout << boost::format("%50s %12s\n") % "PIP threshold: " % pipThreshold;
    cout << boost::format("%50s %12s\n") % "PEP threshold: " % pepThreshold;
    cout << boost::format("%50s %12s\n") % "Number of 1-SNP credible sets: " % numSingleSnpCS;
    cout << boost::format("%50s %12s\n") % "Number of multi-SNP credible sets: " % (numCS-numSingleSnpCS);
    cout << boost::format("%50s %12s\n") % "Total number of SNPs in credible sets: " % sumCSsize;
    cout << boost::format("%50s %12.1f\n") % "Average credible set size: " % (sumCSsize/float(numCS));
    cout << boost::format("%50s %12.1f\n") % "Estimated total number of causal variants: " %nnz;
    cout << boost::format("%50s %12.4f\n") % "Estimated power: " % (cumsumPip/nnz);
    cout << boost::format("%50s %12.1f\n") % "Estimated number of identified causal variants: " % cumsumPip;
    cout << boost::format("%50s %12.4f\n") % "Estimated proportion of variance explained: " % cumsumPropVar;
}


void GCTB::calcCredibleSets(Data &data, McmcSamples &snpEffects, const float pipThreshold, const float pepThreshold, const string &title){
        
    string filename1 = title + ".lcs";
    string filename2 = title + ".lcsRes";
    string filename3 = title + ".gcs";
    string filename4 = title + ".gcsRes";
    ofstream out1(filename1.c_str());
    ofstream out2(filename2.c_str());
    ofstream out3(filename3.c_str());
    ofstream out4(filename4.c_str());
    
    // set the unconverged SNP effects to be zero
    string unconvergedSnpFile = title + ".skepticalSNPs";
    ifstream in(unconvergedSnpFile.c_str());
    if (in) data.readUnconvergedSnplist(unconvergedSnpFile);
    
    //data.filterSnpByGelmanRubinStat(1.2);
    
    for (unsigned j=0; j<data.numSnps; ++j) {
        SnpInfo *snpj = data.snpInfoVec[j];
        if (snpj->unconverged) snpEffects.datMatSp.col(j) *= 0;
    }
    
    // get total genetic variance over MCMC iterations
    VectorXf totalVar;
    totalVar.setZero(snpEffects.nrow);
    
    for (int k=0; k<snpEffects.datMatSp.outerSize(); ++k) {
        for (SpMat::InnerIterator it(snpEffects.datMatSp,k); it; ++it) {
            totalVar(it.row()) += it.value() * it.value();
        }
    }
    
    // calculate variance explained by each SNP
    for (unsigned j=0; j<data.numSnps; ++j) {
        SnpInfo *snpj = data.snpInfoVec[j];
        VectorXf betaj = snpEffects.datMatSp.col(j);
        snpj->varExplained = (betaj.array().square()/totalVar.array()).mean();  // per-SNP variance explained is the mean of MCMC samples of variance explained
    }
    
    // Calculate local credible sets per LD block
    //    unsigned numLDBlocks = data.numLDBlocks;
    //    vector<int> numSnpInRegion(numLDBlocks);
    //
    //    for(int i = 0; i < numLDBlocks;i++){
    //        LDBlockInfo *block = data.ldBlockInfoVec[i];
    //        numSnpInRegion[i] = block->numSnpInBlock;
    //    }
    
    vector<CredibleSetInfo*> csInfoVec;
    
    // sort each SNP's LD friends by their PIP
    map<SnpInfo*, vector<SnpInfo*> >::iterator it, end = data.LDmap.end();
    for (it=data.LDmap.begin(); it!=data.LDmap.end(); ++it) {
        std::sort(it->second.begin(), it->second.end(), &GCTB::comparePIP);
    }
    
    
    // sort SNPs by PIP
    vector<SnpInfo*> snpInfoVecSorted = data.snpInfoVec;
    std::sort(snpInfoVecSorted.begin(), snpInfoVecSorted.end(), &GCTB::comparePIP);
    
//    for(unsigned j=0; j < 20; j++){
//        SnpInfo *snpj = snpInfoVecSorted[j];
//        cout << snpj->index << " " << snpj->pip << endl;
//    }
    
    // construct CS for each SNP
    unsigned numCS = 0;
    unsigned sumCSsize = 0;
    for(unsigned j=0; j < data.numSnps; j++){
        SnpInfo *snpj = snpInfoVecSorted[j];
        float cumPIP = snpj->pip;
        vector<SnpInfo*> cs;
        cs.push_back(snpj);
        vector<SnpInfo*> &LDfriend = data.LDmap[snpj];
        unsigned numLDfriends = LDfriend.size();
        bool csValid = false;
        if (cumPIP > pipThreshold) {
            csValid = true;
        } else {
            for (unsigned k=0; k<numLDfriends; ++k) {
                SnpInfo *snpk = LDfriend[k];
                if (k==j) continue;
                if (!snpk->inCS) {
                    cumPIP += snpk->pip;
                    cs.push_back(snpk);
                }
                if (cumPIP > pipThreshold) {
                    csValid = true;
                    break;
                }
            }
        }
        if (csValid) {
            // calculate posterior probability of SNP-based heritability enrichment
            VectorXf csPGVmcmc;
            csPGVmcmc.setZero(snpEffects.nrow);
            float csPGV = 0.0;
            unsigned csSize = cs.size();
            for (unsigned k=0; k<csSize; ++k) {
                SnpInfo *snpk = cs[k];
                csPGV += snpk->varExplained;
                VectorXf betak = snpEffects.datMatSp.col(snpk->index - 1);
                VectorXf vark = betak.array().square();
                csPGVmcmc += vark;
            }
            csPGVmcmc = csPGVmcmc.array()/totalVar.array();
            unsigned numPositives = 0;
            float average = float(csSize)/data.numSnps;
            unsigned numMcmcSamples = csPGVmcmc.size();
            for (unsigned k=0; k<numMcmcSamples; ++k) {
                if (csPGVmcmc[k] > average) ++numPositives;
            }
            float csPEP = numPositives/float(numMcmcSamples);
            float csPGVenrich = csPGV*data.numSnps/float(csSize);
                        
            if (csPEP > pepThreshold) {
                CredibleSetInfo *csInfo = new CredibleSetInfo(++numCS, pipThreshold, cumPIP, csPGV, cs);
                csInfo->windGenVarEnrich = csPGVenrich;
                csInfo->windGenVarEnrichPP = csPEP;
                csInfoVec.push_back(csInfo);
                sumCSsize += csInfo->size;
                
                for (unsigned k=0; k<csSize; ++k) {
                    SnpInfo *snpk = cs[k];
                    snpk->inCS = true;
                }
            }
        }
        
        //if(!(j%10000)) cout << " Computed credible sets for SNP " << j << " numCS " << numCS << " sumCSsize " << sumCSsize << "\r" << flush;
    }
    
    
    unsigned numCSwithGRres = 0;
    for (unsigned i=0; i<numCS; ++i) {
        CredibleSetInfo *cs = csInfoVec[i];
        cs->getNumUnconvgSNPs(1.2);
        if (cs->numUnconvgSNPs != -1) ++numCSwithGRres;
    }
        
    unsigned numGenes = data.geneInfoVec.size();
    
    out1 << boost::format("%12s %12s %12s %12s %12s %12s ")
    % "CS"
    % "Size"
    % "PIP"
    % "PGV"
    % "PGVenrich"
    % "PEP";
        
    if (numCSwithGRres) out1 << boost::format("%14s ") % "NumUnconvgSNPs";
   out1 << boost::format("%12s ") % "SNP";
    if (numGenes) {
        GeneInfo *gene = data.geneInfoVec[0];
        out1 << boost::format("%16s %16s ")
        % ("ENSGID_" + gene->genomeBuild)
        % ("GeneName_" + gene->genomeBuild);
    }
    out1 << endl;
    
    for (unsigned i=0; i<numCS; ++i) {
        CredibleSetInfo *cs = csInfoVec[i];
        out1 << boost::format("%12s %12s %12.6f %12.6f %12.6f %12.6f ")
        % (i+1)
        % cs->size
        % cs->sumPIP
        % cs->propVar
        % cs->windGenVarEnrich
        % cs->windGenVarEnrichPP;
 
        if (numCSwithGRres) {
            out1 << boost::format("%14s ") % cs->numUnconvgSNPs;
        }

        for (unsigned k=0; k<cs->size; ++k) {
            SnpInfo *snpk = cs->snpVec[k];
            if (k==0) out1 << setw(8) << " " << snpk->ID;
            else out1 << "," << snpk->ID;
        }
        
        if (numGenes) {
            vector<GeneInfo*> geneVec;
            for (unsigned g=0; g<numGenes; ++g) {
                GeneInfo *gene = data.geneInfoVec[g];
                if (gene->containAllSnps(cs->snpVec)) {
                    geneVec.push_back(gene);
                }
            }
            for (unsigned g=0; g<geneVec.size(); ++g) {
                GeneInfo *gene = geneVec[g];
                if (g==0) out1 << setw(8) << " " << gene->ensgid;
                else out1 << "," << gene->ensgid;
            }
            for (unsigned g=0; g<geneVec.size(); ++g) {
                GeneInfo *gene = geneVec[g];
                if (g==0) out1 << setw(8) << " " << gene->name;
                else out1 << "," << gene->name;
            }
        }
        
        out1 << endl;
    }
    out1.close();
    
    
    // summarise the window CS results
    // calculate the power, CS size, and prop hsq for all local credible sets
    
    VectorXd snpPip(data.numSnps);
    for (unsigned i=0; i<data.numSnps; ++i) {
        snpPip[i] = data.snpInfoVec[i]->pip;
    }
    
    double nnz = data.numSnps * snpPip.mean();
    
    unsigned numCSsingleton = 0;
    double cumsumPip = 0.0;
    double cumsumPropVar = 0.0;
    for (unsigned i=0; i<csInfoVec.size(); ++i) {
        CredibleSetInfo *cs = csInfoVec[i];
        if (cs->size == 1) ++numCSsingleton;
        cumsumPip += cs->sumPIP;
        cumsumPropVar += cs->propVar;
    }
    
    out2 << boost::format("%50s %12s\n") % "PIP threshold: " % pipThreshold;
    out2 << boost::format("%50s %12s\n") % "PEP threshold: " % pepThreshold;
    out2 << boost::format("%50s %12s\n") % "Number of 1-SNP credible sets: " % numCSsingleton;
    out2 << boost::format("%50s %12s\n") % "Number of multi-SNP credible sets: " % (numCS-numCSsingleton);
    out2 << boost::format("%50s %12s\n") % "Total number of SNPs in credible sets: " % sumCSsize;
    out2 << boost::format("%50s %12.1f\n") % "Average credible set size: " % (sumCSsize/float(numCS));
    out2 << boost::format("%50s %12.1f\n") % "Estimated total number of causal variants: " %nnz;
    out2 << boost::format("%50s %12.4f\n") % "Estimated power: " % (cumsumPip/nnz);
    out2 << boost::format("%50s %12.1f\n") % "Estimated number of identified causal variants: " % cumsumPip;
    out2 << boost::format("%50s %12.4f\n") % "Estimated proportion of variance explained: " % cumsumPropVar;
    
    out2.close();
    

    VectorXf threshold_vec(11);
    threshold_vec << 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9;

    out3 << boost::format("%12s %12s %12s\n") % "Alpha" % "SNP" % "PIP";
    
    std::sort(data.snpInfoVec.begin(), data.snpInfoVec.end(), &GCTB::comparePIP);
    
    for (unsigned k=0; k<threshold_vec.size(); ++k) {
        double cumPip = 0.0;
        for (unsigned j=0; j<data.numSnps; ++j){
            SnpInfo *snp = data.snpInfoVec[j];
            cumPip += snp->pip;
            out3 << boost::format("%12.2f %12s %12s\n") % threshold_vec[k]  % snp->ID % snp->pip;
            if (cumPip > threshold_vec[k]*nnz) break;
        }
    }
    out3.close();
    
    
    out4 << boost::format("%8s %12s %12s\n") % "Threshold" % "CS_size" % "Prop_hsq";
    
    threshold_vec.resize(12);
    threshold_vec << 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1;
    for (unsigned k=0; k<threshold_vec.size(); ++k) {
        double cumPip = 0.0;
        double propVar = 0.0;
        int cs_size = 0;
        for (unsigned j=0; j<data.numSnps; ++j){
            SnpInfo *snp = data.snpInfoVec[j];
            cumPip += snp->pip;
            propVar += snp->varExplained;
            ++cs_size;
            if (cumPip > threshold_vec[k]*nnz) {
                //cout << threshold_vec[k] << " " << cs_size << " " << propVar << " " << cumPip << " " << nnz << " " << threshold_vec[k]*nnz << " " << snpPip.sum() << " " << snpPip[0] << endl;
                out4 << boost::format("%8s %12s %12.6f\n")
                % threshold_vec[k]
                % cs_size
                % propVar;
                break;
            }
            if (j == (data.numSnps-1)) {
                out4 << boost::format("%8s %12s %12.6f\n")
                % threshold_vec[k]
                % data.numSnps
                % 1.0;
            }
        }
    }
    out4.close();
    
    
    cout << "The estimated total number of causal variants is " << nnz << "." << endl;
    cout << "Identified " << numCS << " local credible sets (including " << numCSsingleton << " single-SNP local credible sets)." << endl;
    
    cout << "Output local credible set results into [" + filename1 + "]." << endl;
    cout << "Output local credible set result summary into [" + filename2 + "]." << endl;
    cout << "Output global credible set result into [" + filename3 + "]." << endl;
    cout << "Output global credible set result summary into [" + filename4 + "]." << endl;
    
    cout << endl << "Summary:" << endl;
    cout << boost::format("%50s %12s\n") % "PIP threshold: " % pipThreshold;
    cout << boost::format("%50s %12s\n") % "PEP threshold: " % pepThreshold;
    cout << boost::format("%50s %12s\n") % "Number of 1-SNP credible sets: " % numCSsingleton;
    cout << boost::format("%50s %12s\n") % "Number of multi-SNP credible sets: " % (numCS-numCSsingleton);
    cout << boost::format("%50s %12s\n") % "Total number of SNPs in credible sets: " % sumCSsize;
    cout << boost::format("%50s %12.1f\n") % "Average credible set size: " % (sumCSsize/float(numCS));
    cout << boost::format("%50s %12.1f\n") % "Estimated total number of causal variants: " %nnz;
    cout << boost::format("%50s %12.4f\n") % "Estimated power: " % (cumsumPip/nnz);
    cout << boost::format("%50s %12.1f\n") % "Estimated number of identified causal variants: " % cumsumPip;
    cout << boost::format("%50s %12.4f\n") % "Estimated proportion of variance explained: " % cumsumPropVar;
}


void GCTB::clearGenotypes(Data &data){
    data.X.resize(0,0);
}

void GCTB::stratify(Data &data, const string &ldmatrixFile, const bool multiLDmat, const string &geneticMapFile, const float genMapN, const string &snpResFile, const string &mcmcSampleFile, const string &annotationFile, const bool transpose, const string &continuousAnnoFile, const unsigned flank, const string &eQTLFile, const string &gwasSummaryFile, const float pValueThreshold, const bool imputeN, const string &filename, const string &bayesType, unsigned chainLength, unsigned burnin, const unsigned thin, const unsigned outputFreq){
    if (multiLDmat)
        data.readMultiLDmatInfoFile(ldmatrixFile);
    else
        data.readLDmatrixInfoFile(ldmatrixFile + ".info");
    data.inputSnpInfoAndResults(snpResFile, bayesType);
    if (!annotationFile.empty())
        data.readAnnotationFile(annotationFile, transpose, true);
    else
        data.readAnnotationFileFormat2(continuousAnnoFile, flank*1000, eQTLFile);
    data.readGwasSummaryFile(gwasSummaryFile, 1, 0, 0, pValueThreshold, imputeN, true);
    data.includeMatchedSnp();
    if (geneticMapFile.empty()) {
        if (multiLDmat)
            data.readMultiLDmatBinFile(ldmatrixFile);
        else
            data.readLDmatrixBinFile(ldmatrixFile + ".bin");
    } else {
        if (multiLDmat)
            data.readMultiLDmatBinFileAndShrink(ldmatrixFile, genMapN);
        else
            data.readLDmatrixBinFileAndShrink(ldmatrixFile + ".bin");
    }
    data.buildSparseMME(false, true);
    data.makeAnnowiseSparseLDM(data.ZPZsp, data.annoInfoVec, data.snpInfoVec);
    
    McmcSamples *snpEffects = inputMcmcSamples(mcmcSampleFile, "SnpEffects", "bin");
    McmcSamples *hsq = inputMcmcSamples(mcmcSampleFile, "hsq", "txt");

    Model *model;
    
    if (bayesType == "S") {
        model = new PostHocStratifyS(data, false, *snpEffects, *hsq, thin, hsq->mean()[0]);
    }
    else if (bayesType == "SMix") {
        McmcSamples *deltaS = inputMcmcSamples(mcmcSampleFile, "DeltaS", "bin");
        model = new PostHocStratifySMix(data, false, *snpEffects, *hsq, *deltaS, thin, hsq->mean()[0]);
    }
    else
        throw(" Error: Wrong bayes type: " + bayesType + " in the annotation-stratified summary-data-based Bayesian analysis.");

    
    if (chainLength > snpEffects->nrow) chainLength = snpEffects->nrow;
    if (burnin > chainLength) burnin = 0.2*chainLength;
    
    runMcmc(*model, 1, chainLength, burnin, thin, outputFreq, filename, false, false);
}

void GCTB::solveSnpEffectsByConjugateGradientMethod(Data &data, const float lambda, const string &filename) const {
    cout << "\nSolving SNP effects by conjugate gradient method ..." << endl;
    cout << "  Lambda = " << lambda << endl;
    //SpMat L(data.numIncdSnps, data.numIncdSnps);
    //SpMat C(data.numIncdSnps, data.numIncdSnps);
    SpMat C(data.numIncdSnps, data.numIncdSnps);
    //L.reserve(data.windSize);
    
    vector<Triplet<float> > tripletList;
    tripletList.reserve(data.windSize.cast<double>().sum());
    
    float val = 0.0;
    for (unsigned i=0; i<data.numIncdSnps; ++i) {
        SnpInfo *snpi = data.incdSnpInfoVec[i];
        if (!(i % 100000)) cout << "  making sparse LD matrix for SNP " << i << " " << data.windSize[i] << " " << snpi->windSize << " " << data.ZPZsp[i].size() << " " << data.ZPZsp[i].nonZeros() << endl;
        for (SparseVector<float>::InnerIterator it(data.ZPZsp[i]); it; ++it) {
            //if (it.index() > i) break;
            //L.insert(i, it.index()) = it.value();
            val = it.value();
            if (it.index() == i) val += lambda;  // adding lambda to diagonals
            tripletList.push_back(Triplet<float>(i, it.index(), val));
        }
    }
    //C = L.transpose().triangularView<Upper>();
    C.setFromTriplets(tripletList.begin(), tripletList.end());
    C.makeCompressed();
    tripletList.clear();
    
    cout << "Running conjugate gradient algorithm ..." << endl;
    
    Gadget::Timer timer;
    timer.setTime();

    ConjugateGradient<SparseMatrix<float, Eigen::ColMajor, long long>, Lower|Upper> cg;
    
    cout << "  preconditioning ..." << endl;
    
    cg.compute(C);
    
    cout << "  solving ..." << endl;
    
    VectorXf sol(data.numIncdSnps);
    sol = cg.solve(data.ZPy);
    
    timer.getTime();
    
    cout << "#iterations:     " << cg.iterations() << endl;
    cout << "estimated error: " << cg.error()      << endl;
    cout << "time used:       " << timer.format(timer.getElapse()) << endl;
    
    ofstream out(filename.c_str());
    out << boost::format("%6s %20s %6s %12s %6s %6s %12s %12s\n")
    % "Index"
    % "Name"
    % "Chrom"
    % "Position"
    % "A1"
    % "A2"
    % "A1Frq"
    % "A1Sol";
    for (unsigned i=0, idx=0; i<data.numSnps; ++i) {
        SnpInfo *snp = data.snpInfoVec[i];
        if(!data.fullSnpFlag[i]) continue;
        out << boost::format("%6s %20s %6s %12s %6s %6s %12.6f %12.6f\n")
        % (idx+1)
        % snp->ID
        % snp->chrom
        % snp->physPos
        % (snp->flipped ? snp->a2 : snp->a1)
        % (snp->flipped ? snp->a1 : snp->a2)
        % (snp->flipped ? 1.0-snp->af : snp->af)
        % (snp->flipped ? -sol[idx] : sol[idx]);
        ++idx;
    }
    out.close();
}

void GCTB::pip2p(const Data &data, const VectorXf &pip, const float propNull, VectorXf &pval){
    VectorXf pipSrt = pip;
    std::sort(pipSrt.data(), pipSrt.data() + pipSrt.size(), greater<float>());
    float cumsum = 0.0;
    float pvali = 0.0;
    float numNull = data.numIncdSnps*propNull;
    map<float, float> pip2pMap;
    for (unsigned i=0; i<data.numIncdSnps; ++i){
        cumsum += pipSrt[i];
        pvali = (i+1 - cumsum) / numNull;
        pip2pMap[pipSrt[i]] = pvali;
    }
    for (unsigned i=0; i<data.numIncdSnps; ++i){
        pval[i] = pip2pMap[pip[i]];
    }
}

float GCTB::tuneEigenCutoff(Data &data, const Options &opt){
    cout << "\nFinding the best eigen cutoff from [" << opt.eigenCutoff.transpose() << "] based on pseudo summary data validation." << endl;
    
    Gadget::Timer timer;
    timer.setTime();
    
    VectorXf nGWASblock = data.nGWASblock;
    data.nGWASblock = data.pseudoGwasNtrnBlock;
    
    unsigned size = opt.eigenCutoff.size();
    VectorXf cor(size);
    VectorXf rel(size);
    
    cout << boost::format("%10s %25s %20s\n") % "Cutoff" % "Prediction accuracy (r)" % "Relative accuracy";
    
    for (unsigned i=0; i<size; ++i) {
        float cutoff = opt.eigenCutoff[i];

        data.readEigenMatrixBinaryFileAndMakeWandQ(opt.eigenMatrixFile, cutoff, data.pseudoGwasEffectTrn, data.pseudoGwasNtrnBlock, false, false);
        //data.readEigenMatrixBinaryFile(opt.eigenMatrixFile, cutoff);
        //data.constructWandQ(data.pseudoGwasEffectTrn, data.pseudoGwasNtrn);
        
        data.initVariances(opt.heritability, opt.propVarRandom);
        bool nDistAuto = false;
        bool print = false;
        Model *modeli = new ApproxBayesR(data, data.lowRankModel, data.varGenotypic, data.varResidual, opt.pis, opt.piPar, opt.gamma, opt.estimatePi, opt.noscale, opt.hsqPercModel, opt.robustMode, opt.algorithm, print);
        
        vector<McmcSamples*> mcmcSampleVeci;
        MCMC mcmc;

        unsigned chainLength = 150;
        unsigned burnin = 100;
        unsigned thin = 1;
        mcmcSampleVeci = mcmc.run(*modeli, 1, chainLength, burnin, thin, print, opt.outputFreq, opt.title, print, print);

        VectorXf betaMean;
        for (unsigned i=0; i<mcmcSampleVeci.size(); ++i) {
            McmcSamples *mcmcSamples = mcmcSampleVeci[i];
            if (mcmcSamples->label == "SnpEffects") {
                betaMean = mcmcSamples->posteriorMean;
            }
        }
        
        // compute prediction accuracy
        cor[i] = betaMean.dot(data.b_val) / sqrt(betaMean.squaredNorm() * data.varPhenotypic);
        rel[i] = cor[i]/cor[0];
        
        cout << boost::format("%10s %25s %20s\n") % cutoff % cor[i] % rel[i];

    }
    
    data.nGWASblock = nGWASblock;
    
    int bestCutoff_index;
    cor.maxCoeff(&bestCutoff_index);
    float bestCutoff = opt.eigenCutoff[bestCutoff_index];

    if (cor[0] < 0) {
        if (rel.maxCoeff() > -1.25)
            bestCutoff = opt.eigenCutoff[0];
    } else {
        if (rel.maxCoeff() < 1.25)
            bestCutoff = opt.eigenCutoff[0];
    }
    
    timer.getTime();

    if (bestCutoff == opt.eigenCutoff.minCoeff()) {
        cout << "==============================================" << endl;
        cout << "Warning: the best eigen cutoff is the minimum value in the tuning set. We suggest expand the tuning set by including lower candidate values, e.g. --ldm-eigen-cutoff 0.995,0.9,0.8,0.7,0.6  (time used: " << timer.format(timer.getElapse()) << ")." << endl;
        cout << "==============================================" << endl;
    } else {
        cout << bestCutoff << " is selected to be the eigen cutoff to continue the analysis (time used: " << timer.format(timer.getElapse()) << ")."  << endl;
    }
    
    return bestCutoff;
}


