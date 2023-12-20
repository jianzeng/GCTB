//
//  gctb.cpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "gctb.hpp"

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
                        const bool sampleOverlap, const bool imputeN, const bool noscale, const bool readLDMfromTxtFile, const bool imputeSummary, const unsigned includeBlock){
    data.readEigenMatrix(eigenMatrixFile, eigenCutoff);
    if (!includeSnpFile.empty()) data.includeSnp(includeSnpFile);
    if (!excludeSnpFile.empty()) data.excludeSnp(excludeSnpFile);
    if (includeChr) data.includeChr(includeChr);
    if (includeBlock) data.includeBlock(includeBlock);
    if (excludeAmbiguousSNP) data.excludeAmbiguousSNP();
    if (!excludeRegionFile.empty()) data.excludeRegion(excludeRegionFile);
    if (excludeMHC) data.excludeMHC();
    if (!annotationFile.empty())
        data.readAnnotationFile(annotationFile, transpose, true);
    else if (!continuousAnnoFile.empty())
        data.readAnnotationFileFormat2(continuousAnnoFile, flank*1000, eQTLFile);
    if (!ldscoreFile.empty()) data.readLDscoreFile(ldscoreFile);
    if (!gwasSummaryFile.empty()) {
        bool removeOutlierN = imputeSummary;
        data.readGwasSummaryFile(gwasSummaryFile, afDiff, mafmin, mafmax, pValueThreshold, imputeN, removeOutlierN);
        if (imputeSummary) {
            data.readEigenMatrixBinaryFile(eigenMatrixFile, eigenCutoff);
            data.impG(includeBlock);
            return;
        }
        data.includeMatchedSnp();
    }
    

    /// partition ld into blocks
//    if(!ldBlockInfoFile.empty()) data.readLDBlockInfoFile(ldBlockInfoFile);
        
    if(!gwasSummaryFile.empty()) data.buildMMEeigen(eigenMatrixFile, sampleOverlap, eigenCutoff, noscale);
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

Model* GCTB::buildModel(Data &data, const string &bedFile, const string &gwasFile, const string &bayesType, const unsigned windowWidth,
                        const float heritability, const float propVarRandom, const float pi, const float piAlpha, const float piBeta, const bool estimatePi, const bool noscale,
                        const VectorXf &pis, const VectorXf &piPar, const VectorXf &gamma, const bool estimateSigmaSq,
                        const float phi, const float kappa, const string &algorithm, const unsigned snpFittedPerWindow,
                        const float varS, const vector<float> &S, const float overdispersion, const bool estimatePS,
                        const float icrsq, const float spouseCorrelation, const bool diagnosticMode, const bool originalModel, const bool perSnpGV, const bool robustMode, const bool nDistAuto){
    data.initVariances(heritability, propVarRandom);
//    if (!bedFile.empty()) {   // TMP_JZ
//        unsigned n_gwas = data.numKeptInds;
//        data.readFamFile(bedFile + ".fam");
//        data.numKeptInds = data.numInds;
//        data.readBedFile(noscale, bedFile + ".bed");
//        data.numKeptInds = n_gwas;
//    }
    if (!gwasFile.empty()) {
        if (data.numAnnos) {
            if (bayesType == "S")
                return new StratApproxBayesS(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, S, algorithm, robustMode);
            else if (bayesType == "RC")
                return new ApproxBayesRC(data, data.lowRankModel, data.varGenotypic, data.varResidual, pis, piPar, gamma, estimatePi, estimateSigmaSq, noscale, originalModel, perSnpGV, overdispersion, estimatePS, spouseCorrelation, diagnosticMode, robustMode, algorithm, nDistAuto);
            else
                throw(" Error: Wrong bayes type: " + bayesType + " in the annotation-stratified summary-data-based Bayesian analysis.");
        }
        else {
            if (bayesType == "C")
                return new ApproxBayesC(data, data.lowRankModel, data.varGenotypic, data.varResidual, data.varRandom, pi, piAlpha, piBeta, estimatePi, noscale, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, diagnosticMode, robustMode);
            else if (bayesType == "B")
            return new ApproxBayesB(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, noscale, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, diagnosticMode, robustMode);
            else if (bayesType == "S")
                return new ApproxBayesS(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, S, algorithm, diagnosticMode, robustMode);
            else if (bayesType == "ST")
                return new ApproxBayesST(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, overdispersion, estimatePS, varS, S, true);
            else if (bayesType == "T")
                return new ApproxBayesST(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, overdispersion, estimatePS, varS, S, false);
            else if (bayesType == "SMix")
                return new ApproxBayesSMix(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, overdispersion, estimatePS, varS, S);
            else if (bayesType == "R")
                return new ApproxBayesR(data, data.lowRankModel, data.varGenotypic, data.varResidual, pis, piPar, gamma, estimatePi, estimateSigmaSq, noscale, originalModel, overdispersion, estimatePS, spouseCorrelation, diagnosticMode, robustMode, algorithm);
            else if (bayesType == "Kap")
                return new ApproxBayesKappa(data, data.lowRankModel, data.varGenotypic, data.varResidual, pis, piPar, gamma, estimatePi, noscale, originalModel, icrsq, kappa);
            else if (bayesType == "RS")
                return new ApproxBayesRS(data, data.lowRankModel, data.varGenotypic, data.varResidual, pis, piPar, gamma, estimatePi, varS, S, algorithm, noscale, originalModel, overdispersion, estimatePS, spouseCorrelation, diagnosticMode, robustMode, algorithm);
            else
                throw(" Error: Wrong bayes type: " + bayesType + " in the summary-data-based Bayesian analysis.");
        }
    }
    if (data.numAnnos) {
        if (bayesType == "RC") {
            data.readBedFile(noscale, bedFile + ".bed");
            return new BayesRC(data, data.varGenotypic, data.varResidual, data.varRandom, pis, piPar, gamma, estimatePi, noscale, originalModel, algorithm);
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
        return new BayesR(data, data.varGenotypic, data.varResidual, data.varRandom, pis, piPar, gamma, estimatePi, noscale, originalModel, algorithm);
    }
    else if (bayesType == "S") {
        data.readBedFile(noscale, bedFile + ".bed");
        return new BayesS(data, data.varGenotypic, data.varResidual, data.varRandom, pi, piAlpha, piBeta, estimatePi, varS, S, algorithm);
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
        return new BayesRS(data, data.varGenotypic, data.varResidual, data.varRandom, pis, piPar, gamma, estimatePi, varS, S, noscale, originalModel, algorithm);
    }
    else if (bayesType == "Cap") {
        //data.readBedFile(bedFile + ".bed");
        data.buildSparseMME(bedFile + ".bed", windowWidth);
        return new ApproxBayesC(data, data.lowRankModel, data.varGenotypic, data.varResidual, data.varRandom, pi, piAlpha, piBeta, estimatePi, noscale, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, diagnosticMode, robustMode);
    }
    else if (bayesType == "Sap") {
        data.buildSparseMME(bedFile + ".bed", windowWidth);
        return new ApproxBayesS(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, S, algorithm, diagnosticMode, robustMode);
    }
    else {
        throw(" Error: Wrong bayes type: " + bayesType);
    }
}

vector<McmcSamples*> GCTB::runMcmc(Model &model, const unsigned chainLength, const unsigned burnin, const unsigned thin, const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior){
    MCMC mcmc;
    return mcmc.run(model, chainLength, burnin, thin, true, outputFreq, title, writeBinPosterior, writeTxtPosterior);
}

vector<McmcSamples*> GCTB::multi_chain_mcmc(Data &data, const string &bayesType, const unsigned windowWidth, const float heritability, const float propVarRandom, const float pi, const float piAlpha, const float piBeta, const bool estimatePi, const VectorXf &pis, const VectorXf &gamma, const float phi, const float kappa, const string &algorithm, const unsigned snpFittedPerWindow, const float varS, const vector<float> &S, const float overdispersion, const bool estimatePS, const float icrsq, const float spouseCorrelation, const bool diagnosticMode, const bool robustMode, const unsigned numChains, const unsigned chainLength, const unsigned burnin, const unsigned thin, const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior){
    
    data.initVariances(heritability, propVarRandom);

    vector<Model*> modelVec(numChains);
    
    for (unsigned i=0; i<numChains; ++i) {
        if (data.numAnnos) {
            if (bayesType == "S")
                modelVec[i] = new StratApproxBayesS(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, S, algorithm, robustMode, true, !i);
            else
                throw(" Error: " + bayesType + " is not available in the multi-chain annotation-stratified Bayesian analysis.");
        }
        else {
            if (bayesType == "C")
                modelVec[i] = new ApproxBayesC(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, diagnosticMode, robustMode, true, !i);
            else if (bayesType == "S")
                modelVec[i] = new ApproxBayesS(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, S, algorithm, diagnosticMode, robustMode, true, !i);
            else if (bayesType == "ST")
                modelVec[i] = new ApproxBayesST(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, overdispersion, estimatePS, varS, S, true, true, !i);
            else if (bayesType == "T")
                modelVec[i] = new ApproxBayesST(data, data.lowRankModel, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, overdispersion, estimatePS, varS, S, false, true, !i);
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
        mcmcSampleVecChain[i] = mcmc.run(*modelVec[i], chainLength, burnin, thin, print, outputFreq, title, (writeBinPosterior && print), (writeTxtPosterior && print));
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

void GCTB::outputResults(const Data &data, const vector<McmcSamples*> &mcmcSampleVec, const string &bayesType, const bool noscale, const string &filename){
    vector<McmcSamples*> mcmcSamplesPar;
    for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
        McmcSamples *mcmcSamples = mcmcSampleVec[i];
        if (mcmcSamples->label == "SnpEffects") {
            //mcmcSamples->readDataBin(mcmcSamples->filename);
//            data.outputSnpResults(mcmcSamples->posteriorMean, mcmcSamples->posteriorSqrMean, mcmcSamples->pip, noscale, filename + ".snpRes");
            McmcSamples *pip = NULL;
            for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
                if (mcmcSampleVec[i]->label == "PIP") {
                    pip = mcmcSampleVec[i];
                    break;
                }
            }
            if (pip == NULL) {
                //data.outputSnpResults(mcmcSamples->posteriorMean, mcmcSamples->posteriorSqrMean, mcmcSamples->pip, noscale, filename + ".snpRes");
                data.outputSnpResults(mcmcSamples->posteriorMean, mcmcSamples->posteriorSqrMean, mcmcSamples->lastSample, mcmcSamples->pip, noscale, filename + ".snpRes");
            } else
                data.outputSnpResults(mcmcSamples->posteriorMean, mcmcSamples->posteriorSqrMean, mcmcSamples->lastSample, pip->posteriorMean, noscale, filename + ".snpRes");
        }
        else if (mcmcSamples->label == "CovEffects") {
            if (mcmcSamples->datMat.size()) data.outputFixedEffects(mcmcSamples->datMat, filename + ".covRes");
        }
        else if (mcmcSamples->label == "RandCovEffects") {
            data.outputRandomEffects(mcmcSamples->datMat, filename + ".randCovRes");
        }
        else if (mcmcSamples->label == "WindowDelta") {
            //mcmcSamples->readDataBin(mcmcSamples->filename);
            data.outputWindowResults(mcmcSamples->posteriorMean, filename + ".window");
        } else {
            mcmcSamplesPar.push_back(mcmcSamples);
        }
    }
    if (bayesType == "SMix") {
        McmcSamples *snpEffects = NULL;
        McmcSamples *delta = NULL;
        for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
            if (mcmcSampleVec[i]->label == "SnpEffects") snpEffects = mcmcSampleVec[i];
            if (mcmcSampleVec[i]->label == "DeltaS") delta = mcmcSampleVec[i];
        }
        string newfilename = filename + ".snpRes";
        ofstream out(newfilename.c_str());
        out << boost::format("%6s %20s %6s %12s %8s %12s %12s %8s %8s\n")
        % "Id"
        % "Name"
        % "Chrom"
        % "Position"
        % "GeneFrq"
        % "Effect"
        % "SE"
        % "PIP"
        % "PiS";
        for (unsigned i=0; i<data.numIncdSnps; ++i) {
            SnpInfo *snp = data.incdSnpInfoVec[i];
            out << boost::format("%6s %20s %6s %12s %8.3f %12.6f %12.6f %8.3f %8.3f\n")
            % (i+1)
            % snp->ID
            % snp->chrom
            % snp->physPos
            % snp->af
            % snpEffects->posteriorMean[i]
            % sqrt(snpEffects->posteriorSqrMean[i]-snpEffects->posteriorMean[i]*snpEffects->posteriorMean[i])
            % snpEffects->pip[i]
            % delta->posteriorMean[i];
        }
        out.close();
    }
    if (bayesType == "RC") {
        McmcSamples *snpEffects = NULL;
        vector<McmcSamples*> deltaPiVec;
        for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
            if (mcmcSampleVec[i]->label == "SnpEffects") snpEffects = mcmcSampleVec[i];
            if (mcmcSampleVec[i]->label.substr(0, 7) == "DeltaPi") deltaPiVec.push_back(mcmcSampleVec[i]);
        }
        string newfilename = filename + ".snpRes";
        ofstream out(newfilename.c_str());
        out << boost::format("%6s %20s %6s %12s %6s %6s %12s %12s %12s")
        % "Id"
        % "Name"
        % "Chrom"
        % "Position"
        % "A1"
        % "A2"
        % "A1Frq"
        % "A1Effect"
        % "SE";
        for (unsigned i=0; i<deltaPiVec.size(); ++i) {
            out << boost::format(" %12s") % deltaPiVec[i]->label.substr(5);
        }
        out << boost::format(" %14s %14s") % "PIP" % "Pvalue";
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
            float effect = (snp->flipped ? - snpEffects->posteriorMean[idx] : snpEffects->posteriorMean[idx]);
            float se = sqrt(snpEffects->posteriorSqrMean[idx]-snpEffects->posteriorMean[idx]*snpEffects->posteriorMean[idx]);
            out << boost::format("%6s %20s %6s %12s %6s %6s %12.6f %12.6f %12.6f")
            % (i+1)
            % snp->ID
            % snp->chrom
            % snp->physPos
            % (snp->flipped ? snp->a2 : snp->a1)
            % (snp->flipped ? snp->a1 : snp->a2)
            % (snp->flipped ? 1.0-snp->af : snp->af)
            % (noscale ? effect : effect/sqrt2pq)
            % (noscale ? se : se/sqrt2pq);
            for (unsigned j = 0; j < deltaPiVec.size(); ++j) {
                out << boost::format(" %12.6f") % deltaPiVec[j]->posteriorMean[idx];
            }
            out << " " << setw(14) << (pip_vec[idx]) << " " << setw(14) << pval[idx];
            out << endl;
            ++idx;
        }
        out.close();
    }
}

McmcSamples* GCTB::inputMcmcSamples(const string &mcmcSampleFile, const string &label, const string &fileformat){
    cout << "reading MCMC samples for " << label << endl;
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

    data.inputSnpResultsOnly(snpResFile);
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
    % "Id"
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

void GCTB::calcCredibleSets(Data &data, const string &snpResFile, McmcSamples &snpEffects, const float csThreshold, const int windowWidth, const string &title){
    string alphaStr = to_string(int(csThreshold*100));
    string windowWidthStr = to_string(int(windowWidth/1000));
    
    string filename1 = title + "." + windowWidthStr + "kb_" + alphaStr + "_CS";
    string filename2 = title + ".genomewide_" + alphaStr + "_CS";
    string filename3 = title + ".genomewide_CS_summary";
    ofstream out1(filename1.c_str());
    ofstream out2(filename2.c_str());
    ofstream out3(filename3.c_str());

    data.inputSnpResultsOnly(snpResFile);
    data.getNonoverlapWindowInfo(windowWidth);

    // calculate variance explained by each SNP
    VectorXf totalVar;
    totalVar.setZero(snpEffects.nrow);

    for (int k=0; k<snpEffects.datMatSp.outerSize(); ++k) {
        for (SpMat::InnerIterator it(snpEffects.datMatSp,k); it; ++it) {
            totalVar(it.row()) += it.value() * it.value();
        }
    }
    for (unsigned j=0; j<data.numSnps; ++j) {
        SnpInfo *snpj = data.snpInfoVec[j];
        VectorXf betaj = snpEffects.datMatSp.col(j);
        snpj->varExplained = (betaj.array().square()/totalVar.array()).mean();
    }
    
    float vg = 0;
    for (unsigned j=0; j<data.numSnps; ++j){
        SnpInfo *snp = data.snpInfoVec[j];
        vg += 2.0*snp->af*(1.0-snp->af)*snp->effect*snp->effect;
    }

    VectorXf windowVar;
    VectorXf windowVarImproper;
    windowVar.setZero(data.numWindows);
    windowVarImproper.setZero(data.numWindows);
    for (unsigned i=0; i<data.numWindows; ++i) {
        for (unsigned j=0; j<data.windSize[i]; ++j){
            SnpInfo *snp = data.snpInfoVec[data.windStart[i]+j];
            windowVar[i] += snp->varExplained;
            windowVarImproper[i] += 2.0*snp->af*(1.0-snp->af)*snp->effect*snp->effect;
        }
    }
    windowVarImproper /= vg;
        
    out1 << boost::format("%6s %12s %12s %12s %12s %12s %12s\n")
    % "Window"
    % "TotalNumSnps"
    % "PropVar"
    % "PropVarImproper"
    % (alphaStr + "_CS_size")
    % (alphaStr + "_CS_PropVar")
    % (alphaStr + "_CS_SNPs");
    
    unsigned nCS10CV = 0;
    
    map<int, vector<SnpInfo*> > credibleSet;
    for (unsigned i=0; i<data.numWindows; ++i) {
        //cout << "i " << i << " " << data.windSize[i] << " " << data.windStart[i] << endl;
        vector<SnpInfo*> snpveci(data.windSize[i]);
        for (unsigned j=0; j<data.windSize[i]; ++j) {
            snpveci[j] = data.snpInfoVec[data.windStart[i]+j];
        }
        //cout << "snpveci.size " << snpveci.size() << endl;
        std::sort(snpveci.begin(), snpveci.end(), &GCTB::comparePIP);
        float cumPip = 0.0;
        float propVar = 0.0;
        for (unsigned j=0; j<data.windSize[i]; ++j){
            cumPip += snpveci[j]->pip;
            propVar += snpveci[j]->varExplained;
            credibleSet[i].push_back(snpveci[j]);
            //cout << i << " " << j << " " << snpveci[j]->pip << endl;

            if (cumPip > csThreshold) {
                out1 << boost::format("%6s %12s %12s %12s %12s %12s ")
                % (i+1)
                % data.windSize[i]
                % windowVar[i]
                % windowVarImproper[i]
                % credibleSet[i].size()
                % propVar;
                vector<SnpInfo*>::iterator it, begin = credibleSet[i].begin();
                for (it=begin; it!=credibleSet[i].end(); ++it) {
                    if (it == begin) out1 << "\t" << (*it)->ID;
                    else out1 << "," << (*it)->ID;
                }
                out1 << endl;
                if (credibleSet[i].size() <= 10) ++nCS10CV;
                break;
            }
        }
    }
    out1.close();
    
    out2 << boost::format("%12s %12s %12s\n") % "SNP" % "PIP" % "PropVar";
        
    VectorXf snpPip(data.numSnps);
    for (unsigned i=0; i<data.numSnps; ++i) {
        snpPip[i] = data.snpInfoVec[i]->pip;
    }

    float nnz = data.numSnps * snpPip.mean();
    
    cout << "The estimated total number of causal variants is " << nnz << "." << endl;
    cout << "There are " << nCS10CV << " 100kb windows with credible set size up to 10." << endl;
        
    std::sort(data.snpInfoVec.begin(), data.snpInfoVec.end(), &GCTB::comparePIP);

    float cumPip = 0.0;
    for (unsigned j=0; j<data.numSnps; ++j){
        SnpInfo *snp = data.snpInfoVec[j];
        cumPip += snp->pip;
        out2 << boost::format("%12s %12s %12s\n") % snp->ID % snp->pip % snp->varExplained;
        if (cumPip > csThreshold*nnz) break;
    }
    out2.close();
    
    out3 << boost::format("%8s %12s %12s\n") % "Threshold" % "CS_size" % "Prop_hsq";
        
    VectorXf threshold_vec(12);
    threshold_vec << 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1;
    for (unsigned k=0; k<threshold_vec.size(); ++k) {
        float cumPip = 0.0;
        float propVar = 0.0;
        int cs_size = 0;
        for (unsigned j=0; j<data.numSnps; ++j){
            SnpInfo *snp = data.snpInfoVec[j];
            cumPip += snp->pip;
            propVar += snp->varExplained;
            ++cs_size;
            if (cumPip > threshold_vec[k]*nnz) {
                out3 << boost::format("%8s %12s %12.6f\n")
                % threshold_vec[k]
                % cs_size
                % propVar;
                break;
            }
            if (j == (data.numSnps-1)) {
                out3 << boost::format("%8s %12s %12.6f\n")
                % threshold_vec[k]
                % data.numSnps
                % 1.0;
            }
        }
    }
    
    cout << "Output window credible set results into [" + filename1 + "]." << endl;
    cout << "Output genome-wide credible set result into [" + filename2 + "]." << endl;
    cout << "Output genome-wide credible set result summary into [" + filename3 + "]." << endl;
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
    
    runMcmc(*model, chainLength, burnin, thin, outputFreq, filename, false, false);
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
    % "Id"
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
    
    unsigned numKeptInds = data.numKeptInds;    
    data.numKeptInds = data.pseudoGwasNtrn;
    
    unsigned size = opt.eigenCutoff.size();
    VectorXf cor(size);
    VectorXf rel(size);
    
    cout << boost::format("%10s %25s %20s\n") % "Cutoff" % "Prediction accuracy (r)" % "Relative accuracy";
    
    for (unsigned i=0; i<size; ++i) {
        float cutoff = opt.eigenCutoff[i];

        data.readEigenMatrixBinaryFileAndMakeWandQ(opt.eigenMatrixFile, cutoff, data.pseudoGwasEffectTrn, data.pseudoGwasNtrn, false, false);
        //data.readEigenMatrixBinaryFile(opt.eigenMatrixFile, cutoff);
        //data.constructWandQ(data.pseudoGwasEffectTrn, data.pseudoGwasNtrn);

        data.initVariances(opt.heritability, opt.propVarRandom);
        bool print = false;
        Model *modeli = new ApproxBayesR(data, data.lowRankModel, data.varGenotypic, data.varResidual, opt.pis, opt.piPar, opt.gamma, opt.estimatePi, opt.estimateSigmaSq, opt.noscale, opt.originalModel, opt.overdispersion, opt.estimatePS, opt.spouseCorrelation, opt.diagnosticMode, opt.robustMode, opt.algorithm, print);
        
        vector<McmcSamples*> mcmcSampleVeci;
        MCMC mcmc;

        unsigned chainLength = 150;
        unsigned burnin = 100;
        unsigned thin = 1;
        mcmcSampleVeci = mcmc.run(*modeli, chainLength, burnin, thin, print, opt.outputFreq, opt.title, print, print);

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
    
    data.numKeptInds = numKeptInds;
    
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

//void GCTB::autoDetermineNumComponents(Data &data, Options &opt, const VectorXf &pis, const VectorXf &gamma, const unsigned numiters){
//    vector<McmcSamples*> mcmcSampleVec = gctb.runMcmc(*model, opt.chainLength, opt.burnin, opt.thin,
//                                                  opt.outputFreq, opt.title, opt.writeBinPosterior, opt.writeTxtPosterior);
//
//}

