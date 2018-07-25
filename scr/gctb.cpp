//
//  gctb.cpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "gctb.hpp"

void GCTB::inputIndInfo(Data &data, const string &bedFile, const string &phenotypeFile, const string &keepIndFile, const unsigned keepIndMax, const unsigned mphen, const string &covariateFile){
    data.readFamFile(bedFile + ".fam");
    data.readPhenotypeFile(phenotypeFile, mphen);
    data.readCovariateFile(covariateFile);
    data.keepMatchedInd(keepIndFile, keepIndMax);
}

void GCTB::inputSnpInfo(Data &data, const string &bedFile, const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile, const unsigned includeChr, const string &skeletonSnpFile, const string &geneticMapFile, const bool readGenotypes){
    data.readBimFile(bedFile + ".bim");
    if (!includeSnpFile.empty()) data.includeSnp(includeSnpFile);
    if (!excludeSnpFile.empty()) data.excludeSnp(excludeSnpFile);
    if (includeChr) data.includeChr(includeChr);
    if (!excludeRegionFile.empty()) data.excludeRegion(excludeRegionFile);
    if (!skeletonSnpFile.empty()) data.includeSkeletonSnp(skeletonSnpFile);
    if (!geneticMapFile.empty()) data.readGeneticMapFile(geneticMapFile);
    data.includeMatchedSnp();
    if (readGenotypes) data.readBedFile(bedFile + ".bed");
}

void GCTB::inputSnpInfo(Data &data, const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile, const string &gwasSummaryFile, const string &ldmatrixFile, const unsigned includeChr, const string &skeletonSnpFile, const string &geneticMapFile, const string &annotationFile, const bool multiLDmat, const bool excludeMHC, const float afDiff){
    if (multiLDmat)
        data.readMultiLDmatInfoFile(ldmatrixFile);
    else
        data.readLDmatrixInfoFile(ldmatrixFile + ".info");
    if (!includeSnpFile.empty()) data.includeSnp(includeSnpFile);
    if (!excludeSnpFile.empty()) data.excludeSnp(excludeSnpFile);
    data.includeChr(includeChr);
    if (!excludeRegionFile.empty()) data.excludeRegion(excludeRegionFile);
    if (excludeMHC) data.excludeMHC();
    if (!skeletonSnpFile.empty()) data.includeSkeletonSnp(skeletonSnpFile);
    if (!geneticMapFile.empty()) data.readGeneticMapFile(geneticMapFile);
    if (!annotationFile.empty()) data.readAnnotationFile(annotationFile);
    if (!gwasSummaryFile.empty()) data.readGwasSummaryFile(gwasSummaryFile, afDiff);
    data.includeMatchedSnp();
    if (geneticMapFile.empty()) {
        if (multiLDmat)
            data.readMultiLDmatBinFile(ldmatrixFile);
        else
            data.readLDmatrixBinFile(ldmatrixFile + ".bin");
    } else {
        if (multiLDmat)
            data.readMultiLDmatBinFileAndShrink(ldmatrixFile);
        else
            data.readLDmatrixBinFileAndShrink(ldmatrixFile + ".bin");
    }
    if (!gwasSummaryFile.empty()) data.buildSparseMME();
}

void GCTB::inputSnpInfo(Data &data, const string &bedFile, const string &gwasSummaryFile, const float afDiff){
    data.readFamFile(bedFile + ".fam");
    data.readBimFile(bedFile + ".bim");

    data.keptIndInfoVec = data.makeKeptIndInfoVec(data.indInfoVec);
    data.numKeptInds =  (unsigned) data.keptIndInfoVec.size();
    
    data.readGwasSummaryFile(gwasSummaryFile, afDiff);
    data.includeMatchedSnp();
    data.readBedFile(bedFile + ".bed");
    data.buildSparseMME();
}

Model* GCTB::buildModel(Data &data, const string &bedFile, const string &gwasFile, const string &bayesType, const unsigned windowWidth,
                        const float heritability, const float pi, const bool estimatePi, const VectorXf &pis, const VectorXf &gamma,
                        const float phi, const float kappa, const string &algorithm, const unsigned snpFittedPerWindow,
                        const float varS, const vector<float> &S, const float overdispersion, const bool estimatePS,
                        const float icrsq, const float spouseCorrelation, const bool diagnosticMode){
    data.initVariances(heritability);
    if (!gwasFile.empty()) {
        if (data.numAnnos) {
            if (bayesType == "S")
                return new StratApproxBayesS(data, data.varGenotypic, data.varResidual, pi, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, S, algorithm);
            else
                throw(" Error: Wrong bayes type: " + bayesType + " in the annotation-stratified summary-data-based Bayesian analysis.");
        }
        else {
            if (bayesType == "C")
                return new ApproxBayesC(data, data.varGenotypic, data.varResidual, pi, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, diagnosticMode);
            else if (bayesType == "S")
                return new ApproxBayesS(data, data.varGenotypic, data.varResidual, pi, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, S, algorithm, diagnosticMode);
            else if (bayesType == "R")
                return new ApproxBayesR(data, data.varGenotypic, data.varResidual, pis, gamma, estimatePi, icrsq, spouseCorrelation);
            else if (bayesType == "Kap")
                return new ApproxBayesKappa(data, data.varGenotypic, data.varResidual, pis, gamma, estimatePi, icrsq, spouseCorrelation, kappa);
            else
                throw(" Error: Wrong bayes type: " + bayesType + " in the summary-data-based Bayesian analysis.");
        }
    }
    if (bayesType == "B") {
        data.readBedFile(bedFile + ".bed");
        return new BayesB(data, data.varGenotypic, data.varResidual, pi, estimatePi);
    }
    if (bayesType == "C") {
        data.readBedFile(bedFile + ".bed");
        return new BayesC(data, data.varGenotypic, data.varResidual, pi, estimatePi, algorithm);
    } 
    if (bayesType == "R") {
        data.readBedFile(bedFile + ".bed");
        return new BayesR(data, data.varGenotypic, data.varResidual, pis, gamma, estimatePi, algorithm);
    }
    else if (bayesType == "S") {
        data.readBedFile(bedFile + ".bed");
        return new BayesS(data, data.varGenotypic, data.varResidual, pi, estimatePi, varS, S, algorithm);
    }
    else if (bayesType == "N") {
        data.readBedFile(bedFile + ".bed");
        data.getNonoverlapWindowInfo(windowWidth);
        return new BayesN(data, data.varGenotypic, data.varResidual, pi, estimatePi, snpFittedPerWindow);
    }
    else if (bayesType == "NS") {
        data.readBedFile(bedFile + ".bed");
        data.getNonoverlapWindowInfo(windowWidth);
        return new BayesNS(data, data.varGenotypic, data.varResidual, pi, estimatePi, varS, S, snpFittedPerWindow, algorithm);
    }
    else if (bayesType == "Cap") {
        //data.readBedFile(bedFile + ".bed");
        data.buildSparseMME(bedFile + ".bed", windowWidth);
        return new ApproxBayesC(data, data.varGenotypic, data.varResidual, pi, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, diagnosticMode);
    }
    else if (bayesType == "Sap") {
        data.buildSparseMME(bedFile + ".bed", windowWidth);
        return new ApproxBayesS(data, data.varGenotypic, data.varResidual, pi, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, S, algorithm, diagnosticMode);
    }
    else {
        throw(" Error: Wrong bayes type: " + bayesType);
    }
}

vector<McmcSamples*> GCTB::runMcmc(Model &model, const unsigned chainLength, const unsigned burnin, const unsigned thin, const unsigned outputFreq, const string &title, const bool writeBinPosterior){
    MCMC mcmc;
    return mcmc.run(model, chainLength, burnin, thin, outputFreq, title, writeBinPosterior);
}

void GCTB::saveMcmcSamples(const vector<McmcSamples*> &mcmcSampleVec, const string &filename){
    for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
        McmcSamples *mcmcSamples = mcmcSampleVec[i];
        if (mcmcSamples->label == "SnpEffects" )  continue;
        if (mcmcSamples->label == "WindowDelta") continue;
        mcmcSamples->writeDataTxt(filename);
    }
}

void GCTB::outputResults(const Data &data, const vector<McmcSamples*> &mcmcSampleVec, const string &filename){
    vector<McmcSamples*> mcmcSamplesPar;
    for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
        McmcSamples *mcmcSamples = mcmcSampleVec[i];
        if (mcmcSamples->label == "SnpEffects") {
            //mcmcSamples->readDataBin(mcmcSamples->filename);
            data.outputSnpResults(mcmcSamples->posteriorMean, mcmcSamples->posteriorSqrMean, mcmcSamples->pip, filename + ".snpRes");
        }
        else if (mcmcSamples->label == "CovEffects") {
            data.outputFixedEffects(mcmcSamples->datMat, filename + ".covRes");
        }
        else if (mcmcSamples->label == "WindowDelta") {
            //mcmcSamples->readDataBin(mcmcSamples->filename);
            data.outputWindowResults(mcmcSamples->posteriorMean, filename + ".window");
        } else {
            mcmcSamplesPar.push_back(mcmcSamples);
        }
    }
}

McmcSamples* GCTB::inputMcmcSamples(const string &mcmcSampleFile, const string &label, const string &fileformat){
    if (myMPI::rank==0) cout << "reading MCMC samples for " << label << endl;
    McmcSamples *mcmcSamples = new McmcSamples(label);
    if (fileformat == "bin") mcmcSamples->readDataBin(mcmcSampleFile + "." + label);
    if (fileformat == "txt") mcmcSamples->readDataTxt(mcmcSampleFile + "." + label);
    return mcmcSamples;
}

void GCTB::estimateHsq(const Data &data, const McmcSamples &snpEffects, const string &filename){
    Heritability hsq;
    float phenVar = Gadget::calcVariance(data.y);
    hsq.getEstimate(data, snpEffects, phenVar);
    hsq.writeRes(filename);
    hsq.writeMcmcSamples(filename);
}

void GCTB::inputSnpResults(Data &data, const string &snpResFile){
    
}

void GCTB::predict(const Data &data, const string &filename){
    Predict pred;
    pred.getAccuracy(data, filename + ".predRes");
    pred.writeRes(data, filename + ".ghat");
}

void GCTB::clearGenotypes(Data &data){
    data.X.resize(0,0);
}
