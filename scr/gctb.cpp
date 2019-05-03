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
void GCTB::inputSnpInfo(Data &data, const string &bedFile, const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile, const unsigned includeChr, const bool excludeAmbiguousSNP, const string &skeletonSnpFile, const string &geneticMapFile, const float mafmin, const float mafmax, const bool noscale, const bool readGenotypes){
    data.readBimFile(bedFile + ".bim");
    if (!includeSnpFile.empty()) data.includeSnp(includeSnpFile);
    if (!excludeSnpFile.empty()) data.excludeSnp(excludeSnpFile);
    if (includeChr) data.includeChr(includeChr);
    if (excludeAmbiguousSNP) data.excludeAmbiguousSNP();
//    if (mafmin || mafmax) data.excludeSNPwithMaf(mafmin, mafmax);  // need to read in genotype data first
    if (!excludeRegionFile.empty()) data.excludeRegion(excludeRegionFile);
    if (!skeletonSnpFile.empty()) data.includeSkeletonSnp(skeletonSnpFile);
    if (!geneticMapFile.empty()) data.readGeneticMapFile(geneticMapFile);
    data.includeMatchedSnp();
    if (readGenotypes) data.readBedFile(noscale, bedFile + ".bed");
}

void GCTB::inputSnpInfo(Data &data, const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile, const string &gwasSummaryFile, const string &ldmatrixFile, const unsigned includeChr, const bool excludeAmbiguousSNP, const string &skeletonSnpFile, const string &geneticMapFile, const float genMapN, const string &annotationFile, const bool transpose, const string &continuousAnnoFile, const unsigned flank, const string &eQTLFile, const string &ldscoreFile, const bool multiLDmat, const bool excludeMHC, const float afDiff, const float mafmin, const float mafmax, const bool sampleOverlap, const bool imputeN, const string &bayesType, const bool noscale){
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
    if (!gwasSummaryFile.empty()) data.readGwasSummaryFile(gwasSummaryFile, afDiff, mafmin, mafmax, imputeN);
    data.includeMatchedSnp();
 //   if (geneticMapFile.empty()) {
        if (multiLDmat)
            data.readMultiLDmatBinFile(ldmatrixFile);
        else
            data.readLDmatrixBinFile(ldmatrixFile + ".bin");
//    } else {
//        if (multiLDmat)
//            data.readMultiLDmatBinFileAndShrink(ldmatrixFile, genMapN);
//        else
//            data.readLDmatrixBinFileAndShrink(ldmatrixFile + ".bin");
//    }
    if (!gwasSummaryFile.empty()) data.buildSparseMME(sampleOverlap, bayesType, noscale);
}

void GCTB::inputSnpInfo(Data &data, const string &bedFile, const string &gwasSummaryFile, const float afDiff, const float mafmin, const float mafmax, const bool sampleOverlap, const bool imputeN, const string &bayesType, const bool noscale){
    data.readFamFile(bedFile + ".fam");
    data.readBimFile(bedFile + ".bim");

    data.keptIndInfoVec = data.makeKeptIndInfoVec(data.indInfoVec);
    data.numKeptInds =  (unsigned) data.keptIndInfoVec.size();
    
    data.readGwasSummaryFile(gwasSummaryFile, afDiff, mafmin, mafmax, imputeN);
    data.includeMatchedSnp();
    data.readBedFile(noscale, bedFile + ".bed");
    data.buildSparseMME(sampleOverlap, bayesType, noscale);
}

Model* GCTB::buildModel(Data &data, const string &bedFile, const string &gwasFile, const string &bayesType, const unsigned windowWidth,
                        const float heritability, const float pi, const float piAlpha, const float piBeta, const bool estimatePi, const bool noscale, const VectorXf &pis, const VectorXf &gamma,
                        const float phi, const float kappa, const string &algorithm, const unsigned snpFittedPerWindow,
                        const float varS, const vector<float> &S, const float overdispersion, const bool estimatePS,
                        const float icrsq, const float spouseCorrelation, const bool diagnosticMode){
    data.initVariances(heritability);
    if (!gwasFile.empty()) {
        if (data.numAnnos) {
            if (bayesType == "S")
                return new StratApproxBayesS(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, S, algorithm);
            else
                throw(" Error: Wrong bayes type: " + bayesType + " in the annotation-stratified summary-data-based Bayesian analysis.");
        }
        else {
            if (bayesType == "C")
                return new ApproxBayesC(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, noscale, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, diagnosticMode);
            else if (bayesType == "S")
                return new ApproxBayesS(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, S, algorithm, diagnosticMode);
            else if (bayesType == "ST")
                return new ApproxBayesST(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, overdispersion, estimatePS, varS, S, true);
            else if (bayesType == "T")
                return new ApproxBayesST(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, overdispersion, estimatePS, varS, S, false);
            else if (bayesType == "SMix")
                return new ApproxBayesSMix(data, data.varGenotypic, data.varResidual, pi, overdispersion, estimatePS, varS, S);
            else if (bayesType == "R")
                return new ApproxBayesR(data, data.varGenotypic, data.varResidual, pis, piAlpha, piBeta, gamma, estimatePi, noscale, overdispersion, estimatePS, spouseCorrelation, diagnosticMode);
            else if (bayesType == "Kap")
                return new ApproxBayesKappa(data, data.varGenotypic, data.varResidual, pis, piAlpha, piBeta, gamma, estimatePi, noscale, icrsq, kappa);
            else
                throw(" Error: Wrong bayes type: " + bayesType + " in the summary-data-based Bayesian analysis.");
        }
    }
    if (bayesType == "B") {
        data.readBedFile(noscale, bedFile + ".bed");
        return new BayesB(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, noscale);
    }
    if (bayesType == "C") {
         data.readBedFile(noscale, bedFile + ".bed");
        return new BayesC(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, noscale, algorithm);
    } 
    if (bayesType == "R") {
        data.readBedFile(noscale, bedFile + ".bed");
        return new BayesR(data, data.varGenotypic, data.varResidual, pis, piAlpha, piBeta, gamma, estimatePi, noscale, algorithm);
    }
    else if (bayesType == "S") {
        data.readBedFile(noscale, bedFile + ".bed");
        return new BayesS(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, varS, S, algorithm);
    }
    else if (bayesType == "SMix") {
        data.readBedFile(noscale, bedFile + ".bed");
        return new BayesSMix(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, varS, S, algorithm);
    }
    else if (bayesType == "N") {
        data.readBedFile(noscale, bedFile + ".bed");
        data.getNonoverlapWindowInfo(windowWidth);
        return new BayesN(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, noscale, snpFittedPerWindow);
    }
    else if (bayesType == "NS") {
        data.readBedFile(noscale, bedFile + ".bed");
        data.getNonoverlapWindowInfo(windowWidth);
        return new BayesNS(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, varS, S, snpFittedPerWindow, algorithm);
    }
    else if (bayesType == "Cap") {
        //data.readBedFile(bedFile + ".bed");
        data.buildSparseMME(bedFile + ".bed", windowWidth);
        return new ApproxBayesC(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, noscale, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, diagnosticMode);
    }
    else if (bayesType == "Sap") {
        data.buildSparseMME(bedFile + ".bed", windowWidth);
        return new ApproxBayesS(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, S, algorithm, diagnosticMode);
    }
    else {
        throw(" Error: Wrong bayes type: " + bayesType);
    }
}

vector<McmcSamples*> GCTB::runMcmc(Model &model, const unsigned chainLength, const unsigned burnin, const unsigned thin, const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior){
    MCMC mcmc;
    return mcmc.run(model, chainLength, burnin, thin, true, outputFreq, title, writeBinPosterior, writeTxtPosterior);
}

vector<McmcSamples*> GCTB::multi_chain_mcmc(Data &data, const string &bayesType, const unsigned windowWidth, const float heritability, const float pi, const float piAlpha, const float piBeta, const bool estimatePi, const VectorXf &pis, const VectorXf &gamma, const float phi, const float kappa, const string &algorithm, const unsigned snpFittedPerWindow, const float varS, const vector<float> &S, const float overdispersion, const bool estimatePS, const float icrsq, const float spouseCorrelation, const bool diagnosticMode, const unsigned numChains, const unsigned chainLength, const unsigned burnin, const unsigned thin, const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior){
    
    data.initVariances(heritability);

    vector<Model*> modelVec(numChains);
    
    for (unsigned i=0; i<numChains; ++i) {
        if (data.numAnnos) {
            if (bayesType == "S")
                modelVec[i] = new StratApproxBayesS(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, S, algorithm, true, !i);
            else
                throw(" Error: " + bayesType + " is not available in the multi-chain annotation-stratified Bayesian analysis.");
        }
        else {
            if (bayesType == "C")
                modelVec[i] = new ApproxBayesC(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, diagnosticMode, true, !i);
            else if (bayesType == "S")
                modelVec[i] = new ApproxBayesS(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, S, algorithm, diagnosticMode, true, !i);
            else if (bayesType == "ST")
                modelVec[i] = new ApproxBayesST(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, overdispersion, estimatePS, varS, S, true, true, !i);
            else if (bayesType == "T")
                modelVec[i] = new ApproxBayesST(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, overdispersion, estimatePS, varS, S, false, true, !i);
            else
                throw(" Error: " + bayesType + " is not available in the multi-chain Bayesian analysis.");
        }
    }
    
    vector<vector<McmcSamples*> > mcmcSampleVecChain;
    mcmcSampleVecChain.resize(numChains);
    
//    if (myMPI::rank==0)
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
            data.outputSnpResults(mcmcSamples->posteriorMean, mcmcSamples->posteriorSqrMean, mcmcSamples->pip, noscale, filename + ".snpRes");
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
    if (bayesType == "SMix") {
        McmcSamples *snpEffects = NULL;
        McmcSamples *delta = NULL;
        for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
            if (mcmcSampleVec[i]->label == "SnpEffects") snpEffects = mcmcSampleVec[i];
            if (mcmcSampleVec[i]->label == "DeltaS") delta = mcmcSampleVec[i];
        }
//        if (myMPI::rank) return;
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
}

McmcSamples* GCTB::inputMcmcSamples(const string &mcmcSampleFile, const string &label, const string &fileformat){
//    if (myMPI::rank==0)
        cout << "reading MCMC samples for " << label << endl;
    McmcSamples *mcmcSamples = new McmcSamples(label);
    if (fileformat == "bin") mcmcSamples->readDataBin(mcmcSampleFile + "." + label);
//    if (fileformat == "txt") mcmcSamples->readDataTxt(mcmcSampleFile + "." + label);
    if (fileformat == "txt") mcmcSamples->readDataTxt(mcmcSampleFile + ".Par", label);
    return mcmcSamples;
}

void GCTB::estimateHsq(const Data &data, const McmcSamples &snpEffects, const string &filename){
    Heritability hsq;
    float phenVar = Gadget::calcVariance(data.y);
    hsq.getEstimate(data, snpEffects, phenVar);
    hsq.writeRes(filename);
    hsq.writeMcmcSamples(filename);
}

void GCTB::predict(const Data &data, const string &filename){
    Predict pred;
    pred.getAccuracy(data, filename + ".predRes");
    pred.writeRes(data, filename + ".ghat");
}

void GCTB::clearGenotypes(Data &data){
    data.X.resize(0,0);
}

void GCTB::stratify(Data &data, const string &ldmatrixFile, const bool multiLDmat, const string &geneticMapFile, const float genMapN, const string &snpResFile, const string &mcmcSampleFile, const string &annotationFile, const bool transpose, const string &continuousAnnoFile, const unsigned flank, const string &eQTLFile, const string &gwasSummaryFile, const bool imputeN, const string &filename, const string &bayesType, unsigned chainLength, unsigned burnin, const unsigned thin, const unsigned outputFreq){
    if (multiLDmat)
        data.readMultiLDmatInfoFile(ldmatrixFile);
    else
        data.readLDmatrixInfoFile(ldmatrixFile + ".info");
    data.inputSnpInfoAndResults(snpResFile, bayesType);
    if (!annotationFile.empty())
        data.readAnnotationFile(annotationFile, transpose, true);
    else
        data.readAnnotationFileFormat2(continuousAnnoFile, flank*1000, eQTLFile);
    data.readGwasSummaryFile(gwasSummaryFile, 1, 0, 0, imputeN);
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
    data.buildSparseMME(false, bayesType, true);
    data.makeAnnowiseSparseLDM(data.ZPZsp, data.annoInfoVec, data.snpInfoVec);
    
    McmcSamples *snpEffects = inputMcmcSamples(mcmcSampleFile, "SnpEffects", "bin");
    McmcSamples *hsq = inputMcmcSamples(mcmcSampleFile, "hsq", "txt");

    Model *model;
    
    if (bayesType == "S") {
        model = new PostHocStratifyS(data, *snpEffects, *hsq, thin, hsq->mean()[0]);
    }
    else if (bayesType == "SMix") {
        McmcSamples *deltaS = inputMcmcSamples(mcmcSampleFile, "DeltaS", "bin");
        model = new PostHocStratifySMix(data, *snpEffects, *hsq, *deltaS, thin, hsq->mean()[0]);
    }
    else
        throw(" Error: Wrong bayes type: " + bayesType + " in the annotation-stratified summary-data-based Bayesian analysis.");

    
    if (chainLength > snpEffects->nrow) chainLength = snpEffects->nrow;
    if (burnin > chainLength) burnin = 0.2*chainLength;
    
    runMcmc(*model, chainLength, burnin, thin, outputFreq, filename, false, false);
}

