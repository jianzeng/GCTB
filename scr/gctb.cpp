//
//  gctb.cpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "gctb.hpp"

void GCTB::inputIndInfo(Data &data, const string &bedFile, const string &phenotypeFile, const string &keepIndFile, const unsigned keepIndMax, const unsigned mphen, const string &covariateFile, const string &residualDiagFile){
    data.readFamFile(bedFile + ".fam");
    data.readPhenotypeFile(phenotypeFile, mphen);
    data.readCovariateFile(covariateFile);
    data.readResidualDiagFile(residualDiagFile);
    data.keepMatchedInd(keepIndFile, keepIndMax);
}

void GCTB::inputSnpInfo(Data &data, const string &bedFile, const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile, const unsigned includeChr, const bool excludeAmbiguousSNP, const string &skeletonSnpFile, const string &geneticMapFile,  const string &annotationFile, const bool transpose, const string &continuousAnnoFile, const unsigned flank, const string &eQTLFile, const float mafmin, const float mafmax, const bool noscale, const bool readGenotypes){
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
    if (!gwasSummaryFile.empty()) data.readGwasSummaryFile(gwasSummaryFile, afDiff, mafmin, mafmax, pValueThreshold, imputeN);
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
    if(opt.bLDMeigen) data.truncBlockEigen(opt.eigCutMethod, opt.eigThreshold); 
    if (!gwasSummaryFile.empty()) data.buildSparseMME(sampleOverlap, noscale);
    if (!windowFile.empty()) data.binSnpByWindowID();
}

void GCTB::inputSnpInfo(Data &data, const string &bedFile, const string &gwasSummaryFile, const float afDiff, const float mafmin, const float mafmax, const float pValueThreshold, const bool sampleOverlap, const bool imputeN, const bool noscale){
    data.readFamFile(bedFile + ".fam");
    data.readBimFile(bedFile + ".bim");

    data.keptIndInfoVec = data.makeKeptIndInfoVec(data.indInfoVec);
    data.numKeptInds =  (unsigned) data.keptIndInfoVec.size();
    
    data.readGwasSummaryFile(gwasSummaryFile, afDiff, mafmin, mafmax, pValueThreshold, imputeN);
    data.includeMatchedSnp();
    data.readBedFile(noscale, bedFile + ".bed");
    data.buildSparseMME(sampleOverlap, noscale);
}

Model* GCTB::buildModel(Data &data, const string &bedFile, const string &gwasFile, const string &bayesType, const unsigned windowWidth,
                        const float heritability, const float pi, const float piAlpha, const float piBeta, const bool estimatePi, const bool noscale,
                        const VectorXf &pis, const VectorXf &piPar, const VectorXf &gamma, const bool estimateSigmaSq,
                        const float phi, const float kappa, const string &algorithm, const unsigned snpFittedPerWindow,
                        const float varS, const vector<float> &S, const float overdispersion, const bool estimatePS,
                        const float icrsq, const float spouseCorrelation, const bool diagnosticMode, const bool originalModel, const bool perSnpGV, const bool robustMode){
    data.initVariances(heritability);
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
                return new StratApproxBayesS(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, S, algorithm, robustMode);
            else if (bayesType == "RC")
                return new ApproxBayesRC(data, data.varGenotypic, data.varResidual, pis, piPar, gamma, estimatePi, estimateSigmaSq, noscale, originalModel, perSnpGV, overdispersion, estimatePS, spouseCorrelation, diagnosticMode, robustMode, algorithm);
            else
                throw(" Error: Wrong bayes type: " + bayesType + " in the annotation-stratified summary-data-based Bayesian analysis.");
        }
        else {
            if (bayesType == "C")
                return new ApproxBayesC(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, noscale, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, diagnosticMode, robustMode);
            else if (bayesType == "B")
            return new ApproxBayesB(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, noscale, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, diagnosticMode, robustMode);
            else if (bayesType == "S")
                return new ApproxBayesS(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, S, algorithm, diagnosticMode, robustMode);
            else if (bayesType == "ST")
                return new ApproxBayesST(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, overdispersion, estimatePS, varS, S, true);
            else if (bayesType == "T")
                return new ApproxBayesST(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, overdispersion, estimatePS, varS, S, false);
            else if (bayesType == "SMix")
                return new ApproxBayesSMix(data, data.varGenotypic, data.varResidual, pi, overdispersion, estimatePS, varS, S);
            else if (bayesType == "R")
                return new ApproxBayesR(data, data.varGenotypic, data.varResidual, pis, piPar, gamma, estimatePi, estimateSigmaSq, noscale, originalModel, overdispersion, estimatePS, spouseCorrelation, diagnosticMode, robustMode, algorithm);
            else if(bayesType == "Reigen")
                return new ApproxBayesReigen(data, data.varGenotypic, data.varResidual, pis, piPar, gamma, estimatePi, estimateSigmaSq, noscale, originalModel, algorithm, false, true);
            else if (bayesType == "Kap")
                return new ApproxBayesKappa(data, data.varGenotypic, data.varResidual, pis, piPar, gamma, estimatePi, noscale, originalModel, icrsq, kappa);
            else if (bayesType == "RS")
                return new ApproxBayesRS(data, data.varGenotypic, data.varResidual, pis, piPar, gamma, estimatePi, varS, S, algorithm, noscale, originalModel, overdispersion, estimatePS, spouseCorrelation, diagnosticMode, robustMode, algorithm);
            else
                throw(" Error: Wrong bayes type: " + bayesType + " in the summary-data-based Bayesian analysis.");
        }
    }
    if (data.numAnnos) {
        if (bayesType == "RC") {
            data.readBedFile(noscale, bedFile + ".bed");
            return new BayesRC(data, data.varGenotypic, data.varResidual, pis, piPar, gamma, estimatePi, noscale, originalModel, algorithm);
        }
        else
            throw(" Error: Wrong bayes type: " + bayesType + " in the annotation-stratified Bayesian analysis.");
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
        return new BayesR(data, data.varGenotypic, data.varResidual, pis, piPar, gamma, estimatePi, noscale, originalModel, algorithm);
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
    else if (bayesType == "RS") {
        data.readBedFile(noscale, bedFile + ".bed");
        return new BayesRS(data, data.varGenotypic, data.varResidual, pis, piPar, gamma, estimatePi, varS, S, noscale, originalModel, algorithm);
    }
    else if (bayesType == "Cap") {
        //data.readBedFile(bedFile + ".bed");
        data.buildSparseMME(bedFile + ".bed", windowWidth);
        return new ApproxBayesC(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, noscale, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, diagnosticMode, robustMode);
    }
    else if (bayesType == "Sap") {
        data.buildSparseMME(bedFile + ".bed", windowWidth);
        return new ApproxBayesS(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, S, algorithm, diagnosticMode, robustMode);
    }
    else {
        throw(" Error: Wrong bayes type: " + bayesType);
    }
}

vector<McmcSamples*> GCTB::runMcmc(Model &model, const unsigned chainLength, const unsigned burnin, const unsigned thin, const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior){
    MCMC mcmc;
    return mcmc.run(model, chainLength, burnin, thin, true, outputFreq, title, writeBinPosterior, writeTxtPosterior);
}

vector<McmcSamples*> GCTB::multi_chain_mcmc(Data &data, const string &bayesType, const unsigned windowWidth, const float heritability, const float pi, const float piAlpha, const float piBeta, const bool estimatePi, const VectorXf &pis, const VectorXf &gamma, const float phi, const float kappa, const string &algorithm, const unsigned snpFittedPerWindow, const float varS, const vector<float> &S, const float overdispersion, const bool estimatePS, const float icrsq, const float spouseCorrelation, const bool diagnosticMode, const bool robustMode, const unsigned numChains, const unsigned chainLength, const unsigned burnin, const unsigned thin, const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior){
    
    data.initVariances(heritability);

    vector<Model*> modelVec(numChains);
    
    for (unsigned i=0; i<numChains; ++i) {
        if (data.numAnnos) {
            if (bayesType == "S")
                modelVec[i] = new StratApproxBayesS(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, S, algorithm, robustMode, true, !i);
            else
                throw(" Error: " + bayesType + " is not available in the multi-chain annotation-stratified Bayesian analysis.");
        }
        else {
            if (bayesType == "C")
                modelVec[i] = new ApproxBayesC(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, diagnosticMode, robustMode, true, !i);
            else if (bayesType == "S")
                modelVec[i] = new ApproxBayesS(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, varS, S, algorithm, diagnosticMode, robustMode, true, !i);
            else if (bayesType == "ST")
                modelVec[i] = new ApproxBayesST(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, overdispersion, estimatePS, varS, S, true, true, !i);
            else if (bayesType == "T")
                modelVec[i] = new ApproxBayesST(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, overdispersion, estimatePS, varS, S, false, true, !i);
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
            data.outputSnpResults(mcmcSamples->posteriorMean, mcmcSamples->posteriorSqrMean, mcmcSamples->lastSample, mcmcSamples->pip, noscale, filename + ".snpRes");
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
        McmcSamples *deltaPi2 = NULL;
        McmcSamples *deltaPi3 = NULL;
        McmcSamples *deltaPi4 = NULL;
        for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
            if (mcmcSampleVec[i]->label == "SnpEffects") snpEffects = mcmcSampleVec[i];
            if (mcmcSampleVec[i]->label == "DeltaPi2") deltaPi2 = mcmcSampleVec[i];
            if (mcmcSampleVec[i]->label == "DeltaPi3") deltaPi3 = mcmcSampleVec[i];
            if (mcmcSampleVec[i]->label == "DeltaPi4") deltaPi4 = mcmcSampleVec[i];
        }
        string newfilename = filename + ".snpRes";
        ofstream out(newfilename.c_str());
        out << boost::format("%6s %20s %6s %12s %6s %6s %12s %12s %12s %8s %8s %8s %8s\n")
        % "Id"
        % "Name"
        % "Chrom"
        % "Position"
        % "A1"
        % "A2"
        % "A1Frq"
        % "A1Effect"
        % "SE"
        % "PIP"
        % "Pi2"
        % "Pi3"
        % "Pi4";
        for (unsigned i=0, idx=0; i<data.numSnps; ++i) {
            SnpInfo *snp = data.snpInfoVec[i];
            if(!data.fullSnpFlag[i]) continue;
            float sqrt2pq = sqrt(2.0*snp->af*(1.0-snp->af));
            float effect = (snp->flipped ? - snpEffects->posteriorMean[idx] : snpEffects->posteriorMean[idx]);
            float se = sqrt(snpEffects->posteriorSqrMean[idx]-snpEffects->posteriorMean[idx]*snpEffects->posteriorMean[idx]);
            out << boost::format("%6s %20s %6s %12s %6s %6s %12.6f %12.6f %12.6f %8.3f %8.3f %8.3f %8.3f\n")
            % (i+1)
            % snp->ID
            % snp->chrom
            % snp->physPos
            % (snp->flipped ? snp->a2 : snp->a1)
            % (snp->flipped ? snp->a1 : snp->a2)
            % (snp->flipped ? 1.0-snp->af : snp->af)
            % (noscale ? effect : effect/sqrt2pq)
            % (noscale ? se : se/sqrt2pq)
            % snpEffects->pip[i]
            % deltaPi2->posteriorMean[i]
            % deltaPi3->posteriorMean[i]
            % deltaPi4->posteriorMean[i];
            ++idx;
        }
        out.close();
    }

}

McmcSamples* GCTB::inputMcmcSamples(const string &mcmcSampleFile, const string &label, const string &fileformat){
    cout << "reading MCMC samples for " << label << endl;
    McmcSamples *mcmcSamples = new McmcSamples(label);
    if (fileformat == "bin") mcmcSamples->readDataBin(mcmcSampleFile + "." + label);
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
    data.readGwasSummaryFile(gwasSummaryFile, 1, 0, 0, pValueThreshold, imputeN);
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

