//
//  options.cpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "options.hpp"
#include <limits>

void Options::inputOptions(const int argc, const char* argv[]){
    stringstream ss;
    for (unsigned i=1; i<argc; ++i) {
        if (!strcmp(argv[i], "--inp-file")) {
            optionFile = argv[++i];
            readFile(optionFile);
            return;
        } else {
            if (i==1) ss << "\nOptions:\n\n";
        }
        if (!strcmp(argv[i], "--bayes")) {
            analysisType = "Bayes";
            bayesType = argv[++i];
            ss << "--bayes " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--sbayes")) {
            analysisType = "SBayes";
            bayesType = argv[++i];
            ss << "--sbayes " << argv[i] << "\n";
        }else if(!strcmp(argv[i], "--eig-cutoff")){
            eigCutMethod = argv[++i];
            eigThreshold = stof(argv[++i]);
            cout << eigCutMethod << " " << argv[i] << std::endl;
            if(!strncmp(argv[i], "--", 2)){
                throw("--eig-cutoff need two paramters: method (value, count, percent) threshold");
            }
            ss << "--eig-cutoff " << eigCutMethod << " " << eigThreshold << std::endl;
        }else if (!strcmp(argv[i], "--cg")) {
            analysisType = "ConjugateGradient";
            ss << "--cg " << "\n";
        }
        else if (!strcmp(argv[i], "--estimate-hsq")) {
            analysisType = "hsq";
            ss << "--estimate-hsq " << "\n";
        }
        else if (!strcmp(argv[i], "--estimate-pi")) {
            analysisType = "Pi";
            ss << "--estimate-pi " << "\n";
        }
//        else if (!strcmp(argv[i], "--make-ldm")) {
//            analysisType = "LDmatrix";
//            ss << "--make-ldm " << "\n";
//        }
        else if (!strcmp(argv[i], "--make-full-ldm")) {
            analysisType = "LDmatrix";
            outLDmatType = "full";
            ss << "--make-full-ldm " << "\n";
        }
        else if (!strcmp(argv[i], "--make-band-ldm")) {
            analysisType = "LDmatrix";
            outLDmatType = "band";
            ss << "--make-band-ldm " << "\n";
        }
        else if (!strcmp(argv[i], "--make-shrunk-ldm")) {
            analysisType = "LDmatrix";
            outLDmatType = "shrunk";
            ss << "--make-shrunk-ldm " << "\n";
        }
        else if (!strcmp(argv[i], "--make-sparse-ldm")) {
            analysisType = "LDmatrix";
            outLDmatType = "sparse";
            ss << "--make-sparse-ldm " << "\n";
        }
        else if (!strcmp(argv[i], "--make-sparse-shrunk-ldm")) {
            analysisType = "LDmatrix";
            outLDmatType = "sparseshrunk";
            ss << "--make-sparse-shrunk-ldm " << "\n";
        }
        else if (!strcmp(argv[i], "--make-block-ldm")) {
            analysisType = "LDmatrix";
            outLDmatType = "block";
            ss << "--make-block-ldm " << "\n";
        }
        else if (!strcmp(argv[i], "--xci")) {
            analysisType = "XCI";
            bayesType = argv[++i];
            ss << "--xci " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--alg")) {
            algorithm = argv[++i];
            ss << "--alg " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--predict")) {
            analysisType = "Predict";
            ss << "--predict " << "\n";
        }
        else if (!strcmp(argv[i], "--bfile")) {
            bedFile = argv[++i];
            ss << "--bfile " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--pheno")) {
            phenotypeFile = argv[++i];
            ss << "--pheno " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--gen-map")) { // Genetic map file option
            geneticMapFile = argv[++i];
            ss << "--gen-map " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--annot")) {
            annotationFile = argv[++i];
            ss << "--annot " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--cont-annot")) {
            continuousAnnoFile = argv[++i];
            ss << "--cont-annot " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--ldsc")) {
            ldscoreFile = argv[++i];
            ss << "--ldsc " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--covar")) {
            covariateFile = argv[++i];
            ss << "--covar " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--random-covar")) {
            randomCovariateFile = argv[++i];
            ss << "--random-covar " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--res-diag")) {
            residualDiagFile = argv[++i];
            ss << "--res-diag " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--mpheno")) {
            mphen = atoi(argv[++i]);
            ss << "--mpheno " << mphen << "\n";
        }
        else if (!strcmp(argv[i], "--keep")) {
            keepIndFile = argv[++i];
            ss << "--keep " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--keep-max")) {
            keepIndMax = atoi(argv[++i]);
            ss << "--keep-max " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--extract")) {
            includeSnpFile = argv[++i];
            ss << "--extract " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--exclude")) {
            excludeSnpFile = argv[++i];
            ss << "--exclude " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--exclude-region")) {
            excludeRegionFile = argv[++i];
            ss << "--exclude-region " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--mcmc-samples")) {
            mcmcSampleFile = argv[++i];
            ss << "--mcmc-samples " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--gwas-summary")) {
            gwasSummaryFile = argv[++i];
            ss << "--gwas-summary " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--ldm")) {
            ldmatrixFile = argv[++i];
            ss << "--ldm " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--mldm")) {
            ldmatrixFile = argv[++i];
            multiLDmat = true;
            ss << "--mldm " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--make-ldm-eigen")) {
            analysisType = "LDmatrixEigen";
            ss << "--make-ldm-eigen " << "\n";
        }
        else if (!strcmp(argv[i], "--ldm-eigen")) {
            eigenMatrixFile = argv[++i];
            ss << "--ldm-eigen " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--ldm-eigen-cutoff")) {
            Gadget::Tokenizer strvec;
            strvec.getTokens(argv[++i], " ,");
            eigenCutoff.resize(strvec.size());
            for (unsigned j=0; j<strvec.size(); ++j) eigenCutoff[j] = stof(strvec[j]);
            std::sort(eigenCutoff.data(), eigenCutoff.data()+eigenCutoff.size());
            ss << "--ldm-eigen-cutoff " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--block-info")) {
            ldBlockInfoFile = argv[++i];
            ss << "--block-info " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--block")) {
            includeBlock = atoi(argv[++i]);
            ss << "--block " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--merge-block-ldm-info")) {
            analysisType = "LDmatrix";
            mergeLdm = true;
            outLDmatType = "block";
            ss << "--merge-block-ldm-info " << "\n";
        }
        else if (!strcmp(argv[i], "--impute-summary")) {
            analysisType = "ImputeSumStats";
            imputeSummary = true;
            ss << "--impute-summary " << "\n";
        }
        else if (!strcmp(argv[i], "--merge-block-gwas-summary")) {
            analysisType = "MergeGwasSummary";
            outLDmatType = "block";
            ss << "--merge-block-gwas-summary " << "\n";
        }
        else if (!strcmp(argv[i], "--per-snp-gv")) {
            perSnpGV = true;
            ss << "--per-snp-gv " << "\n";
        }
        else if (!strcmp(argv[i], "--snp-res")) {
            snpResFile = argv[++i];
            ss << "--snp-res " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--wind")) {
            windowWidth = unsigned(atof(argv[++i]) * Megabase);
            ss << "--wind " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--wind-file")) {
            windowFile = argv[++i];
            ss << "--wind-file " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--pi")) {
            Gadget::Tokenizer strvec;
            strvec.getTokens(argv[++i], " ,");
            if (strvec.size() != 1 && (bayesType != "R" && bayesType != "Kap" && bayesType != "RS" && bayesType != "RC"))
                throw("Error: When NOT using Bayes R or its variants option you can only specify one mixture proportion parameter.");
            if (strvec.size() == 1) {
                for (unsigned j=0; j<strvec.size(); ++j) pi = stof(strvec[j]);
            } else {
                pis.resize(strvec.size());
                for (unsigned j=0; j<strvec.size(); ++j) pis[j] = stof(strvec[j]);
            }
            ss << "--pi " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--pi-alpha")) {
            piAlpha = atof(argv[++i]);
            ss << "--pi-alpha " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--pi-beta")) {
            piBeta = atof(argv[++i]);
            ss << "--pi-beta " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--hsq")) {
            heritability = atof(argv[++i]);
            ss << "--hsq " << argv[i] << "\n";
        }
//        else if (!strcmp(argv[i], "--varg")) {
//            varGenotypic = atof(argv[++i]);
//            ss << "--varg " << argv[i] << "\n";
//        }
//        else if (!strcmp(argv[i], "--vare")) {
//            varResidual = atof(argv[++i]);
//            ss << "--vare " << argv[i] << "\n";
//        }
        else if (!strcmp(argv[i], "--var-random")) {
            propVarRandom = atof(argv[++i]);
            ss << "--var-random " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--num-chains")) {
            numChains = atoi(argv[++i]);
            ss << "--num-chains " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--chain-length")) {
            chainLength = atoi(argv[++i]);
            ss << "--chain-length " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--burn-in")) {
            burnin = atoi(argv[++i]);
            ss << "--burn-in " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--out-freq")) {
            outputFreq = atoi(argv[++i]);
            ss << "--out-freq " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--seed")) {
            seed = atoi(argv[++i]);
            ss << "--seed " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--wind-nnz")) {
            snpFittedPerWindow = atoi(argv[++i]);
            ss << "--wind-nnz " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--out")) {
            title = argv[++i];
            ss << "--out " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--no-mcmc-bin")) {
            writeBinPosterior = false;
            ss << "--no-mcmc-bin " << "\n";
        }
        else if (!strcmp(argv[i], "--no-mcmc-txt")) {
            writeTxtPosterior = false;
            ss << "--no-mcmc-txt " << "\n";
        }
        else if (!strcmp(argv[i], "--thin")) {
            thin = atoi(argv[++i]);
            ss << "--thin " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--fix-sigma2")) {
            estimateSigmaSq = false;
            ss << "--fix-sigma2 " << "\n";
        }
        else if (!strcmp(argv[i], "--fix-pi")) {
            estimatePi = false;
            ss << "--fix-pi " << "\n";
        }
        else if (!strcmp(argv[i], "--fix-pi-ndc")) {
            estimatePiNDC = false;
            ss << "--fix-pi-ndc " << "\n";
        }
        else if (!strcmp(argv[i], "--fix-pi-gxe")) {
            estimatePiGxE = false;
            ss << "--fix-pi-gxe " << "\n";
        }
        else if (!strcmp(argv[i], "--pi-par")) {
            Gadget::Tokenizer strvec;
            strvec.getTokens(argv[++i], " ,");
            piPar.resize(strvec.size());
            for (unsigned j=0; j<strvec.size(); ++j) {
                piPar[j] = stof(strvec[j]);
            }
            ss << "--pi-par " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--pi-ndc-par")) {
            Gadget::Tokenizer strvec;
            strvec.getTokens(argv[++i], " ,");
            if (strvec.size()>2) throw("Error: --pi-ndc-par only allow two parameters!");
            for (unsigned j=0; j<strvec.size(); ++j) {
                piNDCpar[j] = stof(strvec[j]);
            }
            ss << "--pi-ndc-par " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--varS")) {
            varS = atof(argv[++i]);
            ss << "--varS " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--kappa")) {
            kappa = atof(argv[++i]);
            ss << "--kappa " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--S")) {
            Gadget::Tokenizer strvec;
            strvec.getTokens(argv[++i], " ,");
            S.resize(strvec.size());
            for (unsigned j=0; j<strvec.size(); ++j) {
                S[j] = stof(strvec[j]);
            }
            ss << "--S " << argv[i] << "\n";
        }
        // else if (!strcmp(argv[i], "--set-pis")) {
        //     Gadget::Tokenizer strvec;
        //     strvec.getTokens(argv[++i], " ,");
        //     pis.resize(strvec.size());
        //     for (unsigned j=0; j<strvec.size(); ++j) {
        //         pis[j] = stof(strvec[j]);
        //     }
        //     ss << "--set-pis " << argv[i] << "\n";
        // }
        else if (!strcmp(argv[i], "--gamma")) {
            Gadget::Tokenizer strvec;
            strvec.getTokens(argv[++i], " ,");
            gamma.resize(strvec.size());
            for (unsigned j=0; j<strvec.size(); ++j) {
                gamma[j] = stof(strvec[j]);
            }
            ss << "--gamma " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--thread")) {
            numThread = atoi(argv[++i]);
            ss << "--thread " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--chr")) {
            includeChr = atoi(argv[++i]);
            ss << "--chr " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--ld")) {
            LDthreshold = atof(argv[++i]);
            ss << "--ld " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--rsq")) {
            rsqThreshold = atof(argv[++i]);
            ss << "--rsq " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--p-value")) {
            pValueThreshold = atof(argv[++i]);
            ss << "--p-value " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--bin-snp")) {
            analysisType = "LDmatrix";
            binSnp = true;
            ss << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--unscale-genotype")) {
            noscale = true;
            ss << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--snp")) {
            snpRange = argv[++i];
            ss << "--snp " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--part")) {
            partParam = argv[++i];
            ss << "--part " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--ne")) {
            effpopNE = atof(argv[++i]);
            ss << "--ne " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--genmap-n")) {
            genMapN = atof(argv[++i]);
            ss << "--genmap-n " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--shrunk-cutoff")) {
            cutOff = atof(argv[++i]);
            ss << "--shrunk-cutoff " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--multi-thread-eigen")) {
            multiThreadEigen = true;
            ss << "--multi-thread-eigen " << "\n";
        }
        else if (!strcmp(argv[i], "--chisq")) {
            chisqThreshold = atof(argv[++i]);
            ss << "--chisq " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--pi-ndc")) {
            piNDC = atof(argv[++i]);
            ss << "--pi-ndc " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--pi-gxe")) {
            piGxE = atof(argv[++i]);
            ss << "--pi-gxe " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--write-ldm-txt")) {
            writeLdmTxt = true;
            ss << "--write-ldm-txt " << "\n";
        }
        else if (!strcmp(argv[i], "--read-ldm-txt")) {
            readLdmTxt = true;
            ss << "--read-ldm-txt " << "\n";
        }
        else if (!strcmp(argv[i], "--exclude-mhc")) {
            excludeMHC = true;
            ss << "--exclude-mhc " << "\n";
        }
        else if (!strcmp(argv[i], "--phi")) {
            phi = atof(argv[++i]);
            ss << "--phi " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--overdispersion")) {
            overdispersion = atof(argv[++i]);
            ss << "--overdispersion " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--skeleton-snp")) {
            skeletonSnpFile = argv[++i];
            ss << "--skeleton-snp " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--direct-prune")) {
            directPrune = true;
            ss << "--direct-prune " << "\n";
        }
        else if (!strcmp(argv[i], "--estimate-ps")) {
            estimatePS = true;
            ss << "--estimate-ps " << "\n";
        }
        else if (!strcmp(argv[i], "--inter-chr-rsq")) {
            icrsq = atof(argv[++i]);
            ss << "--inter-chr-rsq " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--spouse-corr")) {
            spouseCorrelation = atof(argv[++i]);
            ss << "--spouse-corr " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--diagnostic-mode")) {
            diagnosticMode = true;
            ss << "--diagnostic-mode " << "\n";
        }
        else if (!strcmp(argv[i], "--filter-af-diff")) {
            afDiff = atof(argv[++i]);
            ss << "--filter-af-diff " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--jackknife")) {
            jackknife = true;
            ss << "--jackknife " << "\n";
        }
        else if (!strcmp(argv[i], "--ambiguous-snp")) {
            excludeAmbiguousSNP = true;
            ss << "--ambiguous-snp " << "\n";
        }
        else if (!strcmp(argv[i], "--impute-n")) {
            imputeN = true;
            ss << "--impute-n " << "\n";
        }
        else if (!strcmp(argv[i], "--maf")) {
            mafmin = atof(argv[++i]);
            ss << "--maf " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--max-maf")) {
            mafmax = atof(argv[++i]);
            ss << "--max-maf " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--transpose")) {
            transpose = true;
            ss << "--transpose " << "\n";
        }
        else if (!strcmp(argv[i], "--overlap")) {
            sampleOverlap = true;
            ss << "--overlap " << "\n";
        }
        else if (!strcmp(argv[i], "--stratify")) {
            analysisType = "Stratify";
            bayesType = argv[++i];
            ss << "--stratify " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--flank")) {
            flank = atof(argv[++i]);
            ss << "--flank " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--eqtl")) {
            eQTLFile = argv[++i];
            ss << "--eqtl " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--simu")) {
            simuMode = true;
            ss << "--simu " << "\n";
        }
        else if (!strcmp(argv[i], "--original-model")) {
            originalModel = true;
            ss << "--original-model " << "\n";
        }
        else if (!strcmp(argv[i], "--lambda")) {
            lambda = atof(argv[++i]);
            ss << "--lambda " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--two-stage-model")) {
            twoStageModel = true;
            ss << "--two-stage-model " << "\n";
        }
        else if (!strcmp(argv[i], "--robust")) {
            robustMode = true;
            ss << "--robust " << "\n";
        }
        else {
            stringstream errmsg;
            errmsg << "\nError: invalid option \"" << argv[i] << "\".\n";
            throw (errmsg.str());
        }
    }
    // Error throwing for Bayes R specific options
    // cout << "diff " << std::abs(pis.sum() - 1.0) << endl;
    float epsilon = std::numeric_limits<float>::epsilon();
    if (std::abs(pis.sum() - 1.0) > epsilon) 
    {
        throw("Error: When using Bayes R option and --pi the pis must sum to 1. Please adjust.");
    }
    if (pis.size() != gamma.size()) 
    {
        throw("Error: Length of mixing proportions vector " + to_string(static_cast<long long>(pis.size())) +
              " does not match length of variance scaling factors vector " + 
              to_string(static_cast<long long>(gamma.size())) +
              ". \n" + 
              "When using Bayes R option please specify starting mixing proportions and variance scaling factors." +
              "\n" + 
              "The flags for these are --pi and --gamma.");
    }
    
    // BayesS type of model do not allow scaled genotypes
    if (bayesType == "S" || bayesType == "ST" || bayesType == "T" || bayesType == "SMix" || bayesType == "RS") {
        noscale = true;
    }
   
    if (analysisType == "hsq") noscale = true;
    
    if (chainLength < burnin) {
        throw("Error: Chain length is smaller than burn-in!");
    }
    
    cout << ss.str() << endl;
    
    setThread();
}

void Options::readFile(const string &file){  // input options from file
    optionFile = file;
    stringstream ss;
    ss << "\nOptions:\n\n";
    ss << boost::format("%20s %-1s %-20s\n") %"optionFile" %":" %file;
    makeTitle();
    
    ifstream in(file.c_str());
    if (!in) throw ("Error: can not open the file [" + file + "] to read.");
    
    string key, value;
    while (in >> key >> value) {
        if (key == "bedFile") {
            bedFile = value;
        } else if (key == "phenotypeFile") {
            phenotypeFile = value;
        } else if (key == "mpheno") {
            mphen = stoi(value);
        } else if (key == "bedFile") {
            bedFile = value;
        } else if (key == "geneticMapFile") {
            geneticMapFile = value;
        } else if (key == "annotationFile") {
            annotationFile = value;
        } else if (key == "covariateFile") {
            covariateFile = value;
        } else if (key == "analysisType") {
            analysisType = value;
        } else if (key == "bayesType") {
            bayesType = value;
        } else if (key == "algorithm") {
            algorithm = value;
        } else if (key == "keepIndFile") {
            keepIndFile = value;
        } else if (key == "keepIndMax") {
            keepIndMax = stoi(value);
        } else if (key == "includeSnpFile") {
            includeSnpFile = value;
        } else if (key == "excludeSnpFile") {
            excludeSnpFile = value;
        } else if (key == "mcmcSampleFile") {
            mcmcSampleFile = value;
        } else if (key == "gwasSummaryFile") {
            gwasSummaryFile = value;
        } else if (key == "LDmatrixFile") {
            ldmatrixFile = value;
        } else if (key == "multiLDmatrixFile") {
            ldmatrixFile = value;
            multiLDmat = true;
        } else if (key == "snpResFile") {
            snpResFile = value;
        } else if (key == "windowWidth") {
            windowWidth = unsigned(stof(value) * Megabase);
        } else if (key == "pi") {
            pi = stof(value);
        } else if (key == "heritability") {
            heritability = stof(value);
//        } else if (key == "varGenotypic") {
//            varGenotypic = stof(value);
//        } else if (key == "varResidual") {
//            varResidual = stof(value);
        } else if (key == "chainLength") {
            chainLength = stoi(value);
        } else if (key == "burnin") {
            burnin = stoi(value);
        } else if (key == "outputFreq") {
            outputFreq = stoi(value);
        } else if (key == "seed") {
            seed = stoi(value);
        } else if (key == "snpFittedPerWindow") {
            snpFittedPerWindow = stoi(value);
        } else if (key == "writeBinPosterior" && value == "No") {
            writeBinPosterior = false;
        } else if (key == "writeTxtPosterior" && value == "No") {
            writeTxtPosterior = false;
        } else if (key == "thin") {
            thin = stoi(value);
        } else if (key == "estimatePi" && value == "No") {
            estimatePi = false;
        } else if (key == "estimatePiNDC" && value == "No") {
            estimatePiNDC = false;
        } else if (key == "estimatePiGxE" && value == "No") {
            estimatePiGxE = false;
        } else if (key == "outputResults" && value == "No") {
            outputResults = false;
        } else if (key == "varS") {
            varS = stof(value);
        } else if (key == "S") {
            Gadget::Tokenizer strvec;
            strvec.getTokens(value, " ,");
            S.resize(strvec.size());
            for (unsigned j=0; j<strvec.size(); ++j) {
                S[j] = stof(strvec[j]);
            }
        } else if (key == "numThread") {
            numThread = stoi(value);
        } else if (key == "includeChr") {
            includeChr = stoi(value);
        } else if (key == "LDthreshold") {
            LDthreshold = stof(value);
        } else if (key == "chisqThreshold") {
            chisqThreshold = stof(value);
        } else if (key == "snpRange") {
            snpRange = value;
        } else if (key == "multiThreadEigen") {
            multiThreadEigen = true;
        } else if (key == "outLDmatType") {
            outLDmatType = value;
        } else if (key == "piNDC") {
            piNDC = stof(value);
        } else if (key == "piGxE") {
            piGxE = stof(value);
        } else if (key == "writeLDmatTxt") {
            writeLdmTxt = true;
        } else if (key == "excludeMHC") {
            excludeMHC = true;
        } else if (key == "phi") {
            phi = stof(value);
        } else if (key == "overdispersion") {
            overdispersion = stof(value);
        } else if (key == "skeletonSnpFile") {
            skeletonSnpFile = value;
        } else if (key == "directPrune") {
            directPrune = true;
        } else if (key == "simu") {
            simuMode = true;
        } else if (key.substr(0,2) == "//" ||
                   key.substr(0,1) == "#") {
            continue;
        } else {
            throw("\nError: invalid option " + key + " " + value + "\n");
        }
        ss << boost::format("%20s %-1s %-20s\n") %key %":" %value;
    }
    in.close();
    
    cout << ss.str() << endl;

    setThread();
}

void Options::makeTitle(void){
    title = optionFile;
    size_t pos = optionFile.rfind('.');
    if (pos != string::npos) {
        title = optionFile.substr(0,pos);
    }
}

void Options::setThread(void){
    omp_set_num_threads(numThread);
    if (numThread == 1) return;
    if (multiThreadEigen) {
        Eigen::initParallel();
        Eigen::setNbThreads(numThread);
        cout << "Eigen library is using " << Eigen::nbThreads( ) << " threads." << endl;
    }
#pragma omp parallel
    printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
}
