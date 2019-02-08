//
//  xci.cpp
//  gctb
//
//  Created by Jian Zeng on 27/10/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "xci.hpp"

using namespace std;

void XCI::sortIndBySex(vector<IndInfo*> &indInfoVec){
    vector<IndInfo*> male, female;
    long numInds = indInfoVec.size();
    IndInfo *ind;
    for (unsigned i=0; i<numInds; ++i) {
        ind = indInfoVec[i];
        if (ind->sex == 1) male.push_back(ind);
        else if (ind->sex == 2) female.push_back(ind);
    }
    indInfoVec.resize(0);
    indInfoVec.reserve(male.size() + female.size());
    indInfoVec.insert(indInfoVec.end(), male.begin(), male.end());
    indInfoVec.insert(indInfoVec.end(), female.begin(), female.end());
}

void XCI::restoreFamFileOrder(vector<IndInfo*> &indInfoVec){
    long numInds = indInfoVec.size();
    vector<IndInfo*> vec(numInds);
    IndInfo *ind;
    for (unsigned i=0; i<numInds; ++i) {
        ind = indInfoVec[i];
        vec[ind->famFileOrder] = ind;
    }
    indInfoVec = vec;
}

void XCI::inputIndInfo(Data &data, const string &bedFile, const string &phenotypeFile, const string &keepIndFile,
                       const unsigned keepIndMax, const unsigned mphen, const string &covariateFile){
    data.readFamFile(bedFile + ".fam");
    data.readPhenotypeFile(phenotypeFile, mphen);
    data.readCovariateFile(covariateFile);
    sortIndBySex(data.indInfoVec);
    data.keepMatchedInd(keepIndFile, keepIndMax);
    
    numKeptMales   = 0;
    numKeptFemales = 0;
    IndInfo *ind;
    for (unsigned i=0; i<data.numKeptInds; ++i){
        ind = data.keptIndInfoVec[i];
        if (ind->sex == 1) ++numKeptMales;
        else if (ind->sex == 2) ++numKeptFemales;
    }
    
    // MPI
//    unsigned numKeptMales_all, numKeptFemales_all;
//    MPI_Allreduce(&numKeptMales, &numKeptMales_all, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
//    MPI_Allreduce(&numKeptFemales, &numKeptFemales_all, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
//   if (myMPI::rank==0)
//        cout << "Matched " << numKeptMales_all << " males and " << numKeptFemales_all << " females." << endl;
    cout << "Matched " << numKeptMales << " males and " << numKeptFemales << " females." << endl;
    
    restoreFamFileOrder(data.indInfoVec);
}

void XCI::inputSnpInfo(Data &data, const string &bedFile, const string &includeSnpFile, const string &excludeSnpFile,
                       const unsigned includeChr, const bool readGenotypes){
    data.readBimFile(bedFile + ".bim");
    if (!includeSnpFile.empty()) data.includeSnp(includeSnpFile);
    if (!excludeSnpFile.empty()) data.excludeSnp(excludeSnpFile);
    data.includeChr(includeChr);
    data.includeMatchedSnp();
    if (readGenotypes) readBedFile(data, bedFile + ".bed");  // XCI method: (1) adjust column mean separately in males and females, (2) compute snp2pq from females only
}

void XCI::readBedFile(Data &data, const string &bedFile){
    // features: (1) adjust column mean separately in males and females, (2) compute snp2pq from females only
    unsigned i = 0, j = 0;
    
    if (data.numIncdSnps == 0) throw ("Error: No SNP is retained for analysis.");
    if (data.numKeptInds == 0) throw ("Error: No individual is retained for analysis.");
    
    data.Z.resize(data.numKeptInds, data.numIncdSnps);
    data.ZPZdiag.resize(data.numIncdSnps);
    data.snp2pq.resize(data.numIncdSnps);
    
    // Read bed file
    FILE *in = fopen(bedFile.c_str(), "rb");
    if (!in) throw ("Error: can not open the file [" + bedFile + "] to read.");
//    if (myMPI::rank==0)
        cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
    char header[3];
    fread(header, sizeof(header), 1, in);
    if (!in || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
        cerr << "Error: Incorrect first three bytes of bed file: " << bedFile << endl;
        exit(1);
    }
    
//    unsigned numKeptMales_all;
//    unsigned numKeptFemales_all;
//    MPI_Allreduce(&numKeptMales, &numKeptMales_all, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
//    MPI_Allreduce(&numKeptFemales, &numKeptFemales_all, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    
    // Read genotypes
    SnpInfo *snpInfo = NULL;
    IndInfo *indInfo = NULL;
    unsigned snp = 0;
    unsigned nmiss_male=0;
//    unsigned nmiss_male_all;
    unsigned nmiss_female=0;
//    unsigned nmiss_female_all;
    float sum_male=0.0, mean_male;
//    float sum_male_all=0.0, mean_male_all;
    float sum_female=0.0, mean_female;
//    flaot sum_female_all=0.0, mean_female_all;
    
    const int bedToGeno[4] = {2, -9, 1, 0};
    unsigned size = (data.numInds+3)>>2;
    int genoValue;
    unsigned long long skip = 0;
    
    for (j = 0, snp = 0; j < data.numSnps; j++) {  // code adopted from BOLT-LMM with modification
        snpInfo = data.snpInfoVec[j];
        sum_male = 0.0;
        sum_female = 0.0;
        nmiss_male = 0;
        nmiss_female = 0;
        
        if (!snpInfo->included) {
            //            in.ignore(size);
            skip += size;
            continue;
        }
        
        if (skip) fseek(in, skip, SEEK_CUR);
        skip = 0;
        
        char *bedLineIn = new char[size];
        fread(bedLineIn, 1, size, in);
        
//        int wrong_coding = 0;
//        
//        set<int> maleGeno;
//        set<int> femaleGeno;
        
        for (i = 0; i < data.numInds; i++) {
            indInfo = data.indInfoVec[i];
            if (!indInfo->kept) continue;
            genoValue = bedToGeno[(bedLineIn[i>>2]>>((i&3)<<1))&3];
            if (indInfo->sex == 1) {
                //cout << "Male " << genoValue << endl;
                //if (genoValue == 1) wrong_coding = 1;
                //else if (genoValue == 2) genoValue = 1;
                if (genoValue == 2) {
                    genoValue = 1;
//                    wrong_coding = 1;
                } else if (genoValue == 1) {
                    genoValue = -9;
                }
                
//                maleGeno.insert(genoValue);
                
                if (genoValue == -9) ++nmiss_male;   // missing genotype
                else sum_male += genoValue;
            } else {
                //cout << "Female " << genoValue << endl;
                if (genoValue == -9) ++nmiss_female;   // missing genotype
                else sum_female += genoValue;
                
//                femaleGeno.insert(genoValue);
                
            }
            data.Z(indInfo->index, snp) = genoValue;
        }
        delete[] bedLineIn;
        
//        if (myMPI::rank==0) if (!(snp % 1000)) cout << "snp " << snp << " male " << data.Z.col(snp).transpose().head(10) << " ... female " << data.Z.col(snp).transpose().tail(10) << " ..." << endl;
//        if (myMPI::rank==0) {
//            cout << "snp " << snp << " male ";
//            for (set<int>::iterator it=maleGeno.begin(); it!=maleGeno.end(); ++it) {
//                cout << *it << " ";
//            }
//            cout << " female ";
//            for (set<int>::iterator it=maleGeno.begin(); it!=maleGeno.end(); ++it) {
//                cout << *it << " ";
//            }
//            cout << endl;
//        }
//        
//        int sum_wrong_coding;
//        MPI_Allreduce(&wrong_coding, &sum_wrong_coding, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
//        if (myMPI::rank==0)
//            if (sum_wrong_coding) cout << "Warning: SNP " << snpInfo->ID << " is coded as 0/2 in male X chromosome!" << endl;

//        MPI_Allreduce(&sum_male, &sum_male_all, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
//        MPI_Allreduce(&nmiss_male, &nmiss_male_all, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
//        MPI_Allreduce(&sum_female, &sum_female_all, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
//        MPI_Allreduce(&nmiss_female, &nmiss_female_all, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
        
        // fill missing values with the mean
//        mean_male_all = sum_male_all/float(numKeptMales_all - nmiss_male_all);
//        mean_female_all = sum_female_all/float(numKeptFemales_all - nmiss_female_all);
        mean_male = sum_male/float(numKeptMales - nmiss_male);
        mean_female = sum_female/float(numKeptFemales - nmiss_female);
        if (nmiss_male) {
            for (i=0; i<numKeptMales; ++i) {
//                if (data.Z(i,snp) == -9) data.Z(i,snp) = mean_male_all;
                if (data.Z(i,snp) == -9) data.Z(i,snp) = mean_male;
            }
        }
        if (nmiss_female) {
            for (i=numKeptMales; i<data.numKeptInds; ++i) {
//                if (data.Z(i,snp) == -9) data.Z(i,snp) = mean_female_all;
                if (data.Z(i,snp) == -9) data.Z(i,snp) = mean_female;
            }
        }
        
        // compute allele frequency
//        snpInfo->af = 0.5f*mean_female_all;
        snpInfo->af = 0.5f*mean_female;
        data.snp2pq[snp] = 2.0f*snpInfo->af*(1.0f-snpInfo->af);
        
        //cout << "snp " << snp << "     " << Z.col(snp).sum() << endl;
        
//        data.Z.col(snp).head(numKeptMales).array() -= mean_male_all; // center column by 2p rather than the real mean
//        data.Z.col(snp).tail(numKeptFemales).array() -= mean_female_all; // center column by 2p rather than the real mean
        data.Z.col(snp).head(numKeptMales).array() -= mean_male; // center column by 2p rather than the real mean
        data.Z.col(snp).tail(numKeptFemales).array() -= mean_female; // center column by 2p rather than the real mean

        if (++snp == data.numIncdSnps) break;
    }
    fclose(in);
    
    
    // standardize genotypes
//    VectorXf my_ZPZdiag = data.Z.colwise().squaredNorm();
//    
//    MPI_Allreduce(&my_ZPZdiag[0], &data.ZPZdiag[0], data.numIncdSnps, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    
//    if (myMPI::rank==0)
//        cout << "Genotype data for " << numKeptMales_all + numKeptFemales_all << " individuals (" << numKeptMales_all << " males and " << numKeptFemales_all << " females) and " << data.numIncdSnps << " SNPs are included from [" + bedFile + "]." << endl;
    cout << "Genotype data for " << numKeptMales + numKeptFemales << " individuals (" << numKeptMales << " males and " << numKeptFemales << " females) and " << data.numIncdSnps << " SNPs are included from [" + bedFile + "]." << endl;

}


Model* XCI::buildModel(Data &data, const string &bayesType, const float heritability, const float pi, const float piAlpha, const float piBeta, const bool estimatePi, const float piNDC){
    data.initVariances(heritability);
    if (bayesType == "B") {
        return new BayesBXCI(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, piNDC, numKeptMales, numKeptFemales);
    }
    if (bayesType == "C") {
        return new BayesCXCI(data, data.varGenotypic, data.varResidual, pi, piAlpha, piBeta, estimatePi, piNDC, numKeptMales, numKeptFemales);
    }
    else {
        throw(" Error: Wrong bayes type: " + bayesType);
    }

}

void XCI::simu(Data &data, const unsigned numQTL, const float heritability, const float probNDC, const bool removeQTL, const string &title){
    vector<unsigned> indices(data.numIncdSnps);
    std::iota(indices.begin(), indices.end(), 0);
    std::random_shuffle(indices.begin(), indices.end());
    vector<SnpInfo*> QTLvec(numQTL);
    MatrixXf Q(data.numKeptInds, numQTL);
    VectorXf alpha(numQTL);
    VectorXf isNDC(numQTL);
    SnpInfo *qtl;
    unsigned numFDC = 0;
    for (unsigned j=0; j<numQTL; ++j) {
        qtl = data.incdSnpInfoVec[indices[j]];
        qtl->isQTL = true;
        QTLvec[j] = qtl;
        Q.col(j) = data.Z.col(qtl->index);
        alpha[j] = Stat::snorm();
        if (Stat::ranf() < 1.0f - probNDC) {
            for (unsigned i=numKeptMales; i<data.numKeptInds; ++i) {
                if ( Q(i,j) == 1 ) Q(i,j) = Stat::ranf() < 0.5 ? 1 : 0;
                else Q(i,j) *= 0.5f;
            }
            //Q.col(j).tail(numKeptFemales) *= 0.5f;
            isNDC[j] = 0;
            ++numFDC;
        } else {
            isNDC[j] = 1;
        }
    }
    
    VectorXf g = Q*alpha;
    
    // calculate genetic variance with MPI
    float sumg = g.sum();
    float ssg = g.squaredNorm();
    unsigned ng = data.numKeptInds;
//    float my_sumg = g.sum();
//    float my_ssg  = g.squaredNorm();
//    float sumg, ssg;
//    unsigned ng;
//    MPI_Allreduce(&my_sumg, &sumg, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
//    MPI_Allreduce(&my_ssg, &ssg, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
//    MPI_Allreduce(&data.numKeptInds, &ng, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    
    float genVar = ssg/float(ng) - sumg/float(ng)*sumg/float(ng);
    float resVar = genVar*(1.0-heritability)/heritability;
    float resSD  = sqrt(resVar);
    
//    my_sumg = g.head(numKeptMales).sum();
//    my_ssg  = g.head(numKeptMales).squaredNorm();
//    MPI_Allreduce(&my_sumg, &sumg, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
//    MPI_Allreduce(&my_ssg, &ssg, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
//    MPI_Allreduce(&numKeptMales, &ng, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    sumg = g.head(numKeptMales).sum();
    ssg = g.head(numKeptMales).squaredNorm();
    float genVarMale = ssg/float(ng) - sumg*sumg/float(ng*ng);

//    my_sumg = g.tail(numKeptFemales).sum();
//    my_ssg  = g.tail(numKeptFemales).squaredNorm();
//    MPI_Allreduce(&my_sumg, &sumg, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
//    MPI_Allreduce(&my_ssg, &ssg, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
//    MPI_Allreduce(&numKeptFemales, &ng, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    sumg = g.tail(numKeptFemales).sum();
    ssg  = g.tail(numKeptFemales).squaredNorm();
    float genVarFemale = ssg/float(ng) - sumg*sumg/float(ng*ng);

    for (unsigned i=0; i<data.numKeptInds; ++i) {
        data.y[i] = g[i] + Stat::snorm()*resSD;
    }
//    float my_ypy = (data.y.array()-data.y.mean()).square().sum();
    data.ypy = (data.y.array()-data.y.mean()).square().sum();
    
//    MPI_Allreduce(&my_ypy, &data.ypy, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    
    data.varGenotypic = genVar;
    data.varResidual  = resVar;
    
    if (removeQTL) {
        MatrixXf Ztmp(data.numKeptInds, data.numIncdSnps-numQTL);
        VectorXf ZPZdiagTmp(data.numIncdSnps-numQTL);
        VectorXf snp2pqTmp(data.numIncdSnps-numQTL);
        vector<string> snpEffectNamesTmp(data.numIncdSnps-numQTL);
        for (unsigned j=0, k=0; j<data.numIncdSnps; ++j) {
            if (data.incdSnpInfoVec[j]->isQTL) continue;
            Ztmp.col(k) = data.Z.col(j);
            ZPZdiagTmp[k] = data.ZPZdiag[j];
            snp2pqTmp[k] = data.snp2pq[j];
            snpEffectNamesTmp[k] = data.snpEffectNames[j];
            ++k;
        }
        data.Z = Ztmp;
        data.ZPZdiag = ZPZdiagTmp;
        data.snp2pq = snp2pqTmp;
        data.snpEffectNames = snpEffectNamesTmp;
        data.numIncdSnps -= numQTL;
    }
    
    string filename = title + ".QTLinfo";
    ofstream out(filename.c_str());
    out << boost::format("%6s %20s %6s %12s %8s %12s %6s\n")
    % "Id"
    % "Name"
    % "Chrom"
    % "Position"
    % "GeneFrq"
    % "Effect"
    % "EscapeXCI";
    for (unsigned j=0; j<numQTL; ++j) {
        qtl = data.incdSnpInfoVec[indices[j]];
        out << boost::format("%6s %20s %6s %12s %8s %12s %6s\n")
        % (j+1)
        % qtl->ID
        % qtl->chrom
        % qtl->physPos
        % qtl->af
        % alpha[j]
        % isNDC[j];
    }
    out.close();
    
//    if (!myMPI::rank) {
        cout << "\nSimulated " << numQTL << " QTL with " << numQTL - numFDC << " QTL escaped from XCI." << endl;
        cout << "Simulated genotypic variance: " << genVar << " (male: " << genVarMale << "; female: " << genVarFemale << ")" <<  endl;
        cout << "Simulated residual  variance: " << resVar << endl;
        if (removeQTL) cout << "QTL removed from the analysis." << endl;
        cout << "Saved simulated QTL info to [" << filename << "]." << endl;
//    }
    
}

void XCI::outputResults(const Data &data, const vector<McmcSamples*> &mcmcSampleVec, const string &title){
    McmcSamples *snpEffects = NULL;
    McmcSamples *gamma = NULL;
    for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
        if (mcmcSampleVec[i]->label == "SnpEffects") snpEffects = mcmcSampleVec[i];
        if (mcmcSampleVec[i]->label == "Gamma") gamma = mcmcSampleVec[i];
    }
//    if (myMPI::rank) return;
    string filename = title + ".snpRes";
    ofstream out(filename.c_str());
    out << boost::format("%6s %20s %6s %12s %8s %12s %12s %8s %8s\n")
    % "Id"
    % "Name"
    % "Chrom"
    % "Position"
    % "GeneFrq"
    % "Effect"
    % "SE"
    % "PIP"
    % "PrNDC";
    for (unsigned i=0, idx=0; i<data.numSnps; ++i) {
        SnpInfo *snp = data.snpInfoVec[i];
        if(!data.fullSnpFlag[i]) continue;
        if(snp->isQTL) continue;
        out << boost::format("%6s %20s %6s %12s %8.3f %12.6f %12.6f %8.3f %8.3f\n")
        % (idx+1)
        % snp->ID
        % snp->chrom
        % snp->physPos
        % snp->af
        % snpEffects->posteriorMean[idx]
        % sqrt(snpEffects->posteriorSqrMean[idx]-snpEffects->posteriorMean[idx]*snpEffects->posteriorMean[idx])
        % snpEffects->pip[idx]
        % gamma->posteriorMean[idx];
        ++idx;
    }
    out.close();
}


void BayesCXCI::SnpEffects::sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z, const VectorXf &ZPZdiag,
                                         const VectorXf &ZPZdiagMale, const VectorXf &ZPZdiagFemale,
                                         const unsigned nmale, const unsigned nfemale, const float p,
                                         const float sigmaSq, const float pi, const float varem, const float varef,
                                         VectorXf &gamma, VectorXf &ghatm, VectorXf &ghatf){
    // sample beta, delta, gamma jointly
    // f(beta, delta, gamma) propto f(beta | delta, gamma) f(delta | gamma) f(gamma)
    
    //cout << "check pi " << pi << " p " << p << " sigmaSq " << sigmaSq << endl;
    
    sumSq = 0.0;
    numNonZeros = 0;
    
    ghatm.setZero(ycorrm.size());
    ghatf.setZero(ycorrf.size());

    Vector2f rhsFemale;     // 0: FDC, 1: NDC, corresponding to gamma
    Vector2f my_rhs, rhs;   // 0: FDC, 1: NDC, corresponding to gamma
    Vector2f invLhs;        // 0: FDC, 1: NDC, corresponding to gamma
    Vector2f uhat;          // 0: FDC, 1: NDC, corresponding to gamma
    Vector2f logGamma;      // 0: FDC, 1: NDC, corresponding to gamma
    
    float oldSample, sample;
    float logDelta0, logDelta1, probDelta1;
    float logPi = log(pi);
    float logPiComp = log(1.0-pi);
    float logSigmaSq = log(sigmaSq);
    float invVarem = 1.0f/varem;
    float invVaref = 1.0f/varef;
    float invSigmaSq = 1.0f/sigmaSq;
    float logP = log(p);
    float logPcomp = log(1.0f-p);
    float rhsMale;
    float probGamma1;
    float sampleGamma;
    
    for (unsigned i=0; i<size; ++i) {
        oldSample = values[i];
        ycorrm += Z.col(i).head(nmale) * oldSample;
        if (gamma[i])
            ycorrf += Z.col(i).tail(nfemale) * oldSample;
        else
            ycorrf += Z.col(i).tail(nfemale) * oldSample * 0.5f;
        
        rhsMale = Z.col(i).head(nmale).dot(ycorrm) * invVarem;
        rhsFemale[1] = Z.col(i).tail(nfemale).dot(ycorrf) * invVaref;
        rhsFemale[0] = rhsFemale[1] * 0.5f;
        
        rhs[1] = rhsMale + rhsFemale[1];
        rhs[0] = rhsMale + rhsFemale[0];

//        my_rhs[1] = rhsMale + rhsFemale[1];
//        my_rhs[0] = rhsMale + rhsFemale[0];
//        
//        MPI_Allreduce(&my_rhs[0], &rhs[0], 2, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

        invLhs[1] = 1.0f/(ZPZdiagMale[i]*invVarem + ZPZdiagFemale[i]*invVaref       + invSigmaSq);
        invLhs[0] = 1.0f/(ZPZdiagMale[i]*invVarem + ZPZdiagFemale[i]*0.25f*invVaref + invSigmaSq);
        
        uhat.array() = invLhs.array()*rhs.array();
        
        //sample gamma
        
        logGamma[1] = 0.5f*rhs[1]*uhat[1] + logf(sqrt(invLhs[1])*pi + expf(0.5f*(logSigmaSq-rhs[1]*uhat[1]))*(1.0f-pi)) + logP;
        logGamma[0] = 0.5f*rhs[0]*uhat[0] + logf(sqrt(invLhs[0])*pi + expf(0.5f*(logSigmaSq-rhs[0]*uhat[0]))*(1.0f-pi)) + logPcomp;
        probGamma1 = 1.0f/(1.0f + expf(logGamma[0] - logGamma[1]));
        
        sampleGamma = bernoulli.sample(probGamma1);
        gamma[i] = sampleGamma;
        
        
        // sample delta
        
        logDelta1 = 0.5*(logf(invLhs[sampleGamma]) + uhat[sampleGamma]*rhs[sampleGamma]) + logPi;
        logDelta0 = 0.5*logSigmaSq + logPiComp;
        probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
        
        if (bernoulli.sample(probDelta1)) {
            
            // sample effect
            
            sample = normal.sample(uhat[sampleGamma], invLhs[sampleGamma]);
            values[i] = sample;
            sumSq += sample * sample;
            ++numNonZeros;
            
            ycorrm -= Z.col(i).head(nmale) * sample;
            ghatm += Z.col(i).head(nmale) * sample;
            
            if (gamma[i]) {
                ycorrf -= Z.col(i).tail(nfemale) * sample;
                ghatf += Z.col(i).tail(nfemale) * sample;
            } else {
                ycorrf -= Z.col(i).tail(nfemale) * sample * 0.5f;
                ghatf += Z.col(i).tail(nfemale) * sample * 0.5f;
            }
        }
        else {
            values[i] = 0.0;
        }
    }
}

//void BayesCXCI::SnpEffects::sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z, const VectorXf &ZPZdiag,
//                                         const VectorXf &ZPZdiagMale, const VectorXf &ZPZdiagFemale,
//                                         const unsigned nmale, const unsigned nfemale, const float p,
//                                         const float sigmaSq, const float pi, const float varem, const float varef,
//                                         VectorXf &gamma, VectorXf &ghatm, VectorXf &ghatf){
//    // sample beta, delta, gamma jointly
//    // f(beta, delta, gamma) propto f(beta | delta, gamma) f(delta | gamma) f(gamma)
//    
//    //cout << "check pi " << pi << " p " << p << " sigmaSq " << sigmaSq << endl;
//    
//    sumSq = 0.0;
//    numNonZeros = 0;
//    
//    ghatm.setZero(ycorrm.size());
//    ghatf.setZero(ycorrf.size());
//    
//    Vector2f rhsFemale;     // 0: FDC, 1: NDC, corresponding to gamma
//    Vector2f my_rhs, rhs;   // 0: FDC, 1: NDC, corresponding to gamma
//    Vector2f invLhs;        // 0: FDC, 1: NDC, corresponding to gamma
//    Vector2f uhat;          // 0: FDC, 1: NDC, corresponding to gamma
//    Vector2f logGamma;      // 0: FDC, 1: NDC, corresponding to gamma
//    
//    float oldSample, sample;
//    float oldGamma;
//    float logDelta0, logDelta1, probDelta1;
//    float logPi = log(pi);
//    float logPiComp = log(1.0-pi);
//    float logSigmaSq = log(sigmaSq);
//    float invVarem = 1.0f/varem;
//    float invVaref = 1.0f/varef;
//    float invSigmaSq = 1.0f/sigmaSq;
//    float logP = log(p);
//    float logPcomp = log(1.0f-p);
//    float rhsMale;
//    float probGamma1;
//    float sampleGamma;
//    
//    for (unsigned i=0; i<size; ++i) {
//        oldSample = values[i];
//        oldGamma = gamma[i];
//
//                ycorrm += Z.col(i).head(nmale) * oldSample;
//                if (gamma[i])
//                    ycorrf += Z.col(i).tail(nfemale) * oldSample;
//                else
//                    ycorrf += Z.col(i).tail(nfemale) * oldSample * 0.5f;
//
//        rhsMale = Z.col(i).head(nmale).dot(ycorrm) * invVarem;
//        rhsFemale[1] = Z.col(i).tail(nfemale).dot(ycorrf) * invVaref;
//        rhsFemale[0] = rhsFemale[1] * 0.5f;
//        
//        my_rhs[1] = rhsMale + rhsFemale[1];
//        my_rhs[0] = rhsMale + rhsFemale[0];
//        
//        MPI_Allreduce(&my_rhs[0], &rhs[0], 2, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
//        
//        invLhs[1] = 1.0f/(ZPZdiagMale[i]*invVarem + ZPZdiagFemale[i]*invVaref       + invSigmaSq);
//        invLhs[0] = 1.0f/(ZPZdiagMale[i]*invVarem + ZPZdiagFemale[i]*0.25f*invVaref + invSigmaSq);
//        
//        uhat.array() = invLhs.array()*rhs.array();
//        
//        //sample gamma
//        
//        logGamma[1] = 0.5f*rhs[1]*uhat[1] + logf(sqrt(invLhs[1])*pi + expf(0.5f*(logSigmaSq-rhs[1]*uhat[1]))*(1.0f-pi)) + logP;
//        logGamma[0] = 0.5f*rhs[0]*uhat[0] + logf(sqrt(invLhs[0])*pi + expf(0.5f*(logSigmaSq-rhs[0]*uhat[0]))*(1.0f-pi)) + logPcomp;
//        probGamma1 = 1.0f/(1.0f + expf(logGamma[0] - logGamma[1]));
//        
//        sampleGamma = bernoulli.sample(probGamma1);
//        gamma[i] = sampleGamma;
//        
//        
//        // sample delta
//        
//        logDelta1 = 0.5*(logf(invLhs[sampleGamma]) + uhat[sampleGamma]*rhs[sampleGamma]) + logPi;
//        logDelta0 = 0.5*logSigmaSq + logPiComp;
//        probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
//        
//        if (bernoulli.sample(probDelta1)) {
//            
//            // sample effect
//            
//            sample = normal.sample(uhat[sampleGamma], invLhs[sampleGamma]);
//            values[i] = sample;
//            sumSq += sample * sample;
//            ++numNonZeros;
//            
//            ycorrm -= Z.col(i).head(nmale) * sample;
////            ycorrm += Z.col(i).head(nmale) * (oldSample - sample);
//            ghatm += Z.col(i).head(nmale) * sample;
//            
//            if (gamma[i]) {
////                if (oldGamma)
////                    ycorrf += Z.col(i).tail(nfemale) * (oldSample - sample);
////                else
////                    ycorrf += Z.col(i).tail(nfemale) * (0.5f*oldSample - sample);
//
//                ycorrf -= Z.col(i).tail(nfemale) * sample;
//                ghatf += Z.col(i).tail(nfemale) * sample;
//            } else {
////                if (oldGamma)
////                    ycorrf += Z.col(i).tail(nfemale) * (oldSample - 0.5f*sample);
////                else
////                    ycorrf += Z.col(i).tail(nfemale) * 0.5f*(oldSample - sample);
//                
//                ycorrf -= Z.col(i).tail(nfemale) * sample * 0.5f;
//                ghatf += Z.col(i).tail(nfemale) * sample * 0.5f;
//            }
//        }
//        else {
//            if (oldSample) {
////                ycorrm += Z.col(i).head(nmale) * oldSample;
////                if (gamma[i]) {
////                    if (oldGamma)
////                        ycorrf += Z.col(i).tail(nfemale) * oldSample;
////                    else
////                        ycorrf -= Z.col(i).tail(nfemale) * 0.5f*oldSample;
////                } else {
////                    if (oldGamma)
////                        ycorrf += Z.col(i).tail(nfemale) * 0.5f*oldSample;
////                    else
////                        ycorrf += Z.col(i).tail(nfemale) * 0.5f*oldSample;
////                }
//            }
//            values[i] = 0.0;
//        }
//    }
//}

//void BayesCXCI::SnpEffects::sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z, const VectorXf &ZPZdiag,
//                                         const VectorXf &ZPZdiagMale, const VectorXf &ZPZdiagFemale,
//                                         const unsigned nmale, const unsigned nfemale, const float p,
//                                         const float sigmaSq, const float pi, const float varem, const float varef,
//                                         VectorXf &gamma, VectorXf &ghatm, VectorXf &ghatf){
//    // sample beta, delta, gamma from their full conditionals
//    
//    //cout << "check pi " << pi << " p " << p << " sigmaSq " << sigmaSq << endl;
//    
//    sumSq = 0.0;
//    numNonZeros = 0;
//    
//    ghatm.setZero(ycorrm.size());
//    ghatf.setZero(ycorrf.size());
//    
//    float oldSample, oldGamma;
//    float logDelta0minus1, probDelta1;
//    float logGamma0minus1, probGamma1;
//    float logPi = log(pi);
//    float logPiComp = log(1.0-pi);
//    float invVarem = 1.0f/varem;
//    float invVaref = 1.0f/varef;
//    float invSigmaSq = 1.0f/sigmaSq;
//    float logP = log(p);
//    float logPcomp = log(1.0f-p);
//    float rhsMale, rhsFemale, rhs;
//    float lhsMale, lhsFemale, lhs;
//    float invLhs, uhat;
//    
//    for (unsigned i=0; i<size; ++i) {
//        oldSample = values[i];
//        oldGamma = gamma[i];
//        
//        ycorrm += Z.col(i).head(nmale) * oldSample;
//        if (gamma[i])
//            ycorrf += Z.col(i).tail(nfemale) * oldSample;
//        else
//            ycorrf += Z.col(i).tail(nfemale) * oldSample * 0.5f;
//       
//        rhsMale = Z.col(i).head(nmale).dot(ycorrm) * invVarem;
//        rhsFemale = Z.col(i).tail(nfemale).dot(ycorrf) * invVaref;
//        lhsMale = ZPZdiagMale[i]*invVarem;
//        lhsFemale = ZPZdiagFemale[i]*invVaref;
//        if (gamma[i]) {
//            rhs = rhsMale + rhsFemale;
//            lhs = lhsMale + lhsFemale;
//        } else {
//            rhs = rhsMale + 0.5f*rhsFemale;
//            lhs = lhsMale + 0.25f*lhsFemale;
//        }
//        
//        // sample delta
//        
//        logDelta0minus1 = 2.0f*rhs*values[i] - lhs*values[i]*values[i] + logPiComp - logPi;
//        probDelta1 = 1.0f/(1.0f + expf(logDelta0minus1));
//        
//        if (bernoulli.sample(probDelta1)) {
//            
//            // sample gamma
//            
//            
//            logGamma0minus1 = rhsFemale*values[i] - 0.75f*lhsFemale*values[i]*values[i] + logPcomp - logP;
//            probGamma1 = 1.0f/(1.0f + expf(logGamma0minus1));
//    
//            if (bernoulli.sample(probGamma1)) {
//                
//                // sample beta
//               
//                rhs = rhsMale + rhsFemale;
//                invLhs = 1.0f/(lhsMale + lhsFemale + invSigmaSq);
//                uhat = invLhs*rhs;
//                values[i] = normal.sample(uhat, invLhs);
//                gamma[i] = 1.0;
//                
////                cout << i << " male " << rhsMale/(lhsMale + invSigmaSq) << " female " << rhsFemale/(lhsFemale + invSigmaSq) << endl;
//                
////                if (oldGamma)
////                    ycorrf += Z.col(i).tail(nfemale) * (oldSample - values[i]);
////                else
////                    ycorrf += Z.col(i).tail(nfemale) * (0.5f*oldSample - values[i]);
//                
//                ycorrf -= Z.col(i).tail(nfemale) * values[i];
//                
//                ghatf += Z.col(i).tail(nfemale) * values[i];
//                
//            } else {
//                
//                rhs = rhsMale + 0.5f*rhsFemale;
//                invLhs = 1.0f/(lhsMale + 0.25f*lhsFemale + invSigmaSq);
//                uhat = invLhs*rhs;
//                values[i] = normal.sample(uhat, invLhs);
//                gamma[i] = 0.0;
// 
////                cout << i << " male " << rhsMale/(lhsMale + invSigmaSq) << " female " << rhsFemale/(lhsFemale + invSigmaSq) << endl;
//
////                if (oldGamma)
////                    ycorrf += Z.col(i).tail(nfemale) * (oldSample - 0.5f*values[i]);
////                else
////                    ycorrf += Z.col(i).tail(nfemale) * 0.5f*(oldSample - values[i]);
//                
//                ycorrf -= Z.col(i).tail(nfemale) * values[i] * 0.5f;
//
//                ghatf += Z.col(i).tail(nfemale) * 0.5f*values[i];
//
//            }
//            
////            ycorrm += Z.col(i).head(nmale) * (oldSample - values[i]);
//            ycorrm -= Z.col(i).head(nmale) * values[i];
//            
//            ghatm += Z.col(i).head(nmale) * values[i];
//            sumSq += values[i] * values[i];
//            ++numNonZeros;
//            
//        } else {
//            
////            if (oldSample) {
////                ycorrm += Z.col(i).head(nmale) * oldSample;
////                if (oldGamma)
////                    ycorrf += Z.col(i).tail(nfemale) * oldSample;
////                else
////                    ycorrf += Z.col(i).tail(nfemale) * 0.5f*oldSample;
////            }
//            gamma[i] = bernoulli.sample(p);
//            values[i] = 0.0;
//        }
//    }
//    
//}

void BayesCXCI::ProbNDC::sampleFromFC(const unsigned numSnps, const unsigned numNDC){
    //cout << numSnps << " " << numNDC << endl;
    float alphaTilde = numNDC + alpha;
    float betaTilde  = numSnps - numNDC + beta;
    value = Beta::sample(alphaTilde, betaTilde);
}

void BayesCXCI::Rounding::computeYcorr(const VectorXf &y, const MatrixXf &X, const MatrixXf &Z,
                                           const VectorXf &gamma, const unsigned int nmale, const unsigned int nfemale,
                                           const VectorXf &fixedEffects, const VectorXf &snpEffects, VectorXf &ycorrm, VectorXf &ycorrf){
    if (count++ % 100) return;
    VectorXf oldYcorrm = ycorrm;
    VectorXf oldYcorrf = ycorrf;
    VectorXf ycorr = y - X*fixedEffects;
    for (unsigned i=0; i<snpEffects.size(); ++i) {
        if (snpEffects[i]) {
            if (gamma[i]) {
                ycorr -= Z.col(i)*snpEffects[i];
            } else {
                ycorr.head(nmale) -= Z.col(i).head(nmale)*snpEffects[i];
                ycorr.tail(nfemale) -= Z.col(i).tail(nfemale)*snpEffects[i]*0.5f;
            }
        }
    }
    ycorrm = ycorr.head(nmale);
    ycorrf = ycorr.tail(nfemale);
    float ss = (ycorrm - oldYcorrm).squaredNorm() + (ycorrf - oldYcorrf).squaredNorm();
//    float my_ss = (ycorrm - oldYcorrm).squaredNorm() + (ycorrf - oldYcorrf).squaredNorm();
//    float ss;
//    MPI_Allreduce(&my_ss, &ss, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    value = sqrt(ss);
}

void BayesCXCI::sampleUnknowns(){
    ycorr.head(nmale) = ycorrm;
    ycorr.tail(nfemale) = ycorrf;
    
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    
    ycorrm = ycorr.head(nmale);
    ycorrf = ycorr.tail(nfemale);
    
    unsigned cnt=0;
    do {
        snpEffects.sampleFromFC(ycorrm, ycorrf, data.Z, data.ZPZdiag, ZPZdiagMale, ZPZdiagFemale,
                                nmale, nfemale, piNDC.value, sigmaSq.value, pi.value, varem.value, varef.value, gamma.values, ghatm, ghatf);
        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros == 0);
    
    piNDC.sampleFromFC(snpEffects.size, gamma.values.sum());
    sigmaSq.sampleFromFC(snpEffects.sumSq, snpEffects.numNonZeros);
    if(estimatePi) pi.sampleFromFC(snpEffects.size, snpEffects.numNonZeros);
    
    varem.sampleFromFC(ycorrm);
    varef.sampleFromFC(ycorrf);
    vargm.compute(ghatm);
    vargf.compute(ghatf);
    hsqm.compute(vargm.value, varem.value);
    hsqf.compute(vargf.value, varef.value);
    
//    vare.sampleFromFC(ycorr);
//    varg.compute(ghat);
//    hsq.compute(varg.value, vare.value);
    
    rounding.computeYcorr(data.y, data.X, data.Z, gamma.values, nmale, nfemale, fixedEffects.values, snpEffects.values, ycorrm, ycorrf);
    nnzSnp.getValue(snpEffects.numNonZeros);
    
//    static unsigned iter = 0;
//    if (++iter < 5000) {
//        genVarPrior += (varg.value - genVarPrior)/iter;
//        piPrior += (pi.value - piPrior)/iter;
//        scale.compute(genVarPrior, piPrior, sigmaSq.scale);
//    }
//    scale.compute(varg.value, pi.value, sigmaSq.scale);
}


void BayesBXCI::SnpEffects::sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag,
                                         const VectorXf &ZPZdiagMale, const VectorXf &ZPZdiagFemale,
                                         const unsigned int nmale, const unsigned int nfemale,
                                         const float p, const VectorXf &sigmaSq, const float pi,
                                         const float vare, VectorXf &gamma, VectorXf &ghat) {
    numNonZeros = 0;
    
    ghat.setZero(ycorr.size());
    
    Vector2f rhsFemale;     // 0: FDC, 1: NDC, corresponding to gamma
    Vector2f my_rhs, rhs;   // 0: FDC, 1: NDC, corresponding to gamma
    Vector2f invLhs;        // 0: FDC, 1: NDC, corresponding to gamma
    Vector2f uhat;          // 0: FDC, 1: NDC, corresponding to gamma
    Vector2f logGamma;      // 0: FDC, 1: NDC, corresponding to gamma
    
    float oldSample, sample;
    float logDelta0, logDelta1, probDelta1;
    float logPi = log(pi);
    float logPiComp = log(1.0-pi);
    float invVare = 1.0f/vare;
    float logP = log(p);
    float logPcomp = log(1.0f-p);
    float rhsMale;
    float probGamma1;
    float sampleGamma;
    float beta;
    
    for (unsigned i=0; i<size; ++i) {
        oldSample = values[i];
        ycorr.head(nmale) += Z.col(i).head(nmale) * oldSample;
        if (gamma[i])
            ycorr.tail(nfemale) += Z.col(i).tail(nfemale) * oldSample;
        else
            ycorr.tail(nfemale) += Z.col(i).tail(nfemale) * oldSample * 0.5f;
        
        rhsMale = Z.col(i).head(nmale).dot(ycorr.head(nmale)) * invVare;
        rhsFemale[1] = Z.col(i).tail(nfemale).dot(ycorr.tail(nfemale)) * invVare;
        rhsFemale[0] = rhsFemale[1] * 0.5f;

        rhs[1] = rhsMale + rhsFemale[1];
        rhs[0] = rhsMale + rhsFemale[0];
        
//        my_rhs[1] = rhsMale + rhsFemale[1];
//        my_rhs[0] = rhsMale + rhsFemale[0];
//        
//        MPI_Allreduce(&my_rhs[0], &rhs[0], 2, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        
        invLhs[1] = 1.0f/((ZPZdiagMale[i] + ZPZdiagFemale[i]      )*invVare + 1.0f/sigmaSq[i]);
        invLhs[0] = 1.0f/((ZPZdiagMale[i] + ZPZdiagFemale[i]*0.25f)*invVare + 1.0f/sigmaSq[i]);
        
        uhat.array() = invLhs.array()*rhs.array();
        
        //sample gamma
        
        logGamma[1] = 0.5f*rhs[1]*uhat[1] + logf(sqrt(invLhs[1])*pi + expf(0.5f*(logf(sigmaSq[i])-rhs[1]*uhat[1]))*(1.0f-pi)) + logP;
        logGamma[0] = 0.5f*rhs[0]*uhat[0] + logf(sqrt(invLhs[0])*pi + expf(0.5f*(logf(sigmaSq[i])-rhs[0]*uhat[0]))*(1.0f-pi)) + logPcomp;
        probGamma1 = 1.0f/(1.0f + expf(logGamma[0] - logGamma[1]));
        sampleGamma = bernoulli.sample(probGamma1);
        gamma[i] = sampleGamma;
        
        
        // sample delta
        
        logDelta1 = 0.5*(logf(invLhs[sampleGamma]) + uhat[sampleGamma]*rhs[sampleGamma]) + logPi;
        logDelta0 = 0.5*logf(sigmaSq[i]) + logPiComp;
        probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
        
        if (bernoulli.sample(probDelta1)) {
            
            // sample effect
            
            sample = normal.sample(uhat[sampleGamma], invLhs[sampleGamma]);
            values[i] = sample;
            betaSq[i] = sample * sample;
            ++numNonZeros;
            
            ycorr.head(nmale) -= Z.col(i).head(nmale) * sample;
            ghat .head(nmale) += Z.col(i).head(nmale) * sample;
            
            if (gamma[i]) {
                ycorr.tail(nfemale) -= Z.col(i).tail(nfemale) * sample;
                ghat .tail(nfemale) += Z.col(i).tail(nfemale) * sample;
            } else {
                ycorr.tail(nfemale) -= Z.col(i).tail(nfemale) * sample * 0.5f;
                ghat .tail(nfemale) += Z.col(i).tail(nfemale) * sample * 0.5f;
            }
        }
        else {
            beta = normal.sample(0, sigmaSq[i]);
            betaSq[i] = beta*beta;
            values[i] = 0.0;
        }
    }
}

void BayesBXCI::sampleUnknowns(){
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    unsigned cnt=0;
    do {
        snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, ZPZdiagMale, ZPZdiagFemale,
                                nmale, nfemale, piNDC.value, sigmaSq.values, pi.value, vare.value, gamma.values, ghat);
        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros == 0);
    piNDC.sampleFromFC(snpEffects.size, gamma.values.sum());
    sigmaSq.sampleFromFC(snpEffects.betaSq);
    if(estimatePi) pi.sampleFromFC(snpEffects.size, snpEffects.numNonZeros);
    vare.sampleFromFC(ycorr);
    
    varg.compute(ghat);
    hsq.compute(varg.value, vare.value);
    
    rounding.computeYcorr(data.y, data.X, data.Z, gamma.values, nmale, nfemale, fixedEffects.values, snpEffects.values, ycorrm, ycorrf);
    nnzSnp.getValue(snpEffects.numNonZeros);
    
//    static unsigned iter = 0;
//    if (++iter < 5000) {
//        genVarPrior += (varg.value - genVarPrior)/iter;
//        piPrior += (pi.value - piPrior)/iter;
//        scale.compute(genVarPrior, piPrior, sigmaSq.scale);
//    }
    scale.compute(varg.value, pi.value, sigmaSq.scale);
}

