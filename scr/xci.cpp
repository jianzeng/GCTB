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
        else {
            cout << "Error: Individual " + ind->famID + " " + ind->indID + " have sex " + to_string(static_cast<long long>(ind->sex)) << endl;
        }
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
                       const unsigned keepIndMax, const unsigned mphen, const string &covariateFile, const bool femaleOnly){
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
    unsigned numKeptMales_all, numKeptFemales_all;
    MPI_Allreduce(&numKeptMales, &numKeptMales_all, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&numKeptFemales, &numKeptFemales_all, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
   if (myMPI::rank==0)
        cout << "Matched " << numKeptMales_all << " males and " << numKeptFemales_all << " females." << endl;
    
    restoreFamFileOrder(data.indInfoVec);
}

void XCI::inputSnpInfo(Data &data, const string &bedFile, const string &includeSnpFile, const string &excludeSnpFile,
                       const unsigned includeChr, const bool readGenotypes){
    data.readBimFile(bedFile + ".bim");
    if (!includeSnpFile.empty()) data.includeSnp(includeSnpFile);
    if (!excludeSnpFile.empty()) data.excludeSnp(excludeSnpFile);
    data.includeChr(includeChr);
    data.includeMatchedSnp();
    if (readGenotypes) readBedFile(data, bedFile + ".bed");  // XCI method: (1) adjust column mean separately in males and females, (2) compute snp2pq from males only
}

void XCI::readBedFile(Data &data, const string &bedFile){
    // features: (1) adjust column mean separately in males and females, (2) compute snp2pq from males only because the true 2pq in females is affacted by dosage model
    unsigned i = 0, j = 0;
    
    if (data.numIncdSnps == 0) throw ("Error: No SNP is retained for analysis.");
    if (data.numKeptInds == 0) throw ("Error: No individual is retained for analysis.");
    
    data.Z.resize(data.numKeptInds, data.numIncdSnps);
    data.ZPZdiag.resize(data.numIncdSnps);
    data.snp2pq.resize(data.numIncdSnps);
    
    // Read bed file
    FILE *in = fopen(bedFile.c_str(), "rb");
    if (!in) throw ("Error: can not open the file [" + bedFile + "] to read.");
    if (myMPI::rank==0)
        cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
    char header[3];
    fread(header, sizeof(header), 1, in);
    if (!in || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
        cerr << "Error: Incorrect first three bytes of bed file: " << bedFile << endl;
        exit(1);
    }
    
    unsigned numKeptMales_all;
    unsigned numKeptFemales_all;
    MPI_Allreduce(&numKeptMales, &numKeptMales_all, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&numKeptFemales, &numKeptFemales_all, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    
    // Read genotypes
    SnpInfo *snpInfo = NULL;
    IndInfo *indInfo = NULL;
    unsigned snp = 0;
    unsigned nmiss_male=0, nmiss_male_all;
    unsigned nmiss_female=0, nmiss_female_all;
    float sum_male=0.0, sum_male_all=0.0, mean_male_all;
    float sum_female=0.0, sum_female_all=0.0, mean_female_all;
    
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

        MPI_Allreduce(&sum_male, &sum_male_all, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&nmiss_male, &nmiss_male_all, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&sum_female, &sum_female_all, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&nmiss_female, &nmiss_female_all, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
        
        // fill missing values with the mean
        mean_male_all = sum_male_all/float(numKeptMales_all - nmiss_male_all);
        mean_female_all = sum_female_all/float(numKeptFemales_all - nmiss_female_all);
        if (nmiss_male) {
            for (i=0; i<numKeptMales; ++i) {
                if (data.Z(i,snp) == -9) data.Z(i,snp) = mean_male_all;
            }
        }
        if (nmiss_female) {
            for (i=numKeptMales; i<data.numKeptInds; ++i) {
                if (data.Z(i,snp) == -9) data.Z(i,snp) = mean_female_all;
            }
        }
        
        // compute allele frequency
//        if (numKeptMales_all) {
//            snpInfo->af = mean_male_all;
//            data.snp2pq[snp] = snpInfo->af*(1.0f-snpInfo->af);
//        } else {
            snpInfo->af = 0.5f*mean_female_all;
            data.snp2pq[snp] = 2.0f*snpInfo->af*(1.0f-snpInfo->af);
//        }
        
        //cout << "snp " << snp << "     " << Z.col(snp).sum() << endl;
        
        // center columns
        data.Z.col(snp).head(numKeptMales).array() -= mean_male_all; // center column by 2p rather than the real mean
        data.Z.col(snp).tail(numKeptFemales).array() -= mean_female_all; // center column by 2p rather than the real mean
        
        if (++snp == data.numIncdSnps) break;
    }
    fclose(in);
    
    // standardize genotypes
//    VectorXf my_ZPZdiag = data.Z.colwise().squaredNorm();
//    
//    MPI_Allreduce(&my_ZPZdiag[0], &data.ZPZdiag[0], data.numIncdSnps, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    
    if (myMPI::rank==0)
        cout << "Genotype data for " << numKeptMales_all + numKeptFemales_all << " individuals (" << numKeptMales_all << " males and " << numKeptFemales_all << " females) and " << data.numIncdSnps << " SNPs are included from [" + bedFile + "]." << endl;
}


Model* XCI::buildModel(Data &data, const string &bayesType, const float heritability, const float pi, const VectorXf &piPar, const bool estimatePi, const float piNDC, const Vector2f &piNDCpar, const bool estimatePiNDC, const float piGxE, const bool estimatePiGxE, const unsigned windowWidth){
    data.initVariances(heritability);
    bool noscale = true;
    if (bayesType == "B") {
        return new BayesBXCI(data, data.varGenotypic, data.varResidual, pi, piPar[0], piPar[1], estimatePi, piNDC, piNDCpar, estimatePiNDC, numKeptMales, numKeptFemales, noscale);
    }
    if (bayesType == "C") {
        return new BayesCXCI(data, data.varGenotypic, data.varResidual, pi, piPar[0], piPar[1], estimatePi, piNDC, piNDCpar, estimatePiNDC, numKeptMales, numKeptFemales, noscale);
    }
    if (bayesType == "Cgxs") {
        Vector3f pis;
        pis << 1.0-pi, pi*(1.0-piGxE), pi*piGxE;  // prob of zero effect, sex-shared effect, sex-specific effect
        return new BayesCXCIgxs(data, data.varGenotypic, data.varResidual, pis, piPar, estimatePi, piNDC, piNDCpar, estimatePiNDC, estimatePiGxE, numKeptMales, numKeptFemales, noscale);
    }
    if (bayesType == "Ngxs") {
        Vector3f pis;
        pis << 1.0-pi, pi*(1.0-piGxE), pi*piGxE;  // prob of zero effect, sex-shared effect, sex-specific effect
        data.getNonoverlapWindowInfo(windowWidth);
        return new BayesNXCIgxs(data, data.varGenotypic, data.varResidual, pis, piPar, estimatePi, piNDC, piNDCpar, estimatePiNDC, estimatePiGxE, numKeptMales, numKeptFemales, noscale);
    }
    if (bayesType == "Xgxs") {
        Vector3f piBeta, piDosage;
        piBeta << 1.0-pi, pi*(1.0-piGxE), pi*piGxE;  // prob of zero effect, sex-shared effect, sex-specific effect
//        piDosage << 1.0-pi, pi*piNDC, pi*(1.0-piNDC);
        piDosage << 0.33,0.33,0.33;
        return new BayesXgxs(data, data.varGenotypic, data.varResidual, piBeta, piPar, piDosage, piNDCpar, estimatePi, estimatePiGxE, numKeptMales, numKeptFemales, noscale);
    }
    else {
        throw(" Error: Wrong bayes type: " + bayesType);
    }
}

Model* XCI::buildModelStageOne(Data &data, const string &bayesType, const float heritability, const float pi, const VectorXf &piPar, const bool estimatePi, const float piNDC, const Vector2f &piNDCpar, const bool estimatePiNDC){
    if (!myMPI::rank) {
        cout << "Running a two-stage model: " << endl;
        cout << "Stage One: estimating dosage compensation model with females only data." << endl;
    }
    data.initVariances(heritability);
    bool noscale = true;
    return new BayesCXCI(data, data.varGenotypic, data.varResidual, pi, piPar[0], piPar[1], estimatePi, piNDC, piNDCpar, estimatePiNDC, 0, numKeptFemales, noscale);
}

Model* XCI::buildModelStageTwo(Data &data, const string &bayesType, const float heritability, const float pi, const VectorXf &piPar, const bool estimatePi, const float piNDC, const Vector2f &piNDCpar, const bool estimatePiNDC, const string &snpResFile, const float piGxE, const bool estimatePiGxE){
    if (!myMPI::rank) {
        cout << endl << "State Two: estimating genotype-by-sex effects given the estimated dosage compensation probabilities." << endl;
    }
    bool noscale = true;
    VectorXf snpPiNDC(data.numIncdSnps);
    readSnpPiNDC(snpPiNDC, data, snpResFile);
    Vector3f pis;
    pis << 1.0-pi, pi*(1.0-piGxE), pi*piGxE;  // prob of zero effect, sex-shared effect, sex-specific effect
    return new BayesCXCIgxs(data, data.varGenotypic, data.varResidual, pis, piPar, estimatePi, piNDC, piNDCpar, estimatePiNDC, snpPiNDC, estimatePiGxE, numKeptMales, numKeptFemales, noscale);
}


void XCI::simu(Data &data, const float pi, const float heritability, const float probNDC, const float probGxS, const bool removeQTL, const string &title, const int seed){
    if (myMPI::rank==0) cout << "Simulation start ..." << endl;
    unsigned numQTL = pi*data.numIncdSnps;
    vector<unsigned> indices(data.numIncdSnps);
    std::iota(indices.begin(), indices.end(), 0);
    std::random_shuffle(indices.begin(), indices.end());
    vector<SnpInfo*> QTLvec(numQTL);
    MatrixXf Q(data.numKeptInds, numQTL);
    VectorXf alphaMale(numQTL);
    VectorXf alphaFemale(numQTL);
    VectorXf isNDC(numQTL);
    SnpInfo *qtl;
    unsigned numFDC = 0;
    
    Stat::engine.seed(seed);
    
    for (unsigned j=0; j<numQTL; ++j) {
        qtl = data.incdSnpInfoVec[indices[j]];
        qtl->isQTL = true;
        QTLvec[j] = qtl;
        Q.col(j) = data.Z.col(qtl->index);
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
        alphaMale[j] = Stat::snorm();
        if (Stat::ranf() < probGxS) {
            alphaFemale[j] = Stat::snorm();
        } else {
            alphaFemale[j] = alphaMale[j];
        }
//        cout << j << " " << alphaMale[j] << " " << alphaFemale[j] << endl;
    }
    
//    VectorXf g = Q*alpha;
    VectorXf g(data.numKeptInds);
    g.head(numKeptMales) = Q.block(0, 0, numKeptMales, numQTL)*alphaMale;
    g.tail(numKeptFemales) = Q.block(numKeptMales, 0, numKeptFemales, numQTL)*alphaFemale;
    
//    cout << "alpha " << alphaMale << " " << alphaMale << endl;
//    cout << Q.block(0, 0, 5, numQTL) << endl;
//    cout << g.head(10) << endl;
    
    // calculate genetic variance with MPI
    // set heritability based on the genetic variance in females
    float my_sum, my_ss, sum, ss;
    unsigned n;
    my_sum = g.head(numKeptMales).sum();
    my_ss  = g.head(numKeptMales).squaredNorm();
    MPI_Allreduce(&my_sum, &sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&my_ss, &ss, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&numKeptMales, &n, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    float genVarMale = ss/float(n) - sum*sum/float(n*n);
    
    my_sum = g.tail(numKeptFemales).sum();
    my_ss  = g.tail(numKeptFemales).squaredNorm();
    MPI_Allreduce(&my_sum, &sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&my_ss, &ss, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&numKeptFemales, &n, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    float genVarFemale = ss/float(n) - sum*sum/float(n*n);

    my_sum = g.sum();
    my_ss  = g.squaredNorm();
    MPI_Allreduce(&my_sum, &sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&my_ss, &ss, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&data.numKeptInds, &n, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    float genVarAll = ss/float(n) - sum*sum/float(n*n);
    
    float resVar = genVarFemale*(1.0-heritability)/heritability;
    float resSD  = sqrt(resVar);

    for (unsigned i=0; i<data.numKeptInds; ++i) {
        data.y[i] = g[i] + Stat::snorm()*resSD;
    }
    float my_ypy = (data.y.array()-data.y.mean()).square().sum();
    
    MPI_Allreduce(&my_ypy, &data.ypy, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    
    my_sum = data.y.head(numKeptMales).sum();
    my_ss  = data.y.head(numKeptMales).squaredNorm();
    MPI_Allreduce(&my_sum, &sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&my_ss, &ss, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&numKeptMales, &n, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    float phenVarMale = ss/float(n) - sum*sum/float(n*n);
    
    my_sum = data.y.tail(numKeptFemales).sum();
    my_ss  = data.y.tail(numKeptFemales).squaredNorm();
    MPI_Allreduce(&my_sum, &sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&my_ss, &ss, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&numKeptFemales, &n, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    float phenVarFemale = ss/float(n) - sum*sum/float(n*n);
    
    my_sum = data.y.sum();
    my_ss  = data.y.squaredNorm();
    MPI_Allreduce(&my_sum, &sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&my_ss, &ss, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&data.numKeptInds, &n, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    float phenVarAll = ss/float(n) - sum*sum/float(n*n);
    
    data.varGenotypic = genVarFemale;
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
    if (!myMPI::rank) {
        ofstream out(filename.c_str());
        out << boost::format("%6s %20s %6s %12s %8s %12s %12s %6s\n")
        % "Id"
        % "Name"
        % "Chrom"
        % "Position"
        % "GeneFrq"
        % "EffectMale"
        % "EffectFemale"
        % "EscapeXCI";
        for (unsigned j=0; j<numQTL; ++j) {
            qtl = data.incdSnpInfoVec[indices[j]];
            out << boost::format("%6s %20s %6s %12s %8s %12s %12s %6s\n")
            % (j+1)
            % qtl->ID
            % qtl->chrom
            % qtl->physPos
            % qtl->af
            % alphaMale[j]
            % alphaFemale[j]
            % isNDC[j];
        }
        out.close();
    }
    
    string phenfilename = title + ".phen.thread" + to_string(static_cast<long long>(myMPI::rank));
    ofstream out2(phenfilename.c_str());
    for (unsigned i=0; i<data.numKeptInds; ++i) {
        IndInfo *ind = data.keptIndInfoVec[i];
        unsigned idx = ind->index;
        out2 << boost::format("%12s %12s %12.6f\n")
        % ind->famID
        % ind->indID
        % data.y[idx];
    }
    out2.close();
    
    if (!myMPI::rank) {
        cout << "\nSimulated " << numQTL << " QTL with " << probNDC*100 << "% escaped from XCI and " << probGxS*100 << "% with genotype-by-sex effects." << endl;
        cout << "Simulated genotypic variance: " << genVarAll << " (male: " << genVarMale << "; female: " << genVarFemale << ")" <<  endl;
        cout << "Simulated residual  variance: " << resVar << endl;
        cout << "Simulated heritability: " << genVarAll/phenVarAll << " (male: " << genVarMale/phenVarMale << "; female: " << genVarFemale/phenVarFemale << ")" <<  endl;
        if (removeQTL) cout << "QTL removed from the analysis." << endl;
        cout << "Saved simulated QTL info to [" << filename << "]." << endl;
    }
}

void XCI::outputResults(const Data &data, const vector<McmcSamples*> &mcmcSampleVec, const string &bayesType, const string &title){
    if (myMPI::rank) return;
    if (bayesType == "Xgxs") {
        McmcSamples *snpEffectsMale = NULL;
        McmcSamples *snpEffectsFemale = NULL;
        McmcSamples *deltaNDC = NULL;
        McmcSamples *deltaFDC = NULL;
        McmcSamples *deltaGxS = NULL;
        for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
            if (mcmcSampleVec[i]->label == "SnpEffectsMale") snpEffectsMale = mcmcSampleVec[i];
            if (mcmcSampleVec[i]->label == "SnpEffectsFemale") snpEffectsFemale = mcmcSampleVec[i];
            if (mcmcSampleVec[i]->label == "DeltaNDC") deltaNDC = mcmcSampleVec[i];
            if (mcmcSampleVec[i]->label == "DeltaFDC") deltaFDC = mcmcSampleVec[i];
            if (mcmcSampleVec[i]->label == "DeltaGxS") deltaGxS = mcmcSampleVec[i];
        }
        
        string filename = title + ".snpRes";
        ofstream out(filename.c_str());
        out << boost::format("%6s %20s %6s %12s %8s %12s %12s %12s %12s %8s %8s %8s %8s\n")
        % "Id"
        % "Name"
        % "Chrom"
        % "Position"
        % "GeneFrq"
        % "Effect_Male"
        % "SE_Male"
        % "Effect_Female"
        % "SE_Female"
        % "PIP"
        % "PrNDC"
        % "PrFDC"
        % "PrGxS";
        for (unsigned i=0, idx=0; i<data.numSnps; ++i) {
            SnpInfo *snp = data.snpInfoVec[i];
            if(!data.fullSnpFlag[i]) continue;
            //            if(snp->isQTL) continue;
            out << boost::format("%6s %20s %6s %12s %8.3f %12.6f %12.6f %12.6f %12.6f %8.3f %8.3f %8.3f %8.3f\n")
            % (idx+1)
            % snp->ID
            % snp->chrom
            % snp->physPos
            % snp->af
            % snpEffectsMale->posteriorMean[idx]
            % sqrt(snpEffectsMale->posteriorSqrMean[idx]-snpEffectsMale->posteriorMean[idx]*snpEffectsMale->posteriorMean[idx])
            % snpEffectsFemale->posteriorMean[idx]
            % sqrt(snpEffectsFemale->posteriorSqrMean[idx]-snpEffectsFemale->posteriorMean[idx]*snpEffectsFemale->posteriorMean[idx])
            % snpEffectsMale->pip[idx]
            % deltaNDC->posteriorMean[idx]
            % deltaFDC->posteriorMean[idx]
            % deltaGxS->posteriorMean[idx];
            ++idx;
        }
        out.close();
    }
    else if (bayesType == "Cgxs" || bayesType == "Ngxs") {
        McmcSamples *snpEffectsMale = NULL;
        McmcSamples *snpEffectsFemale = NULL;
        McmcSamples *deltaNDC = NULL;
        McmcSamples *deltaGxS = NULL;
        for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
            if (mcmcSampleVec[i]->label == "SnpEffectsMale") snpEffectsMale = mcmcSampleVec[i];
            if (mcmcSampleVec[i]->label == "SnpEffectsFemale") snpEffectsFemale = mcmcSampleVec[i];
            if (mcmcSampleVec[i]->label == "DeltaNDC") deltaNDC = mcmcSampleVec[i];
            if (mcmcSampleVec[i]->label == "DeltaGxS") deltaGxS = mcmcSampleVec[i];
        }

        string filename = title + ".snpRes";
        ofstream out(filename.c_str());
        out << boost::format("%6s %20s %6s %12s %8s %12s %12s %12s %12s %8s %8s %8s")
        % "Id"
        % "Name"
        % "Chrom"
        % "Position"
        % "GeneFrq"
        % "Effect_Male"
        % "SE_Male"
        % "Effect_Female"
        % "SE_Female"
        % "PIP"
        % "PrNDC"
        % "PrGxS";
        if (bayesType == "Ngxs") boost::format("%6s") % "Window";
        out << "\n";
        for (unsigned i=0, idx=0; i<data.numSnps; ++i) {
            SnpInfo *snp = data.snpInfoVec[i];
            if(!data.fullSnpFlag[i]) continue;
//            if(snp->isQTL) continue;
            out << boost::format("%6s %20s %6s %12s %8.3f %12.6f %12.6f %12.6f %12.6f %8.3f %8.3f %8.3f")
            % (idx+1)
            % snp->ID
            % snp->chrom
            % snp->physPos
            % snp->af
            % snpEffectsMale->posteriorMean[idx]
            % sqrt(snpEffectsMale->posteriorSqrMean[idx]-snpEffectsMale->posteriorMean[idx]*snpEffectsMale->posteriorMean[idx])
            % snpEffectsFemale->posteriorMean[idx]
            % sqrt(snpEffectsFemale->posteriorSqrMean[idx]-snpEffectsFemale->posteriorMean[idx]*snpEffectsFemale->posteriorMean[idx])
            % snpEffectsMale->pip[idx]
            % deltaNDC->posteriorMean[idx]
            % deltaGxS->posteriorMean[idx];
            if (bayesType == "Ngxs") boost::format("%6s") % snp->window;
            out << "\n";
            ++idx;
        }
        out.close();
    }
    else {
        McmcSamples *snpEffects = NULL;
        McmcSamples *deltaNDC = NULL;
        for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
            if (mcmcSampleVec[i]->label == "SnpEffects") snpEffects = mcmcSampleVec[i];
            if (mcmcSampleVec[i]->label == "DeltaNDC") deltaNDC = mcmcSampleVec[i];
        }
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
//            if(snp->isQTL) continue;
            out << boost::format("%6s %20s %6s %12s %8.3f %12.6f %12.6f %8.3f %8.3f\n")
            % (idx+1)
            % snp->ID
            % snp->chrom
            % snp->physPos
            % snp->af
            % snpEffects->posteriorMean[idx]
            % sqrt(snpEffects->posteriorSqrMean[idx]-snpEffects->posteriorMean[idx]*snpEffects->posteriorMean[idx])
            % snpEffects->pip[idx]
            % deltaNDC->posteriorMean[idx];
            ++idx;
        }
        out.close();
    }
}

void XCI::readSnpPiNDC(VectorXf &snpPiNDC, const Data &data, const string &snpResFile) {
    ifstream in(snpResFile.c_str());
    if (!in) throw ("Error: can not open the SNP result file [" + snpResFile + "] to read.");
    if (myMPI::rank==0)
    cout << "Reading estimated SNP NDC probabilities from [" + snpResFile + "]." << endl;

    Gadget::Tokenizer colData;
    string inputStr;
    string sep(" \t");
    unsigned line=0;
    unsigned idx = 999;
    
    while (getline(in,inputStr)) {
        colData.getTokens(inputStr, sep);
        if (!line) {
            for (unsigned i=0; i<colData.size(); ++i) {
                if (colData[i] == "PrNDC") {
                    idx = i;
                    break;
                }
            }
            if (idx == 999 && myMPI::rank==0) {
                cout << "Error: No column in [" + snpResFile + "] called PrNDC!" << endl;
            }
        } else {
            snpPiNDC[line-1] = atof(colData[idx].c_str());
        }
        ++line;
    }
    in.close();
    if (myMPI::rank==0)
    cout << "The mean of the estimated PrNDC across " << data.numIncdSnps << " SNPs: " << snpPiNDC.mean() << endl;
}


void BayesCXCI::FixedEffects::sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &X, const unsigned nmale, const unsigned nfemale,
                                           const VectorXf &XPXdiagMale, const VectorXf &XPXdiagFemale, const float varem, const float varef){
    float invVarem = 1.0f/varem;
    float invVaref = 1.0f/varef;
    float rhs;
    for (unsigned i=0; i<size; ++i) {
        if (!XPXdiagMale[i] && !XPXdiagFemale[i]) continue;
        float oldSample = values[i];
        float my_rhs = X.col(i).head(nmale).dot(ycorrm)*invVarem + X.col(i).tail(nfemale).dot(ycorrf)*invVaref;
        MPI_Allreduce(&my_rhs, &rhs, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        rhs += XPXdiagMale[i]*oldSample*invVarem + XPXdiagFemale[i]*oldSample*invVaref;
        float invLhs = 1.0f/(XPXdiagMale[i]*invVarem + XPXdiagFemale[i]*invVaref);
        float bhat = invLhs*rhs;
        values[i] = Normal::sample(bhat, invLhs);
        ycorrm += X.col(i).head(nmale) * (oldSample - values[i]);
        ycorrf += X.col(i).tail(nfemale) * (oldSample - values[i]);
    }
}


void BayesCXCI::SnpEffects::sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z, const VectorXf &ZPZdiag,
                                         const VectorXf &ZPZdiagMale, const VectorXf &ZPZdiagFemale,
                                         const VectorXf &ZPZdiagMaleRank, const VectorXf &ZPZdiagFemaleRank,
                                         const unsigned nmale, const unsigned nfemale, const float p,
                                         const float sigmaSq, const float pi, const float varem, const float varef,
                                         VectorXf &deltaNDC, VectorXf &ghatm, VectorXf &ghatf){
    // sample beta, delta, deltaNDC jointly
    // f(beta, delta, deltaNDC) propto f(beta | delta, deltaNDC) f(delta | deltaNDC) f(deltaNDC)
    
    //cout << "check pi " << pi << " p " << p << " sigmaSq " << sigmaSq << endl;
    
    sumSq = 0.0;
    numNonZeros = 0;
    
    ghatm.setZero(ycorrm.size());
    ghatf.setZero(ycorrf.size());
    
    Vector2f rhsFemale;     // 0: FDC, 1: NDC, corresponding to deltaNDC
    Vector2f my_rhs, rhs;   // 0: FDC, 1: NDC, corresponding to deltaNDC
    Vector2f invLhs;        // 0: FDC, 1: NDC, corresponding to deltaNDC
    Vector2f uhat;          // 0: FDC, 1: NDC, corresponding to deltaNDC
    Vector2f logDeltaNDC;   // 0: FDC, 1: NDC, corresponding to deltaNDC
    
    Vector2f dmcoef;
    dmcoef << 0.5, 1.0;

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
    float probDeltaNDC1;
    float oldSampleDeltaNDC, sampleDeltaNDC;
    
    for (unsigned i=0; i<size; ++i) {
        oldSample = values[i];
        oldSampleDeltaNDC = deltaNDC[i];

        rhsMale = Z.col(i).head(nmale).dot(ycorrm) + ZPZdiagMaleRank[i]*oldSample;
        rhsFemale[1] = Z.col(i).tail(nfemale).dot(ycorrf) + dmcoef[deltaNDC[i]]*ZPZdiagFemaleRank[i]*oldSample;
        rhsFemale[0] = 0.5f*rhsFemale[1];
        rhsMale *= invVarem;
        rhsFemale *= invVaref;

//        ycorrm += Z.col(i).head(nmale) * oldSample;
//        if (deltaNDC[i])
//            ycorrf += Z.col(i).tail(nfemale) * oldSample;
//        else
//            ycorrf += Z.col(i).tail(nfemale) * oldSample * 0.5f;
//        
//        rhsMale = Z.col(i).head(nmale).dot(ycorrm) * invVarem;
//        rhsFemale[1] = Z.col(i).tail(nfemale).dot(ycorrf) * invVaref;
//        rhsFemale[0] = rhsFemale[1] * 0.5f;
        
        my_rhs[1] = rhsMale + rhsFemale[1];
        my_rhs[0] = rhsMale + rhsFemale[0];
        
        MPI_Allreduce(&my_rhs[0], &rhs[0], 2, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        
        invLhs[1] = 1.0f/(ZPZdiagMale[i]*invVarem + ZPZdiagFemale[i]*invVaref       + invSigmaSq);
        invLhs[0] = 1.0f/(ZPZdiagMale[i]*invVarem + ZPZdiagFemale[i]*0.25f*invVaref + invSigmaSq);
        
        uhat.array() = invLhs.array()*rhs.array();
        
        //sample deltaNDC
        
        logDeltaNDC[1] = 0.5f*rhs[1]*uhat[1] + logf(sqrt(invLhs[1])*pi + expf(0.5f*(logSigmaSq-rhs[1]*uhat[1]))*(1.0f-pi)) + logP;
        logDeltaNDC[0] = 0.5f*rhs[0]*uhat[0] + logf(sqrt(invLhs[0])*pi + expf(0.5f*(logSigmaSq-rhs[0]*uhat[0]))*(1.0f-pi)) + logPcomp;
        probDeltaNDC1 = 1.0f/(1.0f + expf(logDeltaNDC[0] - logDeltaNDC[1]));
        
        sampleDeltaNDC = bernoulli.sample(probDeltaNDC1);
        deltaNDC[i] = sampleDeltaNDC;
        
        
        // sample delta
        
        logDelta1 = 0.5*(logf(invLhs[sampleDeltaNDC]) + uhat[sampleDeltaNDC]*rhs[sampleDeltaNDC]) + logPi;
        logDelta0 = 0.5*logSigmaSq + logPiComp;
        probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
        
        if (bernoulli.sample(probDelta1)) {
            
            // sample effect
            
            sample = normal.sample(uhat[sampleDeltaNDC], invLhs[sampleDeltaNDC]);
            values[i] = sample;
            sumSq += sample * sample;
            ++numNonZeros;
            
            ycorrm += Z.col(i).head(nmale)*(oldSample - sample);
            ycorrf += Z.col(i).tail(nfemale)*(dmcoef[oldSampleDeltaNDC]*oldSample - dmcoef[sampleDeltaNDC]*sample);

            ghatm += Z.col(i).head(nmale) * sample;
            ghatf += Z.col(i).tail(nfemale) * (dmcoef[sampleDeltaNDC]*sample);
            
//            ycorrm -= Z.col(i).head(nmale) * sample;
//            ghatm += Z.col(i).head(nmale) * sample;
//            
//            if (deltaNDC[i]) {
//                ycorrf -= Z.col(i).tail(nfemale) * sample;
//                ghatf += Z.col(i).tail(nfemale) * sample;
//            } else {
//                ycorrf -= Z.col(i).tail(nfemale) * sample * 0.5f;
//                ghatf += Z.col(i).tail(nfemale) * sample * 0.5f;
//            }
        }
        else {
            if (oldSample) {
                ycorrm += Z.col(i).head(nmale)*oldSample;
                ycorrf += Z.col(i).tail(nfemale)*(dmcoef[oldSampleDeltaNDC]*oldSample);
            }
            values[i] = 0.0;
        }
    }
}


void BayesCXCI::ProbNDC::sampleFromFC(const unsigned numSnps, const unsigned numNDC){
    //cout << numSnps << " " << numNDC << endl;
    float alphaTilde = numNDC + alpha;
    float betaTilde  = numSnps - numNDC + beta;
    value = Beta::sample(alphaTilde, betaTilde);
}

void BayesCXCI::Rounding::computeYcorr(const VectorXf &y, const MatrixXf &X, const MatrixXf &Z,
                                           const VectorXf &deltaNDC, const unsigned int nmale, const unsigned int nfemale,
                                           const VectorXf &fixedEffects, const VectorXf &snpEffects, VectorXf &ycorrm, VectorXf &ycorrf){
    if (count++ % 100) return;
    VectorXf oldYcorrm = ycorrm;
    VectorXf oldYcorrf = ycorrf;
    VectorXf ycorr = y - X*fixedEffects;
    for (unsigned i=0; i<snpEffects.size(); ++i) {
        if (snpEffects[i]) {
            if (deltaNDC[i]) {
                ycorr -= Z.col(i)*snpEffects[i];
            } else {
                ycorr.head(nmale) -= Z.col(i).head(nmale)*snpEffects[i];
                ycorr.tail(nfemale) -= Z.col(i).tail(nfemale)*snpEffects[i]*0.5f;
            }
        }
    }
    ycorrm = ycorr.head(nmale);
    ycorrf = ycorr.tail(nfemale);
    float my_ss = (ycorrm - oldYcorrm).squaredNorm() + (ycorrf - oldYcorrf).squaredNorm();
    float ss;
    MPI_Allreduce(&my_ss, &ss, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    value = sqrt(ss);
}

void BayesCXCI::sampleUnknowns(){
    
    fixedEffects.sampleFromFC(ycorrm, ycorrf, data.X, nmale, nfemale, XPXdiagMale, XPXdiagFemale, varem.value, varef.value);
    
//    unsigned cnt=0;
//    do {
        snpEffects.sampleFromFC(ycorrm, ycorrf, data.Z, data.ZPZdiag, ZPZdiagMale, ZPZdiagFemale, ZPZdiagMaleRank, ZPZdiagFemaleRank,
                                nmale, nfemale, piDeltaNDC.value, sigmaSq.value, pi.value, varem.value, varef.value, deltaNDC.values, ghatm, ghatf);
//        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
//    } while (snpEffects.numNonZeros == 0);
    
    if(estimatePiNDC) piDeltaNDC.sampleFromFC(snpEffects.size, deltaNDC.values.sum());
    piNDC.value = deltaNDC.values.sum()/(float)snpEffects.size;
//    piNDC.value = deltaNDC.values.sum()/snpEffects.size;
    sigmaSq.sampleFromFC(snpEffects.sumSq, snpEffects.numNonZeros);
//    sigmaSq.value = 1;
    if(estimatePi) pi.sampleFromFC(snpEffects.size, snpEffects.numNonZeros);
    //pi.value = snpEffects.numNonZeros/snpEffects.size;
    
    varem.sampleFromFC(ycorrm);
    varef.sampleFromFC(ycorrf);
    vargm.compute(ghatm);
    vargf.compute(ghatf);
    hsqm.compute(vargm.value, varem.value);
    hsqf.compute(vargf.value, varef.value);
    
//    vare.sampleFromFC(ycorr);
//    varg.compute(ghat);
//    hsq.compute(varg.value, vare.value);
    
    rounding.computeYcorr(data.y, data.X, data.Z, deltaNDC.values, nmale, nfemale, fixedEffects.values, snpEffects.values, ycorrm, ycorrf);
    nnzSnp.getValue(snpEffects.numNonZeros);
    
//    static unsigned iter = 0;
//    if (++iter < 5000) {
//        genVarPrior += (varg.value - genVarPrior)/iter;
//        piPrior += (pi.value - piPrior)/iter;
//        scale.compute(genVarPrior, piPrior, sigmaSq.scale);
//    }
//    scale.compute(varg.value, pi.value, sigmaSq.scale);
}

void BayesCXCI::getZPZdiag(const Data &data){
    // MPI
    XPXdiagMale.setZero(data.numFixedEffects);
    ZPZdiagMale.setZero(data.numIncdSnps);
    XPXdiagFemale.setZero(data.numFixedEffects);
    ZPZdiagFemale.setZero(data.numIncdSnps);
    VectorXf my_XPXdiagMale   = data.X.block(0, 0, nmale, data.numFixedEffects).colwise().squaredNorm();
    VectorXf my_XPXdiagFemale = data.X.block(nmale, 0, nfemale, data.numFixedEffects).colwise().squaredNorm();
    ZPZdiagMaleRank   = data.Z.block(0, 0, nmale, data.numIncdSnps).colwise().squaredNorm();
    ZPZdiagFemaleRank = data.Z.block(nmale, 0, nfemale, data.numIncdSnps).colwise().squaredNorm();
    MPI_Allreduce(&my_XPXdiagMale[0], &XPXdiagMale[0], data.numFixedEffects, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&my_XPXdiagFemale[0], &XPXdiagFemale[0], data.numFixedEffects, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&ZPZdiagMaleRank[0], &ZPZdiagMale[0], data.numIncdSnps, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&ZPZdiagFemaleRank[0], &ZPZdiagFemale[0], data.numIncdSnps, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
}


void BayesBXCI::SnpEffects::sampleFromFC(VectorXf &ycorr, const MatrixXf &Z, const VectorXf &ZPZdiag,
                                         const VectorXf &ZPZdiagMale, const VectorXf &ZPZdiagFemale,
                                         const unsigned int nmale, const unsigned int nfemale,
                                         const float p, const VectorXf &sigmaSq, const float pi,
                                         const float vare, VectorXf &deltaNDC, VectorXf &ghat) {
    numNonZeros = 0;
    
    ghat.setZero(ycorr.size());
    
    Vector2f rhsFemale;     // 0: FDC, 1: NDC, corresponding to deltaNDC
    Vector2f my_rhs, rhs;   // 0: FDC, 1: NDC, corresponding to deltaNDC
    Vector2f invLhs;        // 0: FDC, 1: NDC, corresponding to deltaNDC
    Vector2f uhat;          // 0: FDC, 1: NDC, corresponding to deltaNDC
    Vector2f logDeltaNDC;      // 0: FDC, 1: NDC, corresponding to deltaNDC
    
    float oldSample, sample;
    float logDelta0, logDelta1, probDelta1;
    float logPi = log(pi);
    float logPiComp = log(1.0-pi);
    float invVare = 1.0f/vare;
    float logP = log(p);
    float logPcomp = log(1.0f-p);
    float rhsMale;
    float probDeltaNDC1;
    float sampleDeltaNDC;
    float beta;
    
    for (unsigned i=0; i<size; ++i) {
        oldSample = values[i];
        ycorr.head(nmale) += Z.col(i).head(nmale) * oldSample;
        if (deltaNDC[i])
            ycorr.tail(nfemale) += Z.col(i).tail(nfemale) * oldSample;
        else
            ycorr.tail(nfemale) += Z.col(i).tail(nfemale) * oldSample * 0.5f;
        
        rhsMale = Z.col(i).head(nmale).dot(ycorr.head(nmale)) * invVare;
        rhsFemale[1] = Z.col(i).tail(nfemale).dot(ycorr.tail(nfemale)) * invVare;
        rhsFemale[0] = rhsFemale[1] * 0.5f;
        
        my_rhs[1] = rhsMale + rhsFemale[1];
        my_rhs[0] = rhsMale + rhsFemale[0];
        
        MPI_Allreduce(&my_rhs[0], &rhs[0], 2, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        
        invLhs[1] = 1.0f/((ZPZdiagMale[i] + ZPZdiagFemale[i]      )*invVare + 1.0f/sigmaSq[i]);
        invLhs[0] = 1.0f/((ZPZdiagMale[i] + ZPZdiagFemale[i]*0.25f)*invVare + 1.0f/sigmaSq[i]);
        
        uhat.array() = invLhs.array()*rhs.array();
        
        //sample deltaNDC
        
        logDeltaNDC[1] = 0.5f*rhs[1]*uhat[1] + logf(sqrt(invLhs[1])*pi + expf(0.5f*(logf(sigmaSq[i])-rhs[1]*uhat[1]))*(1.0f-pi)) + logP;
        logDeltaNDC[0] = 0.5f*rhs[0]*uhat[0] + logf(sqrt(invLhs[0])*pi + expf(0.5f*(logf(sigmaSq[i])-rhs[0]*uhat[0]))*(1.0f-pi)) + logPcomp;
        probDeltaNDC1 = 1.0f/(1.0f + expf(logDeltaNDC[0] - logDeltaNDC[1]));
        sampleDeltaNDC = bernoulli.sample(probDeltaNDC1);
        deltaNDC[i] = sampleDeltaNDC;
        
        
        // sample delta
        
        logDelta1 = 0.5*(logf(invLhs[sampleDeltaNDC]) + uhat[sampleDeltaNDC]*rhs[sampleDeltaNDC]) + logPi;
        logDelta0 = 0.5*logf(sigmaSq[i]) + logPiComp;
        probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
        
        if (bernoulli.sample(probDelta1)) {
            
            // sample effect
            
            sample = normal.sample(uhat[sampleDeltaNDC], invLhs[sampleDeltaNDC]);
            values[i] = sample;
            betaSq[i] = sample * sample;
            ++numNonZeros;
            
            ycorr.head(nmale) -= Z.col(i).head(nmale) * sample;
            ghat .head(nmale) += Z.col(i).head(nmale) * sample;
            
            if (deltaNDC[i]) {
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
    fixedEffects.sampleFromFC(ycorrm, ycorrf, data.X, nmale, nfemale, XPXdiagMale, XPXdiagFemale, varem.value, varef.value);
//    unsigned cnt=0;
//    do {
        snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, ZPZdiagMale, ZPZdiagFemale,
                                nmale, nfemale, piDeltaNDC.value, sigmaSq.values, pi.value, vare.value, deltaNDC.values, ghat);
//        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
//    } while (snpEffects.numNonZeros == 0);
    if(estimatePiNDC) piDeltaNDC.sampleFromFC(snpEffects.size, deltaNDC.values.sum());
    piNDC.value = deltaNDC.values.sum()/(float)snpEffects.size;
    sigmaSq.sampleFromFC(snpEffects.betaSq);
    if(estimatePi) pi.sampleFromFC(snpEffects.size, snpEffects.numNonZeros);
    vare.sampleFromFC(ycorr);
    
    varg.compute(ghat);
    hsq.compute(varg.value, vare.value);
    
    rounding.computeYcorr(data.y, data.X, data.Z, deltaNDC.values, nmale, nfemale, fixedEffects.values, snpEffects.values, ycorrm, ycorrf);
    nnzSnp.getValue(snpEffects.numNonZeros);
    
//    static unsigned iter = 0;
//    if (++iter < 5000) {
//        genVarPrior += (varg.value - genVarPrior)/iter;
//        piPrior += (pi.value - piPrior)/iter;
//        scale.compute(genVarPrior, piPrior, sigmaSq.scale);
//    }
    scale.compute(varg.value, pi.value, sigmaSq.scale);
}


void SBayesCXCI::SnpEffects::sampleFromFC(VectorXf &rcorrm, VectorXf &rcorrf, const MatrixXf &ZPZ, const MatrixXf &ZPZmale, const MatrixXf &ZPZfemale, const float piNDC, const float sigmaSq, const float pi, const float varem, const float varef, VectorXf &deltaNDC, VectorXf &ghatm, VectorXf &ghatf){
    sumSq = 0.0;
    numNonZeros = 0;
    
    Vector2f rhsFemale;     // 0: FDC, 1: NDC, corresponding to deltaNDC
    Vector2f rhs;           // 0: FDC, 1: NDC, corresponding to deltaNDC
    Vector2f invLhs;        // 0: FDC, 1: NDC, corresponding to deltaNDC
    Vector2f uhat;          // 0: FDC, 1: NDC, corresponding to deltaNDC
    Vector2f logDeltaNDC;   // 0: FDC, 1: NDC, corresponding to deltaNDC
    
    Vector2f dmcoef;
    dmcoef << 0.5, 1.0;
    
    float oldSample, sample;
    float logDelta0, logDelta1, probDelta1;
    float logPi = log(pi);
    float logPiComp = log(1.0-pi);
    float logSigmaSq = log(sigmaSq);
    float invVarem = 1.0f/varem;
    float invVaref = 1.0f/varef;
    float invSigmaSq = 1.0f/sigmaSq;
    float logPiNDC = log(piNDC);
    float logPiNDCcomp = log(1.0f-piNDC);
    float rhsMale;
    float probDeltaNDC1;
    float oldSampleDeltaNDC, sampleDeltaNDC;
    
    for (unsigned i=0; i<size; ++i) {
        oldSample = values[i];
        oldSampleDeltaNDC = deltaNDC[i];

        rhsMale = rcorrm[i] + ZPZmale(i,i)*oldSample;
        rhsFemale[1] = rcorrf[i] + dmcoef[deltaNDC[i]]*ZPZfemale(i,i)*oldSample;
        rhsFemale[0] = 0.5f*rhsFemale[1];
        rhsMale *= invVarem;
        rhsFemale *= invVaref;
        
        rhs[1] = rhsMale + rhsFemale[1];
        rhs[0] = rhsMale + rhsFemale[0];
        
        invLhs[1] = 1.0f/(ZPZmale(i,i)*invVarem + ZPZfemale(i,i)*invVaref       + invSigmaSq);
        invLhs[0] = 1.0f/(ZPZmale(i,i)*invVarem + ZPZfemale(i,i)*0.25f*invVaref + invSigmaSq);
        
        uhat.array() = invLhs.array()*rhs.array();
        
        //sample deltaNDC
        
        logDeltaNDC[1] = 0.5f*rhs[1]*uhat[1] + logf(sqrt(invLhs[1])*pi + expf(0.5f*(logSigmaSq-rhs[1]*uhat[1]))*(1.0f-pi)) + logPiNDC;
        logDeltaNDC[0] = 0.5f*rhs[0]*uhat[0] + logf(sqrt(invLhs[0])*pi + expf(0.5f*(logSigmaSq-rhs[0]*uhat[0]))*(1.0f-pi)) + logPiNDCcomp;
        probDeltaNDC1 = 1.0f/(1.0f + expf(logDeltaNDC[0] - logDeltaNDC[1]));
        
        sampleDeltaNDC = bernoulli.sample(probDeltaNDC1);
        deltaNDC[i] = sampleDeltaNDC;
        
        // sample delta
        
        logDelta1 = 0.5*(logf(invLhs[sampleDeltaNDC]) + uhat[sampleDeltaNDC]*rhs[sampleDeltaNDC]) + logPi;
        logDelta0 = 0.5*logSigmaSq + logPiComp;
        probDelta1 = 1.0f/(1.0f + expf(logDelta0-logDelta1));
        
        if (bernoulli.sample(probDelta1)) {
            
            // sample effect
            
            sample = normal.sample(uhat[sampleDeltaNDC], invLhs[sampleDeltaNDC]);
            values[i] = sample;
            sumSq += sample * sample;
            ++numNonZeros;
            
            rcorrm += ZPZmale.col(i)*(oldSample - sample);
            rcorrf += ZPZfemale.col(i).cwiseProduct(deltaNDC)*(dmcoef[oldSampleDeltaNDC]*oldSample - dmcoef[sampleDeltaNDC]*sample);
        }
        else {
            if (oldSample) {
                rcorrm += ZPZmale.col(i)*oldSample;
                rcorrf += ZPZfemale.col(i).cwiseProduct(deltaNDC)*(dmcoef[oldSampleDeltaNDC]*oldSample);
            }
            values[i] = 0.0;
        }
    }
}

void SBayesCXCI::sampleUnknowns() {

}


void BayesCXCIgxs::SnpEffects::sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z, const VectorXf &ZPZdiag,
                                            const VectorXf &ZPZdiagMale, const VectorXf &ZPZdiagFemale,
                                            const VectorXf &ZPZdiagMaleRank, const VectorXf &ZPZdiagFemaleRank,
                                            const unsigned nmale, const unsigned nfemale, const float piNDC,
                                            const float sigmaSq, const Vector3f &pis, const float varem, const float varef,
                                            VectorXf &deltaNDC, VectorXf &deltaGxS, VectorXf &ghatm, VectorXf &ghatf){
    // sample beta, delta, deltaNDC jointly
    // f(beta, delta, deltaNDC) propto f(beta | delta, deltaNDC) f(delta | deltaNDC) f(deltaNDC)
    
    //cout << "check pi " << pi << " p " << p << " sigmaSq " << sigmaSq << endl;
    
    sumSq = 0.0;
    numNonZeros = 0;
    numSnpMixComp.setZero(3);
    deltaGxS.setZero(size);
    
    ghatm.setZero(ycorrm.size());
    ghatf.setZero(ycorrf.size());
    
    Array2f my_rhs, rhs;   // 0: male, 1: female under NDC
    Array2f rhsFemale;     // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f rhsSame;       // 0: FDC, 1: NDC, corresponding to deltaNDC; same effect in males and females
    Array2f invLhsFemale;  // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f invLhsSame;    // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f uhatFemale;    // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f uhatSame;      // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f logDeltaNDC;      // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f quadGxE;       // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f quadSame;      // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f SEuhatFoverM;
    
    Array2f dmcoef;
    dmcoef << 0.5, 1.0;

    Array3f logPis = pis.array().log();  // 0: null effect, 1: same effect, 2: sex-specific effect
    Array3f logDelta;
    Array3f probDelta;
    
    float oldSampleMale, oldSampleFemale;
    float sampleMale, sampleFemale;
    float logSigmaSq = log(sigmaSq);
    float invVarem = 1.0f/varem;
    float invVaref = 1.0f/varef;
    float invSigmaSq = 1.0f/sigmaSq;
    float logPiNDC = log(piNDC);
    float logPiNDCcomp = log(1.0f-piNDC);
    float rhsMale;
    float invLhsMale;
    float uhatMale;
    float probDeltaNDC1;
    float oldSampleDeltaNDC, sampleDeltaNDC;
    
    unsigned delta;
    
    for (unsigned i=0; i<size; ++i) {
        oldSampleMale   = values(i,0);
        oldSampleFemale = values(i,1);
        oldSampleDeltaNDC = deltaNDC[i];

        my_rhs[0] = Z.col(i).head(nmale).dot(ycorrm) + ZPZdiagMaleRank[i]*oldSampleMale;
        my_rhs[1] = Z.col(i).tail(nfemale).dot(ycorrf) + dmcoef[deltaNDC[i]]*ZPZdiagFemaleRank[i]*oldSampleFemale;
        my_rhs[0] *= invVarem;
        my_rhs[1] *= invVaref;

//        ycorrm += Z.col(i).head(nmale) * oldSampleMale;
//        if (deltaNDC[i])
//            ycorrf += Z.col(i).tail(nfemale) * oldSampleFemale;
//        else
//            ycorrf += Z.col(i).tail(nfemale) * oldSampleFemale * 0.5f;
//        
//        my_rhs[0] = Z.col(i).head(nmale).dot(ycorrm) * invVarem;
//        my_rhs[1] = Z.col(i).tail(nfemale).dot(ycorrf) * invVaref;
        
//        cout << myMPI::rank << " snp " << i << " my_rhs " << my_rhs.transpose() << endl;
        
        MPI_Allreduce(&my_rhs[0], &rhs[0], 2, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        
        //cout << "snp " << i << " rhs " << rhs.transpose() << " " << deltaNDC[i] << " " << oldSampleMale << " " << oldSampleFemale << " ycorrm " << ycorrm.squaredNorm() << " ycorrf " << ycorrf.squaredNorm() << " dmcoef[deltaNDC[i]] " << dmcoef[deltaNDC[i]] << " ZPZdiagFemale[i] " << ZPZdiagFemale[i] << " invVarem " << invVarem << " invVaref " << invVaref << endl;
//        cout << myMPI::rank << " snp " << i << " rhs " << rhs.transpose() << " " << Z.col(i).head(nmale).squaredNorm() << " " << Z.col(i).tail(nfemale).squaredNorm() << " nmale " << nmale << " nfemale " << nfemale << " ZPZdiagMale[i] " << ZPZdiagMale[i] << " ZPZdiagFemale[i] " << ZPZdiagFemale[i] << endl;
        
        //float tmp;
        //cin >> tmp;
        
        rhsMale        = rhs[0];
        rhsFemale[0]   = rhs[1]*0.5f;
        rhsFemale[1]   = rhs[1];
        rhsSame        = rhsMale + rhsFemale;
        
        invLhsMale      = 1.0f/(ZPZdiagMale[i]*invVarem + invSigmaSq);
        invLhsFemale[0] = 1.0f/(ZPZdiagFemale[i]*0.25f*invVaref + invSigmaSq);
        invLhsFemale[1] = 1.0f/(ZPZdiagFemale[i]*invVaref       + invSigmaSq);
        invLhsSame[0]   = 1.0f/(ZPZdiagMale[i]*invVarem + ZPZdiagFemale[i]*0.25f*invVaref + invSigmaSq);
        invLhsSame[1]   = 1.0f/(ZPZdiagMale[i]*invVarem + ZPZdiagFemale[i]*invVaref       + invSigmaSq);
        
//        cout << ZPZdiagMale[i] << " " << ZPZdiagFemale[i] << " " << invLhsMale << " " << invLhsFemale.transpose() << " " << invVarem << " " << invVaref << endl;
        
        uhatMale   = invLhsMale  *rhsMale;
        uhatFemale = invLhsFemale*rhsFemale;
        uhatSame   = invLhsSame  *rhsSame;
        
        quadGxE  = uhatMale*rhsMale + uhatFemale*rhsFemale;
        quadSame = uhatSame*rhsSame;
        
//        if (uhatFemale[1]/uhatMale > 0.45 && uhatFemale[1]/uhatMale < 0.55) {
//            quadGxE.setZero(2);
//        }
//        else if (uhatFemale[0]/uhatMale > 1.95 && uhatFemale[0]/uhatMale < 2.05) {
//            quadGxE.setZero(2);
//        }
        
        SEuhatFoverM[1] = sqrt(invLhsFemale[1]/(uhatMale*uhatMale) + invLhsMale*uhatFemale[1]*uhatFemale[1]/(uhatMale*uhatMale*uhatMale*uhatMale));
        SEuhatFoverM[0] = sqrt(invLhsFemale[0]/(uhatMale*uhatMale) + invLhsMale*uhatFemale[0]*uhatFemale[0]/(uhatMale*uhatMale*uhatMale*uhatMale));
        if ((uhatFemale[1]/uhatMale - 0.001*SEuhatFoverM[1]) < 0.5 && (uhatFemale[1]/uhatMale + 0.001*SEuhatFoverM[1]) > 0.5) {
            quadGxE.setZero(2);
        } else if ((uhatFemale[0]/uhatMale - 0.001*SEuhatFoverM[0]) < 2.0 && (uhatFemale[0]/uhatMale + 0.001*SEuhatFoverM[0]) > 2.0) {
            quadGxE.setZero(2);
        }

//        cout << "uhatMale " << uhatMale << " uhatFemale " << uhatFemale.transpose() << " uhatSame " << uhatSame.transpose() << endl;
//        cout << "rhsSame " << rhsSame.transpose() << " invLhsSame " << invLhsSame.transpose() << " invVarem " << invVarem << " invVaref " << invVaref << endl;

//        float logLikeNDCSame = -0.5f*logSigmaSq + 0.5f*quadSame[1] + 0.5f*logf(invLhsSame[1]);
//        float logLikeFDCGxE  = -logSigmaSq + 0.5f*quadGxE[0] + 0.5f*logf(invLhsMale*invLhsFemale[0]);
//        
//        if (abs(logLikeNDCSame - logLikeFDCGxE) < 1e-3) {  // when 'NDC + Same effect' and 'FDC + GxE effect' models are not distinguishable, 'NDC + Same effect' model is preferred to avoid identifiability problem.
//            //sample deltaNDC
//            logDeltaNDC[1] = 0.5f*quadSame[1] + logf(sqrt(invLhsSame[1])*(pis[1]+pis[2]) + expf(0.5f*(logSigmaSq-quadSame[1]))*pis[0]) + logPiNDC;
//            logDeltaNDC[0] = 0.5f*quadSame[0] + logf(sqrt(invLhsSame[0])*(pis[1]+pis[2]) + expf(0.5f*(logSigmaSq-quadSame[0]))*pis[0]) + logPiNDCcomp;
//            
//            probDeltaNDC1 = 1.0f/(1.0f + expf(logDeltaNDC[0] - logDeltaNDC[1]));
//            sampleDeltaNDC = bernoulli.sample(probDeltaNDC1);
//            deltaNDC[i] = sampleDeltaNDC;
//
//            // sample delta
//            logDelta[0] = logPis[0];
//            logDelta[1] = 0.5f*(logf(invLhsSame[sampleDeltaNDC]) - logSigmaSq + quadSame[sampleDeltaNDC]) + logPis[1];
//            logDelta[2] = logPis[2];
//
//        }
//        else {
            //sample deltaNDC
        logDeltaNDC[1] = 0.5f*quadSame[1] + logf(sqrt(invLhsSame[1])*pis[1] + expf(0.5f*(quadGxE[1]-quadSame[1]))*sqrt(invLhsMale*invLhsFemale[1]*invSigmaSq)*pis[2] + expf(0.5f*(logSigmaSq-quadSame[1]))*pis[0]) + logPiNDC;
        logDeltaNDC[0] = 0.5f*quadSame[0] + logf(sqrt(invLhsSame[0])*pis[1] + expf(0.5f*(quadGxE[0]-quadSame[0]))*sqrt(invLhsMale*invLhsFemale[0]*invSigmaSq)*pis[2] + expf(0.5f*(logSigmaSq-quadSame[0]))*pis[0]) + logPiNDCcomp;
            
            probDeltaNDC1 = 1.0f/(1.0f + expf(logDeltaNDC[0] - logDeltaNDC[1]));
            sampleDeltaNDC = bernoulli.sample(probDeltaNDC1);
//        cout << quadSame.transpose() <<  " " << invLhsSame.transpose() << " " << quadGxE.transpose() << " " << invLhsMale << " " << invLhsFemale.transpose() << endl;
//        sampleDeltaNDC = 1;
            deltaNDC[i] = sampleDeltaNDC;
//        cout << "logDeltaNDC " << logDeltaNDC.transpose() << endl;
//        cout << sqrt(invLhsSame[1])*pis[1] << " " << (0.5f*(quadGxE[1]-quadSame[1])) << " " << sqrt(invLhsMale*invLhsFemale[1]*invSigmaSq)*pis[2] << " " << expf(0.5f*(logSigmaSq-quadSame[1]))*pis[0] << endl;
//        cout << "probDeltaNDC1 " << probDeltaNDC1 << " sampleDeltaNDC " << sampleDeltaNDC << endl;
//        cout << "quadGxE " << quadGxE.transpose() << " quadSame " << quadSame.transpose() << endl;
        
            // sample delta
            logDelta[0] = logPis[0];
            logDelta[1] = 0.5f*(logf(invLhsSame[sampleDeltaNDC]) - logSigmaSq + quadSame[sampleDeltaNDC]) + logPis[1];
            logDelta[2] = 0.5f*(logf(invLhsMale)+logf(invLhsFemale[sampleDeltaNDC]) - 2.0f*logSigmaSq + quadGxE[sampleDeltaNDC]) + logPis[2];
//        }

        for (unsigned j=0; j<3; ++j) {
            probDelta[j] = 1.0f/(logDelta-logDelta[j]).exp().sum();
        }

        delta = bernoulli.sample(probDelta);
//        cout << logDelta.transpose() << endl;
//        cout << "probDelta " << probDelta.transpose() << " delta " << delta << endl;
//        delta = 1;
        ++numSnpMixComp[delta];

        //        cout << logLikeNDCSame << " " << logLikeFDCGxE << " " << deltaNDC[i] << " " << delta << endl;

        if (delta) {
            // sample effect
            
            if (delta == 1) {  // same effect size in males and females
                values(i,0) = values(i,1) = sampleMale = sampleFemale = normal.sample(uhatSame[sampleDeltaNDC], invLhsSame[sampleDeltaNDC]);
//                cout << sampleMale << " " << sampleFemale << endl;
            } else {  // different effect sizes in males and females
                deltaNDC[i] = sampleDeltaNDC = bernoulli.sample(piNDC);
                values(i,0) = sampleMale = normal.sample(uhatMale, invLhsMale);
                values(i,1) = sampleFemale = normal.sample(uhatFemale[sampleDeltaNDC], invLhsFemale[sampleDeltaNDC]);
                deltaGxS[i] = 1;
            }
            
            sumSq += sampleMale * sampleMale + sampleFemale * sampleFemale;
            ++numNonZeros;

            ycorrm += Z.col(i).head(nmale)*(oldSampleMale - sampleMale);
            ycorrf += Z.col(i).tail(nfemale)*(dmcoef[oldSampleDeltaNDC]*oldSampleFemale - dmcoef[sampleDeltaNDC]*sampleFemale);
            
            ghatm += Z.col(i).head(nmale) * sampleMale;
            ghatf += Z.col(i).tail(nfemale) * (dmcoef[sampleDeltaNDC]*sampleFemale);

//            ycorrm -= Z.col(i).head(nmale) * sampleMale;
//            ghatm += Z.col(i).head(nmale) * sampleMale;
//            
//            if (deltaNDC[i]) {
//                ycorrf -= Z.col(i).tail(nfemale) * sampleFemale;
//                ghatf += Z.col(i).tail(nfemale) * sampleFemale;
//            } else {
//                ycorrf -= Z.col(i).tail(nfemale) * sampleFemale * 0.5f;
//                ghatf += Z.col(i).tail(nfemale) * sampleFemale * 0.5f;
//            }
        }
        else {
            if (oldSampleMale)   ycorrm += Z.col(i).head(nmale)*oldSampleMale;
            if (oldSampleFemale) ycorrf += Z.col(i).tail(nfemale)*(dmcoef[oldSampleDeltaNDC]*oldSampleFemale);
            values(i,0) = values(i,1) = 0.0;
        }
    }
}

void BayesCXCIgxs::SnpEffects::sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z, const VectorXf &ZPZdiag,
                                            const VectorXf &ZPZdiagMale, const VectorXf &ZPZdiagFemale,
                                            const VectorXf &ZPZdiagMaleRank, const VectorXf &ZPZdiagFemaleRank,
                                            const unsigned nmale, const unsigned nfemale,
                                            const VectorXf &snpPiNDC, const VectorXf &logSnpPiNDC, const VectorXf &logSnpPiNDCcomp,
                                            const float sigmaSq, const Vector3f &pis, const float varem, const float varef,
                                            VectorXf &deltaNDC, VectorXf &deltaGxS, VectorXf &ghatm, VectorXf &ghatf){
    // sample beta, delta, deltaNDC jointly
    // f(beta, delta, deltaNDC) propto f(beta | delta, deltaNDC) f(delta | deltaNDC) f(deltaNDC)
    
    //////////////////////////////////////////////
    ////// Given per-snp piNDC ///////////////////
    //////////////////////////////////////////////
    
    //cout << "check pi " << pi << " p " << p << " sigmaSq " << sigmaSq << endl;
    
    sumSq = 0.0;
    numNonZeros = 0;
    numSnpMixComp.setZero(3);
    deltaGxS.setZero(size);
    
    ghatm.setZero(ycorrm.size());
    ghatf.setZero(ycorrf.size());
    
    Array2f my_rhs, rhs;   // 0: male, 1: female under NDC
    Array2f rhsFemale;     // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f rhsSame;       // 0: FDC, 1: NDC, corresponding to deltaNDC; same effect in males and females
    Array2f invLhsFemale;  // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f invLhsSame;    // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f uhatFemale;    // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f uhatSame;      // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f logDeltaNDC;      // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f quadGxE;       // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f quadSame;      // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f SEuhatFoverM;
    
    Array2f dmcoef;
    dmcoef << 0.5, 1.0;
    
    Array3f logPis = pis.array().log();  // 0: null effect, 1: same effect, 2: sex-specific effect
    Array3f logDelta;
    Array3f probDelta;
    
    float oldSampleMale, oldSampleFemale;
    float sampleMale, sampleFemale;
    float logSigmaSq = log(sigmaSq);
    float invVarem = 1.0f/varem;
    float invVaref = 1.0f/varef;
    float invSigmaSq = 1.0f/sigmaSq;
    float rhsMale;
    float invLhsMale;
    float uhatMale;
    float probDeltaNDC1;
    float oldSampleDeltaNDC, sampleDeltaNDC;
    
    unsigned delta;
    
    for (unsigned i=0; i<size; ++i) {
        oldSampleMale   = values(i,0);
        oldSampleFemale = values(i,1);
        oldSampleDeltaNDC = deltaNDC[i];
        
        my_rhs[0] = Z.col(i).head(nmale).dot(ycorrm) + ZPZdiagMaleRank[i]*oldSampleMale;
        my_rhs[1] = Z.col(i).tail(nfemale).dot(ycorrf) + dmcoef[deltaNDC[i]]*ZPZdiagFemaleRank[i]*oldSampleFemale;
        my_rhs[0] *= invVarem;
        my_rhs[1] *= invVaref;
        
        MPI_Allreduce(&my_rhs[0], &rhs[0], 2, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        
        //cout << "snp " << i << " rhs " << rhs.transpose() << " " << deltaNDC[i] << " " << oldSampleMale << " " << oldSampleFemale << " ycorrm " << ycorrm.squaredNorm() << " ycorrf " << ycorrf.squaredNorm() << " dmcoef[deltaNDC[i]] " << dmcoef[deltaNDC[i]] << " ZPZdiagFemale[i] " << ZPZdiagFemale[i] << " invVarem " << invVarem << " invVaref " << invVaref << endl;
        //        cout << myMPI::rank << " snp " << i << " rhs " << rhs.transpose() << " " << Z.col(i).head(nmale).squaredNorm() << " " << Z.col(i).tail(nfemale).squaredNorm() << " nmale " << nmale << " nfemale " << nfemale << " ZPZdiagMale[i] " << ZPZdiagMale[i] << " ZPZdiagFemale[i] " << ZPZdiagFemale[i] << endl;
        
        //float tmp;
        //cin >> tmp;
        
        rhsMale        = rhs[0];
        rhsFemale[0]   = rhs[1]*0.5f;
        rhsFemale[1]   = rhs[1];
        rhsSame        = rhsMale + rhsFemale;
        
        invLhsMale      = 1.0f/(ZPZdiagMale[i]*invVarem + invSigmaSq);
        invLhsFemale[0] = 1.0f/(ZPZdiagFemale[i]*0.25f*invVaref + invSigmaSq);
        invLhsFemale[1] = 1.0f/(ZPZdiagFemale[i]*invVaref       + invSigmaSq);
        invLhsSame[0]   = 1.0f/(ZPZdiagMale[i]*invVarem + ZPZdiagFemale[i]*0.25f*invVaref + invSigmaSq);
        invLhsSame[1]   = 1.0f/(ZPZdiagMale[i]*invVarem + ZPZdiagFemale[i]*invVaref       + invSigmaSq);
        
        //        cout << ZPZdiagMale[i] << " " << ZPZdiagFemale[i] << " " << invLhsMale << " " << invLhsFemale.transpose() << " " << invVarem << " " << invVaref << endl;
        
        uhatMale   = invLhsMale  *rhsMale;
        uhatFemale = invLhsFemale*rhsFemale;
        uhatSame   = invLhsSame  *rhsSame;
        
        quadGxE  = uhatMale*rhsMale + uhatFemale*rhsFemale;
        quadSame = uhatSame*rhsSame;
        
        
//        SEuhatFoverM[1] = sqrt(invLhsFemale[1]/(uhatMale*uhatMale) + invLhsMale*uhatFemale[1]*uhatFemale[1]/(uhatMale*uhatMale*uhatMale*uhatMale));
//        SEuhatFoverM[0] = sqrt(invLhsFemale[0]/(uhatMale*uhatMale) + invLhsMale*uhatFemale[0]*uhatFemale[0]/(uhatMale*uhatMale*uhatMale*uhatMale));
//        if ((uhatFemale[1]/uhatMale - 0.001*SEuhatFoverM[1]) < 0.5 && (uhatFemale[1]/uhatMale + 0.001*SEuhatFoverM[1]) > 0.5) {
//            quadGxE.setZero(2);
//        } else if ((uhatFemale[0]/uhatMale - 0.001*SEuhatFoverM[0]) < 2.0 && (uhatFemale[0]/uhatMale + 0.001*SEuhatFoverM[0]) > 2.0) {
//            quadGxE.setZero(2);
//        }
        
        //sample deltaNDC
        logDeltaNDC[1] = 0.5f*quadSame[1] + logf(sqrt(invLhsSame[1])*pis[1] + expf(0.5f*(quadGxE[1]-quadSame[1]))*sqrt(invLhsMale*invLhsFemale[1]*invSigmaSq)*pis[2] + expf(0.5f*(logSigmaSq-quadSame[1]))*pis[0]) + logSnpPiNDC[i];
        logDeltaNDC[0] = 0.5f*quadSame[0] + logf(sqrt(invLhsSame[0])*pis[1] + expf(0.5f*(quadGxE[0]-quadSame[0]))*sqrt(invLhsMale*invLhsFemale[0]*invSigmaSq)*pis[2] + expf(0.5f*(logSigmaSq-quadSame[0]))*pis[0]) + logSnpPiNDCcomp[i];
        
        probDeltaNDC1 = 1.0f/(1.0f + expf(logDeltaNDC[0] - logDeltaNDC[1]));
        sampleDeltaNDC = bernoulli.sample(probDeltaNDC1);
        //        cout << quadSame.transpose() <<  " " << invLhsSame.transpose() << " " << quadGxE.transpose() << " " << invLhsMale << " " << invLhsFemale.transpose() << endl;
        //        sampleDeltaNDC = 1;
        deltaNDC[i] = sampleDeltaNDC;
        //        cout << "logDeltaNDC " << logDeltaNDC.transpose() << endl;
        //        cout << sqrt(invLhsSame[1])*pis[1] << " " << (0.5f*(quadGxE[1]-quadSame[1])) << " " << sqrt(invLhsMale*invLhsFemale[1]*invSigmaSq)*pis[2] << " " << expf(0.5f*(logSigmaSq-quadSame[1]))*pis[0] << endl;
        //        cout << "probDeltaNDC1 " << probDeltaNDC1 << " sampleDeltaNDC " << sampleDeltaNDC << endl;
        //        cout << "quadGxE " << quadGxE.transpose() << " quadSame " << quadSame.transpose() << endl;
        
        // sample delta
        logDelta[0] = logPis[0];
        logDelta[1] = 0.5f*(logf(invLhsSame[sampleDeltaNDC]) - logSigmaSq + quadSame[sampleDeltaNDC]) + logPis[1];
        logDelta[2] = 0.5f*(logf(invLhsMale)+logf(invLhsFemale[sampleDeltaNDC]) - 2.0f*logSigmaSq + quadGxE[sampleDeltaNDC]) + logPis[2];
        //        }
        
        for (unsigned j=0; j<3; ++j) {
            probDelta[j] = 1.0f/(logDelta-logDelta[j]).exp().sum();
        }
        
        delta = bernoulli.sample(probDelta);
        //        cout << logDelta.transpose() << endl;
        //        cout << "probDelta " << probDelta.transpose() << " delta " << delta << endl;
        //        delta = 1;
        ++numSnpMixComp[delta];
        
        //        cout << logLikeNDCSame << " " << logLikeFDCGxE << " " << deltaNDC[i] << " " << delta << endl;
        
        if (delta) {
            // sample effect
            
            if (delta == 1) {  // same effect size in males and females
                values(i,0) = values(i,1) = sampleMale = sampleFemale = normal.sample(uhatSame[sampleDeltaNDC], invLhsSame[sampleDeltaNDC]);
                //                cout << sampleMale << " " << sampleFemale << endl;
            } else {  // different effect sizes in males and females
                values(i,0) = sampleMale = normal.sample(uhatMale, invLhsMale);
                values(i,1) = sampleFemale = normal.sample(uhatFemale[sampleDeltaNDC], invLhsFemale[sampleDeltaNDC]);
                deltaGxS[i] = 1;
            }
            
            sumSq += sampleMale * sampleMale + sampleFemale * sampleFemale;
            ++numNonZeros;
            
            ycorrm += Z.col(i).head(nmale)*(oldSampleMale - sampleMale);
            ycorrf += Z.col(i).tail(nfemale)*(dmcoef[oldSampleDeltaNDC]*oldSampleFemale - dmcoef[sampleDeltaNDC]*sampleFemale);
            
            ghatm += Z.col(i).head(nmale) * sampleMale;
            ghatf += Z.col(i).tail(nfemale) * (dmcoef[sampleDeltaNDC]*sampleFemale);
            
        }
        else {
            if (oldSampleMale)   ycorrm += Z.col(i).head(nmale)*oldSampleMale;
            if (oldSampleFemale) ycorrf += Z.col(i).tail(nfemale)*(dmcoef[oldSampleDeltaNDC]*oldSampleFemale);
            values(i,0) = values(i,1) = 0.0;
        }
    }
}


void SBayesCXCIgxs::SnpEffects::sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z, const VectorXf &ZPZdiag,
                                            const VectorXf &ZPZdiagMale, const VectorXf &ZPZdiagFemale,
                                            const unsigned nmale, const unsigned nfemale, const float piNDC,
                                            const float sigmaSq, const Vector3f &pis, const float varem, const float varef,
                                            VectorXf &deltaNDC, VectorXf &deltaGxS, VectorXf &ghatm, VectorXf &ghatf){
    ////////////////////////////////////
    /////////// TO BE FINISH ///////////
    ////////////////////////////////////
    
    // sample beta, delta, deltaNDC jointly
    // f(beta, delta, deltaNDC) propto f(beta | delta, deltaNDC) f(delta | deltaNDC) f(deltaNDC)
    
    //cout << "check pi " << pi << " p " << p << " sigmaSq " << sigmaSq << endl;
    
    sumSq = 0.0;
    numNonZeros = 0;
    numSnpMixComp.setZero(3);
    deltaGxS.setZero(size);
    
    ghatm.setZero(ycorrm.size());
    ghatf.setZero(ycorrf.size());
    
    Array2f my_rhs, rhs;   // 0: male, 1: female under NDC
    Array2f rhsFemale;     // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f rhsSame;       // 0: FDC, 1: NDC, corresponding to deltaNDC; same effect in males and females
    Array2f invLhsFemale;  // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f invLhsSame;    // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f uhatFemale;    // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f uhatSame;      // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f logDeltaNDC;      // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f quadGxE;       // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f quadSame;      // 0: FDC, 1: NDC, corresponding to deltaNDC
    
    Array3f logPis = pis.array().log();  // 0: null effect, 1: same effect, 2: sex-specific effect
    Array3f logDelta;
    Array3f probDelta;
    
    float oldSampleMale, oldSampleFemale;
    float sampleMale, sampleFemale;
    float logSigmaSq = log(sigmaSq);
    float invVarem = 1.0f/varem;
    float invVaref = 1.0f/varef;
    float invSigmaSq = 1.0f/sigmaSq;
    float logPiNDC = log(piNDC);
    float logPiNDCcomp = log(1.0f-piNDC);
    float rhsMale;
    float invLhsMale;
    float uhatMale;
    float probDeltaNDC1;
    float sampleDeltaNDC;
    
    unsigned delta;
    
    for (unsigned i=0; i<size; ++i) {
        oldSampleMale   = values(i,0);
        oldSampleFemale = values(i,1);
        ycorrm += Z.col(i).head(nmale) * oldSampleMale;
        if (deltaNDC[i])
            ycorrf += Z.col(i).tail(nfemale) * oldSampleFemale;
        else
            ycorrf += Z.col(i).tail(nfemale) * oldSampleFemale * 0.5f;
        
        my_rhs[0] = Z.col(i).head(nmale).dot(ycorrm) * invVarem;
        my_rhs[1] = Z.col(i).tail(nfemale).dot(ycorrf) * invVaref;
        
        MPI_Allreduce(&my_rhs[0], &rhs[0], 2, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        
        rhsMale        = rhs[0];
        rhsFemale[0]   = rhs[1]*0.5f;
        rhsFemale[1]   = rhs[1];
        rhsSame        = rhsMale + rhsFemale;
        
        invLhsMale      = 1.0f/(ZPZdiagMale[i]*invVarem + invSigmaSq);
        invLhsFemale[0] = 1.0f/(ZPZdiagFemale[i]*0.25f*invVaref + invSigmaSq);
        invLhsFemale[1] = 1.0f/(ZPZdiagFemale[i]*invVaref       + invSigmaSq);
        invLhsSame[0]   = 1.0f/(ZPZdiagMale[i]*invVarem + ZPZdiagFemale[i]*0.25f*invVaref + invSigmaSq);
        invLhsSame[1]   = 1.0f/(ZPZdiagMale[i]*invVarem + ZPZdiagFemale[i]*invVaref       + invSigmaSq);
        
        //        cout << ZPZdiagMale[i] << " " << ZPZdiagFemale[i] << " " << invLhsMale << " " << invLhsFemale.transpose() << " " << invVarem << " " << invVaref << endl;
        
        uhatMale   = invLhsMale  *rhsMale;
        uhatFemale = invLhsFemale*rhsFemale;
        uhatSame   = invLhsSame  *rhsSame;
        
        quadGxE  = uhatMale*rhsMale + uhatFemale*rhsFemale;
        quadSame = uhatSame*rhsSame;
        
        if (uhatFemale[1]/uhatMale > 0.45 && uhatFemale[1]/uhatMale < 0.55) {
            quadGxE.setZero(2);
        }
        else if (uhatFemale[0]/uhatMale > 1.95 && uhatFemale[0]/uhatMale < 2.05) {
            quadGxE.setZero(2);
        }
        
        //        cout << "uhatMale " << uhatMale << " uhatFemale " << uhatFemale.transpose() << " uhatSame " << uhatSame.transpose() << endl;
        //        cout << "rhsSame " << rhsSame.transpose() << " invLhsSame " << invLhsSame.transpose() << " invVarem " << invVarem << " invVaref " << invVaref << endl;
        
        //        float logLikeNDCSame = -0.5f*logSigmaSq + 0.5f*quadSame[1] + 0.5f*logf(invLhsSame[1]);
        //        float logLikeFDCGxE  = -logSigmaSq + 0.5f*quadGxE[0] + 0.5f*logf(invLhsMale*invLhsFemale[0]);
        //
        //        if (abs(logLikeNDCSame - logLikeFDCGxE) < 1e-3) {  // when 'NDC + Same effect' and 'FDC + GxE effect' models are not distinguishable, 'NDC + Same effect' model is preferred to avoid identifiability problem.
        //            //sample deltaNDC
        //            logDeltaNDC[1] = 0.5f*quadSame[1] + logf(sqrt(invLhsSame[1])*(pis[1]+pis[2]) + expf(0.5f*(logSigmaSq-quadSame[1]))*pis[0]) + logPiNDC;
        //            logDeltaNDC[0] = 0.5f*quadSame[0] + logf(sqrt(invLhsSame[0])*(pis[1]+pis[2]) + expf(0.5f*(logSigmaSq-quadSame[0]))*pis[0]) + logPiNDCcomp;
        //
        //            probDeltaNDC1 = 1.0f/(1.0f + expf(logDeltaNDC[0] - logDeltaNDC[1]));
        //            sampleDeltaNDC = bernoulli.sample(probDeltaNDC1);
        //            deltaNDC[i] = sampleDeltaNDC;
        //
        //            // sample delta
        //            logDelta[0] = logPis[0];
        //            logDelta[1] = 0.5f*(logf(invLhsSame[sampleDeltaNDC]) - logSigmaSq + quadSame[sampleDeltaNDC]) + logPis[1];
        //            logDelta[2] = logPis[2];
        //
        //        }
        //        else {
        //sample deltaNDC
        logDeltaNDC[1] = 0.5f*quadSame[1] + logf(sqrt(invLhsSame[1])*pis[1] + expf(0.5f*(quadGxE[1]-quadSame[1]))*sqrt(invLhsMale*invLhsFemale[1]*invSigmaSq)*pis[2] + expf(0.5f*(logSigmaSq-quadSame[1]))*pis[0]) + logPiNDC;
        logDeltaNDC[0] = 0.5f*quadSame[0] + logf(sqrt(invLhsSame[0])*pis[1] + expf(0.5f*(quadGxE[0]-quadSame[0]))*sqrt(invLhsMale*invLhsFemale[0]*invSigmaSq)*pis[2] + expf(0.5f*(logSigmaSq-quadSame[0]))*pis[0]) + logPiNDCcomp;
        
        probDeltaNDC1 = 1.0f/(1.0f + expf(logDeltaNDC[0] - logDeltaNDC[1]));
        sampleDeltaNDC = bernoulli.sample(probDeltaNDC1);
        //        cout << quadSame.transpose() <<  " " << invLhsSame.transpose() << " " << quadGxE.transpose() << " " << invLhsMale << " " << invLhsFemale.transpose() << endl;
        //        sampleDeltaNDC = 1;
        deltaNDC[i] = sampleDeltaNDC;
        //        cout << "logDeltaNDC " << logDeltaNDC.transpose() << endl;
        //        cout << sqrt(invLhsSame[1])*pis[1] << " " << (0.5f*(quadGxE[1]-quadSame[1])) << " " << sqrt(invLhsMale*invLhsFemale[1]*invSigmaSq)*pis[2] << " " << expf(0.5f*(logSigmaSq-quadSame[1]))*pis[0] << endl;
        //        cout << "probDeltaNDC1 " << probDeltaNDC1 << " sampleDeltaNDC " << sampleDeltaNDC << endl;
        //        cout << "quadGxE " << quadGxE.transpose() << " quadSame " << quadSame.transpose() << endl;
        
        // sample delta
        logDelta[0] = logPis[0];
        logDelta[1] = 0.5f*(logf(invLhsSame[sampleDeltaNDC]) - logSigmaSq + quadSame[sampleDeltaNDC]) + logPis[1];
        logDelta[2] = 0.5f*(logf(invLhsMale)+logf(invLhsFemale[sampleDeltaNDC]) - 2.0f*logSigmaSq + quadGxE[sampleDeltaNDC]) + logPis[2];
        //        }
        
        for (unsigned j=0; j<3; ++j) {
            probDelta[j] = 1.0f/(logDelta-logDelta[j]).exp().sum();
        }
        
        delta = bernoulli.sample(probDelta);
        //        cout << logDelta.transpose() << endl;
        //        cout << "probDelta " << probDelta.transpose() << " delta " << delta << endl;
        //        delta = 1;
        ++numSnpMixComp[delta];
        
        //        cout << logLikeNDCSame << " " << logLikeFDCGxE << " " << deltaNDC[i] << " " << delta << endl;
        
        if (delta) {
            // sample effect
            
            if (delta == 1) {  // same effect size in males and females
                values(i,0) = values(i,1) = sampleMale = sampleFemale = normal.sample(uhatSame[sampleDeltaNDC], invLhsSame[sampleDeltaNDC]);
                //                cout << sampleMale << " " << sampleFemale << endl;
            } else {  // different effect sizes in males and females
                values(i,0) = sampleMale = normal.sample(uhatMale, invLhsMale);
                values(i,1) = sampleFemale = normal.sample(uhatFemale[sampleDeltaNDC], invLhsFemale[sampleDeltaNDC]);
                deltaGxS[i] = 1;
            }
            
            sumSq += sampleMale * sampleMale + sampleFemale * sampleFemale;
            ++numNonZeros;
            
            ycorrm -= Z.col(i).head(nmale) * sampleMale;
            ghatm += Z.col(i).head(nmale) * sampleMale;
            
            if (deltaNDC[i]) {
                ycorrf -= Z.col(i).tail(nfemale) * sampleFemale;
                ghatf += Z.col(i).tail(nfemale) * sampleFemale;
            } else {
                ycorrf -= Z.col(i).tail(nfemale) * sampleFemale * 0.5f;
                ghatf += Z.col(i).tail(nfemale) * sampleFemale * 0.5f;
            }
        }
        else {
            values(i,0) = values(i,1) = 0.0;
        }
    }
}


void BayesCXCIgxs::Rounding::computeYcorr(const VectorXf &y, const MatrixXf &X, const MatrixXf &Z,
                                       const VectorXf &deltaNDC, const unsigned int nmale, const unsigned int nfemale,
                                       const VectorXf &fixedEffects, const MatrixXf &snpEffects, VectorXf &ycorrm, VectorXf &ycorrf){
    if (count++ % 100) return;
    VectorXf oldYcorrm = ycorrm;
    VectorXf oldYcorrf = ycorrf;
    VectorXf ycorr = y - X*fixedEffects;
    for (unsigned i=0; i<snpEffects.rows(); ++i) {
        if (snpEffects(i,0)) {
            ycorr.head(nmale) -= Z.col(i).head(nmale)*snpEffects(i,0);
            if (deltaNDC[i]) {
                ycorr.tail(nfemale) -= Z.col(i).tail(nfemale)*snpEffects(i,1);
            } else {
                ycorr.tail(nfemale) -= Z.col(i).tail(nfemale)*snpEffects(i,1)*0.5f;
            }
        }
    }
    ycorrm = ycorr.head(nmale);
    ycorrf = ycorr.tail(nfemale);
    float my_ss = (ycorrm - oldYcorrm).squaredNorm() + (ycorrf - oldYcorrf).squaredNorm();
    float ss;
    MPI_Allreduce(&my_ss, &ss, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    value = sqrt(ss);
}

void BayesCXCIgxs::ProbMixComps::getValues(VectorXf &pis){
    values = pis;
    for (unsigned i=0; i<ndist; ++i) {
        (*this)[i]->value=values[i];
    }
}

void BayesCXCIgxs::sampleUnknowns(){
    fixedEffects.sampleFromFC(ycorrm, ycorrf, data.X, nmale, nfemale, XPXdiagMale, XPXdiagFemale, varem.value, varef.value);

//    unsigned cnt=0;
//    do {
    if (snpPiNDCgiven) { // two-stage model
        snpEffects.sampleFromFC(ycorrm, ycorrf, data.Z, data.ZPZdiag, ZPZdiagMale, ZPZdiagFemale, ZPZdiagMaleRank, ZPZdiagFemaleRank,
                                nmale, nfemale, snpPiNDC, logSnpPiNDC, logSnpPiNDCcomp,
                                sigmaSq.value, pis.values, varem.value, varef.value,
                                deltaNDC.values, deltaGxS.values, ghatm, ghatf);
    } else {
        snpEffects.sampleFromFC(ycorrm, ycorrf, data.Z, data.ZPZdiag, ZPZdiagMale, ZPZdiagFemale, ZPZdiagMaleRank, ZPZdiagFemaleRank,
                                nmale, nfemale, piDeltaNDC.value, sigmaSq.value, pis.values, varem.value, varef.value,
                                deltaNDC.values, deltaGxS.values, ghatm, ghatf);
    }
//        if (++cnt == 100) throw("Error: Zero SNP effect in the model for 100 cycles of sampling");
//    } while (snpEffects.numNonZeros == 0);
    
    if (estimatePiNDC) piDeltaNDC.sampleFromFC(snpEffects.size, deltaNDC.values.sum());
//    piDeltaNDC.sampleFromFC(snpEffects.numNonZeros, deltaNDC.values.sum());
    piNDC.value = deltaNDC.values.sum()/(float)snpEffects.size;
    sigmaSq.sampleFromFC(snpEffects.sumSq, 2.0*snpEffects.numNonZeros);  // both male and female effects contribute to sigmaSq
    if (estimatePi) {
        if (estimatePiGxS) {
            pis.sampleFromFC(snpEffects.numSnpMixComp);
        } else {
            pi.sampleFromFC(snpEffects.size, snpEffects.numNonZeros);
            VectorXf vec(3);
            vec << 1.0-pi.value, pi.value*(1.0-piGxSgiven), pi.value*piGxSgiven;
            pis.getValues(vec);
        }
    }
    pi.value = snpEffects.numNonZeros/(float)snpEffects.size;
    piGxS.value = snpEffects.numNonZeros ? snpEffects.numSnpMixComp[2]/(float)snpEffects.numNonZeros : 0.0;
    
    varem.sampleFromFC(ycorrm);
    varef.sampleFromFC(ycorrf);
    vargm.compute(ghatm);
    vargf.compute(ghatf);
    hsqm.compute(vargm.value, varem.value);
    hsqf.compute(vargf.value, varef.value);
    
    rounding.computeYcorr(data.y, data.X, data.Z, deltaNDC.values, nmale, nfemale, fixedEffects.values, snpEffects.values, ycorrm, ycorrf);
    nnzSnp.getValue(snpEffects.numNonZeros);
    
    snpEffectsMale.values = snpEffects.values.col(0);
    snpEffectsFemale.values = snpEffects.values.col(1);
}


void BayesNXCIgxs::SnpEffects::sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z, const vector<MatrixXf> &ZPZdiagMale,
                                            const vector<MatrixXf> &ZPZdiagFemale, const unsigned nmale, const unsigned nfemale,
                                            const float piNDC, const float sigmaSq, const Vector3f &pis, const float varem, const float varef,
                                            VectorXf &deltaNDC, VectorXf &deltaGxS, VectorXf &ghatm, VectorXf &ghatf){
    sumSq = 0.0;
    numNonZeros = 0;
    numNonZeroWind = 0;
    numSnpMixComp.setZero(3);
    deltaGxS.setZero(size);
    windDeltaGxS.setZero(numWindows);
    
    ghatm.setZero(ycorrm.size());
    ghatf.setZero(ycorrf.size());
    
    Array2f logDeltaNDC;      // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f quadGxE;       // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f quadSame;      // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f detInvLhsFemale;
    Array2f detInvLhsSame;
    
    Array3f logPis = pis.array().log();  // 0: null effect, 1: same effect, 2: sex-specific effect
    Array3f logDelta;
    Array3f probDelta;
    
    VectorXf oldSampleMale, oldSampleFemale;
    VectorXf sampleMale, sampleFemale;
    VectorXf my_rhsMale, rhsMale;
    VectorXf my_rhsFemaleNDC, rhsFemaleNDC;
    VectorXf uhatMale;
    
    MatrixXf rhsFemale;    // col(0): FDC; col(1): NDC
    MatrixXf rhsSame;      // col(0): FDC; col(1): NDC; same effect in males and females
    MatrixXf uhatFemale;   // col(0): FDC; col(1): NDC
    MatrixXf uhatSame;     // col(0): FDC; col(1): NDC

    MatrixXf invLhsMale;
    
    vector<MatrixXf> invLhsFemale;   // 0: FDC, 1: NDC
    vector<MatrixXf> invLhsSame;     // 0: FDC, 1: NDC
    
    invLhsFemale.resize(2);
    invLhsSame.resize(2);
    
    LLT<MatrixXf> llt;

    float logSigmaSq = log(sigmaSq);
    float invVarem = 1.0f/varem;
    float invVaref = 1.0f/varef;
    float invSigmaSq = 1.0f/sigmaSq;
    float logPiNDC = log(piNDC);
    float logPiNDCcomp = log(1.0f-piNDC);
    float probDeltaNDC1;
    float sampleDeltaNDC;
    float detInvLhsMale;
    
    VectorXf invSigmaSqVec;
    
    unsigned delta;
    unsigned starti, sizei;
    
    for (unsigned i=0; i<numWindows; ++i) {
        starti = windStart[i];
        sizei  = windSize[i];
        
//        cout << i << " " << starti << " " << sizei << " " << size << " " << values.cols() << " " << values.rows() << endl;
//        cout << values.block(starti,0,sizei,1) << endl;
        
        oldSampleMale   = values.block(starti,0,sizei,1);
        oldSampleFemale = values.block(starti,1,sizei,1);
        ycorrm += Z.block(0,starti,nmale,sizei) * oldSampleMale;
        if (deltaNDC[i])
            ycorrf += Z.block(nmale,starti,nfemale,sizei) * oldSampleFemale;
        else
            ycorrf += Z.block(nmale,starti,nfemale,sizei) * oldSampleFemale * 0.5f;
        
        my_rhsMale = Z.block(0,starti,nmale,sizei).transpose() * ycorrm * invVarem;
        my_rhsFemaleNDC = Z.block(nmale,starti,nfemale,sizei).transpose() * ycorrf * invVaref;
        
        rhsMale.resize(sizei);
        rhsFemaleNDC.resize(sizei);
        
        MPI_Allreduce(&my_rhsMale[0], &rhsMale[0], sizei, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&my_rhsFemaleNDC[0], &rhsFemaleNDC[0], sizei, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        
        rhsFemale.resize(sizei,2);
        rhsFemale.col(0) = rhsFemaleNDC*0.5f;
        rhsFemale.col(1) = rhsFemaleNDC;
        
        rhsSame.resize(sizei,2);
        rhsSame.col(0) = rhsMale + rhsFemale.col(0);
        rhsSame.col(1) = rhsMale + rhsFemale.col(1);
        
        invSigmaSqVec.setConstant(sizei, invSigmaSq);
        
        invLhsMale  = ZPZdiagMale[i];
        invLhsMale *= invVarem;
        invLhsFemale[0]  = ZPZdiagFemale[i];
        invLhsFemale[0] *= 0.25f*invVaref;
        invLhsFemale[1]  = ZPZdiagFemale[i];
        invLhsFemale[1] *= invVaref;
        invLhsSame[0] = invLhsMale + invLhsFemale[0];
        invLhsSame[1] = invLhsMale + invLhsFemale[1];
        
        invLhsMale += invSigmaSqVec.asDiagonal();
        invLhsMale  = invLhsMale.inverse();
        
        invLhsFemale[0] += invSigmaSqVec.asDiagonal();
        invLhsFemale[1] += invSigmaSqVec.asDiagonal();
        invLhsFemale[0]  = invLhsFemale[0].inverse();
        invLhsFemale[1]  = invLhsFemale[1].inverse();

        invLhsSame[0] += invSigmaSqVec.asDiagonal();
        invLhsSame[1] += invSigmaSqVec.asDiagonal();
        invLhsSame[0]  = invLhsSame[0].inverse();
        invLhsSame[1]  = invLhsSame[1].inverse();
        
        uhatMale   = invLhsMale*rhsMale;
        
        uhatFemale.resize(sizei,2);
        uhatSame.resize(sizei,2);

        uhatFemale.col(0) = invLhsFemale[0]*rhsFemale.col(0);
        uhatFemale.col(1) = invLhsFemale[1]*rhsFemale.col(1);
        
        uhatSame.col(0)   = invLhsSame[0]*rhsSame.col(0);
        uhatSame.col(1)   = invLhsSame[1]*rhsSame.col(1);
        
        quadGxE[0] = uhatMale.dot(rhsMale) + uhatFemale.col(0).dot(rhsFemale.col(0));
        quadGxE[1] = uhatMale.dot(rhsMale) + uhatFemale.col(1).dot(rhsFemale.col(1));
        quadSame[0] = uhatSame.col(0).dot(rhsSame.col(0));
        quadSame[1] = uhatSame.col(1).dot(rhsSame.col(1));
        
        detInvLhsMale = invLhsMale.determinant();
        detInvLhsFemale[0] = invLhsFemale[0].determinant();
        detInvLhsFemale[1] = invLhsFemale[1].determinant();
        detInvLhsSame[0] = invLhsSame[0].determinant();
        detInvLhsSame[1] = invLhsSame[1].determinant();
        
        // sample deltaNDC
        logDeltaNDC[1] = 0.5f*quadSame[1] + logf(sqrt(detInvLhsSame[1])*pis[1] + expf(0.5f*(quadGxE[1]-quadSame[1]))*sqrt(detInvLhsMale*detInvLhsFemale[1]*invSigmaSq/sizei)*pis[2] + expf(0.5f*(sizei*logSigmaSq-quadSame[1]))*pis[0]) + logPiNDC;
        logDeltaNDC[0] = 0.5f*quadSame[0] + logf(sqrt(detInvLhsSame[0])*pis[1] + expf(0.5f*(quadGxE[0]-quadSame[0]))*sqrt(detInvLhsMale*detInvLhsFemale[0]*invSigmaSq/sizei)*pis[2] + expf(0.5f*(sizei*logSigmaSq-quadSame[0]))*pis[0]) + logPiNDCcomp;
        
        probDeltaNDC1 = 1.0f/(1.0f + expf(logDeltaNDC[0] - logDeltaNDC[1]));
        sampleDeltaNDC = bernoulli.sample(probDeltaNDC1);
        windDeltaNDC[i] = sampleDeltaNDC;
        deltaNDC.segment(starti,sizei) = VectorXf::Constant(sizei,sampleDeltaNDC);
        
        // sample delta
        logDelta[0] = logPis[0];
        logDelta[1] = 0.5f*(logf(detInvLhsSame[sampleDeltaNDC]) - sizei*logSigmaSq + quadSame[sampleDeltaNDC]) + logPis[1];
        logDelta[2] = 0.5f*(logf(detInvLhsMale)+logf(detInvLhsFemale[sampleDeltaNDC]) - 2.0f*sizei*logSigmaSq + quadGxE[sampleDeltaNDC]) + logPis[2];
        
        for (unsigned j=0; j<3; ++j) {
            probDelta[j] = 1.0f/(logDelta-logDelta[j]).exp().sum();
        }
        
        delta = bernoulli.sample(probDelta);
        ++numSnpMixComp[delta];
        
        if (delta) {
            // sample effect
            if (delta == 1) {  // same effect size in males and females
                sampleMale.setZero(sizei);
                for (unsigned j=0;j<sizei;j++) sampleMale[j] = Stat::snorm();
                llt.compute(invLhsSame[sampleDeltaNDC]); // cholesky decomposition
                sampleMale = uhatSame.col(sampleDeltaNDC) + llt.matrixL()*sampleMale;
                values.block(starti,0,sizei,1) = values.block(starti,1,sizei,1) = sampleFemale = sampleMale;
            }
            else {  // different effect sizes in males and females
                sampleMale.setZero(sizei);
                for (unsigned j=0;j<sizei;j++) sampleMale[j] = Stat::snorm();
                llt.compute(invLhsMale); // cholesky decomposition
                sampleMale = uhatMale + llt.matrixL()*sampleMale;
                values.block(starti,0,sizei,1) = sampleMale;
                
                sampleFemale.setZero(sizei);
                for (unsigned j=0;j<sizei;j++) sampleFemale[j] = Stat::snorm();
                llt.compute(invLhsFemale[sampleDeltaNDC]); // cholesky decomposition
                sampleFemale = uhatFemale.col(sampleDeltaNDC) + llt.matrixL()*sampleFemale;
                values.block(starti,1,sizei,1) = sampleFemale;
                
                windDeltaGxS[i] = 1;
                deltaGxS.segment(starti,sizei) = VectorXf::Ones(sizei);
            }
            
            sumSq += sampleMale.dot(sampleMale) + sampleFemale.dot(sampleFemale);
            windDeltaEff[i] = 1;
            ++numNonZeroWind;
            numNonZeros += sizei;
            
            ycorrm -= Z.block(0,starti,nmale,sizei) * sampleMale;
            ghatm  += Z.block(0,starti,nmale,sizei) * sampleMale;
            
            if (deltaNDC[i]) {
                ycorrf -= Z.block(nmale,starti,nfemale,sizei) * sampleFemale;
                ghatf  += Z.block(nmale,starti,nfemale,sizei) * sampleFemale;
            } else {
                ycorrf -= Z.block(nmale,starti,nfemale,sizei) * sampleFemale * 0.5f;
                ghatf  += Z.block(nmale,starti,nfemale,sizei) * sampleFemale * 0.5f;
            }
        }
        else {
            values.block(starti,0,sizei,1) = values.block(starti,1,sizei,1) = VectorXf::Zero(sizei);
        }
    }
}

void BayesNXCIgxs::sampleUnknowns(){

    fixedEffects.sampleFromFC(ycorrm, ycorrf, data.X, nmale, nfemale, XPXdiagMale, XPXdiagFemale, varem.value, varef.value);
    
    snpEffects.sampleFromFC(ycorrm, ycorrf, data.Z, ZPZblockDiagMale, ZPZblockDiagFemale,
                            nmale, nfemale, piDeltaNDC.value, sigmaSq.value, pis.values, varem.value, varef.value,
                            deltaNDC.values, deltaGxS.values, ghatm, ghatf);
    
    if(estimatePiNDC) piDeltaNDC.sampleFromFC(snpEffects.numWindows, snpEffects.windDeltaNDC.sum());
    piNDC.value = snpEffects.windDeltaNDC.sum()/(float)snpEffects.numWindows;
    sigmaSq.sampleFromFC(snpEffects.sumSq, 2.0*snpEffects.numNonZeros);  // both male and female effects contribute to sigmaSq
    if (estimatePi) {
        if (estimatePiGxS) {
            pis.sampleFromFC(snpEffects.numSnpMixComp);
        } else {
            pi.sampleFromFC(snpEffects.size, snpEffects.numNonZeros);
            VectorXf vec;
            vec << 1.0-pi.value, pi.value*(1.0-piGxSgiven), pi.value*piGxSgiven;
            pis.getValues(vec);
        }
    }
    pi.value = snpEffects.numNonZeroWind/(float)snpEffects.numWindows;
    piGxS.value = snpEffects.numNonZeroWind ? snpEffects.numSnpMixComp[2]/(float)snpEffects.numNonZeroWind : 0.0;
    
    varem.sampleFromFC(ycorrm);
    varef.sampleFromFC(ycorrf);
    vargm.compute(ghatm);
    vargf.compute(ghatf);
    hsqm.compute(vargm.value, varem.value);
    hsqf.compute(vargf.value, varef.value);
    
    rounding.computeYcorr(data.y, data.X, data.Z, deltaNDC.values, nmale, nfemale, fixedEffects.values, snpEffects.values, ycorrm, ycorrf);
    nnzSnp.getValue(snpEffects.numNonZeros);
    nnzWind.getValue(snpEffects.numNonZeroWind);
    
    snpEffectsMale.values = snpEffects.values.col(0);
    snpEffectsFemale.values = snpEffects.values.col(1);
}

void BayesNXCIgxs::getZPZblockDiag(const Data &data){
    long numWindows = data.windSize.size();
    ZPZblockDiagMale.resize(numWindows);
    ZPZblockDiagFemale.resize(numWindows);
    for (unsigned i=0; i<numWindows; ++i) {
        unsigned starti = data.windStart[i];
        unsigned sizei  = data.windSize[i];
//        cout << i << " " << starti << " " << sizei << " " << endl;
        ZPZblockDiagMale[i] = data.Z.block(0,starti,nmale,sizei).transpose() * data.Z.block(0,starti,nmale,sizei);
        ZPZblockDiagFemale[i] = data.Z.block(nmale,starti,nfemale,sizei).transpose() * data.Z.block(nmale,starti,nfemale,sizei);
    }
}


void BayesXgxs::SnpEffects::sampleFromFC(VectorXf &ycorrm, VectorXf &ycorrf, const MatrixXf &Z,
                                         const VectorXf &ZPZdiagMale, const VectorXf &ZPZdiagFemale, const unsigned int nmale, const unsigned int nfemale,
                                         const float sigmaSq, const Vector3f &piDosage, const Vector3f &piBeta, const float varem, const float varef,
                                         VectorXf &deltaNDC, VectorXf &deltaFDC, VectorXf &deltaGxS, VectorXf &ghatm, VectorXf &ghatf){
    sumSq = 0.0;
    numNonZeros = 0;
    numSnpBetaMix.setZero(3);
    numSnpDosageMix.setZero(3);
    deltaNDC.setZero(size);
    deltaFDC.setZero(size);
    deltaGxS.setZero(size);
    
    ghatm.setZero(ycorrm.size());
    ghatf.setZero(ycorrf.size());
    
    Array2f my_rhs, rhs;   // 0: male, 1: female under NDC
    Array2f rhsFemale;     // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f rhsSame;       // 0: FDC, 1: NDC, corresponding to deltaNDC; same effect in males and females
    Array2f invLhsFemale;  // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f invLhsSame;    // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f uhatFemale;    // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f uhatSame;      // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f quadGxE;       // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f quadSame;      // 0: FDC, 1: NDC, corresponding to deltaNDC
    Array2f logDeltaNDC;      // 0: FDC, 1: NDC, corresponding to deltaNDC
    
    Array3f logPiBeta = piBeta.array().log();  // 0: null effect, 1: same effect, 2: sex-specific effect
    Array3f logPiDosage = piDosage.array().log();
    Array3f logDeltaBeta;
    Array3f logDeltaDosage;
    Array3f probDeltaBeta;
    Array3f probDeltaDosage;
    
    float oldSampleMale, oldSampleFemale;
    float sampleMale, sampleFemale;
    float logSigmaSq = log(sigmaSq);
    float invVarem = 1.0f/varem;
    float invVaref = 1.0f/varef;
    float invSigmaSq = 1.0f/sigmaSq;
    float rhsMale;
    float invLhsMale;
    float uhatMale;
    
    unsigned deltaBeta;
    unsigned deltaDosage;
    unsigned sampleDeltaNDC;
    
    for (unsigned i=0; i<size; ++i) {
        oldSampleMale   = values(i,0);
        oldSampleFemale = values(i,1);
        ycorrm += Z.col(i).head(nmale) * oldSampleMale;
        if (deltaNDC[i])
            ycorrf += Z.col(i).tail(nfemale) * oldSampleFemale;
        else
            ycorrf += Z.col(i).tail(nfemale) * oldSampleFemale * 0.5f;
        
        my_rhs[0] = Z.col(i).head(nmale).dot(ycorrm) * invVarem;
        my_rhs[1] = Z.col(i).tail(nfemale).dot(ycorrf) * invVaref;
        
        MPI_Allreduce(&my_rhs[0], &rhs[0], 2, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        
        rhsMale        = rhs[0];
        rhsFemale[0]   = rhs[1]*0.5f;
        rhsFemale[1]   = rhs[1];
        rhsSame        = rhsMale + rhsFemale;
        
        invLhsMale      = 1.0f/(ZPZdiagMale[i]*invVarem + invSigmaSq);
        invLhsFemale[0] = 1.0f/(ZPZdiagFemale[i]*0.25f*invVaref + invSigmaSq);
        invLhsFemale[1] = 1.0f/(ZPZdiagFemale[i]*invVaref       + invSigmaSq);
        invLhsSame[0]   = 1.0f/(ZPZdiagMale[i]*invVarem + ZPZdiagFemale[i]*0.25f*invVaref + invSigmaSq);
        invLhsSame[1]   = 1.0f/(ZPZdiagMale[i]*invVarem + ZPZdiagFemale[i]*invVaref       + invSigmaSq);
        
        //        cout << ZPZdiagMale[i] << " " << ZPZdiagFemale[i] << " " << invLhsMale << " " << invLhsFemale.transpose() << " " << invVarem << " " << invVaref << endl;
        
        uhatMale   = invLhsMale  *rhsMale;
        uhatFemale = invLhsFemale*rhsFemale;
        uhatSame   = invLhsSame  *rhsSame;
        
        quadGxE  = uhatMale*rhsMale + uhatFemale*rhsFemale;
        quadSame = uhatSame*rhsSame;
        
        if (uhatFemale[1]/uhatMale > 0.45 && uhatFemale[1]/uhatMale < 0.55) {
            quadGxE.setZero(2);
        }
        else if (uhatFemale[0]/uhatMale > 1.95 && uhatFemale[0]/uhatMale < 2.05) {
            quadGxE.setZero(2);
        }
        
        logDeltaNDC[1] = 0.5f*quadSame[1] + logf(sqrt(invLhsSame[1])*piBeta[1] + expf(0.5f*(quadGxE[1]-quadSame[1]))*sqrt(invLhsMale*invLhsFemale[1]*invSigmaSq)*piBeta[2] + expf(0.5f*(logSigmaSq-quadSame[1]))*piBeta[0]);
        logDeltaNDC[0] = 0.5f*quadSame[0] + logf(sqrt(invLhsSame[0])*piBeta[1] + expf(0.5f*(quadGxE[0]-quadSame[0]))*sqrt(invLhsMale*invLhsFemale[0]*invSigmaSq)*piBeta[2] + expf(0.5f*(logSigmaSq-quadSame[0]))*piBeta[0]);

        //        cout << "uhatMale " << uhatMale << " uhatFemale " << uhatFemale.transpose() << " uhatSame " << uhatSame.transpose() << endl;
        //        cout << "rhsSame " << rhsSame.transpose() << " invLhsSame " << invLhsSame.transpose() << " invVarem " << invVarem << " invVaref " << invVaref << endl;
        
        //        logLikeNDCSame = -0.5f*logSigmaSq + 0.5f*quadSame[1] + 0.5f*logf(invLhsSame[1]);
        //        logLikeFDCGxE  = -logSigmaSq + 0.5f*quadGxE[0] + 0.5f*logf(invLhsMale*invLhsFemale[0]);
        //
        
        //sample deltaDosage
        logDeltaDosage[1] = logDeltaNDC[1] + logPiDosage[1];
        logDeltaDosage[2] = logDeltaNDC[0] + logPiDosage[2];
        logDeltaDosage[0] = logf(0.5) + 0.5f*quadSame[1] + logf(sqrt(invLhsSame[1])*piBeta[1] + expf(0.5f*(quadGxE[1]-quadSame[1]))*sqrt(invLhsMale*invLhsFemale[1]*invSigmaSq)*piBeta[2] + expf(0.5f*(logSigmaSq-quadSame[1]))*piBeta[0] + expf(0.5f*(quadSame[0]-quadSame[1]))*sqrt(invLhsSame[0])*piBeta[1] + expf(0.5f*(quadGxE[0]-quadSame[0]-quadSame[1]))*sqrt(invLhsMale*invLhsFemale[0]*invSigmaSq)*piBeta[2] + expf(0.5f*(logSigmaSq-quadSame[0]-quadSame[1]))*piBeta[0]) + logPiDosage[0];
        
        for (unsigned j=0; j<3; ++j) {
            probDeltaDosage[j] = 1.0f/(logDeltaDosage-logDeltaDosage[j]).exp().sum();
        }
        
        cout << "logPiDosage " << logPiDosage.transpose() << endl;
        cout << "snp " << i << " " << logDeltaDosage.transpose() << " " << probDeltaDosage.transpose() << endl;
//        cout << sqrt(invLhsSame[1])*piBeta[1] << " " << expf(0.5f*(quadGxE[1]-quadSame[1]))*sqrt(invLhsMale*invLhsFemale[1]*invSigmaSq)*piBeta[2] << " " << expf(0.5f*(logSigmaSq-quadSame[1]))*piBeta[0] << " " << expf(0.5f*(quadSame[0]-quadSame[1]))*sqrt(invLhsSame[0])*piBeta[1] << " " << expf(0.5f*(quadGxE[0]-quadSame[0]-quadSame[1]))*sqrt(invLhsMale*invLhsFemale[0]*invSigmaSq)*piBeta[2] << " " << expf(0.5f*(logSigmaSq-quadSame[0]-quadSame[1]))*piBeta[0] << endl;
//        cout << quadSame[0] << " " << quadSame[1] << endl;
        
        deltaDosage = bernoulli.sample(probDeltaDosage);
        
        if (deltaDosage == 1) {
            sampleDeltaNDC = bernoulli.sample(probDeltaDosage[1]);
            deltaNDC[i] = 1;
        } else if (deltaDosage == 2) {
            sampleDeltaNDC = bernoulli.sample(1.0-probDeltaDosage[2]);
            deltaFDC[i] = 1;
        } else {
            sampleDeltaNDC = bernoulli.sample(0.5);
        }
        deltaNDC[i] = sampleDeltaNDC;
        ++numSnpDosageMix[deltaDosage];
        
//        probDeltaNDC1 = 1.0f/(1.0f + expf(logDeltaNDC[0] - logDeltaNDC[1]));
//        sampleDeltaNDC = bernoulli.sample(probDeltaNDC1);
        //        cout << quadSame.transpose() <<  " " << invLhsSame.transpose() << " " << quadGxE.transpose() << " " << invLhsMale << " " << invLhsFemale.transpose() << endl;
        //        sampleDeltaNDC = 1;
//        deltaNDC[i] = sampleDeltaNDC;
        //        cout << "logDeltaNDC " << logDeltaNDC.transpose() << endl;
        //        cout << sqrt(invLhsSame[1])*pis[1] << " " << (0.5f*(quadGxE[1]-quadSame[1])) << " " << sqrt(invLhsMale*invLhsFemale[1]*invSigmaSq)*pis[2] << " " << expf(0.5f*(logSigmaSq-quadSame[1]))*pis[0] << endl;
        //        cout << "probDeltaNDC1 " << probDeltaNDC1 << " sampleDeltaNDC " << sampleDeltaNDC << endl;
        //        cout << "quadGxE " << quadGxE.transpose() << " quadSame " << quadSame.transpose() << endl;
        
        // sample delta
        logDeltaBeta[0] = logPiBeta[0];
        logDeltaBeta[1] = 0.5f*(logf(invLhsSame[sampleDeltaNDC]) - logSigmaSq + quadSame[sampleDeltaNDC]) + logPiBeta[1];
        logDeltaBeta[2] = 0.5f*(logf(invLhsMale)+logf(invLhsFemale[sampleDeltaNDC]) - 2.0f*logSigmaSq + quadGxE[sampleDeltaNDC]) + logPiBeta[2];
        
        for (unsigned j=0; j<3; ++j) {
            probDeltaBeta[j] = 1.0f/(logDeltaBeta-logDeltaBeta[j]).exp().sum();
        }
        
        deltaBeta = bernoulli.sample(probDeltaBeta);
        //        cout << logDelta.transpose() << endl;
        //        cout << "probDelta " << probDelta.transpose() << " delta " << delta << endl;
        //        delta = 1;
        ++numSnpBetaMix[deltaBeta];
        
        
        if (deltaBeta) {
            // sample effect
            
            if (deltaBeta == 1) {  // same effect size in males and females
                values(i,0) = values(i,1) = sampleMale = sampleFemale = normal.sample(uhatSame[sampleDeltaNDC], invLhsSame[sampleDeltaNDC]);
                //                cout << sampleMale << " " << sampleFemale << endl;
            } else {  // different effect sizes in males and females
                values(i,0) = sampleMale = normal.sample(uhatMale, invLhsMale);
                values(i,1) = sampleFemale = normal.sample(uhatFemale[sampleDeltaNDC], invLhsFemale[sampleDeltaNDC]);
                deltaGxS[i] = 1;
            }
            
            sumSq += sampleMale * sampleMale + sampleFemale * sampleFemale;
            ++numNonZeros;
            
            ycorrm -= Z.col(i).head(nmale) * sampleMale;
            ghatm += Z.col(i).head(nmale) * sampleMale;
            
            if (deltaNDC[i]) {
                ycorrf -= Z.col(i).tail(nfemale) * sampleFemale;
                ghatf += Z.col(i).tail(nfemale) * sampleFemale;
            } else {
                ycorrf -= Z.col(i).tail(nfemale) * sampleFemale * 0.5f;
                ghatf += Z.col(i).tail(nfemale) * sampleFemale * 0.5f;
            }
        }
        else {
            values(i,0) = values(i,1) = 0.0;
        }
    }
}

void BayesXgxs::sampleUnknowns(){
    
    fixedEffects.sampleFromFC(ycorrm, ycorrf, data.X, nmale, nfemale, XPXdiagMale, XPXdiagFemale, varem.value, varef.value);
    
    snpEffects.sampleFromFC(ycorrm, ycorrf, data.Z, ZPZdiagMale, ZPZdiagFemale, nmale, nfemale, sigmaSq.value,
                            piDosage.values, piBeta.values, varem.value, varef.value,
                            deltaNDC.values, deltaFDC.values, deltaGxS.values, ghatm, ghatf);
    
    sigmaSq.sampleFromFC(snpEffects.sumSq, 2.0*snpEffects.numNonZeros);  // both male and female effects contribute to sigmaSq
    
    if(estimatePi) {
        piDosage.sampleFromFC(snpEffects.numSnpDosageMix);
        piBeta.sampleFromFC(snpEffects.numSnpBetaMix);
    }
    piNDC.value = snpEffects.numSnpDosageMix[1]/(float)snpEffects.size;
    piFDC.value = snpEffects.numSnpDosageMix[2]/(float)snpEffects.size;
    pi.value = snpEffects.numNonZeros/(float)snpEffects.size;
    piGxS.value = snpEffects.numNonZeros ? snpEffects.numSnpBetaMix[2]/(float)snpEffects.numNonZeros : 0.0;
    
    varem.sampleFromFC(ycorrm);
    varef.sampleFromFC(ycorrf);
    vargm.compute(ghatm);
    vargf.compute(ghatf);
    hsqm.compute(vargm.value, varem.value);
    hsqf.compute(vargf.value, varef.value);
    
    rounding.computeYcorr(data.y, data.X, data.Z, deltaNDC.values, nmale, nfemale, fixedEffects.values, snpEffects.values, ycorrm, ycorrf);
    nnzSnp.getValue(snpEffects.numNonZeros);
    
    snpEffectsMale.values = snpEffects.values.col(0);
    snpEffectsFemale.values = snpEffects.values.col(1);
}
