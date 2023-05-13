//
//  data.cpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "data.hpp"

// most read file methods are adopted from GCTA with modification

bool SnpInfo::isProximal(const SnpInfo &snp2, const float genWindow) const {
    return chrom == snp2.chrom && fabs(genPos - snp2.genPos) < genWindow;
}

bool SnpInfo::isProximal(const SnpInfo &snp2, const unsigned physWindow) const {
    return chrom == snp2.chrom && abs(physPos - snp2.physPos) < physWindow;
}

void AnnoInfo::getSnpInfo() {
    vector<SnpInfo*> incdSnpVec;
    unsigned numIncdSnps = 0;
    SnpInfo *snp;
    for (unsigned i=0; i<size; ++i) {
        snp = memberSnpVec[i];
        if (snp->included) {
            incdSnpVec.push_back(snp);
            snp->annoIdx.resize(snp->numAnnos);   // the index of SNP in the annotation
            for (unsigned j=0; j<snp->numAnnos; ++j) {
                if (snp->annoVec[j] == this) {
                    snp->annoIdx[j] = numIncdSnps;
                }
            }
            ++numIncdSnps;
        }
    }
    memberSnpVec = incdSnpVec;
    size = numIncdSnps;
    snp2pq.resize(size);
    for (unsigned i=0; i<size; ++i) {
        snp2pq[i] = memberSnpVec[i]->twopq;
    }
}

void AnnoInfo::print() {
    cout << boost::format("%6s %30s %12s %8.3f\n")
    % (idx+1)
    % label
    % size
    % fraction;
}

void Data::readFamFile(const string &famFile){
    // ignore phenotype column
    ifstream in(famFile.c_str());
    if (!in) throw ("Error: can not open the file [" + famFile + "] to read.");
    cout << "Reading PLINK FAM file from [" + famFile + "]." << endl;
    indInfoVec.clear();
    indInfoMap.clear();
    string fid, pid, dad, mom, sex, phen;
    unsigned idx = 0;
    while (in >> fid >> pid >> dad >> mom >> sex >> phen) {
        IndInfo *ind = new IndInfo(idx++, fid, pid, dad, mom, atoi(sex.c_str()));
        indInfoVec.push_back(ind);
        if (indInfoMap.insert(pair<string, IndInfo*>(ind->catID, ind)).second == false) {
            throw ("Error: Duplicate individual ID found: \"" + fid + "\t" + pid + "\".");
        }
    }
    in.close();
    numInds = (unsigned) indInfoVec.size();
    cout << numInds << " individuals to be included from [" + famFile + "]." << endl;
}

void Data::readBimFile(const string &bimFile) {
    // Read bim file: recombination rate is defined between SNP i and SNP i-1
    ifstream in(bimFile.c_str());
    if (!in) throw ("Error: can not open the file [" + bimFile + "] to read.");
    cout << "Reading PLINK BIM file from [" + bimFile + "]." << endl;
    snpInfoVec.clear();
    snpInfoMap.clear();
    string id, allele1, allele2;
    unsigned chr, physPos;
    float genPos;
    unsigned idx = 0;
    while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2) {
        SnpInfo *snp = new SnpInfo(idx++, id, allele1, allele2, chr, genPos, physPos);
        snpInfoVec.push_back(snp);
        chromosomes.insert(snp->chrom);
        if (snpInfoMap.insert(pair<string, SnpInfo*>(id, snp)).second == false) {
            throw ("Error: Duplicate SNP ID found: \"" + id + "\".");
        }
    }
    in.close();
    numSnps = (unsigned) snpInfoVec.size();
    cout << numSnps << " SNPs to be included from [" + bimFile + "]." << endl;
}

/*void Data::readBedFile(const string &bedFile){
    unsigned i = 0, j = 0, k = 0;
    
    if (numIncdSnps == 0) throw ("Error: No SNP is retained for analysis.");
    if (numKeptInds == 0) throw ("Error: No individual is retained for analysis.");
    
    Z.resize(numKeptInds, numIncdSnps);
    ZPZdiag.resize(numIncdSnps);
    snp2pq.resize(numIncdSnps);
    
    // Read bed file
    char ch[1];
    bitset<8> b;
    unsigned allele1=0, allele2=0;
    ifstream BIT(bedFile.c_str(), ios::binary);
    if (!BIT) throw ("Error: can not open the file [" + bedFile + "] to read.");
    cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
    for (i = 0; i < 3; i++) BIT.read(ch, 1); // skip the first three bytes
    SnpInfo *snpInfo = NULL;
    unsigned snp = 0, ind = 0;
    unsigned nmiss = 0;
    float mean = 0.0;
    for (j = 0, snp = 0; j < numSnps; j++) { // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 10: hetezygote; 01: missing
        snpInfo = snpInfoVec[j];
        mean = 0.0;
        nmiss = 0;
        if (!snpInfo->included) {
            for (i = 0; i < numInds; i += 4) BIT.read(ch, 1);
            continue;
        }
        for (i = 0, ind = 0; i < numInds;) {
            BIT.read(ch, 1);
            if (!BIT) throw ("Error: problem with the BED file ... has the FAM/BIM file been changed?");
            b = ch[0];
            k = 0;
            while (k < 7 && i < numInds) {
                if (!indInfoVec[i]->kept) k += 2;
                else {
                    allele1 = (!b[k++]);
                    allele2 = (!b[k++]);
                    if (allele1 == 0 && allele2 == 1) {  // missing genotype
                        Z(ind++, snp) = -9;
                        ++nmiss;
                    } else {
                        mean += Z(ind++, snp) = allele1 + allele2;
                    }
                }
                i++;
            }
        }

        // fill missing values with the mean
        mean /= float(numKeptInds-nmiss);
        if (nmiss) {
            for (i=0; i<numKeptInds; ++i) {
                if (Z(i,snp) == -9) Z(i,snp) = mean;
            }
        }
        
        // compute allele frequency
        snpInfo->af = 0.5f*mean;
        snp2pq[snp] = 2.0f*snpInfo->af*(1.0f-snpInfo->af);

        //cout << "snp " << snp << "     " << Z.col(snp).sum() << endl;

        if (++snp == numIncdSnps) break;
    }
    BIT.clear();
    BIT.close();
    
    
    // standardize genotypes
    for (i=0; i<numIncdSnps; ++i) {
        Z.col(i).array() -= Z.col(i).mean();
        //Z.col(i).array() /= sqrtf(gadgets::calcVariance(Z.col(i))*numKeptInds);
        ZPZdiag[i] = Z.col(i).squaredNorm();
    }
    
    //cout << "Z" << endl << Z << endl;
    
    cout << "Genotype data for " << numKeptInds << " individuals and " << numIncdSnps << " SNPs are included from [" + bedFile + "]." << endl;
}*/

void Data::readBedFile(const bool noscale, const string &bedFile){
    unsigned i = 0, j = 0;
    
    if (numIncdSnps == 0) throw ("Error: No SNP is retained for analysis.");
    if (numKeptInds == 0) throw ("Error: No individual is retained for analysis.");
    
    Z.resize(numKeptInds, numIncdSnps);
    ZPZdiag.resize(numIncdSnps);
    snp2pq.resize(numIncdSnps);
    
    // Read bed file
    FILE *in = fopen(bedFile.c_str(), "rb");
    if (!in) throw ("Error: can not open the file [" + bedFile + "] to read.");
    cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
    char header[3];
    fread(header, sizeof(header), 1, in);
    if (!in || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
        cerr << "Error: Incorrect first three bytes of bed file: " << bedFile << endl;
        exit(1);
    }
    
    // Read genotypes
    SnpInfo *snpInfo = NULL;
    IndInfo *indInfo = NULL;
    unsigned snp = 0;
    unsigned nmiss=0;
    float sum=0.0, mean=0.0;
    
    const int bedToGeno[4] = {2, -9, 1, 0};
    unsigned size = (numInds+3)>>2;
    int genoValue;
    unsigned long long skip = 0;
    
    for (j = 0, snp = 0; j < numSnps; j++) {  // code adopted from BOLT-LMM with modification
        snpInfo = snpInfoVec[j];
        sum = 0.0;
        nmiss = 0;
        
        if (!snpInfo->included) {
//            in.ignore(size);
            skip += size;
            continue;
        }
        
        if (skip) fseek(in, skip, SEEK_CUR);
        skip = 0;
 
        char *bedLineIn = new char[size];
        fread(bedLineIn, 1, size, in);

        for (i = 0; i < numInds; i++) {
            indInfo = indInfoVec[i];
            if (!indInfo->kept) continue;
            genoValue = bedToGeno[(bedLineIn[i>>2]>>((i&3)<<1))&3];
            
            Z(indInfo->index, snp) = genoValue;
            if (genoValue == -9) ++nmiss;   // missing genotype
            else sum += genoValue;
        }
        delete[] bedLineIn;
                
        // fill missing values with the mean
        mean = sum/float(numKeptInds - nmiss);
        if (nmiss) {
            for (i=0; i<numKeptInds; ++i) {
                if (Z(i,snp) == -9) Z(i,snp) = mean;
            }
        }
        
        // compute allele frequency
        snpInfo->af = 0.5f*mean;
        snp2pq[snp] = snpInfo->twopq = 2.0f*snpInfo->af*(1.0f-snpInfo->af);
        
        //cout << "snp " << snp << "     " << Z.col(snp).sum() << " " << mean_all << endl;
        
        Z.col(snp).array() -= mean; // center column by 2p rather than the real mean
        if (!noscale) Z.col(snp).array() /= sqrt(snp2pq[snp]);  // standardise to have variance one
        
        Z.col(snp).array() *= RinverseSqrt.array();

        if (++snp == numIncdSnps) break;
    }
    fclose(in);
    
    // standardize genotypes
//    VectorXf colsums = Z.colwise().sum();
    //Z.rowwise() -= colsums.transpose()/numKeptInds;  // center
    
    ZPZdiag = Z.colwise().squaredNorm();
    
    cout << "Genotype data for " << numKeptInds << " individuals and " << numIncdSnps << " SNPs are included from [" + bedFile + "]." << endl;
}

void Data::readPhenotypeFile(const string &phenFile, const unsigned mphen) {
    // NA: missing phenotype
    ifstream in(phenFile.c_str());
    if (!in) throw ("Error: can not open the phenotype file [" + phenFile + "] to read.");
    cout << "Reading phenotypes from [" + phenFile + "]." << endl;
    map<string, IndInfo*>::iterator it, end=indInfoMap.end();
    IndInfo *ind = NULL;
    Gadget::Tokenizer colData;
    string inputStr;
    string sep(" \t");
    string id;
    unsigned line=0;
    while (getline(in,inputStr)) {
        colData.getTokens(inputStr, sep);
        id = colData[0] + ":" + colData[1];
        it = indInfoMap.find(id);
        if (it != end && colData[mphen+1] != "NA") {
            ind = it->second;
            ind->phenotype = atof(colData[mphen+1].c_str());
            ++line;
        }
    }
    in.close();
    cout << "Non-missing phenotypes of trait " << mphen << " of " << line << " individuals are included from [" + phenFile + "]." << endl;
}

void Data::readCovariateFile(const string &covarFile){
    if (covarFile.empty()) return;
    ifstream in(covarFile.c_str());
    if (!in) throw ("Error: can not open the file [" + covarFile + "] to read.");
    map<string, IndInfo*>::iterator it, end=indInfoMap.end();
    IndInfo *ind = NULL;
    Gadget::Tokenizer colData;
    string inputStr;
    string sep(" \t");
    string id;
    unsigned line=0;
    unsigned numCovariates=0;
    while (getline(in,inputStr)) {
        colData.getTokens(inputStr, sep);
        if (line==0) {
            numCovariates = (unsigned)colData.size() - 2;
            numFixedEffects = numCovariates + 1;
            fixedEffectNames.resize(numFixedEffects);
            fixedEffectNames[0] = "Intercept";
            for (unsigned i=0; i<numCovariates; ++i) {
                fixedEffectNames[i+1] = colData[i+2];
            }
            ++line;
        }
        id = colData[0] + ":" + colData[1];
        it = indInfoMap.find(id);
        if (it != end) {
            ind = it->second;
            ind->covariates.resize(numCovariates + 1);  // plus intercept
            ind->covariates[0] = 1;
            for (unsigned i=2; i<colData.size(); ++i) {
                ind->covariates[i-1] = atof(colData[i].c_str());
            }
            ++line;
        }
    }
    in.close();
    
    cout << "Read " << numCovariates << " covariates from [" + covarFile + "]." << endl;
}

void Data::readRandomCovariateFile(const string &covarFile){
    if (covarFile.empty()) return;
    ifstream in(covarFile.c_str());
    if (!in) throw ("Error: can not open the file [" + covarFile + "] to read.");
    map<string, IndInfo*>::iterator it, end=indInfoMap.end();
    IndInfo *ind = NULL;
    Gadget::Tokenizer colData;
    string inputStr;
    string sep(" \t");
    string id;
    unsigned line=0;
    while (getline(in,inputStr)) {
        colData.getTokens(inputStr, sep);
        if (line==0) {
            numRandomEffects = (unsigned)colData.size() - 2;
            randomEffectNames.resize(numRandomEffects);
            for (unsigned i=0; i<numRandomEffects; ++i) {
                randomEffectNames[i] = colData[i+2];
            }
            ++line;
        }
        id = colData[0] + ":" + colData[1];
        it = indInfoMap.find(id);
        if (it != end) {
            ind = it->second;
            ind->randomCovariates.resize(numRandomEffects);
            for (unsigned i=0; i<numRandomEffects; ++i) {
                ind->randomCovariates[i] = atof(colData[i+2].c_str());
            }
            ++line;
        }
    }
    in.close();
    
    cout << "Read " << numRandomEffects << " covariates as random effects from [" + covarFile + "]." << endl;
}

void Data::readResidualDiagFile(const string &resDiagFile){
    if (resDiagFile.empty()) return;
    ifstream in(resDiagFile.c_str());
    if (!in) throw ("Error: can not open the file [" + resDiagFile + "] to read.");
    map<string, IndInfo*>::iterator it, end=indInfoMap.end();
    IndInfo *ind = NULL;
    string fid, pid, resDiag;
    string id;
    unsigned line=0;
    while (in >> fid >> pid >> resDiag) {
        id = fid + ":" + pid;
        it = indInfoMap.find(id);
        if (it != end) {
            ind = it->second;
            ind->rinverse = 1.0f/atof(resDiag.c_str());
            ++line;
        }
    }
    in.close();
    weightedRes = true;
    
    cout << "Read residual diagonal values for " << line << " individuals from [" + resDiagFile + "]." << endl;
}

void Data::keepMatchedInd(const string &keepIndFile, const unsigned keepIndMax){  // keepIndFile is optional
    map<string, IndInfo*>::iterator it, end=indInfoMap.end();
    IndInfo *ind = NULL;
    vector<string> keep;
    keep.reserve(numInds);
    unsigned cnt=0;
    
    if (numFixedEffects) {
        for (unsigned i=0; i<numInds; ++i) {
            ind = indInfoVec[i];
            ind->kept = false;
            if (numRandomEffects) {
                if (ind->phenotype!=-9 && ind->covariates.size() && ind->randomCovariates.size()) {
                    if (keepIndMax > cnt++)
                        keep.push_back(ind->catID);
                }
            }
            else {
                if (ind->phenotype!=-9 && ind->covariates.size()) {
                    if (keepIndMax > cnt++)
                        keep.push_back(ind->catID);
                }
            }
        }
    }
    else {
        numFixedEffects = 1;
        fixedEffectNames = {"Intercept"};
        for (unsigned i=0; i<numInds; ++i) {
            ind = indInfoVec[i];
            ind->kept = false;
            if (numRandomEffects) {
                if (ind->phenotype!=-9 && ind->randomCovariates.size()) {
                    ind->covariates.resize(1);
                    ind->covariates << 1;
                    if (keepIndMax > cnt++)
                        keep.push_back(ind->catID);
                }
            }
            else {
                if (ind->phenotype!=-9) {
                    ind->covariates.resize(1);
                    ind->covariates << 1;
                    if (keepIndMax > cnt++)
                        keep.push_back(ind->catID);
                }
            }
        }
    }
    
    if (!keepIndFile.empty()) {
        ifstream in(keepIndFile.c_str());
        if (!in) throw ("Error: can not open the file [" + keepIndFile + "] to read.");
        string fid, pid;
        keep.clear();
        while (in >> fid >> pid) {
            keep.push_back(fid + ":" + pid);
        }
        in.close();
    }
    
    numKeptInds = 0;

    for (unsigned i=0; i<keep.size(); ++i) {
        it = indInfoMap.find(keep[i]);
        if (it == end) {
            Gadget::Tokenizer token;
            token.getTokens(keep[i], ":");
            throw("Error: Individual " + token[0] + " " + token[1] + " from file [" + keepIndFile + "] does not exist!");
        } else {
            ind = it->second;
            if (ind->phenotype != -9) {
                ind->kept = true;
            } else {
                throw("Error: Individual " + ind->famID + " " + ind->indID + " from file [" + keepIndFile + "] does not have phenotype!");
            }
        }
    }
    
    keptIndInfoVec = makeKeptIndInfoVec(indInfoVec);
    numKeptInds =  (unsigned) keptIndInfoVec.size();
    
    RinverseSqrt.setOnes(numKeptInds);
    Rsqrt.setOnes(numKeptInds);
    for (unsigned i=0; i<numKeptInds; ++i) {
        RinverseSqrt[i] = sqrt(keptIndInfoVec[i]->rinverse);
        Rsqrt[i] = 1.0f/RinverseSqrt[i];
    }
    
    y.setZero(numKeptInds);
    for (unsigned i=0; i<numKeptInds; ++i) {
        y[i] = keptIndInfoVec[i]->phenotype;
    }
    y.array() *= RinverseSqrt.array();
    ypy = (y.array()-y.mean()).square().sum();
    varPhenotypic = ypy/numKeptInds;
    
    X.resize(numKeptInds, numFixedEffects);
    for (unsigned i=0; i<numKeptInds; ++i) {
        X.row(i) = keptIndInfoVec[i]->covariates.array() * RinverseSqrt[i];
    }
    XPXdiag = X.colwise().squaredNorm();
    
    if (numRandomEffects) {
        W.resize(numKeptInds, numRandomEffects);
        for (unsigned i=0; i<numKeptInds; ++i) {
            W.row(i) = keptIndInfoVec[i]->randomCovariates.array() * RinverseSqrt[i];
        }
        WPWdiag = W.colwise().squaredNorm();
    }
    
    cout << numKeptInds << " matched individuals are kept." << endl;
}

void Data::initVariances(const float heritability, const float propVarRandom){
//    float varPhenotypic = ypy/numKeptInds;;
    varGenotypic = varPhenotypic * heritability;
    varResidual  = varPhenotypic - varGenotypic;
    varRandom    = varPhenotypic * propVarRandom;
//    cout <<ypy<<" "<<numKeptInds<<" "<<varPhenotypic<<" " <<varGenotypic << " " <<varResidual << endl;
}

void Data::includeSnp(const string &includeSnpFile){
    ifstream in(includeSnpFile.c_str());
    if (!in) throw ("Error: can not open the file [" + includeSnpFile + "] to read.");
    for (unsigned i=0; i<numSnps; ++i) {
        snpInfoVec[i]->included = false;
    }
    map<string, SnpInfo*>::iterator it, end = snpInfoMap.end();
    string id;
    while (in >> id) {
        it = snpInfoMap.find(id);
        if (it != end) {
            it->second->included = true;
        }
    }
    in.close();
}

void Data::excludeSnp(const string &excludeSnpFile){
    ifstream in(excludeSnpFile.c_str());
    if (!in) throw ("Error: can not open the file [" + excludeSnpFile + "] to read.");
    map<string, SnpInfo*>::iterator it, end = snpInfoMap.end();
    string id;
    while (in >> id) {
        it = snpInfoMap.find(id);
        if (it != end) {
            it->second->included = false;
        }
    }
    in.close();
}

void Data::includeChr(const unsigned chr){
    if (!chr) return;
    for (unsigned i=0; i<numSnps; ++i){
        SnpInfo *snp = snpInfoVec[i];
        if (snp->chrom != chr) snp->included = false;
    }
}

void Data::includeBlock(const unsigned block){
    if (!block) return;
    if (block > ldBlockInfoVec.size()) throw("Error: Reqest to include block " + to_string(block) + " but there are only " + to_string(ldBlockInfoVec.size()) + " in total!");
    for (unsigned i=0; i<numLDBlocks; ++i) {
        LDBlockInfo *blockInfo = ldBlockInfoVec[i];
        if (block != i+1) blockInfo->kept = false;
    }
    LDBlockInfo *blockInfo = ldBlockInfoVec[block-1];
    unsigned cnt = 0;
    if (blockInfo->numSnpInBlock) {
        for (unsigned i=0; i<numSnps; ++i){
            SnpInfo *snpInfo = snpInfoVec[i];
            snpInfo->included = false;
        }
        for (unsigned i=0; i<blockInfo->numSnpInBlock; ++i) {
            SnpInfo *snpInfo = blockInfo->snpInfoVec[i];
            snpInfo->included = true;
            ++cnt;
        }
    }
    else {
        for (unsigned i=0; i<numSnps; ++i){
            SnpInfo *snpInfo = snpInfoVec[i];
            if (snpInfo->chrom != blockInfo->chrom) snpInfo->included = false;
            else if (snpInfo->physPos < blockInfo->startPos) snpInfo->included = false;
            else if (snpInfo->physPos > blockInfo->endPos) snpInfo->included = false;
            else ++cnt;
        }
    }
    cout << "Included " << cnt << " SNPs in block " << blockInfo->ID << endl;
}

void Data::includeSkeletonSnp(const string &skeletonSnpFile){
    ifstream in(skeletonSnpFile.c_str());
    if (!in) throw ("Error: can not open the file [" + skeletonSnpFile + "] to read.");
    map<string, SnpInfo*>::iterator it, end = snpInfoMap.end();
    string id;
    SnpInfo *snp;
    numSkeletonSnps = 0;
    while (in >> id) {
        it = snpInfoMap.find(id);
        if (it != end) {
            snp = it->second;
            snp->included = true;
            snp->skeleton = true;
            ++numSkeletonSnps;
        }
    }
    cout << numSkeletonSnps << " skeleton SNPs are included." << endl;
    in.close();
}

void Data::excludeMHC(){
    long cnt = 0;
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        if (snp->chrom == 6) {
            if (snp->physPos > 28e6 && snp->physPos < 34e6) {
                snp->included = false;
                ++cnt;
            }
        }
    }
    cout << cnt << " SNPs in the MHC region (Chr6:28-34Mb) are excluded." << endl;
}

void Data::excludeAmbiguousSNP(){
    long cnt = 0;
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        if ((snp->a1 == "T" && snp->a2 == "A") ||
            (snp->a1 == "A" && snp->a2 == "T") ||
            (snp->a1 == "G" && snp->a2 == "C") ||
            (snp->a1 == "C" && snp->a2 == "G")) {
            snp->included = false;
            ++cnt;
        }
    }
    cout << cnt << " SNPs with ambiguous nucleotides, i.e. A/T or G/C, are excluded." << endl;
}

void Data::excludeSNPwithMaf(const float mafmin, const float mafmax){
    unsigned cntmin=0, cntmax=0;
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        float maf = snp->af < 0.5 ? snp->af : 1.0 - snp->af;
        if (mafmin && maf < mafmin) {
            snp->included = false;
            ++cntmin;
        }
        if (mafmax && maf > mafmax) {
            snp->included = false;
            ++cntmax;
        }
    }
    if (mafmin) cout << cntmin << " SNPs with MAF below " << mafmin << " are excluded." << endl;
    if (mafmax) cout << cntmax << " SNPs with MAF above " << mafmax << " are excluded." << endl;
}

void Data::excludeRegion(const string &excludeRegionFile){
    ifstream in(excludeRegionFile.c_str());
    if (!in) throw ("Error: can not open the file [" + excludeRegionFile + "] to read.");

    map<int, vector<pair<long, long> > > regions;
    int chr;
    long regStart, regEnd;
    while (in >> chr >> regStart >> regEnd) {
        regions[chr].push_back(pair<long, long>(regStart, regEnd));
    }
    
    unsigned cnt = 0;
    map<int, vector<pair<long, long> > >::iterator it, end = regions.end();
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        it = regions.find(snp->chrom);
        if (it != end) {
            long size = it->second.size();
            for (unsigned j=0; j<size; ++j) {
                regStart = it->second[j].first;
                regEnd   = it->second[j].second;
                if (snp->physPos > regStart && snp->physPos < regEnd) {
                    snp->included = false;
                    ++cnt;
                }
            }
        }
    }
    cout << cnt << " SNPs are excluded due to --exclude-region." << endl;
}

void Data::reindexSnp(vector<SnpInfo*> snpInfoVec){
    SnpInfo *snp;
    for (unsigned i=0, idx=0; i<snpInfoVec.size(); ++i) {
        snp = snpInfoVec[i];
        if (snp->included) {
            snp->index = idx++;
        } else {
            snp->index = -9;
        }
    }
}

void Data::includeMatchedSnp(){
    reindexSnp(snpInfoVec);  // reindex for MPI purpose in terms of full snplist
    fullSnpFlag.resize(numSnps);
    for (int i=0; i<numSnps; ++i) fullSnpFlag[i] = snpInfoVec[i]->included; // for output purpose
    
    incdSnpInfoVec = makeIncdSnpInfoVec(snpInfoVec);
    numIncdSnps = (unsigned) incdSnpInfoVec.size();
    reindexSnp(incdSnpInfoVec);
    snp2pq.resize(numIncdSnps);
    
    map<int, vector<SnpInfo*> > chrmap;
    map<int, vector<SnpInfo*> >::iterator it;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        if (chrmap.find(snp->chrom) == chrmap.end()) {
            chrmap[snp->chrom] = *new vector<SnpInfo*>;
        }
        chrmap[snp->chrom].push_back(snp);
    }
    numChroms = (unsigned) chrmap.size();
    chromInfoVec.clear();
    for (it=chrmap.begin(); it!=chrmap.end(); ++it) {
        int id = it->first;
        vector<SnpInfo*> &vec = it->second;
        ChromInfo *chr = new ChromInfo(id, (unsigned)vec.size(), vec[0]->index, vec.back()->index);
        chromInfoVec.push_back(chr);
        //cout << "size chrom " << id << ": " << vec.back()->physPos - vec[0]->physPos << endl;
    }

    cout << numIncdSnps << " SNPs on " << numChroms << " chromosomes are included." << endl;
    
//    if (numAnnos) setAnnoInfoVec();

//    if (numAnnos) {
//        if (myMPI::rank==0) cout << "\nAnnotation info:" << endl;
//        numSnpAnnoVec.resize(numAnnos);
//        for (unsigned i=0; i<numAnnos; ++i) {
//            AnnoInfo *anno = annoInfoVec[i];
//            numSnpAnnoVec[i] = anno->size;
//            anno->getSnpInfo();
//            anno->fraction = float(anno->size)/float(numIncdSnps);
//            anno->print();
//        }
//        if (myMPI::rank==0) cout << endl;
//    }
}

vector<SnpInfo*> Data::makeIncdSnpInfoVec(const vector<SnpInfo*> &snpInfoVec){
    vector<SnpInfo*> includedSnps;
    includedSnps.reserve(numSnps);
    snpEffectNames.reserve(numSnps);
    SnpInfo *snp = NULL;
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        if(snp->included) {
            //snp->index = j++;  // reindex snps
            includedSnps.push_back(snp);
            snpEffectNames.push_back(snp->ID);
        }
    }
    return includedSnps;
}

vector<IndInfo*> Data::makeKeptIndInfoVec(const vector<IndInfo*> &indInfoVec){
    vector<IndInfo*> keptInds;
    keptInds.reserve(numInds);
    IndInfo *ind = NULL;
    for (unsigned i=0, j=0; i<numInds; ++i) {
        ind = indInfoVec[i];
        if(ind->kept) {
            ind->index = j++;  // reindex inds
            keptInds.push_back(ind);
        }
    }
    return keptInds;
}



void Data::computeAlleleFreq(const MatrixXf &Z, vector<SnpInfo*> &incdSnpInfoVec, VectorXf &snp2pq){
    cout << "Computing allele frequencies ..." << endl;
    snp2pq.resize(numIncdSnps);
    SnpInfo *snp = NULL;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snp = incdSnpInfoVec[i];
        snp->af = 0.5f*Z.col(i).mean();
        snp2pq[i] = snp->twopq = 2.0f*snp->af*(1.0f-snp->af);
    }
}

void Data::getWindowInfo(const vector<SnpInfo*> &incdSnpInfoVec, const unsigned windowWidth, VectorXi &windStart, VectorXi &windSize){
    cout << "Creating windows (window width: " + to_string(static_cast<long long>(windowWidth/1e6)) + "Mb) ..." << endl;
    int i=0, j=0;
    windStart.setZero(numIncdSnps);
    windSize.setZero(numIncdSnps);
    SnpInfo *snpi, *snpj;
    for (i=0; i<numIncdSnps; ++i) {
        snpi = incdSnpInfoVec[i];
        snpi->resetWindow();
        for (j=i; j>=0; --j) {
            snpj = incdSnpInfoVec[j];
            if (snpi->isProximal(*snpj, windowWidth/2)) {
                snpi->windStart = snpj->index;
                snpi->windSize++;
            } else break;
        }
        for (j=i+1; j<numIncdSnps; ++j) {
            snpj = incdSnpInfoVec[j];
            if (snpi->isProximal(*snpj, windowWidth/2)) {
                snpi->windSize++;
            } else break;
        }
        if(!(i%10000))
            cout << "SNP " << i << " Window Size " << snpi->windSize << endl;
        if (!snpi->windSize) {
            throw("Error: SNP " + snpi->ID + " has zero SNPs in its window!");
        }
        windStart[i] = snpi->windStart;
        windSize [i] = snpi->windSize;
        snpi->windEnd = snpi->windStart + snpi->windSize - 1;
    }
}

void Data::getNonoverlapWindowInfo(const unsigned windowWidth){
    if (!windowWidth) throw("Error: Did you forget to set window width by --wind [Mb]?");
    unsigned window = 0;
    unsigned currChr = incdSnpInfoVec[0]->chrom;
    unsigned long startPos = incdSnpInfoVec[0]->physPos;
    vector<int> windStartVec = {0};
    SnpInfo *snp;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snp = incdSnpInfoVec[i];
        if (snp->physPos - startPos > windowWidth || snp->chrom > currChr) {
            currChr = snp->chrom;
            startPos = snp->physPos;
            windStartVec.push_back(i);
            ++window;
        }
        snp->window = window;
    }
    
    numWindows = windStartVec.size();

    makeWindows = true;
    
    windStart = VectorXi::Map(&windStartVec[0], numWindows);
    windSize.setZero(numWindows);

    for (unsigned i=0; i<numWindows; ++i) {
        if (i != numWindows-1)
            windSize[i] = windStart[i+1] - windStart[i];
        else
            windSize[i] = numIncdSnps - windStart[i];
    }
    
    cout << "Created " << numWindows << " non-overlapping " << windowWidth/1e3 << "kb windows with average size of " << windSize.sum()/float(numWindows) << " SNPs." << endl;
}


void Data::buildSparseMME(const string &bedFile, const unsigned windowWidth){
    cout << "Building sparse MME ..." << endl;
    
    getWindowInfo(incdSnpInfoVec, windowWidth, windStart, windSize);
    
    //cout << "windStart " << windStart.transpose() << endl;
    //cout << "windSize " << windSize.transpose() << endl;
    
    if (numIncdSnps == 0) throw ("Error: No SNP is retained for analysis.");
    if (numKeptInds == 0) throw ("Error: No individual is retained for analysis.");
    
    ZPZ.resize(numIncdSnps);
    for (unsigned i=0; i<numIncdSnps; ++i) {
        ZPZ[i].resize(windSize[i]);
    }
    ZPZdiag.resize(numIncdSnps);
    ZPX.resize(numIncdSnps, numFixedEffects);
    ZPy.resize(numIncdSnps);
    
    Gadget::Timer timer;
    timer.setTime();
    
    const int bedToGeno[4] = {2, -9, 1, 0};
    
#pragma omp parallel for
    for (unsigned chr=0; chr<numChroms; ++chr) {
        
        // Read bed file
        VectorXf genotypes(numKeptInds);
        ifstream in(bedFile.c_str(), ios::binary);
        if (!in) throw ("Error: can not open the file [" + bedFile + "] to read.");
        if (chr==0)
            cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
        char header[3];
        in.read((char *) header, 3);
        if (!in || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
            cerr << "Error: Incorrect first three bytes of bed file: " << bedFile << endl;
            exit(1);
        }

        ChromInfo *chrinfo = chromInfoVec[chr];
        
        unsigned start = chrinfo->startSnpIdx;
        unsigned end = chrinfo->endSnpIdx;
        unsigned lastWindStart = chrinfo->startSnpIdx;
        
        //cout << "thread " << omp_get_thread_num() << " chr " << chr << " snp " << start << "-" << chrinfo->endSnpIdx << endl;

        IndInfo *indi = NULL;
        SnpInfo *snpj = NULL;
        SnpInfo *snpk = NULL;

        int genoValue;
        unsigned i, j, k;
        unsigned inc; // index of included SNP
        
        for (j = 0, inc = start; j < numSnps; j++) {

            unsigned size = (numInds+3)>>2;
            
            snpj = snpInfoVec[j];

            if (snpj->index < start || !snpj->included) {
                in.ignore(size);
                continue;
            }
            
            char *bedLineIn = new char[size];
            in.read((char *)bedLineIn, size);

            if(!(inc%1000)) {
                cout << " thread " << omp_get_thread_num() << " read snp " << inc << " windStart " << snpj->windStart << " windSize " << snpj->windSize << endl;
            }

            float mean = 0.0;
            unsigned nmiss = 0;

            for (i = 0; i < numInds; i++) {
                indi = indInfoVec[i];
                if (!indi->kept) continue;
                genoValue = bedToGeno[(bedLineIn[i>>2]>>((i&3)<<1))&3];
                genotypes[indi->index] = genoValue;
                if (genoValue == -9) ++nmiss;   // missing genotype
                else mean += genoValue;
            }
            delete[] bedLineIn;
            
            // fill missing values with the mean
            mean /= float(numKeptInds-nmiss);
            if (nmiss) {
                for (i=0; i<numKeptInds; ++i) {
                    if (genotypes[i] == -9) genotypes[i] = mean;
                }
            }
            
            // compute allele frequency
            snpj->af = 0.5f*mean;
            snp2pq[inc] = snpj->twopq = 2.0f*snpj->af*(1.0f-snpj->af);
            
            // center genotypes
            genotypes.array() -= genotypes.mean();
            snpj->genotypes = genotypes;
            
            // compute Zj'Z[j] with Z[j] for genotype matrix of SNPs in the window of SNP j
            ZPZdiag[inc] = ZPZ[inc][inc - snpj->windStart] = genotypes.squaredNorm();
            for (k = snpj->windStart; k<inc; ++k) {
                snpk = incdSnpInfoVec[k];
                ZPZ[inc][k - snpj->windStart] = ZPZ[k][inc - snpk->windStart] = genotypes.dot(snpk->genotypes);
            }
            
            // release memory for genotypes of anterior SNPs of the window
            if (lastWindStart != snpj->windStart) {
                for (k=lastWindStart; k<snpj->windStart; ++k) {
                    incdSnpInfoVec[k]->genotypes.resize(0);
                }
                lastWindStart = snpj->windStart;
            }
            
            // compute Zj'X
            ZPX.row(inc) = genotypes.transpose()*X;
            
            // compute Zj'y
            ZPy[inc] = genotypes.dot(y);
            
            if (inc++ == end) break;
        }

        in.close();
    }

    n.setConstant(numIncdSnps, numKeptInds);
    tss.setConstant(numIncdSnps, ypy);

    timer.getTime();
    
    cout << "Average window size " << windSize.sum()/numIncdSnps << endl;
    cout << "Genotype data for " << numKeptInds << " individuals and " << numIncdSnps << " SNPs are included from [" + bedFile + "]." << endl;
    cout << "Construction of sparse MME completed (time used: " << timer.format(timer.getElapse()) << ")" << endl;
    
//    for (unsigned i=0; i<ZPZ.size(); ++i) {
//        cout << i << " " << ZPZ[i].transpose() << endl;
//    }
    
//    cout << "ZPZdiag " << ZPZdiag.transpose() << endl;
//    cout << "ZPZ.back() " << ZPZ.back().transpose() << endl;
//    cout << "ZPy " << ZPy.transpose() << endl;
//
//    string outfile = bedFile + ".ma";
//    ofstream out(outfile.c_str());
//    for (unsigned i=0; i<ZPy.size(); ++i) {
//        out << incdSnpInfoVec[i]->ID << "   " << setprecision(12) << ZPy[i]/ZPZdiag[i] << endl;
//    }
//    out.close();
}

void Data::outputSnpResults(const VectorXf &posteriorMean, const VectorXf &posteriorSqrMean, const VectorXf &pip, const bool noscale, const string &filename) const {
    ofstream out(filename.c_str());
    out << boost::format("%6s %20s %6s %12s %6s %6s %12s %12s %12s %8s")
    % "Id"
    % "Name"
    % "Chrom"
    % "Position"
    % "A1"
    % "A2"
    % "A1Frq"
    % "A1Effect"
    % "SE"
    % "PIP";
    if (makeWindows) out << boost::format("%8s") % "Window";
    out << endl;
    for (unsigned i=0, idx=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        if(!fullSnpFlag[i]) continue;
        //        if(snp->isQTL) continue;)
        float sqrt2pq = sqrt(2.0*snp->af*(1.0-snp->af));
        float effect = (snp->flipped ? -posteriorMean[idx] : posteriorMean[idx]);
        float se = sqrt(posteriorSqrMean[idx]-posteriorMean[idx]*posteriorMean[idx]);
        out << boost::format("%6s %20s %6s %12s %6s %6s %12.6f %12.6f %12.6f %8.8f")
        % (idx+1)
        % snp->ID
        % snp->chrom
        % snp->physPos
        % (snp->flipped ? snp->a2 : snp->a1)
        % (snp->flipped ? snp->a1 : snp->a2)
        % (snp->flipped ? 1.0-snp->af : snp->af)
        % (noscale ? effect : effect/sqrt2pq)
        % (noscale ? se : se/sqrt2pq)
        % pip[idx];
        if (makeWindows) out << boost::format("%8s") % snp->window;
        out << endl;
        ++idx;
    }
    out.close();
}

void Data::outputSnpResults(const VectorXf &posteriorMean, const VectorXf &posteriorSqrMean, const VectorXf &lastSample, const VectorXf &pip, const bool noscale, const string &filename) const {
    ofstream out(filename.c_str());
    out << boost::format("%6s %20s %6s %12s %6s %6s %12s %12s %12s %14s %14s")
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
    % "LastSampleEff";
    if (makeWindows) out << boost::format("%8s") % "Window";
    out << endl;
    for (unsigned i=0, idx=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        if(!fullSnpFlag[i]) continue;
        //        if(snp->isQTL) continue;)
        float sqrt2pq = sqrt(2.0*snp->af*(1.0-snp->af));
        float effect = (snp->flipped ? -posteriorMean[idx] : posteriorMean[idx]);
        float lastBeta = (snp->flipped ? -lastSample[idx] : lastSample[idx]);
        float se = sqrt(posteriorSqrMean[idx]-posteriorMean[idx]*posteriorMean[idx]);
        out << boost::format("%6s %20s %6s %12s %6s %6s %12.6f %12.6f %12.6f %14.8f %14.6f")
        % (idx+1)
        % snp->ID
        % snp->chrom
        % snp->physPos
        % (snp->flipped ? snp->a2 : snp->a1)
        % (snp->flipped ? snp->a1 : snp->a2)
        % (snp->flipped ? 1.0-snp->af : snp->af)
        % (noscale ? effect : effect/sqrt2pq)
        % (noscale ? se : se/sqrt2pq)
        % pip[idx]
        % (noscale ? lastBeta : lastBeta/sqrt2pq);
        if (makeWindows) out << boost::format("%8s") % snp->window;
        out << endl;
        ++idx;
    }
    out.close();
}

void Data::inputSnpResults(const string &snpResFile){
    ifstream in(snpResFile.c_str());
    if (!in) throw ("Error: can not open the SNP result file [" + snpResFile + "] to read.");
    cout << "Reading SNP results from [" + snpResFile + "]." << endl;
    
    SnpInfo *snp;
    map<string, SnpInfo*>::iterator it;
    string name;
    int id, chrom, pos, window;
    float freq, effect, se, pip;
    unsigned line=0, match=0;
    string header;
    getline(in, header);
    while (in >> id >> name >> chrom >> pos >> freq >> effect >> se >> pip >> window) {
        ++line;
        it = snpInfoMap.find(name);
        if (it == snpInfoMap.end()) continue;
        snp = it->second;
        if (snp->included) {
            snp->effect = effect;
            ++match;
        }
    }
    in.close();
    
    cout << match << " matched SNPs in the SNP result file (in total " << line << " SNPs)." << endl;
}

void Data::inputSnpInfoAndResults(const string &snpResFile, const string &bayesType){
    ifstream in(snpResFile.c_str());
    if (!in) throw ("Error: can not open the SNP result file [" + snpResFile + "] to read.");
    cout << "Reading SNP info and results from [" + snpResFile + "]." << endl;
    
    SnpInfo *snp;
    map<string, SnpInfo*>::iterator it;
    string name;
    int id, chrom, pos, window;
    float freq, effect, se, pip;
    unsigned line=0, match=0;
    string header;
    getline(in, header);
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        snp->included = false;
    }
    if (bayesType == "SMix") {
        float piS;
        while (in >> id >> name >> chrom >> pos >> freq >> effect >> se >> pip >> piS) {
            ++line;
            it = snpInfoMap.find(name);
            if (it == snpInfoMap.end()) {
                throw("Error: SNP " + name + " is not in the LD matrix!");
            }
            snp = it->second;
            snp->included = true;
            snp->gwas_af = freq;
            snp->effect = effect;
        }
    }
    else {
        while (in >> id >> name >> chrom >> pos >> freq >> effect >> se >> pip >> window) {
            ++line;
            //        SnpInfo *snp = new SnpInfo(id-1, name, "NA", "NA", chrom, 0, pos);
            it = snpInfoMap.find(name);
            if (it == snpInfoMap.end()) {
                throw("Error: SNP " + name + " is not in the LD matrix!");
            }
            snp = it->second;
            snp->included = true;
            snp->gwas_af = freq;
            snp->effect = effect;
            //        snpInfoVec.push_back(snp);
            //        snpInfoMap.insert(pair<string, SnpInfo*>(name, snp));
            //        chromosomes.insert(snp->chrom);
        }
    }
    in.close();
    
//    numSnps = (unsigned) snpInfoVec.size();
    cout << line << " SNPs in the SNP result file." << endl;
}


void Data::summarizeSnpResults(const SpMat &snpEffects, const string &filename) const {
    cout << "SNP results to be summarized in " << filename << endl;
    unsigned nrow = snpEffects.rows();
    VectorXf effectSum(numIncdSnps), effectMean(numIncdSnps);
    VectorXf pipSum(numIncdSnps), pip(numIncdSnps);  // posterior inclusion probability
    for (unsigned i=0; i<numIncdSnps; ++i) {
        effectSum[i] = snpEffects.col(i).sum();
        pipSum[i] = (VectorXf(snpEffects.col(i)).array()!=0).count();
    }
    effectMean = effectSum/(float)nrow;
    pip = pipSum/(float)nrow;
    
    ofstream out(filename.c_str());
    out << boost::format("%6s %20s %6s %12s %6s %6s %12s %12s %8s %8s\n")
    % "Id"
    % "Name"
    % "Chrom"
    % "Position"
    % "A1"
    % "A2"
    % "A1Frq"
    % "A1Effect"
    % "PIP"
    % "Window";
    for (unsigned i=0, idx=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        if(!fullSnpFlag[i]) continue;
        out << boost::format("%6s %20s %6s %12s %6s %6s %12.6f %12.6f %8.3f %8s\n")
        % (idx+1)
        % snp->ID
        % snp->chrom
        % snp->physPos
        % (snp->flipped ? snp->a2 : snp->a1)
        % (snp->flipped ? snp->a1 : snp->a2)
        % (snp->flipped ? 1.0-snp->af : snp->af)
        % (snp->flipped ? -effectMean[idx] : effectMean[idx])
        % pip[idx]
        % snp->window;
        ++idx;
    }
    out.close();
}

void Data::outputFixedEffects(const MatrixXf &fixedEffects, const string &filename) const {
    ofstream out(filename.c_str());
    long nrow = fixedEffects.rows();
    VectorXf mean = fixedEffects.colwise().mean();
    VectorXf sd = (fixedEffects.rowwise() - mean.transpose()).colwise().squaredNorm().cwiseSqrt()/sqrt(nrow);
    for (unsigned i=0; i<numFixedEffects; ++i) {
        out << boost::format("%20s %12.6f %12.6f\n") % fixedEffectNames[i] %mean[i] %sd[i];
    }
    out.close();
}

void Data::outputRandomEffects(const MatrixXf &randomEffects, const string &filename) const {
    ofstream out(filename.c_str());
    long nrow = randomEffects.rows();
    VectorXf mean = randomEffects.colwise().mean();
    VectorXf sd = (randomEffects.rowwise() - mean.transpose()).colwise().squaredNorm().cwiseSqrt()/sqrt(nrow);
    for (unsigned i=0; i<numRandomEffects; ++i) {
        out << boost::format("%20s %12.6f %12.6f\n") % randomEffectNames[i] %mean[i] %sd[i];
    }
    out.close();
}

void Data::outputWindowResults(const VectorXf &posteriorMean, const string &filename) const {
    ofstream out(filename.c_str());
    out << boost::format("%6s %8s\n") %"Id" %"PIP";
    for (unsigned i=0; i<posteriorMean.size(); ++i) {
        out << boost::format("%6s %8.3f\n")
        % (i+1)
        % posteriorMean[i];
    }
    out.close();
}

void Data::readGwasSummaryFile(const string &gwasFile, const float afDiff, const float mafmin, const float mafmax, const float pValueThreshold, const bool imputeN, const bool removeOutlierN){
    ifstream in(gwasFile.c_str());
    if (!in) throw ("Error: can not open the GWAS summary data file [" + gwasFile + "] to read.");
    cout << "Reading GWAS summary data from [" + gwasFile + "]." << endl;
    
    string header;
    getline(in, header);

    SnpInfo *snp;
    map<string, SnpInfo*>::iterator it;
    string id, allele1, allele2, freq, b, se, pval, n;
    unsigned line=0, match=0;
    unsigned numInconAllele=0, numInconAf=0, numFixed=0, numMafMin=0, numMafMax=0, numOutlierN=0;
    unsigned numPvalPruned=0;
    unsigned numFlip=0;
    bool inconAllele, inconAf, fixed, ismafmin, ismafmax, isPvalPruned;
    float gwas_af;
    while (in >> id >> allele1 >> allele2 >> freq >> b >> se >> pval >> n) {
        ++line;
        it = snpInfoMap.find(id);
        if (it == snpInfoMap.end()) {
            //cout << "cannot find SNP " << id << endl;
            continue;
        }
        snp = it->second;
        if (!snp->included) {
            //cout << "exclude SNP " << id << endl;
            continue;
        }
        inconAllele = inconAf = fixed = ismafmin = ismafmax = isPvalPruned = false;
        if (allele1 == snp->a1 && allele2 == snp->a2) {
            gwas_af = atof(freq.c_str());
            snp->gwas_b  = atof(b.c_str());
            snp->gwas_af = gwas_af != -1 ? gwas_af : snp->af;  // set -1 in gwas summary file if allele frequencies are not available
            snp->gwas_se = atof(se.c_str());
            snp->gwas_n  = atof(n.c_str());
            snp->gwas_pvalue = atof(pval.c_str());
        } else if (allele1 == snp->a2 && allele2 == snp->a1) {
            gwas_af = atof(freq.c_str());
            snp->gwas_b  = -atof(b.c_str());
            snp->gwas_af = gwas_af != -1 ? 1.0 - gwas_af : 1.0 - snp->af;
            snp->gwas_se = atof(se.c_str());
            snp->gwas_n  = atof(n.c_str());
            snp->gwas_pvalue = atof(pval.c_str());
            snp->flipped = true;
//            snp->included = false;
            //cout << snp->index << " " << snp->ID << " " << snp->gwas_af << " " << snp->gwas_b << endl;
            ++numFlip;
        } else {
//            cout << "WARNING: SNP " + id + " has inconsistent allele coding in between the reference and GWAS samples." << endl;
            inconAllele = true;
            ++numInconAllele;
        }
        if (!inconAllele) {
            if (abs(snp->af - snp->gwas_af) > afDiff) {
                inconAf = true;
                ++numInconAf;
            } else if (snp->gwas_af==0 || snp->gwas_af==1) {
                fixed = true;
                ++numFixed;
            } else if (mafmin || mafmax) {
                float maf_ref = snp->af < 0.5 ? snp->af : 1.0 - snp->af;
                float maf_gwas = snp->gwas_af < 0.5 ? snp->gwas_af : 1.0 - snp->gwas_af;
                if (mafmin && (maf_ref < mafmin || maf_gwas < mafmin)) {
                    ismafmin = true;
                    ++numMafMin;
                }
                if (mafmax && (maf_ref > mafmax || maf_gwas > mafmax)) {
                    ismafmax = true;
                    ++numMafMax;
                }
            }
            if (snp->gwas_pvalue > pValueThreshold) {
                isPvalPruned = true;
                ++numPvalPruned;
            }
        }
        if (inconAllele || inconAf || fixed || ismafmin || ismafmax || isPvalPruned) {
            snp->included = false;
            cout << snp->index << " " << snp->ID << endl;
        } else ++match;
    }
    in.close();
        
    if (removeOutlierN){
        unsigned size = 0;
        for (unsigned i=0; i<numSnps; ++i) {
            snp = snpInfoVec[i];
            if (snp->included && snp->gwas_n != -999) ++size;
        }
        ArrayXf perSnpN(size);
        vector<SnpInfo*> snpvec(size);
        for (unsigned i=0, j=0; i<numSnps; ++i) {
            snp = snpInfoVec[i];
            if (snp->included && snp->gwas_n != -999) {
                perSnpN[j] = snp->gwas_n;
                snpvec[j] = snp;
                ++j;
            }
        }
        
        float n_med = Gadget::findMedian(perSnpN);
        float sd = sqrt(Gadget::calcVariance(perSnpN));
        for (unsigned i=0; i<size; ++i) {
            snp = snpvec[i];
            if (perSnpN[i] < n_med - 3*sd || perSnpN[i] > n_med + 3*sd) {
                snp->included = false;
                ++numOutlierN;
            }
        }
        
        match -= numOutlierN;
    }
    
    numIncdSnps = 0;
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        if (!snp->included) continue;
        if (snp->gwas_b == -999) {
            snp->included = false;
        } else {
            ++numIncdSnps;
        }
    }

    if (numFlip) cout << "flipped " << numFlip << " SNPs according to the minor allele in the reference and GWAS samples." << endl;
    if (numInconAllele) cout << "removed " << numInconAllele << " SNPs with inconsistent allele coding in between the reference and GWAS samples." << endl;
    if (numInconAf) cout << "removed " << numInconAf << " SNPs with differences in allele frequency between the reference and GWAS samples > " << afDiff << "." << endl;
    if (numFixed) cout << "removed " << numFixed << " fixed SNPs in the GWAS samples." << endl;
    if (mafmin) cout << "removed " << numMafMin << " SNPs with MAF below " << mafmin << " in either reference and GWAS samples." << endl;
    if (mafmax) cout << "removed " << numMafMax << " SNPs with MAF above " << mafmax << " in either reference and GWAS samples." << endl;
    if (pValueThreshold < 1.0) cout << "removed " << numPvalPruned << " SNPs with GWAS P value greater than " << pValueThreshold << "." << endl;
    if (numOutlierN) cout << "removed " << numOutlierN << " SNPs with per-SNP sample size beyond 3 SD around the median value." << endl;
    cout << match << " matched SNPs in the GWAS summary data (in total " << line << " SNPs)." << endl;
    
    if (imputeN) imputePerSnpSampleSize(snpInfoVec, numIncdSnps, 0);
}

void Data::imputePerSnpSampleSize(vector<SnpInfo*> &snpInfoVec, unsigned &numIncdSnps, float sd) {
    // use input allele frequencies, b_hat and se to impute per-snp N
    // then filter SNPs with N > 3 sd apart from the median value
    ArrayXf n(numIncdSnps);
    ArrayXf p(numIncdSnps);
    ArrayXf bsq(numIncdSnps);
    ArrayXf var(numIncdSnps);
    ArrayXf tpq(numIncdSnps);
    ArrayXf ypy(numIncdSnps);
    SnpInfo *snp;
    unsigned j = 0;
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        if (!snp->included) continue;
        n[j] = snp->gwas_n;
        p[j] = snp->gwas_af;
        bsq[j] = snp->gwas_b*snp->gwas_b;
        var[j] = snp->gwas_se*snp->gwas_se;
        ++j;
    }
    tpq = 2.0*p*(1.0-p);
    ypy = tpq*n.square()*var + tpq*n*bsq;
    float ypy_med = Gadget::findMedian(ypy);
    // Given ypy and n compute 2pq
    //tpq = ypy / (var*n.square() + bsq*n);
    // Given ypy_med and 2pq compute n
    float n_med = Gadget::findMedian(n);
    float vary = ypy_med / n_med;
    n = (vary - tpq*bsq) / (tpq*var);
    // compute sd of n
    float sdOld = sd;
    sd = sqrt(Gadget::calcVariance(n));
    float p_new;
    j = 0;
    numIncdSnps = 0;
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        if (!snp->included) continue;
        if (n[j] < n_med - 3*sd || n[j] > n_med + 3*sd) {
            snp->included = false;
        }
        else {
            snp->gwas_n = n[j];
            //p_new = 0.5 - 0.5*sqrt(1.0-2.0*tpq[j]);
            //snp->gwas_af = p[j] < 0.5 ? p_new : 1.0-p_new;
            ++numIncdSnps;
        }
        ++j;
    }
//    cout << n.mean() << " " << ypy_med << " " << sdOld << " " << sd << " " << numIncdSnps << " " << n.head(10).transpose() << endl;
    if (abs(sd-sdOld) > 0.01) {
        imputePerSnpSampleSize(snpInfoVec, numIncdSnps, sd);
    } else {
        cout << numIncdSnps << " SNPs with per-SNP sample size within 3 sd around the median value of " << n_med << endl;
        string outfile = title + ".imputedPerSnpN";
        ofstream out(outfile.c_str());
        out << boost::format("%15s %12s\n")
        % "ID" % "Imputed_N";
        for (unsigned i=0; i<numSnps; ++i) {
            snp = snpInfoVec[i];
            if (!snp->included) continue;
            out << boost::format("%15s %12s\n")
            % snp->ID
            % snp->gwas_n;
        }
        out.close();
        return;
    }
}

/*
 * Divide the big matrix to small part (From Zhili)
 * @Param totalPart: number of parts would like to divide 
 * @Param curPart:  current part would like to run. 1 based
 * @Param nSNPs:  total number of SNPs
 * @Return string: SNP idx start-end, 1 based
*/
string makeSNPRangeString(int totalPart, int curPart, int nSNPs){
    int nPartSNP = (nSNPs + totalPart - 1) / totalPart;
    int start = nPartSNP * (curPart - 1) + 1;
    int end = nPartSNP * curPart;
    if(end > nSNPs) end = nSNPs;
    return(to_string(start) + "-" + to_string(end));
}

string Data::partLDMatrix(const string &partParam, const string &outfilename, const string &LDmatType){
    Gadget::Tokenizer token;
    token.getTokens(partParam, ",");
    string snpRange = "";
    if(token.size() == 2){
        int nTotalPart = atoi(token[0].c_str());
        int nCurPart = atoi(token[1].c_str());
        cout << "Dividing LD matrix by " << nTotalPart << " parts, current running part " << nCurPart << std::endl;
        if(nTotalPart < nCurPart || nCurPart <= 0){
            throw("--part usage total,curentPart"); 
        }

        snpRange = makeSNPRangeString(nTotalPart, nCurPart, numIncdSnps);
        cout << "  SNP range " << snpRange << endl;

        if(nCurPart == 1){
            ofstream mldmfile((outfilename + ".mldm").c_str());
            for(int i = 1; i <= nTotalPart; i++){
                string snpRange1 = makeSNPRangeString(nTotalPart, i, numIncdSnps);
                string outfilename2 = outfilename + ".snp" + snpRange1 + ".ldm." + LDmatType;
                mldmfile << outfilename2 << endl;
            }
            mldmfile.close();
            cout << " run gctb --mldm " << outfilename << ".mldm --make-full-ldm --out ...  to combine the parted LDM matrix" << std::endl;  
        }
    }
    return snpRange;
}

void Data::makeLDmatrix(const string &bedFile, const string &LDmatType, const float chisqThreshold, const float LDthreshold, const unsigned windowWidth, const string &snpRange, const string &filename, const bool writeLdmTxt){
    
    Gadget::Tokenizer token;
    token.getTokens(snpRange, "-");
    
    unsigned start = 0;
    unsigned end = numIncdSnps;
    
    if (token.size()) {
        start = atoi(token[0].c_str()) - 1;
        end = atoi(token[1].c_str());
        if (end > numIncdSnps) end = numIncdSnps;
    }
    
    unsigned numSnpInRange = end - start;
    
    if (snpRange.empty())
        cout << "Building " + LDmatType + " LD matrix for all SNPs ..." << endl;
    else
        cout << "Building " + LDmatType + " LD matrix for SNPs " << snpRange << " ..." << endl;
    
    if (numIncdSnps == 0) throw ("Error: No SNP is retained for analysis.");
    if (numKeptInds == 0) throw ("Error: No individual is retained for analysis.");
    if (start >= numIncdSnps) throw ("Error: Specified a SNP range of " + snpRange + " but " + to_string(static_cast<long long>(numIncdSnps)) + " SNPs are included.");
    
    Gadget::Timer timer;
    timer.setTime();
    
    unsigned firstWindStart = 0;
    unsigned lastWindEnd = 0;
    if (windowWidth) {
        getWindowInfo(incdSnpInfoVec, windowWidth, windStart, windSize);
    }
    
    // first read in the genotypes of SNPs in the given range
    
    const int bedToGeno[4] = {2, -9, 1, 0};
    unsigned size = (numInds+3)>>2;
    
    MatrixXf ZP(numSnpInRange, numKeptInds);  // SNP x Ind
    D.setZero(numSnpInRange);

    if (numKeptInds < 2) throw("Error: Cannot calculate LD matrix with number of individuals < 2.");
    
    FILE *in1 = fopen(bedFile.c_str(), "rb");
    if (!in1) throw ("Error: can not open the file [" + bedFile + "] to read.");
    cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
    char header[3];
    fread(header, sizeof(header), 1, in1);
    if (!in1 || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
        cerr << "Error: Incorrect first three bytes of bed file: " << bedFile << endl;
        exit(1);
    }

    
    IndInfo *indi = NULL;
    SnpInfo *snpj = NULL;
    SnpInfo *snpk = NULL;
    
    int genoValue;
    unsigned i, j, k;
    unsigned incj, inck; // index of included SNP
    unsigned long long skipj = 0;
    unsigned nmiss;
    float mean;
    
    set<int> chromInRange;
    
    for (j = 0, incj = 0; j < numSnps; j++) {
        snpj = snpInfoVec[j];
        
        if (snpj->index < start || !snpj->included) {
            skipj += size;
            continue;
        }
        
        if (skipj) fseek(in1, skipj, SEEK_CUR);
        skipj = 0;
        
        char *bedLineIn = new char[size];
        fread(bedLineIn, sizeof(char), size, in1);
        
        chromInRange.insert(snpj->chrom);
        
        mean = 0.0;
        nmiss = 0;
        
        for (i = 0; i < numInds; i++) {
            indi = indInfoVec[i];
            if (!indi->kept) continue;
            genoValue = bedToGeno[(bedLineIn[i>>2]>>((i&3)<<1))&3];
            ZP(incj, indi->index) = genoValue;
            if (genoValue == -9) ++nmiss;
            else mean += genoValue;
        }
        delete[] bedLineIn;
        
        // fill missing values with the mean
        snpj->sampleSize = numKeptInds-nmiss;
        mean /= float(snpj->sampleSize);
        if (nmiss) {
            for (i=0; i<numKeptInds; ++i) {
                if (ZP(incj, i) == -9) ZP(incj, i) = mean;
            }
        }
        
        // compute allele frequency
        snpj->af = 0.5f*mean;
        snp2pq[incj] = snpj->twopq = 2.0f*snpj->af*(1.0f-snpj->af);
        
        if (snp2pq[incj]==0) throw ("Error: " + snpj->ID + " is a fixed SNP (MAF=0)!");
        
        // standardize genotypes
        //D[incj] = snp2pq[incj]*snpj->sampleSize;
        D[incj] = Gadget::calcVariance(ZP.row(incj))*numKeptInds;
        
        ZP.row(incj) = (ZP.row(incj).array() - ZP.row(incj).mean())/sqrt(D[incj]);
//        ZP.row(incj) = (ZP.row(incj).array() - mean)/sqrt(D[incj]);
        
        if (windowWidth) {
            if (incj == 0) firstWindStart = snpj->windStart;
            if (incj == numSnpInRange-1) lastWindEnd = snpj->windStart + snpj->windSize;
        }
        
        if (++incj == numSnpInRange) break;
    }
    
    fclose(in1);
    
//    ZP = ZP.colwise() - ZP.rowwise().mean();
//    ZP = ZP.array().colwise() / D.cwiseSqrt().array();
    ZPZdiag = ZP.rowwise().squaredNorm();
    
//    cout << ZP.rowwise().mean() << endl << endl;
//    cout << ZP.block(0, 0, 10, 10) << endl;
    
    // then read in the bed file again to compute Z'Z
    

    MatrixXf denseZPZ;
    denseZPZ.setZero(numSnpInRange, numIncdSnps);
    VectorXf Zk(numKeptInds);
    D.setZero(numIncdSnps);
    
//    float numJackknife = numKeptInds-1;
//    MatrixXf ZPZkCwise;
//    MatrixXf samplVarEmp;
//    samplVarEmp.setZero(numSnpInRange, numIncdSnps);
    
    FILE *in2 = fopen(bedFile.c_str(), "rb");
    fseek(in2, 3, SEEK_SET);
    unsigned long long skipk = 0;
    
    set<int>::iterator setend = chromInRange.end();

    if (numSkeletonSnps) {
        for (k = 0, inck = 0; k < numSnps; k++) {
            snpk = snpInfoVec[k];
            
            if (!snpk->included) {
                skipk += size;
                continue;
            }

            if(!(inck%1000)) cout << " read snp " << inck << "\r" << flush;

            if (chromInRange.find(snpk->chrom) == setend && !snpk->skeleton) {
                skipk += size;
                ++inck;       // ensure the index is correct
                continue;
            }
            
            if (windowWidth) {
                if (inck < firstWindStart) {
                    skipk += size;
                    continue;
                } else if (inck > lastWindEnd) {
                    break;
                }
            }
            
            if (skipk) fseek(in2, skipk, SEEK_CUR);
            skipk = 0;
            
            char *bedLineIn = new char[size];
            fread(bedLineIn, sizeof(char), size, in2);
            
            mean = 0.0;
            nmiss = 0;
            
            for (i = 0; i < numInds; i++) {
                indi = indInfoVec[i];
                if (!indi->kept) continue;
                genoValue = bedToGeno[(bedLineIn[i>>2]>>((i&3)<<1))&3];
                Zk[indi->index] = genoValue;
                if (genoValue == -9) ++nmiss;   // missing genotype
                else mean += genoValue;
            }
            delete[] bedLineIn;
            
            // fill missing values with the mean
            snpk->sampleSize = numKeptInds-nmiss;
            mean /= float(snpk->sampleSize);
            if (nmiss) {
                for (i=0; i<numKeptInds; ++i) {
                    if (Zk[i] == -9) Zk[i] = mean;
                }
            }
            
            // compute allele frequency
            snpk->af = 0.5f*mean;
            snp2pq[inck] = snpk->twopq = 2.0f*snpk->af*(1.0f-snpk->af);
            
            if (snp2pq[inck]==0) throw ("Error: " + snpk->ID + " is a fixed SNP (MAF=0)!");
            
            // standardize genotypes
            //D[inck] = snp2pq[inck]*snpk->sampleSize;
            D[inck] = Gadget::calcVariance(Zk.row(inck))*numKeptInds;

            Zk = (Zk.array() - Zk.mean())/sqrt(D[inck]);
//            Zk = (Zk.array() - mean)/sqrt(D[inck]);
            
            denseZPZ.col(inck) = ZP * Zk;

//            cout << " inck " << inck << " snpk " << k << " chr " << snpk->chrom << " " << ZP*Zk << endl;
            
            ++inck;
        }
    }
    else {
        for (k = 0, inck = 0; k < numSnps; k++) {
            snpk = snpInfoVec[k];
            
            if (!snpk->included) {
                skipk += size;
                continue;
            }
            
            if (windowWidth) {
                if (inck < firstWindStart) {
                    skipk += size;
                    continue;
                } else if (inck > lastWindEnd) {
                    break;
                }
            }
            
            if (skipk) fseek(in2, skipk, SEEK_CUR);
            skipk = 0;
            
            char *bedLineIn = new char[size];
            fread(bedLineIn, sizeof(char), size, in2);
            
            mean = 0.0;
            nmiss = 0;
            
            for (i = 0; i < numInds; i++) {
                indi = indInfoVec[i];
                if (!indi->kept) continue;
                genoValue = bedToGeno[(bedLineIn[i>>2]>>((i&3)<<1))&3];
                Zk[indi->index] = genoValue;
                if (genoValue == -9) ++nmiss;   // missing genotype
                else mean += genoValue;
            }
            delete[] bedLineIn;
            
            // fill missing values with the mean
            snpk->sampleSize = numKeptInds-nmiss;
            mean /= float(snpk->sampleSize);
            if (nmiss) {
                for (i=0; i<numKeptInds; ++i) {
                    if (Zk[i] == -9) Zk[i] = mean;
                }
            }
            
            // compute allele frequency
            snpk->af = 0.5f*mean;
            snp2pq[inck] = snpk->twopq = 2.0f*snpk->af*(1.0f-snpk->af);
            
            if (snp2pq[inck]==0) throw ("Error: " + snpk->ID + " is a fixed SNP (MAF=0)!");
            
            // standardize genotypes
            //D[inck] = snp2pq[inck]*snpk->sampleSize;
            D[inck] = Gadget::calcVariance(Zk)*numKeptInds;

            Zk = (Zk.array() - Zk.mean())/sqrt(D[inck]);
//            Zk = (Zk.array() - mean)/sqrt(D[inck]);
            
            denseZPZ.col(inck) = ZP * Zk;
            
//            // Jackknife estimate of correlation and sampling variance
//            ZPZkCwise = ZP.array().rowwise() * Zk.transpose().array();
//            denseZPZ.col(inck) = ZPZkCwise.rowwise().sum();
//            ZPZkCwise = - (ZPZkCwise.colwise() - denseZPZ.col(inck));
//            ZPZkCwise *= numKeptInds/numJackknife;
//            samplVarEmp.col(inck) = (ZPZkCwise.colwise() - ZPZkCwise.rowwise().mean()).rowwise().squaredNorm() * numJackknife/numKeptInds;
//            // Jackknife end
            
            if(!(inck%1000)) cout << " read snp " << inck << "\r" << flush;
            
            ++inck;
        }
    }

    fclose(in2);
    
    //cout << denseZPZ.block(0, 0, 10, 10) << endl;
    

    // find out per-SNP window position
    
    if (LDmatType == "full") {
        ZPZ.resize(numSnpInRange);
        windStart.setZero(numSnpInRange);
        windSize.setConstant(numSnpInRange, numIncdSnps);
        for (unsigned i=0; i<numSnpInRange; ++i) {
            SnpInfo *snp = incdSnpInfoVec[start+i];
            snp->windStart = 0;
            snp->windSize  = numIncdSnps;
            snp->windEnd   = numIncdSnps-1;
            ZPZ[i] = denseZPZ.row(i);
            snp->ldSamplVar = (1.0 - denseZPZ.row(i).array().square()).square().sum()/numKeptInds;
            snp->ldSum = denseZPZ.row(i).sum();
//            snp->ldSum = samplVarEmp.row(i).sum();
       }
    }
    else if (LDmatType == "band") {
        ZPZ.resize(numSnpInRange);
        if (windowWidth) {  // based on the given window width
            for (unsigned i=0; i<numSnpInRange; ++i) {
                SnpInfo *snp = incdSnpInfoVec[start+i];
                ZPZ[i] = denseZPZ.row(i).segment(snp->windStart, snp->windSize);
                snp->ldSamplVar = (1.0 - ZPZ[i].array().square()).square().sum()/numKeptInds;
                snp->ldSum = ZPZ[i].sum();
            }
        } else {  // based on the given LD threshold
            windStart.setZero(numSnpInRange);
            windSize.setZero(numSnpInRange);
            for (unsigned i=0; i<numSnpInRange; ++i) {
                SnpInfo *snp = incdSnpInfoVec[start+i];
                unsigned windEndi = numIncdSnps;
                for (unsigned j=0; j<numIncdSnps; ++j) {
                    if (abs(denseZPZ(i,j)) > LDthreshold) {
                        windStart[i] = snp->windStart = j;
                        break;
                    }
                }
                for (unsigned j=numIncdSnps; j>0; --j) {
                    if (abs(denseZPZ(i,j-1)) > LDthreshold) {
                        windEndi = j;
                        break;
                    }
                }
                windSize[i] = snp->windSize = windEndi - windStart[i];
                snp->windEnd = windEndi - 1;
                ZPZ[i].resize(windSize[i]);
                VectorXf::Map(&ZPZ[i][0], windSize[i]) = denseZPZ.row(i).segment(windStart[i], windSize[i]);
                snp->ldSamplVar = (1.0 - ZPZ[i].array().square()).square().sum()/numKeptInds;
                snp->ldSum = ZPZ[i].sum();
            }
        }
    }
    else if (LDmatType == "sparse") {
        ZPZsp.resize(numSnpInRange);
        windStart.setZero(numSnpInRange);
        windSize.setZero(numSnpInRange);
        float rsq = 0.0;
        SnpInfo *snpi, *snpj;
        if (numSkeletonSnps) {
            if (LDthreshold) {
                for (unsigned i=0; i<numSnpInRange; ++i) {
                    snpi = incdSnpInfoVec[start+i];
                    snpi->ldSamplVar = 0.0;
                    snpi->ldSum = 0.0;
                    for (unsigned j=0; j<numIncdSnps; ++j) {
                        snpj = incdSnpInfoVec[j];
                        if (snpj->skeleton) {
                            rsq = denseZPZ(i,j)*denseZPZ(i,j);
                            snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/numKeptInds;
                            snpi->ldSum += denseZPZ(i,j);
                        }
                        else {
                            if (abs(denseZPZ(i,j)) < LDthreshold) denseZPZ(i,j) = 0;
                            else {
                                rsq = denseZPZ(i,j)*denseZPZ(i,j);
                                snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/numKeptInds;
                                snpi->ldSum += denseZPZ(i,j);
                            }
                        }
                    }
                    ZPZsp[i] = denseZPZ.row(i).sparseView();
                    SparseVector<float>::InnerIterator it(ZPZsp[i]);
                    windStart[i] = snpi->windStart = it.index();
                    windSize[i] = snpi->windSize = ZPZsp[i].nonZeros();
                    for (; it; ++it) snpi->windEnd = it.index();
                }
            } else {
                for (unsigned i=0; i<numSnpInRange; ++i) {
                    snpi = incdSnpInfoVec[start+i];
                    snpi->ldSamplVar = 0.0;
                    snpi->ldSum = 0.0;
//                    unsigned cnt = 0;
                    for (unsigned j=0; j<numIncdSnps; ++j) {
                        snpj = incdSnpInfoVec[j];
                        if (snpj->skeleton) {
                            rsq = denseZPZ(i,j)*denseZPZ(i,j);
                            snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/numKeptInds;
                            snpi->ldSum += denseZPZ(i,j);
//                            cout << j << " " << denseZPZ(i,j) << endl;
//                            ++cnt;
                        }
                        else {
                            if (i!=j && denseZPZ(i,j)*denseZPZ(i,j)*snpi->sampleSize < chisqThreshold) denseZPZ(i,j) = 0;
                            else {
                                rsq = denseZPZ(i,j)*denseZPZ(i,j);
                                snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/numKeptInds;
                                snpi->ldSum += denseZPZ(i,j);
//                                cout << j << " " << denseZPZ(i,j) << endl;
//                                ++cnt;
                            }
                        }
                    }
//                    cout << "cnt " << cnt << endl;
                    ZPZsp[i] = denseZPZ.row(i).sparseView();
                    SparseVector<float>::InnerIterator it(ZPZsp[i]);
                    windStart[i] = snpi->windStart = it.index();
                    windSize[i] = snpi->windSize = ZPZsp[i].nonZeros();
                    for (; it; ++it) snpi->windEnd = it.index();
                }
            }
        }
        else {
            if (LDthreshold) {
                for (unsigned i=0; i<numSnpInRange; ++i) {
                    snpi = incdSnpInfoVec[start+i];
                    snpi->ldSamplVar = 0.0;
                    snpi->ldSum = 0.0;
                    for (unsigned j=0; j<numIncdSnps; ++j) {
                        snpj = incdSnpInfoVec[j];
                        if (abs(denseZPZ(i,j)) < LDthreshold) denseZPZ(i,j) = 0;
                        else {
                            rsq = denseZPZ(i,j)*denseZPZ(i,j);
                            snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/numKeptInds;
                            snpi->ldSum += denseZPZ(i,j);
                        }
                    }
                    ZPZsp[i] = denseZPZ.row(i).sparseView();
                    SparseVector<float>::InnerIterator it(ZPZsp[i]);
                    windStart[i] = snpi->windStart = it.index();
                    windSize[i] = snpi->windSize = ZPZsp[i].nonZeros();
                    for (; it; ++it) snpi->windEnd = it.index();
                }
            } else {
                for (unsigned i=0; i<numSnpInRange; ++i) {
                    snpi = incdSnpInfoVec[start+i];
                    snpi->ldSamplVar = 0.0;
                    snpi->ldSum = 0.0;
                    for (unsigned j=0; j<numIncdSnps; ++j) {
                        snpj = incdSnpInfoVec[j];
                        if (i!=j && denseZPZ(i,j)*denseZPZ(i,j)*numKeptInds < chisqThreshold) denseZPZ(i,j) = 0;
                        else {
                            rsq = denseZPZ(i,j)*denseZPZ(i,j);
                            snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/numKeptInds;
                            snpi->ldSum += denseZPZ(i,j);
                        }
                    }
                    ZPZsp[i] = denseZPZ.row(i).sparseView();
                    SparseVector<float>::InnerIterator it(ZPZsp[i]);
                    windStart[i] = snpi->windStart = it.index();
                    windSize[i] = snpi->windSize = ZPZsp[i].nonZeros();
                    for (; it; ++it) snpi->windEnd = it.index();
                }
            }
        }
    }
    
    denseZPZ.resize(0,0);
    
    //    cout << denseZPZ.block(0,0,10,10) << endl;
    //    cout << windStart.transpose() << endl;
    //    cout << windSize.transpose() << endl;
    
    
    timer.getTime();
    
    cout << endl;
    displayAverageWindowSize(windSize);
    cout << "LD matrix diagonal mean " << ZPZdiag.mean() << " variance " << Gadget::calcVariance(ZPZdiag) << "." << endl;
    if (ZPZdiag.mean() < 0.8 || ZPZdiag.mean() > 1.2) cout << "ERROR: The mean of LD matrix diagonal values is expected to be close to one. Something is wrong with the LD matrix!" << endl;
    cout << "Genotype data for " << numKeptInds << " individuals and " << numSnpInRange << " SNPs are included from [" + bedFile + "]." << endl;
    cout << "Build of LD matrix completed (time used: " << timer.format(timer.getElapse()) << ")." << endl;
        
    vector<SnpInfo*> snpVecTmp(numSnpInRange);
    for (unsigned i=0; i<numSnpInRange; ++i) {
        snpVecTmp[i] = incdSnpInfoVec[start+i];
    }
    incdSnpInfoVec = snpVecTmp;
    numIncdSnps = numSnpInRange;
    string outfilename = filename;
    if (!snpRange.empty()) outfilename += ".snp" + snpRange;
    outputLDmatrix(LDmatType, outfilename, writeLdmTxt);
}

void Data::outputLDmatrix(const string &LDmatType, const string &filename, const bool writeLdmTxt) const {
    string outfilename = filename + ".ldm." + LDmatType;
    string outfile1 = outfilename + ".info";
    string outfile2 = outfilename + ".bin";
    ofstream out1(outfile1.c_str());
    FILE *out2 = fopen(outfile2.c_str(), "wb");
    ofstream out3;
    string outfile3;
    if (writeLdmTxt) {
        outfile3 = outfilename + ".txt";
        out3.open(outfile3.c_str());
    }
    out1 << boost::format("%6s %15s %10s %15s %6s %6s %12s %10s %10s %10s %10s %15s %10s %12s %12s\n")
    % "Chrom"
    % "ID"
    % "GenPos"
    % "PhysPos"
    % "A1"
    % "A2"
    % "A1Freq"
    % "Index"
    % "WindStart"
    % "WindEnd"
    % "WindSize"
    % "WindWidth"
    % "N"
    % "SamplVar"
    % "LDsum";
    SnpInfo *snp, *windStart, *windEnd;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snp = incdSnpInfoVec[i];
        windStart = incdSnpInfoVec[snp->windStart];
        windEnd = incdSnpInfoVec[snp->windEnd];
        out1 << boost::format("%6s %15s %10s %15s %6s %6s %12f %10s %10s %10s %10s %15s %10s %12.6f %12.6f\n")
        % snp->chrom
        % snp->ID
        % snp->genPos
        % snp->physPos
        % snp->a1
        % snp->a2
        % snp->af
        % snp->index
        % snp->windStart
        % snp->windEnd
        % snp->windSize
        % (windStart->chrom == windEnd->chrom ? windEnd->physPos - windStart->physPos : windStart->chrom-windEnd->chrom)
        % numKeptInds
        % snp->ldSamplVar
        % snp->ldSum;
        if (LDmatType == "sparse") {
            fwrite(ZPZsp[i].innerIndexPtr(), sizeof(unsigned), ZPZsp[i].nonZeros(), out2);
            fwrite(ZPZsp[i].valuePtr(), sizeof(float), ZPZsp[i].nonZeros(), out2);
            if (writeLdmTxt) {
                //out3 << ZPZsp[i].transpose();
                for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                    out3 << boost::format("%-8s %-15s %-8s %-15s %-15s\n")
                    % (i+1)
                    % snp->ID
                    % (it.index()+1)
                    % incdSnpInfoVec[it.index()]->ID
                    % it.value();
                }
            }
        } else {
            fwrite(&ZPZ[i][0], sizeof(float), snp->windSize, out2);
            if (writeLdmTxt) out3 << ZPZ[i].transpose() << endl;
        }
    }
    out1.close();
    fclose(out2);
    
    cout << "Written the SNP info into file [" << outfile1 << "]." << endl;
    cout << "Written the LD matrix into file [" << outfile2 << "]." << endl;
    
    if (writeLdmTxt) {
        out3.close();
        cout << "Written the LD matrix into text file [" << outfile3 << "]." << endl;
    }
}


void Data::displayAverageWindowSize(const VectorXi &windSize){
//    float windSizeMean = 0.0;
//    float windSizeSqMean = 0.0;
//    long size = windSize.size();
//    for (unsigned i=0; i<size; ++i) {
//        windSizeMean += (windSize[i] - windSizeMean)/(i+1);
//        windSizeSqMean += (windSize[i]*windSize[i] - windSizeSqMean)/(i+1);
//    }
//    cout << "Per-SNP window size mean " << windSizeMean << " sd " << windSizeSqMean-windSizeMean*windSizeMean << "." << endl;
    VectorXf windSizeFloat = windSize.cast<float>();
    float windSizeMean = Gadget::calcMean(windSizeFloat);
    float windSizeSD = sqrt(Gadget::calcVariance(windSizeFloat));
    cout << "Per-SNP window size mean " << windSizeMean << " sd " << windSizeSD << "." << endl;
}

void Data::resizeWindow(const vector<SnpInfo *> &incdSnpInfoVec, const VectorXi &windStartOri, const VectorXi &windSizeOri,
                        VectorXi &windStart, VectorXi &windSize){
    bool reindexed = false;
    for (unsigned i=0; i<numSnps; ++i) {
        if (!snpInfoVec[i]->included) {
            reindexed = true;
            break;
        }
    }
    if (reindexed == false) {
        windStart = windStartOri;
        windSize  = windSizeOri;
        return;
    }
    cout << "Resizing per-SNP LD window..." << endl;
    windStart.setZero(numIncdSnps);
    windSize.setZero(numIncdSnps);
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snpi = snpInfoVec[i];
        if (!snpi->included) continue;
        unsigned windEndOri = windStartOri[i] + windSizeOri[i];
        for (unsigned j=windStartOri[i]; j<windEndOri; ++j) {
            SnpInfo *snpj = snpInfoVec[j];
            if (!snpj->included) continue;
            if (!windSize[snpi->index]) {
                windStart[snpi->index] = snpj->index;
            }
            ++windSize[snpi->index];
        }
        snpi->windStart = windStart[snpi->index];
        snpi->windSize  = windSize[snpi->index];
        snpi->windEnd   = snpi->windStart + snpi->windSize - 1;
    }
}

void Data::readLDmatrixInfoFileOld(const string &ldmatrixFile){   // old format: no allele frequency, no header
    ifstream in(ldmatrixFile.c_str());
    if (!in) throw ("Error: can not open the file [" + ldmatrixFile + "] to read.");
    //cout << "Reading SNP info from [" + ldmatrixFile + "]." << endl;
    //snpInfoVec.clear();
    //snpInfoMap.clear();
    string header;
    string id, allele1, allele2;
    unsigned chr, physPos;
    float genPos;
    unsigned idx, windStart, windEnd, windSize, windWidth;
    long sampleSize;
    while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2 >> idx >> windStart >> windEnd >> windSize >> windWidth >> sampleSize) {
        SnpInfo *snp = new SnpInfo(idx, id, allele1, allele2, chr, genPos, physPos);
        snp->windStart = snp->windStartOri = windStart;
        snp->windEnd = snp->windEndOri = windEnd;
        snp->windSize = snp->windSizeOri = windSize;
        snp->sampleSize = sampleSize;
        snpInfoVec.push_back(snp);
        if (snpInfoMap.insert(pair<string, SnpInfo*>(id, snp)).second == false) {
            throw ("Error: Duplicate SNP ID found: \"" + id + "\".");
        }
    }
    in.close();
    numSnps = (unsigned) snpInfoVec.size();
    cout << numSnps << " SNPs to be included from [" + ldmatrixFile + "]." << endl;
}

void Data::readLDmatrixInfoFile(const string &ldmatrixFile){
    ifstream in(ldmatrixFile.c_str());
    if (!in) throw ("Error: can not open the file [" + ldmatrixFile + "] to read.");
    cout << "Reading SNP info from [" + ldmatrixFile + "]." << endl;
    //snpInfoVec.clear();
    //snpInfoMap.clear();
    string header;
    string id, allele1, allele2;
    unsigned chr, physPos;
    float genPos, af, ldSamplVar, ldSum;
    unsigned idx, windStart, windEnd, windSize;
    int windWidth;
    long sampleSize;
    bool skeleton;
    getline(in, header);
    Gadget::Tokenizer token;
    token.getTokens(header, " ");
    
    if (token.size() == 12) {
        in.close();
        readLDmatrixInfoFileOld(ldmatrixFile);
        return;
    }

    if (token.back() == "Skeleton") {
        while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2 >> af >> idx >> windStart >> windEnd >> windSize >> windWidth >> sampleSize >> ldSamplVar >> ldSum >> skeleton) {
            SnpInfo *snp = new SnpInfo(idx, id, allele1, allele2, chr, genPos, physPos);
            snp->af = af;
            snp->twopq = 2.0*af*(1.0-af);
            snp->windStart = snp->windStartOri = windStart;
            snp->windEnd = snp->windEndOri = windEnd;
            snp->windSize = snp->windSizeOri = windSize;
            snp->sampleSize = sampleSize;
            snp->ldSamplVar = ldSamplVar;
            snp->ldSum = ldSum;
            snp->skeleton = skeleton;
            snpInfoVec.push_back(snp);
            if (snpInfoMap.insert(pair<string, SnpInfo*>(id, snp)).second == false) {
                throw ("Error: Duplicate SNP ID found: \"" + id + "\".");
            }
        }
    }
    else if (token.back() == "LDsum") {
        while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2 >> af >> idx >> windStart >> windEnd >> windSize >> windWidth >> sampleSize >> ldSamplVar >> ldSum) {
            SnpInfo *snp = new SnpInfo(idx, id, allele1, allele2, chr, genPos, physPos);
            snp->af = af;
            snp->twopq = 2.0*af*(1.0-af);
            snp->windStart = snp->windStartOri = windStart;
            snp->windEnd = snp->windEndOri = windEnd;
            snp->windSize = snp->windSizeOri = windSize;
            snp->sampleSize = sampleSize;
            snp->ldSamplVar = ldSamplVar;
            snp->ldSum = ldSum;
            snpInfoVec.push_back(snp);
            if (snpInfoMap.insert(pair<string, SnpInfo*>(id, snp)).second == false) {
                throw ("Error: Duplicate SNP ID found: \"" + id + "\".");
            }
        }
    }
    else {
        while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2 >> af >> idx >> windStart >> windEnd >> windSize >> windWidth >> sampleSize >> ldSamplVar) {
            SnpInfo *snp = new SnpInfo(idx, id, allele1, allele2, chr, genPos, physPos);
            snp->af = af;
            snp->twopq = 2.0*af*(1.0-af);
            snp->windStart = snp->windStartOri = windStart;
            snp->windEnd = snp->windEndOri = windEnd;
            snp->windSize = snp->windSizeOri = windSize;
            snp->sampleSize = sampleSize;
            snp->ldSamplVar = ldSamplVar;
            snpInfoVec.push_back(snp);
            if (snpInfoMap.insert(pair<string, SnpInfo*>(id, snp)).second == false) {
                throw ("Error: Duplicate SNP ID found: \"" + id + "\".");
            }
        }
    }
    in.close();
    numSnps = (unsigned) snpInfoVec.size();
    cout << numSnps << " SNPs to be included from [" + ldmatrixFile + "]." << endl;
}

void Data::readLDmatrixBinFile(const string &ldmatrixFile){

    Gadget::Timer timer;
    timer.setTime();
    
    Gadget::Tokenizer token;
    token.getTokens(ldmatrixFile, ".");
    string ldmType = token[token.size()-2];
    sparseLDM = ldmType == "sparse" ? true : false;

    // New part for shrunk part ldm 
    shrunkLDM = ldmType == "shrunk" ? true : false;
    
    VectorXi windStartLDM(numSnps);
    VectorXi windSizeLDM(numSnps);
    
    windStart.resize(numIncdSnps);
    windSize.resize(numIncdSnps);
    
    SnpInfo *snpi, *snpj;
    
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snpi = snpInfoVec[i];
        windStartLDM[i] = snpi->windStart;
        windSizeLDM[i]  = snpi->windSize;
    }

    FILE *in = fopen(ldmatrixFile.c_str(), "rb");
    if (!in) {
        throw("Error: cannot open LD matrix file " + ldmatrixFile);
    }
    
    if (!sparseLDM) resizeWindow(incdSnpInfoVec, windStartLDM, windSizeLDM, windStart, windSize);
    
    if (numIncdSnps == 0) throw ("Error: No SNP is retained for analysis.");
    
    cout << "Reading " + ldmType + " LD matrix from [" + ldmatrixFile + "]..." << endl;
    
    float rsq = 0.0;
    
    if (sparseLDM) {
        ZPZsp.resize(numIncdSnps);
        ZPZdiag.resize(numIncdSnps);
       
        for (unsigned i=0, inci=0; i<numSnps; i++) {
            snpi = snpInfoVec[i];
                        
            unsigned d[windSizeLDM[i]];
            float v[windSizeLDM[i]];
            
            if (!snpi->included) {
                fseek(in, sizeof(d), SEEK_CUR);
                fseek(in, sizeof(v), SEEK_CUR);
                continue;
            }
            
            fread(d, sizeof(d), 1, in);
            fread(v, sizeof(v), 1, in);
            
            ZPZsp[inci].resize(windSizeLDM[i]);
            snpi->ldSamplVar = 0.0;
            snpi->ldSum = 0.0;
            if (!readLDscore) snpi->ldsc = 0.0;
            
            for (unsigned j=0; j<windSizeLDM[i]; ++j) {
                snpj = snpInfoVec[d[j]];
        
                if (snpj->included) {
                    ZPZsp[inci].insertBack(snpj->index) = v[j];
                    rsq = v[j]*v[j];
                    snpi->ldSamplVar += (1.0f-rsq)*(1.0f-rsq)/snpi->sampleSize;
                    snpi->ldSum += v[j];
                    if (!readLDscore) snpi->ldsc += rsq;
                    if (snpj == snpi)
                        ZPZdiag[inci] = v[j];
                }
            }
            SparseVector<float>::InnerIterator it(ZPZsp[inci]);
            windStart[inci] = snpi->windStart = it.index();
            windSize[inci] = snpi->windSize = ZPZsp[inci].nonZeros();
            for (; it; ++it) snpi->windEnd = it.index();
            snpi->numNonZeroLD = snpi->windSize;
            
            if (++inci == numIncdSnps) break;
        }
    }
    else {
        ZPZ.resize(numIncdSnps);
        ZPZdiag.resize(numIncdSnps);
        
        for (unsigned i=0, inci=0; i<numSnps; i++) {
            snpi = snpInfoVec[i];
            
            float v[windSizeLDM[i]];
            
            if (!snpi->included) {
                fseek(in, sizeof(v), SEEK_CUR);
                continue;
            }
            
            fread(v, sizeof(v), 1, in);
            
            ZPZ[inci].resize(windSize[inci]);
            snpi->ldSamplVar = 0.0;
            snpi->ldSum = 0.0;
            if (!readLDscore) snpi->ldsc = 0.0;
            snpi->numNonZeroLD = snpi->windSize;
            
            for (unsigned j=0, incj=0; j<windSizeLDM[i]; ++j) {
                snpj = snpInfoVec[windStartLDM[i]+j];
                if (snpj->included) {
                    ZPZ[inci][incj++] = v[j];
                    rsq = v[j]*v[j];
                    snpi->ldSamplVar += (1.0f-rsq)*(1.0f-rsq)/snpi->sampleSize;
                    snpi->ldSum += v[j];
                    if (!readLDscore) snpi->ldsc += rsq;
                    if (snpj == snpi)
                        ZPZdiag[inci] = v[j];
                }
            }
            
            if (++inci == numIncdSnps) break;
        }
    }
    
    fclose(in);
    
    timer.getTime();
    
//    cout << "Window width " << windowWidth << " Mb." << endl;
    displayAverageWindowSize(windSize);
    cout << "LD matrix diagonal mean " << ZPZdiag.mean() << " sd " << sqrt(Gadget::calcVariance(ZPZdiag)) << "." << endl;
    if (ZPZdiag.mean() < 0.8 || ZPZdiag.mean() > 1.2) throw("ERROR: The mean of LD matrix diagonal values is expected to be close to one. Something is wrong with the LD matrix!");
    cout << "Read LD matrix for " << numIncdSnps << " SNPs (time used: " << timer.format(timer.getElapse()) << ")." << endl;
}

void Data::readLDmatrixBinFileAndShrink(const string &ldmatrixFile){
    
    Gadget::Timer timer;
    timer.setTime();
    
    Gadget::Tokenizer token;
    token.getTokens(ldmatrixFile, ".");
    string ldmType = token[token.size()-2];
    sparseLDM = ldmType == "sparse" ? true : false;
    
    VectorXi windStartLDM(numSnps);
    VectorXi windSizeLDM(numSnps);
    
    windStart.resize(numIncdSnps);
    windSize.resize(numIncdSnps);
    
    SnpInfo *snpi, *snpj;
    
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snpi = snpInfoVec[i];
        windStartLDM[i] = snpi->windStart;
        windSizeLDM[i]  = snpi->windSize;
    }
    
    FILE *in = fopen(ldmatrixFile.c_str(), "rb");
    if (!in) {
        throw("Error: cannot open LD matrix file " + ldmatrixFile);
    }
    
    if (!sparseLDM) resizeWindow(incdSnpInfoVec, windStartLDM, windSizeLDM, windStart, windSize);
    
    if (numIncdSnps == 0) throw ("Error: No SNP is retained for analysis.");
    
    cout << "Reading and shrinking " + ldmType + " LD matrix from [" + ldmatrixFile + "]..." << endl;
    
    float rsq = 0.0;
    
//    float nref = incdSnpInfoVec[0]->sampleSize;
    float nref = 183;  // this is the sample size for genetic map
    float n = 2.0*nref-1.0;
    // Approximation to the harmonic series
    float nsum = log(n) + 0.5772156649 + 1.0 / (2.0 * n) - 1.0 / (12.0 * pow(n, 2.0)) + 1.0 / (120.0 * pow(n, 4.0));
    float theta = (1.0 / nsum) / (2.0 * (nref) + 1.0 / nsum);
    float off  = (1.0 - theta)*(1.0 - theta);
    float diag = 0.5*theta*(1.0 - 0.5*theta);
    
    float shrunkLD;
    float rho;
    float Ne = 11490.672741;
    float shrinkage;
    float sdprod;
    
    if (sparseLDM) {
        ZPZsp.resize(numIncdSnps);
        ZPZdiag.resize(numIncdSnps);
        
        for (unsigned i=0, inci=0; i<numSnps; i++) {
            snpi = snpInfoVec[i];
            
            unsigned d[windSizeLDM[i]];
            float v[windSizeLDM[i]];
            
            if (!snpi->included) {
                fseek(in, sizeof(d), SEEK_CUR);
                fseek(in, sizeof(v), SEEK_CUR);
                continue;
            }
            
            fread(d, sizeof(d), 1, in);
            fread(v, sizeof(v), 1, in);
            
            ZPZsp[inci].resize(windSizeLDM[i]);
            snpi->ldSamplVar = 0.0;
            snpi->ldSum = 0.0;
            if (!readLDscore) snpi->ldsc = 0.0;
            
            for (unsigned j=0; j<windSizeLDM[i]; ++j) {
                snpj = snpInfoVec[d[j]];
                if (snpj->included) {
                    
                    sdprod = sqrt(snpi->twopq*snpj->twopq);
                    shrunkLD = v[j]*sdprod;
                    rho = 4.0 * Ne * abs(snpi->genPos - snpj->genPos)/100.0;
                    shrinkage = exp(-rho / (2.0*nref)) * off;
                    if (shrinkage <= 1e-5) shrinkage = 0.0;
                    shrunkLD *= shrinkage;
                    shrunkLD /= sdprod;
                    if (snpj == snpi) {
                        shrunkLD += diag/sdprod;
                        ZPZdiag[inci] = shrunkLD;
                    }
                    
                    ZPZsp[inci].insertBack(snpj->index) = shrunkLD;
                    rsq = shrunkLD * shrunkLD;
                    snpi->ldSamplVar += shrinkage*shrinkage*(1.0f-rsq)*(1.0f-rsq)/snpi->sampleSize;
                    snpi->ldSum += shrunkLD;
                    if (!readLDscore) snpi->ldsc += rsq;
                }
            }
            SparseVector<float>::InnerIterator it(ZPZsp[inci]);
            windStart[inci] = snpi->windStart = it.index();
            windSize[inci] = snpi->windSize = ZPZsp[inci].nonZeros();
            for (; it; ++it) snpi->windEnd = it.index();
            snpi->numNonZeroLD = snpi->windSize;
            
            if (++inci == numIncdSnps) break;
        }
    }
    else {
        ZPZ.resize(numIncdSnps);
        ZPZdiag.resize(numIncdSnps);
        
        for (unsigned i=0, inci=0; i<numSnps; i++) {
            snpi = snpInfoVec[i];
            
            float v[windSizeLDM[i]];
            
            if (!snpi->included) {
                fseek(in, sizeof(v), SEEK_CUR);
                continue;
            }
            
            fread(v, sizeof(v), 1, in);
            
            ZPZ[inci].resize(windSize[inci]);
            snpi->ldSamplVar = 0.0;
            snpi->ldSum = 0.0;
            if (!readLDscore) snpi->ldsc = 0.0;
            snpi->numNonZeroLD = 0;
            
            for (unsigned j=0, incj=0; j<windSizeLDM[i]; ++j) {
                snpj = snpInfoVec[windStartLDM[i]+j];
                if (snpj->included) {
                    
                    sdprod = sqrt(snpi->twopq*snpj->twopq);
                    shrunkLD = v[j]*sdprod;
                    rho = 4.0 * Ne * abs(snpi->genPos - snpj->genPos)/100.0;
                    shrinkage = exp(-rho / (2.0*nref)) * off;
                    if (shrinkage <= 1e-5) shrinkage = 0.0;
                    else snpi->numNonZeroLD++;
                    shrunkLD *= shrinkage;
                    shrunkLD /= sdprod;
                    if (snpj == snpi) {
                        shrunkLD += diag/sdprod;
                        ZPZdiag[inci] = shrunkLD;
                    }
                    
                    ZPZ[inci][incj++] = shrunkLD;
                    rsq = shrunkLD * shrunkLD;
                    snpi->ldSamplVar += shrinkage*shrinkage*(1.0f-rsq)*(1.0f-rsq)/snpi->sampleSize;
                    snpi->ldSum += shrunkLD;
                    if (!readLDscore) snpi->ldsc += rsq;
                }
            }
            
            if (++inci == numIncdSnps) break;
        }
    }
    
    fclose(in);
    
    timer.getTime();
    
    //    cout << "Window width " << windowWidth << " Mb." << endl;
    displayAverageWindowSize(windSize);
    cout << "LD matrix diagonal mean " << ZPZdiag.mean() << " sd " << sqrt(Gadget::calcVariance(ZPZdiag)) << "." << endl;
    if (ZPZdiag.mean() < 0.8 || ZPZdiag.mean() > 1.2) throw("ERROR: The mean of LD matrix diagonal values is expected to be close to one. Something is wrong with the LD matrix!");
    cout << "Read LD matrix for " << numIncdSnps << " SNPs (time used: " << timer.format(timer.getElapse()) << ")." << endl;
}

void Data::readMultiLDmatInfoFile(const string &mldmatFile){
    ifstream in(mldmatFile.c_str());
    if (!in) throw ("Error: can not open the file [" + mldmatFile + "] to read.");
    cout << "Reading SNP info from [" + mldmatFile + "]..." << endl;
    string inputStr;
    numSnpMldVec.clear();
    while (getline(in, inputStr)) {
        readLDmatrixInfoFile(inputStr+".info");
        numSnpMldVec.push_back(numSnps);
    }
    SnpInfo *snp = snpInfoVec[numSnpMldVec[0]];
    if (snp->index == 0) reindexed = true;
    else reindexed = false;
}

void Data::readMultiLDmatBinFile(const string &mldmatFile){
    //cout << "Hi I'm in here " << endl;
    vector<string> filenameVec;
    ifstream in1(mldmatFile.c_str());
    if (!in1) throw ("Error: can not open the file [" + mldmatFile + "] to read.");
    
    Gadget::Timer timer;
    timer.setTime();
    
    string inputStr;
    string ldmType;
    sparseLDM = true;
    while (getline(in1, inputStr)) {
        filenameVec.push_back(inputStr + ".bin");
        Gadget::Tokenizer token;
        token.getTokens(inputStr, ".");
        ldmType = token[token.size()-1];
        sparseLDM = ldmType == "sparse" ? true : false;
    }

    cout << "Reading " + ldmType + " LD matrices from [" + mldmatFile + "]..." << endl;

    VectorXi windStartLDM(numSnps);
    VectorXi windSizeLDM(numSnps);
    //cout << "I made it here " << endl; 
    for (unsigned j=0, i=0, cnt=0; j<numSnps; ++j) {
        SnpInfo *snp = snpInfoVec[j];
        if (j==numSnpMldVec[i]) {
            if (reindexed)
                cnt = numSnpMldVec[i++];
            else
                cnt = 0;
        }
        snp->windStart += cnt;
        snp->windEnd   += cnt;
        windSizeLDM[j]  = snp->windSize;
        windStartLDM[j] = snp->windStart;
    }
    
    if (sparseLDM) {
        windStart.setZero(numIncdSnps);
        windSize.setZero(numIncdSnps);
        ZPZsp.resize(numIncdSnps);
    }
    else {
        resizeWindow(incdSnpInfoVec, windStartLDM, windSizeLDM, windStart, windSize);
        ZPZ.resize(numIncdSnps);
    }
    ZPZdiag.resize(numIncdSnps);
    
    unsigned starti = 0;
    unsigned incj = 0;
    //cout << "I made it here " << endl;    
    long numFiles = filenameVec.size();
    for (unsigned i=0; i<numFiles; ++i) {
        FILE *in2 = fopen(filenameVec[i].c_str(), "rb");
        if (!in2) {
            throw("Error: cannot open LD matrix file " + filenameVec[i]);
        }
        
        SnpInfo *snpj = NULL;
        SnpInfo *snpk = NULL;
        
        float rsq = 0.0;
        
        if (sparseLDM) {
            for (unsigned j=starti; j<numSnpMldVec[i]; j++) {
                snpj = snpInfoVec[j];
                
                unsigned d[windSizeLDM[j]];
                float v[windSizeLDM[j]];
                
                if (!snpj->included) {
                    fseek(in2, sizeof(d), SEEK_CUR);
                    fseek(in2, sizeof(v), SEEK_CUR);
                    continue;
                }
                
                fread(d, sizeof(d), 1, in2);
                fread(v, sizeof(v), 1, in2);
                
                ZPZsp[incj].resize(windSizeLDM[j]);
                snpj->ldSamplVar = 0.0;
                snpj->ldSum = 0.0;
                if (!readLDscore) snpj->ldsc = 0.0;

                for (unsigned k=0; k<windSizeLDM[j]; ++k) {
                    snpk = snpInfoVec[windStartLDM[j]+d[k]-d[0]];
                    if (snpk->included) {
                        ZPZsp[incj].insertBack(snpk->index) = v[k];
                        rsq = v[k]*v[k];
                        snpj->ldSamplVar += (1.0f-rsq)*(1.0f-rsq)/snpj->sampleSize;
                        snpj->ldSum += v[k];
                        if (!readLDscore) snpj->ldsc += rsq;
                        if (snpk == snpj)
                            ZPZdiag[incj] = v[k];
                    }
                }
                SparseVector<float>::InnerIterator it(ZPZsp[incj]);
                windStart[incj] = snpj->windStart = it.index();
                windSize[incj] = snpj->windSize = ZPZsp[incj].nonZeros();
                snpj->numNonZeroLD = snpj->windSize;
                ++incj;
            }
        }
        else {
            //cout << "I's reading the snps " << endl;
            for (unsigned j=starti; j<numSnpMldVec[i]; j++) {
                snpj = snpInfoVec[j];
                float v[windSizeLDM[j]];
	        //cout << "I'm at SNP " << j << " windStartLDM[j] " << windStartLDM[j] << " numSnpMldVec[i] " << numSnpMldVec[i] << " windSizeLDM[j] " << windSizeLDM[j] << endl;
                if (!snpj->included) {
                    fseek(in2, sizeof(v), SEEK_CUR);
                    continue;
                }
                
                fread(v, sizeof(v), 1, in2);
                //cout << "Here 1 " << endl;                 
                ZPZ[incj].resize(windSize[incj]);
                snpj->ldSamplVar = 0.0;
                snpj->ldSum = 0.0;
                if (!readLDscore) snpj->ldsc = 0.0;
                snpj->numNonZeroLD = snpj->windSize;

                for (unsigned k=0, inck=0; k<windSizeLDM[j]; ++k) {
                    //cout << "I'm at k " << k << endl; 
                    snpk = snpInfoVec[windStartLDM[j]+k];
                    if (snpk->included) {
                        ZPZ[incj][inck++] = v[k];
                        rsq = v[k]*v[k];
                        snpj->ldSamplVar += (1.0f-rsq)*(1.0f-rsq)/snpj->sampleSize;
                        snpj->ldSum += v[k];
                        if (!readLDscore) snpj->ldsc += rsq;
                        if (snpk == snpj)
                            ZPZdiag[incj] = v[k];
                    }
                }
                 //cout << "Here 2 " << endl;
                ++incj;
            }
        }
        
        fclose(in2);
        cout << "Read " + ldmType + " LD matrix for " << numSnpMldVec[i]-starti << " SNPs from [" << filenameVec[i] << "]." << endl;
        
        starti = numSnpMldVec[i];
    }
    
    timer.getTime();
    
    displayAverageWindowSize(windSize);
    cout << "LD matrix diagonal mean " << ZPZdiag.mean() << " sd " << sqrt(Gadget::calcVariance(ZPZdiag)) << "." << endl;
    if (ZPZdiag.mean() < 0.8 || ZPZdiag.mean() > 1.2) throw("ERROR: The mean of LD matrix diagonal values is expected to be close to one. Something is wrong with the LD matrix!");
    cout << "Read LD matrix for " << numIncdSnps << " SNPs (time used: " << timer.format(timer.getElapse()) << ")." << endl;
}

void Data::readMultiLDmatBinFileAndShrink(const string &mldmatFile, const float genMapN){
    //cout << "Hi I'm in here " << endl;
    vector<string> filenameVec;
    ifstream in1(mldmatFile.c_str());
    if (!in1) throw ("Error: can not open the file [" + mldmatFile + "] to read.");
    
    Gadget::Timer timer;
    timer.setTime();
    
    string inputStr;
    string ldmType;
    sparseLDM = true;
    while (getline(in1, inputStr)) {
        filenameVec.push_back(inputStr + ".bin");
        Gadget::Tokenizer token;
        token.getTokens(inputStr, ".");
        ldmType = token[token.size()-1];
        sparseLDM = ldmType == "sparse" ? true : false;
    }
    
    cout << "Reading and shrinking " + ldmType + " LD matrices from [" + mldmatFile + "]..." << endl;
    
    VectorXi windStartLDM(numSnps);
    VectorXi windSizeLDM(numSnps);
    //cout << "I made it here " << endl;
    for (unsigned j=0, i=0, cnt=0; j<numSnps; ++j) {
        SnpInfo *snp = snpInfoVec[j];
        if (j==numSnpMldVec[i]) {
            if (reindexed)
                cnt = numSnpMldVec[i++];
            else
                cnt = 0;
        }
        snp->windStart += cnt;
        snp->windEnd   += cnt;
        windSizeLDM[j]  = snp->windSize;
        windStartLDM[j] = snp->windStart;
    }
    
    if (sparseLDM) {
        windStart.setZero(numIncdSnps);
        windSize.setZero(numIncdSnps);
        ZPZsp.resize(numIncdSnps);
    }
    else {
        resizeWindow(incdSnpInfoVec, windStartLDM, windSizeLDM, windStart, windSize);
        ZPZ.resize(numIncdSnps);
    }
    ZPZdiag.resize(numIncdSnps);
    
    unsigned starti = 0;
    unsigned incj = 0;
    //cout << "I made it here " << endl;
    long numFiles = filenameVec.size();
    for (unsigned i=0; i<numFiles; ++i) {
        FILE *in2 = fopen(filenameVec[i].c_str(), "rb");
        if (!in2) {
            throw("Error: cannot open LD matrix file " + filenameVec[i]);
        }
        
        SnpInfo *snpj = NULL;
        SnpInfo *snpk = NULL;
        
        float rsq = 0.0;
        
        float nref = genMapN;
        float n = 2.0*nref-1.0;
        // Approximation to the harmonic series
        float nsum = log(n) + 0.5772156649 + 1.0 / (2.0 * n) - 1.0 / (12.0 * pow(n, 2.0)) + 1.0 / (120.0 * pow(n, 4.0));
        float theta = (1.0 / nsum) / (2.0 * (nref) + 1.0 / nsum);
        float off  = (1.0 - theta)*(1.0 - theta);
        float diag = 0.5*theta*(1.0 - 0.5*theta);
        
        float shrunkLD;
        float rho;
        float Ne = 11490.672741;
        float shrinkage;
        float sdprod;

        if (sparseLDM) {
            for (unsigned j=starti; j<numSnpMldVec[i]; j++) {
                snpj = snpInfoVec[j];
                
                unsigned d[windSizeLDM[j]];
                float v[windSizeLDM[j]];
                
                if (!snpj->included) {
                    fseek(in2, sizeof(d), SEEK_CUR);
                    fseek(in2, sizeof(v), SEEK_CUR);
                    continue;
                }
                
                fread(d, sizeof(d), 1, in2);
                fread(v, sizeof(v), 1, in2);
                
                ZPZsp[incj].resize(windSizeLDM[j]);
                snpj->ldSamplVar = 0.0;
                snpj->ldSum = 0.0;
                if (!readLDscore) snpj->ldsc = 0.0;
                
                for (unsigned k=0; k<windSizeLDM[j]; ++k) {
                    snpk = snpInfoVec[windStartLDM[j]+d[k]-d[0]];
                    if (snpk->included) {
                        
                        sdprod = sqrt(snpj->twopq*snpk->twopq);
                        shrunkLD = v[k]*sdprod;
                        rho = 4.0 * Ne * abs(snpj->genPos - snpk->genPos)/100.0;
                        shrinkage = exp(-rho / (2.0*nref)) * off;
                        if (shrinkage <= 1e-5) shrinkage = 0.0;
                        shrunkLD *= shrinkage;
                        shrunkLD /= sdprod;
                        if (snpk == snpj) {
                            shrunkLD += diag/sdprod;
                            ZPZdiag[incj] = shrunkLD;
                        }
                        
                        ZPZsp[incj].insertBack(snpk->index) = shrunkLD;
                        rsq = shrunkLD * shrunkLD;
                        snpj->ldSamplVar += shrinkage*shrinkage*(1.0f-rsq)*(1.0f-rsq)/snpj->sampleSize;
                        snpj->ldSum += shrunkLD;
                        if (!readLDscore) snpj->ldsc += rsq;
                    }
                }
                SparseVector<float>::InnerIterator it(ZPZsp[incj]);
                windStart[incj] = snpj->windStart = it.index();
                windSize[incj] = snpj->windSize = ZPZsp[incj].nonZeros();
                snpj->numNonZeroLD = snpj->windSize;
                ++incj;
            }
        }
        else {
            //cout << "I's reading the snps " << endl;
            for (unsigned j=starti; j<numSnpMldVec[i]; j++) {
                snpj = snpInfoVec[j];
                float v[windSizeLDM[j]];
                //cout << "I'm at SNP " << j << " windStartLDM[j] " << windStartLDM[j] << " numSnpMldVec[i] " << numSnpMldVec[i] << " windSizeLDM[j] " << windSizeLDM[j] << endl;
                if (!snpj->included) {
                    fseek(in2, sizeof(v), SEEK_CUR);
                    continue;
                }
                
                fread(v, sizeof(v), 1, in2);
                //cout << "Here 1 " << endl;
                ZPZ[incj].resize(windSize[incj]);
                snpj->ldSamplVar = 0.0;
                snpj->ldSum = 0.0;
                if (!readLDscore) snpj->ldsc = 0.0;
                snpj->numNonZeroLD = 0;
                
                for (unsigned k=0, inck=0; k<windSizeLDM[j]; ++k) {
                    //cout << "I'm at k " << k << endl;
                    snpk = snpInfoVec[windStartLDM[j]+k];
                    if (snpk->included) {
                        
                        sdprod = sqrt(snpj->twopq*snpk->twopq);
                        shrunkLD = v[k]*sdprod;
                        rho = 4.0 * Ne * abs(snpj->genPos - snpk->genPos)/100.0;
                        shrinkage = exp(-rho / (2.0*nref)) * off;
                        if (shrinkage <= 1e-5) shrinkage = 0.0;
                        else snpj->numNonZeroLD++;
                        shrunkLD *= shrinkage;
                        shrunkLD /= sdprod;
                        if (snpk == snpj) {
                            shrunkLD += diag/sdprod;
                            ZPZdiag[incj] = shrunkLD;
                        }
                        
                        ZPZ[incj][inck++] = shrunkLD;
                        rsq = shrunkLD * shrunkLD;
                        snpj->ldSamplVar += shrinkage*shrinkage*(1.0f-rsq)*(1.0f-rsq)/snpj->sampleSize;
                        snpj->ldSum += shrunkLD;
                        if (!readLDscore) snpj->ldsc += rsq;
                    }
                }
                //cout << "Here 2 " << endl;
                ++incj;
            }
        }
        
        fclose(in2);
        cout << "Read " + ldmType + " LD matrix for " << numSnpMldVec[i]-starti << " SNPs from [" << filenameVec[i] << "]." << endl;
        
        starti = numSnpMldVec[i];
    }
    
    timer.getTime();
    
    displayAverageWindowSize(windSize);
    cout << "LD matrix diagonal mean " << ZPZdiag.mean() << " sd " << sqrt(Gadget::calcVariance(ZPZdiag)) << "." << endl;
    if (ZPZdiag.mean() < 0.8 || ZPZdiag.mean() > 1.2) throw("ERROR: The mean of LD matrix diagonal values is expected to be close to one. Something is wrong with the LD matrix!");
    cout << "Read LD matrix for " << numIncdSnps << " SNPs (time used: " << timer.format(timer.getElapse()) << ")." << endl;
}

void Data::resizeLDmatrix(const string &LDmatType, const float chisqThreshold, const unsigned windowWidth, const float LDthreshold, const float effpopNE, const float cutOff, const float genMapN) {
    if (LDmatType == "full") {
//        for (unsigned i=0; i<numIncdSnps; ++i) {
//            SnpInfo *snp = incdSnpInfoVec[i];
//            VectorXf rsq = ZPZ[i].cwiseProduct(ZPZ[i]);
//            VectorXf rsq_adj = rsq - (VectorXf::Ones(snp->windSize) - rsq)/float(snp->sampleSize-2);
//            snp->ldSamplVar = (VectorXf::Ones(snp->windSize) - rsq_adj).squaredNorm()/snp->sampleSize;
//            snp->ldSum = ZPZ[i].sum();
//        }
        return;
    }
    //if (LDmatType == "shrunk") return;  // TMP; to be removed
    snp2pq.resize(numIncdSnps);
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        snp2pq[i] = snp->twopq = 2.0*snp->af*(1.0-snp->af);
    }
    float rsq = 0.0;    
    if (LDmatType == "sparse") {
        if (ZPZsp.size() == 0) {
            cout << "Making a sparse LD matrix by setting the non-significant LD to be zero..." << endl;
            ZPZsp.resize(numIncdSnps);
            SnpInfo *snpi, *snpj;
            if (LDthreshold) {
                for (unsigned i=0; i<numIncdSnps; ++i) {
                    snpi = incdSnpInfoVec[i];
                    ZPZsp[i].resize(snpi->windSize);
                    snpi->ldSamplVar = 0.0;
                    snpi->ldSum = 0.0;
                    for (unsigned j=0; j<snpi->windSize; ++j) {
                        snpj = incdSnpInfoVec[snpi->windStart + j];
                        if (abs(ZPZ[i][j]) > LDthreshold || snpj->skeleton) {
                            ZPZsp[i].insertBack(snpi->windStart + j) = ZPZ[i][j];
                            rsq = ZPZ[i][j]*ZPZ[i][j];
                            snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snpi->sampleSize;
                            snpi->ldSum += ZPZ[i][j];
                        }
                    }
                    SparseVector<float>::InnerIterator it(ZPZsp[i]);
                    windStart[i] = snpi->windStart = it.index();
                    windSize[i] = snpi->windSize = ZPZsp[i].nonZeros();
                    for (; it; ++it) snpi->windEnd = it.index();
                    ZPZ[i].resize(0);
                }
            } else {
                if (windowWidth) {
                    for (unsigned i=0; i<numIncdSnps; ++i) {
                        snpi = incdSnpInfoVec[i];
                        ZPZsp[i].resize(snpi->windSize);
                        snpi->ldSamplVar = 0.0;
                        snpi->ldSum = 0.0;
                        for (unsigned j=0; j<snpi->windSize; ++j) {
                            snpj = incdSnpInfoVec[snpi->windStart + j];
                            if (i==j || ZPZ[i][j]*ZPZ[i][j]*snpi->sampleSize > chisqThreshold ||
                                snpi->isProximal(*incdSnpInfoVec[snpi->windStart + j], windowWidth/2)) {
                                ZPZsp[i].insertBack(snpi->windStart + j) = ZPZ[i][j];
                                rsq = ZPZ[i][j]*ZPZ[i][j];
                                snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snpi->sampleSize;
                                snpi->ldSum += ZPZ[i][j];
                            }
                        }
                        //            ZPZsp[i] = ZPZ[i].sparseView();
                        SparseVector<float>::InnerIterator it(ZPZsp[i]);
                        windStart[i] = snpi->windStart = it.index();
                        windSize[i] = snpi->windSize = ZPZsp[i].nonZeros();
                        for (; it; ++it) snpi->windEnd = it.index();
                        ZPZ[i].resize(0);
                        //cout << i << " windsize " << snp->windSize << " " << ZPZsp[i].size() << endl;
                    }
                } else {
                    cout << "Using a chisq threshold of " << chisqThreshold << endl; 
                    for (unsigned i=0; i<numIncdSnps; ++i) {
                        snpi = incdSnpInfoVec[i];
                        ZPZsp[i].resize(snpi->windSize);
                        snpi->ldSamplVar = 0.0;
                        snpi->ldSum = 0.0;
                        for (unsigned j=0; j<snpi->windSize; ++j) {
                            snpj = incdSnpInfoVec[snpi->windStart + j];
                            if (i==j || ZPZ[i][j]*ZPZ[i][j]*snpi->sampleSize > chisqThreshold) {
                                ZPZsp[i].insertBack(snpi->windStart + j) = ZPZ[i][j];
                                rsq = ZPZ[i][j]*ZPZ[i][j];
                                snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snpi->sampleSize;
                                snpi->ldSum += ZPZ[i][j];
                            }
                        }
                        //            ZPZsp[i] = ZPZ[i].sparseView();
                        SparseVector<float>::InnerIterator it(ZPZsp[i]);
                        windStart[i] = snpi->windStart = it.index();
                        windSize[i] = snpi->windSize = ZPZsp[i].nonZeros();
                        for (; it; ++it) snpi->windEnd = it.index();
                        ZPZ[i].resize(0);
                        //cout << i << " windsize " << snp->windSize << " " << ZPZsp[i].size() << endl;
                    }
                }
            }
        } else {
            cout << "Pruning a sparse LD matrix by chisq threshold of " << chisqThreshold << endl;
            SnpInfo *snpi, *snpj;
            for (unsigned i=0; i<numIncdSnps; ++i) {
                snpi = incdSnpInfoVec[i];
                snpi->ldSamplVar = 0.0;
                snpi->ldSum = 0.0;
                for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                    snpj = incdSnpInfoVec[it.index()];
                    rsq = it.value()*it.value();
                    if (rsq*snpi->sampleSize <= chisqThreshold) it.valueRef() = 0.0;
                    else {
                        snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snpi->sampleSize;
                        snpi->ldSum += it.value();
                    }
                }
                ZPZsp[i].prune(0.0);
                SparseVector<float>::InnerIterator it(ZPZsp[i]);
                windStart[i] = snpi->windStart = it.index();
                windSize[i] = snpi->windSize = ZPZsp[i].nonZeros();
                for (; it; ++it) snpi->windEnd = it.index();
            }
        }
    }
    if (LDmatType == "band") {
        VectorXi windStartOri = windStart;
        VectorXi windEndi;
        windEndi.resize(numIncdSnps);
        VectorXi windSizeOri = windSize;
        VectorXf ZPZiTmp;
        if (windowWidth) {
            cout << "Resizing LD matrix based on a window width of " << windowWidth*1e-6 << " Mb..." << endl;
            getWindowInfo(incdSnpInfoVec, windowWidth, windStart, windSize);
            for (unsigned i=0; i<numIncdSnps; ++i) {
                SnpInfo *snp = incdSnpInfoVec[i];
                windStart[i] = snp->windStart = max(windStart[i], windStartOri[i]);
                windSize[i]  = snp->windSize  = min(windSize[i], windSizeOri[i]);
                ZPZiTmp = ZPZ[i].segment(windStart[i]-windStartOri[i], windSize[i]);
                ZPZ[i] = ZPZiTmp;
                snp->ldSamplVar = (1.0 - ZPZ[i].array().square()).square().sum()/snp->sampleSize;
                snp->ldSum = ZPZ[i].sum();
            }
        } else if (LDthreshold) {
            cout << "Resizing LD matrix based on a LD threshold of " << LDthreshold << "..." << endl;
            for (unsigned i=0; i<numIncdSnps; ++i) {
                windEndi[i] = windSizeOri[i];
                // Gather up the window ends
                for (unsigned j=windSizeOri[i]; j>i; --j) {
                    if (abs(ZPZ[i][j-1]) > LDthreshold) {
                        windEndi[i] = j;
                        break;
                    }
                }
                // Fill the lower triangle with the necessary zeroes
                //if (windEndi != numIncdSnps) {
                for (unsigned j=windEndi[i]; j<numIncdSnps; ++j) {
                      ZPZ[j][i] = 0.0;
                }
                //}
            }
            // Cycle again over the vectors and get the window start
            for (unsigned i=0; i<numIncdSnps; ++i) {
                SnpInfo *snp = incdSnpInfoVec[i];
                // unsigned windEndi = windSizeOri[i];
                // Run over the lower triangle and make equal
                for (unsigned j=0; j<(i-1); ++j) {
                    if (abs(ZPZ[i][j]) > 0.0) {
                        windStart[i] = snp->windStart = windStartOri[i] + j;
                        break;
                    }
                }
                // Resize the vectors of vectors accordingly
                // cout << "Wind start, Wind end " << windStart[i] << ", "<< windEndi[i] << endl;
                windSize[i] = snp->windSize = windStartOri[i] + windEndi[i] - windStart[i];
                ZPZiTmp = ZPZ[i].segment(windStart[i] - windStartOri[i], windSize[i]);
                ZPZ[i] = ZPZiTmp;
            }
        } else {
            cout << "Resizing LD matrix based on a chisquare threshold of " << chisqThreshold << "..." << endl;
            for (unsigned i=0; i<numIncdSnps; ++i) {
                windEndi[i] = windSizeOri[i];
                 SnpInfo *snp = incdSnpInfoVec[i];
                // Gather up the window ends
                for (unsigned j=windSizeOri[i]; j>i; --j) {
                    if (ZPZ[i][j]*ZPZ[i][j]*snp->sampleSize > chisqThreshold) {
                        windEndi[i] = j;
                        break;
                    }
                }
                // Fill the lower triangle with the necessary zeroes
                //if (windEndi != numIncdSnps) {
                for (unsigned j=windEndi[i]; j<numIncdSnps; ++j) {
                      ZPZ[j][i] = 0.0;
                }
                //}
            }
            // Cycle again over the vectors and get the window start
            for (unsigned i=0; i<numIncdSnps; ++i) {
                SnpInfo *snp = incdSnpInfoVec[i];
                // unsigned windEndi = windSizeOri[i];
                // Run over the lower triangle and make equal
                for (unsigned j=0; j<(i-1); ++j) {
                    if (abs(ZPZ[i][j]) > 0.0) {
                        windStart[i] = snp->windStart = windStartOri[i] + j;
                        break;
                    }
                }
                // Resize the vectors of vectors accordingly
                // cout << "Wind start, Wind end " << windStart[i] << ", "<< windEndi[i] << endl;
                windSize[i] = snp->windSize = windStartOri[i] + windEndi[i] - windStart[i];
                ZPZiTmp = ZPZ[i].segment(windStart[i] - windStartOri[i], windSize[i]);
                ZPZ[i] = ZPZiTmp;
                snp->ldSamplVar = (1.0 - ZPZ[i].array().square()).square().sum()/snp->sampleSize;
                snp->ldSum = ZPZ[i].sum();
            }
//        } else {
//            cout << "Resizing LD matrix based on a chisq threshold of " << chisqThreshold << "..." << endl;
//            for (unsigned i=0; i<numIncdSnps; ++i) {
//                SnpInfo *snp = incdSnpInfoVec[i];
//                unsigned windEndi = windSizeOri[i];
//                for (unsigned j=0; j<windSizeOri[i]; ++j) {
//                    if (ZPZ[i][j]*ZPZ[i][j]*snp->sampleSize > chisqThreshold) {
//                        windStart[i] = snp->windStart = windStartOri[i] + j;
//                        break;
//                    }
//                }
//                for (unsigned j=windSizeOri[i]; j>0; --j) {
//                    if (ZPZ[i][j]*ZPZ[i][j]*snp->sampleSize > chisqThreshold) {
//                        windEndi = j;
//                        break;
//                    }
//                }
//                windSize[i] = snp->windSize = windStartOri[i] + windEndi - windStart[i];
//                ZPZiTmp = ZPZ[i].segment(windStart[i] - windStartOri[i], windSize[i]);
//                ZPZ[i] = ZPZiTmp;
//                snp->ldSamplVar = (1.0 - ZPZ[i].array().square()).square().sum()/snp->sampleSize;
//                snp->ldSum = ZPZ[i].sum();
//            }
        }
    }
    if (LDmatType == "shrunk") {
        VectorXi windStartOri = windStart;
        VectorXi windSizeOri = windSize;
        cout << "Resizing LD matrix using shrunk matrix properties ..." << endl;
        // ----------------------------------------------------
        // NEW - Calculate the mutation rate components 
        // ----------------------------------------------------
        // m is the number of individuals in the reference panel for each variant
        // Need to compute theta, which is related to the mutation rate
        VectorXf nmsumi;
        VectorXf thetai;
        VectorXf mi;
        VectorXf gmapi;
        VectorXf sdss;
        //cout << "Snps size " << incdSnpInfoVec.size() << endl;
        nmsumi.resize(numIncdSnps);
        thetai.resize(numIncdSnps);
        sdss.resize(numIncdSnps);
        // mi.resize(numIncdSnps);
        cout << "\nUsing genetic map sample size of " << genMapN << " please alter with --genmap-n if inappropriate." << endl;
        float m = genMapN;
        gmapi.resize(numIncdSnps);
        for (unsigned i=0; i<numIncdSnps; ++i) {
            SnpInfo *snp = incdSnpInfoVec[i];
            //mi[i] = (snp->sampleSize);
            int  n = 2 * m - 1;
//            cout << "snp " << i << " sample size " << snp->sampleSize << endl;
            // Approximation to the harmonic series
            nmsumi[i] = log(n) + 0.5772156649 + 1.0/ (2.0 * n) - 1.0 / (12.0 * pow(n, 2)) + 1.0 / (120.0 * pow(n, 4));
            //cout << nmsumi[i] << endl;
            // Calculate theta
            //thetai[i] = (1.0 / nmsumi[i]) / (2.0 * (snp->sampleSize) + 1.0 / nmsumi[i]);
            thetai[i] = (1.0 / nmsumi[i]) / (2.0 * m + 1.0 / nmsumi[i]);
            //cout <<  thetai[i] << endl;
            // Pull out the standard deviation for each variant
            sdss[i] = sqrt(2.0 * (snp->af) * (1.0 - (snp->af)));
            // cout << "snp " << i << " af " << snp->af << endl;
            // 
            gmapi[i] = snp->genPos;
        }
        long int nmsum;
        float theta;
        // cout << "I made it here  for shrunk sparse " << endl;
        
        // Rescale Z to the covariance scale
        VectorXf sdssSub;
        for (unsigned i=0; i<numIncdSnps; ++i) {
            sdssSub = sdss.segment(windStart[i], windSizeOri[i]);
            ZPZ[i] = (sdss[i] / 2.0) *  sdssSub.array() * ZPZ[i].array();
        }
        // cout << "I made it here  for shrunk sparse " << endl;
        // --------------------------------------------------------
        // Compute the shrinkage value and then shrink the elements
        // --------------------------------------------------------
        // cout << "Made it here 2 " << endl;
        float mapdiffi; 
        float rho;
        float shrinkage;
        float Ne = effpopNE;
        cout << "Using European effective population size Ne=" << Ne << " please alter with --ne if inappropriate. " << endl;
        float cutoff = cutOff;
        for (unsigned i=0; i<numIncdSnps; ++i) {
            if (!(i%1000)) cout << i << " SNPs processed\r";
            SnpInfo *snp = incdSnpInfoVec[i];
            for (unsigned j=i; j<=snp->windEnd; ++j) {
                mapdiffi = abs(gmapi[j] - gmapi[i]);
                rho = 4 * Ne * (mapdiffi / 100);
                shrinkage = exp(-rho / (2 * m));
                if (shrinkage <= cutoff)
                {
                    shrinkage = 0.0;
                }
                // // Multiple each covariance matrix element with the shrinkage value
                float value = ZPZ[i][j];
                if (i != (windStart[i] + j)) value = value * shrinkage;
                // // Complete as SigHAat from Li and Stephens 2003
                value =  value * ((1.0 - thetai[i]) * (1.0 - thetai[(windStart[i] + j)]));
                // If it's the diagonal element add the extra term
                if (i == (windStart[i] + j))
                {
                    value = value + 0.5f * thetai[i] * (1.0 - 0.5f * thetai[i]);
                }  
                ZPZ[j][i] = ZPZ[i][j] = value;
            }
            if(!(i%1000)) cout << " Completed snp " << i << "\r" << flush;
        }
        // Now back to correlation
        for (unsigned i=0; i<numIncdSnps; ++i) {
            sdssSub = sdss.segment(windStart[i], windSizeOri[i]);
            ZPZ[i]  = (2.0 / sdss[i]) *  (1.0 / sdssSub.array()) * ZPZ[i].array();
        }
    }
    if (LDmatType == "sparseshrunk") {
        cout << "Resizing sparse LD matrix using shrunk matrix properties ..." << endl;
        // ----------------------------------------------------
        // NEW - Calculate the mutation rate components 
        // ----------------------------------------------------
        // m is the number of individuals in the reference panel for each variant
        // Need to compute theta, which is related to the mutation rate
        VectorXf nmsumi;
        VectorXf thetai;
        VectorXf mi;
        VectorXf gmapi;
        VectorXf sdss;
        //cout << "Snps size " << incdSnpInfoVec.size() << endl;
        nmsumi.resize(numIncdSnps);
        thetai.resize(numIncdSnps);
        sdss.resize(numIncdSnps);
        // mi.resize(numIncdSnps);
        cout << "\nUsing genetic map sample size of " << genMapN << " please alter with --genmap-n if inappropriate." << endl;
        float m = genMapN;
        gmapi.resize(numIncdSnps);
        for (unsigned i=0; i<numIncdSnps; ++i) {
            SnpInfo *snp = incdSnpInfoVec[i];
            // mi[i] = (snp->sampleSize);
            int  n = 2.0 * m - 1.0;
            // cout << "snp " << i << " sample size " << snp->sampleSize << endl;
            // Approximation to the harmonic series
            nmsumi[i] = log(n) + 0.5772156649 + 1.0 / (2.0 * n) - 1.0 / (12.0 * pow(n, 2)) + 1.0 / (120.0 * pow(n, 4));
            //cout << nmsumi[i] << endl;
            // Calculate theta
//            thetai[i] = (1.0 / nmsumi[i]) / (2.0 * (snp->sampleSize) + 1 / nmsumi[i]);
            thetai[i] = (1.0 / nmsumi[i]) / (2.0 * m + 1 / nmsumi[i]);
            //cout <<  thetai[i] << endl;
            // Pull out the standard deviation for each variant
            sdss[i] = sqrt(2.0 * (snp->af) * (1.0 - (snp->af)));
            //
            gmapi[i] = snp->genPos;
        }
        long int nmsum;
        float theta;
        // Rescale Z to the covariance scale
        float rsq = 0.0; 
        for (unsigned i=0; i<numIncdSnps; ++i) {
            // cout << "Before " << ZPZsp[i] << endl;
            for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                // snpj = incdSnpInfoVec[it.index()];
                rsq = (sdss[i] / 2.0) *  sdss[it.index()] * it.value();
                // cout << "sdss[i] " << sdss[i] << " sdss[it.index()] " << sdss[it.index()] << " it.value() " << it.value() << " it.index() " << it.index() << endl;
                // cout << "rsq " << rsq << endl;
                it.valueRef() = rsq;
            }
            // cout << "After " << ZPZsp[i] << endl;
        }
        // --------------------------------------------------------
        // Compute the shrinkage value and then shrink the elements
        // --------------------------------------------------------
        float mapdiffi; 
        float rho;
        float shrinkage;
        float Ne = effpopNE;
        cout << "Using European effective population size Ne=" << Ne << " please alter with --ne if inappropriate. " << endl;
        float cutoff = cutOff;
        for (unsigned i=0; i<numIncdSnps; ++i) {
            // -----------------------------
            // Shrinkage using sparse matrix
            // -----------------------------
            float ZPZij = 0.0;
            for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                // cout << " j " << j << " (windStart[i] + j) " << (windStart[i] + j) << endl;
                mapdiffi = abs(gmapi[it.index()] - gmapi[i]);
                rho = 4.0 * Ne * (mapdiffi / 100.0);
                shrinkage = exp(-rho / (2 * m)); 
                // cout << "Shrinkage " << shrinkage << endl;
                if (shrinkage <= cutoff)
                {
                    shrinkage = 0.0;
                }
                // Multiple each covariance matrix element with the shrinkage value
                ZPZij = it.value() * shrinkage;
                // Complete as SigHAat from Li and Stephens 2003
                ZPZij =  ZPZij * ((1.0 - thetai[i]) * (1.0 - thetai[it.index()]));
                // If it's the diagonal element add the extra term
                if (i == (it.index()))
                {
                    ZPZij = ZPZij + 0.5f * thetai[i] * (1.0 - 0.5f * thetai[i]);
                } 
                it.valueRef()= ZPZij; 
            }
            // cout << "After " << ZPZsp[i] << endl;
            if(!(i%1000)) cout << " Completed snp " << i << "\r" << flush;
        }
        // // Now back to correlation
        for (unsigned i=0; i<numIncdSnps; ++i) {
            for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                // snpj = incdSnpInfoVec[it.index()];
                rsq = (2.0 / sdss[i]) *  (1.0 / sdss[it.index()]) * it.value();
                // cout << "sdss[i] " << sdss[i] << " sdss[it.index()] " << sdss[it.index()] << " it.value() " << it.value() << " it.index() " << it.index() << endl;
                // cout << "rsq " << rsq << endl;
                it.valueRef() = rsq;
            }
        }
    }
    displayAverageWindowSize(windSize);
}


// =============================================================================================
// Make shrunk matrix start
// =============================================================================================

// =============================================================================================
// Function read the genetic map and pass matched SNPs to main SNP inclusion exclusion tools
// =============================================================================================

void Data::readGeneticMapFile(const string &geneticMapFile){
    ifstream in(geneticMapFile.c_str());
    if (!in) throw ("Error: can not open the file [" + geneticMapFile + "] to read.");
    cout << "Reading genetic map info from [" + geneticMapFile + "]." << endl;
    
    SnpInfo *snp;
    map<string, SnpInfo*>::iterator it;
    string id, gmGenPos, gmPPos;
    unsigned line=0, match=0;
    unsigned incon=0;
    while (in >> id >> gmPPos >> gmGenPos) {
        ++line;
        it = snpInfoMap.find(id);
        if (it == snpInfoMap.end()) continue;
        snp = it->second;
        if (!snp->included) continue;
//        snp->gen_map_ppos = atof(gmPPos.c_str());
        snp->genPos = atof(gmGenPos.c_str());
        ++match;
    }
    in.close();
    
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        if (!snp->included) continue;
        if (snp->genPos == -999 || snp->genPos == 0) {
            //cout << "Who went false snp " << i << endl;
            snp->included = false;
        }
    }
    cout << match << " matched SNPs in the genetic map file (in total " << line << " SNPs)." << endl;
}

void Data::readfreqFile(const string &freqFile){
    ifstream in(freqFile.c_str());
    if (!in) throw ("Error: can not open the file [" + freqFile + "] to read.");
    cout << "Reading allele frequency file from [" + freqFile + "]." << endl;
    
    SnpInfo *snp;
    map<string, SnpInfo*>::iterator it;
    string id, A1, freq;
    unsigned line=0, match=0;
    unsigned incon=0;
    while (in >> id >> A1 >> freq) {
        ++line;
        it = snpInfoMap.find(id);
        if (it == snpInfoMap.end()) continue;
        snp = it->second;
        if (!snp->included) continue;
        snp->af = atof(freq.c_str());
        ++match;
    }
    in.close();
    
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        if (!snp->included) continue;
        if (snp->af == -1) {
            //cout << "Who went false snp " << i << endl;
            snp->included = false;
        }
    }
    cout << match << " matched SNPs in the allele frequency file (in total " << line << " SNPs)." << endl;
}


// =============================================================================================
// Function to build the shrunk matrix 
// =============================================================================================

void Data::makeshrunkLDmatrix(const string &bedFile, const string &LDmatType, const string &snpRange, const string &filename, const bool writeLdmTxt, const float effpopNE, const float cutOff, const float genMapN){
    

    Gadget::Tokenizer token;
    token.getTokens(snpRange, "-");
    
    unsigned start = 0;
    unsigned end = numIncdSnps;
    
    if (token.size()) {
        start = atoi(token[0].c_str()) - 1;
        end = atoi(token[1].c_str());
        if (end > numIncdSnps) end = numIncdSnps;
    }

    
    unsigned numSnpInRange = end - start;
    
    if (snpRange.empty())
        cout << "Building shrunk LD matrix for all SNPs ..." << endl;
    else
        cout << "Building shrunk LD matrix for SNPs " << snpRange << " ..." << endl;
    
    if (numIncdSnps == 0) throw ("Error: No SNP is retained for analysis.");
    if (numKeptInds == 0) throw ("Error: No individual is retained for analysis.");
    if (start >= numIncdSnps) throw ("Error: Specified a SNP range of " + snpRange + " but " + to_string(static_cast<long long>(numIncdSnps)) + " SNPs are included.");
    
    Gadget::Timer timer;
    timer.setTime();
    
    unsigned firstWindStart = 0;
    unsigned lastWindEnd = 0;
    
    // first read in the genotypes of SNPs in the given range
    
    const int bedToGeno[4] = {2, -9, 1, 0};
    unsigned size = (numInds+3)>>2;
    
    MatrixXf ZP(numSnpInRange, numKeptInds);  // SNP x Ind
    D.setZero(numSnpInRange);

    if (numKeptInds < 2) throw("Error: Cannot calculate LD matrix with number of individuals < 2.");

    FILE *in1 = fopen(bedFile.c_str(), "rb");
    if (!in1) throw ("Error: can not open the file [" + bedFile + "] to read.");
    cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
    char header[3];
    fread(header, sizeof(header), 1, in1);
    if (!in1 || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
        cerr << "Error: Incorrect first three bytes of bed file: " << bedFile << endl;
        exit(1);
    }

    
    IndInfo *indi = NULL;
    SnpInfo *snpj = NULL;
    SnpInfo *snpk = NULL;
    
    int genoValue;
    unsigned i, j, k;
    unsigned incj, inck; // index of included SNP
    unsigned long long skipj = 0;
    unsigned nmiss;
    float mean;
    
    for (j = 0, incj = 0; j < numSnps; j++) {
        snpj = snpInfoVec[j];
        // cout << "Genetic map position " << snpj->gen_map_pos << endl;
        if (snpj->index < start || !snpj->included) {
            skipj += size;
            continue;
        }
        
        if (skipj) fseek(in1, skipj, SEEK_CUR);
        skipj = 0;
        
        char *bedLineIn = new char[size];
        fread(bedLineIn, sizeof(char), size, in1);
        
        mean = 0.0;
        nmiss = 0;
        
        for (i = 0; i < numInds; i++) {
            indi = indInfoVec[i];
            if (!indi->kept) continue;
            genoValue = bedToGeno[(bedLineIn[i>>2]>>((i&3)<<1))&3];
            ZP(incj, indi->index) = genoValue;
            if (genoValue == -9) ++nmiss;
            else mean += genoValue;
        }
        delete[] bedLineIn;
        
        // fill missing values with the mean
        snpj->sampleSize = numKeptInds-nmiss;
        mean /= float(snpj->sampleSize);
        if (nmiss) {
            for (i=0; i<numKeptInds; ++i) {
                if (ZP(incj, i) == -9) ZP(incj, i) = mean;
            }
        }
        
        // compute allele frequency
        snpj->af = 0.5f*mean;
        snp2pq[incj] = 2.0f*snpj->af*(1.0f-snpj->af);
        
        if (snp2pq[incj]==0) throw ("Error: " + snpj->ID + " is a fixed SNP (MAF=0)!");
        
        // standardize genotypes
        //D[incj] = snp2pq[incj]*snpj->sampleSize;
        D[incj] = Gadget::calcVariance(ZP.row(incj))*numKeptInds;
        
        ZP.row(incj) = (ZP.row(incj).array() - mean)/sqrt(D[incj]);
        
        if (++incj == numSnpInRange) break;
    }
    
    fclose(in1);
    

    ZPZdiag = ZP.rowwise().squaredNorm(); 
    
    // then read in the bed file again to compute Z'Z
  
    MatrixXf denseZPZ;
    denseZPZ.setZero(numSnpInRange, numIncdSnps);
    VectorXf Zk(numKeptInds);
    D.setZero(numIncdSnps);

    FILE *in2 = fopen(bedFile.c_str(), "rb");
    fseek(in2, 3, SEEK_SET);
    unsigned long long skipk = 0;

    for (k = 0, inck = 0; k < numSnps; k++) {
        snpk = snpInfoVec[k];
        
        if (!snpk->included) {
            skipk += size;
            continue;
        }
        
        if (skipk) fseek(in2, skipk, SEEK_CUR);
        skipk = 0;
        
        char *bedLineIn = new char[size];
        fread(bedLineIn, sizeof(char), size, in2);
        
        mean = 0.0;
        nmiss = 0;
        
        for (i = 0; i < numInds; i++) {
            indi = indInfoVec[i];
            if (!indi->kept) continue;
            genoValue = bedToGeno[(bedLineIn[i>>2]>>((i&3)<<1))&3];
            Zk[indi->index] = genoValue;
            if (genoValue == -9) ++nmiss;   // missing genotype
            else mean += genoValue;
        }
        delete[] bedLineIn;
        
        // fill missing values with the mean
        snpk->sampleSize = numKeptInds-nmiss;
        mean /= float(snpk->sampleSize);
        if (nmiss) {
            for (i=0; i<numKeptInds; ++i) {
                if (Zk[i] == -9) Zk[i] = mean;
            }
        }
        
        // compute allele frequency
        snpk->af = 0.5f*mean;
        snp2pq[inck] = 2.0f*snpk->af*(1.0f-snpk->af);
        
        if (snp2pq[inck]==0) throw ("Error: " + snpk->ID + " is a fixed SNP (MAF=0)!");
        
        // standardize genotypes
        //D[inck] = snp2pq[inck]*snpk->sampleSize;
        D[inck] = Gadget::calcVariance(Zk)*numKeptInds;
        
        Zk = (Zk.array() - mean)/sqrt(D[inck]);
        
        denseZPZ.col(inck) = ZP * Zk;
        
        if(!(inck%1000)) cout << " read snp " << inck << "\r" << flush;
        
        ++inck;
    }

    fclose(in2);
    
    // ----------------------------------------------------
    // NEW - Calculate the mutation rate components 
    // ----------------------------------------------------
    // m is the number of individuals in the reference panel for each variant
    // Need to compute theta, which is related to the mutation rate
    VectorXf nmsumi(numIncdSnps);
    VectorXf thetai(numIncdSnps);
    VectorXf mi(numIncdSnps);
    VectorXf gmapi(numIncdSnps);
    VectorXf sdss(numIncdSnps);
    cout << "\nUsing genetic map sample size of " << genMapN << " please alter with --genmap-n if inappropriate." << endl;
    float m = genMapN;
    for (unsigned i=0; i<numIncdSnps; ++i) {
            SnpInfo *snp = incdSnpInfoVec[i];
            //mi[i] = (snp->sampleSize);
            int  n = 2.0 * m - 1.0;
            // Approximation to the harmonic series
            nmsumi[i] = log(n) + 0.5772156649 + 1.0 / (2.0 * n) - 1.0 / (12.0 * pow(n, 2.0)) + 1.0 / (120.0 * pow(n, 4.0));
            //cout << nmsumi[i] << endl;
            // Calculate theta
            //thetai[i] = (1.0 / nmsumi[i]) / (2.0 * (snp->sampleSize) + 1.0 / nmsumi[i]);
            thetai[i] = (1.0 / nmsumi[i]) / (2.0 * m + 1.0 / nmsumi[i]);
            //cout <<  thetai[i] << endl;
            // Pull out the standard deviation for each variant
            sdss[i] = sqrt(2.0 * (snp->af) * (1.0 - (snp->af)));
            // 
            gmapi[i] = snp->genPos;
    }
    long int nmsum;
    float theta;
     

    // Rescale Z to the covariance scale
    for (unsigned i=0; i<numIncdSnps; ++i) {
        denseZPZ.col(i) = (sdss[i] / 2.0) *  sdss.array() * denseZPZ.col(i).array();
    }
    // --------------------------------------------------------
    // Compute the shrinkage value and then shrink the elements
    // --------------------------------------------------------
    float mapdiffi; 
    float rho;
    float shrinkage;
    float Ne = effpopNE;
    cout << "\nUsing European effective population size Ne=" << Ne << " please alter with --ne if inappropriate." << endl;
    float cutoff = cutOff;
    for (unsigned i=0; i<numSnpInRange; ++i) {
        for (unsigned j=0; j<numIncdSnps; ++j) {
            mapdiffi = abs(gmapi[j] - gmapi[start+i]);
            rho = 4.0 * Ne * (mapdiffi / 100.0);
            shrinkage = exp(-rho / (2 * m)); 
            // if (i <=10 && j <= 10)
            // { 
            //     cout << "Snp " << i << " " << j << " Mapdiff " << mapdiffi << " shrinkage " << shrinkage << endl;
            //     cout << "Start position " << start  << " start plus i " << start + i << " gmap start plus i " << gmapi[start +i] << endl;
            //     // cout << gmapi[i] << " " << gmapi[j] << endl;
            //     // cout << mapdiffi << endl;
            //     // cout << shrinkage << endl;
            //     // cout << mi[i] << endl;
            // }
            if (shrinkage <= cutoff) {
                shrinkage = 0.0;
            }
            // Multiple each covariance matrix element with the shrinkage value
            denseZPZ(i, j) *= shrinkage;
            // Complete as SigHAat from Li and Stephens 2003
            denseZPZ(i, j) *= (1.0 - thetai[start+i]) * (1.0 - thetai[j]);
            // If it's the diagonal element add the extra term
            if (start+i == j)
            {
               denseZPZ(i, j) = denseZPZ(i, j) + 0.5f * thetai[start+i] * (1 - 0.5f * thetai[start+i]);
            }
            // Make the upper triangle equal to the lower triangle
            //denseZPZ(j, i) =  denseZPZ(i, j);
        }
    }
    // Now back to correlation
    for (unsigned i=0; i<numIncdSnps; ++i) {
        denseZPZ.col(i) = (2.0 / sdss[i]) *  (1.0 / sdss.array()) * denseZPZ.col(i).array();
    }
    // --------------------------------
    // Copy to vector of vectors format
    // --------------------------------
    ZPZ.resize(numSnpInRange);
    windStart.setZero(numSnpInRange);
    windSize.setConstant(numSnpInRange, numIncdSnps);
    // cout << "Num snp range " << numSnpInRange << endl;
    for (unsigned i=0; i<numSnpInRange; ++i) {
        SnpInfo *snp = incdSnpInfoVec[start+i];
        snp->windStart = 0;
        snp->windSize  = numIncdSnps;
        snp->windEnd   = numIncdSnps-1;
        ZPZ[i] = denseZPZ.row(i);
        snp->ldSum = denseZPZ.row(i).sum();
    }
        // cout << ZPZ[0] << endl;
    //}
    // else if (LDmatType == "band") {
    //     ZPZ.resize(numSnpInRange);
    //     if (windowWidth) {  // based on the given window width
    //         for (unsigned i=0; i<numSnpInRange; ++i) {
    //             SnpInfo *snp = incdSnpInfoVec[start+i];
    //             ZPZ[i] = denseZPZ.row(i).segment(snp->windStart, snp->windSize);
    //         }
    //     } else {  // based on the given LD threshold
    //         windStart.setZero(numSnpInRange);
    //         windSize.setZero(numSnpInRange);
    //         for (unsigned i=0; i<numSnpInRange; ++i) {
    //             SnpInfo *snp = incdSnpInfoVec[start+i];
    //             unsigned windEndi = numIncdSnps;
    //             for (unsigned j=0; j<numIncdSnps; ++j) {
    //                 if (abs(denseZPZ(i,j)) > LDthreshold) {
    //                     windStart[i] = snp->windStart = j;
    //                     break;
    //                 }
    //             }
    //             for (unsigned j=numIncdSnps; j>0; --j) {
    //                 if (abs(denseZPZ(i,j-1)) > LDthreshold) {
    //                     windEndi = j;
    //                     break;
    //                 }
    //             }
    //             windSize[i] = snp->windSize = windEndi - windStart[i];
    //             snp->windEnd = windEndi - 1;
    //             ZPZ[i].resize(windSize[i]);
    //             VectorXf::Map(&ZPZ[i][0], windSize[i]) = denseZPZ.row(i).segment(windStart[i], windSize[i]);
    //         }
    //     }
    // }
    // else if (LDmatType == "sparse") {
    //     ZPZsp.resize(numSnpInRange);
    //     windStart.setZero(numSnpInRange);
    //     windSize.setZero(numSnpInRange);
    //     if (LDthreshold) {
    //         for (unsigned i=0; i<numSnpInRange; ++i) {
    //             SnpInfo *snp = incdSnpInfoVec[start+i];
    //             for (unsigned j=0; j<numIncdSnps; ++j) {
    //                 if (abs(denseZPZ(i,j)) < LDthreshold) denseZPZ(i,j) = 0;
    //             }
    //             ZPZsp[i] = denseZPZ.row(i).sparseView();
    //             SparseVector<float>::InnerIterator it(ZPZsp[i]);
    //             windStart[i] = snp->windStart = it.index();
    //             windSize[i] = snp->windSize = ZPZsp[i].nonZeros();
    //             for (; it; ++it) snp->windEnd = it.index();
    //         }
    //     } else {
    //         for (unsigned i=0; i<numSnpInRange; ++i) {
    //             SnpInfo *snp = incdSnpInfoVec[start+i];
    //             for (unsigned j=0; j<numIncdSnps; ++j) {
    //                 if (denseZPZ(i,j)*denseZPZ(i,j)*snp->sampleSize < chisqThreshold) denseZPZ(i,j) = 0;
    //             }
    //             ZPZsp[i] = denseZPZ.row(i).sparseView();
    //             SparseVector<float>::InnerIterator it(ZPZsp[i]);
    //             windStart[i] = snp->windStart = it.index();
    //             windSize[i] = snp->windSize = ZPZsp[i].nonZeros();
    //             for (; it; ++it) snp->windEnd = it.index();
    //         }
    //     }
    // }
    
    denseZPZ.resize(0,0);
    
    //    cout << denseZPZ.block(0,0,10,10) << endl;
    //    cout << windStart.transpose() << endl;
    //    cout << windSize.transpose() << endl;
    
    
    timer.getTime();
    
    
    cout << endl;
    displayAverageWindowSize(windSize);
    cout << "LD matrix diagonal mean " << ZPZdiag.mean() << " variance " << Gadget::calcVariance(ZPZdiag) << "." << endl;
    if (ZPZdiag.mean() < 0.8 || ZPZdiag.mean() > 1.2) cout << "ERROR: The mean of LD matrix diagonal values is expected to be close to one. Something is wrong with the LD matrix!" << endl;
cout << "Genotype data for " << numKeptInds << " individuals and " << numSnpInRange << " SNPs are included from [" + bedFile + "]." << endl;
    cout << "Build of LD matrix completed (time used: " << timer.format(timer.getElapse()) << ")." << endl;
    
    
    vector<SnpInfo*> snpVecTmp(numSnpInRange);
    for (unsigned i=0; i<numSnpInRange; ++i) {
        snpVecTmp[i] = incdSnpInfoVec[start+i];
    }
    incdSnpInfoVec = snpVecTmp;
    numIncdSnps = numSnpInRange;
    string outfilename = filename;
    if (!snpRange.empty()) outfilename += ".snp" + snpRange;
    outputLDmatrix(LDmatType, outfilename, writeLdmTxt);
}

// =============================================================================================
// Make shrunk matrix end
// =============================================================================================

void Data::buildSparseMME(const bool sampleOverlap, const bool noscale){
    VectorXf Dref = snp2pq*numKeptInds;
    snp2pq.resize(numIncdSnps);
    D.resize(numIncdSnps);
//    ZPZdiag.resize(numIncdSnps);
    ZPy.resize(numIncdSnps);
    b.resize(numIncdSnps);
    n.resize(numIncdSnps);
    se.resize(numIncdSnps);
    tss.resize(numIncdSnps);
    SnpInfo *snp;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snp = incdSnpInfoVec[i];
        snp->af = snp->gwas_af;
        snp2pq[i] = snp->twopq = 2.0f*snp->gwas_af*(1.0f-snp->gwas_af);
        if(snp2pq[i]==0) cout << "Error: SNP " << snp->ID << " af " << snp->af << " has 2pq = 0." << endl;
        D[i] = snp2pq[i]*snp->gwas_n;
        b[i] = snp->gwas_b;
        n[i] = snp->gwas_n;
        se[i]= snp->gwas_se;
        tss[i] = D[i]*(n[i]*se[i]*se[i] + b[i]*b[i]);
//        D[i] = 1.0/(se[i]*se[i]+b[i]*b[i]/snp->gwas_n);  // NEW!
//        snp2pq[i] = snp->twopq = D[i]/snp->gwas_n;       // NEW!
    }
    //b.array() -= b.mean();  // DO NOT CENTER b
    
    // estimate phenotypic variance based on the input allele frequencies in GWAS
    //ypy = (D.array()*(n.array()*se.array().square()+b.array().square())).mean();
    VectorXf ypySrt = D.array()*(n.array()*se.array().square()+b.array().square());
    VectorXf varpSrt = ypySrt.array()/n.array();
    std::sort(ypySrt.data(), ypySrt.data() + ypySrt.size());
    std::sort(varpSrt.data(), varpSrt.data() + varpSrt.size());
    ypy = ypySrt[ypySrt.size()/2];  // median
    varPhenotypic = varpSrt[varpSrt.size()/2];
    
    //numKeptInds = n.mean();
    
    VectorXf nSrt = n;
    std::sort(nSrt.data(), nSrt.data() + nSrt.size());
    numKeptInds = nSrt[nSrt.size()/2]; // median
    
        // NEW
        // compute D and snp2pq based on n, se and b, assuming varp = 1
        // these quantities are used in sbayes, as they are more reliable than input allele frequencies
        for (unsigned i=0; i<numIncdSnps; ++i) {
            snp = incdSnpInfoVec[i];
            D[i] = varPhenotypic/(se[i]*se[i]+b[i]*b[i]/snp->gwas_n);  // NEW!
            snp2pq[i] = snp->twopq = D[i]/snp->gwas_n;       // NEW!
            tss[i] = D[i]*(n[i]*se[i]*se[i] + b[i]*b[i]);
            // Need to adjust R and C models X'X matrix depending scale of genotypes or not
            if (noscale == true) {
                D[i] = snp2pq[i]*snp->gwas_n;
            } else {
                D[i] = snp->gwas_n;
            }
        }
        //ypy = numKeptInds;
        // NEW END
    
    if (ZPZ.size() || ZPZsp.size()) {
        if (sparseLDM == true) {
            for (unsigned i=0; i<numIncdSnps; ++i) {
                snp = incdSnpInfoVec[i];
                //cout << i << " " << ZPZsp[i].nonZeros() << " " << D.size() << endl;
                for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                    //cout << it.index() << " ";
                    it.valueRef() *= sqrt(D[i]*D[it.index()]);
                }
            }
        } else {
            for (unsigned i=0; i<numIncdSnps; ++i) {
                snp = incdSnpInfoVec[i];
                for (unsigned j=0; j<snp->windSize; ++j) {
                    ZPZ[i][j] *= sqrt(D[i]*D[snp->windStart+j]);
                }
            }
        }

        // sum of sampling variance of LD for each SNP with all other SNPs
        // for significant LD, the sampling variance is proportional to the (ratio of ref and gwas n) + 1
        // for insignificant LD, the sampling variance is 1 over gwas n
        LDsamplVar.resize(numIncdSnps);
        LDscore.resize(numIncdSnps);
        for (unsigned i=0; i<numIncdSnps; ++i) {
            snp = incdSnpInfoVec[i];
            if (sampleOverlap) {
                LDsamplVar[i] = 0;
            } else {
                LDsamplVar[i]  = (snp->gwas_n + snp->sampleSize)/float(numIncdSnps)*snp->ldSamplVar;
            }
            LDsamplVar[i] += (numIncdSnps - snp->numNonZeroLD)/float(numIncdSnps);
            //if (sampleOverlap) LDsamplVar[i] = 0;
            LDscore[i] = snp->ldsc; //*snp->gwas_n;
        }

        ZPZdiag.array() *= D.array();
    }
    else {
        Dratio = D.array()/Dref.array();
        ZPZdiag.array() *= Dratio.array();
        DratioSqrt = Dratio.array().sqrt();
        for (unsigned i=0; i<numIncdSnps; ++i) {
            Z.col(i) *= DratioSqrt[i];
        }
    }
    
    
    if (noscale) {
        ZPy = ZPZdiag.cwiseProduct(b);
    } else {
        ZPy = ZPZdiag.cwiseProduct(b).cwiseProduct(snp2pq.array().sqrt().matrix());
    }
    chisq = ZPy.cwiseProduct(b);
    
//    cout << "!!!!!!!!" << endl << ZPZdiag.transpose() << endl << endl;
    
//        ofstream out("ldsc.txt");
//        for (unsigned i=0; i<numIncdSnps; ++i) {
//            snp = incdSnpInfoVec[i];
//            out << chisq[i] << "\t" << LDscore[i] << "\t" << LDsamplVar[i] << "\t" << n[i] << "\t" << n[i]*(numIncdSnps+snp->windSize)/float(numIncdSnps) << endl;
//        }
//        out.close();

    //    cout << "ZPZdiag " << ZPZdiag.transpose() << endl;
    //    cout << "ZPZ.back() " << ZPZ.back().transpose() << endl;
    //    cout << "ZPZ.front() " << ZPZ.front().transpose() << endl;
    //    cout << "ZPy " << ZPy.head(100).transpose() << endl;
    //    cout << "b.mean() " << b.mean() << endl;
    
    // estimate ypy
    // ypy = tss.mean();
    // ypy = (D.array()*(n.array()*se.array().square()+b.array().square())).mean();
    numKeptInds = n.mean();
    
    //cout << ZPZ.size() << " " << ZPy.size() << " " << ypy << endl;
    //    cout << ZPy << endl;
    //    for (unsigned i=0; i<numIncdSnps; ++i) {
    //        cout << D[i] << "\t" << ZPZdiag[i] << endl;
    //    }
    //
    //    cout << ZPZ << endl;
    
    // no fixed effects
//    numFixedEffects = 0;
//    fixedEffectNames.resize(0);
//    XPX.resize(0,0);
//    ZPX.resize(0,0);
//    XPy.resize(0);
    numFixedEffects = 1;
    fixedEffectNames.resize(1);
    fixedEffectNames[0] = "Intercept";
    ZPX.resize(numIncdSnps,1);
    if (ZPZ.size()) {
        for (unsigned i=0; i<numIncdSnps; ++i) {
            ZPX(i,0) = ZPZ[i].sum();
        }
    } else if (ZPZsp.size()) {
        for (unsigned i=0; i<numIncdSnps; ++i) {
            ZPX(i,0) = ZPZsp[i].sum();
        }
    } else {
        throw("Error: either dense or sparse ldm does not exist!");
    }
    XPX.resize(1,1);
    XPX << ZPX.col(0).sum();
    XPXdiag.resize(1);
    XPXdiag << XPX(0,0);
    XPy.resize(1,1);
    XPy << ZPy.sum();
    
//    cout << "ZPZ" << endl;
//    for (unsigned i=0; i<numIncdSnps; ++i) {
//        cout << ZPZ[i].transpose() << endl;
//    }
//    
//    cout << "ZPZdiag" << endl;
//    cout << ZPZdiag.transpose() << endl;
//    
//    cout << "D" << endl << D.transpose() << endl << endl;
//    
//    cout << "n" << endl << n.transpose() << endl << endl;
//    
//    cout << "2pq" << endl << snp2pq.transpose() << endl << endl;
//    
//    cout << "b" << endl << b.transpose() << endl << endl;
//
//    cout << "ZPy" << endl;
//    cout << ZPy.transpose() << endl;
    
    // data summary
    cout << "\nData summary:" << endl;
    cout << boost::format("%40s %8s %8s\n") %"" %"mean" %"sd";
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP Phenotypic variance" %Gadget::calcMean(varpSrt) %sqrt(Gadget::calcVariance(varpSrt));
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP heterozygosity" %Gadget::calcMean(snp2pq) %sqrt(Gadget::calcVariance(snp2pq));
    cout << boost::format("%40s %8.0f %8.0f\n") %"GWAS SNP sample size" %Gadget::calcMean(n) %sqrt(Gadget::calcVariance(n));
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP effect" %Gadget::calcMean(b) %sqrt(Gadget::calcVariance(b));
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP SE" %Gadget::calcMean(se) %sqrt(Gadget::calcVariance(se));
    cout << boost::format("%40s %8.3f %8.3f\n") %"MME left-hand-side diagonals" %Gadget::calcMean(ZPZdiag) %sqrt(Gadget::calcVariance(ZPZdiag));
    cout << boost::format("%40s %8.3f %8.3f\n") %"MME right-hand-side" %Gadget::calcMean(ZPy) %sqrt(Gadget::calcVariance(ZPy));
    cout << boost::format("%40s %8.3f %8.3f\n") %"LD sampling variance" %Gadget::calcMean(LDsamplVar) %sqrt(Gadget::calcVariance(LDsamplVar));
    cout << boost::format("%40s %8.3f %8.3f\n") %"LD score" %Gadget::calcMean(LDscore) %sqrt(Gadget::calcVariance(LDscore));
//    cout << "\n  Median of per-SNP phenotypic variance: " << varp << endl;
    
//    ofstream out("tmp.txt");
//    out << "refZPZdiag\t gwasZPZdiag\t b\t ZPy\t refsnp2pq\t gwassnp2pq n" << endl;
//    for (unsigned i=0; i<numIncdSnps; ++i) {
//        out << refZPZdiag[i] << "\t" << ZPZdiag[i] << "\t" << b[i] << "\t" << ZPy[i] << "\t" << refsnp2pq[i] << "\t" << snp2pq[i] << "\t" << n[i] << endl;
//    }
//    out.close();
    
    if (numAnnos) setAnnoInfoVec();

}

void Data::setAnnoInfoVec() {
    bool print = numAnnos > 50 ? false : true;
    if (print) cout << "\nAnnotation info:" << endl;
    numSnpAnnoVec.resize(numAnnos);
    for (unsigned i=0; i<numAnnos; ++i) {
        AnnoInfo *anno = annoInfoVec[i];
        anno->idx = i;
        anno->getSnpInfo();
        anno->fraction = float(anno->size)/float(numIncdSnps);
        if (print) anno->print();
        numSnpAnnoVec[i] = anno->size;
    }
    if (print) cout << endl;
    
    // set up snp-annot pair names
    snpAnnoPairNames.clear();
    numAnnoPerSnpVec.resize(numIncdSnps);
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        numAnnoPerSnpVec[i] = snp->numAnnos;
        for (unsigned j=0; j<snp->numAnnos; ++j) {
            string snpAnnoPair = snp->ID + "\t" + snp->annoVec[j]->label;
            snpAnnoPairNames.push_back(snpAnnoPair);
        }
    }
    
    // set up annotation by annotation matrix APA
    annoMat.setZero(numIncdSnps, numAnnos);
//    annoMat.setZero(numIncdSnps, numAnnos+1); // first is the intercept
//    annoMat.col(0) = VectorXf::Ones(numIncdSnps);
    APA.setZero(numAnnos, numAnnos);
//    APA.setZero(numAnnos+1, numAnnos+1);
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        for (unsigned j=0; j<snp->numAnnos; ++j) {
            unsigned annoIdx = snp->annoVec[j]->idx;
//            annoMat(i,annoIdx+1) = 1;
            annoMat(i,annoIdx) = 1;
        }
    }
        
    // TMP
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        annoMat.row(i) = snp->annoValues;
    }
    annoSD.setZero(numAnnos);
    for (unsigned i=0; i<numAnnos; ++i) {
        AnnoInfo *anno = annoInfoVec[i];
        if (annoMat.col(i).sum() == anno->size) annoSD[i] = 1; // binary annotation
        else annoSD[i] = sqrt(Gadget::calcVariance(annoMat.col(i))); // quantitative annotation
    }
    
    
    annoMean.setZero(numAnnos);
    for (unsigned i=0; i<numAnnos; ++i) {
        annoMean[i] = annoMat.col(i).mean();
        //if (i) annoMat.col(i).array() -= annoMean[i];  // center the annotation matrix
    }
    APA = annoMat.transpose()*annoMat;    
}


void Data::outputSnpEffectSamples(const SpMat &snpEffects, const unsigned burnin, const unsigned outputFreq, const string&snpResFile, const string &filename) const {
    cout << "writing SNP effect samples into " << filename << endl;
    unsigned nrow = snpEffects.rows();
    vector<string> snpName;
    vector<float> sample;

    ifstream in(snpResFile.c_str());
    if (!in) throw ("Error: can not open the snpRes file [" + snpResFile + "] to read.");

    Gadget::Tokenizer colData;
    string inputStr;
    string sep(" \t");
    string id;
    unsigned line=0;
    while (getline(in,inputStr)) {
        ++line;
        if (line==1) continue;
        colData.getTokens(inputStr, sep);
        snpName.push_back(colData[1]);
    }
    in.close();
    long numSnps = snpName.size();
    
    ofstream out(filename.c_str());
    out << boost::format("%6s %20s %8s\n")
    % "Iteration"
    % "Name"
    % "Sample";
    
    cout << "Size of mcmc samples " << snpEffects.rows() << " " << snpEffects.cols() << endl;
    
    unsigned idx=0;
    for (unsigned iter=0; iter<nrow; ++iter) {
        if (iter < burnin) continue;
        if (!(iter % outputFreq)) {
            ++idx;
            for (unsigned j=0; j<numSnps; ++j) {
                if (snpEffects.coeff(iter, j)) {
                    out << boost::format("%6s %20s %8s\n")
                    % idx
                    % snpName[j]
                    % snpEffects.coeff(iter, j);
                }
            }
        }
    }
    
    out.close();
}

void Data::directPruneLDmatrix(const string &ldmatrixFile, const string &outLDmatType, const float chisqThreshold, const string &title, const bool writeLdmTxt){
    readLDmatrixInfoFile(ldmatrixFile + ".info");
    Gadget::Tokenizer token;
    token.getTokens(ldmatrixFile, ".");
    string ldmType = token[token.size()-1];
    if (ldmType != "full") {
        throw("Error: --direct-prune only prunes a full LD matrix to a sparse matrix!");
    }
    cout << "Direct pruning " + ldmType + " LD matrix from [" + ldmatrixFile + "]..." << endl;

    FILE *in = fopen((ldmatrixFile+".bin").c_str(), "rb");
    if (!in) {
        throw("Error: cannot open LD matrix file " + ldmatrixFile);
    }
    SnpInfo *snp;
    VectorXf vec;
    float rsq = 0.0;

    string outfilename = title + ".ldm." + outLDmatType;
    string outfile1 = outfilename + ".info";
    string outfile2 = outfilename + ".bin";
    ofstream out1(outfile1.c_str());
    FILE *out2 = fopen(outfile2.c_str(), "wb");
    ofstream out3;
    string outfile3;
    if (writeLdmTxt) {
        outfile3 = outfilename + ".txt";
        out3.open(outfile3.c_str());
    }
    out1 << boost::format("%6s %15s %10s %15s %6s %6s %12s %10s %10s %10s %10s %15s %10s %12s %12s\n")
    % "Chrom"
    % "ID"
    % "GenPos"
    % "PhysPos"
    % "A1"
    % "A2"
    % "A1Freq"
    % "Index"
    % "WindStart"
    % "WindEnd"
    % "WindSize"
    % "WindWidth"
    % "N"
    % "SamplVar"
    % "LDsum";

    SnpInfo *windStart, *windEnd;
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        snp->ldSamplVar = 0.0;
        snp->ldSum = 0.0;
        
        vec.resize(snp->windSize);
        fread(vec.data(), sizeof(float), snp->windSize, in);
        
//        if (!(i%100)) cout << " SNP " << i << "\r";
        
        for (unsigned j=0; j<snp->windSize; ++j) {
            rsq = vec[j]*vec[j];
            if (i!=j && rsq*snp->sampleSize < chisqThreshold) {
                vec[j] = 0;
            } else {
                snp->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snp->sampleSize;
                snp->ldSum += vec[j];
            }
        }
        
        SparseVector<float> spvec = vec.sparseView();
        SparseVector<float>::InnerIterator it(spvec);
        snp->windStart = it.index();
        snp->windSize = spvec.nonZeros();
        for (; it; ++it) snp->windEnd = it.index();
        
        out1 << boost::format("%6s %15s %10s %15s %6s %6s %12f %10s %10s %10s %10s %15s %10s %12.6f %12.6f\n")
        % snp->chrom
        % snp->ID
        % snp->genPos
        % snp->physPos
        % snp->a1
        % snp->a2
        % snp->af
        % snp->index
        % snp->windStart
        % snp->windEnd
        % snp->windSize
        % -9
        % snp->sampleSize
        % snp->ldSamplVar
        % snp->ldSum;
        fwrite(spvec.innerIndexPtr(), sizeof(unsigned), snp->windSize, out2);
        fwrite(spvec.valuePtr(), sizeof(float), snp->windSize, out2);
        if (writeLdmTxt) out3 << vec.transpose() << endl;
    }
    out1.close();
    fclose(out2);
    
    cout << "Written the SNP info into file [" << outfile1 << "]." << endl;
    cout << "Written the LD matrix into file [" << outfile2 << "]." << endl;
    
    if (writeLdmTxt) {
        out3.close();
        cout << "Written the LD matrix into text file [" << outfile3 << "]." << endl;
    }
}

void Data::addLDmatrixInfo(const string &ldmatrixFile) {
    ifstream in(ldmatrixFile.c_str());
    if (!in) throw ("Error: can not open the file [" + ldmatrixFile + "] to read.");
    cout << "Reading SNP info from [" + ldmatrixFile + "]." << endl;
    string header;
    string id, allele1, allele2;
    unsigned chr, physPos;
    float genPos, af, ldSamplVar, ldSum;
    unsigned idx, windStart, windEnd, windSize, windWidth;
    long sampleSize;
    bool skeleton;
    getline(in, header);
    Gadget::Tokenizer token;
    token.getTokens(header, " ");
    map<string, SnpInfo*>::iterator it, end = snpInfoMap.end();
    unsigned count = 0;
    if (token.back() == "Skeleton") {
        while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2 >> af >> idx >> windStart >> windEnd >> windSize >> windWidth >> sampleSize >> ldSamplVar >> ldSum >> skeleton) {
            it = snpInfoMap.find(id);
            if (it == end) throw ("Error: SNP " + id + " not found in .bim file.");
            SnpInfo *snp = it->second;
            if (snp->index != idx) throw ("Error: The index of SNP " + id + " is different from that in .bim file.");
            snp->windStart = windStart;
            snp->windEnd = windEnd;
            snp->windSize = windSize;
            snp->sampleSize = sampleSize;
            snp->ldSamplVar = ldSamplVar;
            snp->ldSum = ldSum;
            snp->skeleton = skeleton;
            ++count;
        }
    }
    else if (token.back() == "LDsum") {
        while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2 >> af >> idx >> windStart >> windEnd >> windSize >> windWidth >> sampleSize >> ldSamplVar >> ldSum) {
            it = snpInfoMap.find(id);
            if (it == end) throw ("Error: SNP " + id + " not found in .bim file.");
            SnpInfo *snp = it->second;
            if (snp->index != idx) throw ("Error: The index of SNP " + id + " is different from that in .bim file.");
            snp->windStart = windStart;
            snp->windEnd = windEnd;
            snp->windSize = windSize;
            snp->sampleSize = sampleSize;
            snp->ldSamplVar = ldSamplVar;
            snp->ldSum = ldSum;
            ++count;
        }
    }
    else {
        while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2 >> af >> idx >> windStart >> windEnd >> windSize >> windWidth >> sampleSize >> ldSamplVar) {
            it = snpInfoMap.find(id);
            if (it == end) throw ("Error: SNP " + id + " not found in .bim file.");
            SnpInfo *snp = it->second;
            if (snp->index != idx) throw ("Error: The index of SNP " + id + " is different from that in .bim file.");
            snp->windStart = windStart;
            snp->windEnd = windEnd;
            snp->windSize = windSize;
            snp->sampleSize = sampleSize;
            snp->ldSamplVar = ldSamplVar;
            ++count;
        }
    }
    in.close();
    cout << count << " SNPs to be included from [" + ldmatrixFile + "]." << endl;

}

void Data::jackknifeLDmatrix(const string &ldmatrixFile, const string &outLDmatType, const string &title, const bool writeLdmTxt) {
    // Jackknife estimate of correlation and sampling variance
    Gadget::Tokenizer token;
    token.getTokens(ldmatrixFile, ".");
    if (token[token.size()-1] != "sparse") {
        throw("Error: --jackknife method only works for sparse LD matrix at the moment!");
    }

    addLDmatrixInfo(ldmatrixFile + ".info");

    FILE *in = fopen((ldmatrixFile + ".bin").c_str(), "rb");
    if (!in) {
        throw("Error: cannot open LD matrix file " + ldmatrixFile + ".bin");
    }

    if (numIncdSnps == 0) throw ("Error: No SNP is retained for analysis.");
    
    cout << "Reading and jackknifing sparse LD matrix from [" + ldmatrixFile + ".bin]..." << endl;
    
    SnpInfo *snpi, *snpj;
    SnpInfo *windStart, *windEnd;
    SparseVector<float> ZPZspvec;
    VectorXf ZiCwiseZj(numKeptInds);
    VectorXf ones;
    ones.setOnes(numKeptInds);
    float numJackknife = numKeptInds-1;
    float samplVarEmp;
    
    string outfilename = title + ".ldm." + outLDmatType;
    string outfile = outfilename + ".info";
    ofstream out(outfile.c_str());
    out << boost::format("%6s %15s %10s %15s %6s %6s %12s %10s %10s %10s %10s %15s %10s %12s %12s\n")
    % "Chrom"
    % "ID"
    % "GenPos"
    % "PhysPos"
    % "A1"
    % "A2"
    % "A1Freq"
    % "Index"
    % "WindStart"
    % "WindEnd"
    % "WindSize"
    % "WindWidth"
    % "N"
    % "SamplVar"
    % "LDsum";
    
    
    // standardize genotype matrix
    for (unsigned i=0; i<numIncdSnps; ++i) {
        Z.col(i) /= sqrt(numKeptInds*snp2pq[i]);
    }

    
    for (unsigned i=0, inci=0; i<numSnps; i++) {
        snpi = snpInfoVec[i];
        
        unsigned d[snpi->windSize];
        float v[snpi->windSize];
        
        if (!snpi->included) {
            fseek(in, sizeof(d), SEEK_CUR);
            fseek(in, sizeof(v), SEEK_CUR);
            continue;
        }
        
        fread(d, sizeof(d), 1, in);
        fread(v, sizeof(v), 1, in);
        
        ZPZspvec.resize(snpi->windSize);
//        snpi->ldSamplVar = 0.0;
        snpi->ldSum = 0.0;
        
        for (unsigned j=0; j<snpi->windSize; ++j) {
            snpj = snpInfoVec[d[j]];
            if (snpj->included) {
                ZiCwiseZj  = v[j]*ones - Z.col(snpi->index).cwiseProduct(Z.col(snpj->index));
                ZiCwiseZj *= numKeptInds/numJackknife;
                samplVarEmp = (ZiCwiseZj.array() - ZiCwiseZj.mean()).square().sum() * numJackknife/numKeptInds;
                snpi->ldSum += samplVarEmp;
            }
        }

        windStart = incdSnpInfoVec[snpi->windStart];
        windEnd = incdSnpInfoVec[snpi->windEnd];
        out << boost::format("%6s %15s %10s %15s %6s %6s %12f %10s %10s %10s %10s %15s %10s %12.6f %12.6f\n")
        % snpi->chrom
        % snpi->ID
        % snpi->genPos
        % snpi->physPos
        % snpi->a1
        % snpi->a2
        % snpi->af
        % snpi->index
        % snpi->windStart
        % snpi->windEnd
        % snpi->windSize
        % (windStart->chrom == windEnd->chrom ? windEnd->physPos - windStart->physPos : windStart->chrom-windEnd->chrom)
        % snpi->sampleSize
        % snpi->ldSamplVar
        % snpi->ldSum;

        if (++inci == numIncdSnps) break;
    }

    out.close();
    cout << "Written the SNP info into file [" << outfile << "]." << endl;

}


void Data::readAnnotationFile(const string &annoFile, const bool transpose, const bool allowMultiAnno) {
    ifstream in(annoFile.c_str());
    if (!in) throw ("Error: can not open the annotation file [" + annoFile + "] to read.");
    cout << "Reading SNP annotation from [" + annoFile + "]." << endl;
    
    map<string, SnpInfo*>::iterator it, end=snpInfoMap.end();
    SnpInfo *snp = NULL;
    Gadget::Tokenizer header;
    Gadget::Tokenizer colData;
    string inputStr;
    string sep(" \t");
    
    getline(in,inputStr);
    header.getTokens(inputStr, sep);
    long size=header.size();
    long line=0;

    // transpose the data matrix if columns are SNPs rather than annotation categories
    if (transpose) {
        vector<vector<string> > datmat;
        datmat.push_back(header);
        while (getline(in,inputStr)) {
            colData.getTokens(inputStr, sep);
            if (colData.size() != size) {
                throw("Error: the annotation file does not have consistent number of columns!");
            }
            datmat.push_back(colData);
            datmat[++line].resize(size);
            for (unsigned j=0; j<size; ++j) {
                datmat[line][j] = colData[j];
            }
        }
        line = datmat.size();
        string outfile = Gadget::getFileName(annoFile) + ".transpose" + Gadget::getFileSuffix(annoFile);
        ofstream out(outfile.c_str());
        cout << "Writing the tranposed SNP annotation to [" + outfile + "]." << endl;
        for (unsigned j=0; j<size; ++j) {
            for (unsigned i=0; i<line; ++i) {
                if (i==0) out << boost::format("%16s ") %datmat[i][j];
                else out << boost::format("%8s ") %datmat[i][j];
            }
            out << "\n";
        }
        readAnnotationFile(outfile, false, allowMultiAnno);
        return;
    }
    ///////////////////////////
    
    annoInfoVec.clear();
    annoNames.clear();
    for (unsigned i=1; i<size; ++i) {
        AnnoInfo *anno = new AnnoInfo(i-1, header[i]);
        annoInfoVec.push_back(anno);
        annoNames.push_back(header[i]);
    }
    numAnnos = annoInfoVec.size();

    string id;
    while (getline(in,inputStr)) {
        colData.getTokens(inputStr, sep);
        id = colData[0];
        it = snpInfoMap.find(id);
        if (it != end) {
            snp = it->second;
            snp->annoValues.setZero(size-1);
            for (unsigned j=1; j<size; ++j) {
                if (atof(colData[j].c_str())) {
                    snp->annoVec.push_back(annoInfoVec[j-1]);
                    snp->annoValues[j-1] = atof(colData[j].c_str());  // NEW LINE
                    snp->numAnnos++;
                }
            }
            ++line;
        }
    }
    in.close();
    
    unsigned numMultiAnno = 0;
    
    if (allowMultiAnno) {
        for (unsigned i=0; i<numSnps; ++i) {
            snp = snpInfoVec[i];
            if (!snp->included) continue;
            if (!snp->numAnnos) {
                snp->included = false;   // remove SNPs without any annotation from the model
            } else {
                for (unsigned j=0; j<snp->numAnnos; ++j) {
                    AnnoInfo *anno = snp->annoVec[j];
                    anno->memberSnpVec.push_back(snp);
                    anno->size++;
                }
                if (snp->numAnnos > 1) ++numMultiAnno;
            }
        }
    }
    else {
        for (unsigned i=0; i<numSnps; ++i) {
            snp = snpInfoVec[i];
            if (!snp->included) continue;
            if (!snp->numAnnos || snp->numAnnos>1) {
                snp->included = false;
            } else {
                for (unsigned j=0; j<snp->numAnnos; ++j) {
                    AnnoInfo *anno = snp->annoVec[j];
                    anno->memberSnpVec.push_back(snp);
                    anno->size++;
                }
                if (snp->annoVec.size() > 1) ++numMultiAnno;
            }
        }
    }
    
    cout << line << " matched SNPs in the annotation file (" << numAnnos << " annotations and " << numMultiAnno << " SNPs have more than one annotation)." << endl;
}

void Data::makeWindowAnno(const string &annoFile, const float windowWidth){
    cout << "Making window annotations ..." << endl;
    VectorXf annoSum;
    annoSum.setZero(numAnnos);
    unsigned currChr = 0;
    unsigned windStart = 0;
    unsigned windEnd = 0;
    SnpInfo *snpi;
    SnpInfo *snpj;
    float bound = 0.5*windowWidth;
    MatrixXf windAnnoMat(numIncdSnps, numAnnos);
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snpi = incdSnpInfoVec[i];
        if (snpi->chrom != currChr) {
            currChr = snpi->chrom;
            windStart = i;
            windEnd = i;
            annoSum.setZero(numAnnos);
        }
        for (unsigned j=windStart; j<numIncdSnps; ++j) {
            snpj = incdSnpInfoVec[j];
            if (snpi->physPos - snpj->physPos > bound) {
                annoSum -= annoMat.row(j);
                windStart++;
            } else {
                break;
            }
        }
        for (unsigned j=windEnd; j<numIncdSnps; ++j) {
            snpj = incdSnpInfoVec[j];
            if (snpj->physPos - snpi->physPos < bound) {
                annoSum += annoMat.row(j);
                windEnd++;
            } else {
                break;
            }
        }
        windAnnoMat.row(i) = annoSum/float(windEnd - windStart);
        
//        if (snpi->ID == "rs7588213") {
//            cout << "chrom " << snpi->chrom << endl;
//            cout << "prev_chrom " << incdSnpInfoVec[i-1]->chrom << endl;
//            cout << "windStart " << windStart << endl;
//            cout << "windEnd " << windEnd << endl;
//            cout << "BP " << snpi->physPos << " windStartBP " << incdSnpInfoVec[windStart]->physPos << " windEndBP " << incdSnpInfoVec[windEnd]->physPos << endl;
//            cout << "annoSum " << annoSum << endl;
//        }
    }
    
    string outfile = Gadget::getFileName(annoFile) + ".windowAnno" + Gadget::getFileSuffix(annoFile);
    ofstream out(outfile.c_str());
    cout << "Writing the window annotation to [" + outfile + "]." << endl;
    out << boost::format("%16s ") %"SNP";
    for (unsigned k=0; k<numAnnos; ++k) {
        out << boost::format("%16s ") % annoNames[k];
    }
    out << endl;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        out << boost::format("%16s ") %snp->ID;
        out << windAnnoMat.row(i) << endl;
    }
    throw ("Finished!");
}

void Data::readAnnotationFileFormat2(const string &continuousAnnoFile, const unsigned flank, const string &eQTLFile) {
    ifstream in;
    in.open(continuousAnnoFile.c_str());
    if (!in) throw ("Error: can not open the annotation file [" + continuousAnnoFile + "] to read.");
    cout << "Reading SNP annotation from [" + continuousAnnoFile + "]." << endl;
    
    SnpInfo *snp = NULL;
    AnnoInfo *anno = NULL;
    unsigned chrom;
    unsigned startBP;
    unsigned endBP;
    string name;

    unsigned line = 0;
    map<unsigned, vector<AnnoInfo*> > annoChrMap;
    map<string, AnnoInfo*> annoInfoMap;
    vector<AnnoInfo*> annoVec;
    while (in >> chrom >> startBP >> endBP >> name) {
        anno = new AnnoInfo(line++, name);
        anno->chrom = chrom;
        anno->startBP = startBP;
        anno->endBP = endBP;
        annoChrMap[chrom].push_back(anno);
        annoInfoMap[name] = anno;
        annoVec.push_back(anno);
    }
    in.close();
    in.clear();
    
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        if (!snp->included) continue;
        vector<AnnoInfo*> &annovec = annoChrMap[snp->chrom];
        for (unsigned j=0; j<annovec.size(); ++j) {
            anno = annoChrMap[snp->chrom][j];
            if (snp->physPos >= (anno->startBP - flank) && snp->physPos <= (anno->endBP + flank)) {
                snp->annoVec.push_back(anno);
                snp->annoMap[anno->idx] = anno;
                snp->numAnnos++;
                anno->memberSnpVec.push_back(snp);
                anno->memberSnpMap[snp->index] = snp;
                anno->size++;
            }
        }
    }
    
    if (!eQTLFile.empty()) {
        ifstream in;
        in.open(eQTLFile.c_str());
        if (!in) throw ("Error: can not open the eQTL file [" + eQTLFile + "] to read.");
        cout << "Reading eQTL info from [" + eQTLFile + "]." << endl;
        Gadget::Tokenizer header;
        Gadget::Tokenizer colData;
        string inputStr;
        string sep(" \t");
        getline(in, inputStr);
        header.getTokens(inputStr, sep);
        unsigned snpIDidx = header.getIndex("SNP");
        unsigned chrIDidx = header.getIndex("Chr");
        unsigned geneIDidx = header.getIndex("Gene");
        string snpID, geneID;
        unsigned chrID;
        map<string, SnpInfo*>::iterator itsnp, endsnp = snpInfoMap.end();
        map<string, AnnoInfo*>::iterator itgene, endgene = annoInfoMap.end();
        while (getline(in, inputStr)) {
            colData.getTokens(inputStr, sep);
            snpID = colData[snpIDidx];
            geneID = colData[geneIDidx];
            chrID = atoi(colData[chrIDidx].c_str());
            itsnp = snpInfoMap.find(snpID);
            if (itsnp != endsnp) {
                snp = itsnp->second;
                if (!snp->included) continue;
                itgene = annoInfoMap.find(geneID);
                if (itgene != endgene) {
                    anno = itgene->second;
                    if (snp->annoMap.insert(pair<int, AnnoInfo*>(anno->idx, anno)).second) {
                        snp->annoVec.push_back(anno);
                        snp->numAnnos++;
                    }
                    if (anno->memberSnpMap.insert(pair<int, SnpInfo*>(snp->index, snp)).second) {
                        anno->memberSnpVec.push_back(snp);
                        anno->size++;
                    }
                }
            }
        }
        in.close();
        in.clear();
    }

    unsigned numUnannoSnp = 0;
    unsigned numAnnoSnp = 0;
    unsigned numMultiAnnoSnp = 0;
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        if (!snp->included) continue;
        if (snp->numAnnos == 0) ++numUnannoSnp;
        else {
            ++numAnnoSnp;
            if (snp->numAnnos > 1) ++numMultiAnnoSnp;
        }
    }
    
    
    annoInfoVec.clear();
    annoNames.clear();
    numAnnos = 0;
    unsigned meanAnnoSize = 0;
    for (unsigned i=0; i<line; ++i) {
        anno = annoVec[i];
        if (anno->size) {
            annoInfoVec.push_back(anno);
            annoNames.push_back(anno->label);
            ++numAnnos;
            meanAnnoSize += anno->size;
       }
    }
    meanAnnoSize /= numAnnos;
    
    if (flank) cout << "add " << flank << " kb to both sides of each annotation." << endl;
    cout << numAnnos << " nonempty annotations (" << line << " annotations in total)." << endl;
    cout << numAnnoSnp << " annotated SNPs (" << numMultiAnnoSnp << " SNPs with more than one annotation) and " << numUnannoSnp << " unannotated SNPs." << endl;
    cout << meanAnnoSize << " SNPs per annotation on average." << endl;
}

void Data::readLDscoreFile(const string &ldscoreFile) {
    ifstream in(ldscoreFile.c_str());
    if (!in) throw ("Error: can not open the LD score file [" + ldscoreFile + "] to read.");
    cout << "Reading SNP LD scores from [" + ldscoreFile + "]." << endl;

    readLDscore = true;
    
    SnpInfo *snp;
    map<string, SnpInfo*>::iterator it;
    string id;
    float ldsc;
    unsigned line = 0, match = 0;
    while(in >> id >> ldsc) {
        ++line;
        it = snpInfoMap.find(id);
        if (it == snpInfoMap.end()) continue;
        snp = it->second;
        if (snp->included) {
            snp->ldsc = ldsc;
            ++match;
        }
    }
    in.close();
    
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        if (!snp->included) continue;
        if (!snp->ldsc) {
            snp->included = false;
        }
    }
    
    cout << "Read LD scores for " << match << " matched SNPs (in total " << line << " SNPs)." << endl;
}

void Data::makeAnnowiseSparseLDM(const vector<SparseVector<float> > &ZPZsp, const vector<AnnoInfo *> &annoInfoVec, const vector<SnpInfo*> &snpInfoVec) {
    // for memory efficiency, only store the upper triangular off-diagonals
    cout << "\nMaking annotation-wise sparse LD matrix ..." << endl;
    annowiseZPZsp.resize(numAnnos);
    annowiseZPZdiag.resize(numAnnos);
    
    long chunkSize = numAnnos/omp_get_max_threads();
#pragma omp parallel for schedule(dynamic, chunkSize)
    for (unsigned i=0; i<numAnnos; ++i) {
        AnnoInfo *anno = annoInfoVec[i];
        annowiseZPZsp[i].resize(anno->size, anno->size);
        annowiseZPZdiag[i].resize(anno->size);
        SnpInfo *snpj, *snpk;
        unsigned r, c;
        float v;
        for (unsigned j=0; j<anno->size; ++j) {
//            if (!(j%10000))
//                cout << "  annotaion " << std::setw(2) << std::left << i+1 << " snp " << std::setw(6) << std::left << j << "\r";
            snpj = anno->memberSnpVec[j];
            r = j;
            c = 0;
            VectorXf dense;
            dense.setZero(numIncdSnps);
            for (SparseVector<float>::InnerIterator it(ZPZsp[snpj->index]); it; ++it) {
                dense[it.index()] = it.value();
            }
            for (unsigned k=j+1; k<anno->size; ++k) {
                snpk = anno->memberSnpVec[k];
                if (snpj->chrom != snpk->chrom) continue;
                //                v = ZPZsp[snpj->index].coeff(snpk->index);
                v = dense[snpk->index];
                if (v) {
                    annowiseZPZsp[i].insert(c++, r) = v;
                }
            }
            annowiseZPZdiag[i][j] = dense[snpj->index];
//            if (j==anno->size-1)
//                cout << "  annotaion " << std::setw(2) << std::left << i+1 << " snp " << std::setw(6) << std::left << j+1 << " nonzeros " << annowiseZPZsp[i].nonZeros() << "\r";
        }
//        if (myMPI::rank==0) cout << endl;
        if (!(i%100))
            cout << "  annotaion " << std::setw(2) << std::left << i+1 << " size " << std::setw(6) << std::left << anno->size << " nonzeros " << annowiseZPZsp[i].nonZeros() << endl;
        annowiseZPZsp[i].makeCompressed();
    }
}

void Data::getZPZspmat(){ // get ZPZspmat from ZPZsp which is a vector of sparse vectors
    ZPZspmat.resize(numIncdSnps, numIncdSnps);
    vector<Triplet<float> > tripletList;
    tripletList.reserve(windSize.cast<double>().sum());
    float val = 0.0;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snpi = incdSnpInfoVec[i];
        if (!(i % 100000)) cout << "  making sparse LD matrix for SNP " << i << " " << windSize[i] << " " << snpi->windSize << " " << ZPZsp[i].size() << " " << ZPZsp[i].nonZeros() << endl;
        for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
            tripletList.push_back(Triplet<float>(i, it.index(), it.value()));
            //if (it.index() > numIncdSnps - 1) cout << i << " " << it.index() << " " << it.value() << endl;
        }
    }
    ZPZspmat.setFromTriplets(tripletList.begin(), tripletList.end());
    ZPZspmat.makeCompressed();
    tripletList.clear();
}

void Data::getZPZmat(){ // get ZPZmat from ZPZ which is a vector of vectors
    ZPZmat.setZero(numIncdSnps, numIncdSnps);
    for (unsigned i=0; i<numIncdSnps; ++i) {
        ZPZmat.col(i) = ZPZ[i];
    }
}

void Data::binSnpByLDrsq(const float rsqThreshold, const string &title){
// only work for full LD matrix at the moment
    
    // initialise a unique window ID for each SNP
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        snp->window = i;
    }
    
    // merge SNP window ID if two SNPs are in high LD
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        VectorXf rsq = ZPZ[i].cwiseProduct(ZPZ[i]);
        for (unsigned j=0; j<numIncdSnps; ++j) {
            if (rsq[j] > rsqThreshold) incdSnpInfoVec[j]->window = snp->window;
        }
    }

    // reorder the window ID
    map<int, vector<SnpInfo*> > uniqueWinIDmap;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        uniqueWinIDmap[snp->window].push_back(snp);
    }
    map<int, vector<SnpInfo*> >::iterator it, end = uniqueWinIDmap.end();
    int newWinID = 1;
    for (it = uniqueWinIDmap.begin(); it!=end; ++it) {
        for (unsigned i=0; i<it->second.size(); ++i) {
            it->second[i]->window = newWinID;
        }
        ++newWinID;
    }
    
    // output SNP window ID info
    string outfile = title + ".window";
    ofstream out(outfile.c_str());
    out << boost::format("%18s %8s\n") % "SNP" % "Window";
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        out << boost::format("%18s %8s\n") % snp->ID % snp->window;
    }
    out.close();
    
    cout << "Based on LD Rsq threshold of " << rsqThreshold << ", " << numIncdSnps << " SNPs are grouped into " << newWinID << " bins." << endl;
}

void Data::readWindowFile(const string &windowFile){
    ifstream in(windowFile.c_str());
    if (!in) throw ("Error: can not open the SNP window file [" + windowFile + "] to read.");
    cout << "Reading SNP window info from [" + windowFile + "]." << endl;
    
    SnpInfo *snp;
    map<string, SnpInfo*>::iterator it;
    string snpID;
    int windowID;
    set<int> uniqueWinID;
    unsigned line=0, match=0;
    string header;
    getline(in, header);
    while (in >> snpID >> windowID) {
        ++line;
        it = snpInfoMap.find(snpID);
        if (it == snpInfoMap.end()) continue;
        snp = it->second;
        snp->window = windowID;
        ++match;
        uniqueWinID.insert(windowID);
    }
    
    numWindows = uniqueWinID.size();

    cout << match << " matched SNPs in the SNP window file (in total " << line << " SNPs) are binned into " << numWindows << " windows." << endl;
}

void Data::binSnpByWindowID(){
    // reorder the window ID in case the original window ID are not continous
    map<int, int> orderedWinIDmap;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        orderedWinIDmap[snp->window] = 0;
    }
    map<int, int>::iterator it, end = orderedWinIDmap.end();
    int newWinID = 0;
    for (it = orderedWinIDmap.begin(); it!=end; ++it) {
        it->second = newWinID++;
    }

    // set up window-snp map
    windowSnpIdxVec.resize(numWindows);
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        newWinID = orderedWinIDmap[snp->window];
//        cout << i << " " << snp->window << " " << newWinID << " " << windowSnpIdxVec.size() << endl;
        windowSnpIdxVec[newWinID].push_back(snp->index);
    }
    
    windSize.resize(numWindows);
    for (unsigned i=0; i<numWindows; ++i) {
        windSize[i] = windowSnpIdxVec[i].size();
    }
    
    cout << "Average window size for each SNP is " << windSize.mean() << endl;
}

void Data::filterSnpByLDrsq(const float rsqThreshold){
    cout << "Filtering SNPs by rsq ... " << rsqThreshold << endl;
    unsigned numExcdSnpsBeforeFilter = 0;
    unsigned numExcdSnpsAfterFilter = 0;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snpi = incdSnpInfoVec[i];
        if (!snpi->included) ++numExcdSnpsBeforeFilter;
    }
    if (sparseLDM) {
        for (unsigned i=0; i<numIncdSnps; ++i) {
            SnpInfo *snpi = incdSnpInfoVec[i];
            if (!snpi->included) continue;
            float rsq = 0.0;
            for (SparseVector<float>::InnerIterator it(ZPZsp[i]); it; ++it) {
                SnpInfo *snpj = incdSnpInfoVec[it.index()];
                rsq = it.value()*it.value();
                if (snpi!=snpj && rsq > rsqThreshold) {
                    if (snpi->gwas_pvalue < snpj->gwas_pvalue) snpj->included = false;
                    else snpi->included = false;
                }
            }
        }
    }
    else {
        for (unsigned i=0; i<numIncdSnps; ++i) {
            SnpInfo *snpi = incdSnpInfoVec[i];
            if (!snpi->included) continue;
            VectorXf rsq = ZPZ[i].cwiseProduct(ZPZ[i]);
            for (unsigned j=0; j<numIncdSnps; ++j) {
                SnpInfo *snpj = incdSnpInfoVec[j];
                if (!snpj->included) continue;
                if (snpi!=snpj && rsq[j] > rsqThreshold) {
                    if (snpi->gwas_pvalue < snpj->gwas_pvalue) snpj->included = false;
                    else snpi->included = false;
                }
            }
        }
    }
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snpi = incdSnpInfoVec[i];
        if (!snpi->included) ++numExcdSnpsAfterFilter;
    }
    
    // restore the LD window information for preparing reading LD matrix again
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        snp->windStart = snp->windStartOri;
        snp->windSize = snp->windSizeOri;
        snp->windEnd = snp->windEndOri;
    }

    cout << "Removed " << numExcdSnpsAfterFilter - numExcdSnpsBeforeFilter << " SNPs with LD R-square above " << rsqThreshold << "." << endl;
}

void Data::readLDmatrixTxtFile(const string &ldmatrixFile) {
    Gadget::Timer timer;
    timer.setTime();
    
    Gadget::Tokenizer token;
    token.getTokens(ldmatrixFile, ".");
    string ldmType = token[token.size()-2];
    if (ldmType != "full") throw("Error: can only read full matrix from text file.");
    
    VectorXi windStartLDM(numSnps);
    VectorXi windSizeLDM(numSnps);
    
    windStart.resize(numIncdSnps);
    windSize.resize(numIncdSnps);
    
    SnpInfo *snpi, *snpj;
    
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snpi = snpInfoVec[i];
        windStartLDM[i] = snpi->windStart;
        windSizeLDM[i]  = snpi->windSize;
    }

    ifstream in(ldmatrixFile.c_str());
    if (!in) {
        throw("Error: cannot open LD matrix file " + ldmatrixFile);
    }

    resizeWindow(incdSnpInfoVec, windStartLDM, windSizeLDM, windStart, windSize);
    
    if (numIncdSnps == 0) throw ("Error: No SNP is retained for analysis.");

    cout << "Reading " + ldmType + " LD matrix from [" + ldmatrixFile + "]..." << endl;

    float rsq = 0.0;

    ZPZ.resize(numIncdSnps);
    ZPZdiag.resize(numIncdSnps);
    Gadget::Tokenizer colData;
    string inputStr;
    string sep(" \t,");
    
    for (unsigned i=0, inci=0; i<numSnps; i++) {
        snpi = snpInfoVec[i];
        getline(in,inputStr);
        
        if (!snpi->included) continue;
        
        colData.getTokens(inputStr, sep);
        
        ZPZ[inci].resize(windSize[inci]);
        snpi->ldSamplVar = 0.0;
        snpi->ldSum = 0.0;
        if (!readLDscore) snpi->ldsc = 0.0;
        snpi->numNonZeroLD = snpi->windSize;
        
        for (unsigned j=0, incj=0; j<windSizeLDM[i]; ++j) {
            unsigned idx = windStartLDM[i] + j;
            snpj = snpInfoVec[idx];
            if (snpj->included) {
                float ld = atof(colData[idx].c_str());
                ZPZ[inci][incj++] = ld;
                rsq = ld*ld;
                snpi->ldSamplVar += (1.0f-rsq)*(1.0f-rsq)/snpi->sampleSize;
                snpi->ldSum += ld;
                if (!readLDscore) snpi->ldsc += rsq;
                if (snpj == snpi)
                    ZPZdiag[inci] = ld;
            }
        }
        
        if (++inci == numIncdSnps) break;
    }
    
    in.close();
    
    timer.getTime();
    
    //    cout << "Window width " << windowWidth << " Mb." << endl;
    displayAverageWindowSize(windSize);
    cout << "LD matrix diagonal mean " << ZPZdiag.mean() << " sd " << sqrt(Gadget::calcVariance(ZPZdiag)) << "." << endl;
    if (ZPZdiag.mean() < 0.8 || ZPZdiag.mean() > 1.2) throw("ERROR: The mean of LD matrix diagonal values is expected to be close to one. Something is wrong with the LD matrix!");
    cout << "Read LD matrix for " << numIncdSnps << " SNPs (time used: " << timer.format(timer.getElapse()) << ")." << endl;
}

