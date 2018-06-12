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

void Data::readFamFile(const string &famFile){
    // ignore phenotype column
    ifstream in(famFile.c_str());
    if (!in) throw ("Error: can not open the file [" + famFile + "] to read.");
    if (myMPI::rank==0)
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
    if (myMPI::rank==0)
        cout << numInds << " individuals to be included from [" + famFile + "]." << endl;
}

void Data::readBimFile(const string &bimFile) {
    // Read bim file: recombination rate is defined between SNP i and SNP i-1
    ifstream in(bimFile.c_str());
    if (!in) throw ("Error: can not open the file [" + bimFile + "] to read.");
    if (myMPI::rank==0)
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
    if (myMPI::rank==0)
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

void Data::readBedFile(const string &bedFile){
    unsigned i = 0, j = 0;
    
    if (numIncdSnps == 0) throw ("Error: No SNP is retained for analysis.");
    if (numKeptInds == 0) throw ("Error: No individual is retained for analysis.");
    
    Z.resize(numKeptInds, numIncdSnps);
    ZPZdiag.resize(numIncdSnps);
    snp2pq.resize(numIncdSnps);
    
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

    unsigned numKeptInds_all;
    MPI_Allreduce(&numKeptInds, &numKeptInds_all, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    
    // Read genotypes
    SnpInfo *snpInfo = NULL;
    IndInfo *indInfo = NULL;
    unsigned snp = 0;
    unsigned nmiss=0, nmiss_all;
    float sum=0.0, sum_all=0.0, mean_all;
    
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
        
        MPI_Allreduce(&sum, &sum_all, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&nmiss, &nmiss_all, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
        
        // fill missing values with the mean
        mean_all = sum_all/float(numKeptInds_all - nmiss_all);
        if (nmiss) {
            for (i=0; i<numKeptInds; ++i) {
                if (Z(i,snp) == -9) Z(i,snp) = mean_all;
            }
        }
        
        // compute allele frequency
        snpInfo->af = 0.5f*mean_all;
        snp2pq[snp] = snpInfo->twopq = 2.0f*snpInfo->af*(1.0f-snpInfo->af);
        
        //cout << "snp " << snp << "     " << Z.col(snp).sum() << endl;
        
        Z.col(snp).array() -= mean_all; // center column by 2p rather than the real mean

        if (++snp == numIncdSnps) break;
    }
    fclose(in);
    
    
    // standardize genotypes
    VectorXf colsums = Z.colwise().sum();
    VectorXf colsums_all(numIncdSnps);
    
    MPI_Allreduce(&colsums[0], &colsums_all[0], numIncdSnps, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    
    //Z.rowwise() -= colsums_all.transpose()/numKeptInds_all;  // center
    VectorXf my_ZPZdiag = Z.colwise().squaredNorm();
    
    MPI_Allreduce(&my_ZPZdiag[0], &ZPZdiag[0], numIncdSnps, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    
    if (myMPI::rank==0)
        cout << "Genotype data for " << numKeptInds_all << " individuals and " << numIncdSnps << " SNPs are included from [" + bedFile + "]." << endl;
}

void Data::readPhenotypeFile(const string &phenFile, const unsigned mphen) {
    // NA: missing phenotype
    ifstream in(phenFile.c_str());
    if (!in) throw ("Error: can not open the phenotype file [" + phenFile + "] to read.");
    if (myMPI::rank==0)
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
    if (myMPI::rank==0)
        cout << "Non-missing phenotypes of trait " << mphen << " of " << line << " individuals are included from [" + phenFile + "]." << endl;
}

void Data::keepMatchedInd(const string &keepIndFile, const unsigned keepIndMax){  // keepIndFile is optional
    map<string, IndInfo*>::iterator it, end=indInfoMap.end();
    IndInfo *ind = NULL;
    vector<string> keep;
    keep.reserve(numInds);
    unsigned cnt=0;
    for (unsigned i=0; i<numInds; ++i) {
        ind = indInfoVec[i];
        ind->kept = false;
        if (ind->phenotype!=-9) {
            if (keepIndMax > cnt++)
                keep.push_back(ind->catID);
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
    
    unsigned numKeptInds_all = 0;

    if (myMPI::partition == "byrow") {
        unsigned total_size = (unsigned) keep.size();
        unsigned batch_size = total_size/myMPI::clusterSize;
        unsigned my_start = myMPI::rank*batch_size;
        unsigned my_end = (myMPI::rank+1)==myMPI::clusterSize ? total_size : my_start + batch_size;
        unsigned my_size = my_end - my_start;
                
        myMPI::iStart = my_start;
        myMPI::iSize  = my_size;
        
        vector<string>::const_iterator first = keep.begin() + my_start;
        vector<string>::const_iterator last  = keep.begin() + my_end;
        vector<string> my_keep(first, last);
        
        for (unsigned i=0; i<my_size; ++i) {
            it = indInfoMap.find(my_keep[i]);
            if (it == end) {
                Gadget::Tokenizer token;
                token.getTokens(my_keep[i], ":");
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
        
        y.setZero(numKeptInds);
        for (unsigned i=0; i<numKeptInds; ++i) {
            y[i] = keptIndInfoVec[i]->phenotype;
        }
        float my_ypy = (y.array()-y.mean()).square().sum();
        
        MPI_Allreduce(&my_ypy, &ypy, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        
        MPI_Allreduce(&numKeptInds, &numKeptInds_all, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
        
        if (myMPI::rank==0) {
            cout << numKeptInds_all << " matched individuals are kept." << endl;
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        printf("%d individuals assigned to processor %s at rank %d.\n", numKeptInds, myMPI::processorName, myMPI::rank);
        //cout << numKeptInds << " individuals assigned to processor " << myMPI::processorName << " at rank " << myMPI::rank << "." << endl;
    }
    else {
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
        numKeptInds_all = numKeptInds;
        
        y.setZero(numKeptInds);
        for (unsigned i=0; i<numKeptInds; ++i) {
            y[i] = keptIndInfoVec[i]->phenotype;
        }
        ypy = (y.array()-y.mean()).square().sum();
        
        if (myMPI::rank==0) {
            cout << numKeptInds << " matched individuals are kept." << endl;
        }
    }
}

void Data::initVariances(const float heritability){
    float varPhenotypic;
    if (myMPI::partition == "byrow") {
        unsigned numKeptInds_all = 0;
        MPI_Allreduce(&numKeptInds, &numKeptInds_all, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
        varPhenotypic = ypy/numKeptInds_all;
    } else {
        varPhenotypic = ypy/numKeptInds;
    }
    varGenotypic = varPhenotypic * heritability;
    varResidual  = varPhenotypic - varGenotypic;
    //cout <<varPhenotypic<<" " <<varGenotypic << " " <<varResidual << endl;
}

void Data::readCovariateFile(const string &covarFile){
    if (!covarFile.empty()) {
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
                for (unsigned i=0; i<numCovariates; ++i)
                    fixedEffectNames[i+1] = colData[i+2];
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
        
        if (myMPI::rank==0)
            cout << "Read " << numCovariates << " covariates from [" + covarFile + "]." << endl;
        
        X.resize(numKeptInds, numFixedEffects);
        for (unsigned i=0; i<numKeptInds; ++i) {
            ind = keptIndInfoVec[i];
            if (ind->covariates.size() < numFixedEffects) {
                cout << "Error: Individual " + ind->famID + " " + ind->indID + " has missing covariate(s)!" << endl;
            }
        }
        for (unsigned i=0; i<numKeptInds; ++i) {
            ind = keptIndInfoVec[i];
            X.row(i) = ind->covariates;
        }
        VectorXf my_XPXdiag = X.colwise().squaredNorm();
        XPXdiag.setZero(numFixedEffects);
        MPI_Allreduce(&my_XPXdiag[0], &XPXdiag[0], numFixedEffects, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    }
    else {
        // only intercept for now
        numFixedEffects = 1;
        fixedEffectNames = {"Intercept"};
        X.setOnes(numKeptInds,1);
        XPX.resize(1,1);
        XPXdiag.resize(1);
        XPy.resize(1);
        
        if (myMPI::partition == "byrow") {
            unsigned numKeptInds_all;
            MPI_Allreduce(&numKeptInds, &numKeptInds_all, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
            
            XPX << numKeptInds_all;
            XPXdiag << numKeptInds_all;
            
            float sum = y.sum();
            unsigned sum_all;
            MPI_Allreduce(&sum, &sum_all, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            XPy << sum_all;
        }
        else {
            XPX << numKeptInds;
            XPXdiag << numKeptInds;
            XPy << y.sum();
        }
    }
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
    if (myMPI::rank==0) {
        cout << numSkeletonSnps << " skeleton SNPs are included." << endl;
    }
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
    if (myMPI::rank==0) {
        cout << cnt << " SNPs in the MHC region (Chr6:28-34Mb) are excluded." << endl;
    }
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
    
    
    if (myMPI::partition == "bycol") {  // MPI by chromosome
        vector<int> chromVec;
        for (set<int>::iterator it=chromosomes.begin(); it!=chromosomes.end(); ++it) {
            chromVec.push_back(*it);
        }
        unsigned my_chrom = (myMPI::rank+1)==myMPI::clusterSize ? 99 : chromVec[myMPI::rank];

        SnpInfo *snp;
        for (unsigned i=0; i<numSnps; ++i) {  // each core takes one chromosome, if #chrom > #core, the last core takes all rest chroms
            snp = snpInfoVec[i];
            if (my_chrom != 99) {
                if (snp->chrom != my_chrom) snp->included = false;
            } else {
                if (snp->chrom < chromVec[myMPI::rank]) snp->included = false;
            }
        }
        
        incdSnpInfoVec = makeIncdSnpInfoVec(snpInfoVec);
        numIncdSnps = (unsigned) incdSnpInfoVec.size();
        snp2pq.resize(numIncdSnps);
        
        // setup MPI for collecting SNP info
        myMPI::iSize = numIncdSnps;
        myMPI::iStart = incdSnpInfoVec[0]->index;
        myMPI::srcounts.resize(myMPI::clusterSize);
        myMPI::displs.resize(myMPI::clusterSize);
        
        MPI_Allgather(&myMPI::iSize, 1, MPI_INT, &myMPI::srcounts[0], 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(&myMPI::iStart, 1, MPI_INT, &myMPI::displs[0], 1, MPI_INT, MPI_COMM_WORLD);

        
        reindexSnp(incdSnpInfoVec);  // reindex based on the snplist for each core

        //cout << "rank " << myMPI::rank << " " << myMPI::iStart << " " << myMPI::iSize << endl;
        
        unsigned numIncdSnps_all;
        MPI_Reduce(&numIncdSnps, &numIncdSnps_all, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
        if (myMPI::rank==0) cout << numIncdSnps_all << " SNPs are included." << endl;
        MPI_Barrier(MPI_COMM_WORLD);
        cout << numIncdSnps << " SNPs assigned to processor " << myMPI::processorName << " at rank " << myMPI::rank << "." << endl;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    else {
        incdSnpInfoVec = makeIncdSnpInfoVec(snpInfoVec);
        numIncdSnps = (unsigned) incdSnpInfoVec.size();
        reindexSnp(incdSnpInfoVec);
        snp2pq.resize(numIncdSnps);
    }
    
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

    if (myMPI::rank==0) cout << numIncdSnps << " SNPs on " << numChroms << " chromosomes are included." << endl;
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
    if (myMPI::rank==0)
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
    if (myMPI::rank==0)
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
        if(!(i%10000) && myMPI::rank==0)
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
    
    long numberWindows = windStartVec.size();

    windStart = VectorXi::Map(&windStartVec[0], numberWindows);
    windSize.setZero(numberWindows);

    for (unsigned i=0; i<numberWindows; ++i) {
        if (i != numberWindows-1)
            windSize[i] = windStart[i+1] - windStart[i];
        else
            windSize[i] = numIncdSnps - windStart[i];
    }
    
    if (myMPI::rank==0)
        cout << "Created " << numberWindows << " non-overlapping " << windowWidth/1e3 << "kb windows with average size of " << windSize.sum()/float(numberWindows) << " SNPs." << endl;
}


void Data::buildSparseMME(const string &bedFile, const unsigned windowWidth){
    if (myMPI::rank==0)
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

            if(!(inc%1000) && myMPI::rank==0) {
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

    if (myMPI::rank==0) {
        cout << "Average window size " << windSize.sum()/numIncdSnps << endl;
        cout << "Genotype data for " << numKeptInds << " individuals and " << numIncdSnps << " SNPs are included from [" + bedFile + "]." << endl;
        cout << "Construction of sparse MME completed (time used: " << timer.format(timer.getElapse()) << ")" << endl;
    }
    
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


void Data::outputSnpResults(const VectorXf &posteriorMean, const VectorXf &posteriorSqrMean, const VectorXf &pip, const string &filename) const {
    if (myMPI::rank) return;
    ofstream out(filename.c_str());
    out << boost::format("%6s %20s %6s %12s %6s %8s %12s %12s %8s %8s\n")
    % "ID"
    % "Name"
    % "Chrom"
    % "Position"
    % "Allele"
    % "Freq"
    % "Effect"
    % "SE"
    % "PIP"
    % "Window";
    for (unsigned i=0, idx=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        if(!fullSnpFlag[i]) continue;
        if(snp->isQTL) continue;
        out << boost::format("%6s %20s %6s %12s %6s %8.6f %12.6f %12.6f %8.3f %8s\n")
        % (idx+1)
        % snp->ID
        % snp->chrom
        % snp->physPos
        % snp->a2
        % snp->af
        % posteriorMean[idx]
        % sqrt(posteriorSqrMean[idx]-posteriorMean[idx]*posteriorMean[idx])
        % pip[idx]
        % snp->window;
        ++idx;
    }
    out.close();
}

void Data::inputSnpResults(const string &snpResFile){
    ifstream in(snpResFile.c_str());
    if (!in) throw ("Error: can not open the SNP result file [" + snpResFile + "] to read.");
    if (myMPI::rank==0)
        cout << "Reading SNP results from [" + snpResFile + "]." << endl;
    
    SnpInfo *snp;
    map<string, SnpInfo*>::iterator it;
    string name, a2;
    int id, chrom, pos, window;
    float freq, effect, se, pip;
    unsigned line=0, match=0;
    string header;
    getline(in, header);
    while (in >> id >> name >> chrom >> pos >> a2 >> freq >> effect >> se >> pip >> window) {
        ++line;
        it = snpInfoMap.find(name);
        if (it == snpInfoMap.end()) continue;
        snp = it->second;
        if (snp->included) {
            if (snp->a2 == a2)
                snp->effect = effect;
            else
                snp->effect = -effect;
            ++match;
        }
    }
    in.close();
    
    if (myMPI::rank==0) {
        cout << match << " matched SNPs in the SNP result file (in total " << line << " SNPs)." << endl;
    }
}


void Data::summarizeSnpResults(const SparseMatrix<float> &snpEffects, const string &filename) const {
    if (myMPI::rank==0) {
        cout << "SNP results to be summarized in " << filename << endl;
    }
    unsigned nrow = snpEffects.rows();
    VectorXf effectSum(numIncdSnps), effectMean(numIncdSnps);
    VectorXf pipSum(numIncdSnps), pip(numIncdSnps);  // posterior inclusion probability
    for (unsigned i=0; i<numIncdSnps; ++i) {
        effectSum[i] = snpEffects.col(i).sum();
        pipSum[i] = (VectorXf(snpEffects.col(i)).array()!=0).count();
    }
    MPI_Allreduce(MPI_IN_PLACE, &nrow, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&effectSum[0], &effectMean[0], numIncdSnps, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&pipSum[0], &pip[0], numIncdSnps, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    effectMean /= (float)nrow;
    pip /= (float)nrow;
    
    if (myMPI::rank) return;
    
    ofstream out(filename.c_str());
    out << boost::format("%6s %20s %6s %12s %6s %8s %12s %8s %8s\n")
    % "Id"
    % "Name"
    % "Chrom"
    % "Position"
    % "Allele"
    % "Freq"
    % "Effect"
    % "PIP"
    % "Window";
    for (unsigned i=0, idx=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        if(!fullSnpFlag[i]) continue;
        out << boost::format("%6s %20s %6s %12s %6s %8.3f %12.6f %8.3f %8s\n")
        % (idx+1)
        % snp->ID
        % snp->chrom
        % snp->physPos
        % snp->a2
        % snp->af
        % effectMean[idx]
        % pip[idx]
        % snp->window;
        ++idx;
    }
    out.close();
}

void Data::outputFixedEffects(const MatrixXf &fixedEffects, const string &filename) const {
    if (myMPI::rank) return;
    ofstream out(filename.c_str());
    long nrow = fixedEffects.rows();
    VectorXf mean = fixedEffects.colwise().mean();
    VectorXf sd = (fixedEffects.rowwise() - mean.transpose()).colwise().squaredNorm().cwiseSqrt()/sqrt(nrow);
    for (unsigned i=0; i<numFixedEffects; ++i) {
        out << boost::format("%20s %12.6f %12.6f\n") % fixedEffectNames[i] %mean[i] %sd[i];
    }
    out.close();
}

void Data::outputWindowResults(const VectorXf &posteriorMean, const string &filename) const {
    if (myMPI::rank) return;
    ofstream out(filename.c_str());
    out << boost::format("%6s %8s\n") %"Id" %"PIP";
    for (unsigned i=0; i<posteriorMean.size(); ++i) {
        out << boost::format("%6s %8.3f\n")
        % (i+1)
        % posteriorMean[i];
    }
    out.close();
}

void Data::readGwasSummaryFile(const string &gwasFile){
    ifstream in(gwasFile.c_str());
    if (!in) throw ("Error: can not open the GWAS summary data file [" + gwasFile + "] to read.");
    if (myMPI::rank==0)
        cout << "Reading GWAS summary data from [" + gwasFile + "]." << endl;
    
    SnpInfo *snp;
    map<string, SnpInfo*>::iterator it;
    string id, allele1, allele2, freq, b, se, pval, n;
    unsigned line=0, match=0;
    unsigned incon=0;
    while (in >> id >> allele1 >> allele2 >> freq >> b >> se >> pval >> n) {
        ++line;
        it = snpInfoMap.find(id);
        if (it == snpInfoMap.end()) continue;
        snp = it->second;
        if (!snp->included) continue;
        if (allele1 == snp->a1 && allele2 == snp->a2) {
            snp->gwas_b  = atof(b.c_str());
            snp->gwas_af = atof(freq.c_str());
            snp->gwas_se = atof(se.c_str());
            snp->gwas_n  = atof(n.c_str());
            ++match;
        } else if (allele1 == snp->a2 && allele2 == snp->a1) {
            snp->gwas_b  = -atof(b.c_str());
            snp->gwas_af = 1.0-atof(freq.c_str());
            snp->gwas_se = atof(se.c_str());
            snp->gwas_n  = atof(n.c_str());
            ++match;
        } else {
            cout << "WARNING: SNP " + id + " has inconsistent allele coding in between the reference and GWAS samples." << endl;
            snp->included = false;
            ++incon;
        }
        if (snp->gwas_af==0 || snp->gwas_af==1) throw ("Error: SNP " + id + " is a fixed SNP!");
    }
    in.close();
    
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        if (!snp->included) continue;
        if (snp->gwas_b == -999) {
            snp->included = false;
        }
    }

    if (myMPI::rank==0) {
        if (incon) cout << "removed " << incon << " SNPs with inconsistent allele coding in between the reference and GWAS samples." << endl;
        cout << match << " matched SNPs in the GWAS summary data (in total " << line << " SNPs)." << endl;
    }

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
        
        if (snp2pq[incj]==0) throw ("Error: " + snpj->ID + " is a fixed SNP!");
        
        // standardize genotypes
        D[incj] = snp2pq[incj]*snpj->sampleSize;
        
        ZP.row(incj) = (ZP.row(incj).array() - mean)/sqrt(D[incj]);
        
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

            if(!(inck%1000) && myMPI::rank==0) cout << " read snp " << inck << "\r" << flush;

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
            
            if (snp2pq[inck]==0) throw ("Error: " + snpk->ID + " is a fixed SNP!");
            
            // standardize genotypes
            D[inck] = snp2pq[inck]*snpk->sampleSize;
            
            Zk = (Zk.array() - mean)/sqrt(D[inck]);
            
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
            
            if (snp2pq[inck]==0) throw ("Error: " + snpk->ID + " is a fixed SNP!");
            
            // standardize genotypes
            D[inck] = snp2pq[inck]*snpk->sampleSize;
            
            Zk = (Zk.array() - mean)/sqrt(D[inck]);
            
            denseZPZ.col(inck) = ZP * Zk;
            
            if(!(inck%1000) && myMPI::rank==0) cout << " read snp " << inck << "\r" << flush;
            
            ++inck;
        }
    }

    fclose(in2);
    
//    cout << denseZPZ << endl;
    
    
    VectorXf snp2pqTmp(numIncdSnps);
    for (unsigned i=0, j=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        if (snp->included) {
            snp2pqTmp[j++] = snp2pq[i];
        }
    }
    snp2pq = snp2pqTmp;

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
            snp->ldSamplVar = (1.0 - denseZPZ.row(i).array().square()).square().sum()/snp->sampleSize;
            snp->ldSum = denseZPZ.row(i).sum();
       }
    }
    else if (LDmatType == "band") {
        ZPZ.resize(numSnpInRange);
        if (windowWidth) {  // based on the given window width
            for (unsigned i=0; i<numSnpInRange; ++i) {
                SnpInfo *snp = incdSnpInfoVec[start+i];
                ZPZ[i] = denseZPZ.row(i).segment(snp->windStart, snp->windSize);
                snp->ldSamplVar = (1.0 - ZPZ[i].array().square()).square().sum()/snp->sampleSize;
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
                snp->ldSamplVar = (1.0 - ZPZ[i].array().square()).square().sum()/snp->sampleSize;
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
                            snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snpi->sampleSize;
                            snpi->ldSum += denseZPZ(i,j);
                        }
                        else {
                            if (abs(denseZPZ(i,j)) < LDthreshold) denseZPZ(i,j) = 0;
                            else {
                                rsq = denseZPZ(i,j)*denseZPZ(i,j);
                                snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snpi->sampleSize;
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
                            snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snpi->sampleSize;
                            snpi->ldSum += denseZPZ(i,j);
//                            cout << j << " " << denseZPZ(i,j) << endl;
//                            ++cnt;
                        }
                        else {
                            if (denseZPZ(i,j)*denseZPZ(i,j)*snpi->sampleSize < chisqThreshold) denseZPZ(i,j) = 0;
                            else {
                                rsq = denseZPZ(i,j)*denseZPZ(i,j);
                                snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snpi->sampleSize;
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
                            snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snpi->sampleSize;
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
                        if (denseZPZ(i,j)*denseZPZ(i,j)*snpi->sampleSize < chisqThreshold) denseZPZ(i,j) = 0;
                        else {
                            rsq = denseZPZ(i,j)*denseZPZ(i,j);
                            snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snpi->sampleSize;
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
    % "A2Freq"
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
        % snp->sampleSize
        % snp->ldSamplVar
        % snp->ldSum;
        if (LDmatType == "sparse") {
            fwrite(ZPZsp[i].innerIndexPtr(), sizeof(unsigned), ZPZsp[i].nonZeros(), out2);
            fwrite(ZPZsp[i].valuePtr(), sizeof(float), ZPZsp[i].nonZeros(), out2);
            if (writeLdmTxt) out3 << ZPZsp[i].transpose() << endl;
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
    float windSizeMean = 0.0;
    float windSizeSqMean = 0.0;
    long size = windSize.size();
    for (unsigned i=0; i<size; ++i) {
        windSizeMean += (windSize[i] - windSizeMean)/(i+1);
        windSizeSqMean += (windSize[i]*windSize[i] - windSizeSqMean)/(i+1);
    }
    cout << "Per-SNP window size mean " << windSizeMean << " sd " << windSizeSqMean-windSizeMean*windSizeMean << "." << endl;
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
    unsigned idx, windStart, windEnd, windSize, windWidth;
    long sampleSize;
    bool skeleton;
    getline(in, header);
    Gadget::Tokenizer token;
    token.getTokens(header, " ");

    if (token.back() == "Skeleton") {
        while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2 >> af >> idx >> windStart >> windEnd >> windSize >> windWidth >> sampleSize >> ldSamplVar >> ldSum >> skeleton) {
            SnpInfo *snp = new SnpInfo(idx, id, allele1, allele2, chr, genPos, physPos);
            snp->af = af;
            snp->twopq = 2.0*af*(1.0-af);
            snp->windStart = windStart;
            snp->windEnd = windEnd;
            snp->windSize = windSize;
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
            snp->windStart = windStart;
            snp->windEnd = windEnd;
            snp->windSize = windSize;
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
            snp->windStart = windStart;
            snp->windEnd = windEnd;
            snp->windSize = windSize;
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
            snpi->ldsc = 0.0;
            
            for (unsigned j=0; j<windSizeLDM[i]; ++j) {
                snpj = snpInfoVec[d[j]];
                if (snpj->included) {
                    ZPZsp[inci].insertBack(snpj->index) = v[j];
                    rsq = v[j]*v[j];
                    snpi->ldSamplVar += (1.0f-rsq)*(1.0f-rsq)/snpi->sampleSize;
                    snpi->ldSum += v[j];
                    snpi->ldsc += rsq;
                    if (snpj == snpi)
                        ZPZdiag[inci] = v[j];
                }
            }
            SparseVector<float>::InnerIterator it(ZPZsp[inci]);
            windStart[inci] = snpi->windStart = it.index();
            windSize[inci] = snpi->windSize = ZPZsp[inci].nonZeros();
            for (; it; ++it) snpi->windEnd = it.index();
            
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
            snpi->ldsc = 0.0;
            
            for (unsigned j=0, incj=0; j<windSizeLDM[i]; ++j) {
                snpj = snpInfoVec[windStartLDM[i]+j];
                if (snpj->included) {
                    ZPZ[inci][incj++] = v[j];
                    rsq = v[j]*v[j];
                    snpi->ldSamplVar += (1.0f-rsq)*(1.0f-rsq)/snpi->sampleSize;
                    snpi->ldSum += v[j];
                    snpi->ldsc += rsq;
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
    cout << "LD matrix diagnal mean " << ZPZdiag.mean() << " sd " << sqrt(Gadget::calcVariance(ZPZdiag)) << "." << endl;
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
                snpj->ldsc = 0.0;

                for (unsigned k=0; k<windSizeLDM[j]; ++k) {
                    snpk = snpInfoVec[windStartLDM[j]+d[k]-d[0]];
                    if (snpk->included) {
                        ZPZsp[incj].insertBack(snpk->index) = v[k];
                        rsq = v[k]*v[k];
                        snpj->ldSamplVar += (1.0f-rsq)*(1.0f-rsq)/snpj->sampleSize;
                        snpj->ldSum += v[k];
                        snpj->ldsc += rsq;
                        if (snpk == snpj)
                            ZPZdiag[incj] = v[k];
                    }
                }
                SparseVector<float>::InnerIterator it(ZPZsp[incj]);
                windStart[incj] = snpj->windStart = it.index();
                windSize[incj] = snpj->windSize = ZPZsp[incj].nonZeros();
                ++incj;
            }
        }
        else {
            for (unsigned j=starti; j<numSnpMldVec[i]; j++) {
                snpj = snpInfoVec[j];
                float v[windSizeLDM[j]];
                
                if (!snpj->included) {
                    fseek(in2, sizeof(v), SEEK_CUR);
                    continue;
                }
                
                fread(v, sizeof(v), 1, in2);
                
                ZPZ[incj].resize(windSize[incj]);
                snpj->ldSamplVar = 0.0;
                snpj->ldSum = 0.0;
                snpj->ldsc = 0.0;

                for (unsigned k=0, inck=0; k<windSizeLDM[j]; ++k) {
                    snpk = snpInfoVec[windStartLDM[j]+k];
                    if (snpk->included) {
                        ZPZ[incj][inck++] = v[k];
                        rsq = v[k]*v[k];
                        snpj->ldSamplVar += (1.0f-rsq)*(1.0f-rsq)/snpj->sampleSize;
                        snpj->ldSum += v[k];
                        snpj->ldsc += rsq;
                        if (snpk == snpj)
                            ZPZdiag[incj] = v[k];
                    }
                }
                ++incj;
            }
        }
        
        fclose(in2);
        cout << "Read " + ldmType + " LD matrix for " << numSnpMldVec[i]-starti << " SNPs from [" << filenameVec[i] << "]." << endl;
        
        starti = numSnpMldVec[i];
    }
    
    timer.getTime();
    
    displayAverageWindowSize(windSize);
    cout << "LD matrix diagnal mean " << ZPZdiag.mean() << " sd " << sqrt(Gadget::calcVariance(ZPZdiag)) << "." << endl;
    cout << "Read LD matrix for " << numIncdSnps << " SNPs (time used: " << timer.format(timer.getElapse()) << ")." << endl;
}

void Data::resizeLDmatrix(const string &LDmatType, const float chisqThreshold, const unsigned windowWidth, const float LDthreshold) {
    if (LDmatType == "full") return;
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
                            if (ZPZ[i][j]*ZPZ[i][j]*snpi->sampleSize > chisqThreshold ||
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
                    for (unsigned i=0; i<numIncdSnps; ++i) {
                        snpi = incdSnpInfoVec[i];
                        ZPZsp[i].resize(snpi->windSize);
                        snpi->ldSamplVar = 0.0;
                        snpi->ldSum = 0.0;
                        for (unsigned j=0; j<snpi->windSize; ++j) {
                            snpj = incdSnpInfoVec[snpi->windStart + j];
                            if (ZPZ[i][j]*ZPZ[i][j]*snpi->sampleSize > chisqThreshold) {
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
                    if (rsq*snpi->sampleSize < chisqThreshold) it.valueRef() = 0.0;
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
                SnpInfo *snp = incdSnpInfoVec[i];
                unsigned windEndi = windSizeOri[i];
                for (unsigned j=0; j<windSizeOri[i]; ++j) {
                    if (abs(ZPZ[i][j]) > LDthreshold) {
                        windStart[i] = snp->windStart = windStartOri[i] + j;
                        break;
                    }
                }
                for (unsigned j=windSizeOri[i]; j>0; --j) {
                    if (abs(ZPZ[i][j-1]) > LDthreshold) {
                        windEndi = j;
                        break;
                    }
                }
                windSize[i] = snp->windSize = windStartOri[i] + windEndi - windStart[i];
                ZPZiTmp = ZPZ[i].segment(windStart[i] - windStartOri[i], windSize[i]);
                ZPZ[i] = ZPZiTmp;
                snp->ldSamplVar = (1.0 - ZPZ[i].array().square()).square().sum()/snp->sampleSize;
                snp->ldSum = ZPZ[i].sum();
            }
        } else {
            cout << "Resizing LD matrix based on a chisq threshold of " << chisqThreshold << "..." << endl;
            for (unsigned i=0; i<numIncdSnps; ++i) {
                SnpInfo *snp = incdSnpInfoVec[i];
                unsigned windEndi = windSizeOri[i];
                for (unsigned j=0; j<windSizeOri[i]; ++j) {
                    if (ZPZ[i][j]*ZPZ[i][j]*snp->sampleSize > chisqThreshold) {
                        windStart[i] = snp->windStart = windStartOri[i] + j;
                        break;
                    }
                }
                for (unsigned j=windSizeOri[i]; j>0; --j) {
                    if (ZPZ[i][j]*ZPZ[i][j]*snp->sampleSize > chisqThreshold) {
                        windEndi = j;
                        break;
                    }
                }
                windSize[i] = snp->windSize = windStartOri[i] + windEndi - windStart[i];
                ZPZiTmp = ZPZ[i].segment(windStart[i] - windStartOri[i], windSize[i]);
                ZPZ[i] = ZPZiTmp;
                snp->ldSamplVar = (1.0 - ZPZ[i].array().square()).square().sum()/snp->sampleSize;
                snp->ldSum = ZPZ[i].sum();
            }
        }
    }
    displayAverageWindowSize(windSize);
}

void Data::buildSparseMME(){
    VectorXf refZPZdiag = ZPZdiag;
    VectorXf refsnp2pq = snp2pq;
    
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
    }
    
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
            LDsamplVar[i]  = (snp->gwas_n + snp->sampleSize)/float(numIncdSnps)*snp->ldSamplVar;
            LDsamplVar[i] += (numIncdSnps - snp->windSize)/float(numIncdSnps);
            LDscore[i] = snp->ldsc*snp->gwas_n;
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
    
    
    ZPy = ZPZdiag.cwiseProduct(b);
    chisq = ZPy.cwiseProduct(b);
    
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
    //ypy = (D.array()*(n.array()*se.array().square()+b.array().square())).mean();
    VectorXf ypySrt = D.array()*(n.array()*se.array().square()+b.array().square());
    VectorXf varpSrt = ypySrt.array()/n.array();
    std::sort(ypySrt.data(), ypySrt.data() + ypySrt.size());
    std::sort(varpSrt.data(), varpSrt.data() + varpSrt.size());
    ypy = ypySrt[ypySrt.size()/2];  // median
    float varp = varpSrt[varpSrt.size()/2];

    //numKeptInds = n.mean();
    
    VectorXf nSrt = n;
    std::sort(nSrt.data(), nSrt.data() + nSrt.size());
    numKeptInds = nSrt[nSrt.size()/2]; // median
    
    //cout << ZPZ.size() << " " << ZPy.size() << " " << ypy << endl;
    //    cout << ZPy << endl;
    //    for (unsigned i=0; i<numIncdSnps; ++i) {
    //        cout << D[i] << "\t" << ZPZdiag[i] << endl;
    //    }
    //
    //    cout << ZPZ << endl;
    
    // no fixed effects
    numFixedEffects = 0;
    fixedEffectNames.resize(0);
    XPX.resize(0,0);
    ZPX.resize(0,0);
    XPy.resize(0);
    
    // data summary
    cout << "\nData summary:" << endl;
    cout << boost::format("%40s %8s %8s\n") %"" %"mean" %"sd";
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP heterozygosity" %snp2pq.mean() %sqrt(Gadget::calcVariance(snp2pq));
    cout << boost::format("%40s %8.0f %8.0f\n") %"GWAS SNP sample size" %n.mean() %sqrt(Gadget::calcVariance(n));
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP effect" %b.mean() %sqrt(Gadget::calcVariance(b));
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP SE" %se.mean() %sqrt(Gadget::calcVariance(se));
    cout << boost::format("%40s %8.3f %8.3f\n") %"MME left-hand-side diagonals" %ZPZdiag.mean() %sqrt(Gadget::calcVariance(ZPZdiag));
    cout << boost::format("%40s %8.3f %8.3f\n") %"MME right-hand-side" %ZPy.mean() %sqrt(Gadget::calcVariance(ZPy));
    cout << boost::format("%40s %8.3f %8.3f\n") %"LD sampling variance" %LDsamplVar.mean() %sqrt(Gadget::calcVariance(LDsamplVar));
    cout << "\n  Median of per-SNP phenotypic variance: " << varp << endl;
    
//    ofstream out("tmp.txt");
//    out << "refZPZdiag\t gwasZPZdiag\t b\t ZPy\t refsnp2pq\t gwassnp2pq n" << endl;
//    for (unsigned i=0; i<numIncdSnps; ++i) {
//        out << refZPZdiag[i] << "\t" << ZPZdiag[i] << "\t" << b[i] << "\t" << ZPy[i] << "\t" << refsnp2pq[i] << "\t" << snp2pq[i] << "\t" << n[i] << endl;
//    }
//    out.close();
}


void Data::outputSnpEffectSamples(const SparseMatrix<float> &snpEffects, const unsigned burnin, const unsigned outputFreq, const string&snpResFile, const string &filename) const {
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


