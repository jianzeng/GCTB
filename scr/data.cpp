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
    ifstream in(bedFile.c_str(), ios::binary);
    if (!in) throw ("Error: can not open the file [" + bedFile + "] to read.");
    if (myMPI::rank==0)
        cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
    char header[3];
    in.read((char *) header, 3);
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
    int genoValue;
    for (j = 0, snp = 0; j < numSnps; j++) {  // code adopted from BOLT-LMM with modification
        snpInfo = snpInfoVec[j];
        sum = 0.0;
        nmiss = 0;
        
        unsigned size = (numInds+3)>>2;
        
        if (!snpInfo->included) {
            in.ignore(size);
            continue;
        }
 
        char *bedLineIn = new char[size];
        in.read((char *)bedLineIn, size);

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
        snp2pq[snp] = 2.0f*snpInfo->af*(1.0f-snpInfo->af);
        
        //cout << "snp " << snp << "     " << Z.col(snp).sum() << endl;

        if (++snp == numIncdSnps) break;
    }
    in.close();
    
    
    // standardize genotypes
    VectorXf colsums = Z.colwise().sum();
    VectorXf colsums_all(numIncdSnps);
    
    MPI_Allreduce(&colsums[0], &colsums_all[0], numIncdSnps, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    
    Z.rowwise() -= colsums_all.transpose()/numKeptInds_all;  // center
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
            X.row(i) = keptIndInfoVec[i]->covariates;
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
        snp2pq[i] = 2.0f*snp->af*(1.0f-snp->af);
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
    }
}

void Data::getNonoverlapWindowInfo(const unsigned windowWidth){
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
            snp2pq[inc] = 2.0f*snpj->af*(1.0f-snpj->af);
            
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
    out << boost::format("%6s %20s %6s %12s %8s %12s %12s %8s %8s\n")
    % "Id"
    % "Name"
    % "Chrom"
    % "Position"
    % "GeneFrq"
    % "Effect"
    % "SE"
    % "PIP"
    % "Window";
    for (unsigned i=0, idx=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        if(!fullSnpFlag[i]) continue;
        if(snp->isQTL) continue;
        out << boost::format("%6s %20s %6s %12s %8.6f %12.6f %12.6f %8.3f %8s\n")
        % (idx+1)
        % snp->ID
        % snp->chrom
        % snp->physPos
        % snp->af
        % posteriorMean[idx]
        % sqrt(posteriorSqrMean[idx]-posteriorMean[idx]*posteriorMean[idx])
        % pip[idx]
        % snp->window;
        ++idx;
    }
    out.close();
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
    out << boost::format("%6s %20s %6s %12s %8s %12s %8s %8s\n")
    % "Id"
    % "Name"
    % "Chrom"
    % "Position"
    % "GeneFrq"
    % "Effect"
    % "PIP"
    % "Window";
    for (unsigned i=0, idx=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        if(!fullSnpFlag[i]) continue;
        out << boost::format("%6s %20s %6s %12s %8.3f %12.6f %8.3f %8s\n")
        % (idx+1)
        % snp->ID
        % snp->chrom
        % snp->physPos
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

void Data::makeLDmatrix(const string &bedFile, const unsigned windowWidth, const string &filename){

    cout << "Building LD matrix ..." << endl;
    
    getWindowInfo(incdSnpInfoVec, windowWidth, windStart, windSize);
    
    if (numIncdSnps == 0) throw ("Error: No SNP is retained for analysis.");
    if (numKeptInds == 0) throw ("Error: No individual is retained for analysis.");
    
    ZPZ.resize(numIncdSnps);
    ZPZdiag.resize(numIncdSnps);
    for (unsigned i=0; i<numIncdSnps; ++i) {
        ZPZ[i].resize(windSize[i]);
    }
    D.setZero(numIncdSnps);
    
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
        
        cout << "thread " << omp_get_thread_num() << " of " << omp_get_num_threads() << " is processing chrom " << chrinfo->id << " (" << chrinfo->size << " SNPs)." << endl;
        
        //cout << "chrom " << chrinfo->id << " start " << start << " end " << end << " lastWindStart " << lastWindStart << endl;
        
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
            
            //if (snpj->ID=="rs10165221") cout << snpj->ID << " af " << snpj->af << " z " << genotypes.head(10).transpose() << endl;

            // compute allele frequency
            snpj->af = 0.5f*mean;
            snp2pq[inc] = 2.0f*snpj->af*(1.0f-snpj->af);
            
            if (snp2pq[inc]==0) throw ("Error: " + snpj->ID + " is a fixed SNP!");
            
            // standardize genotypes
            D[inc] = snp2pq[inc]*(numKeptInds-nmiss);
            genotypes.array() -= genotypes.mean();
            genotypes.array() /= sqrtf(D[inc]);
            snpj->genotypes = genotypes;
            
            // compute Zj'Z[j] with Z[j] for genotype matrix of SNPs in the window of SNP j
            ZPZdiag[inc] = ZPZ[inc][inc - snpj->windStart] = genotypes.squaredNorm();
            for (k = snpj->windStart; k<inc; ++k) {
                snpk = incdSnpInfoVec[k];
                ZPZ[inc][k - snpj->windStart] = ZPZ[k][inc - snpk->windStart] = genotypes.dot(snpk->genotypes);
            }
            
            //if (snpj->ID=="rs10165221") cout << snpj->ID << " af " << snpj->af << " z " << genotypes.head(10).transpose() << endl;
            
            // release memory for genotypes of anterior SNPs of the window
            if (lastWindStart != snpj->windStart) {
                for (k=lastWindStart; k<snpj->windStart; ++k) {
                    incdSnpInfoVec[k]->genotypes.resize(0);
                }
                lastWindStart = snpj->windStart;
            }
            
            if (++inc == end) break;
        }
        
        in.close();
    }
    
    timer.getTime();
    
    cout << "Average window size " << windSize.sum()/numIncdSnps << "." << endl;
    cout << "LD matrix diagonal mean " << ZPZdiag.mean() << " variance " << (ZPZdiag.array()-ZPZdiag.mean()).square().sum()/numIncdSnps << "." << endl;
    cout << "Genotype data for " << numKeptInds << " individuals and " << numIncdSnps << " SNPs are included from [" + bedFile + "]." << endl;
    cout << "Build of LD matrix completed (time used: " << timer.format(timer.getElapse()) << ")." << endl;
    
    
    string outfilename = filename + ".ldm.w" + to_string(static_cast<long long>(windowWidth/1e6)) + "mb";
    string outfile1 = outfilename + ".info";
    ofstream out1(outfile1.c_str());
    SnpInfo *snp;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snp = incdSnpInfoVec[i];
        out1 << boost::format("%6s %15s %6s %15s %6s %6s %6s %6s %6s\n")
        %snp->chrom
        %snp->ID
        %snp->genPos
        %snp->physPos
        %snp->a1
        %snp->a2
        %snp->index
        %snp->windStart
        %snp->windSize;
    }
    out1.close();
    
    
    string outfile2 = outfilename + ".bin";
    FILE *out2 = fopen(outfile2.c_str(), "wb");
    if (!out2) {
        throw("Error: cannot open file " + outfile2);
    }
    
//    float ww = float(windowWidth)/1e6f;
//    fwrite(&ww, sizeof(float), 1, out2);
    
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        fwrite(&ZPZ[i][0], sizeof(float), snp->windSize, out2);
        //cout << "chr " << snp->chrom << " SNP " << snp->ID << " ZPZ[" << i << "] "; //<< ZPZ[i].transpose() << endl;
        //if(i==51) cout << ZPZ[i].transpose() << endl;
//        for (unsigned j=0; j<snp->windSize; ++j) {
//            SnpInfo *snp2 = incdSnpInfoVec[snp->windStart+j];
//            cout << snp2->ID << " " << ZPZ[i][j] << " ";
//            fwrite(&ZPZ[i][j], sizeof(float), 1, out2);
//        }
        //cout << endl;
    }
    fclose(out2);
    
    cout << "Written the LD matrix into file [" << outfile1 << "] ..." << endl;
    cout << "Written the SNP info into file [" << outfile2 << "] ..." << endl;
}

void Data::makeLDmatrix(const string &bedFile, const float LDthreshold, const string &snpRange, const string &filename){
    
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
        cout << "Building LD matrix for all SNPs ..." << endl;
    else
        cout << "Building LD matrix for SNPs " << snpRange << " ..." << endl;
    
    if (numIncdSnps == 0) throw ("Error: No SNP is retained for analysis.");
    if (numKeptInds == 0) throw ("Error: No individual is retained for analysis.");
    if (start >= numIncdSnps) throw ("Error: Specified a SNP range of " + snpRange + " but " + to_string(numIncdSnps) + " SNPs are included.");

    Gadget::Timer timer;
    timer.setTime();
    
    
    const int bedToGeno[4] = {2, -9, 1, 0};
    
    // first read in the genotypes of SNPs in the given range
    
    MatrixXf ZP(numSnpInRange, numKeptInds);  // SNP x Ind
    D.setZero(numSnpInRange);
    ifstream in(bedFile.c_str(), ios::binary);
    if (!in) throw ("Error: can not open the file [" + bedFile + "] to read.");
    cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
    char header[3];
    in.read((char *) header, 3);
    if (!in || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
        cerr << "Error: Incorrect first three bytes of bed file: " << bedFile << endl;
        exit(1);
    }
    
    IndInfo *indi = NULL;
    SnpInfo *snpj = NULL;
    
    int genoValue;
    unsigned i, j;
    unsigned inc; // index of included SNP
    
    for (j = 0, inc = 0; j < numSnps; j++) {
        
        unsigned size = (numInds+3)>>2;
        
        snpj = snpInfoVec[j];
        
        if (snpj->index < start || !snpj->included) {
            in.ignore(size);
            continue;
        }
        
        char *bedLineIn = new char[size];
        in.read((char *)bedLineIn, size);
        
        float mean = 0.0;
        unsigned nmiss = 0;
        
        for (i = 0; i < numInds; i++) {
            indi = indInfoVec[i];
            if (!indi->kept) continue;
            genoValue = bedToGeno[(bedLineIn[i>>2]>>((i&3)<<1))&3];
            ZP(inc, indi->index) = genoValue;
            if (genoValue == -9) ++nmiss;   // missing genotype
            else mean += genoValue;
        }
        delete[] bedLineIn;
        
        // fill missing values with the mean
        mean /= float(numKeptInds-nmiss);
        if (nmiss) {
            for (i=0; i<numKeptInds; ++i) {
                if (ZP(inc, i) == -9) ZP(inc, i) = mean;
            }
        }
        
        // compute allele frequency
        snpj->af = 0.5f*mean;
        snp2pq[inc] = 2.0f*snpj->af*(1.0f-snpj->af);
        
        if (snp2pq[inc]==0) throw ("Error: " + snpj->ID + " is a fixed SNP!");
        
        // standardize genotypes
        D[inc] = snp2pq[inc]*(numKeptInds-nmiss);
        
        if (++inc == numSnpInRange) break;
    }
    
    in.close();
    
    ZP = ZP.colwise() - ZP.rowwise().mean();
    ZP = ZP.array().colwise() / D.cwiseSqrt().array();
    
    ZPZdiag = ZP.rowwise().squaredNorm();
    
    
    // then read in the bed file again to compute Z'Z
    
    in.open(bedFile.c_str(), ios::binary);
    in.seekg(0, ios::beg);
    in.read((char *) header, 3);
    
    MatrixXf denseZPZ(numSnpInRange, numIncdSnps);
    VectorXf genotypes(numKeptInds);
    D.setZero(numIncdSnps);
    
    for (j = 0, inc = 0; j < numSnps; j++) {
        
        unsigned size = (numInds+3)>>2;
        
        snpj = snpInfoVec[j];
        
        if (!snpj->included) {
            in.ignore(size);
            continue;
        }
        
        char *bedLineIn = new char[size];
        in.read((char *)bedLineIn, size);
        
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
        snp2pq[inc] = 2.0f*snpj->af*(1.0f-snpj->af);
        
        if (snp2pq[inc]==0) throw ("Error: " + snpj->ID + " is a fixed SNP!");
        
        // standardize genotypes
        D[inc] = snp2pq[inc]*(numKeptInds-nmiss);
        genotypes.array() -= genotypes.mean();
        genotypes.array() /= sqrtf(D[inc]);
        
        denseZPZ.col(inc) = ZP * genotypes;
        
        if(!(inc%1000) && myMPI::rank==0) cout << " read snp " << inc << "\r" << flush;

        ++inc;
    }
    
    in.close();

    
    // find out per-SNP window position
    ZPZ.resize(numSnpInRange);
    windStart.setZero(numSnpInRange);
    windSize.setZero(numSnpInRange);
    
//    cout << denseZPZ.block(0,0,10,10) << endl;
//    cout << LDthreshold << endl;
    
    for (int i=0; i<numSnpInRange; ++i) {
        SnpInfo *snp = incdSnpInfoVec[start+i];
        unsigned windEndi = numIncdSnps;
        for (int j=0; j<numIncdSnps; ++j) {
            if (abs(denseZPZ(i,j)) > LDthreshold) {
                windStart[i] = snp->windStart = j;
                break;
            }
        }
        for (int j=numIncdSnps; j>0; --j) {
            if (abs(denseZPZ(i,j-1)) > LDthreshold) {
                windEndi = j;
                break;
            }
        }
        windSize[i] = snp->windSize = windEndi - windStart[i];
        ZPZ[i].resize(windSize[i]);
        VectorXf::Map(&ZPZ[i][0], windSize[i]) = denseZPZ.row(i).segment(windStart[i], windSize[i]);
    }
    
//    cout << windStart.transpose() << endl;
//    cout << windSize.transpose() << endl;

    
    timer.getTime();
    
    
    cout << "Window size mean " << windSize.mean() << " sd " << unsigned(sqrt(Gadget::calcVariance(windSize.cast<float>()))) << "." << endl;
    cout << "LD matrix diagonal mean " << ZPZdiag.mean() << " variance " << Gadget::calcVariance(ZPZdiag) << "." << endl;
    cout << "Genotype data for " << numKeptInds << " individuals and " << numSnpInRange << " SNPs are included from [" + bedFile + "]." << endl;
    cout << "Build of LD matrix completed (time used: " << timer.format(timer.getElapse()) << ")." << endl;
    
    
    string outfilename = filename + ".ldm";
    if (!snpRange.empty()) outfilename += ".snp" + snpRange;
    string outfile1 = outfilename + ".info";
    string outfile2 = outfilename + ".bin";
//    string outfile3 = outfilename + ".txt";
    ofstream out1(outfile1.c_str());
    FILE *out2 = fopen(outfile2.c_str(), "wb");
//    ofstream out3(outfile3.c_str());
    if (!out2) {
        throw("Error: cannot open file " + outfile2);
    }
    SnpInfo *snp;
    for (unsigned i=0; i<numSnpInRange; ++i) {
        snp = incdSnpInfoVec[i+start];
        out1 << boost::format("%6s %15s %6s %15s %6s %6s %6s %6s %6s\n")
        %snp->chrom
        %snp->ID
        %snp->genPos
        %snp->physPos
        %snp->a1
        %snp->a2
        %snp->index
        %snp->windStart
        %snp->windSize;
        fwrite(&ZPZ[i][0], sizeof(float), snp->windSize, out2);
//        out3 << ZPZ[i].transpose() << endl;
    }
    out1.close();
    fclose(out2);
//    out3.close();
    
    cout << "Written the LD matrix into file [" << outfile1 << "]." << endl;
    cout << "Written the SNP info into file [" << outfile2 << "]." << endl;
//    cout << "Written the LD matrix into text file [" << outfile3 << "]." << endl;

}

void Data::readLDmatrixInfoFile(const string &ldmatrixFile){
    ifstream in(ldmatrixFile.c_str());
    if (!in) throw ("Error: can not open the file [" + ldmatrixFile + "] to read.");
    cout << "Reading SNP info from [" + ldmatrixFile + "]." << endl;
    //snpInfoVec.clear();
    //snpInfoMap.clear();
    string id, allele1, allele2;
    unsigned chr, physPos;
    float genPos;
    unsigned idx, windStart, windSize;
    while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2 >> idx >> windStart >> windSize) {
        SnpInfo *snp = new SnpInfo(idx, id, allele1, allele2, chr, genPos, physPos);
        snp->windStart = windStart;
        snp->windSize = windSize;
        snpInfoVec.push_back(snp);
        if (snpInfoMap.insert(pair<string, SnpInfo*>(id, snp)).second == false) {
            throw ("Error: Duplicate SNP ID found: \"" + id + "\".");
        }
    }
    in.close();
    numSnps = (unsigned) snpInfoVec.size();
    cout << numSnps << " SNPs to be included from [" + ldmatrixFile + "]." << endl;
}


void Data::resizeWindow(const vector<SnpInfo *> &incdSnpInfoVec, const VectorXi &windStartOri, const VectorXi &windSizeOri,
                        VectorXi &windStart, VectorXi &windSize){
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
    }
}

void Data::readLDmatrixBinFile(const string &ldmatrixFile){
    
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
//    float windowWidth;
//    fread(&windowWidth, sizeof(float), 1, in);
//
//    getWindowInfo(incdSnpInfoVec, windowWidth*1e6, windStart, windSize);
    
    resizeWindow(incdSnpInfoVec, windStartLDM, windSizeLDM, windStart, windSize);
    
    if (numIncdSnps == 0) throw ("Error: No SNP is retained for analysis.");
    
    cout << "Reading LD matrix from [" + ldmatrixFile + "]..." << endl;
    
    ZPZ.resize(numIncdSnps);
    ZPZdiag.resize(numIncdSnps);
    for (unsigned i=0; i<numIncdSnps; ++i) {
        ZPZ[i].resize(windSize[i]);
    }

    Gadget::Timer timer;
    timer.setTime();
    
    for (unsigned i=0, inci=0; i<numSnps; i++) {
        snpi = snpInfoVec[i];
        
        float v[windSizeLDM[i]];
        
        if (!snpi->included) {
            fseek(in, sizeof(v), SEEK_CUR);
            continue;
        }
        
        fread(v, sizeof(v), 1, in);

        for (unsigned j=0, incj=0; j<windSizeLDM[i]; ++j) {
            snpj = snpInfoVec[windStartLDM[i]+j];
            if (snpj->included) {
                ZPZ[inci][incj++] = v[j];
                if (snpj == snpi)
                    ZPZdiag[inci] = v[j];
            }
        }
        
        if (++inci == numIncdSnps) break;
    }
    
    fclose(in);
    
    timer.getTime();
    
//    cout << "Window width " << windowWidth << " Mb." << endl;
    cout << "Window size mean " << windSize.mean() << " sd " << unsigned(sqrt(Gadget::calcVariance(windSize.cast<float>()))) << "." << endl;
    cout << "LD matrix diagnal mean " << ZPZdiag.mean() << " sd " << sqrt(Gadget::calcVariance(ZPZdiag)) << "." << endl;
    cout << "Read LD matrix for " << numIncdSnps << " SNPs (time used: " << timer.format(timer.getElapse()) << ")." << endl;
}

void Data::buildSparseMME(){
    snp2pq.resize(numIncdSnps);
    D.resize(numIncdSnps);
    ZPZdiag.resize(numIncdSnps);
    ZPy.resize(numIncdSnps);
    b.resize(numIncdSnps);
    n.resize(numIncdSnps);
    se.resize(numIncdSnps);
    tss.resize(numIncdSnps);
    SnpInfo *snp;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snp = incdSnpInfoVec[i];
        snp->af = snp->gwas_af;
        snp2pq[i] = 2.0f*snp->gwas_af*(1.0f-snp->gwas_af);
        D[i] = snp2pq[i]*snp->gwas_n;
        b[i] = snp->gwas_b;
        n[i] = snp->gwas_n;
        se[i]= snp->gwas_se;
        tss[i] = D[i]*(n[i]*se[i]*se[i] + b[i]*b[i]);
    }
    
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snp = incdSnpInfoVec[i];
        for (unsigned j=0; j<snp->windSize; ++j) {
            ZPZ[i][j] *= sqrt(D[i]*D[snp->windStart+j]);
        }
        //ZPZdiag[i] = ZPZ[i][i-snp->windStart];
    }
    ZPZdiag = D;
    
    //b.array() -= b.mean();
    ZPy.array() = D.array()*b.array();
    
//    cout << "ZPZdiag " << ZPZdiag.transpose() << endl;
//    cout << "ZPZ.back() " << ZPZ.back().transpose() << endl;
//    cout << "ZPZ.front() " << ZPZ.front().transpose() << endl;
//    cout << "ZPy " << ZPy.head(100).transpose() << endl;
//    cout << "b.mean() " << b.mean() << endl;
    
    // estimate ypy
    ypy = (D.array()*(n.array()*se.array().square()+b.array().square())).mean();
    numKeptInds = n.mean();
    
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
    cout << boost::format("%30s %8s %8s\n") %"" %"mean" %"sd";
    cout << boost::format("%30s %8.3f %8.3f\n") %"GWAS SNP heterozygosity" %snp2pq.mean() %sqrt(Gadget::calcVariance(snp2pq));
    cout << boost::format("%30s %8.0f %8.0f\n") %"GWAS SNP sample size" %n.mean() %sqrt(Gadget::calcVariance(n));
    cout << boost::format("%30s %8.3f %8.3f\n") %"GWAS SNP effect" %b.mean() %sqrt(Gadget::calcVariance(b));
    cout << boost::format("%30s %8.3f %8.3f\n") %"GWAS SNP SE" %se.mean() %sqrt(Gadget::calcVariance(se));
    cout << boost::format("%30s %8.3f %8.3f\n") %"MME left-hand-side diagonals" %ZPZdiag.mean() %sqrt(Gadget::calcVariance(ZPZdiag));
    cout << boost::format("%30s %8.3f %8.3f\n") %"MME right-hand-side" %ZPy.mean() %sqrt(Gadget::calcVariance(ZPy));
    cout << boost::format("%30s %8.3f\n") %"Total sum of square" %ypy;
}

void Data::readMultiLDmatInfoFile(const string &mldmatFile){
    ifstream in(mldmatFile.c_str());
    if (!in) throw ("Error: can not open the file [" + mldmatFile + "] to read.");
    cout << "Reading SNP info from [" + mldmatFile + "]..." << endl;
    string inputStr;
    numSnpMldVec.clear();
    while (getline(in, inputStr)) {
//        vector<SnpInfo*> vec;
        readLDmatrixInfoFile(inputStr+".info");
        numSnpMldVec.push_back(numSnps);
//        mldmVec.push_back(vec);
    }
}

void Data::readMultiLDmatBinFile(const string &mldmatFile){
    vector<string> filenameVec;
    ifstream in1(mldmatFile.c_str());
    if (!in1) throw ("Error: can not open the file [" + mldmatFile + "] to read.");
    cout << "Reading LD matrices from [" + mldmatFile + "]..." << endl;
    
    Gadget::Timer timer;
    timer.setTime();
    
    string inputStr;
    while (getline(in1, inputStr)) {
        filenameVec.push_back(inputStr + ".bin");
    }
    
    VectorXi windStartLDM(numSnps);
    VectorXi windSizeLDM(numSnps);
    
    for (unsigned j=0, i=0, cnt=0; j<numSnps; ++j) {
        SnpInfo *snp = snpInfoVec[j];
        if (j==numSnpMldVec[i]) {
            if (snp->windStart == 0)
                cnt = numSnpMldVec[i++];
            else
                cnt = 0;
        }
        windStartLDM[j] = snp->windStart + cnt;
        windSizeLDM[j] = snp->windSize + cnt;
    }
    
    resizeWindow(incdSnpInfoVec, windStartLDM, windSizeLDM, windStart, windSize);
    
    ZPZ.resize(numIncdSnps);
    ZPZdiag.resize(numIncdSnps);
    for (unsigned j=0; j<numIncdSnps; ++j) {
        ZPZ[j].resize(windSize[j]);
    }
    
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
        
        for (unsigned j=starti; j<numSnpMldVec[i]; j++) {
            snpj = snpInfoVec[j];
            float v[windSizeLDM[j]];
            
            if (!snpj->included) {
                fseek(in2, sizeof(v), SEEK_CUR);
                continue;
            }
            
            fread(v, sizeof(v), 1, in2);
            
            for (unsigned k=0, inck=0; k<windSizeLDM[j]; ++k) {
                snpk = snpInfoVec[windStartLDM[j]+k];
                if (snpk->included) {
                    ZPZ[incj][inck++] = v[k];
                    if (snpk == snpj)
                        ZPZdiag[incj] = v[k];
                }
            }
            ++incj;
        }
        
        fclose(in2);
        cout << "Read LD matrix for " << numSnpMldVec[i]-starti << " SNPs from [" << filenameVec[i] << "]." << endl;

        starti = numSnpMldVec[i];
}
    
    timer.getTime();
    
    cout << "Window size mean " << windSize.mean() << " sd " << unsigned(sqrt(Gadget::calcVariance(windSize.cast<float>()))) << "." << endl;
    cout << "LD matrix diagnal mean " << ZPZdiag.mean() << " sd " << sqrt(Gadget::calcVariance(ZPZdiag)) << "." << endl;
    cout << "Read LD matrix for " << numIncdSnps << " SNPs (time used: " << timer.format(timer.getElapse()) << ")." << endl;
}


void Data::outputLDmatrix(const string &filename) const {
    string outfilename = filename + ".ldm";
    string outfile1 = outfilename + ".info";
    string outfile2 = outfilename + ".bin";
    ofstream out1(outfile1.c_str());
    FILE *out2 = fopen(outfile2.c_str(), "wb");
    if (!out2) {
        throw("Error: cannot open file " + outfile2);
    }
    SnpInfo *snp;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snp = incdSnpInfoVec[i];
        out1 << boost::format("%6s %15s %6s %15s %6s %6s %6s %6s %6s\n")
        %snp->chrom
        %snp->ID
        %snp->genPos
        %snp->physPos
        %snp->a1
        %snp->a2
        %snp->index
        %snp->windStart
        %snp->windSize;
        fwrite(&ZPZ[i][0], sizeof(float), snp->windSize, out2);
    }
    out1.close();
    fclose(out2);
    
    cout << "Written the LD matrix into file [" << outfile1 << "]." << endl;
    cout << "Written the SNP info into file [" << outfile2 << "]." << endl;
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
