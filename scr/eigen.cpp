//
//  eigen.cpp
//  gctb
//
//  Created by Shouye Liu and Jian Zeng on 05/01/2023.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "data.hpp"

///////////////////////////////////////////////////////////////////////////////////////
////////    Step 1. perform eigen-decomposition for ld blocks                   ///////
///////////////////////////////////////////////////////////////////////////////////////


void Data::readLDBlockInfoFile(const string &ldBlockInfoFile){
    // Read bim file: recombination rate is defined between SNP i and SNP i-1
    ifstream in(ldBlockInfoFile.c_str());
    if (!in) throw ("Error: can not open the file [" + ldBlockInfoFile + "] to read.");
    cout << "Reading ld block info from file [" + ldBlockInfoFile + "]." << endl;
    ldBlockInfoVec.clear();
    ldBlockInfoMap.clear();
    string header;
    string id;
    int  chr,start, stop,pdist;
    float start_cm, stop_cm,gdist;
    int idx = 0;
    getline(in, header);
    while (in >>chr>>id>>start>>stop>>start_cm>>stop_cm>>pdist>>gdist) {
        LDBlockInfo *ld = new LDBlockInfo(idx++, id, chr);
        ld->startPos    = start;
        ld->endPos     = stop;
        ld->start_cm    = start_cm;
        ld->stop_cm     = stop_cm;
        ld->pdist       = pdist;
        ld->gdist       = gdist;
        ldBlockInfoVec.push_back(ld);
        //chromosomes.insert(ld->chr);
        if (ldBlockInfoMap.insert(pair<string, LDBlockInfo*>(id, ld)).second == false) {
            throw ("Error: Duplicate LD block ID found: \"" + id + "\".");
        }
    }
    in.close();
    numLDBlocks = (unsigned) ldBlockInfoVec.size();
    cout << numLDBlocks << " LD Blocks to be included from [" + ldBlockInfoFile + "]." << endl;
}

void Data::eigenDecomposition( const MatrixXf &X, const float &prop, VectorXf &eigenValAdjusted, MatrixXf &eigenVecAdjusted, VectorXf &cumsumNonNeg){
    // VectorXf cumsumNonNeg; // cumulative sums of non-negative values
    float sumNonNeg = 0.0;

    SelfAdjointEigenSolver<MatrixXf> eigensolver(X);
    VectorXf eigenVal = eigensolver.eigenvalues();
    MatrixXf eigenVec = eigensolver.eigenvectors();
    int revIdx = eigenVal.size();
    cumsumNonNeg.resize(revIdx);
    cumsumNonNeg.setZero();
    revIdx = revIdx -1;
    if(eigenVal(revIdx) < 0) cout << "Error, all eigenvector are negative" << endl;
    cumsumNonNeg(revIdx) = eigenVal(revIdx);
    sumNonNeg = eigenVal(revIdx);
    revIdx = revIdx -1;
    
    while( eigenVal(revIdx) > 1e-10 ){
        sumNonNeg = sumNonNeg + eigenVal(revIdx);
        cumsumNonNeg(revIdx) = eigenVal(revIdx) + cumsumNonNeg(revIdx + 1);
        revIdx =revIdx - 1;
        if(revIdx < 0) break;
    }
    // cout << "revIdx: " << revIdx << endl;
    // cout << "size: " << eigenVal.size()  << " eigenVal: " << eigenVal << endl;
    // cout << "cumsumNoNeg: " << cumsumNonNeg << endl;
    cumsumNonNeg = cumsumNonNeg/sumNonNeg;
    bool haveValue = false;
    // cout << "cumsumNonNeg: " << cumsumNonNeg << endl;
    // cout << "revIdx : " << revIdx << endl;
    // cout << "revIdx: "  << revIdx  << endl;
    for (revIdx = revIdx + 1; revIdx < eigenVal.size(); revIdx ++ ){
        // cout << "revIdx: " << revIdx << " cumsumNonNeg: " << cumsumNonNeg(revIdx) << endl;
        if(prop >= cumsumNonNeg(revIdx) ){
            revIdx = revIdx -1;
            haveValue = true;
            break;
        }
    }
    // cout << "revIdx : " << revIdx << endl;
    if(!haveValue) revIdx = eigenVal.size() - 1;
    // cout << "cumsumNonNeg: " << cumsumNonNeg.size() << endl;
    // cout << "cumsumNoNeg: " << cumsumNonNeg << endl;
    // cout << "revIdx: "  << revIdx  << endl;

    eigenVecAdjusted = eigenVec.rightCols(eigenVal.size() - revIdx);
    eigenValAdjusted = eigenVal.tail(eigenVal.size() - revIdx);
    // cout << "eigenValAdjusted size: " << eigenValAdjusted.size() << " eigenValue eventually: " << eigenValAdjusted << endl;
    //eigenvalueNum = eigenVal.size() - revIdx;
    // cout << endl;
}

MatrixXf Data::generateLDmatrixPerBlock(const string &bedFile, const vector<string> &snplists){
    int numSnpInRange = snplists.size();
    IndInfo *indi = NULL;
    SnpInfo *snpj = NULL;
    SnpInfo *snpk = NULL;

    snpj = snpInfoMap.at(snplists[0]); // start
    snpk = snpInfoMap.at(snplists[numSnpInRange - 1]); // end;
    unsigned start = snpj->index;
    unsigned end = snpk->index;
    
    if (numIncdSnps == 0) throw ("Error: No SNP is retained for analysis.");
    if (numKeptInds == 0) throw ("Error: No individual is retained for analysis.");
    // if (start >= numIncdSnps) throw ("Error: Specified a SNP range of " + snpRange + " but " + to_string(static_cast<long long>(numIncdSnps)) + " SNPs are included.");
    
    // Gadget::Timer timer;
    // timer.setTime();
    
    //////////////////////////////////////////////////////
    // Step 1. read in the genotypes of SNPs in the given range
    //////////////////////////////////////////////////////
    
    const int bedToGeno[4] = {2, -9, 1, 0};
    unsigned size = (numInds+3)>>2;
    
    MatrixXf ZP(numSnpInRange, numKeptInds);  // SNP x Ind
    VectorXf Dtmp;
    Dtmp.setZero(numSnpInRange);

    if (numKeptInds < 2) throw("Error: Cannot calculate LD matrix with number of individuals < 2.");
    
    FILE *in1 = fopen(bedFile.c_str(), "rb");
    if (!in1) throw ("Error: can not open the file [" + bedFile + "] to read.");
    // cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
    char header[3];
    fread(header, sizeof(header), 1, in1);
    if (!in1 || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
        cerr << "Error: Incorrect first three bytes of bed file: " << bedFile << endl;
        exit(1);
    }

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
        
        Dtmp[incj] = Gadget::calcVariance(ZP.row(incj))*numKeptInds;
        
        ZP.row(incj) = (ZP.row(incj).array() - ZP.row(incj).mean())/sqrt(Dtmp[incj]);
        
        if (++incj == numSnpInRange) break;
    }
    
    fclose(in1);
    
    ZPZdiag = ZP.rowwise().squaredNorm();
    
    //////////////////////////////////////////////////////
    // Step 2. read in the bed file again to compute Z'Z
    //////////////////////////////////////////////////////
    
    MatrixXf denseZPZ;
    denseZPZ.setZero(numSnpInRange, numSnpInRange);
    VectorXf Zk(numKeptInds);
    Dtmp.setZero(numIncdSnps);
    
    FILE *in2 = fopen(bedFile.c_str(), "rb");
    fseek(in2, 3, SEEK_SET);
    unsigned long long skipk = 0;
    
    set<int>::iterator setend = chromInRange.end();

    if (numSkeletonSnps) {
        for (k = 0, inck = 0; k < numSnps; k++) {
            snpk = snpInfoVec[k];
            
            // if (!snpk->included) {
            if (snpk->index < start || !snpk->included) {
                skipk += size;
                continue;
            }

            if (chromInRange.find(snpk->chrom) == setend && !snpk->skeleton) {
                skipk += size;
                ++inck;       // ensure the index is correct
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
            snp2pq[inck] = snpk->twopq = 2.0f*snpk->af*(1.0f-snpk->af);
            
            if (snp2pq[inck]==0) throw ("Error: " + snpk->ID + " is a fixed SNP!");
            Dtmp[inck] = Gadget::calcVariance(Zk.row(inck))*numKeptInds;
            Zk = (Zk.array() - Zk.mean())/sqrt(Dtmp[inck]);
            denseZPZ.col(inck) = ZP * Zk;

            //++inck;
            if (++inck == numSnpInRange) break;
        }
    }
    else {
        for (k = 0, inck = 0; k < numSnps; k++) {
            snpk = snpInfoVec[k];
            
            // if (!snpk->included) {
            if (snpk->index < start || !snpk->included) {
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
            snp2pq[inck] = snpk->twopq = 2.0f*snpk->af*(1.0f-snpk->af);
            
            if (snp2pq[inck]==0) throw ("Error: " + snpk->ID + " is a fixed SNP!");

            Dtmp[inck] = Gadget::calcVariance(Zk)*numKeptInds;

            Zk = (Zk.array() - Zk.mean())/sqrt(Dtmp[inck]);
            
            denseZPZ.col(inck) = ZP * Zk;

            if (++inck == numSnpInRange) break;
        }
    }

    fclose(in2);
    // timer.getTime();
    return denseZPZ;
}

void Data::getEigenDataFromFullLDM(const string &filename, const float eigenCutoff = 0.995){
    
    string outfilename = filename + ".eigen.bin";
    FILE *out3 = fopen(outfilename.c_str(), "wb");

    for (unsigned blk=0; blk < numKeptLDBlocks; blk++){
        
        MatrixXf eigenVec;
        VectorXf eigenVal, cumsumNonNegPerLD;
        // cout << "rval: " << rval << endl;
        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        eigenDecomposition(ZPZmat, eigenCutoff, eigenVal, eigenVec, cumsumNonNegPerLD);
        // cout << "rval: " << rval.row(0) << endl;
        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
        // save svd matrix
        
        int32_t numEigenValue = eigenVal.size();
        int32_t numSnpInBlock = keptLdBlockInfoVec[blk]->numSnpInBlock;  // TMP
        float eigenValueSum = eigenVal.sum();

        cout << " Generate Eigen decomposition result for LD block " << blk << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << endl;
        
        // save summary
        // 1, the number of SNPs in the block
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, out3);
        // 2, the number of eigenvalues at with the given cutoff
        fwrite(&numEigenValue, sizeof(int32_t), 1, out3);
        // 3. sum of the selected eigenvalues
        fwrite(&eigenValueSum, sizeof(float), 1, out3);
        //4. eigenvalue cutoff based on the proportion of variance explained in LD
        fwrite(&eigenCutoff, sizeof(float), 1, out3);
        // 5. save the selected eigenvalues
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, out3);
        // 6. save a series of proportions of variance in LD, please note, the order of
        // these proportions are the opposite of the order of eigenValues
        fwrite(cumsumNonNegPerLD.data(), sizeof(float), numSnpInBlock, out3);
        // 7. save eigen vector;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, out3);
        
    }
    
    fclose(out3);

    //outputEigenDataForLDM(filename);
    
}

void Data::makeFullLdmForLdBlocks(const string &bedFile, const string &ldBlockInfoFile, const string &filename, const bool writeLdmTxt, int ldBlockRegionWind){
    int i,j;
    vector<locus_bp> snpVec;
    SnpInfo *snp;
    
    map<int, string>  chrEndSnp;
    for (i = 1; i < numIncdSnps; i++) {
        snp = incdSnpInfoVec[i];
        if(incdSnpInfoVec[i]->chrom != incdSnpInfoVec[i-1]->chrom){
            chrEndSnp.insert(pair<int, string>(incdSnpInfoVec[i - 1]->chrom,incdSnpInfoVec[i - 1]->ID ));
        }
    }
    chrEndSnp.insert(pair<int, string>(incdSnpInfoVec[numIncdSnps - 1]->chrom,incdSnpInfoVec[numIncdSnps - 1]->ID ));
    //Step 1.2  Read block file
    readLDBlockInfoFile(ldBlockInfoFile);
    /////////////////////////////////////////
    // Step 2. Map snps to blocks
    /////////////////////////////////////////
    vector<string> block2snp_1(numLDBlocks), block2snp_2(numLDBlocks);
    map<string,int> keptLdBlock2AllLdBlcokMap;
    vector<locus_bp>::iterator iter;
    map<int, string>::iterator chrIter;
    LDBlockInfo *ldblock;
    for (i = 0; i < numIncdSnps ; i++) {
        snp = incdSnpInfoVec[i];
        snpVec.push_back(locus_bp(snp->ID, snp->chrom, snp->physPos ));
    }
#pragma omp parallel for private(iter, chrIter)
    for (i = 0; i < numLDBlocks; i++) {
        // find lowest snp_name in the block
        ldblock = ldBlockInfoVec[i];
        
        iter = find_if(snpVec.begin(), snpVec.end(), locus_bp( ldblock->ID ,ldblock->chrom, ldblock->startPos - ldBlockRegionWind));
        if (iter != snpVec.end()) block2snp_1[i] = iter->locusName;
        else block2snp_1[i] = "NA";
    }
#pragma omp parallel for private(iter, chrIter)
    for (i = 0; i < numLDBlocks; i++) {
        ldblock = ldBlockInfoVec[i];
        if (block2snp_1[i] == "NA") {
            block2snp_2[i] = "NA";
            continue;
        }
        iter = find_if(snpVec.begin(), snpVec.end(), locus_bp(ldblock->ID, ldblock->chrom, ldblock->endPos + ldBlockRegionWind));
        if (iter != snpVec.end()){
            if (iter->bp ==  ldblock->endPos + ldBlockRegionWind){
                block2snp_2[i] = iter->locusName;
            }else {
                if(iter!=snpVec.begin()){
                    iter--;
                    block2snp_2[i] = iter->locusName;
                }
                else block2snp_2[i] = "NA";
            }
        }
        else {
            chrIter = chrEndSnp.find(ldblock->chrom);
            if (chrIter == chrEndSnp.end()) block2snp_2[i] = "NA";
            else block2snp_2[i] = chrIter->second;
        }
    }
    int mapped = 0;
    for (i = 0; i < numLDBlocks; i++) {
        ldblock = ldBlockInfoVec[i];
        if (block2snp_1[i] != "NA" && block2snp_2[i] != "NA")
        {
            mapped++;
            // ldblock->kept = true;
            keptLdBlock2AllLdBlcokMap.insert(pair<string, int>(ldblock->ID, i));
        } else {
            ldblock->kept = false;
        }
    }
    if (mapped < 1) throw(0, "No SNP can be mapped to the provided LD block list. Please check the input data regarding chromosome and bp.");
    else cout << mapped << " LD blocks have at least one SNP." << endl;
    
    keptLdBlockInfoVec = makeKeptLDBlockInfoVec(ldBlockInfoVec);
    numKeptLDBlocks = (unsigned) keptLdBlockInfoVec.size();
    
    map<string, int>::iterator iter1, iter2;
    map<string, int> snp_name_map;
    VectorXf snpNumInldblock(numKeptLDBlocks);
    VectorXf eigenNumInldblock(numKeptLDBlocks);
    vector<VectorXf> cumsumNonNeg(numKeptLDBlocks);
    
    for (i = 0; i < numIncdSnps; i++) {
        snp = incdSnpInfoVec[i];
        snp_name_map.insert(pair<string,int>(snp->ID, i));
        
    }
    // eigenvector and eigenvalue
    //eigenValLdBlock.resize(numKeptLDBlocks);
    //eigenVecLdBlock.resize(numKeptLDBlocks);
    string LDmatType = "full";
    //    string outfilename = filename + "." + LDmatType + ".eigen" + ".bin";
    //    FILE *out = fopen(outfilename.c_str(), "wb");
    
    string outfilename = filename + ".ldm." + LDmatType;
    //    string outfile1 = outfilename + ".info";
    string outfile2 = outfilename + ".bin";
    //    ofstream out1(outfile1.c_str());
    FILE *out2 = fopen(outfile2.c_str(), "wb");
    ofstream out3;
    string outfile3;
    if (writeLdmTxt) {
        outfile3 = outfilename + ".txt";
        out3.open(outfile3.c_str());
    }
    
    bool readBedBool = true;
    for (i = 0; i < numKeptLDBlocks; i++) {
        ldblock = keptLdBlockInfoVec[i];
        // cout << "ldblock id: " << ldblock->ID << endl;
        iter1 = snp_name_map.find(block2snp_1[keptLdBlock2AllLdBlcokMap.at(ldblock->ID ) ]);
        iter2 = snp_name_map.find(block2snp_2[keptLdBlock2AllLdBlcokMap.at(ldblock->ID ) ]);
        bool skip = false;
        if (iter1 == snp_name_map.end() || iter2 == snp_name_map.end() || iter1->second >= iter2->second) ldblock->kept = false;
        snpNumInldblock[i] = iter2->second - iter1->second + 1;
        // cout << "ldblock->kept: " << ldblock->kept << endl;
        if(!ldblock->kept) continue;
        vector<int> snp_indx;
        for (j = iter1->second; j <= iter2->second; j++) {
            snp_indx.push_back(j);
            ldblock->gwasSnpNameVecInBlock.push_back(incdSnpInfoVec[j]->ID);
            // cout << incdSnpInfoVec[j]->ID << " ";
        }
        if(readBedBool) {
            cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
            readBedBool = false;
        }
        //        MatrixXf eigenVec;
        //        VectorXf eigenVal,cumsumNonNegPerLD;
        MatrixXf rval = generateLDmatrixPerBlock(bedFile, ldblock->gwasSnpNameVecInBlock);
        //        // cout << "rval: " << rval << endl;
        //        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        //        eigenDecomposition(rval, eigenCutoff,eigenVal, eigenVec,cumsumNonNegPerLD);
        //        // cout << "rval: " << rval.row(0) << endl;
        //        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
        //        // save svd matrix
        //        int32_t numEigenValue = eigenVal.size();
                int32_t numSnpInBlock = ldblock->gwasSnpNameVecInBlock.size();
        //        eigenNumInldblock[i] = numEigenValue;
        //        float eigenValueSum = eigenVal.sum();
        //        // save summary
        //        // 1, nrow of eigenVecGene[i]
        //        fwrite(&numSnpInBlock, sizeof(int32_t), 1, out3);
        //        //        cout << "rval: " << rval << endl;
        //        // cout << "eigenVec: "  << eigenVec << endl;
        //        // cout << "";
        //        // 2, ncol of eigenVecGene[i]
        //        fwrite(&numEigenValue, sizeof(int32_t), 1, out3);
        //        // 3. sum of eigen values
        //        fwrite(&eigenValueSum, sizeof(float), 1, out3);
        //        //4. eigenCutoff
        //        fwrite(&eigenCutoff, sizeof(float), 1, out3);
        //        // 5. save eigen value
        //        fwrite(eigenVal.data(), sizeof(float), numEigenValue, out3);
        //        // 6. save a series of proportions of variance in LD, please note, the order of
        //        // these proportions are the opposite of the order of eigenValues
        //        fwrite(cumsumNonNegPerLD.data(), sizeof(float), numSnpInBlock, out3);
        //        // 7. save eigen vector;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numSnpInBlock;
        fwrite(rval.data(), sizeof(float), nElements, out2);
        //        cout << " Generate and save Eigen decomposition result for LD block " << i << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << "\r" << flush;
        
        if (writeLdmTxt) {
            for (unsigned ii=0; ii<numSnpInBlock; ++ii){
                for (unsigned jj=0; jj<numSnpInBlock; ++jj) {
                    out3 << ldblock->ID << "\t" << ldblock->gwasSnpNameVecInBlock[ii] << "\t" << ldblock->gwasSnpNameVecInBlock[jj] << "\t" << rval(ii,jj) << endl;
                }
            }
        }
    }
    
    fclose(out2);
    
    cout << "Written the LD matrix into file [" << outfile2 << "]." << endl;

    if (writeLdmTxt) {
        out3.close();
        cout << "Written the LD matrix into text file [" << outfile3 << "]." << endl;
    }
    
    // cout << "size of eigenValuegene: " << eigenValGene.size() << endl;
    outputLdDataForBlockLDM(filename);

    
}


void Data::getEigenDataForLDBlock(const string &bedFile, const string &ldBlockInfoFile, int ldBlockRegionWind, const string &filename, const float eigenCutoff){
    int i,j;
    vector<locus_bp> snpVec;
    SnpInfo *snp;

    map<int, string>  chrEndSnp;
    for (i = 1; i < numIncdSnps; i++) {
        snp = incdSnpInfoVec[i];
        if(incdSnpInfoVec[i]->chrom != incdSnpInfoVec[i-1]->chrom){
            chrEndSnp.insert(pair<int, string>(incdSnpInfoVec[i - 1]->chrom,incdSnpInfoVec[i - 1]->ID ));
        }
    }
    chrEndSnp.insert(pair<int, string>(incdSnpInfoVec[numIncdSnps - 1]->chrom,incdSnpInfoVec[numIncdSnps - 1]->ID ));
    //Step 1.2  Read block file
    readLDBlockInfoFile(ldBlockInfoFile);
    /////////////////////////////////////////
    // Step 2. Map snps to blocks
    /////////////////////////////////////////
    vector<string> block2snp_1(numLDBlocks), block2snp_2(numLDBlocks);
    map<string,int> keptLdBlock2AllLdBlcokMap;
    vector<locus_bp>::iterator iter;
    map<int, string>::iterator chrIter;
    LDBlockInfo *ldblock;
    for (i = 0; i < numIncdSnps ; i++) {
        snp = incdSnpInfoVec[i];
        snpVec.push_back(locus_bp(snp->ID, snp->chrom, snp->physPos ));
    }
#pragma omp parallel for private(iter, chrIter)
    for (i = 0; i < numLDBlocks; i++) {
        // find lowest snp_name in the block
        ldblock = ldBlockInfoVec[i];

        iter = find_if(snpVec.begin(), snpVec.end(), locus_bp( ldblock->ID ,ldblock->chrom, ldblock->startPos - ldBlockRegionWind));
        if (iter != snpVec.end()) block2snp_1[i] = iter->locusName;
        else block2snp_1[i] = "NA";
    }
#pragma omp parallel for private(iter, chrIter)
    for (i = 0; i < numLDBlocks; i++) {
        ldblock = ldBlockInfoVec[i];
        if (block2snp_1[i] == "NA") {
            block2snp_2[i] = "NA";
            continue;
        }
        iter = find_if(snpVec.begin(), snpVec.end(), locus_bp(ldblock->ID, ldblock->chrom, ldblock->endPos + ldBlockRegionWind));
        if (iter != snpVec.end()){
            if (iter->bp ==  ldblock->endPos + ldBlockRegionWind){
                block2snp_2[i] = iter->locusName;
            }else {
                if(iter!=snpVec.begin()){
                    iter--;
                    block2snp_2[i] = iter->locusName;
                }
                else block2snp_2[i] = "NA";
            }
        }
        else {
            chrIter = chrEndSnp.find(ldblock->chrom);
            if (chrIter == chrEndSnp.end()) block2snp_2[i] = "NA";
            else block2snp_2[i] = chrIter->second;
        }
    }
    int mapped = 0;
    for (i = 0; i < numLDBlocks; i++) {
        ldblock = ldBlockInfoVec[i];
        if (block2snp_1[i] != "NA" && block2snp_2[i] != "NA")
        {
            mapped++;
            // ldblock->kept = true;
            keptLdBlock2AllLdBlcokMap.insert(pair<string, int>(ldblock->ID, i));
        } else {
            ldblock->kept = false;
        }
    }
    if (mapped < 1) throw(0, "No SNP can be mapped to the provided LD block list. Please check the input data regarding chromosome and bp.");
    else cout << mapped << " LD blocks have at least one SNP." << endl;

    keptLdBlockInfoVec = makeKeptLDBlockInfoVec(ldBlockInfoVec);
    numKeptLDBlocks = (unsigned) keptLdBlockInfoVec.size();

    map<string, int>::iterator iter1, iter2;
    map<string, int> snp_name_map;
    VectorXf snpNumInldblock(numKeptLDBlocks);
    VectorXf eigenNumInldblock(numKeptLDBlocks);
    vector<VectorXf> cumsumNonNeg(numKeptLDBlocks);
    
    for (i = 0; i < numIncdSnps; i++) {
        snp = incdSnpInfoVec[i];
        snp_name_map.insert(pair<string,int>(snp->ID, i));
        
    }
      // eigenvector and eigenvalue
    //eigenValLdBlock.resize(numKeptLDBlocks);
    //eigenVecLdBlock.resize(numKeptLDBlocks);
    string LDmatType = "ldblock";
    string outfilename = filename + "." + LDmatType + ".eigen" + ".bin";
    FILE *out3 = fopen(outfilename.c_str(), "wb");

    bool readBedBool = true;
    for (i = 0; i < numKeptLDBlocks; i++) {
        ldblock = keptLdBlockInfoVec[i];
        // cout << "ldblock id: " << ldblock->ID << endl;
        iter1 = snp_name_map.find(block2snp_1[keptLdBlock2AllLdBlcokMap.at(ldblock->ID ) ]);
        iter2 = snp_name_map.find(block2snp_2[keptLdBlock2AllLdBlcokMap.at(ldblock->ID ) ]);
        bool skip = false;
        if (iter1 == snp_name_map.end() || iter2 == snp_name_map.end() || iter1->second >= iter2->second) ldblock->kept = false;
        snpNumInldblock[i] = iter2->second - iter1->second + 1;
        // cout << "ldblock->kept: " << ldblock->kept << endl;
        if(!ldblock->kept) continue;
        vector<int> snp_indx;
        for (j = iter1->second; j <= iter2->second; j++) {
            snp_indx.push_back(j);
            ldblock->gwasSnpNameVecInBlock.push_back(incdSnpInfoVec[j]->ID);
           // cout << incdSnpInfoVec[j]->ID << " ";
        }
        if(readBedBool) {
            cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
            readBedBool = false;
        }
        MatrixXf eigenVec;
        VectorXf eigenVal,cumsumNonNegPerLD;
        MatrixXf rval = generateLDmatrixPerBlock(bedFile, ldblock->gwasSnpNameVecInBlock);
        // cout << "rval: " << rval << endl;
        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        eigenDecomposition(rval, eigenCutoff,eigenVal, eigenVec,cumsumNonNegPerLD);
        // cout << "rval: " << rval.row(0) << endl;
        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
        // save svd matrix
        int32_t numEigenValue = eigenVal.size();
        int32_t numSnpInBlock = ldblock->gwasSnpNameVecInBlock.size();
        eigenNumInldblock[i] = numEigenValue;
        float eigenValueSum = eigenVal.sum();
        // save summary
        // 1, nrow of eigenVecGene[i]
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, out3);
        //        cout << "rval: " << rval << endl;
        // cout << "eigenVec: "  << eigenVec << endl;
        // cout << "";
        // 2, ncol of eigenVecGene[i]
        fwrite(&numEigenValue, sizeof(int32_t), 1, out3);
        // 3. sum of eigen values
        fwrite(&eigenValueSum, sizeof(float), 1, out3);
        //4. eigenCutoff
        fwrite(&eigenCutoff, sizeof(float), 1, out3);
        // 5. save eigen value
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, out3);
        // 6. save a series of proportions of variance in LD, please note, the order of
        // these proportions are the opposite of the order of eigenValues
        fwrite(cumsumNonNegPerLD.data(), sizeof(float), numSnpInBlock, out3);
        // 7. save eigen vector;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, out3);
        cout << " Generate and save Eigen decomposition result for LD block " << i << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << "\r" << flush;
    }
    fclose(out3);
    // cout << "size of eigenValuegene: " << eigenValGene.size() << endl;
    //outputEigenDataForLDM(filename);
    
    cout << "To explain " << eigenCutoff*100 << "% variance in LD, on average " << int(eigenNumInldblock.mean()) << " eigenvalues are selected across LD blocks (mean number of SNPs is " << int(snpNumInldblock.mean()) << ")." << endl;

}

void Data::outputLdDataForBlockLDM(const string &filename) const {
    
    string outfilename;
        outfilename = filename + ".ldm.full";
    
    string outfile1 = outfilename + ".info";
    string outfile2 = outfilename + ".blk.info";

    // write snp info
    ofstream out1(outfile1.c_str());
    out1 << boost::format("%6s %15s %10s %15s %6s %6s %12s %10s\n")
    % "Chrom"
    % "ID"
    % "GenPos"
    % "PhysPos"
    % "A1"
    % "A2"
    % "A1Freq"
    % "N";
    SnpInfo *snp, *windStart, *windEnd;
    for (unsigned i=0; i < numIncdSnps ; ++i) {
        snp = incdSnpInfoVec[i];
        out1 << boost::format("%6s %15s %10s %15s %6s %6s %12f %10s \n")
        % snp->chrom
        % snp->ID
        % snp->genPos
        % snp->physPos
        % snp->a1
        % snp->a2
        % snp->af
        % numKeptInds;
    }
    out1.close();

    // svd matrix for ld blocks here.
    ofstream out2(outfile2.c_str());
    out2 << boost::format("%6s %15s %10s %15s %15s %15s\n")
    % "Chrom"
    % "ldBlockID"
    % "start"
    % "end"
    % "snpInLdBlock"
    % "NumSnpInLdBlock";
    LDBlockInfo * ldblock;
    for (unsigned i=0; i < numKeptLDBlocks ; ++i) {
        ldblock = keptLdBlockInfoVec[i];
        cout << "ldblock id: " << ldblock->ID << endl;
        for(unsigned j = 0; j < ldblock->gwasSnpNameVecInBlock.size() ; ++j){
            out2 << boost::format("%6s %15s %10s %15s %15s %15s \n")
                % ldblock->chrom
                % ldblock->ID
                % ldblock->startPos
                % ldblock->endPos
                % ldblock->gwasSnpNameVecInBlock[j]
                % ldblock->gwasSnpNameVecInBlock.size();
        }
    }
    out2.close();
    
    cout << "Written the LD matrix SNP info into file [" << outfile1 << "]." << endl;
    cout << "Written the LD matrix block info into file [" << outfile2 << "]." << endl;
    //cout << "Written SVD matrices from per LD BLOCK matrix into file [" << outfile3 << "](Set proportion of variance in LD as " +to_string(eigenCutoff)+ ")." << endl;
}


///////////// Step 2.1 read LD matrix of ld blocks info (5)
void Data::readBlockLDMblockInfoFile(const string &infoFile){
    // Read bim file: recombination rate is defined between SNP i and SNP i-1
    ifstream in(infoFile.c_str());
    if (!in) throw ("Error: can not open the file [" + infoFile + "] to read.");
    cout << "Reading ld ldm info from file [" + infoFile + "]." << endl;
    ldBlockInfoVec.clear();
    ldBlockInfoMap.clear();
    map<string,int> ld2snpMap;
        
    string header;
    string id;
    int  chr, blockStart, blockEnd, snpNum;
    int idx = 0;
    int snpCount =  1;
    string snpName;
    LDBlockInfo *ldblock;
    getline(in, header);
    while (in >>chr>>id >>blockStart>>blockEnd>>snpName >> snpNum ) {
        if (ld2snpMap.insert(pair<string, int>(id + "_" + snpName , snpCount)).second == false) {
            throw ("Error: Duplicate LDBlock-SNP pair found: \"" + id + "_" + snpName + "\".");
        } else{
            if(snpCount == 1){
                ldblock = new LDBlockInfo(idx++, id, chr);
                ldblock->startPos = blockStart;
                ldblock->endPos   = blockEnd;
            }
            if(snpCount < snpNum){
                ldblock->gwasSnpNameVecInBlock.push_back(snpName);
            }
            if(snpCount == snpNum){
                ldblock->gwasSnpNameVecInBlock.push_back(snpName);
                ldblock->numSnpInBlock = snpNum;
                ldBlockInfoVec.push_back(ldblock);
                if (ldBlockInfoMap.insert(pair<string, LDBlockInfo*>(id, ldblock)).second == false) {
                        throw ("Error: Duplicate LD block ID found: \"" + id + "\".");}
                snpCount = 1;
                continue;
            }
            snpCount++;
        }
    }
    in.close();
    numLDBlocks = (unsigned) ldBlockInfoVec.size();
    cout << numLDBlocks << " LD Blocks to be included from [" + infoFile + "]." << endl;
}

void Data::readBlockLDMsnpInfoFile(const string &snpInfoFile){
    ifstream in(snpInfoFile.c_str());
    if (!in) throw ("Error: can not open the file [" + snpInfoFile + "] to read.");
    cout << "Reading ld snp info from file [" + snpInfoFile + "]." << endl;
    snpInfoVec.clear();
    snpInfoMap.clear();
    string header;
    string id, allele1, allele2;
    int chr, physPos,ld_n;
    float genPos;
    float allele1Freq;
    int idx = 0;
    getline(in, header);
    while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2>>allele1Freq>>ld_n) {
        SnpInfo *snp = new SnpInfo(idx++, id, allele1, allele2, chr, genPos, physPos);
        snp->af = allele1Freq;
        snp->ld_n = ld_n;
        snpInfoVec.push_back(snp);
        chromosomes.insert(snp->chrom);
        if (snpInfoMap.insert(pair<string, SnpInfo*>(id, snp)).second == false) {
            throw ("Error: Duplicate SNP ID found: \"" + id + "\".");
        }
    }
    in.close();
    numSnps = (unsigned) snpInfoVec.size();
    cout << numSnps << " SNPs to be included from [" + snpInfoFile + "]." << endl;
}

void Data::readBlockLDMbinaryAndDoEigenDecomposition(const string &binFile, const string &filename, const float &eigenCutoff, const bool writeLdmTxt){
    FILE *fp = fopen(binFile.c_str(), "rb");
    if(!fp){throw ("Error: can not open the file [" + binFile + "] to read.");}
    vector<int> numSnpInRegion;
    LDBlockInfo * ldblock;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        ldblock = ldBlockInfoVec[i];
        numSnpInRegion[i] = ldblock->numSnpInBlock;
    }
        
    string outfilename = filename + ".eigen.bin";
    FILE *out = fopen(outfilename.c_str(), "wb");
    
    ofstream out2;
    string outfile2;
    if (writeLdmTxt) {
        outfile2 = filename + ".eigen.txt";
        out2.open(outfile2.c_str());
    }

    keptLdBlockInfoVec = ldBlockInfoVec;

    for(int i = 0; i < numLDBlocks; i++){
        
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        int32_t blockSize = block->numSnpInBlock;
        
        MatrixXf ldm(blockSize, blockSize);
        uint64_t nElements = (uint64_t)blockSize * (uint64_t)blockSize;
                
        if(fread(ldm.data(), sizeof(float), nElements, fp) != nElements){
            cout << "fread(U.data(), sizeof(float), nElements, fp): " << fread(ldm.data(), sizeof(float), nElements, fp) << endl;
            cout << "nEle: " << nElements << " ldm.size: " << ldm.size() <<  " ldm.col: " << ldm.cols() << " row: " << ldm.rows() << endl;
            throw("In LD block " + to_string(i) + ",size error in " + binFile);
            // cout << "Read " << svdLDfile << " error (U)" << endl;
            // throw("read file error");
        }
        
        MatrixXf eigenVec;
        VectorXf eigenVal, cumsumNonNegPerLD;
        // cout << "rval: " << rval << endl;
        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        eigenDecomposition(ldm, eigenCutoff, eigenVal, eigenVec, cumsumNonNegPerLD);
        // cout << "rval: " << rval.row(0) << endl;
        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
        // save svd matrix

        
        int32_t numEigenValue = eigenVal.size();
        int32_t numSnpInBlock = blockSize;
        float eigenValueSum = eigenVal.sum();

        cout << " Generate Eigen decomposition result for LD block " << i << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << endl;
        
        // save summary
        // 1, the number of SNPs in the block
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, out);
        // 2, the number of eigenvalues at with the given cutoff
        fwrite(&numEigenValue, sizeof(int32_t), 1, out);
        // 3. sum of the selected eigenvalues
        fwrite(&eigenValueSum, sizeof(float), 1, out);
        //4. eigenvalue cutoff based on the proportion of variance explained in LD
        fwrite(&eigenCutoff, sizeof(float), 1, out);
        // 5. save the selected eigenvalues
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, out);
        // 6. save a series of proportions of variance in LD, please note, the order of
        // these proportions are the opposite of the order of eigenValues
        fwrite(cumsumNonNegPerLD.data(), sizeof(float), numSnpInBlock, out);
        // 7. save eigen vector;
        nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, out);
        
        if (writeLdmTxt) {
            out2 << "Block " << block->ID << endl;
            out2 << "numSnpInBlock " << numSnpInBlock << endl;
            out2 << "numEigenValue " << numEigenValue << endl;
            out2 << "eigenValueSum " << eigenValueSum << endl;
            out2 << "eigenCutoff " << eigenCutoff << endl;
            out2 << "eigenVal\n" << eigenVal.transpose() << endl;
            out2 << "cumsumNonNegPerLD\n" << cumsumNonNegPerLD.transpose() << endl;
            out2 << "eigenVec\n" << eigenVec << endl;
            out2 << endl;
        }
        
    }

    fclose(out);

    if (writeLdmTxt) {
        out2.close();
        cout << "Written the eigen data for block LD matrix into text file [" << outfile2 << "]." << endl;
    }

}

void Data::readEigenMatrixBinaryFile(const string &eigenBinFile, const float eigenCutoff){
    FILE *fp = fopen(eigenBinFile.c_str(), "rb");
    if(!fp){throw ("Error: can not open the file [" + eigenBinFile + "] to read.");}
    vector<int>numSnpInRegion;
    LDBlockInfo * ldblock;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        ldblock = ldBlockInfoVec[i];
        numSnpInRegion[i] = ldblock->numSnpInBlock;
    }
    eigenValLdBlock.resize(numLDBlocks);
    eigenVecLdBlock.resize(numLDBlocks);
    
    for(int i = 0; i < numLDBlocks; i++){
        int32_t cur_m = 0;
        int32_t cur_k = 0;
        float sumLambda = 0;
        float svdVarProp =0;
        // 1. marker number
        if(fread(&cur_m, sizeof(int32_t), 1, fp) != 1){
            throw("Read " + eigenBinFile + " error (m)");
        }
        if(cur_m != numSnpInRegion[i]){
            throw("In LD block " + to_string(i) + ", inconsistent marker number to marker information in " + eigenBinFile);
        }
        // 2. ncol of eigenVec (number of eigenvalues)
        if(fread(&cur_k, sizeof(int32_t), 1, fp) != 1){
            throw("In LD block " + to_string(i) + ", error about number of eigenvalues in  " + eigenBinFile);
            // cout << "Read " << eigenBinFile << " error (k)" << endl;
            // throw("read file error");
        }
        // 3. sum of eigen values
        if(fread(&sumLambda, sizeof(float), 1, fp) != 1){
            throw("In LD block " + to_string(i) + ", error about sumLambda in " + eigenBinFile);
            // cout << "Read " << eigenBinFile << " error sumLambda" << endl;
            // throw("read file error");
        }
        // 4. eigenCutoff
        if(fread(&svdVarProp, sizeof(float), 1, fp) != 1){
            throw("In LD block " + to_string(i) + ", error about svdVarProp used in " + eigenBinFile);
            // cout << "Read " << eigenBinFile << " error svdVarProp" << endl;
            // throw("read file error");
        }
        // 5. eigen values
        VectorXf lambda(cur_k);
        if(fread(lambda.data(), sizeof(float), cur_k, fp) != cur_k){
            throw("In LD block " + to_string(i) + ",size error about eigenvalues in " + eigenBinFile);
            // cout << "Read " << eigenBinFile << " error (lambda)" << endl;
            // throw("read file error");
        }
        // 6. save a series of proportions of variance in LD
        VectorXf cumsumNonNeg(cur_m);
        if(fread(cumsumNonNeg.data(), sizeof(float), cur_m, fp) != cur_m){
            throw("In LD block " + to_string(i) + ",size error cumsum of non-negative eigenvalues in " + eigenBinFile);
            // cout << "Read " << eigenBinFile << " error (cumsumNonNeg)" << endl;
            // throw("read file error");
        }
        //7. read eigen vector
        MatrixXf U(cur_m, cur_k);
        uint64_t nElements = (uint64_t)cur_m * (uint64_t)cur_k;
        if(fread(U.data(), sizeof(float), nElements, fp) != nElements){
            cout << "fread(U.data(), sizeof(float), nElements, fp): " << fread(U.data(), sizeof(float), nElements, fp) << endl;
            cout << "nEle: " << nElements << " U.size: " << U.size() <<  " U.col: " << U.cols() << " row: " << U.rows() << endl;
            throw("In LD block " + to_string(i) + ",size error about eigenvectors in " + eigenBinFile);
            // cout << "Read " << eigenBinFile << " error (U)" << endl;
            // throw("read file error");
        }
        bool haveValue = false;
        int revIdx = 0;
        if(svdVarProp != eigenCutoff & i == 0){
            cout << "Warning: current proportion of variance in LD block is set as " + to_string(eigenCutoff)+ ". But the proportion of variance is set as "<< to_string(svdVarProp) + " in "  + eigenBinFile + ".\n";
            // throw("");
        }
        // cout << "lambda: " << lambda << endl;
        // cout << "U: " << U << endl;
        // eigenVecLdBlock[i] = U;
        // eigenValLdBlock[i] = lambda;
        for (;revIdx < cur_m; revIdx ++ ){
            if(cumsumNonNeg(revIdx) < 1e-10 || cumsumNonNeg(revIdx) == 0 ) {continue;}
            if(eigenCutoff >= cumsumNonNeg(revIdx) ){
                revIdx = revIdx -1;
                haveValue = true;
                break;
            }
        }
        if(!haveValue) revIdx = cur_m - 1;
        eigenVecLdBlock[i] = U.rightCols(cur_m - revIdx);
        eigenValLdBlock[i] = lambda.tail(cur_m - revIdx);
        // cout << "U: " <<  U << endl;
    }
}

void Data::readBlockLDmatrixAndDoEigenDecomposition(const string &LDmatrixFile, const float eigenCutoff, const bool writeLdmTxt){
    readBlockLDMblockInfoFile(LDmatrixFile + ".blk.info");
    readBlockLDMsnpInfoFile(LDmatrixFile + ".info");
    readBlockLDMbinaryAndDoEigenDecomposition(LDmatrixFile + ".bin", LDmatrixFile, eigenCutoff, writeLdmTxt);
}

void Data::readEigenMatrix(const string &eigenMatrixFile, const float eigenCutoff){
    readBlockLDMblockInfoFile(eigenMatrixFile + ".blk.info");
    readBlockLDMsnpInfoFile(eigenMatrixFile + ".info");
    readEigenMatrixBinaryFile(eigenMatrixFile + ".eigen.bin", eigenCutoff);
}

vector<LDBlockInfo*> Data::makeKeptLDBlockInfoVec(const vector<LDBlockInfo*> &ldBlockInfoVec){
    vector<LDBlockInfo*> keptLDBlock;
    ldblockNames.clear();
    LDBlockInfo * ldblock = NULL;
    for (unsigned i=0, j=0; i< numLDBlocks; ++i) {
        ldblock = ldBlockInfoVec[i];
        if(ldblock->kept) {
            ldblock->index = j++;  // reindex inds
            keptLDBlock.push_back(ldblock);
            ldblockNames.push_back(ldblock->ID);
        }
    }
    return keptLDBlock;
}

void Data::buildMMEeigen(const bool sampleOverlap, const bool noscale){
    includeMatchedBlocks();
    // update gwas information
    if (numIncdSnps!=0) constructWandQ(noscale);
    if (numAnnos) setAnnoInfoVec();
    lowRankModel = true;
}

void Data::includeMatchedBlocks(){
    // this step is to construct gwasSnp2geneVec
    //cout << "Construct various maps." << endl;
    SnpInfo * snp;
    LDBlockInfo * ldblock;
    
    for(unsigned j = 0; j < numLDBlocks; j++){
        ldblock = ldBlockInfoVec[j];
        ldblock->block2GwasSnpVec.clear();
        for(int i = 0; i < numIncdSnps; i++){
            snp = incdSnpInfoVec[i];
            // find if given snp belongs to ld block;
            if (std::find(ldblock->gwasSnpNameVecInBlock.begin(), ldblock->gwasSnpNameVecInBlock.end(), snp->ID) != ldblock->gwasSnpNameVecInBlock.end()) {
                ldblock->block2GwasSnpVec.push_back(i);
                ldblock->memberSnpVec.push_back(snp);
            }
        }
        if(ldblock->block2GwasSnpVec.size() == 0){
            ldblock->kept = false;
        } else {
            ldblock->startSnpIdx = ldblock->block2GwasSnpVec[0];
            ldblock->endSnpIdx = ldblock->block2GwasSnpVec[ldblock->numSnpInBlock-1];
        }
    }
    
    keptLdBlockInfoVec = makeKeptLDBlockInfoVec(ldBlockInfoVec);
    numKeptLDBlocks = (unsigned) keptLdBlockInfoVec.size();

    ldblock2gwasSnpMap.clear();
    // Construct map from ld to snp
    for(unsigned i = 0; i < numKeptLDBlocks; i++){
        ldblock = keptLdBlockInfoVec[i];
        ldblock2gwasSnpMap.insert(pair<int, vector<int> > (i,ldblock->block2GwasSnpVec));
    }
    
    cout << numKeptLDBlocks << " LD blocks are included." << endl;
}

void Data::constructWandQ(const bool noscale){
    VectorXf nMinusOne;
    snp2pq.resize(numIncdSnps);
    D.resize(numIncdSnps);
   // ZPZdiag.resize(numIncdSnps);
    ZPy.resize(numIncdSnps);
    b.resize(numIncdSnps);
    n.resize(numIncdSnps);
    nMinusOne.resize(numIncdSnps);
    se.resize(numIncdSnps);
    tss.resize(numIncdSnps);
    SnpInfo *snp;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snp = incdSnpInfoVec[i];
        snp->af = snp->gwas_af;
        snp2pq[i] = snp->twopq = 2.0f*snp->gwas_af*(1.0f-snp->gwas_af);
        if(snp2pq[i]==0) cout << "Error: SNP " << snp->ID << " af " << snp->af << " has 2pq = 0." << endl;
        D[i] = snp2pq[i]*snp->gwas_n;
        b[i] = snp->gwas_b * sqrt(snp2pq[i]); // scale the marginal effect so that it's in per genotype SD unit
        n[i] = snp->gwas_n;
        nMinusOne[i] = snp->gwas_n - 1;
        se[i]= snp->gwas_se;
        tss[i] = D[i]*(n[i]*se[i]*se[i] + b[i]*b[i]);
        //  D[i] = 1.0/(se[i]*se[i]+b[i]*b[i]/snp->gwas_n);  // NEW!
        //  snp2pq[i] = snp->twopq = D[i]/snp->gwas_n;       // NEW!
    }
    cout << endl;
    LDBlockInfo * ldblock;
    wcorrBlocks.resize(numKeptLDBlocks);
    // Qblocks.resize(numKeptLDBlocks);
    numSnpsBlock.resize(numKeptLDBlocks);
    numEigenvalBlock.resize(numKeptLDBlocks);
    Qblocks.clear();
    VectorXf sqrtLambda;
    // save gwas marginal effect into block
    gwasMarginEffectInBlock.resize(numKeptLDBlocks);
    for (unsigned i = 0; i < numKeptLDBlocks; i++){
        ldblock = keptLdBlockInfoVec[i];
        gwasMarginEffectInBlock[i] = b(ldblock->block2GwasSnpVec);
        // calculate wbcorr and Qblocks
        sqrtLambda = eigenValLdBlock[i].array().sqrt();
        //cout << eigenVecLdBlock[i].transpose().rows() << " " << eigenVecLdBlock[i].transpose().cols() << " " << gwasMarginEffectInBlock[i].size() << endl;
        wcorrBlocks[i] = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (eigenVecLdBlock[i].transpose() * gwasMarginEffectInBlock[i] );
        // cout << "eigenVecLdBlock[i]: " << eigenVecLdBlock[i] << endl;
        // cout << "wcorrBlocks[i]: " << wcorrBlocks[i] << endl;
        // cout << "sqrtLambda: " << sqrtLambda << endl;
        // cout << gwasMarginEffectInBlock[i] << endl;
        MatrixXf tmpQblocks = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
        MatrixDat matrixDat = MatrixDat(ldblock->gwasSnpNameVecInBlock,tmpQblocks );
        // cout << "Qblock: " << endl;
        // cout << matrixDat.values << endl;
        Qblocks.push_back(matrixDat);
        numSnpsBlock[i] = Qblocks[i].ncol;
        numEigenvalBlock[i] = Qblocks[i].nrow;
    }
    
    //b.array() -= b.mean();  // DO NOT CENTER b
    // estimate phenotypic variance based on the input allele frequencies in GWAS

    //  Vp_buf = h_buf * N_buf * se_buf * se_buf + h_buf * b_buf * b_buf * N_buf / (N_buf - 1.0);
    //VectorXf ypySrt = D.array()*n.array()*se.array().square() + D.array() * b.array().square() * n.array() / nMinusOne.array();
    VectorXf ypySrt = D.array()*(n.array()*se.array().square()+b.array().square());
    VectorXf varpSrt = ypySrt.array()/n.array();
    std::sort(ypySrt.data(), ypySrt.data() + ypySrt.size());
    std::sort(varpSrt.data(), varpSrt.data() + varpSrt.size());
    ypy = ypySrt[ypySrt.size()/2];  // median
    varPhenotypic = varpSrt[varpSrt.size()/2];
    //cout << "varPhenotypic: " << varPhenotypic << endl;
    VectorXf nSrt = n;
    std::sort(nSrt.data(), nSrt.data() + nSrt.size());
    numKeptInds = nSrt[nSrt.size()/2]; // median
    
    nGWASblock.resize(numKeptLDBlocks);
    for (unsigned i=0; i<numKeptLDBlocks; ++i) {
        nGWASblock[i] = numKeptInds;
    }

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
    // cout << endl << "snp2pq: " << endl << snp2pq << endl;
    //ypy = numKeptInds;
    // NEW END
    // data summary
    cout << "\nData summary:" << endl;
    cout << boost::format("%40s %8s %8s\n") %"" %"mean" %"sd";
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP Phenotypic variance" %Gadget::calcMean(varpSrt) %sqrt(Gadget::calcVariance(varpSrt));
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP heterozygosity" %Gadget::calcMean(snp2pq) %sqrt(Gadget::calcVariance(snp2pq));
    cout << boost::format("%40s %8.0f %8.0f\n") %"GWAS SNP sample size" %Gadget::calcMean(n) %sqrt(Gadget::calcVariance(n));
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP effect (in genotype SD unit)" %Gadget::calcMean(b) %sqrt(Gadget::calcVariance(b));
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP SE" %Gadget::calcMean(se) %sqrt(Gadget::calcVariance(se));
    cout << boost::format("%40s %8.3f %8.3f\n") %"LD block size" %Gadget::calcMean(numSnpsBlock) %sqrt(Gadget::calcVariance(numSnpsBlock));
    cout << boost::format("%40s %8.3f %8.3f\n") %"LD block rank" %Gadget::calcMean(numEigenvalBlock) %sqrt(Gadget::calcVariance(numEigenvalBlock));
}
