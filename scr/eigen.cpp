//
//  eigen.cpp
//  gctb
//
//  Created by Shouye Liu and Jian Zeng on 05/01/2023.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "data.hpp"
#include <sys/stat.h>
#include <dirent.h>
#include <cmath>
#include <stdexcept>

namespace {
/** Anchor SNPs for --impute-summary: must have finite beta and positive finite SE (z = b/SE). */
inline bool gwasUsableAnchorSnp(const SnpInfo *snp) {
    return snp->included && std::isfinite(snp->gwas_b) && std::isfinite(snp->gwas_se) && snp->gwas_se > 0.0;
}

[[noreturn]] inline void throwEigenReadErr(const string &msg) {
    cerr << msg << endl;
    throw std::runtime_error(msg);
}
}

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
    int  chr,start, stop;
    int idx = 0;
    getline(in, header);
    while (in >> id >> chr >> start >> stop) {
        LDBlockInfo *ld = new LDBlockInfo(idx++, id, chr);
        ld->startPos    = start;
        ld->endPos     = stop;
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

void Data::eigenDecomposition( const MatrixXf &X, const float &prop, VectorXf &eigenValAdjusted, MatrixXf &eigenVecAdjusted, float &sumPosEigVal){
    // VectorXf cumsumNonNeg; // cumulative sums of non-negative values
    sumPosEigVal = 0.0;

    SelfAdjointEigenSolver<MatrixXf> eigensolver(X);
    VectorXf eigenVal = eigensolver.eigenvalues();
    MatrixXf eigenVec = eigensolver.eigenvectors();
    int n = (int)eigenVal.size();
    if (n == 0) {
        eigenValAdjusted.resize(0);
        eigenVecAdjusted.resize(0, 0);
        return;
    }
    // 1x1 matrix (e.g. sub-LD with one retained SNP): one eigenvalue explains all variance; the general loop uses revIdx-1 and breaks for n==1.
    if (n == 1) {
        if (eigenVal[0] < 0) cout << "Error, all eigenvector are negative" << endl;
        sumPosEigVal = (eigenVal[0] > 1e-10f) ? eigenVal[0] : 0.f;
        eigenValAdjusted.resize(1);
        eigenValAdjusted[0] = eigenVal[0];
        eigenVecAdjusted = eigenVec;
        return;
    }

    int revIdx = n;
    VectorXf cumsumNonNeg(revIdx);
    cumsumNonNeg.setZero();
    revIdx = revIdx -1;
    if(eigenVal[revIdx] < 0) cout << "Error, all eigenvector are negative" << endl;
    cumsumNonNeg[revIdx] = eigenVal[revIdx];
    sumPosEigVal = eigenVal[revIdx];
    revIdx = revIdx -1;
    
    int numPosEigVal = 0;
    while (revIdx >= 0 && eigenVal[revIdx] > 1e-10 ){
        sumPosEigVal = sumPosEigVal + eigenVal[revIdx];
        cumsumNonNeg[revIdx] = eigenVal[revIdx] + cumsumNonNeg[revIdx + 1];
        ++numPosEigVal;
        revIdx = revIdx - 1;
        if(revIdx < 0) break;
    }
    // cout << "revIdx: " << revIdx << endl;
    // cout << "size: " << eigenVal.size()  << " eigenVal: " << eigenVal << endl;
    // cout << "cumsumNoNeg: " << cumsumNonNeg << endl;
    cumsumNonNeg = cumsumNonNeg/sumPosEigVal;
    bool haveValue = false;
    // cout << "cumsumNonNeg: " << cumsumNonNeg << endl;
    // cout << "revIdx : " << revIdx << endl;
    // cout << "revIdx: "  << revIdx  << endl;
    for (revIdx = revIdx + 1; revIdx < eigenVal.size(); revIdx ++ ){
        // cout << "revIdx: " << revIdx << " cumsumNonNeg: " << cumsumNonNeg(revIdx) << endl;
        if(prop >= cumsumNonNeg[revIdx] ){
            revIdx = revIdx -1;
            haveValue = true;
            break;
        }
    }
    // cout << "revIdx : " << revIdx << endl;
    //if(!haveValue) revIdx = eigenVal.size() - 1;
    if(!haveValue) revIdx = numPosEigVal - 1;
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
        
        if (snp2pq[incj]==0) throw ("Error: " + snpj->ID + " is a fixed SNP (MAF=0)!");
        
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
            
            if (snp2pq[inck]==0) throw ("Error: " + snpk->ID + " is a fixed SNP (MAF=0)!");
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
            
            if (snp2pq[inck]==0) throw ("Error: " + snpk->ID + " is a fixed SNP (MAF=0)!");

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
        VectorXf eigenVal;
        float sumPosEigVal;
        // cout << "rval: " << rval << endl;
        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        eigenDecomposition(ZPZmat, eigenCutoff, eigenVal, eigenVec, sumPosEigVal);
        // cout << "rval: " << rval.row(0) << endl;
        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
        // save svd matrix
        
        int32_t numEigenValue = eigenVal.size();
        int32_t numSnpInBlock = keptLdBlockInfoVec[blk]->numSnpInBlock;  // TMP

        cout << " Generate Eigen decomposition result for LD block " << blk << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << endl;
        
        // save summary
        // 1, the number of SNPs in the block
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, out3);
        // 2, the number of eigenvalues at with the given cutoff
        fwrite(&numEigenValue, sizeof(int32_t), 1, out3);
        // 3. sum of all the positive eigenvalues
        fwrite(&sumPosEigVal, sizeof(float), 1, out3);
        // 4. eigenvalue cutoff based on the proportion of variance explained in LD
        fwrite(&eigenCutoff, sizeof(float), 1, out3);
        // 5. the selected eigenvalues
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, out3);
        // 6. the selected eigenvector;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, out3);
        
    }
    
    fclose(out3);

    //outputEigenDataForLDM(filename);
    
}

void Data::mapSnpsToBlocks(){
    //unsigned snpIdx = 0;
    int mapped = 0;
    ChromInfo* chr;
    LDBlockInfo *block;
    SnpInfo *snp;
    vector<ChromInfo*>::iterator it = chromInfoVec.begin(), end = chromInfoVec.end();
    unsigned chrCur = (*it)->id;
//    int chrStartSnpIdx = (*it)->startSnpIdx;
//    int chrEndSnpIdx = (*it)->endSnpIdx;
    
    for (unsigned i=0; i<numLDBlocks; ++i) {
        block = ldBlockInfoVec[i];
        block->snpNameVec.clear();
        block->snpInfoVec.clear();
        
//        cout << i << " " << numLDBlocks << " " << block->chrom << " " << chrCur << endl;
        
        if (block->chrom > chrCur) {  // move to a new chromosome
            if (it != end) ++it;
            if (it == end) {
                block->kept = false;
                continue;
            }
//            chrStartSnpIdx = (*it)->startSnpIdx;
//            chrEndSnpIdx = (*it)->endSnpIdx;
            chrCur = (*it)->id;
        }
        if (block->chrom != chrCur) {
            block->kept = false;
            continue;
        }
        
//        cout << i << " " << numLDBlocks << " " << block->chrom << " " << block->startPos << " " << block->endPos << endl;

        for (unsigned snpIdx = 0; snpIdx < numSnps; ++snpIdx) {
            //cout << "snpIdx " << snpIdx << " chrEndSnpIdx " << chrEndSnpIdx << endl;
            snp = snpInfoVec[snpIdx];
            if (!snp->included) continue;
            if (snp->chrom != chrCur) continue;
            if (snp->physPos < block->startPos) continue;
            else if (snp->physPos > block->endPos) break;
            block->snpNameVec.push_back(snp->ID);
            block->snpInfoVec.push_back(snp);
            snp->block = block->ID;
            snp->blockIdx = block->index;
//            cout << snpIdx << " " << snp->ID << " " << snp->chrom << " " << snp->physPos << " " << block->startPos << " " << block->endPos << " " << snp->included << endl;
        }
        block->numSnpInBlock = block->snpInfoVec.size();
        if (block->numSnpInBlock) {
            block->startSnpIdx = block->snpInfoVec[0]->index;
            block->endSnpIdx = block->snpInfoVec[block->numSnpInBlock-1]->index;
            block->kept = true;
            ++mapped;
        } else {
            block->kept = false;
        }
        //cout << i << " " << numLDBlocks << " " << block->numSnpInBlock << endl;

    }
        
    if (mapped < 1) throw(0, "No SNP can be mapped to the provided LD block list. Please check the input data regarding chromosome and bp.");
    else cout << mapped << " LD block(s) are retained." << endl;
}

void Data::makeBlockLDmatrix(const string &bedFile, const string &LDmatType, const unsigned block, const string &dirname, const bool writeLdmTxt, int ldBlockRegionWind){
    cout << "Making block LD matricies ..." << endl;
    
    struct stat sb;
    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        string create_cmd = "mkdir " + dirname;
        system(create_cmd.c_str());
        cout << "Created folder [" << dirname << "] to store LD matrices." << endl;
    }
    
    mapSnpsToBlocks();
    
    keptLdBlockInfoVec = makeKeptLDBlockInfoVec(ldBlockInfoVec);
    numKeptLDBlocks = (unsigned) keptLdBlockInfoVec.size();
    
#pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < numKeptLDBlocks; i++) {
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];

        string outBinfile = dirname + "/block" + ldblock->ID + ".ldm.bin";
        FILE *outbin = fopen(outBinfile.c_str(), "wb");
        ofstream outtxt;
        string outTxtfile;
        if (writeLdmTxt) {
            outTxtfile = dirname + "/block" + ldblock->ID + ".ldm.txt";
            outtxt.open(outTxtfile.c_str());
        }
        string outSnpfile = dirname + "/block" + ldblock->ID + ".snp.info";
        string outldmfile = dirname + "/block" + ldblock->ID + ".ldm.info";
        string pwldfilename = dirname + "/block" + ldblock->ID + ".rsq0.5.pwld";  // pairwise LD file

        if(!i) cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
            
        MatrixXf rval = generateLDmatrixPerBlock(bedFile, ldblock->snpNameVec);
        
//        unsigned numSnpInBlock = ldblock->numSnpInBlock;
//        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numSnpInBlock;
//        fwrite(rval.data(), sizeof(float), nElements, outbin);
        
        // save the lower trangular matrix only
        unsigned numSnpInBlock = ldblock->numSnpInBlock;
        for (unsigned row = 0; row < numSnpInBlock; ++row) {
            for (unsigned col = 0; col <= row; ++col) {
                float value = rval(row, col);
                fwrite(&value, sizeof(float), 1, outbin);
            }
        }
        
        if (writeLdmTxt) {
//            for (unsigned ii=0; ii<numSnpInBlock; ++ii){
//                for (unsigned jj=0; jj<numSnpInBlock; ++jj) {
//                    outtxt << ldblock->ID << "\t" << ldblock->snpNameVec[ii] << "\t" << ldblock->snpNameVec[jj] << "\t" << rval(ii,jj) << endl;
//                }
//            }
            for (unsigned row = 0; row < numSnpInBlock; ++row) {
                for (unsigned col = 0; col <= row; ++col) {
                    outtxt << rval(row, col) << "\t";
                }
                outtxt << endl;
            }
        }
        
        fclose(outbin);
        if (writeLdmTxt) outtxt.close();
        
        outputBlockLDmatrixInfo(*ldblock, outSnpfile, outldmfile);
        
        outputLDfriends(rval, ldblock, dirname);
        
        
        if(!(i%1)) cout << " computed block " << ldblock->ID << "\r" << flush;

        if (block) {
            cout << "Written the LD matrix into file [" << outBinfile << "]." << endl;
            if (writeLdmTxt) cout << "Written the LD matrix into file [" << outTxtfile << "]." << endl;
            cout << "Written the LD matrix SNP info into file [" << outSnpfile << "]." << endl;
            cout << "Written the LD matrix ldm info into file [" << outldmfile << "]." << endl;
            cout << "Written the high pairwise LD correlations (rsq>0.5) into [" << pwldfilename << "]." << endl;
        }
    }
    
    if (!block) {
        cout << "Written the LD matrix into folder [" << dirname << "/block*.ldm.bin]." << endl;
        if (writeLdmTxt) cout << "Written the LD matrix into text file [" << dirname << "/block*.ldm.txt]." << endl;
        
        //if (chromInfoVec.size() >= 22) {  // genome-wide build of LD matrices
        
            mergeLdmInfo(LDmatType, dirname, true);
        
        //} else {
//            cout << "Written the LD matrix into folder [" << dirname << "/block*.snp.info]." << endl;
//            cout << "Written the LD matrix into folder [" << dirname << "/block*.ldm.info]." << endl;
//            cout << "Written the high pairwise LD correlations (rsq>0.5) into [" << dirname << "/block*.rsq0.5.pwld]." << endl;
        //}
    }
}

void Data::makeBlockLDmatrixFromEigen(const string &eigenDirname, const float eigenCutoff, const string &LDmatType, const unsigned block, const string &dirname, const bool writeLdmTxt, int ldBlockRegionWind){
    (void)ldBlockRegionWind;
    cout << "Making block LD matrices from eigen decomposition in [" << eigenDirname << "] ..." << endl;
    
    struct stat sb;
    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        string create_cmd = "mkdir " + dirname;
        system(create_cmd.c_str());
        cout << "Created folder [" << dirname << "] to store LD matrices." << endl;
    }
    
    readBlockLdmInfoFile(eigenDirname, block);
    readBlockLdmSnpInfoFile(eigenDirname, block);
    readEigenMatrixBinaryFile(eigenDirname, eigenCutoff, false, ".");
    
#pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < numLDBlocks; i++) {
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (!ldblock->kept) continue;
        
        string outBinfile = dirname + "/block" + ldblock->ID + ".ldm.bin";
        FILE *outbin = fopen(outBinfile.c_str(), "wb");
        ofstream outtxt;
        string outTxtfile;
        if (writeLdmTxt) {
            outTxtfile = dirname + "/block" + ldblock->ID + ".ldm.txt";
            outtxt.open(outTxtfile.c_str());
        }
        string outSnpfile = dirname + "/block" + ldblock->ID + ".snp.info";
        string outldmfile = dirname + "/block" + ldblock->ID + ".ldm.info";
        string pwldfilename = dirname + "/block" + ldblock->ID + ".rsq0.5.pwld";  // pairwise LD file
        
        MatrixXf rval = eigenVecLdBlock[i] * eigenValLdBlock[i].asDiagonal() * eigenVecLdBlock[i].transpose();
        
        unsigned numSnpInBlock = ldblock->numSnpInBlock;
        for (unsigned row = 0; row < numSnpInBlock; ++row) {
            for (unsigned col = 0; col <= row; ++col) {
                float value = rval(row, col);
                fwrite(&value, sizeof(float), 1, outbin);
            }
        }
        
        if (writeLdmTxt) {
            for (unsigned row = 0; row < numSnpInBlock; ++row) {
                for (unsigned col = 0; col <= row; ++col) {
                    outtxt << rval(row, col) << "\t";
                }
                outtxt << endl;
            }
        }
        
        fclose(outbin);
        if (writeLdmTxt) outtxt.close();
        
        outputBlockLDmatrixInfo(*ldblock, outSnpfile, outldmfile);
        
        outputLDfriends(rval, ldblock, dirname);
        
        
        if(!(i%1)) cout << " computed block " << ldblock->ID << "\r" << flush;

        if (block) {
            cout << "Written the LD matrix into file [" << outBinfile << "]." << endl;
            if (writeLdmTxt) cout << "Written the LD matrix into file [" << outTxtfile << "]." << endl;
            cout << "Written the LD matrix SNP info into file [" << outSnpfile << "]." << endl;
            cout << "Written the LD matrix ldm info into file [" << outldmfile << "]." << endl;
            cout << "Written the high pairwise LD correlations (rsq>0.5) into [" << pwldfilename << "]." << endl;
        }
    }
    
    if (!block) {
        cout << "Written the LD matrix into folder [" << dirname << "/block*.ldm.bin]." << endl;
        if (writeLdmTxt) cout << "Written the LD matrix into text file [" << dirname << "/block*.ldm.txt]." << endl;
        
        mergeLdmInfo(LDmatType, dirname, true);
    }
}


void Data::outputBlockLDmatrixInfo(const LDBlockInfo &block, const string &outSnpfile, const string &outldmfile) const {
    // write snp info
    ofstream out1(outSnpfile.c_str());
    out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %12s %10s %10s\n")
    % "Chrom"
    % "ID"
    % "Index"
    % "GenPos"
    % "PhysPos"
    % "A1"
    % "A2"
    % "A1Freq"
    % "N"
    % "Block";
    SnpInfo *snp;
    for (unsigned i=0; i < block.numSnpInBlock; ++i) {
        snp = block.snpInfoVec[i];
        out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %12f %10s %10s\n")
        % snp->chrom
        % snp->ID
        % snp->index
        % snp->genPos
        % snp->physPos
        % snp->a1
        % snp->a2
        % snp->af
        % numKeptInds
        % snp->block;
    }
    out1.close();
    
    // svd matrix for ld blocks here.
    ofstream out2(outldmfile.c_str());
    out2 << boost::format("%10s %6s %15s %15s %15s %15s %12s\n")
    % "Block"
    % "Chrom"
    % "StartSnpIdx"
    % "StartSnpID"
    % "EndSnpIdx"
    % "EndSnpID"
    % "NumSnps";
    out2 << boost::format("%10s %6s %15s %15s %15s %15s %12s\n")
    % block.ID
    % block.chrom
    % block.startSnpIdx
    % block.snpInfoVec[0]->ID
    % block.endSnpIdx
    % block.snpInfoVec[block.numSnpInBlock-1]->ID
    % block.numSnpInBlock;
    out2.close();
}

//void Data::outputBlockLDmatrixInfo(const unsigned block, const string &dirname) const {
//    struct stat sb;
//    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
//        throw("Error: Folder " + dirname + " does not exist!");
//    }
//
//    string outSnpfile;
//    string outldmfile;
//
//    if (block) {
//        LDBlockInfo *blockinfo = keptLdBlockInfoVec[0];
//        outSnpfile = dirname + "/block" + blockinfo->ID + ".snp.info";
//        outldmfile = dirname + "/block" + blockinfo->ID + ".ldm.info";
//    } else {
//        outSnpfile = dirname + "/snp.info";
//        outldmfile = dirname + "/ldm.info";
//    }
//
//    // write snp info
//    ofstream out1(outSnpfile.c_str());
//    out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %12s %10s %10s\n")
//    % "Chrom"
//    % "ID"
//    % "Index"
//    % "GenPos"
//    % "PhysPos"
//    % "A1"
//    % "A2"
//    % "A1Freq"
//    % "N"
//    % "Block";
//    SnpInfo *snp;
//    for (unsigned i=0; i < numIncdSnps ; ++i) {
//        snp = incdSnpInfoVec[i];
//        out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %12f %10s %10s\n")
//        % snp->chrom
//        % snp->ID
//        % snp->index
//        % snp->genPos
//        % snp->physPos
//        % snp->a1
//        % snp->a2
//        % snp->af
//        % numKeptInds
//        % snp->block;
//    }
//    out1.close();
//
//    // svd matrix for ld blocks here.
//    ofstream out2(outldmfile.c_str());
//    out2 << boost::format("%10s %6s %15s %15s %15s %15s %12s\n")
//    % "Block"
//    % "Chrom"
//    % "StartSnpIdx"
//    % "StartSnpID"
//    % "EndSnpIdx"
//    % "EndSnpID"
//    % "NumSnps";
//    LDBlockInfo * ldblock;
//    for (unsigned i=0; i < numKeptLDBlocks ; ++i) {
//        ldblock = keptLdBlockInfoVec[i];
//        //cout << "ldblock id: " << ldblock->ID << endl;
//        out2 << boost::format("%10s %6s %15s %15s %15s %15s %12s\n")
//        % ldblock->ID
//        % ldblock->chrom
//        % ldblock->startSnpIdx
//        % incdSnpInfoVec[ldblock->startSnpIdx]->ID
//        % ldblock->endSnpIdx
//        % incdSnpInfoVec[ldblock->endSnpIdx]->ID
//        % ldblock->snpNameVec.size();
//    }
//    out2.close();
//
//    cout << "Written the LD matrix SNP info into file [" << outSnpfile << "]." << endl;
//    cout << "Written the LD matrix ldm info into file [" << outldmfile << "]." << endl;
//}


//void Data::makeBlockLDmatrix(const string &bedFile, const string &LDmatType, const unsigned block, const string &dirname, const bool writeLdmTxt, int ldBlockRegionWind){
//    int i,j;
//    vector<locus_bp> snpVec;
//    SnpInfo *snp;
//
//    map<int, string>  chrEndSnp;
//    for (i = 1; i < numIncdSnps; i++) {
//        snp = incdSnpInfoVec[i];
//        if(incdSnpInfoVec[i]->chrom != incdSnpInfoVec[i-1]->chrom){
//            chrEndSnp.insert(pair<int, string>(incdSnpInfoVec[i - 1]->chrom,incdSnpInfoVec[i - 1]->ID ));
//        }
//    }
//    chrEndSnp.insert(pair<int, string>(incdSnpInfoVec[numIncdSnps - 1]->chrom,incdSnpInfoVec[numIncdSnps - 1]->ID ));
//    //Step 1.2  Read block file
//    //readLDBlockInfoFile(ldBlockInfoFile);
//
//
//    /////////////////////////////////////////
//    // Step 2. Map snps to blocks
//    /////////////////////////////////////////
//    vector<string> block2snp_1(numLDBlocks), block2snp_2(numLDBlocks);
//    map<string,int> keptLdBlock2AllLdBlcokMap;
//    vector<locus_bp>::iterator iter;
//    map<int, string>::iterator chrIter;
//    LDBlockInfo *ldblock;
//    for (i = 0; i < numIncdSnps ; i++) {
//        snp = incdSnpInfoVec[i];
//        snpVec.push_back(locus_bp(snp->ID, snp->chrom, snp->physPos ));
//    }
//#pragma omp parallel for private(iter, chrIter)
//    for (i = 0; i < numLDBlocks; i++) {
//        // find lowest snp_name in the block
//        ldblock = ldBlockInfoVec[i];
//
//        iter = find_if(snpVec.begin(), snpVec.end(), locus_bp( ldblock->ID ,ldblock->chrom, ldblock->startPos - ldBlockRegionWind));
//        if (iter != snpVec.end()) block2snp_1[i] = iter->locusName;
//        else block2snp_1[i] = "NA";
//    }
//#pragma omp parallel for private(iter, chrIter)
//    for (i = 0; i < numLDBlocks; i++) {
//        ldblock = ldBlockInfoVec[i];
//        if (block2snp_1[i] == "NA") {
//            block2snp_2[i] = "NA";
//            continue;
//        }
//        iter = find_if(snpVec.begin(), snpVec.end(), locus_bp(ldblock->ID, ldblock->chrom, ldblock->endPos + ldBlockRegionWind));
//        if (iter != snpVec.end()){
//            if (iter->bp ==  ldblock->endPos + ldBlockRegionWind){
//                block2snp_2[i] = iter->locusName;
//            }else {
//                if(iter!=snpVec.begin()){
//                    iter--;
//                    block2snp_2[i] = iter->locusName;
//                }
//                else block2snp_2[i] = "NA";
//            }
//        }
//        else {
//            chrIter = chrEndSnp.find(ldblock->chrom);
//            if (chrIter == chrEndSnp.end()) block2snp_2[i] = "NA";
//            else block2snp_2[i] = chrIter->second;
//        }
//    }
//    int mapped = 0;
//    for (i = 0; i < numLDBlocks; i++) {
//        ldblock = ldBlockInfoVec[i];
//        if (block2snp_1[i] != "NA" && block2snp_2[i] != "NA")
//        {
//            mapped++;
//            // ldblock->kept = true;
//            keptLdBlock2AllLdBlcokMap.insert(pair<string, int>(ldblock->ID, i));
//        } else {
//            ldblock->kept = false;
//        }
//    }
//
//    if (mapped < 1) throw(0, "No SNP can be mapped to the provided LD block list. Please check the input data regarding chromosome and bp.");
//    else cout << mapped << " LD block(s) are retained." << endl;
//
//    keptLdBlockInfoVec = makeKeptLDBlockInfoVec(ldBlockInfoVec);
//    numKeptLDBlocks = (unsigned) keptLdBlockInfoVec.size();
//
//    map<string, int>::iterator iter1, iter2;
//    map<string, int> snp_name_map;
//    VectorXf snpNumInldblock(numKeptLDBlocks);
//    VectorXf eigenNumInldblock(numKeptLDBlocks);
//    vector<VectorXf> cumsumNonNeg(numKeptLDBlocks);
//
//    for (i = 0; i < numIncdSnps; i++) {
//        snp = incdSnpInfoVec[i];
//        snp_name_map.insert(pair<string,int>(snp->ID, i));
//
//    }
//    // eigenvector and eigenvalue
//    //eigenValLdBlock.resize(numKeptLDBlocks);
//    //eigenVecLdBlock.resize(numKeptLDBlocks);
//
//    //    string outfilename = filename + "." + LDmatType + ".eigen" + ".bin";
//    //    FILE *out = fopen(outfilename.c_str(), "wb");
//
//    struct stat sb;
//    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
//        // Folder doesn't exist, create it
//        string create_cmd = "mkdir " + dirname;
//        system(create_cmd.c_str());
//        cout << "Created folder [" << dirname << "] to store results." << endl;
//
//    }
//
//    string outBinfile;
//    string outTxtfile;
//    bool readBedBool = true;
//    for (i = 0; i < numKeptLDBlocks; i++) {
//        ldblock = keptLdBlockInfoVec[i];
//
//        outBinfile = dirname + "/block" + ldblock->ID + ".ldm.bin";
//        FILE *outbin = fopen(outBinfile.c_str(), "wb");
//        ofstream outtxt;
//        outTxtfile;
//        if (writeLdmTxt) {
//            outTxtfile = dirname + "/block" + ldblock->ID + ".ldm.txt";
//            outtxt.open(outTxtfile.c_str());
//        }
//
//        // cout << "ldblock id: " << ldblock->ID << endl;
//        iter1 = snp_name_map.find(block2snp_1[keptLdBlock2AllLdBlcokMap.at(ldblock->ID ) ]);
//        iter2 = snp_name_map.find(block2snp_2[keptLdBlock2AllLdBlcokMap.at(ldblock->ID ) ]);
//
//        bool skip = false;
//        if (iter1 == snp_name_map.end() || iter2 == snp_name_map.end() || iter1->second >= iter2->second) ldblock->kept = false;
//        snpNumInldblock[i] = iter2->second - iter1->second + 1;
//        // cout << "ldblock->kept: " << ldblock->kept << endl;
//        if(!ldblock->kept) continue;
//        vector<int> snp_indx;
//        SnpInfo *snp;
//        for (j = iter1->second; j <= iter2->second; j++) {
//            snp_indx.push_back(j);
//            snp = incdSnpInfoVec[j];
//            ldblock->snpNameVec.push_back(snp->ID);
//            snp->block = ldblock->ID;
//            // cout << incdSnpInfoVec[j]->ID << " ";
//        }
//        if(readBedBool) {
//            cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
//            readBedBool = false;
//        }
//        //        MatrixXf eigenVec;
//        //        VectorXf eigenVal,cumsumNonNegPerLD;
//        MatrixXf rval = generateLDmatrixPerBlock(bedFile, ldblock->snpNameVec);
//        //        // cout << "rval: " << rval << endl;
//        //        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
//        //        eigenDecomposition(rval, eigenCutoff,eigenVal, eigenVec,cumsumNonNegPerLD);
//        //        // cout << "rval: " << rval.row(0) << endl;
//        //        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
//        //        // save svd matrix
//        //        int32_t numEigenValue = eigenVal.size();
//        int32_t numSnpInBlock = ldblock->snpNameVec.size();
//
//        ldblock->startSnpIdx = snp_indx[0];
//        ldblock->endSnpIdx = snp_indx[snp_indx.size()-1];
//
//        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numSnpInBlock;
//        fwrite(rval.data(), sizeof(float), nElements, outbin);
//        //        cout << " Generate and save Eigen decomposition result for LD block " << i << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << "\r" << flush;
//
//        if (writeLdmTxt) {
//            for (unsigned ii=0; ii<numSnpInBlock; ++ii){
//                for (unsigned jj=0; jj<numSnpInBlock; ++jj) {
//                    outtxt << ldblock->ID << "\t" << ldblock->snpNameVec[ii] << "\t" << ldblock->snpNameVec[jj] << "\t" << rval(ii,jj) << endl;
//                }
//            }
//        }
//
//        fclose(outbin);
//        if (writeLdmTxt) outtxt.close();
//    }
//
//    if (block) {
//        cout << "Written the LD matrix into file [" << outBinfile << "]." << endl;
//        if (writeLdmTxt) cout << "Written the LD matrix into file [" << outTxtfile << "]." << endl;
//    }
//    else {
//        cout << "Written the LD matrix into folder [" << dirname << "/block*.ldm.bin]." << endl;
//        if (writeLdmTxt) cout << "Written the LD matrix into text file [" << dirname << "/block*.ldm.txt]." << endl;
//    }
//
//    outputBlockLDmatrixInfo(block, dirname);
//
//}


void Data::impG(const unsigned block, double diag_mod){
    VectorXi numImpSnp;
    VectorXi numValidTyp;  // included in GWAS with non-zero SE (anchors for LD imputation)
    numImpSnp.setZero(numLDBlocks);
    numValidTyp.setZero(numLDBlocks);
    for (unsigned i = 0; i < numLDBlocks; i++ ){
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (!ldblock->kept) continue;
        for (unsigned j=0; j<ldblock->numSnpInBlock; ++j) {
            SnpInfo *snp = ldblock->snpInfoVec[j];
            if (gwasUsableAnchorSnp(snp)) {
                ++numValidTyp[i];
            } else {
                // Not in GWAS, or present but SE missing/zero/non-finite — impute from LD with anchors
                ++numImpSnp[i];
            }
        }
    }
    unsigned totalNumImpSnp = numImpSnp.sum();
    
    cout << boost::format("%12s %12s %12s %12s %12s\n") % "Block" % "TotalSNPs" % "Anchors" % "ToImpute" % "PctImpute";
    for (unsigned i = 0; i < numLDBlocks; i++ ){
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        cout << boost::format("%12s %12s %12s %12s %12.3f\n") % (i+1) % ldblock->numSnpInBlock % numValidTyp[i] % numImpSnp[i] % (float(numImpSnp[i])/float(ldblock->numSnpInBlock));
        if (numValidTyp[i] == 0) {
            cout << "  Warning: Block " << (i+1) << " has no anchor SNP (finite beta and SE>0); block skipped." << endl;
            ldblock->kept = false;
        }
    }
        
    cout << "Imputing summary statistics for " << to_string(totalNumImpSnp) << " SNPs (not in GWAS and/or SE=0), using LD with SNPs that have valid SE..." << endl;

    Gadget::Timer timer;
    timer.setTime();
    
#pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < numLDBlocks; i++ ){
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (!ldblock->kept) continue;
        
        if (numImpSnp[i] && numValidTyp[i]) {
            
            Stat::Normal normal;
            
            /// Step 1. construct LD
            MatrixXf LDPerBlock = eigenVecLdBlock[i] * eigenValLdBlock[i].asDiagonal() * eigenVecLdBlock[i].transpose();
            
            LDPerBlock.diagonal().array() += (float)diag_mod;
            /// Step 2. Construct the LD correlation matrix among anchor SNPs (LDtt) and between imputation targets and anchors (LDit).
            // Anchors: gwasUsableAnchorSnp (finite b, finite SE>0).  Else imputed.
            VectorXi typedSnpIdx(numValidTyp[i]);
            VectorXi untypedSnpIdx(numImpSnp[i]);
            VectorXf zTypSnp(numValidTyp[i]);
            VectorXf nTypSnp(numValidTyp[i]);
            VectorXf varyTypSnp(numValidTyp[i]);
            for(unsigned j=0, idxTyp=0, idxImp=0; j < ldblock->numSnpInBlock; j++){
                SnpInfo *snp = ldblock->snpInfoVec[j];
                if (gwasUsableAnchorSnp(snp)) {
                    typedSnpIdx[idxTyp] = j;
                    zTypSnp[idxTyp] = float(snp->gwas_b / snp->gwas_se);
                    nTypSnp[idxTyp] = snp->gwas_n;
                    float hetj = 2.0 * snp->gwas_af * (1.0 - snp->gwas_af);
                    varyTypSnp[idxTyp] = hetj * (snp->gwas_n * snp->gwas_se * snp->gwas_se + snp->gwas_b * snp->gwas_b);
                    ++idxTyp;
                } else {
                    untypedSnpIdx[idxImp] = j;
                    ++idxImp;
                }
            }
            // Step 2.2 construct LDtt and LDit and Ztt.
            MatrixXf LDtt = LDPerBlock(typedSnpIdx,typedSnpIdx);
            MatrixXf LDit = LDPerBlock(untypedSnpIdx,typedSnpIdx);
            // Step 2.3 //  The Z score for the missing SNPs;
            VectorXf LDi_Z = LDtt.ldlt().solve(zTypSnp);
            VectorXf zImpSnp = LDit * LDi_Z;
            // Step 3. re-calcualte beta and se
            // if snp is missing use median to replace N
            std::sort(nTypSnp.data(), nTypSnp.data() + nTypSnp.size());
            float nMedian = nTypSnp[nTypSnp.size()/2];  // median
            // calcuate median of phenotypic variance
            std::sort(varyTypSnp.data(), varyTypSnp.data() + varyTypSnp.size());
            float varyMedian = varyTypSnp[varyTypSnp.size()/2];  // median
            // begin impute
            for(unsigned j = 0; j < numImpSnp[i]; j++){
                SnpInfo *snp = ldblock->snpInfoVec[untypedSnpIdx[j]];
                if (snp->keepSumstatIntact) continue;
                float base = float(sqrt(2.0 * snp->af *(1.0 - snp->af) * (nMedian + zImpSnp[j] * zImpSnp[j])));
                if (!(base > 1e-12f) || !std::isfinite(base) || !std::isfinite(zImpSnp[j])) {
                    cout << "Warning: cannot impute SNP " << snp->ID << " (MAF at 0/1 in reference or non-finite Z); excluding from analysis." << endl;
                    snp->included = false;
                    snp->gwas_b = -999;
                    snp->gwas_se = -999;
                    continue;
                }
                snp->gwas_b = zImpSnp[j] * sqrt(varyMedian)/base;
                snp->gwas_se = sqrt(varyMedian) / base;
                snp->gwas_n = nMedian;
                snp->gwas_af = snp->af;
                {
                    double denom = std::max(snp->gwas_se, 1e-20);
                    double zpv = std::abs(snp->gwas_b / denom);
                    if (!std::isfinite(zpv)) zpv = 0.0;
                    snp->gwas_pvalue = float(2*(1.0-normal.cdf_01(zpv)));
                }
                snp->included = true;
                //cout << "b " << snp->gwas_b << " se " << snp->gwas_se << " z " << snp->gwas_b/snp->gwas_se << " p " << snp->gwas_pvalue << endl;
            }
        }
        
//        if(!(i%10)) cout << " imputed block " << i << "\r" << flush;
        cout << " imputed block " << i << endl;

        if (block) {
            string outfile = title + ".block" + ldblock->ID + ".imputed.ma";
            ofstream out(outfile.c_str());
            out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n") % "SNP" % "A1" % "A2" % "freq" % "b" % "se" % "p" % "N";
            for (unsigned i=0; i<ldblock->numSnpInBlock; ++i) {
                SnpInfo *snp = ldblock->snpInfoVec[i];
                out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n")
                % snp->ID
                % snp->a1
                % snp->a2
                % snp->gwas_af
                % snp->gwas_b
                % snp->gwas_se
                % snp->gwas_pvalue
                % snp->gwas_n;
            }
            out.close();

            timer.getTime();
            cout << "Imputation of summary statistics is completed (time used: " << timer.format(timer.getElapse()) << ")." << endl;
            cout << "Summary statistics of all SNPs are saved into file [" + outfile + "]." << endl;

        }

    }

    if (block) return;
    
    string outfile = title + ".imputed.ma";
    ofstream out(outfile.c_str());
    out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n") % "SNP" % "A1" % "A2" % "freq" % "b" % "se" % "p" % "N";
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n")
        % snp->ID
        % snp->a1
        % snp->a2
        % snp->gwas_af
        % snp->gwas_b
        % snp->gwas_se
        % snp->gwas_pvalue
        % snp->gwas_n;
    }
    out.close();

    timer.getTime();
    cout << "Imputation of summary statistics is completed (time used: " << timer.format(timer.getElapse()) << ")." << endl;
    cout << "Summary statistics of all SNPs are saved into file [" + outfile + "]." << endl;
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
    for (i = 0; i < numIncdSnps ; i++) {
        snp = incdSnpInfoVec[i];
        snpVec.push_back(locus_bp(snp->ID, snp->chrom, snp->physPos ));
    }
#pragma omp parallel for private(iter, chrIter)
    for (i = 0; i < numLDBlocks; i++) {
        // find lowest snp_name in the block
        LDBlockInfo *ldblock = ldBlockInfoVec[i];

        iter = find_if(snpVec.begin(), snpVec.end(), locus_bp( ldblock->ID ,ldblock->chrom, ldblock->startPos - ldBlockRegionWind));
        if (iter != snpVec.end()) block2snp_1[i] = iter->locusName;
        else block2snp_1[i] = "NA";
    }
#pragma omp parallel for private(iter, chrIter)
    for (i = 0; i < numLDBlocks; i++) {
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
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
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
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
        SnpInfo *snp = incdSnpInfoVec[i];
        snp_name_map.insert(pair<string,int>(snp->ID, i));
        
    }
      // eigenvector and eigenvalue
    //eigenValLdBlock.resize(numKeptLDBlocks);
    //eigenVecLdBlock.resize(numKeptLDBlocks);
    string LDmatType = "block";
    string outfilename = filename + "." + LDmatType + ".eigen" + ".bin";
    FILE *out3 = fopen(outfilename.c_str(), "wb");

    bool readBedBool = true;
    for (i = 0; i < numKeptLDBlocks; i++) {
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
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
            ldblock->snpNameVec.push_back(incdSnpInfoVec[j]->ID);
           // cout << incdSnpInfoVec[j]->ID << " ";
        }
        if(readBedBool) {
            cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
            readBedBool = false;
        }
        MatrixXf eigenVec;
        VectorXf eigenVal;
        float sumPosEigVal = 0;
        MatrixXf rval = generateLDmatrixPerBlock(bedFile, ldblock->snpNameVec);
        // cout << "rval: " << rval << endl;
        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        eigenDecomposition(rval, eigenCutoff,eigenVal, eigenVec, sumPosEigVal);
        // cout << "rval: " << rval.row(0) << endl;
        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
        // save svd matrix
        int32_t numEigenValue = eigenVal.size();
        int32_t numSnpInBlock = ldblock->snpNameVec.size();
        eigenNumInldblock[i] = numEigenValue;
        // save summary
        // 1, nrow of eigenVecGene[i]
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, out3);
        //        cout << "rval: " << rval << endl;
        // cout << "eigenVec: "  << eigenVec << endl;
        // cout << "";
        // 2, ncol of eigenVecGene[i]
        fwrite(&numEigenValue, sizeof(int32_t), 1, out3);
        // 3. sum of all positive eigenvalues
        fwrite(&sumPosEigVal, sizeof(float), 1, out3);
        // 4. eigenCutoff
        fwrite(&eigenCutoff, sizeof(float), 1, out3);
        // 5. eigenvalues
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, out3);
        // 6. eigenvectors;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, out3);
        cout << " Generate and save Eigen decomposition result for LD block " << i << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << "\r" << flush;
    }
    fclose(out3);
    // cout << "size of eigenValuegene: " << eigenValGene.size() << endl;
    //outputEigenDataForLDM(filename);
    
    cout << "To explain " << eigenCutoff*100 << "% variance in LD, on average " << int(eigenNumInldblock.mean()) << " eigenvalues are selected across LD blocks (mean number of SNPs is " << int(snpNumInldblock.mean()) << ")." << endl;

}


void Data::readBlockLdmInfoFile(const string &dirname, const unsigned block){
    // Read bim file: recombination rate is defined between SNP i and SNP i-1
    string infoFile;
    if (block) {
        infoFile = dirname + "/block" + to_string(block) + ".ldm.info";
        ifstream in(infoFile.c_str());
        if (!in) infoFile = dirname + "/ldm.info";
    } else {
        infoFile = dirname + "/ldm.info";
    }
    ifstream in(infoFile.c_str());
    if (!in) throw ("Error: can not open the file [" + infoFile + "] to read.");
    cout << "Reading LDM info from file [" + infoFile + "]." << endl;
    ldBlockInfoVec.clear();
    ldBlockInfoMap.clear();
        
    string header;
    int id;
    int chr, blockStart, blockEnd, snpNum;
    int idx = 0;
    int snpCount =  1;
    string snpName;
    string startSnpID, endSnpID;
    LDBlockInfo *ldblock;
    getline(in, header);
    while (in >> id >> chr >> blockStart >> startSnpID >> blockEnd >> endSnpID >> snpNum) {
        if (block && id != block) continue;

        string idstr = to_string(id);
        ldblock = new LDBlockInfo(idx++, idstr, chr);
        ldblock->startSnpIdx = blockStart;
        ldblock->endSnpIdx   = blockEnd;
        ldblock->numSnpInBlock = snpNum;
                
        ldBlockInfoVec.push_back(ldblock);
        ldBlockInfoMap[idstr] = ldblock;

    }
    in.close();
    numLDBlocks = (unsigned) ldBlockInfoVec.size();
    cout << numLDBlocks << " LD Blocks to be included from [" + infoFile + "]." << endl;
}

void Data::readBlockLdmSnpInfoFile(const string &dirname, const unsigned block){
    string snpInfoFile;
    if (block) {
        snpInfoFile = dirname + "/block" + to_string(block) + ".snp.info";
        ifstream in(snpInfoFile.c_str());
        if (!in) snpInfoFile = dirname + "/snp.info";
    } else {
        snpInfoFile = dirname + "/snp.info";
    }
    ifstream in(snpInfoFile.c_str());
    if (!in) throw ("Error: can not open the file [" + snpInfoFile + "] to read.");
    cout << "Reading LDM SNP info from file [" + snpInfoFile + "]." << endl;
    snpInfoVec.clear();
    snpInfoMap.clear();
    
    map<string,int> ld2snpMap;
    string header;
    string id, allele1, allele2;
    int chr, physPos,ld_n;
    float genPos;
    float allele1Freq;
    int index;
    int blockID;
    LDBlockInfo *ldblock;
    getline(in, header);
    while (in >> chr >> id >> index >> genPos >> physPos >> allele1 >> allele2 >> allele1Freq >> ld_n >> blockID) {
        if (block && blockID != block) continue;
        
        string blockIDstr = to_string(blockID);

        SnpInfo *snp = new SnpInfo(index, id, allele1, allele2, chr, genPos, physPos);
        
//        cout << idx << " " << id << " " << blockID << endl;
        
        snp->af = allele1Freq;
        snp->ld_n = ld_n;
        snp->block = blockIDstr;
        snpInfoVec.push_back(snp);
        chromosomes.insert(snp->chrom);
        
        ldblock = ldBlockInfoMap[blockIDstr];
        
        if (ld2snpMap.insert(pair<string, int>(blockIDstr + "_" + id, index)).second == false) {
            throw ("Error: Duplicate LDBlock-SNP pair found: \"" + blockIDstr + "_" + id + "\".");
        } else{
            ldblock->snpNameVec.push_back(id);
            ldblock->snpInfoVec.push_back(snp);
        }

        if (snpInfoMap.insert(pair<string, SnpInfo*>(id, snp)).second == false) {
            throw ("Error: Duplicate SNP ID found: \"" + id + "\".");
        }
    }
    in.close();
    numSnps = (unsigned) snpInfoVec.size();
    
    for (unsigned i=0; i<numLDBlocks; ++i) {
        LDBlockInfo *b = ldBlockInfoVec[i];
        if (block) {
            // Filtered snp.info: snpInfoVec is only this block's SNPs; startSnpIdx/endSnpIdx in ldm.info are global indices.
            if (b->snpInfoVec.empty()) {
                throw ("Error: no SNPs read for LD block " + b->ID + "; check snp.info block IDs match ldm.info.");
            }
            b->startPos = b->snpInfoVec.front()->physPos;
            b->endPos = b->snpInfoVec.back()->physPos;
        } else {
            if ((unsigned) b->startSnpIdx >= snpInfoVec.size() || (unsigned) b->endSnpIdx >= snpInfoVec.size()) {
                throw ("Error: LD block " + b->ID + " startSnpIdx/endSnpIdx out of range for snp.info (size " + to_string(snpInfoVec.size()) + ").");
            }
            b->startPos = snpInfoVec[b->startSnpIdx]->physPos;
            b->endPos = snpInfoVec[b->endSnpIdx]->physPos;
        }
    }
    
    numKeptInds = ld_n;
    
    cout << numSnps << " SNPs to be included from [" + snpInfoFile + "]." << endl;
}

void Data::readBlockLdmBinaryAndMakeItSparse(const string &dirname, const unsigned block, const float chisqThreshold, const bool writeLdmTxt){
    struct stat sb;
    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        throw("Error: cannot find the folder [" + dirname + "]");
    }
    
    vector<int> numSnpInRegion;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        LDBlockInfo *block = ldBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }
        
    keptLdBlockInfoVec = ldBlockInfoVec;
    
    cout << "Making sparse matrices using a chisq threshold of " << chisqThreshold << endl;

#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numLDBlocks; i++){
        //if (block && i != block - 1) continue;
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[i];

        string infile = dirname + "/block" + blockInfo->ID + ".ldm.bin";
        FILE *fp = fopen(infile.c_str(), "rb");
        if(!fp){throw ("Error: can not open the file [" + infile + "] to read.");}

        string outBinfile = dirname + "/block" + blockInfo->ID + ".sparse.bin";
//        FILE *outbin = fopen(outBinfile.c_str(), "wb");

        string outTxtfile;
        ofstream outtxt;
        if (writeLdmTxt) {
            outTxtfile = dirname + "/block" + blockInfo->ID + ".sparse.txt";
//            outtxt.open(outTxtfile.c_str());
        }
        
        int32_t blockSize = blockInfo->numSnpInBlock;
        
//        MatrixXf ldm(blockSize, blockSize);
//        uint64_t nElements = (uint64_t)blockSize * (uint64_t)blockSize;
//                
//        if(fread(ldm.data(), sizeof(float), nElements, fp) != nElements){
//            cout << "fread(U.data(), sizeof(float), nElements, fp): " << fread(ldm.data(), sizeof(float), nElements, fp) << endl;
//            cout << "nEle: " << nElements << " ldm.size: " << ldm.size() <<  " ldm.col: " << ldm.cols() << " row: " << ldm.rows() << endl;
//            throw("In LD block " + blockInfo->ID + ",size error in " + outBinfile);
//            // cout << "Read " << svdLDfile << " error (U)" << endl;
//            // throw("read file error");
//        }
        
        MatrixXf ldm(blockSize, blockSize);
        readBlockLDmatrix(dirname, blockInfo->ID, blockSize, ldm);

        
//        cout << " set zero " << blockInfo->ID << endl;
        SparseMatrix<float> ldmSp = ldm.sparseView();
        ldmSp.makeCompressed();
        
//        cout << " computed block " << blockInfo->ID << endl;

        // output sparse LDM
        Gadget::writeSparseMatrixBinary(ldmSp, outBinfile);
        
//        cout << " output block " << blockInfo->ID << endl;

        if (writeLdmTxt) {
            Gadget::writeSparseMatrixToText(ldmSp, outTxtfile);
        }
        
        
        fclose(fp);
//        fclose(outbin);
//        if (writeLdmTxt) outtxt.close();

        if(!(i%1)) cout << " computed block " << blockInfo->ID << "\r" << flush;
        
        if (block) {
            cout << "Written the sparse matrix data for block LD matrix into file [" << outBinfile << "]." << endl;
            if (writeLdmTxt) cout << "Written the sparse matrix data for block LD matrix into file [" << outTxtfile << "]." << endl;
        }
    }

    if (!block) {
        cout << "Written the sparse matrix data for block LD matrix into file [" << dirname << "/block*.sparse.bin]." << endl;
        if (writeLdmTxt) cout << "Written the sparse matrix data for block LD matrix into file [" << dirname << "/block*.sparse.txt]." << endl;
    }

}

void Data::readBlockLdmBinaryAndDoEigenDecomposition(const string &dirname, const unsigned block, const float eigenCutoff, const bool writeLdmTxt){
    struct stat sb;
    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        throw("Error: cannot find the folder [" + dirname + "]");
    }
    
    vector<int> numSnpInRegion;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        LDBlockInfo *block = ldBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }
        
    keptLdBlockInfoVec = ldBlockInfoVec;
    
#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numLDBlocks; i++){
        //if (block && i != block - 1) continue;
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[i];

        string outBinfile = dirname + "/block" + blockInfo->ID + ".eigen.bin";
        FILE *outbin = fopen(outBinfile.c_str(), "wb");

        string outTxtfile;
        ofstream outtxt;
        if (writeLdmTxt) {
            outTxtfile = dirname + "/block" + blockInfo->ID + ".eigen.txt";
            outtxt.open(outTxtfile.c_str());
        }
        
        
//        MatrixXf ldm(blockSize, blockSize);
//        uint64_t nElements = (uint64_t)blockSize * (uint64_t)blockSize;
//                
//        if(fread(ldm.data(), sizeof(float), nElements, fp) != nElements){
//            cout << "fread(ldm.data(), sizeof(float), nElements, fp): " << fread(ldm.data(), sizeof(float), nElements, fp) << endl;
//            cout << "nEle: " << nElements << " ldm.size: " << ldm.size() <<  " ldm.col: " << ldm.cols() << " row: " << ldm.rows() << endl;
//            throw("In LD block " + blockInfo->ID + ",size error in " + outBinfile);
//            // cout << "Read " << svdLDfile << " error (U)" << endl;
//            // throw("read file error");
//        }
        
        MatrixXf ldm;
        readBlockLDmatrix(dirname, blockInfo->ID, blockInfo->numSnpInBlock, ldm);
                
        
        MatrixXf eigenVec;
        VectorXf eigenVal;
        float sumPosEigVal;  // sum of all positive eigenvalues
        // cout << "rval: " << rval << endl;
        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        eigenDecomposition(ldm, eigenCutoff, eigenVal, eigenVec, sumPosEigVal);
        // cout << "rval: " << rval.row(0) << endl;
        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
        // save svd matrix

        blockInfo->sumPosEigVal = sumPosEigVal;
        
        int32_t numEigenValue = eigenVal.size();
        int32_t numSnpInBlock = blockInfo->numSnpInBlock;

        //cout << " Generate Eigen decomposition result for LD block " << i << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << endl;
        
        // save summary
        // 1, the number of SNPs in the block
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, outbin);
        // 2, the number of eigenvalues at with the given cutoff
        fwrite(&numEigenValue, sizeof(int32_t), 1, outbin);
        // 3. sum of all the positive eigenvalues
        fwrite(&sumPosEigVal, sizeof(float), 1, outbin);
        // 4. eigenvalue cutoff based on the proportion of variance explained in LD
        fwrite(&eigenCutoff, sizeof(float), 1, outbin);
        // 5. the selected eigenvalues
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, outbin);
        // 6. the selected eigenvector;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, outbin);
        
        if (writeLdmTxt) {
            outtxt << "Block " << blockInfo->ID << endl;
            outtxt << "numSnps " << blockInfo->numSnpInBlock << endl;
            outtxt << "numEigenvalues " << numEigenValue << endl;
            outtxt << "SumPositiveEigenvalues " << sumPosEigVal << endl;
            outtxt << "EigenCutoff " << eigenCutoff << endl;
            outtxt << "Eigenvalues\n" << eigenVal.transpose() << endl;
            outtxt << "Eigenvectors\n" << eigenVec << endl;
            outtxt << endl;
        }
        
        fclose(outbin);
        if (writeLdmTxt) outtxt.close();

        if(!(i%1)) cout << " computed block " << blockInfo->ID << "\r" << flush;
        
        if (block) {
            cout << "Written the eigen data for block LD matrix into file [" << outBinfile << "]." << endl;
            if (writeLdmTxt) cout << "Written the eigen data for block LD matrix into file [" << outTxtfile << "]." << endl;
        }
    }

    if (!block) {
        cout << "Written the eigen data for block LD matrix into file [" << dirname << "/block*.eigen.bin]." << endl;
        if (writeLdmTxt) cout << "Written the eigen data for block LD matrix into file [" << dirname << "/block*.eigen.txt]." << endl;
    }

}

void Data::readSparseBlockLdmBinaryAndDoEigenDecomposition(const string &dirname, const unsigned block, const float eigenCutoff, const bool writeLdmTxt){
    struct stat sb;
    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        throw("Error: cannot find the folder [" + dirname + "]");
    }
    
    vector<int> numSnpInRegion;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        LDBlockInfo *block = ldBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }
        
    keptLdBlockInfoVec = ldBlockInfoVec;
    
#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numLDBlocks; i++){
        //if (block && i != block - 1) continue;
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[i];

        string infile = dirname + "/block" + blockInfo->ID + ".sparse.bin";
        FILE *fp = fopen(infile.c_str(), "rb");
        if(!fp){throw ("Error: can not open the file [" + infile + "] to read.");}

        string outBinfile = dirname + "/block" + blockInfo->ID + ".eigen.bin";
        FILE *outbin = fopen(outBinfile.c_str(), "wb");

        string outTxtfile;
        ofstream outtxt;
        if (writeLdmTxt) {
            outTxtfile = dirname + "/block" + blockInfo->ID + ".eigen.txt";
            outtxt.open(outTxtfile.c_str());
        }
        
        int32_t blockSize = blockInfo->numSnpInBlock;
        
//        MatrixXf ldm(blockSize, blockSize);
        uint64_t nElements = (uint64_t)blockSize * (uint64_t)blockSize;
                
        MatrixXf ldm = MatrixXf(Gadget::readSparseMatrixBinary(infile));
        
//        if(fread(ldm.data(), sizeof(float), nElements, fp) != nElements){
//            cout << "fread(U.data(), sizeof(float), nElements, fp): " << fread(ldm.data(), sizeof(float), nElements, fp) << endl;
//            cout << "nEle: " << nElements << " ldm.size: " << ldm.size() <<  " ldm.col: " << ldm.cols() << " row: " << ldm.rows() << endl;
//            throw("In LD block " + blockInfo->ID + ",size error in " + outBinfile);
//            // cout << "Read " << svdLDfile << " error (U)" << endl;
//            // throw("read file error");
//        }
        
        MatrixXf eigenVec;
        VectorXf eigenVal;
        float sumPosEigVal;  // sum of all positive eigenvalues
        // cout << "rval: " << rval << endl;
        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        eigenDecomposition(ldm, eigenCutoff, eigenVal, eigenVec, sumPosEigVal);
        // cout << "rval: " << rval.row(0) << endl;
        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
        // save svd matrix

        blockInfo->sumPosEigVal = sumPosEigVal;
        
        int32_t numEigenValue = eigenVal.size();
        int32_t numSnpInBlock = blockSize;

        //cout << " Generate Eigen decomposition result for LD block " << i << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << endl;
        
        // save summary
        // 1, the number of SNPs in the block
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, outbin);
        // 2, the number of eigenvalues at with the given cutoff
        fwrite(&numEigenValue, sizeof(int32_t), 1, outbin);
        // 3. sum of all the positive eigenvalues
        fwrite(&sumPosEigVal, sizeof(float), 1, outbin);
        // 4. eigenvalue cutoff based on the proportion of variance explained in LD
        fwrite(&eigenCutoff, sizeof(float), 1, outbin);
        // 5. the selected eigenvalues
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, outbin);
        // 6. the selected eigenvector;
        nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, outbin);
        
        if (writeLdmTxt) {
            outtxt << "Block " << blockInfo->ID << endl;
            outtxt << "numSnps " << numSnpInBlock << endl;
            outtxt << "numEigenvalues " << numEigenValue << endl;
            outtxt << "SumPositiveEigenvalues " << sumPosEigVal << endl;
            outtxt << "EigenCutoff " << eigenCutoff << endl;
            outtxt << "Eigenvalues\n" << eigenVal.transpose() << endl;
            outtxt << "Eigenvectors\n" << eigenVec << endl;
            outtxt << endl;
        }
        
        fclose(fp);
        fclose(outbin);
        if (writeLdmTxt) outtxt.close();

        if(!(i%1)) cout << " computed block " << blockInfo->ID << "\r" << flush;
        
        if (block) {
            cout << "Written the eigen data for block LD matrix into file [" << outBinfile << "]." << endl;
            if (writeLdmTxt) cout << "Written the eigen data for block LD matrix into file [" << outTxtfile << "]." << endl;
        }
    }

    if (!block) {
        cout << "Written the eigen data for block LD matrix into file [" << dirname << "/block*.eigen.bin]." << endl;
        if (writeLdmTxt) cout << "Written the eigen data for block LD matrix into file [" << dirname << "/block*.eigen.txt]." << endl;
    }

}

void Data::readEigenMatrixBinaryFile(const string &dirname, const float eigenCutoff, const bool writeLdmTxt, const string &outputDir){
    if (!Gadget::directoryExist(dirname)) {
        throwEigenReadErr("Error: cannot find the folder [" + dirname + "]");
    }
    
    vector<int>numSnpInRegion;
    LDBlockInfo * block;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        block = ldBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }
    eigenValLdBlock.resize(numLDBlocks);
    eigenVecLdBlock.resize(numLDBlocks);
    cout << "Reading per-block eigen binaries under [" << dirname << "] (" << numLDBlocks << " blocks, serial) ..." << endl;
    // Serial loop: exceptions must not be thrown from OpenMP parallel regions (undefined behavior with libgomp;
    // often terminates with std::string and no message on stderr).
    for(int i = 0; i < (int)numLDBlocks; i++){
        LDBlockInfo * block;
        block = ldBlockInfoVec[i];
        
        if (!block->kept) continue;
        
        int32_t cur_m = 0;
        int32_t cur_k = 0;
        float sumPosEigVal = 0;
        float oldEigenCutoff =0;
        
        string infile = dirname + "/block" + block->ID + ".eigen.bin";
        FILE *fp = fopen(infile.c_str(), "rb");
        if(!fp){
            throwEigenReadErr("Error: can not open the file [" + infile + "] to read.");
        }

        // 1. marker number
        if(fread(&cur_m, sizeof(int32_t), 1, fp) != 1){
            fclose(fp);
            throwEigenReadErr("Read " + infile + " error (m)");
        }
                
        if(cur_m != numSnpInRegion[i]){
            fclose(fp);
            throwEigenReadErr("In LD block " + block->ID + ", inconsistent marker number to marker information in " + infile);
        }
        // 2. ncol of eigenVec (number of eigenvalues)
        if(fread(&cur_k, sizeof(int32_t), 1, fp) != 1){
            fclose(fp);
            throwEigenReadErr("In LD block " + block->ID + ", error about number of eigenvalues in  " + infile);
        }
        // 3. sum of all positive eigenvalues
        if(fread(&sumPosEigVal, sizeof(float), 1, fp) != 1){
            fclose(fp);
            throwEigenReadErr("In LD block " + block->ID + ", error about the sum of positive eigenvalues in " + infile);
        }
        // 4. eigenCutoff
        if(fread(&oldEigenCutoff, sizeof(float), 1, fp) != 1){
            fclose(fp);
            throwEigenReadErr("In LD block " + block->ID + ", error about eigen cutoff used in " + infile);
        }
        // 5. eigenvalues
        VectorXf lambda(cur_k);
        if(fread(lambda.data(), sizeof(float), cur_k, fp) != (size_t)cur_k){
            fclose(fp);
            throwEigenReadErr("In LD block " + block->ID + ",size error about eigenvalues in " + infile);
        }
        // 6. eigenvector
        MatrixXf U(cur_m, cur_k);
        uint64_t nElements = (uint64_t)cur_m * (uint64_t)cur_k;
        size_t nReadU = fread(U.data(), sizeof(float), nElements, fp);
        if(nReadU != nElements){
            fclose(fp);
            cerr << "Expected " << nElements << " eigen matrix floats; got " << nReadU
                 << " (rows " << cur_m << ", cols " << cur_k << ")." << endl;
            throwEigenReadErr("In LD block " + block->ID + ",size error about eigenvectors in " + infile);
        }
        
        fclose(fp);
        
        bool haveValue = false;
        int revIdx = 0;
        if(oldEigenCutoff < eigenCutoff & i == 0){
            cout << "Warning: current proportion of variance in LD block is set as " + to_string(eigenCutoff)+ ". But the proportion of variance is set as "<< to_string(oldEigenCutoff) + " in "  + infile + ".\n";
            // throw("");
        }
        // cout << "lambda: " << lambda << endl;
        // cout << "U: " << U << endl;
        // eigenVecLdBlock[i] = U;
        // eigenValLdBlock[i] = lambda;
        
        if (eigenCutoff < oldEigenCutoff) {
            truncateEigenMatrix(sumPosEigVal, eigenCutoff, lambda, U, eigenValLdBlock[i], eigenVecLdBlock[i]);
        } else {
            eigenValLdBlock[i] = lambda;
            eigenVecLdBlock[i] = U;
        }
        block->sumPosEigVal = sumPosEigVal;
        block->eigenvalues = lambda;
        
        string outTxtfile;
        ofstream outtxt;
        if (writeLdmTxt) {
            Gadget::createDirectory(outputDir);
            outTxtfile = outputDir + "/block" + block->ID + ".eigen.txt";
            outtxt.open(outTxtfile.c_str());
            outtxt << "Block " << block->ID << endl;
            outtxt << "numSnps " << cur_m << endl;
            outtxt << "numEigenvalues " << cur_k << endl;
            outtxt << "SumPositiveEigenvalues " << sumPosEigVal << endl;
            outtxt << "EigenCutoff " << eigenCutoff << endl;
            outtxt << "Eigenvalues\n" << lambda.transpose() << endl;
            outtxt << "Eigenvectors\n" << U << endl;
            outtxt.close();
            cout << "Output eigen matrix for block " << block->ID << " in file [" << outTxtfile << "]." << endl;
        }
    }
    
    lowRankModel = true;
}

bool Data::blockLdmBinFileExists(const string &dirname, const string &blockID) const {
    string path = dirname + "/block" + blockID + ".ldm.bin";
    struct stat sb;
    return stat(path.c_str(), &sb) == 0 && S_ISREG(sb.st_mode);
}

bool Data::blockEigenBinFileExists(const string &dirname, const string &blockID) const {
    string path = dirname + "/block" + blockID + ".eigen.bin";
    struct stat sb;
    return stat(path.c_str(), &sb) == 0 && S_ISREG(sb.st_mode);
}

bool Data::readBlockEigenBinContents(
    const string &dirname,
    const string &blockID,
    int32_t expect_m,
    int32_t &cur_m,
    int32_t &cur_k,
    float &sumPosEigVal,
    float &oldEigenCutoff,
    VectorXf &lambda,
    MatrixXf &U
) {
    string infile = dirname + "/block" + blockID + ".eigen.bin";
    FILE *fp = fopen(infile.c_str(), "rb");
    if (!fp) return false;
    if (fread(&cur_m, sizeof(int32_t), 1, fp) != 1) {
        fclose(fp);
        return false;
    }
    if (cur_m != expect_m) {
        fclose(fp);
        return false;
    }
    if (fread(&cur_k, sizeof(int32_t), 1, fp) != 1) {
        fclose(fp);
        return false;
    }
    if (fread(&sumPosEigVal, sizeof(float), 1, fp) != 1) {
        fclose(fp);
        return false;
    }
    if (fread(&oldEigenCutoff, sizeof(float), 1, fp) != 1) {
        fclose(fp);
        return false;
    }
    lambda.resize(cur_k);
    if (fread(lambda.data(), sizeof(float), cur_k, fp) != (size_t)cur_k) {
        fclose(fp);
        return false;
    }
    U.resize(cur_m, cur_k);
    uint64_t nElements = (uint64_t)cur_m * (uint64_t)cur_k;
    if (fread(U.data(), sizeof(float), nElements, fp) != nElements) {
        fclose(fp);
        return false;
    }
    fclose(fp);
    return true;
}

bool Data::buildSubmatrixLdFromEigenFactors(
    const float eigenCutoff,
    float sumPosEigVal,
    float oldEigenCutoff,
    const VectorXf &lambda,
    const MatrixXf &U,
    const vector<unsigned> &kept,
    MatrixXf &RssOut
) {
    unsigned m2 = (unsigned)kept.size();
    if (m2 == 0) return false;

    VectorXf lambdaUse;
    MatrixXf Uuse;
    if (eigenCutoff < oldEigenCutoff) {
        truncateEigenMatrix(sumPosEigVal, eigenCutoff, lambda, U, lambdaUse, Uuse);
    } else {
        lambdaUse = lambda;
        Uuse = U;
    }
    int k2 = (int)lambdaUse.size();
    if (k2 <= 0) return false;

    MatrixXf U_K((int)m2, k2);
    for (unsigned a = 0; a < m2; ++a)
        U_K.row((int)a) = Uuse.row((int)kept[a]);

    RssOut = U_K * lambdaUse.asDiagonal() * U_K.transpose();
    return true;
}

bool Data::trySubLdEigenFromFullLdm(
    LDBlockInfo *block,
    const string &dirname,
    const float eigenCutoff,
    const VectorXf &gwasBlock,
    VectorXf &eigenValOut,
    MatrixXf &eigenVecOut,
    VectorXf &wcorrOut,
    MatrixXf &Qout,
    float &sumPosEigValOut,
    VectorXi &remapOut,
    vector<unsigned> &keptLocalOut,
    bool *outUsedEigenRecon
) {
    if (outUsedEigenRecon) *outUsedEigenRecon = false;

    int m = block->numSnpInBlock;
    if (m <= 1) return false;
    if (!recomputeEigen) return false;
    unsigned n_skip = 0;
    for (int j = 0; j < m; ++j)
        if (block->snpInfoVec[j]->skip) ++n_skip;
    if (n_skip == 0 || n_skip >= (unsigned)m) return false;

    vector<unsigned> kept;
    kept.reserve(m - n_skip);
    for (int j = 0; j < m; ++j)
        if (!block->snpInfoVec[j]->skip) kept.push_back((unsigned)j);

    unsigned m2 = (unsigned)kept.size();
    if (m2 == 0) return false;

    MatrixXf Rss(m2, m2);
    bool haveRss = false;

    if (blockLdmBinFileExists(dirname, block->ID)) {
        try {
            MatrixXf L;
            readBlockLDmatrix(dirname, block->ID, m, L);
            for (unsigned a = 0; a < m2; ++a)
                for (unsigned b = 0; b < m2; ++b)
                    Rss(a, b) = L((int)kept[a], (int)kept[b]);
            haveRss = true;
        } catch (...) {
            haveRss = false;
        }
    }

    if (!haveRss) {
        int32_t cur_m = 0, cur_k = 0;
        float sumPosEigVal = 0, oldEigenCutoff = 0;
        VectorXf lambda;
        MatrixXf U;
        if (!readBlockEigenBinContents(dirname, block->ID, m, cur_m, cur_k, sumPosEigVal, oldEigenCutoff, lambda, U))
            return false;
        if (!buildSubmatrixLdFromEigenFactors(eigenCutoff, sumPosEigVal, oldEigenCutoff, lambda, U, kept, Rss))
            return false;
        if (outUsedEigenRecon) *outUsedEigenRecon = true;
    }

    eigenDecomposition(Rss, eigenCutoff, eigenValOut, eigenVecOut, sumPosEigValOut);

    remapOut.resize(m);
    int lc = 0;
    for (int j = 0; j < m; ++j) {
        if (block->snpInfoVec[j]->skip) remapOut[j] = -1;
        else remapOut[j] = lc++;
    }
    keptLocalOut.swap(kept);

    VectorXf sqrtLambda = eigenValOut.array().sqrt();
    VectorXf b_kept(m2);
    for (unsigned a = 0; a < m2; ++a)
        b_kept[a] = gwasBlock[(int)keptLocalOut[a]];
    wcorrOut = (1.0f / sqrtLambda.array()).matrix().asDiagonal() * (eigenVecOut.transpose() * b_kept);
    Qout = sqrtLambda.asDiagonal() * eigenVecOut.transpose();
    return true;
}

bool Data::trySubLdEigenFromFullLdmBivariate(
    LDBlockInfo *block,
    const string &dirname,
    const float eigenCutoff,
    const VectorXf &gwas1,
    const VectorXf &gwas2,
    VectorXf &wcorrOut,
    MatrixXf &Qout,
    VectorXf &eigenValOut,
    MatrixXf &eigenVecOut,
    float &sumPosEigValOut,
    VectorXi &remapOut,
    vector<unsigned> &keptLocalOut,
    bool *outUsedEigenRecon
) {
    if (outUsedEigenRecon) *outUsedEigenRecon = false;

    int m = block->numSnpInBlock;
    if (m <= 1) return false;
    if (!recomputeEigen) return false;
    unsigned n_skip = 0;
    for (int j = 0; j < m; ++j)
        if (block->snpInfoVec[j]->skip) ++n_skip;
    if (n_skip == 0 || n_skip >= (unsigned)m) return false;

    vector<unsigned> kept;
    kept.reserve(m - n_skip);
    for (int j = 0; j < m; ++j)
        if (!block->snpInfoVec[j]->skip) kept.push_back((unsigned)j);

    unsigned m2 = (unsigned)kept.size();
    if (m2 == 0) return false;

    MatrixXf Rss(m2, m2);
    bool haveRss = false;

    if (blockLdmBinFileExists(dirname, block->ID)) {
        try {
            MatrixXf L;
            readBlockLDmatrix(dirname, block->ID, m, L);
            for (unsigned a = 0; a < m2; ++a)
                for (unsigned b = 0; b < m2; ++b)
                    Rss(a, b) = L((int)kept[a], (int)kept[b]);
            haveRss = true;
        } catch (...) {
            haveRss = false;
        }
    }

    if (!haveRss) {
        int32_t cur_m = 0, cur_k = 0;
        float sumPosEigVal = 0, oldEigenCutoff = 0;
        VectorXf lambda;
        MatrixXf U;
        if (!readBlockEigenBinContents(dirname, block->ID, m, cur_m, cur_k, sumPosEigVal, oldEigenCutoff, lambda, U))
            return false;
        if (!buildSubmatrixLdFromEigenFactors(eigenCutoff, sumPosEigVal, oldEigenCutoff, lambda, U, kept, Rss))
            return false;
        if (outUsedEigenRecon) *outUsedEigenRecon = true;
    }

    eigenDecomposition(Rss, eigenCutoff, eigenValOut, eigenVecOut, sumPosEigValOut);

    remapOut.resize(m);
    int lc = 0;
    for (int j = 0; j < m; ++j) {
        if (block->snpInfoVec[j]->skip) remapOut[j] = -1;
        else remapOut[j] = lc++;
    }
    keptLocalOut.swap(kept);

    VectorXf sqrtLambda = eigenValOut.array().sqrt();
    VectorXf b1_kept(m2), b2_kept(m2);
    for (unsigned a = 0; a < m2; ++a) {
        b1_kept[a] = gwas1[(int)keptLocalOut[a]];
        b2_kept[a] = gwas2[(int)keptLocalOut[a]];
    }
    VectorXf w1 = (1.0f / sqrtLambda.array()).matrix().asDiagonal() * (eigenVecOut.transpose() * b1_kept);
    VectorXf w2 = (1.0f / sqrtLambda.array()).matrix().asDiagonal() * (eigenVecOut.transpose() * b2_kept);
    wcorrOut.resize(w1.size() + w2.size());
    wcorrOut.head(w1.size()) = w1;
    wcorrOut.tail(w2.size()) = w2;
    Qout = sqrtLambda.asDiagonal() * eigenVecOut.transpose();
    return true;
}

void Data::readEigenMatrixBinaryFileAndMakeWandQ(const string &dirname, const float eigenCutoff, const vector<VectorXf> &GWASeffects, const VectorXf &nGWASblock, const bool noscale, const bool makePseudoSummary){
    if (!wcorrBlocks.size()) {  // only print for the first time reading the data
        cout << "Reading eigenvectors from binary file and making W and Q matrices..." << endl;
    }
    if (!Gadget::directoryExist(dirname)) {
        throw("Error: cannot find the folder [" + dirname + "]");
    }
    
    vector<int>numSnpInRegion(numKeptLDBlocks);
    
    for(int i = 0; i < numKeptLDBlocks;i++){
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }

    unsigned subLdPlanned = 0;
    vector<string> subLdMissingLdBlockIds;
    for (int i = 0; i < numKeptLDBlocks; ++i) {
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        int m = block->numSnpInBlock;
        if (m <= 1) continue;
        if (!recomputeEigen) continue;
        unsigned n_skip = 0;
        for (int j = 0; j < m; ++j)
            if (block->snpInfoVec[j]->skip) ++n_skip;
        if (n_skip == 0 || n_skip >= (unsigned)m) continue;
        if (blockLdmBinFileExists(dirname, block->ID) || blockEigenBinFileExists(dirname, block->ID))
            ++subLdPlanned;
        else
            subLdMissingLdBlockIds.push_back(block->ID);
    }
    if (!subLdMissingLdBlockIds.empty() && !wcorrBlocks.size()) {
        cout << "Warning: --recompute-eigen with --skip: partial skip in " << subLdMissingLdBlockIds.size()
             << " LD block(s); need block<ID>.ldm.bin or block<ID>.eigen.bin in [" << dirname << "]. Missing for block ID(s):";
        for (size_t t = 0; t < subLdMissingLdBlockIds.size(); ++t)
            cout << " " << subLdMissingLdBlockIds[t];
        cout << ". Falling back to full-block block*.eigen.bin for those blocks." << endl;
    }
    if (subLdPlanned && !wcorrBlocks.size()) {
        cout << "For " << subLdPlanned << " LD block(s), the LD submatrix of retained SNPs will be reconstructed (--recompute-eigen) and the GWAS effects of skipped SNPs will not be used." << endl;
    }

    eigenValLdBlock.resize(numLDBlocks);
    eigenVecLdBlock.resize(numLDBlocks);
    wcorrBlocks.resize(numKeptLDBlocks);
    numSnpsBlock.resize(numKeptLDBlocks);
    numEigenvalBlock.resize(numKeptLDBlocks);
    Qblocks.resize(numKeptLDBlocks);

    
    //Constructing pseudo summary statistics for training and validation data sets, with 90% sample size for training and 10% for validation
    
    VectorXf n_trn(numKeptLDBlocks);
    VectorXf n_val(numKeptLDBlocks);
    if (makePseudoSummary) {
        pseudoGwasNtrnBlock.resize(numKeptLDBlocks);
        for(int i = 0; i < numKeptLDBlocks;i++){
            n_trn[i] = 0.9*nGWASblock[i];
            n_val[i] = nGWASblock[i] - n_trn[i];
            pseudoGwasNtrnBlock[i] = n_trn[i];
        }

        pseudoGwasEffectTrn.resize(numKeptLDBlocks);
        pseudoGwasEffectVal.resize(numKeptLDBlocks);
        b_val.setZero(numIncdSnps);
    }
    
#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        block->eigenColRemap.resize(0);
        block->subLdKeptLocalIdx.clear();

        VectorXf eigenValSub;
        MatrixXf eigenVecSub;
        VectorXf wcorrSub;
        MatrixXf Qsub;
        float sumPosSub = 0;
        VectorXi remapSub;
        vector<unsigned> keptSub;

        bool didSubLd = trySubLdEigenFromFullLdm(block, dirname, eigenCutoff, GWASeffects[i], eigenValSub, eigenVecSub, wcorrSub, Qsub, sumPosSub, remapSub, keptSub, nullptr);

        if (didSubLd) {
#pragma omp critical
            {
                cout << "LD block " << block->ID << ": " << block->numSnpInBlock - (int)keptSub.size()
                     << " SNPs skipped, " << keptSub.size() << " retained." << endl;
            }
            eigenValLdBlock[i] = eigenValSub;
            eigenVecLdBlock[i] = eigenVecSub;
            wcorrBlocks[i] = wcorrSub;
            Qblocks[i] = Qsub;
            block->sumPosEigVal = sumPosSub;
            block->eigenvalues = eigenValSub;
            block->eigenColRemap = remapSub;
            block->subLdKeptLocalIdx = keptSub;
            numSnpsBlock[i] = Qblocks[i].cols();
            numEigenvalBlock[i] = Qblocks[i].rows();

            if (makePseudoSummary) {
                long size = eigenValLdBlock[i].size();
                VectorXf rnd(size);
                for (long j = 0; j < size; ++j) rnd[j] = Stat::snorm();
                VectorXf noise_kept = eigenVecLdBlock[i] * (eigenValLdBlock[i].array().sqrt().matrix().asDiagonal() * rnd);
                VectorXf noise_full = VectorXf::Zero(block->numSnpInBlock);
                for (unsigned a = 0; a < keptSub.size(); ++a)
                    noise_full[(int)keptSub[a]] = noise_kept[a];
                pseudoGwasEffectTrn[i] = gwasEffectInBlock[i] + sqrt(1.0f/n_trn[i] - 1.0f/nGWASblock[i]) * noise_full;
                pseudoGwasEffectVal[i] = nGWASblock[i]/n_val[i] * gwasEffectInBlock[i] - n_trn[i]/n_val[i] * pseudoGwasEffectTrn[i];
                b_val.segment(block->startSnpIdx, block->numSnpInBlock) = pseudoGwasEffectVal[i];
            }
            eigenVecLdBlock[i].resize(0,0);
            continue;
        }

        int32_t cur_m = 0;
        int32_t cur_k = 0;
        float sumPosEigVal = 0;
        float oldEigenCutoff =0;
        
        string infile = dirname + "/block" + block->ID + ".eigen.bin";
        FILE *fp = fopen(infile.c_str(), "rb");
        if(!fp){cout << "Error: can not open the file [" + infile + "] to read." << endl;}
        if(!fp){throw ("Error: can not open the file [" + infile + "] to read.");}

        // 1. marker number
        if(fread(&cur_m, sizeof(int32_t), 1, fp) != 1){
            throw("Read " + infile + " error (m)");
        }

        if(cur_m != numSnpInRegion[i]){
            throw("In LD block " + block->ID + ", inconsistent marker number to marker information in " + infile);
        }

        // 2. ncol of eigenVec (number of eigenvalues)
        if(fread(&cur_k, sizeof(int32_t), 1, fp) != 1){
            throw("In LD block " + block->ID + ", error about number of eigenvalues in  " + infile);
            // cout << "Read " << eigenBinFile << " error (k)" << endl;
            // throw("read file error");
        }

        // 3. sum of all positive eigenvalues
        if(fread(&sumPosEigVal, sizeof(float), 1, fp) != 1){
            throw("In LD block " + block->ID + ", error about the sum of positive eigenvalues in " + infile);
            // cout << "Read " << eigenBinFile << " error sumLambda" << endl;
            // throw("read file error");
        }

        // 4. eigenCutoff
        if(fread(&oldEigenCutoff, sizeof(float), 1, fp) != 1){
            throw("In LD block " + block->ID + ", error about eigen cutoff used in " + infile);
            // cout << "Read " << eigenBinFile << " error svdVarProp" << endl;
            // throw("read file error");
        }

        // 5. eigenvalues
        VectorXf lambda(cur_k);
        if(fread(lambda.data(), sizeof(float), cur_k, fp) != cur_k){
            throw("In LD block " + block->ID + ",size error about eigenvalues in " + infile);
            // cout << "Read " << eigenBinFile << " error (lambda)" << endl;
            // throw("read file error");
        }

        // 6. eigenvector
        MatrixXf U(cur_m, cur_k);
        uint64_t nElements = (uint64_t)cur_m * (uint64_t)cur_k;
        if(fread(U.data(), sizeof(float), nElements, fp) != nElements){
            cout << "fread(U.data(), sizeof(float), nElements, fp): " << fread(U.data(), sizeof(float), nElements, fp) << endl;
            cout << "nEle: " << nElements << " U.size: " << U.size() <<  " U.col: " << U.cols() << " row: " << U.rows() << endl;
            throw("In LD block " + block->ID + ",size error about eigenvectors in " + infile);
            // cout << "Read " << eigenBinFile << " error (U)" << endl;
            // throw("read file error");
        }
        
        fclose(fp);   // this is important!
        
        bool haveValue = false;
        int revIdx = 0;
        if(oldEigenCutoff < eigenCutoff & i == 0){
            cout << "Warning: current proportion of variance in LD block is set as " + to_string(eigenCutoff)+ ". But the proportion of variance is set as "<< to_string(oldEigenCutoff) + " in "  + infile + ".\n";
            // throw("");
        }

        // cout << "lambda: " << lambda << endl;
        // cout << "U: " << U << endl;
        // eigenVecLdBlock[i] = U;
        // eigenValLdBlock[i] = lambda;
        
        if (eigenCutoff < oldEigenCutoff) {
            truncateEigenMatrix(sumPosEigVal, eigenCutoff, lambda, U, eigenValLdBlock[i], eigenVecLdBlock[i]);
        } else {
            eigenValLdBlock[i] = lambda;
            eigenVecLdBlock[i] = U;
        }
        block->sumPosEigVal = sumPosEigVal;
        block->eigenvalues = lambda;
                        
        // make w and Q
        VectorXf sqrtLambda = eigenValLdBlock[i].array().sqrt();
        wcorrBlocks[i] = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (eigenVecLdBlock[i].transpose() * GWASeffects[i] );
        //MatrixXf tmpQblocks = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
        //MatrixDat matrixDat = MatrixDat(block->snpNameVec, tmpQblocks);
        Qblocks[i] = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
        
        // if (noscale) {
        //     VectorXf Dsqrt(block->numSnpInBlock);
        //     for (unsigned j=0; j<block->numSnpInBlock; ++j) {
        //         SnpInfo *snp = block->snpInfoVec[j];
        //         Dsqrt[j] = sqrt(snp->twopq)*snp->gwas_scalar;
        //         //Dsqrt[j] = snp->gwas_scalar;
        //     }
        //     Qblocks[i] = Qblocks[i] * Dsqrt.asDiagonal();
        // }

        numSnpsBlock[i] = Qblocks[i].cols();
        numEigenvalBlock[i] = Qblocks[i].rows();
                
        // make pseudo summary data
        if (makePseudoSummary) {
            long size = eigenValLdBlock[i].size();
            VectorXf rnd(size);
            for (unsigned j=0; j<size; ++j) {
                rnd[j] = Stat::snorm();
            }
            
            pseudoGwasEffectTrn[i] = gwasEffectInBlock[i] + sqrt(1.0/n_trn[i] - 1.0/nGWASblock[i]) * eigenVecLdBlock[i] * (eigenValLdBlock[i].array().sqrt().matrix().asDiagonal() * rnd);

            pseudoGwasEffectVal[i] = nGWASblock[i]/n_val[i] * gwasEffectInBlock[i] - n_trn[i]/n_val[i] * pseudoGwasEffectTrn[i];
            b_val.segment(block->startSnpIdx, block->numSnpInBlock) = pseudoGwasEffectVal[i];
        }
        
        eigenVecLdBlock[i].resize(0,0);
    }
}

void Data::truncateEigenMatrix(const float sumPosEigVal, const float eigenCutoff, const VectorXf &oriEigenVal, const MatrixXf &oriEigenVec, VectorXf &newEigenVal, MatrixXf &newEigenVec){
    int revIdx = (int)oriEigenVal.size();
    VectorXf cumsumNonNeg(revIdx);
    cumsumNonNeg.setZero();
    revIdx = revIdx -1;
    if(oriEigenVal(revIdx) < 0) cout << "Error, all eigenvector are negative" << endl;
    cumsumNonNeg(revIdx) = oriEigenVal(revIdx);
    revIdx = revIdx -1;
    
    while (revIdx >= 0 && oriEigenVal(revIdx) > 1e-10 ){
        cumsumNonNeg(revIdx) = oriEigenVal(revIdx) + cumsumNonNeg(revIdx + 1);
        revIdx =revIdx - 1;
        if(revIdx < 0) break;
    }
    // cout << "revIdx: " << revIdx << endl;
    // cout << "size: " << eigenVal.size()  << " eigenVal: " << eigenVal << endl;
    // cout << "cumsumNoNeg: " << cumsumNonNeg << endl;
    cumsumNonNeg = cumsumNonNeg/sumPosEigVal;  // calcualte the cumulative proportion of variance explained
    bool haveValue = false;
    // cout << "cumsumNonNeg: " << cumsumNonNeg << endl;
    // cout << "revIdx : " << revIdx << endl;
    // cout << "revIdx: "  << revIdx  << endl;
    for (revIdx = revIdx + 1; revIdx < oriEigenVal.size(); revIdx ++ ){
        // cout << "revIdx: " << revIdx << " cumsumNonNeg: " << cumsumNonNeg(revIdx) << endl;
        if(eigenCutoff >= cumsumNonNeg(revIdx) ){
            revIdx = revIdx -1;
            haveValue = true;
            break;
        }
    }
    // cout << "revIdx : " << revIdx << endl;
    if(!haveValue) revIdx = oriEigenVal.size() - 1;
    // cout << "cumsumNonNeg: " << cumsumNonNeg.size() << endl;
    // cout << "cumsumNoNeg: " << cumsumNonNeg << endl;
    // cout << "revIdx: "  << revIdx  << endl;

    newEigenVec = oriEigenVec.rightCols(oriEigenVal.size() - revIdx);
    newEigenVal = oriEigenVal.tail(oriEigenVal.size() - revIdx);
    // cout << "eigenValAdjusted size: " << eigenValAdjusted.size() << " eigenValue eventually: " << eigenValAdjusted << endl;
    //eigenvalueNum = eigenVal.size() - revIdx;
    // cout << endl;

}

void Data::readBlockLDmatrixAndMakeItSparse(const string &dirname, const unsigned block, const float chisqThreshold, const bool writeLdmTxt){
    readBlockLdmInfoFile(dirname, block);
    readBlockLdmSnpInfoFile(dirname, block);
    readBlockLdmBinaryAndMakeItSparse(dirname, block, chisqThreshold, writeLdmTxt);
}

void Data::readBlockLDmatrixAndDoEigenDecomposition(const string &dirname, const unsigned block, const float eigenCutoff, const bool writeLdmTxt){
    readBlockLdmInfoFile(dirname, block);
    readBlockLdmSnpInfoFile(dirname, block);
    readBlockLdmBinaryAndDoEigenDecomposition(dirname, block, eigenCutoff, writeLdmTxt);
}

void Data::readSparseBlockLDmatrixAndDoEigenDecomposition(const string &dirname, const unsigned block, const float eigenCutoff, const bool writeLdmTxt){
    readBlockLdmInfoFile(dirname, block);
    readBlockLdmSnpInfoFile(dirname, block);
    readSparseBlockLdmBinaryAndDoEigenDecomposition(dirname, block, eigenCutoff, writeLdmTxt);
}

void Data::readEigenMatrix(const string &dirname, const float eigenCutoff, const bool readBinary, const bool writeLdmTxt, const string &outputDir){
    cout << "Reading LD matrix eigen-decomposition data..." << endl;
    //Gadget::Timer timer;
    //timer.setTime();
    
    readBlockLdmInfoFile(dirname);
    readBlockLdmSnpInfoFile(dirname);
    if (readBinary) readEigenMatrixBinaryFile(dirname, eigenCutoff, writeLdmTxt, outputDir);
    
    //timer.getTime();
    //cout << "Read LD data completed (time used: " << timer.format(timer.getElapse()) << ")." << endl;
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

void Data::buildMMEeigen(const string &dirname, const bool sampleOverlap, const float eigenCutoff, const bool noscale){
    includeMatchedBlocks();
    
//    for (unsigned i=0; i<numIncdSnps; ++i) {
//        SnpInfo *snp = incdSnpInfoVec[i];
//        if (snp->gwas_b == -999) {
//            throw("Error: SNP " + snp->ID + " in the LD reference has no summary data. Run --impute-summary first.");
//        }
//    }
    unsigned nmiss = 0;
    for (unsigned i=0; i<numKeptLDBlocks; ++i) {
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        for (unsigned j=0; j<block->numSnpInBlock; ++j) {
            SnpInfo *snp = block->snpInfoVec[j];
            if (snp->gwas_b == -999) {
//                cout << "Error: SNP " + snp->ID + " in the LD reference has no summary data." << endl;
                ++nmiss;
            }
        }
    }
    if (nmiss) {
        string msg = "Error: " + to_string(nmiss) + " SNPs in the LD reference have no summary data after imputation. Run --impute-summary first, or check GWAS SE (finite, >0) and blocks with no anchor SNPs.";
        cerr << msg << endl;
        throw std::runtime_error(msg);
    }
    
    scaleGwasEffects();
    
    ldBlockInfoMap.clear();
//    snpInfoMap.clear();  // free memory from snpInfoMap which is no longer needed
    
    readEigenMatrixBinaryFileAndMakeWandQ(dirname, eigenCutoff, gwasEffectInBlock, nGWASblock, noscale, true);
    //constructPseudoSummaryData();  // for finding the best eigen cutoff by pseudo validation
    //if (numIncdSnps!=0) constructWandQ(gwasEffectInBlock, numKeptInds);
    //if (numIncdSnps!=0) constructWandQ(eigenCutoff, noscale);

    cout << "\nData summary:" << endl;
    cout << boost::format("%40s %8s %8s\n") %"" %"mean" %"sd";
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP Phenotypic variance" %Gadget::calcMean(varySnp) %sqrt(Gadget::calcVariance(varySnp));
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP heterozygosity" %Gadget::calcMean(snp2pq) %sqrt(Gadget::calcVariance(snp2pq));
    cout << boost::format("%40s %8.0f %8.0f\n") %"GWAS SNP sample size" %Gadget::calcMean(n) %sqrt(Gadget::calcVariance(n));
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP effect (in genotype SD unit)" %Gadget::calcMean(b) %sqrt(Gadget::calcVariance(b));
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP SE" %Gadget::calcMean(se) %sqrt(Gadget::calcVariance(se));
    cout << boost::format("%40s %8.0f %8.0f\n") %"LD block size" %Gadget::calcMean(numSnpsBlock) %sqrt(Gadget::calcVariance(numSnpsBlock));
    cout << boost::format("%40s %8.0f %8.0f\n") %"LD block rank" %Gadget::calcMean(numEigenvalBlock) %sqrt(Gadget::calcVariance(numEigenvalBlock));
    cout << endl;

    if (numAnnos) setAnnoInfoVec();
    lowRankModel = true;

}

void Data::buildMMEeigenBivariate(const string &dirname, const bool sampleOverlap, const float eigenCutoff, const bool noscale){
    includeMatchedBlocks();
    
    // Check for missing data for both traits
    unsigned nmiss1 = 0, nmiss2 = 0;
    for (unsigned i=0; i<numKeptLDBlocks; ++i) {
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        for (unsigned j=0; j<block->numSnpInBlock; ++j) {
            SnpInfo *snp = block->snpInfoVec[j];
            if (snp->gwas_b == -999) ++nmiss1;
            if (snp->gwas_b2 == -999) ++nmiss2;
        }
    }
    if (nmiss1) {
        string msg = "Error: " + to_string(nmiss1) + " SNPs in the LD reference have no summary data for trait 1. Run --impute-summary first, or check GWAS SE (finite, >0) and blocks with no anchor SNPs.";
        cerr << msg << endl;
        throw std::runtime_error(msg);
    }
    if (nmiss2) {
        string msg = "Error: " + to_string(nmiss2) + " SNPs in the LD reference have no summary data for trait 2. Run --impute-summary first, or check GWAS SE (finite, >0) and blocks with no anchor SNPs.";
        cerr << msg << endl;
        throw std::runtime_error(msg);
    }
    
    scaleBivariateGwasEffects();
    
    ldBlockInfoMap.clear();
    
    // Create gwasEffectInBlock for both traits
    vector<VectorXf> gwasEffectInBlock1, gwasEffectInBlock2;
    vector<Vector2f> nGWASblockBivariate;
    
    if (numKeptLDBlocks) {
        gwasEffectInBlock1.resize(numKeptLDBlocks);
        gwasEffectInBlock2.resize(numKeptLDBlocks);
        nGWASblockBivariate.resize(numKeptLDBlocks);
        
        for (unsigned i = 0; i < numKeptLDBlocks; i++){
            LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
            gwasEffectInBlock1[i] = b(ldblock->block2GwasSnpVec);
            gwasEffectInBlock2[i] = b2(ldblock->block2GwasSnpVec);
            
            VectorXf nBlock1 = n(ldblock->block2GwasSnpVec);
            VectorXf nBlock2 = n2(ldblock->block2GwasSnpVec);
            VectorXf nBlock1Srt = nBlock1;
            VectorXf nBlock2Srt = nBlock2;
            std::sort(nBlock1Srt.data(), nBlock1Srt.data() + nBlock1Srt.size());
            std::sort(nBlock2Srt.data(), nBlock2Srt.data() + nBlock2Srt.size());
            nGWASblockBivariate[i] << nBlock1Srt[nBlock1Srt.size()/2], nBlock2Srt[nBlock2Srt.size()/2]; // median for each trait
        }
    }
    
    readEigenMatrixBinaryFileAndMakeWandQBivariate(dirname, eigenCutoff, gwasEffectInBlock1, gwasEffectInBlock2, nGWASblockBivariate, noscale);
    
    // Store nGWASblockBivariate for later use (member variable)
    this->nGWASblockBivariate = nGWASblockBivariate;
    // Also set nGWASblock to trait 1 values for compatibility
    nGWASblock.resize(nGWASblockBivariate.size());
    for (unsigned i = 0; i < nGWASblockBivariate.size(); ++i) {
        nGWASblock[i] = nGWASblockBivariate[i][0];
    }

    cout << "\nBivariate data summary:" << endl;
    cout << boost::format("%40s %8s %8s\n") %"" %"mean" %"sd";
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP Phenotypic variance (Trait 1)" %Gadget::calcMean(varySnp)  %sqrt(Gadget::calcVariance(varySnp));
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP Phenotypic variance (Trait 2)" %Gadget::calcMean(varySnp2) %sqrt(Gadget::calcVariance(varySnp2));
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP heterozygosity"               %Gadget::calcMean(snp2pq)   %sqrt(Gadget::calcVariance(snp2pq));
    cout << boost::format("%40s %8.0f %8.0f\n") %"GWAS SNP sample size (Trait 1)"        %Gadget::calcMean(n)        %sqrt(Gadget::calcVariance(n));
    cout << boost::format("%40s %8.0f %8.0f\n") %"GWAS SNP sample size (Trait 2)"        %Gadget::calcMean(n2)       %sqrt(Gadget::calcVariance(n2));
    cout << boost::format("%40s %8.0f %8.0f\n") %"LD block size"                         %Gadget::calcMean(numSnpsBlock)    %sqrt(Gadget::calcVariance(numSnpsBlock));
    cout << boost::format("%40s %8.0f %8.0f\n") %"LD block rank"                         %Gadget::calcMean(numEigenvalBlock) %sqrt(Gadget::calcVariance(numEigenvalBlock));
    cout << endl;

    if (numAnnos) setAnnoInfoVec();
    lowRankModel = true;
}

void Data::includeMatchedBlocks(){
    // this step is to construct gwasSnp2geneVec
//    cout << "Matching blocks..." << endl;
    SnpInfo * snp;
    LDBlockInfo * ldblock;
    
    for (unsigned i=0; i<numLDBlocks; ++i){
        ldblock = ldBlockInfoVec[i];
        ldblock->block2GwasSnpVec.clear();
    }
    for (unsigned j=0; j<numIncdSnps; ++j){
        snp = incdSnpInfoVec[j];
        ldblock = ldBlockInfoMap[snp->block];
//        cout << j << " " << snp->ID << " " << snp->block << endl;
        ldblock->block2GwasSnpVec.push_back(j);
        ldblock->snpInfoVec.push_back(snp);
    }
    for (unsigned i=0; i<numLDBlocks; ++i){
        ldblock = ldBlockInfoVec[i];
        if(ldblock->block2GwasSnpVec.size() == 0){
            ldblock->kept = false;
        } else {
            ldblock->startSnpIdx = ldblock->block2GwasSnpVec[0];
            ldblock->endSnpIdx = ldblock->block2GwasSnpVec[ldblock->numSnpInBlock-1];
            ldblock->kept = true;
        }
    }
        
    keptLdBlockInfoVec = makeKeptLDBlockInfoVec(ldBlockInfoVec);
    numKeptLDBlocks = (unsigned) keptLdBlockInfoVec.size();

    ldblock2gwasSnpMap.clear();
//    cout << "Construct map from ld to snp" << endl;
    for(unsigned i = 0; i < numKeptLDBlocks; i++){
        ldblock = keptLdBlockInfoVec[i];
        ldblock2gwasSnpMap.insert(pair<int, vector<int> > (i,ldblock->block2GwasSnpVec));
    }
    
    for (unsigned j=0; j<numIncdSnps; ++j){  // update each SNP's LD block index
        snp = incdSnpInfoVec[j];
        ldblock = ldBlockInfoMap[snp->block];
        snp->blockIdx = ldblock->index;
    }
    
    cout << numKeptLDBlocks << " LD blocks are included." << endl;
}

//void Data::constructWandQ(const float eigenCutoff, const bool noscale){
//    VectorXf nMinusOne;
//    snp2pq.resize(numIncdSnps);
//    D.resize(numIncdSnps);
//   // ZPZdiag.resize(numIncdSnps);
//    ZPy.resize(numIncdSnps);
//    b.resize(numIncdSnps);
//    n.resize(numIncdSnps);
//    nMinusOne.resize(numIncdSnps);
//    se.resize(numIncdSnps);
//    tss.resize(numIncdSnps);
//    SnpInfo *snp;
//    for (unsigned i=0; i<numIncdSnps; ++i) {
//        snp = incdSnpInfoVec[i];
//        snp->af = snp->gwas_af;
//        snp2pq[i] = snp->twopq = 2.0f*snp->gwas_af*(1.0f-snp->gwas_af);
//        if(snp2pq[i]==0) cout << "Error: SNP " << snp->ID << " af " << snp->af << " has 2pq = 0." << endl;
//        D[i] = snp2pq[i]*snp->gwas_n;
//        b[i] = snp->gwas_b * sqrt(snp2pq[i]); // scale the marginal effect so that it's in per genotype SD unit
//        n[i] = snp->gwas_n;
//        nMinusOne[i] = snp->gwas_n - 1;
//        se[i]= snp->gwas_se * sqrt(snp2pq[i]);
//        tss[i] = D[i]*(n[i]*se[i]*se[i] + b[i]*b[i]);
//        ZPy[i] = n[i]*b[i];
//        //  D[i] = 1.0/(se[i]*se[i]+b[i]*b[i]/snp->gwas_n);  // NEW!
//        //  snp2pq[i] = snp->twopq = D[i]/snp->gwas_n;       // NEW!
//    }
//    cout << endl;
//    LDBlockInfo * ldblock;
//    wcorrBlocks.resize(numKeptLDBlocks);
//    // Qblocks.resize(numKeptLDBlocks);
//    numSnpsBlock.resize(numKeptLDBlocks);
//    numEigenvalBlock.resize(numKeptLDBlocks);
//    Qblocks.clear();
//    VectorXf sqrtLambda;
//    // save gwas marginal effect into block
//    gwasEffectInBlock.resize(numKeptLDBlocks);
//    for (unsigned i = 0; i < numKeptLDBlocks; i++){
//        ldblock = keptLdBlockInfoVec[i];
//        gwasEffectInBlock[i] = b(ldblock->block2GwasSnpVec);
//        // calculate wbcorr and Qblocks
//        sqrtLambda = eigenValLdBlock[i].array().sqrt();
//        //cout << eigenVecLdBlock[i].transpose().rows() << " " << eigenVecLdBlock[i].transpose().cols() << " " << gwasEffectInBlock[i].size() << endl;
//        wcorrBlocks[i] = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (eigenVecLdBlock[i].transpose() * gwasEffectInBlock[i] );
//        // cout << "eigenVecLdBlock[i]: " << eigenVecLdBlock[i] << endl;
//        // cout << "wcorrBlocks[i]: " << wcorrBlocks[i] << endl;
//        // cout << "sqrtLambda: " << sqrtLambda << endl;
//        // cout << gwasEffectInBlock[i] << endl;
//        MatrixXf tmpQblocks = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
//        MatrixDat matrixDat = MatrixDat(ldblock->snpNameVec,tmpQblocks );
//        // cout << "Qblock: " << endl;
//        // cout << matrixDat.values << endl;
//        Qblocks.push_back(matrixDat);
//        numSnpsBlock[i] = Qblocks[i].ncol;
//        numEigenvalBlock[i] = Qblocks[i].nrow;
//    }
//    
//    //b.array() -= b.mean();  // DO NOT CENTER b
//    // estimate phenotypic variance based on the input allele frequencies in GWAS
//
//    //  Vp_buf = h_buf * N_buf * se_buf * se_buf + h_buf * b_buf * b_buf * N_buf / (N_buf - 1.0);
//    //VectorXf ypySrt = D.array()*n.array()*se.array().square() + D.array() * b.array().square() * n.array() / nMinusOne.array();
//    VectorXf ypySrt = D.array()*(n.array()*se.array().square()+b.array().square());
//    VectorXf varpSrt = ypySrt.array()/n.array();
//    std::sort(ypySrt.data(), ypySrt.data() + ypySrt.size());
//    std::sort(varpSrt.data(), varpSrt.data() + varpSrt.size());
//    ypy = ypySrt[ypySrt.size()/2];  // median
//    varPhenotypic = varpSrt[varpSrt.size()/2];
//    //cout << "varPhenotypic: " << varPhenotypic << endl;
//    VectorXf nSrt = n;
//    std::sort(nSrt.data(), nSrt.data() + nSrt.size());
//    numKeptInds = nSrt[nSrt.size()/2]; // median
//    
//    nGWASblock.resize(numKeptLDBlocks);
//    for (unsigned i=0; i<numKeptLDBlocks; ++i) {
//        nGWASblock[i] = numKeptInds;
//    }
//
//    for (unsigned i=0; i<numIncdSnps; ++i) {
//        snp = incdSnpInfoVec[i];
//        D[i] = varPhenotypic/(se[i]*se[i]+b[i]*b[i]/snp->gwas_n);  // NEW!
//        snp2pq[i] = snp->twopq = D[i]/snp->gwas_n;       // NEW!
//        tss[i] = D[i]*(n[i]*se[i]*se[i] + b[i]*b[i]);
//        // Need to adjust R and C models X'X matrix depending scale of genotypes or not
//        if (noscale == true) {
//            D[i] = snp2pq[i]*snp->gwas_n;
//        } else {
//            D[i] = snp->gwas_n;
//        }
//    }
//    // cout << endl << "snp2pq: " << endl << snp2pq << endl;
//    //ypy = numKeptInds;
//    // NEW END
//    // data summary
//    cout << "\nData summary:" << endl;
//    cout << boost::format("%40s %8s %8s\n") %"" %"mean" %"sd";
//    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP Phenotypic variance" %Gadget::calcMean(varpSrt) %sqrt(Gadget::calcVariance(varpSrt));
//    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP heterozygosity" %Gadget::calcMean(snp2pq) %sqrt(Gadget::calcVariance(snp2pq));
//    cout << boost::format("%40s %8.0f %8.0f\n") %"GWAS SNP sample size" %Gadget::calcMean(n) %sqrt(Gadget::calcVariance(n));
//    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP effect (in genotype SD unit)" %Gadget::calcMean(b) %sqrt(Gadget::calcVariance(b));
//    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP SE" %Gadget::calcMean(se) %sqrt(Gadget::calcVariance(se));
//    cout << boost::format("%40s %8.3f %8.3f\n") %"LD block size" %Gadget::calcMean(numSnpsBlock) %sqrt(Gadget::calcVariance(numSnpsBlock));
//    cout << boost::format("%40s %8.3f %8.3f\n") %"LD block rank" %Gadget::calcMean(numEigenvalBlock) %sqrt(Gadget::calcVariance(numEigenvalBlock));
//}


void Data::mergeLdmInfo(const string &outLDmatType, const string &dirname, const bool print) {
    if (outLDmatType != "block") throw("Error: --merge-ldm-info only works for block LD matrices at the moment!");
    
    string dir_path = dirname; // Replace with your folder path
    string search_str = "snp.info";
    DIR* dirp = opendir(dir_path.c_str());
    
    if (dirp == NULL) {
        throw("Error opening directory [" + dirname + "]");
    }
    
    // find out all ldm in the folder
    vector<string> file_list;
    
    dirent* dp;
    while ((dp = readdir(dirp)) != NULL) {
        string file_name = dp->d_name;
        if (file_name.find(search_str) != string::npos) {
            file_list.push_back(file_name);
        }
    }
    
    closedir(dirp);
    
    if (file_list.size() == 1) {
        if (print) throw ("Error: there is only one info file in folder [" + dirname + "], so no need to merge.");
        else return;
    }

    search_str = "block";

    vector<string> filtered_list;
    for (vector<string>::iterator it = file_list.begin(); it != file_list.end(); ++it) {
        if (it->find("block") != string::npos) {
            filtered_list.push_back(*it);
        }
    }
    file_list.swap(filtered_list); // replace original with filtered

    
    set<unsigned> blockIdxSet;
    blockIdxSet.clear();

    for (vector<string>::iterator it = file_list.begin(); it != file_list.end(); ++it) {
        size_t block_pos = it->find(search_str);
        size_t dot_pos = it->find_first_of(".", block_pos);
        if (block_pos != string::npos && dot_pos != string::npos) {
            string block_num_str = it->substr(block_pos + search_str.size(), dot_pos - block_pos - search_str.size());
            int block_num = atoi(block_num_str.c_str());
            blockIdxSet.insert(block_num);
        }
    }
        
    unsigned nldm = blockIdxSet.size();
    
    cout << "Merging " << nldm << " LD matrices info files" << endl;
    
    set<unsigned>::iterator it = blockIdxSet.begin();
    
    string outSnpInfoFile = dirname + "/snp.info";
    string outldmInfoFile = dirname + "/ldm.info";
    string outpwldFile    = dirname + "/rsq0.5.pwld";

    // If a combined rsq0.5.pwld already exists (e.g. built across blocks elsewhere), merge only
    // snp.info / ldm.info from per-block files and keep the existing pwld (do not require block*.rsq0.5.pwld).
    bool useExistingCombinedPwld = false;
    {
        ifstream ex(outpwldFile.c_str(), ios::ate | ios::binary);
        if (ex.is_open() && static_cast<long long>(ex.tellg()) > 0) {
            useExistingCombinedPwld = true;
        }
    }
    if (useExistingCombinedPwld) {
        cout << "Found existing combined pairwise LD file [" + outpwldFile + "]; merging snp/ldm info only and keeping this pwld file." << endl;
    }

    ofstream out1(outSnpInfoFile.c_str());
    out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %12s %10s %10s\n")
    % "Chrom"
    % "ID"
    % "Index"
    % "GenPos"
    % "PhysPos"
    % "A1"
    % "A2"
    % "A1Freq"
    % "N"
    % "Block";
    
    ofstream out2(outldmInfoFile.c_str());
    out2 << boost::format("%10s %6s %15s %15s %15s %15s %12s\n")
    % "Block"
    % "Chrom"
    % "StartSnpIdx"
    % "StartSnpID"
    % "EndSnpIdx"
    % "EndSnpID"
    % "NumSnps";

    ofstream out3;
    if (!useExistingCombinedPwld) {
        out3.open(outpwldFile.c_str());
        out3 << boost::format("%12s %12s %12s\n")
        % "SNP1"
        % "SNP2"
        % "LDcorrelation";
    }

    unsigned snpIdx = 0;
    unsigned ldmIdx = 0;
    
    for (unsigned i=0; i<nldm; ++i) {
        
        string snpInfoFile = dirname + "/block" + to_string(*it) + ".snp.info";
        string ldmInfoFile = dirname + "/block" + to_string(*it) + ".ldm.info";
        string pwldFile = dirname + "/block" + to_string(*it) + ".rsq0.5.pwld";

        // read snp info file
        ifstream in1(snpInfoFile.c_str());
        if (!in1) throw ("Error: can not open the file [" + snpInfoFile + "] to read.");
        //cout << "Reading SNP info from file [" + snpInfoFile + "]." << endl;
        string header;
        string id, allele1, allele2;
        int chr, physPos,ld_n;
        float genPos;
        float allele1Freq;
        int idx = 0;
        int index;
        string blockID;
        map<string, int> snpID2index;
        getline(in1, header);
        while (in1 >> chr >> id >> index >> genPos >> physPos >> allele1 >> allele2 >> allele1Freq >> ld_n >> blockID) {
            out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %12f %10s %10s\n")
            % chr
            % id
            % snpIdx
            % genPos
            % physPos
            % allele1
            % allele2
            % allele1Freq
            % ld_n
            % blockID;
            snpID2index[id] = snpIdx;
            ++snpIdx;
        }
        in1.close();
        
        
        // read ldm info file
        ifstream in2(ldmInfoFile.c_str());
        if (!in2) throw ("Error: can not open the file [" + ldmInfoFile + "] to read.");
        //cout << "Reading LDM info from file [" + ldmInfoFile + "]." << endl;
        
        int blockStart, blockEnd, snpNum;
        string startSnpID, endSnpID;
        getline(in2, header);
        while (in2 >> id >> chr >> blockStart >> startSnpID >> blockEnd >> endSnpID >> snpNum) {
            out2 << boost::format("%10s %6s %15s %15s %15s %15s %12s\n")
            % id
            % chr
            % snpID2index[startSnpID]
            % startSnpID
            % snpID2index[endSnpID]
            % endSnpID
            % snpNum;
            
            ++ldmIdx;
        }
        in2.close();
        
        
        // read per-block pwld only when building merged rsq0.5.pwld (not when keeping existing combined file)
        if (!useExistingCombinedPwld) {
            ifstream in3(pwldFile.c_str());
            if (!in3) {
                cout << "Warning: can not open the file [" + pwldFile + "] to read; skipping pairwise LD for this block." << endl;
            } else {
                float ldcor;
                string snp1ID, snp2ID;
                getline(in3, header);
                while (in3 >> snp1ID >> snp2ID >> ldcor) {
                    out3 << boost::format("%12s %12s %12.6f\n")
                    % snp1ID
                    % snp2ID
                    % ldcor;
                }
                in3.close();
                remove(pwldFile.c_str());
            }
        }

        ++it;
        
        remove(snpInfoFile.c_str());
        remove(ldmInfoFile.c_str());
    }
    
    out1.close();
    out2.close();
    if (!useExistingCombinedPwld) {
        out3.close();
    }

    cout << "Written " << snpIdx << " SNPs info into file [" + outSnpInfoFile + "]." << endl;
    cout << "Written " << ldmIdx << " LDMs info into file [" + outldmInfoFile + "]." << endl;
    if (useExistingCombinedPwld) {
        cout << "Left pairwise LD file unchanged: [" + outpwldFile + "]." << endl;
    } else {
        cout << "Written pairwise LD correlations (rsq>0.5) into file [" + outpwldFile + "]." << endl;
    }
}

void Data::mergeBlockGwasSummary(const string &gwasSummaryFile, const string &title) {
    string dir_path = "."; // Replace with your folder path
    string search_str1 = gwasSummaryFile + ".block";
    string search_str2 = ".ma";
    DIR* dirp = opendir(dir_path.c_str());
    
    if (dirp == NULL) {
        throw("Error opening directory [" + dir_path + "]");
    }
    
    // find out all .ma files in the folder
    vector<string> file_list;
    
    dirent* dp;
    while ((dp = readdir(dirp)) != NULL) {
        string file_name = dp->d_name;
        if (file_name.find(search_str1) != string::npos && file_name.find(search_str2) != string::npos) {
            file_list.push_back(file_name);
        }
    }
    
    closedir(dirp);
    
    map<unsigned, string> blockIdxMap;
    blockIdxMap.clear();
    
    for (vector<string>::iterator it = file_list.begin(); it != file_list.end(); ++it) {
        size_t block_pos = it->find(search_str1);
        size_t dot_pos = it->find_first_of(".", block_pos);
        if (block_pos != string::npos && dot_pos != string::npos) {
            string block_num_str = it->substr(block_pos + search_str1.size(), dot_pos - block_pos - search_str1.size());
            int block_num = atoi(block_num_str.c_str());
            blockIdxMap[block_num] = *it;
        }
    }
    
    if (blockIdxMap.size() == 0) {
        throw ("Error: there is no info file to merge in folder [" + dir_path + "].");
    }
    
    unsigned nBlk = blockIdxMap.size();
    
    cout << "Merging GWAS summary statistics files across " + to_string(nBlk) + " blocks..." << endl;
    
    map<unsigned, string>::iterator it = blockIdxMap.begin();
    
    string outMaFile = dir_path + "/" + title + ".ma";
    
    ofstream out(outMaFile.c_str());
    out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n") % "SNP" % "A1" % "A2" % "freq" % "b" % "se" % "p" % "N";
    
    unsigned snpIdx = 0;
    
    for (unsigned i=0; i<nBlk; ++i) {
        
        string mafile = it->second;
        
        // read snp info file
        ifstream in(mafile.c_str());
        if (!in) throw ("Error: can not open the file [" + mafile + "] to read.");
        cout << "Reading summary statistics from file [" + mafile + "]." << endl;

        string header;
        getline(in, header);

        string id, allele1, allele2, freq, b, se, pval, n;
        while (in >> id >> allele1 >> allele2 >> freq >> b >> se >> pval >> n) {
            out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n")
            % id
            % allele1
            % allele2
            % freq
            % b
            % se
            % pval
            % n;
            ++snpIdx;
        }
        in.close();
        
        ++it;
    }
    
    out.close();
    
    cout << "Written " << snpIdx << " SNPs info into file [" + outMaFile + "]." << endl;
}

void Data::constructPseudoSummaryData(){
    cout << "Constructing pseudo summary statistics for training and validation data sets, with 90% sample size for training and 10% for validation." << endl;
    
    pseudoGwasEffectTrn.resize(numKeptLDBlocks);
    pseudoGwasEffectVal.resize(numKeptLDBlocks);
    
    VectorXf n_trn(numKeptLDBlocks);
    VectorXf n_val(numKeptLDBlocks);
    pseudoGwasNtrnBlock.resize(numKeptLDBlocks);
    b_val.resize(numIncdSnps);

    for (unsigned i=0; i<numKeptLDBlocks; ++i) {
        LDBlockInfo* block = keptLdBlockInfoVec[i];
        
        n_trn[i] = 0.9*nGWASblock[i];
        n_val[i] = nGWASblock[i] - n_trn[i];
        pseudoGwasNtrnBlock[i] = n_trn[i];

        long size = eigenValLdBlock[i].size();
        VectorXf rnd(size);
        for (unsigned j=0; j<size; ++j) {
            rnd[j] = Stat::snorm();
        }
        
        pseudoGwasEffectTrn[i] = gwasEffectInBlock[i] + sqrt(1.0/n_trn[i] - 1.0/nGWASblock[i]) * eigenVecLdBlock[i] * (eigenValLdBlock[i].array().sqrt().matrix().asDiagonal() * rnd);

        pseudoGwasEffectVal[i] = nGWASblock[i]/n_val[i] * gwasEffectInBlock[i] - n_trn[i]/n_val[i] * pseudoGwasEffectTrn[i];
        b_val.segment(block->startSnpIdx, block->numSnpInBlock) = pseudoGwasEffectVal[i];
    }
}

void Data::readEigenMatrixBinaryFileAndMakeWandQBivariate(const string &dirname, const float eigenCutoff, const vector<VectorXf> &GWASeffects1, const vector<VectorXf> &GWASeffects2, const vector<Vector2f> &nGWASblock, const bool noscale){
    if (!wcorrBlocks.size()) {  // only print for the first time reading the data
        cout << "Reading eigenvectors from binary file and making W and Q matrices for bivariate analysis..." << endl;
    }
    if (!Gadget::directoryExist(dirname)) {
        throw("Error: cannot find the folder [" + dirname + "]");
    }
    
    vector<int>numSnpInRegion(numKeptLDBlocks);
    
    for(int i = 0; i < numKeptLDBlocks;i++){
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }

    unsigned subLdPlannedBiv = 0;
    vector<string> subLdMissingLdBlockIdsBiv;
    for (int i = 0; i < numKeptLDBlocks; ++i) {
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        int m = block->numSnpInBlock;
        if (m <= 1) continue;
        if (!recomputeEigen) continue;
        unsigned n_skip = 0;
        for (int j = 0; j < m; ++j)
            if (block->snpInfoVec[j]->skip) ++n_skip;
        if (n_skip == 0 || n_skip >= (unsigned)m) continue;
        if (blockLdmBinFileExists(dirname, block->ID) || blockEigenBinFileExists(dirname, block->ID))
            ++subLdPlannedBiv;
        else
            subLdMissingLdBlockIdsBiv.push_back(block->ID);
    }
    if (!subLdMissingLdBlockIdsBiv.empty() && !wcorrBlocks.size()) {
        cout << "Warning (bivariate): --recompute-eigen with --skip: partial skip in " << subLdMissingLdBlockIdsBiv.size()
             << " LD block(s); need block<ID>.ldm.bin or block<ID>.eigen.bin in [" << dirname << "]. Missing for block ID(s):";
        for (size_t t = 0; t < subLdMissingLdBlockIdsBiv.size(); ++t)
            cout << " " << subLdMissingLdBlockIdsBiv[t];
        cout << ". Falling back to full-block block*.eigen.bin for those blocks." << endl;
    }
    if (subLdPlannedBiv && !wcorrBlocks.size()) {
        cout << "For " << subLdPlannedBiv << " LD block(s) (bivariate), the LD submatrix of retained SNPs will be reconstructed (--recompute-eigen) and the GWAS effects of skipped SNPs will not be used." << endl;
    }

    eigenValLdBlock.resize(numLDBlocks);
    eigenVecLdBlock.resize(numLDBlocks);
    wcorrBlocks.resize(numKeptLDBlocks);
    numSnpsBlock.resize(numKeptLDBlocks);
    numEigenvalBlock.resize(numKeptLDBlocks);
    Qblocks.resize(numKeptLDBlocks);

#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        block->eigenColRemap.resize(0);
        block->subLdKeptLocalIdx.clear();

        VectorXf eigenValSub;
        MatrixXf eigenVecSub;
        VectorXf wcorrSub;
        MatrixXf Qsub;
        float sumPosSub = 0;
        VectorXi remapSub;
        vector<unsigned> keptSub;

        if (trySubLdEigenFromFullLdmBivariate(block, dirname, eigenCutoff, GWASeffects1[i], GWASeffects2[i], wcorrSub, Qsub, eigenValSub, eigenVecSub, sumPosSub, remapSub, keptSub, nullptr)) {
#pragma omp critical
            {
                cout << "LD block " << block->ID << " (bivariate): " << block->numSnpInBlock - (int)keptSub.size()
                     << " SNPs skipped, " << keptSub.size() << " retained." << endl;
            }
            eigenValLdBlock[i] = eigenValSub;
            eigenVecLdBlock[i] = eigenVecSub;
            wcorrBlocks[i] = wcorrSub;
            Qblocks[i] = Qsub;
            block->sumPosEigVal = sumPosSub;
            block->eigenvalues = eigenValSub;
            block->eigenColRemap = remapSub;
            block->subLdKeptLocalIdx = keptSub;
            numSnpsBlock[i] = Qblocks[i].cols();
            numEigenvalBlock[i] = Qblocks[i].rows();
            eigenVecLdBlock[i].resize(0,0);
            continue;
        }

        int32_t cur_m = 0;
        int32_t cur_k = 0;
        float sumPosEigVal = 0;
        float oldEigenCutoff =0;
        
        string infile = dirname + "/block" + block->ID + ".eigen.bin";
        FILE *fp = fopen(infile.c_str(), "rb");
        if(!fp){cout << "Error: can not open the file [" + infile + "] to read." << endl;}
        if(!fp){throw ("Error: can not open the file [" + infile + "] to read.");}

        // 1. marker number
        if(fread(&cur_m, sizeof(int32_t), 1, fp) != 1){
            throw("Read " + infile + " error (m)");
        }

        if(cur_m != numSnpInRegion[i]){
            throw("In LD block " + block->ID + ", inconsistent marker number to marker information in " + infile);
        }

        // 2. ncol of eigenVec (number of eigenvalues)
        if(fread(&cur_k, sizeof(int32_t), 1, fp) != 1){
            throw("In LD block " + block->ID + ", error about number of eigenvalues in  " + infile);
        }

        // 3. sum of all positive eigenvalues
        if(fread(&sumPosEigVal, sizeof(float), 1, fp) != 1){
            throw("In LD block " + block->ID + ", error about the sum of positive eigenvalues in " + infile);
        }

        // 4. eigenCutoff
        if(fread(&oldEigenCutoff, sizeof(float), 1, fp) != 1){
            throw("In LD block " + block->ID + ", error about eigen cutoff used in " + infile);
        }

        // 5. eigenvalues
        VectorXf lambda(cur_k);
        if(fread(lambda.data(), sizeof(float), cur_k, fp) != cur_k){
            throw("In LD block " + block->ID + ",size error about eigenvalues in " + infile);
        }

        // 6. eigenvector
        MatrixXf U(cur_m, cur_k);
        uint64_t nElements = (uint64_t)cur_m * (uint64_t)cur_k;
        if(fread(U.data(), sizeof(float), nElements, fp) != nElements){
            throw("In LD block " + block->ID + ",size error about eigenvectors in " + infile);
        }
        
        fclose(fp);
        
        if(oldEigenCutoff < eigenCutoff & i == 0){
            cout << "Warning: current proportion of variance in LD block is set as " + to_string(eigenCutoff)+ ". But the proportion of variance is set as "<< to_string(oldEigenCutoff) + " in "  + infile + ".\n";
        }
        
        if (eigenCutoff < oldEigenCutoff) {
            truncateEigenMatrix(sumPosEigVal, eigenCutoff, lambda, U, eigenValLdBlock[i], eigenVecLdBlock[i]);
        } else {
            eigenValLdBlock[i] = lambda;
            eigenVecLdBlock[i] = U;
        }
        block->sumPosEigVal = sumPosEigVal;
        block->eigenvalues = lambda;
                        
        // make w and Q for bivariate: concatenate w1 and w2
        VectorXf sqrtLambda = eigenValLdBlock[i].array().sqrt();
        VectorXf w1 = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (eigenVecLdBlock[i].transpose() * GWASeffects1[i]);
        VectorXf w2 = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (eigenVecLdBlock[i].transpose() * GWASeffects2[i]);
        
        // Concatenate w1 and w2: [w1; w2]
        wcorrBlocks[i].resize(w1.size() + w2.size());
        wcorrBlocks[i].head(w1.size()) = w1;
        wcorrBlocks[i].tail(w2.size()) = w2;
        
        Qblocks[i] = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
        
        numSnpsBlock[i] = Qblocks[i].cols();
        numEigenvalBlock[i] = Qblocks[i].rows();
    }
}

void Data::constructWandQ(const vector<VectorXf> &GWASeffects, const float nGWAS, const bool noscale) {
    wcorrBlocks.resize(numKeptLDBlocks);
    numSnpsBlock.resize(numKeptLDBlocks);
    numEigenvalBlock.resize(numKeptLDBlocks);
    Qblocks.clear();

    for (unsigned i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
        // calculate wbcorr and Qblocks
        VectorXf sqrtLambda = eigenValLdBlock[i].array().sqrt();
        //cout << eigenVecLdBlock[i].transpose().rows() << " " << eigenVecLdBlock[i].transpose().cols() << " " << gwasEffectInBlock[i].size() << endl;
        wcorrBlocks[i] = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (eigenVecLdBlock[i].transpose() * GWASeffects[i] );
        // cout << "eigenVecLdBlock[i]: " << eigenVecLdBlock[i] << endl;
        // cout << "wcorrBlocks[i]: " << wcorrBlocks[i] << endl;
        // cout << "sqrtLambda: " << sqrtLambda << endl;
        // cout << gwasEffectInBlock[i] << endl;
//        MatrixXf tmpQblocks = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
//        MatrixDat matrixDat = MatrixDat(ldblock->snpNameVec, tmpQblocks);
//        // cout << "Qblock: " << endl;
//        // cout << matrixDat.values << endl;
//        Qblocks.push_back(matrixDat);
        Qblocks[i] = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
        
        if (noscale) {
            VectorXf Dsqrt(ldblock->numSnpInBlock);
            for (unsigned j=0; j<ldblock->numSnpInBlock; ++j) {
                SnpInfo *snp = ldblock->snpInfoVec[j];
                //Dsqrt[j] = sqrt(snp->twopq);
                Dsqrt[j] = snp->gwas_scalar;
            }
            Qblocks[i] = Qblocks[i] * Dsqrt.asDiagonal();
        }
        
        numSnpsBlock[i] = Qblocks[i].cols();
        numEigenvalBlock[i] = Qblocks[i].rows();
        
        eigenVecLdBlock[i].resize(0,0);
    }
    
    nGWASblock.resize(numKeptLDBlocks);
    for (unsigned i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
        nGWASblock[i] = nGWAS;
    }
}

void Data::outputWandQ(const string &dirname) {
    // output w and Q in text format
    if (!Gadget::directoryExist(dirname)){
        Gadget::createDirectory(dirname);
        cout << "  Created directory [" << dirname << "] to store w and Q matrices.\n\n";
    }
#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[i];
        
        string outTxtfile1 = dirname + "/block" + blockInfo->ID + ".w.txt";
        string outTxtfile2 = dirname + "/block" + blockInfo->ID + ".Q.txt";
        ofstream outtxt1(outTxtfile1.c_str());
        ofstream outtxt2(outTxtfile2.c_str());

        outtxt1 << wcorrBlocks[i] << endl;
        outtxt2 << Qblocks[i] << endl;
        
        outtxt1.close();
        outtxt2.close();
    }
}

void Data::scaleGwasEffects(){
    snp2pq.resize(numIncdSnps);
    b.resize(numIncdSnps);
    n.resize(numIncdSnps);
    se.resize(numIncdSnps);
    scalar.resize(numIncdSnps);
    tss.resize(numIncdSnps); // only used in SBayesC
    ZPy.resize(numIncdSnps);
    const char *gwasScalarLabel = "estimated by (1/sqrt(n*SE^2+b^2) from summary stats";
    if (gwasScalarMode == GwasScalarMode::Gwas2pq)
        gwasScalarLabel = "sqrt(2pq) based on allele frequency in GWAS)";
    else if (gwasScalarMode == GwasScalarMode::Ref2pq)
        gwasScalarLabel = "(sqrt(2pq) based on allele frequency in LD reference";
    cout << "GWAS effect scaling: " << gwasScalarLabel << "." << endl;

    SnpInfo *snp;
    unsigned ref2pqFallback = 0;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snp = incdSnpInfoVec[i];
        const float refAf = snp->af;
        snp->af = snp->gwas_af;
        snp2pq[i] = snp->twopq = 2.0f*snp->gwas_af*(1.0f-snp->gwas_af);
        if(snp2pq[i]==0) cout << "Error: SNP " << snp->ID << " af " << snp->af << " has 2pq = 0." << endl;
        b[i] = snp->gwas_b;
        n[i] = snp->gwas_n;
        se[i]= snp->gwas_se;
        double sc = 0.0;
        if (gwasScalarMode == GwasScalarMode::LeastSquares) {
            const double denom = double(n[i])*double(se[i])*double(se[i]) + double(b[i])*double(b[i]);
            if (!(denom > 0.0))
                throw("Error: SNP " + snp->ID + " has non-positive n*SE^2+b^2 for --gwas-scalar ls.");
            sc = 1.0/std::sqrt(denom);
        } else if (gwasScalarMode == GwasScalarMode::Ref2pq) {
            if (refAf > 0.0f && refAf < 1.0f)
                sc = std::sqrt(double(2.0f*refAf*(1.0f-refAf)));
            else {
                sc = std::sqrt(double(snp2pq[i]));
                ++ref2pqFallback;
            }
        } else {
            sc = std::sqrt(double(snp2pq[i]));
        }
        snp->gwas_scalar = sc;
        scalar[i] = snp->gwas_scalar;
//        b[i] = snp->gwas_b * sqrt(snp2pq[i]); // scale the marginal effect so that it's in per genotype SD unit
//        n[i] = snp->gwas_n;
//        se[i]= snp->gwas_se * sqrt(snp2pq[i]);
//        tss[i] = n[i]*(n[i]*se[i]*se[i] + b[i]*b[i]);
//        ZPy[i] = n[i]*b[i];
    }
    if (ref2pqFallback)
        cout << "Warning: --gwas-scalar ref2pq used GWAS sqrt(2pq) for " << ref2pqFallback << " SNPs (LD reference AF missing or not in (0,1))." << endl;

    // estimate phenotypic variance
    varySnp = snp2pq.array()*(n.array()*se.array().square()+b.array().square());
    VectorXf varpSrt = varySnp;
    std::sort(varpSrt.data(), varpSrt.data() + varpSrt.size());
    float obsVarPhenotypic = varpSrt[varpSrt.size()/2];
    varPhenotypic = 1.0;
    //cout << "varPhenotypic: " << varPhenotypic << endl;
    VectorXf nSrt = n;
    std::sort(nSrt.data(), nSrt.data() + nSrt.size());
    numKeptInds = nSrt[nSrt.size()/2]; // median
    
    // estimate per-SNP 2pq using the estimated phenotypic variance
    // for (unsigned i=0; i<numIncdSnps; ++i) {
    //     snp = incdSnpInfoVec[i];
    //     snp2pq[i] = snp->twopq = obsVarPhenotypic/(snp->gwas_n*se[i]*se[i]+b[i]*b[i]);       // NEW!
    // }
    
    // scale GWAS effects
    b.array() *= scalar.array();
    se.array() *= scalar.array();
    tss.array() = n.array()*(n.array()*se.array().square() + b.array().square());
    ZPy.array() = n.array()*b.array();
    
    // map to blocks
    if (numKeptLDBlocks) {
        gwasEffectInBlock.resize(numKeptLDBlocks);
        gwasPerSnpNinBlock.resize(numKeptLDBlocks);
        nGWASblock.resize(numKeptLDBlocks);
        for (unsigned i = 0; i < numKeptLDBlocks; i++){
            LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
            gwasEffectInBlock[i] = b(ldblock->block2GwasSnpVec);
            gwasPerSnpNinBlock[i] = n(ldblock->block2GwasSnpVec);
//            nGWASblock[i] = gwasPerSnpNinBlock[i].mean();
            VectorXf nBlockSrt = gwasPerSnpNinBlock[i];
            std::sort(nBlockSrt.data(), nBlockSrt.data() + nBlockSrt.size());
            nGWASblock[i] = nBlockSrt[nBlockSrt.size()/2]; // median
//            nGWASblock[i] = numKeptInds;
        }
    }
    
    // calculate ypy (total sum of squares) for summary statistics
    // ypy = mean phenotypic variance which is one after scaling * number of individuals
    ypy = numKeptInds;


    // output GWAS data
}

void Data::scaleBivariateGwasEffects(){
    // Scale GWAS effects for both traits (--gwas-scalar matches univariate Data::scaleGwasEffects)
    snp2pq.resize(numIncdSnps);
    b.resize(numIncdSnps);
    n.resize(numIncdSnps);
    se.resize(numIncdSnps);
    scalar.resize(numIncdSnps);
    b2.resize(numIncdSnps);
    se2.resize(numIncdSnps);
    scalar2.resize(numIncdSnps);
    n2.resize(numIncdSnps);

    const char *gwasScalarLabel = "ls (1/sqrt(n*SE^2+b^2) from summary stats)";
    if (gwasScalarMode == GwasScalarMode::Gwas2pq)
        gwasScalarLabel = "gwas2pq (sqrt(2pq) from GWAS AF)";
    else if (gwasScalarMode == GwasScalarMode::Ref2pq)
        gwasScalarLabel = "ref2pq (sqrt(2pq) from LD reference AF before GWAS AF overwrite)";
    cout << "GWAS effect scaling (bivariate): " << gwasScalarLabel << "." << endl;
    
    SnpInfo *snp;
    unsigned ref2pqFallback = 0;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snp = incdSnpInfoVec[i];
        const float refAf = snp->af;
        snp->af = snp->gwas_af;  // trait 1 GWAS AF on SnpInfo::af (same convention as univariate)
        snp2pq[i] = snp->twopq = 2.0f*snp->gwas_af*(1.0f-snp->gwas_af);
        if(snp2pq[i]==0) cout << "Error: SNP " << snp->ID << " af " << snp->af << " has 2pq = 0." << endl;
        
        b[i] = snp->gwas_b;
        n[i] = snp->gwas_n;
        se[i] = snp->gwas_se;
        b2[i] = snp->gwas_b2;
        n2[i] = snp->gwas_n2;
        se2[i] = snp->gwas_se2;

        double sc1 = 0.0, sc2 = 0.0;
        if (gwasScalarMode == GwasScalarMode::LeastSquares) {
            const double denom1 = double(n[i])*double(se[i])*double(se[i]) + double(b[i])*double(b[i]);
            if (!(denom1 > 0.0))
                throw("Error: SNP " + snp->ID + " trait 1: non-positive n*SE^2+b^2 for --gwas-scalar ls.");
            sc1 = 1.0/std::sqrt(denom1);
            const double denom2 = double(n2[i])*double(se2[i])*double(se2[i]) + double(b2[i])*double(b2[i]);
            if (!(denom2 > 0.0))
                throw("Error: SNP " + snp->ID + " trait 2: non-positive n*SE^2+b^2 for --gwas-scalar ls.");
            sc2 = 1.0/std::sqrt(denom2);
        } else if (gwasScalarMode == GwasScalarMode::Ref2pq) {
            if (refAf > 0.0f && refAf < 1.0f) {
                const double t = double(2.0f*refAf*(1.0f-refAf));
                sc1 = sc2 = std::sqrt(t);
            } else {
                sc1 = std::sqrt(double(snp2pq[i]));
                const float tpq2 = 2.0f*snp->gwas_af2*(1.0f-snp->gwas_af2);
                if (tpq2 <= 0.f)
                    throw("Error: SNP " + snp->ID + " trait 2 GWAS AF gives 2pq = 0; cannot complete ref2pq fallback.");
                sc2 = std::sqrt(double(tpq2));
                ++ref2pqFallback;
            }
        } else { // Gwas2pq
            sc1 = std::sqrt(double(snp2pq[i]));
            const float tpq2 = 2.0f*snp->gwas_af2*(1.0f-snp->gwas_af2);
            if (tpq2 <= 0.f)
                throw("Error: SNP " + snp->ID + " trait 2 GWAS AF gives 2pq = 0 for --gwas-scalar gwas2pq.");
            sc2 = std::sqrt(double(tpq2));
        }
        snp->gwas_scalar = sc1;
        scalar[i] = snp->gwas_scalar;
        snp->gwas_scalar2 = scalar2[i] = sc2;
    }
    if (ref2pqFallback)
        cout << "Warning: bivariate --gwas-scalar ref2pq used per-trait GWAS sqrt(2pq) for " << ref2pqFallback << " SNPs (LD reference AF missing or not in (0,1))." << endl;

    // estimate phenotypic variance for both traits (before scaling)
    varySnp  = snp2pq.array()*(n.array() *se.array().square() +b.array().square());
    varySnp2 = snp2pq.array()*(n2.array()*se2.array().square()+b2.array().square());
    VectorXf varpSrt = varySnp;
    std::sort(varpSrt.data(), varpSrt.data() + varpSrt.size());
    float obsVarPhenotypic = varpSrt[varpSrt.size()/2];
    varPhenotypic = 1.0;
    VectorXf nSrt = n;
    std::sort(nSrt.data(), nSrt.data() + nSrt.size());
    numKeptInds = nSrt[nSrt.size()/2]; // median sample size trait 1
    VectorXf n2Srt = n2;
    std::sort(n2Srt.data(), n2Srt.data() + n2Srt.size());
    numKeptInds2 = n2Srt[n2Srt.size()/2]; // median sample size trait 2
    
    // Implied 2pq from trait 1 summary stats (legacy bivariate LS path only)
    if (gwasScalarMode == GwasScalarMode::LeastSquares) {
        for (unsigned i=0; i<numIncdSnps; ++i) {
            snp = incdSnpInfoVec[i];
            snp2pq[i] = snp->twopq = obsVarPhenotypic/(snp->gwas_n*se[i]*se[i]+b[i]*b[i]);
        }
    }
    
    // scale GWAS effects for both traits
    b.array() *= scalar.array();
    se.array() *= scalar.array();
    b2.array() *= scalar2.array();
    se2.array() *= scalar2.array();
    
    // calculate ypy (total sum of squares) for summary statistics
    ypy = numKeptInds;
    ypy2 = numKeptInds2;
//    string outfile = "ma.txt";
//    ofstream out(outfile.c_str());
//    for (unsigned i=0; i<numIncdSnps; ++i) {
//        snp = incdSnpInfoVec[i];
//        out << snp->ID << " " << snp->a1 << " " << snp->a2 << " " << snp->af << " " << snp->gwas_b << " " << snp->gwas_se << " " << snp->gwas_pvalue << " " << snp->gwas_n << " " << snp2pq[i]*n[i] << " " << varySnp[i] << " " << b[i] << endl;
//    }
//    out.close();

}


void Data::resizeBlockLDmatrix(const string &inDirname, const string &outLDmatType, const string &includeSnpFile, const string &outDirname, const bool writeLdmTxt){
    readBlockLdmInfoFile(inDirname);
    readBlockLdmSnpInfoFile(inDirname);

    includeSnp(includeSnpFile);
    //includeMatchedSnp();
    //mapSnpsToBlocks();
    //includeMatchedBlocks();
    
    
//    ldBlockInfoVec = keptLdBlockInfoVec;
//    numLDBlocks = numKeptLDBlocks;
    
    struct stat sb;
    if (stat(inDirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        throw("Error: cannot find the folder [" + inDirname + "]");
    }
    
    if (stat(outDirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        string create_cmd = "mkdir " + outDirname;
        system(create_cmd.c_str());
        cout << "Created folder [" << outDirname << "] to store resized LD matrices." << endl;
    }

    vector<int> numSnpInRegion;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        LDBlockInfo *block = ldBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }
        
    keptLdBlockInfoVec = ldBlockInfoVec;

    cout << "Resizing block LD matrices... " << endl;

#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numLDBlocks; i++){
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[i];
                
        int32_t blockSize = blockInfo->numSnpInBlock;
        
        unsigned newBlockSize = 0;
        for (unsigned j=0; j<blockSize; ++j) {
            SnpInfo *snpj = blockInfo->snpInfoVec[j];
            if (snpj->included) {
                ++newBlockSize;
            }
        }
        
        if (!newBlockSize) continue;
        
        
        if (outLDmatType == "block") {
            
            MatrixXf ldm(blockSize, blockSize);
            readBlockLDmatrix(inDirname, blockInfo->ID, blockSize, ldm);

            //        cout << " read block " << blockInfo->ID << endl;
                        
            vector<string> newSnpNameVec;
            vector<SnpInfo*> newSnpInfoVec;
            
            MatrixXf ldm_resized(newBlockSize, newBlockSize);
            
            unsigned jj=0;
            for (unsigned j=0; j<blockSize; ++j) {
                SnpInfo *snpj = blockInfo->snpInfoVec[j];
                if (!snpj->included) continue;
                newSnpNameVec.push_back(snpj->ID);
                newSnpInfoVec.push_back(snpj);
                unsigned kk=0;
                for (unsigned k=0; k<blockSize; ++k) {
                    SnpInfo *snpk = blockInfo->snpInfoVec[k];
                    if (!snpk->included) continue;
                    ldm_resized(jj,kk) = ldm(j,k);
                    ++kk;
                }
                ++jj;
            }
            
            blockInfo->numSnpInBlock = newBlockSize;
            blockInfo->snpNameVec = newSnpNameVec;
            blockInfo->snpInfoVec = newSnpInfoVec;
            
            string outBinfile = outDirname + "/block" + blockInfo->ID + ".ldm.bin";
            FILE *outbin = fopen(outBinfile.c_str(), "wb");
            uint64_t nElements = (uint64_t) newBlockSize * (uint64_t) newBlockSize;
            fwrite(ldm_resized.data(), sizeof(float), nElements, outbin);
            fclose(outbin);
            
            if (writeLdmTxt) {
                string outTxtfile = outDirname + "/block" + blockInfo->ID + ".ldm.txt";
                ofstream outtxt;
                outtxt.open(outTxtfile.c_str());
                for (unsigned ii=0; ii<newBlockSize; ++ii){
                    for (unsigned jj=0; jj<newBlockSize; ++jj) {
                        outtxt << blockInfo->ID << "\t" << blockInfo->snpNameVec[ii] << "\t" << blockInfo->snpNameVec[jj] << "\t" << ldm_resized(ii,jj) << endl;
                    }
                }
                outtxt.close();
            }
            
        }
        else if (outLDmatType == "blockSparse") {
            string infile = inDirname + "/block" + blockInfo->ID + ".sparse.bin";
            FILE *fp = fopen(infile.c_str(), "rb");
            if(!fp){throw ("Error: can not open the file [" + infile + "] to read.");}
            
            SparseMatrix<float> ldmSp = Gadget::readSparseMatrixBinary(infile);
            
            //        cout << " read block " << blockInfo->ID << endl;
            
            vector<string> newSnpNameVec;
            vector<SnpInfo*> newSnpInfoVec;
            
            vector<Triplet<float> > tripletList;
            tripletList.reserve(newBlockSize * newBlockSize);  // preallocate for speed
            
            unsigned jj=0;
            for (unsigned j=0; j<blockSize; ++j) {
                SnpInfo *snpj = blockInfo->snpInfoVec[j];
                if (!snpj->included) continue;
                newSnpNameVec.push_back(snpj->ID);
                newSnpInfoVec.push_back(snpj);
                unsigned kk=0;
                for (unsigned k=0; k<blockSize; ++k) {
                    SnpInfo *snpk = blockInfo->snpInfoVec[k];
                    if (!snpk->included) continue;
                    float val = ldmSp.coeff(j, k);  // safely access sparse matrix element
                    if (val != 0.0f) {
                        tripletList.emplace_back(jj, kk, val);
                    }
                    ++kk;
                }
                ++jj;
            }
            
            SparseMatrix<float> ldmSp_resized(newBlockSize, newBlockSize);
            ldmSp_resized.setFromTriplets(tripletList.begin(), tripletList.end());
            
            blockInfo->numSnpInBlock = newBlockSize;
            blockInfo->snpNameVec = newSnpNameVec;
            blockInfo->snpInfoVec = newSnpInfoVec;
            
            string outBinfile = outDirname + "/block" + blockInfo->ID + ".sparse.bin";
            Gadget::writeSparseMatrixBinary(ldmSp_resized, outBinfile);
            
            if (writeLdmTxt) {
                string outTxtfile = outDirname + "/block" + blockInfo->ID + ".sparse.txt";
                Gadget::writeSparseMatrixToText(ldmSp_resized, outTxtfile);
            }
        }
            
        string outSnpfile = outDirname + "/block" + blockInfo->ID + ".snp.info";
        string outldmfile = outDirname + "/block" + blockInfo->ID + ".ldm.info";
        
        outputBlockLDmatrixInfo(*blockInfo, outSnpfile, outldmfile);
        
        if(!(i%1)) cout << " computed block " << blockInfo->ID << "\r" << flush;
    }
    
    if (outLDmatType == "block") {
        cout << "Written the LD matrix data into file [" << outDirname << "/block*.ldm.bin]." << endl;
        if (writeLdmTxt) cout << "Written the LD matrix data into file [" << outDirname << "/block*.ldm.txt]." << endl;
    } else if (outLDmatType == "blockSparse") {
        cout << "Written the LD matrix data into file [" << outDirname << "/block*.sparse.bin]." << endl;
        if (writeLdmTxt) cout << "Written the LD matrix data into file [" << outDirname << "/block*.sparse.txt]." << endl;
    }
}

void Data::outputBlockLDmatrixTxt(const string &dirname, const unsigned int block){
    readBlockLdmInfoFile(dirname, block);
    readBlockLdmSnpInfoFile(dirname, block);

    struct stat sb;
    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        throw("Error: cannot find the folder [" + dirname + "]");
    }

    vector<int> numSnpInRegion;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        LDBlockInfo *block = ldBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }
            
#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numLDBlocks; i++){
//        if (block && i != block - 1) continue;
        
        LDBlockInfo *blockInfo = ldBlockInfoVec[i];
        
        string infile = dirname + "/block" + blockInfo->ID + ".ldm.bin";
        FILE *fp = fopen(infile.c_str(), "rb");
        if(!fp){throw ("Error: can not open the file [" + infile + "] to read.");}
        
        string outTxtfile = dirname + "/block" + blockInfo->ID + ".ldm.txt";
        ofstream outtxt;
        outtxt.open(outTxtfile.c_str());
        
        int32_t blockSize = blockInfo->numSnpInBlock;
        
        MatrixXf ldm(blockSize, blockSize);
        readBlockLDmatrix(dirname, blockInfo->ID, blockSize, ldm);
        
//        uint64_t nElements = (uint64_t)blockSize * (uint64_t)blockSize;
//        
//        if(fread(ldm.data(), sizeof(float), nElements, fp) != nElements){
//            cout << "fread(U.data(), sizeof(float), nElements, fp): " << fread(ldm.data(), sizeof(float), nElements, fp) << endl;
//            cout << "nEle: " << nElements << " ldm.size: " << ldm.size() <<  " ldm.col: " << ldm.cols() << " row: " << ldm.rows() << endl;
//            throw("In LD block " + blockInfo->ID + ",size error in " + infile);
//            // cout << "Read " << svdLDfile << " error (U)" << endl;
//            // throw("read file error");
//        }
        
//        for (unsigned ii=0; ii<blockSize; ++ii){
//            for (unsigned jj=0; jj<blockSize; ++jj) {
//                outtxt << blockInfo->ID << "\t" << blockInfo->snpNameVec[ii] << "\t" << blockInfo->snpNameVec[jj] << "\t" << ldm(ii,jj) << endl;
//            }
//        }

        //        outtxt << ldm << endl;
        
        for (unsigned row = 0; row < blockSize; ++row) {
            for (unsigned col = 0; col <= row; ++col) {
                outtxt << ldm(row, col) << "\t";
            }
            outtxt << endl;
        }

        outtxt.close();

    }
    
    if (block) {
        cout << "Written the LD matrix data into file [" << dirname << "/block" << block << ".ldm.txt]." << endl;
    } else {
        cout << "Written the LD matrix data into file [" << dirname << "/block*.ldm.txt]." << endl;
    }

}

void Data::resizeBlockLDmatrixAndDoEigenDecomposition(const string &inDirname, const float eigenCutoff, const float rsqThreshold, const string &outDirname, const bool writeLdmTxt){
    //readBlockLdmInfoFile(inDirname + "/ldm.info");
    //readBlockLdmSnpInfoFile(inDirname + "/snp.info");

    struct stat sb;
    if (stat(inDirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        throw("Error: cannot find the folder [" + inDirname + "]");
    }
    
    if (stat(outDirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        string create_cmd = "mkdir " + outDirname;
        system(create_cmd.c_str());
        cout << "Created folder [" << outDirname << "] to store resized LD matrices." << endl;
    }

    vector<int> numSnpInRegion;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        LDBlockInfo *block = ldBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }


    keptLdBlockInfoVec.resize(0);

    cout << "Resizing block LD matrices... " << endl;

#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numLDBlocks; i++){
        
        LDBlockInfo *blockInfo = ldBlockInfoVec[i];
                
        int32_t blockSize = blockInfo->numSnpInBlock;
        
        unsigned newBlockSize = 0;
        for (unsigned j=0; j<blockSize; ++j) {
            SnpInfo *snpj = blockInfo->snpInfoVec[j];
            if (snpj->included) {
                ++newBlockSize;
            }
        }
        
        if (!newBlockSize) {
            blockInfo->kept = false;
            continue;
        } else {
            blockInfo->kept = true;
            keptLdBlockInfoVec.push_back(blockInfo);
        }
        
        MatrixXf ldm(blockSize, blockSize);
        readBlockLDmatrix(inDirname, blockInfo->ID, blockSize, ldm);
        
        //        cout << " read block " << blockInfo->ID << endl;
        
        vector<string> newSnpNameVec;
        vector<SnpInfo*> newSnpInfoVec;
        
        MatrixXf ldm_resized(newBlockSize, newBlockSize);
        
        unsigned jj=0;
        for (unsigned j=0; j<blockSize; ++j) {
            SnpInfo *snpj = blockInfo->snpInfoVec[j];
            if (!snpj->included) continue;
            newSnpNameVec.push_back(snpj->ID);
            newSnpInfoVec.push_back(snpj);
            unsigned kk=0;
            for (unsigned k=0; k<blockSize; ++k) {
                SnpInfo *snpk = blockInfo->snpInfoVec[k];
                if (!snpk->included) continue;
                ldm_resized(jj,kk) = ldm(j,k);
                ++kk;
            }
            ++jj;
        }
        
        blockInfo->numSnpInBlock = newBlockSize;
        blockInfo->snpNameVec = newSnpNameVec;
        blockInfo->snpInfoVec = newSnpInfoVec;
        
        MatrixXf eigenVec;
        VectorXf eigenVal;
        float sumPosEigVal;  // sum of all positive eigenvalues
        // cout << "rval: " << rval << endl;
        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        eigenDecomposition(ldm_resized, eigenCutoff, eigenVal, eigenVec, sumPosEigVal);

        blockInfo->sumPosEigVal = sumPosEigVal;
        
        int32_t numEigenValue = eigenVal.size();
        int32_t numSnpInBlock = newBlockSize;


        
        string outBinfile = outDirname + "/block" + blockInfo->ID + ".eigen.bin";
        FILE *outbin = fopen(outBinfile.c_str(), "wb");
                
        // save summary
        // 1, the number of SNPs in the block
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, outbin);
        // 2, the number of eigenvalues at with the given cutoff
        fwrite(&numEigenValue, sizeof(int32_t), 1, outbin);
        // 3. sum of all the positive eigenvalues
        fwrite(&sumPosEigVal, sizeof(float), 1, outbin);
        // 4. eigenvalue cutoff based on the proportion of variance explained in LD
        fwrite(&eigenCutoff, sizeof(float), 1, outbin);
        // 5. the selected eigenvalues
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, outbin);
        // 6. the selected eigenvector;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, outbin);
        
        fclose(outbin);


        if (writeLdmTxt) {
            string outTxtfile = outDirname + "/block" + blockInfo->ID + ".eigen.txt";
            ofstream outtxt;
            outtxt.open(outTxtfile.c_str());
            outtxt << "Block " << blockInfo->ID << endl;
            outtxt << "numSnps " << numSnpInBlock << endl;
            outtxt << "numEigenvalues " << numEigenValue << endl;
            outtxt << "SumPositiveEigenvalues " << sumPosEigVal << endl;
            outtxt << "EigenCutoff " << eigenCutoff << endl;
            outtxt << "Eigenvalues\n" << eigenVal.transpose() << endl;
            outtxt << "Eigenvectors\n" << eigenVec << endl;
            outtxt << endl;
            outtxt.close();
        }
        
            
        string outSnpfile = outDirname + "/block" + blockInfo->ID + ".snp.info";
        string outldmfile = outDirname + "/block" + blockInfo->ID + ".ldm.info";
        
        outputBlockLDmatrixInfo(*blockInfo, outSnpfile, outldmfile);
        
        outputLDfriends(ldm_resized, blockInfo, outDirname);
        
        if(!(i%1)) cout << " computed block " << blockInfo->ID << "\r" << flush;
    }
    
    
    numKeptLDBlocks = keptLdBlockInfoVec.size();
    
    if (numKeptLDBlocks == 1) {
        cout << "Written the LD matrix data into file [" << outDirname << "/block" << keptLdBlockInfoVec[0]->ID << ".eigen.bin]." << endl;
        if (writeLdmTxt) cout << "Written the LD matrix data into file [" << outDirname << "/block" << keptLdBlockInfoVec[0]->ID << ".eigen.txt]." << endl;
        cout << "Written the high pairwise LD correlations (rsq>0.5) into [" << outDirname << "/block" << keptLdBlockInfoVec[0]->ID << ".rsq0.5.pwld]." << endl;
    } else {
        cout << "Written the LD matrix data into file [" << outDirname << "/block*.eigen.bin]." << endl;
        if (writeLdmTxt) cout << "Written the LD matrix data into file [" << outDirname << "/block*.eigen.txt]." << endl;
        mergeLdmInfo("block", outDirname, true);
    }
}


void Data::readBlockLDmatrix(const string &dirname, const string &blockID, const int32_t blockSize, MatrixXf &ldm) {
    string infile = dirname + "/block" + blockID + ".ldm.bin";

    // Get file size
    struct stat stat_buf;
    if (stat(infile.c_str(), &stat_buf) != 0) {
        throw("Error: cannot stat the file [" + infile + "]");
    }
    size_t file_bytes = stat_buf.st_size;
    size_t n_floats = file_bytes / sizeof(float);

    // Expected sizes
    size_t n_full = (size_t)blockSize * blockSize;
    size_t n_lower_tri = (size_t)blockSize * (blockSize + 1) / 2;

    FILE *fp = fopen(infile.c_str(), "rb");
    if (!fp) {
        throw("Error: cannot open the file [" + infile + "] to read.");
    }

    ldm.setZero(blockSize, blockSize);

    if (n_floats == n_full) {
        // Full square matrix
        if (fread(ldm.data(), sizeof(float), n_full, fp) != n_full) {
            throw("Error reading full matrix from " + infile);
        }
    } else if (n_floats == n_lower_tri) {
        // Lower triangular matrix
        vector<float> buffer(n_lower_tri);
        if (fread(buffer.data(), sizeof(float), n_lower_tri, fp) != n_lower_tri) {
            throw("Error reading lower-triangular matrix from " + infile);
        }

        size_t idx = 0;
        for (int i = 0; i < blockSize; ++i) {
            for (int j = 0; j <= i; ++j) {
                ldm(i, j) = buffer[idx];
                if (i != j) ldm(j, i) = buffer[idx];
                ++idx;
            }
        }
    } else {
        fclose(fp);
        throw("Unrecognized matrix format in " + infile +
              ": got " + to_string(n_floats) + " elements, but expected " +
              to_string(n_full) + " (full) or " + to_string(n_lower_tri) + " (lower-triangular).");
    }

    fclose(fp);
}

void Data::convertToBlockTriangularMatrix(const string &inDirname, const bool writeLdmTxt, const string &outDirname){
    readBlockLdmInfoFile(inDirname);
    readBlockLdmSnpInfoFile(inDirname);

    struct stat sb;
    if (stat(outDirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        string create_cmd = "mkdir " + outDirname;
        system(create_cmd.c_str());
        cout << "Created folder [" << outDirname << "] to store LD matrices." << endl;
    }

#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numLDBlocks; i++){
        LDBlockInfo *blockInfo = ldBlockInfoVec[i];
                
        int32_t blockSize = blockInfo->numSnpInBlock;
        
        MatrixXf ldm(blockSize, blockSize);
        readBlockLDmatrix(inDirname, blockInfo->ID, blockSize, ldm);

        string outBinfile = outDirname + "/block" + blockInfo->ID + ".ldm.bin";
        FILE *outbin = fopen(outBinfile.c_str(), "wb");
        ofstream outtxt;
        string outTxtfile;
        if (writeLdmTxt) {
            outTxtfile = outDirname + "/block" + blockInfo->ID + ".ldm.txt";
            outtxt.open(outTxtfile.c_str());
        }

        // save the lower trangular matrix only
        unsigned numSnpInBlock = blockInfo->numSnpInBlock;
        for (unsigned row = 0; row < numSnpInBlock; ++row) {
            for (unsigned col = 0; col <= row; ++col) {
                float value = ldm(row, col);
                fwrite(&value, sizeof(float), 1, outbin);
            }
        }
        
        if (writeLdmTxt) {
            for (unsigned row = 0; row < numSnpInBlock; ++row) {
                for (unsigned col = 0; col <= row; ++col) {
                    outtxt << ldm(row, col) << "\t";
                }
                outtxt << endl;
            }
        }
        
        fclose(outbin);
        if (writeLdmTxt) outtxt.close();
        
        string outSnpfile = outDirname + "/block" + blockInfo->ID + ".snp.info";
        string outldmfile = outDirname + "/block" + blockInfo->ID + ".ldm.info";

        outputBlockLDmatrixInfo(*blockInfo, outSnpfile, outldmfile);
        
        outputLDfriends(ldm, blockInfo, outDirname);
       
        if(!(i%1)) cout << " computed block " << blockInfo->ID << "\r" << flush;

    }
    
    cout << "Written the LD matrix into folder [" << outDirname << "/block*.ldm.bin]." << endl;
    if (writeLdmTxt) cout << "Written the LD matrix into text file [" << outDirname << "/block*.ldm.txt]." << endl;
    
    mergeLdmInfo("block", outDirname, true);
}


void Data::outputLDfriends(const MatrixXf &ldm, const LDBlockInfo *blockInfo, const string &outDirname){
    // find LD friends for each SNP with LD rsq > threshold
    float rsqThreshold = 0.5;
    int blockSize = ldm.rows();
    vector<vector<int> > ldSnpIdxBlk(blockSize);
    vector<vector<float> > ldcorBlk(blockSize);
    for(unsigned j=0; j < blockSize; j++){
        SnpInfo *snp = blockInfo->snpInfoVec[j];
        for(unsigned k=0; k < j; k++){
            float rsq = ldm(j,k)*ldm(j,k);
            if(rsq > rsqThreshold){
                ldSnpIdxBlk[j].push_back(k);
                ldcorBlk[j].push_back(ldm(j,k));
            }
        }
    }
    
    string pwldfilename = outDirname + "/block" + blockInfo->ID + ".rsq0.5.pwld";  // pairwise LD file
    ofstream out(pwldfilename.c_str());

    out << boost::format("%12s %12s %12s\n")
    % "SNP1"
    % "SNP2"
    % "LDcorrelation";
    
    for (unsigned j=0; j<blockInfo->numSnpInBlock; ++j) {
        unsigned numLDfrd = ldSnpIdxBlk[j].size();
        if (!numLDfrd) continue;
        SnpInfo *snpj = blockInfo->snpInfoVec[j];
        for (unsigned k=0; k<numLDfrd; ++k) {
            SnpInfo *snpk = blockInfo->snpInfoVec[ldSnpIdxBlk[j][k]];
            
            out << boost::format("%12s %12s %12.6f\n")
            % snpj->ID
            % snpk->ID
            % ldcorBlk[j][k];
        }
    }
    out.close();

}
