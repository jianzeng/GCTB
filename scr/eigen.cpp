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
    int revIdx = eigenVal.size();
    VectorXf cumsumNonNeg(revIdx);
    cumsumNonNeg.setZero();
    revIdx = revIdx -1;
    if(eigenVal(revIdx) < 0) cout << "Error, all eigenvector are negative" << endl;
    cumsumNonNeg(revIdx) = eigenVal(revIdx);
    sumPosEigVal = eigenVal(revIdx);
    revIdx = revIdx -1;
    
    while( eigenVal(revIdx) > 1e-10 ){
        sumPosEigVal = sumPosEigVal + eigenVal(revIdx);
        cumsumNonNeg(revIdx) = eigenVal(revIdx) + cumsumNonNeg(revIdx + 1);
        revIdx =revIdx - 1;
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
    unsigned snpIdx = 0;
    int mapped = 0;
    ChromInfo* chr;
    LDBlockInfo *block;
    SnpInfo *snp;
    vector<ChromInfo*>::iterator it = chromInfoVec.begin(), end = chromInfoVec.end();
    unsigned chrCur = (*it)->id;
    int chrEndSnpIdx = (*it)->endSnpIdx;
    
    for (unsigned i=0; i<numLDBlocks; ++i) {
        block = ldBlockInfoVec[i];
        block->snpNameVec.clear();
        block->snpInfoVec.clear();
        
        //cout << i << " " << numLDBlocks << " " << block->chrom << endl;
        
        if (block->chrom > chrCur) {  // move to a new chromosome
            if (it != end) ++it;
            if (it == end) {
                block->kept = false;
                continue;
            }
            chrEndSnpIdx = (*it)->endSnpIdx;
            chrCur = (*it)->id;
        }
        if (block->chrom != chrCur) {
            block->kept = false;
            continue;
        }
        
        //cout << i << " " << numLDBlocks << " " << block->chrom << " " << block->startPos << " " << block->endPos << endl;

        for (; snpIdx <= chrEndSnpIdx; ++snpIdx) {
            //cout << "snpIdx " << snpIdx << " chrEndSnpIdx " << chrEndSnpIdx << endl;
            snp = snpInfoVec[snpIdx];
            //cout << snp->ID << " " << snp->chrom << " " << snp->physPos << " " << block->startPos << " " << block->endPos << endl;
            if (!snp->included) continue;
            if (snp->physPos >= block->startPos && snp->physPos < block->endPos) {
                block->snpNameVec.push_back(snp->ID);
                block->snpInfoVec.push_back(snp);
                snp->block = block->ID;
            } else if (snp->physPos >= block->endPos) {
                break;
            }
        }
        block->numSnpInBlock = block->snpInfoVec.size();
        //cout << "block->numSnpInBlock " << block->numSnpInBlock << endl;
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

        if(!i) cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
            
        MatrixXf rval = generateLDmatrixPerBlock(bedFile, ldblock->snpNameVec);
        
        unsigned numSnpInBlock = ldblock->numSnpInBlock;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numSnpInBlock;
        fwrite(rval.data(), sizeof(float), nElements, outbin);
        
        if (writeLdmTxt) {
            for (unsigned ii=0; ii<numSnpInBlock; ++ii){
                for (unsigned jj=0; jj<numSnpInBlock; ++jj) {
                    outtxt << ldblock->ID << "\t" << ldblock->snpNameVec[ii] << "\t" << ldblock->snpNameVec[jj] << "\t" << rval(ii,jj) << endl;
                }
            }
        }
        
        fclose(outbin);
        if (writeLdmTxt) outtxt.close();
        
        outputBlockLDmatrixInfo(*ldblock, outSnpfile, outldmfile);
        
        if(!(i%1)) cout << " computed block " << ldblock->ID << "\r" << flush;

        if (block) {
            cout << "Written the LD matrix into file [" << outBinfile << "]." << endl;
            if (writeLdmTxt) cout << "Written the LD matrix into file [" << outTxtfile << "]." << endl;
            cout << "Written the LD matrix SNP info into file [" << outSnpfile << "]." << endl;
            cout << "Written the LD matrix ldm info into file [" << outldmfile << "]." << endl;
        }
    }
    
    if (!block) {
        cout << "Written the LD matrix into folder [" << dirname << "/block*.ldm.bin]." << endl;
        if (writeLdmTxt) cout << "Written the LD matrix into text file [" << dirname << "/block*.ldm.txt]." << endl;
        
        if (chromInfoVec.size() >= 22) {  // genome-wide build of LD matrices
            mergeLdmInfo(LDmatType, dirname);
        } else {
            cout << "Written the LD matrix into folder [" << dirname << "/block*.snp.info]." << endl;
            cout << "Written the LD matrix into folder [" << dirname << "/block*.ldm.info]." << endl;
        }
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
    % incdSnpInfoVec[block.startSnpIdx]->ID
    % block.endSnpIdx
    % incdSnpInfoVec[block.endSnpIdx]->ID
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
    VectorXi numTypSnp;
    numImpSnp.setZero(numLDBlocks);
    numTypSnp.setZero(numLDBlocks);
    for (unsigned i = 0; i < numLDBlocks; i++ ){
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (!ldblock->kept) continue;
        for (unsigned j=0; j<ldblock->numSnpInBlock; ++j) {
            SnpInfo *snp = ldblock->snpInfoVec[j];
            if (snp->included) {
                ++numTypSnp[i];
            } else {
                ++numImpSnp[i];
            }
        }
    }
    unsigned totalNumImpSnp = numImpSnp.sum();
    
    if (!totalNumImpSnp) {
        cout << "\nNo SNP needs to be imputed for summary statistics!" << endl;
        return;
    }
    
    cout << "Imputing summary statistics for " << to_string(totalNumImpSnp) << " SNPs in the LD reference but not in the GWAS data file..." << endl;

    Gadget::Timer timer;
    timer.setTime();
    
#pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < numLDBlocks; i++ ){
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (!ldblock->kept) continue;
        
        Stat::Normal normal;

        /// Step 1. construct LD 
        MatrixXf LDPerBlock = eigenVecLdBlock[i] * eigenValLdBlock[i].asDiagonal() * eigenVecLdBlock[i].transpose();

        LDPerBlock.diagonal().array() += (float)diag_mod;
        /// Step 2. Construct the LD correlation matrix among the typed SNPs(LDtt) and the LD correlation matrix among the missing SNPs and typed SNPs (LDit).
        // Step 2.1 divide SNPs into typed and untyped SNPs
        VectorXi typedSnpIdx(numTypSnp[i]);
        VectorXi untypedSnpIdx(numImpSnp[i]);
        VectorXf zTypSnp(numTypSnp[i]);
        VectorXf nTypSnp(numTypSnp[i]);
        VectorXf varyTypSnp(numTypSnp[i]);
        for(unsigned j=0, idxTyp=0, idxImp=0; j < ldblock->numSnpInBlock; j++){
            SnpInfo *snp = ldblock->snpInfoVec[j];
            if(snp->included){
                // typed snp
                typedSnpIdx[idxTyp] = j;
                zTypSnp[idxTyp] = snp->gwas_b / snp->gwas_se;
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
        VectorXf LDi_Z = LDtt.colPivHouseholderQr().solve(zTypSnp);
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
            float base = sqrt(2.0 * snp->af *(1.0 - snp->af) * (nMedian + zImpSnp[j] * zImpSnp[j]));
            snp->gwas_b = zImpSnp[j] * sqrt(varyMedian)/base;
            snp->gwas_se = sqrt(varyMedian) / base;
            snp->gwas_n = nMedian;
            snp->gwas_af = snp->af;
            snp->gwas_pvalue = 2*(1.0-normal.cdf_01(abs(snp->gwas_b/snp->gwas_se)));
            snp->included = true;
            //cout << "b " << snp->gwas_b << " se " << snp->gwas_se << " z " << snp->gwas_b/snp->gwas_se << " p " << snp->gwas_pvalue << endl;
        }
        
        if(!(i%10)) cout << " imputed block " << i << "\r" << flush;
        
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
                % snp->af
                % snp->gwas_b
                % snp->gwas_se
                % snp->gwas_pvalue
                % snp->gwas_n;
            }
            out.close();

            timer.getTime();
            cout << "Imputation of summary statistics is completed (time used: " << timer.format(timer.getElapse()) << ")." << endl;
            cout << "Summary statistics of all SNPs are save into file [" + outfile + "]." << endl;

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
        % snp->af
        % snp->gwas_b
        % snp->gwas_se
        % snp->gwas_pvalue
        % snp->gwas_n;
    }
    out.close();

    timer.getTime();
    cout << "Imputation of summary statistics is completed (time used: " << timer.format(timer.getElapse()) << ")." << endl;
    cout << "Summary statistics of all SNPs are save into file [" + outfile + "]." << endl;
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


void Data::readBlockLdmInfoFile(const string &infoFile){
    // Read bim file: recombination rate is defined between SNP i and SNP i-1
    ifstream in(infoFile.c_str());
    if (!in) throw ("Error: can not open the file [" + infoFile + "] to read.");
    cout << "Reading LDM info from file [" + infoFile + "]." << endl;
    ldBlockInfoVec.clear();
    ldBlockInfoMap.clear();
        
    string header;
    string id;
    int  chr, blockStart, blockEnd, snpNum;
    int idx = 0;
    int snpCount =  1;
    string snpName;
    string startSnpID, endSnpID;
    LDBlockInfo *ldblock;
    getline(in, header);
    while (in >> id >> chr >> blockStart >> startSnpID >> blockEnd >> endSnpID >> snpNum) {
        
        ldblock = new LDBlockInfo(idx++, id, chr);
        ldblock->startSnpIdx = blockStart;
        ldblock->endSnpIdx   = blockEnd;
        ldblock->numSnpInBlock = snpNum;
        
        ldBlockInfoVec.push_back(ldblock);
        ldBlockInfoMap[id] = ldblock;

    }
    in.close();
    numLDBlocks = (unsigned) ldBlockInfoVec.size();
    cout << numLDBlocks << " LD Blocks to be included from [" + infoFile + "]." << endl;
}

void Data::readBlockLdmSnpInfoFile(const string &snpInfoFile){
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
    int idx = 0;
    int index;
    string blockID;
    LDBlockInfo *ldblock;
    getline(in, header);
    while (in >> chr >> id >> index >> genPos >> physPos >> allele1 >> allele2 >> allele1Freq >> ld_n >> blockID) {
        SnpInfo *snp = new SnpInfo(idx++, id, allele1, allele2, chr, genPos, physPos);
        
//        cout << idx << " " << id << " " << blockID << endl;
        
        snp->af = allele1Freq;
        snp->ld_n = ld_n;
        snp->block = blockID;
        snpInfoVec.push_back(snp);
        chromosomes.insert(snp->chrom);
        
        ldblock = ldBlockInfoMap[blockID];
        
        if (ld2snpMap.insert(pair<string, int>(blockID + "_" + id, idx)).second == false) {
            throw ("Error: Duplicate LDBlock-SNP pair found: \"" + blockID + "_" + id + "\".");
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
    cout << numSnps << " SNPs to be included from [" + snpInfoFile + "]." << endl;
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
        if (block && i != block - 1) continue;
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[i];

        string infile = dirname + "/block" + blockInfo->ID + ".ldm.bin";
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
        
        MatrixXf ldm(blockSize, blockSize);
        uint64_t nElements = (uint64_t)blockSize * (uint64_t)blockSize;
                
        if(fread(ldm.data(), sizeof(float), nElements, fp) != nElements){
            cout << "fread(U.data(), sizeof(float), nElements, fp): " << fread(ldm.data(), sizeof(float), nElements, fp) << endl;
            cout << "nEle: " << nElements << " ldm.size: " << ldm.size() <<  " ldm.col: " << ldm.cols() << " row: " << ldm.rows() << endl;
            throw("In LD block " + blockInfo->ID + ",size error in " + outBinfile);
            // cout << "Read " << svdLDfile << " error (U)" << endl;
            // throw("read file error");
        }
        
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

void Data::readEigenMatrixBinaryFile(const string &dirname, const float eigenCutoff){
    if (!Gadget::directoryExist(dirname)) {
        throw("Error: cannot find the folder [" + dirname + "]");
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
        
#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numLDBlocks; i++){
        LDBlockInfo * block;
        block = ldBlockInfoVec[i];
        
        if (!block->kept) continue;
        
        int32_t cur_m = 0;
        int32_t cur_k = 0;
        float sumPosEigVal = 0;
        float oldEigenCutoff =0;
        
        string infile = dirname + "/block" + block->ID + ".eigen.bin";
        FILE *fp = fopen(infile.c_str(), "rb");
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
    }
}

void Data::readEigenMatrixBinaryFileAndMakeWandQ(const string &dirname, const float eigenCutoff, const vector<VectorXf> &GWASeffects, const float nGWAS, const bool makePseudoSummary){
    if (!Gadget::directoryExist(dirname)) {
        throw("Error: cannot find the folder [" + dirname + "]");
    }
    
    vector<int>numSnpInRegion(numKeptLDBlocks);
    
    for(int i = 0; i < numKeptLDBlocks;i++){
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }
    eigenValLdBlock.resize(numLDBlocks);
    eigenVecLdBlock.resize(numLDBlocks);
    wcorrBlocks.resize(numKeptLDBlocks);
    numSnpsBlock.resize(numKeptLDBlocks);
    numEigenvalBlock.resize(numKeptLDBlocks);
    Qblocks.resize(numKeptLDBlocks);

    
    //Constructing pseudo summary statistics for training and validation data sets, with 90% sample size for training and 10% for validation
    
    float n_trn, n_val;
    if (makePseudoSummary) {
        n_trn = 0.9*float(numKeptInds);
        n_val = numKeptInds - n_trn;
        pseudoGwasNtrn = n_trn;

        pseudoGwasEffectTrn.resize(numKeptLDBlocks);
        pseudoGwasEffectVal.resize(numKeptLDBlocks);
        b_val.setZero(numIncdSnps);
    }

#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        int32_t cur_m = 0;
        int32_t cur_k = 0;
        float sumPosEigVal = 0;
        float oldEigenCutoff =0;
        
        string infile = dirname + "/block" + block->ID + ".eigen.bin";
        FILE *fp = fopen(infile.c_str(), "rb");
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
        numSnpsBlock[i] = Qblocks[i].cols();
        numEigenvalBlock[i] = Qblocks[i].rows();
        
        // make pseudo summary data
        if (makePseudoSummary) {
            long size = eigenValLdBlock[i].size();
            VectorXf rnd(size);
            for (unsigned j=0; j<size; ++j) {
                rnd[j] = Stat::snorm();
            }
            
            pseudoGwasEffectTrn[i] = gwasEffectInBlock[i] + sqrt(1.0/n_trn - 1.0/nGWASblock[i]) * eigenVecLdBlock[i] * (eigenValLdBlock[i].array().sqrt().matrix().asDiagonal() * rnd);

            pseudoGwasEffectVal[i] = nGWASblock[i]/n_val * gwasEffectInBlock[i] - n_trn/n_val * pseudoGwasEffectTrn[i];
            b_val.segment(block->startSnpIdx, block->numSnpInBlock) = pseudoGwasEffectVal[i];
        }
        
        eigenVecLdBlock[i].resize(0,0);
    }

    nGWASblock.resize(numKeptLDBlocks);
    for (unsigned i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
        nGWASblock[i] = nGWAS;
    }
}

void Data::truncateEigenMatrix(const float sumPosEigVal, const float eigenCutoff, const VectorXf &oriEigenVal, const MatrixXf &oriEigenVec, VectorXf &newEigenVal, MatrixXf &newEigenVec){
    int revIdx = oriEigenVal.size();
    VectorXf cumsumNonNeg(revIdx);
    cumsumNonNeg.setZero();
    revIdx = revIdx -1;
    if(oriEigenVal(revIdx) < 0) cout << "Error, all eigenvector are negative" << endl;
    cumsumNonNeg(revIdx) = oriEigenVal(revIdx);
    revIdx = revIdx -1;
    
    while(oriEigenVal(revIdx) > 1e-10 ){
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

void Data::readBlockLDmatrixAndDoEigenDecomposition(const string &dirname, const unsigned block, const float eigenCutoff, const bool writeLdmTxt){
    readBlockLdmInfoFile(dirname + "/ldm.info");
    readBlockLdmSnpInfoFile(dirname + "/snp.info");
    readBlockLdmBinaryAndDoEigenDecomposition(dirname, block, eigenCutoff, writeLdmTxt);
}

void Data::readEigenMatrix(const string &dirname, const float eigenCutoff){
    cout << "Reading LD matrix eigen-decomposition data..." << endl;
    //Gadget::Timer timer;
    //timer.setTime();
    
    readBlockLdmInfoFile(dirname + "/ldm.info");
    readBlockLdmSnpInfoFile(dirname + "/snp.info");
    //readEigenMatrixBinaryFile(dirname, eigenCutoff);
    
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
    
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        if (snp->gwas_b == -999) {
            throw("Error: SNP " + snp->ID + " in the LD reference has no summary data. Run --impute-summary first.");
        }
    }
    
    scaleGwasEffects();
    readEigenMatrixBinaryFileAndMakeWandQ(dirname, eigenCutoff, gwasEffectInBlock, numKeptInds, true);
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
    cout << boost::format("%40s %8.3f %8.3f\n") %"LD block size" %Gadget::calcMean(numSnpsBlock) %sqrt(Gadget::calcVariance(numSnpsBlock));
    cout << boost::format("%40s %8.3f %8.3f\n") %"LD block rank" %Gadget::calcMean(numEigenvalBlock) %sqrt(Gadget::calcVariance(numEigenvalBlock));

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


void Data::mergeLdmInfo(const string &outLDmatType, const string &dirname) {
    
    if (outLDmatType != "block") throw("Error: --merge-ldm-info only works for block LD matrices at the moment!");
    
    string dir_path = dirname; // Replace with your folder path
    string search_str = outLDmatType;
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
    
    if (blockIdxSet.size() == 0) {
        throw ("Error: there is no info file to merge in folder [" + dirname + "].");
    }
    
    unsigned nldm = blockIdxSet.size();
    
    set<unsigned>::iterator it = blockIdxSet.begin();
    
    string outSnpInfoFile = dirname + "/snp.info";
    string outldmInfoFile = dirname + "/ldm.info";
    
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
    
    
    
    unsigned snpIdx = 0;
    unsigned ldmIdx = 0;
    
    for (unsigned i=0; i<nldm; ++i) {
        
        string snpInfoFile = dirname + "/block" + to_string(*it) + ".snp.info";
        string ldmInfoFile = dirname + "/block" + to_string(*it) + ".ldm.info";
        
        // read snp info file
        ifstream in1(snpInfoFile.c_str());
        if (!in1) throw ("Error: can not open the file [" + snpInfoFile + "] to read.");
        cout << "Reading SNP info from file [" + snpInfoFile + "]." << endl;
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
        cout << "Reading LDM info from file [" + ldmInfoFile + "]." << endl;
        
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
        
        ++it;
        
        remove(snpInfoFile.c_str());
        remove(ldmInfoFile.c_str());
    }
    
    out1.close();
    out2.close();
    
    cout << "Written " << snpIdx << " SNPs info into file [" + outSnpInfoFile + "]." << endl;
    cout << "Written " << ldmIdx << " LDMs info into file [" + outldmInfoFile + "]." << endl;
    
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
    
    float n_trn = 0.9*float(numKeptInds);
    float n_val = numKeptInds - n_trn;
    pseudoGwasNtrn = n_trn;
    b_val.resize(numIncdSnps);

    for (unsigned i=0; i<numKeptLDBlocks; ++i) {
        LDBlockInfo* block = keptLdBlockInfoVec[i];
        
        long size = eigenValLdBlock[i].size();
        VectorXf rnd(size);
        for (unsigned j=0; j<size; ++j) {
            rnd[j] = Stat::snorm();
        }
        
        pseudoGwasEffectTrn[i] = gwasEffectInBlock[i] + sqrt(1.0/n_trn - 1.0/nGWASblock[i]) * eigenVecLdBlock[i] * (eigenValLdBlock[i].array().sqrt().matrix().asDiagonal() * rnd);

        pseudoGwasEffectVal[i] = nGWASblock[i]/n_val * gwasEffectInBlock[i] - n_trn/n_val * pseudoGwasEffectTrn[i];
        b_val.segment(block->startSnpIdx, block->numSnpInBlock) = pseudoGwasEffectVal[i];
    }
    
    
}

void Data::constructWandQ(const vector<VectorXf> &GWASeffects, const float nGWAS) {
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

void Data::scaleGwasEffects(){
    snp2pq.resize(numIncdSnps);
    b.resize(numIncdSnps);
    n.resize(numIncdSnps);
    se.resize(numIncdSnps);
    tss.resize(numIncdSnps); // only used in SBayesC
    ZPy.resize(numIncdSnps);
    SnpInfo *snp;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snp = incdSnpInfoVec[i];
        snp->af = snp->gwas_af;
        snp2pq[i] = snp->twopq = 2.0f*snp->gwas_af*(1.0f-snp->gwas_af);
        if(snp2pq[i]==0) cout << "Error: SNP " << snp->ID << " af " << snp->af << " has 2pq = 0." << endl;
        b[i] = snp->gwas_b * sqrt(snp2pq[i]); // scale the marginal effect so that it's in per genotype SD unit
        n[i] = snp->gwas_n;
        se[i]= snp->gwas_se * sqrt(snp2pq[i]);
        tss[i] = n[i]*(n[i]*se[i]*se[i] + b[i]*b[i]);
        ZPy[i] = n[i]*b[i];
    }

    // estimate phenotypic variance
    varySnp = (n.array()*(n.array()*se.array().square()+b.array().square()))/n.array();
    VectorXf varpSrt = varySnp;
    std::sort(varpSrt.data(), varpSrt.data() + varpSrt.size());
    varPhenotypic = varpSrt[varpSrt.size()/2];
    //cout << "varPhenotypic: " << varPhenotypic << endl;
    VectorXf nSrt = n;
    std::sort(nSrt.data(), nSrt.data() + nSrt.size());
    numKeptInds = nSrt[nSrt.size()/2]; // median
    
    // map to blocks
    gwasEffectInBlock.resize(numKeptLDBlocks);
    nGWASblock.resize(numKeptLDBlocks);
    for (unsigned i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
        gwasEffectInBlock[i] = b(ldblock->block2GwasSnpVec);
        nGWASblock[i] = numKeptInds;
    }
}


