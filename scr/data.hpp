//
//  data.hpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#ifndef data_hpp
#define data_hpp

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <set>
#include <bitset>
#include <iomanip>     
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <boost/format.hpp>
#include <mpi.h>
#include <omp.h>
#include "gadgets.hpp"
#include "mympi.hpp"

using namespace std;
using namespace Eigen;

class AnnoInfo;

class SnpInfo {
public:
    const string ID;
    const string a1; // the referece allele
    const string a2; // the coded allele
    const int chrom;
    float genPos;
    const int physPos;
    
    int index;
    int window;
    int windStart;  // for window surrounding the SNP
    int windSize;   // for window surrounding the SNP
    int windEnd;
    float af;       // allele frequency
    float twopq;
    bool included;  // flag for inclusion in panel
    bool isQTL;     // for simulation
    bool recoded;   // swap A1 and A2: use A2 as the reference allele and A1 as the coded allele
    bool skeleton;  // skeleton snp for sbayes
    long sampleSize;
    
    VectorXf genotypes; // temporary storage of genotypes of individuals used for building sparse Z'Z
    
    vector<AnnoInfo*> annoPtr;
    vector<unsigned> annoIdx;
    
    float effect;   // estimated effect
    
    // GWAS summary statistics
    float gwas_b;
    float gwas_se;
    float gwas_n;
    float gwas_af;

    float ldSamplVar;    // sum of sampling variance of LD with other SNPs for summary-bayes method
    float ldSum;         // sum of LD with other SNPs
    float ldsc;          // LD score: sum of r^2
    
    int numNonZeroLD;   // may be different from windSize in shrunk ldm
    unsigned numAnnos;

    SnpInfo(const int idx, const string &id, const string &allele1, const string &allele2,
            const int chr, const float gpos, const int ppos)
    : ID(id), index(idx), a1(allele1), a2(allele2), chrom(chr), genPos(gpos), physPos(ppos) {
        window = 0;
        windStart = -1;
        windSize  = 0;
        windEnd   = -1;
        af = -1;
        twopq = -1;
        included = true;
        isQTL = false;
        recoded = false;
        skeleton = false;
        sampleSize = 0;
        effect = 0;
        gwas_b  = -999;
        gwas_se = -999;
        gwas_n  = -999;
        gwas_af = -999;
        ldSamplVar = 0.0;
        ldSum = 0.0;
        ldsc = 0.0;
        numNonZeroLD = 0;
        numAnnos = 0;
    };
    
    void resetWindow(void) {windStart = -1; windSize = 0;};
    bool isProximal(const SnpInfo &snp2, const float genWindow) const;
    bool isProximal(const SnpInfo &snp2, const unsigned physWindow) const;
};

class ChromInfo {
public:
    const int id;
    const unsigned size;
    const int startSnpIdx;
    const int endSnpIdx;
    
    ChromInfo(const int id, const unsigned size, const int startSnp, const int endSnp): id(id), size(size), startSnpIdx(startSnp), endSnpIdx(endSnp){}
};

class AnnoInfo {  // annotation info for SNPs
public:
    int idx;
    const string label;
    unsigned size;
    float fraction;   // fraction of all SNPs in this annotation
    
    unsigned chrom;   // for continuous annotation
    unsigned startBP; // for continuous annotation
    unsigned endBP;   // for continuous annotation
    
    vector<SnpInfo*> memberSnpVec;
    VectorXf snp2pq;
    
    AnnoInfo(const int idx, const string &lab): idx(idx), label(lab){
        size = 0;
        chrom = 0;
        startBP = 0;
        endBP = 0;
    }
    
    void getSnpInfo(void);
    void print(void);
};


class IndInfo {
public:
    const string famID;
    const string indID;
    const string catID;    // catenated family and individual ID
    const string fatherID;
    const string motherID;
    const int famFileOrder; // original fam file order
    const int sex;  // 1: male, 2: female
    
    int index;
    bool kept;
    
    float phenotype;
    
    VectorXf covariates;  // covariates for fixed effects
    
    IndInfo(const int idx, const string &fid, const string &pid, const string &dad, const string &mom, const int sex)
    : famID(fid), indID(pid), catID(fid+":"+pid), fatherID(dad), motherID(mom), index(idx), famFileOrder(idx), sex(sex) {
        phenotype = -9;
        kept = true;
    }
};

class Data {
public:
    MatrixXf X;              // coefficient matrix for fixed effects
    MatrixXf Z;              // coefficient matrix for SNP effects
    VectorXf D;              // 2pqn
    VectorXf y;              // phenotypes
    
    //SparseMatrix<float> ZPZ; // sparse Z'Z because LE is assumed for distant SNPs
    vector<VectorXf> ZPZ;
    vector<SparseVector<float> > ZPZsp;
    SparseMatrix<float> ZPZinv;

    MatrixXf XPX;            // X'X the MME lhs
    MatrixXf ZPX;            // Z'X the covariance matrix of SNPs and fixed effects
    VectorXf XPXdiag;        // X'X diagonal
    VectorXf ZPZdiag;        // Z'Z diagonal
    VectorXf XPy;            // X'y the MME rhs for fixed effects
    VectorXf ZPy;            // Z'y the MME rhs for snp effects
    
    VectorXf snp2pq;         // 2pq of SNPs
    VectorXf se;             // se from GWAS summary data
    VectorXf tss;            // total ss (ypy) for every SNP
    VectorXf b;              // beta from GWAS summary data
    VectorXf n;              // sample size for each SNP in GWAS
    VectorXf Dratio;         // GWAS ZPZdiag over reference ZPZdiag for each SNP
    VectorXf DratioSqrt;     // square root of GWAS ZPZdiag over reference ZPZdiag for each SNP
    VectorXf chisq;          // GWAS chi square statistics = D*b^2
    
    VectorXi windStart;      // leading snp position for each window
    VectorXi windSize;       // number of snps in each window
    
    VectorXf LDsamplVar;     // sum of sampling variance of LD for each SNP with all other SNPs; this is for summary-bayes methods
    VectorXf LDscore;        // sum of r^2 over SNPs in significant LD
    
    float ypy;               // y'y the total sum of squares adjusted for the mean
    float varGenotypic;
    float varResidual;
    
    bool reindexed;
    bool sparseLDM;
    bool shrunkLDM;
    bool readLDscore;
    
    vector<SnpInfo*> snpInfoVec;
    vector<IndInfo*> indInfoVec;

    vector<AnnoInfo*> annoInfoVec;
    vector<string> annoNames;
    
    map<string, SnpInfo*> snpInfoMap;
    map<string, IndInfo*> indInfoMap;


    vector<SnpInfo*> incdSnpInfoVec;
    vector<IndInfo*> keptIndInfoVec;
    
    vector<string> fixedEffectNames;
    vector<string> snpEffectNames;
    
    set<int> chromosomes;
    vector<ChromInfo*> chromInfoVec;
    
    vector<bool> fullSnpFlag;
    
    vector<unsigned> numSnpMldVec;
    vector<unsigned> numSnpAnnoVec;
    
    vector<SparseMatrix<float> > annowiseZPZsp;
    vector<VectorXf> annowiseZPZdiag;
    
    unsigned numFixedEffects;
    unsigned numSnps;
    unsigned numInds;
    unsigned numIncdSnps;
    unsigned numKeptInds;
    unsigned numChroms;
    unsigned numSkeletonSnps;
    unsigned numAnnos;
    
    string label;
    
    Data(){
        numFixedEffects = 0;
        numSnps = 0;
        numInds = 0;
        numIncdSnps = 0;
        numKeptInds = 0;
        numChroms = 0;
        numSkeletonSnps = 0;
        numAnnos = 0;
        
        reindexed = false;
        sparseLDM = false;
        readLDscore = false;
    }
    
    void readFamFile(const string &famFile);
    void readBimFile(const string &bimFile);
    void readBedFile(const string &bedFile);
    void readPhenotypeFile(const string &phenFile, const unsigned mphen);
    void readCovariateFile(const string &covarFile);
    void readGwasSummaryFile(const string &gwasFile, const float afDiff, const float mafmin, const float mafmax);
    void readLDmatrixInfoFile(const string &ldmatrixFile);
    void readLDmatrixBinFile(const string &ldmatrixFile);
    void readGeneticMapFile(const string &freqFile);
    void readfreqFile(const string &geneticMapFile);
    void keepMatchedInd(const string &keepIndFile, const unsigned keepIndMax);
    void includeSnp(const string &includeSnpFile);
    void excludeSnp(const string &excludeSnpFile);
    void includeChr(const unsigned chr);
    void excludeMHC(void);
    void excludeAmbiguousSNP(void);
    void excludeSNPwithMaf(const float mafmin, const float mafmax);
    void excludeRegion(const string &excludeRegionFile);
    void includeSkeletonSnp(const string &skeletonSnpFile);

    void includeMatchedSnp(void);
    vector<SnpInfo*> makeIncdSnpInfoVec(const vector<SnpInfo*> &snpInfoVec);
    vector<IndInfo*> makeKeptIndInfoVec(const vector<IndInfo*> &indInfoVec);
    void getWindowInfo(const vector<SnpInfo*> &incdSnpInfoVec, const unsigned windowWidth, VectorXi &windStart, VectorXi &windSize);
    void getNonoverlapWindowInfo(const unsigned windowWidth);
    void buildSparseMME(const string &bedFile, const unsigned windowWidth);
//    void makeLDmatrix(const string &bedFile, const unsigned windowWidth, const string &filename);
    void makeLDmatrix(const string &bedFile, const string &LDmatType, const float chisqThreshold, const float LDthreshold, const unsigned windowWidth,
                      const string &snpRange, const string &filename, const bool writeLdmTxt);
    void makeshrunkLDmatrix(const string &bedFile, const string &LDmatType, const string &snpRange, const string &filename, const bool writeLdmTxt, const float effpopNE, const float cutOff);
    void resizeWindow(const vector<SnpInfo*> &incdSnpInfoVec, const VectorXi &windStartOri, const VectorXi &windSizeOri,
                      VectorXi &windStartNew, VectorXi &windSizeNew);
    void computeAlleleFreq(const MatrixXf &Z, vector<SnpInfo*> &incdSnpInfoVec, VectorXf &snp2pq);
    void reindexSnp(vector<SnpInfo*> snpInfoVec);
    void initVariances(const float heritability);
    
    void outputSnpResults(const VectorXf &posteriorMean, const VectorXf &posteriorSqrMean, const VectorXf &pip, const string &filename) const;
    void outputFixedEffects(const MatrixXf &fixedEffects, const string &filename) const;
    void outputWindowResults(const VectorXf &posteriorMean, const string &filename) const;
    void summarizeSnpResults(const SparseMatrix<float> &snpEffects, const string &filename) const;
    void buildSparseMME(const bool sampleOverlap);
    void readMultiLDmatInfoFile(const string &mldmatFile);
    void readMultiLDmatBinFile(const string &mldmatFile);
    void outputSnpEffectSamples(const SparseMatrix<float> &snpEffects, const unsigned burnin, const unsigned outputFreq, const string &snpResFile, const string &filename) const;
    void resizeLDmatrix(const string &LDmatType, const float chisqThreshold, const unsigned windowWidth, const float LDthreshold, const float effpopNE, const float cutOff);
    void outputLDmatrix(const string &LDmatType, const string &filename, const bool writeLdmTxt) const;
    void displayAverageWindowSize(const VectorXi &windSize);
    
    void inputSnpResults(const string &snpResFile);
    void inputSnpInfoAndResults(const string &snpResFile);
    void readLDmatrixBinFileAndShrink(const string &ldmatrixFile);
    void readMultiLDmatBinFileAndShrink(const string &mldmatFile);
    void directPruneLDmatrix(const string &ldmatrixFile, const string &outLDmatType, const float chisqThreshold, const string &title, const bool writeLdmTxt);
    void jackknifeLDmatrix(const string &ldmatrixFile, const string &outLDmatType, const string &title, const bool writeLdmTxt);
    void addLDmatrixInfo(const string &ldmatrixFile);
    
    void readAnnotationFile(const string &annotationFile, const bool transpose, const bool allowMultiAnno);
    void readAnnotationFileFormat2(const string &continuousAnnoFile, const unsigned flank); // for continuous annotations
    void setAnnoInfoVec(void);
    void readLDscoreFile(const string &ldscFile);
    void makeAnnowiseSparseLDM(const vector<SparseVector<float> > &ZPZsp, const vector<AnnoInfo *> &annoInfoVec, const vector<SnpInfo*> &snpInfoVec);    
};

#endif /* data_hpp */
