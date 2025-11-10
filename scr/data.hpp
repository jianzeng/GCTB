//
//  data.hpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#ifndef data_hpp
#define data_hpp
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

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
#include <omp.h>
#include <cstdio>
#include "gadgets.hpp"
#include "stat.hpp"

using namespace std;
using namespace Eigen;

typedef SparseMatrix<float, Eigen::ColMajor, long long> SpMat;

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
    int windStartOri;  // original value from .info file
    int windSizeOri;
    int windEndOri;
    float af;       // allele frequency
    float twopq;
    bool included;  // flag for inclusion in panel
    bool isQTL;     // for simulation
    bool iseQTL;
    bool recoded;   // swap A1 and A2: use A2 as the reference allele and A1 as the coded allele
    bool skeleton;  // skeleton snp for sbayes
    bool flipped;   // A1 A2 alleles are flipped in between gwas and LD ref samples
    bool unconverged;  // fail to converge
    bool inCS;  // in a credible set
    bool skip;  // skip sampling its effect
    long sampleSize;
    
    string block;
    int blockIdx;
    
    VectorXf genotypes; // temporary storage of genotypes of individuals used for building sparse Z'Z
    
    vector<AnnoInfo*> annoVec;
    vector<unsigned> annoIdx;   // the index of SNP in the annotation
    map<int, AnnoInfo*> annoMap;
    
    VectorXf annoValues;
    
    float effect;   // estimated effect
    float pip;
    float varExplained;
    float GelmanRubinR;  // Gelman–Rubin's R statistic for convergence diagnostic
    
    // GWAS summary statistics
    double gwas_b;
    double gwas_se;
    float gwas_n;
    float gwas_af;
    double gwas_pvalue;
    double gwas_scalar;

    float ldSamplVar;    // sum of sampling variance of LD with other SNPs for summary-bayes method
    float ldSum;         // sum of LD with other SNPs
    float ldsc;          // LD score: sum of r^2
    
    int ld_n;  // LD reference sample size
    
    int numNonZeroLD;   // may be different from windSize in shrunk ldm
    unsigned numAnnos;

    SnpInfo(const int idx, const string &id, const string &allele1, const string &allele2,
            const int chr, const float gpos, const int ppos)
    : ID(id), index(idx), a1(allele1), a2(allele2), chrom(chr), genPos(gpos), physPos(ppos) {
        window = 0;
        windStart = -1;
        windSize  = 0;
        windEnd   = -1;
        windStartOri = -1;
        windSizeOri = 0;
        windEndOri = -1;
        af = -1;
        twopq = -1;
        included = true;
        isQTL = false;
        iseQTL = false;
        recoded = false;
        skeleton = false;
        flipped = false;
        unconverged = false;
        inCS = false;
        skip = false;
        sampleSize = 0;
        effect = 0;
        varExplained = 0;
        gwas_b  = -999;
        gwas_se = -999;
        gwas_n  = -999;
        gwas_af = -1;
        gwas_pvalue = 1.0;
        gwas_scalar = 0.0;
        ldSamplVar = 0.0;
        ldSum = 0.0;
        ldsc = 0.0;
        numNonZeroLD = 0;
        numAnnos = 0;
        ld_n = -999;
        block = "NA";
        blockIdx = -1;
        GelmanRubinR = -1;
    };
    
    void resetWindow(void) {windStart = -1; windSize = 0;};
    bool isProximal(const SnpInfo &snp2, const float genWindow) const;
    bool isProximal(const SnpInfo &snp2, const unsigned physWindow) const;
};

class LDBlockInfo {
public:
    const string ID;
    const int chrom;
    int index;
    // block
    int startPos;
    int endPos;
    float start_cm;
    float stop_cm;
    int pdist;
    float gdist;
    bool kept;

    // the following variables aim to read svd-ld matrix from SBayesRC-Eigen
    int startSnpIdx;
    int endSnpIdx;
    
    //int idxStart;
    //int idxEnd;
    int preBlock;
    int postBlock;
    //
    vector<string> snpNameVec;
    vector<SnpInfo*> snpInfoVec;
    vector<int> block2GwasSnpVec; // store snps that belong to this block;
    int numSnpInBlock;
    
    VectorXf eigenvalues;
    float sumPosEigVal; // sum of all positive eigenvalues

    LDBlockInfo(const int idx, const string id, const int chr) : index(idx), ID(id), chrom(chr)
    {
        // block info
        startPos = -999;
        endPos = -999;
        start_cm = -999;
        stop_cm = -999;
        pdist = -999;
        gdist = -999;
        // svd'ed ld info
        startSnpIdx = -999;
        endSnpIdx = -999;
        //idxStart = -999;
        //idxEnd = -999;
        preBlock = -999;
        postBlock = -999;
        numSnpInBlock = 0;
        kept = true;
        sumPosEigVal = 0;
    }
};

class locus_bp {
public:
    string locusName;
    int chr;
    int bp;

    locus_bp(string locusNameBuf, int chrBuf, int bpBuf)
    {
        locusName = locusNameBuf;
        chr = chrBuf;
        bp = bpBuf;
    }

    bool operator()(const locus_bp &other)
    {
        return (chr == other.chr && bp <= other.bp);
    }
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
    bool isBinary;
    
    float fraction;   // fraction of all SNPs in this annotation
    float mean;
    float sd;         // standard deviation
    float sum;
    float ssq;        // sum of squares
    
    unsigned chrom;   // for continuous annotation
    unsigned startBP; // for continuous annotation
    unsigned endBP;   // for continuous annotation
    
    vector<SnpInfo*> memberSnpVec;
    map<int, SnpInfo*> memberSnpMap;
    VectorXf snp2pq;
    
    AnnoInfo(const int idx, const string &lab): idx(idx), label(lab){
        size = 0;
        chrom = 0;
        startBP = 0;
        endBP = 0;
        mean = 0.0;
        sd = 0.0;
        sum = 0.0;
        ssq = 0.0;
        isBinary = true;
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
    float rinverse;
    
    VectorXf covariates;  // covariates for fixed effects
    VectorXf randomCovariates;
    
    IndInfo(const int idx, const string &fid, const string &pid, const string &dad, const string &mom, const int sex)
    : famID(fid), indID(pid), catID(fid+":"+pid), fatherID(dad), motherID(mom), index(idx), famFileOrder(idx), sex(sex) {
        phenotype = -9;
        rinverse = 1;
        kept = true;
    }
};

class WindowInfo {
public:
    int index;
    int size;
    int start;
    int end;
    
    double propGenVar;
    double genVarEnrich;
    double genVarEnrichPP;
    
    VectorXf propGenVarMcmc;
    
    vector<SnpInfo*> snpVec;
    
    WindowInfo(const int index, const vector<SnpInfo*> &snpVec): index(index), snpVec(snpVec) {
        size = snpVec.size();
        start = snpVec[0]->index;
        end = snpVec[size-1]->index;
        propGenVar = 0.0;
        genVarEnrich = 0.0;
        genVarEnrichPP = 0.0;
    }
    
    void calcVarEnrichPP(float numWindows);
};

class CredibleSetInfo {
public:
    int index;
    int size;
    float threshold;
    double sumPIP;
    double propVar;
        
    VectorXf PVEmcmc;

    vector<SnpInfo*> snpVec;
    
    int windSize;
    double windPropGenVar;
    double windGenVarEnrich;
    double windGenVarEnrichPP;
    
    int numUnconvgSNPs;  // -1: not tested
    
    CredibleSetInfo(const int index, const float threshold, const float sumPIP, const float propVar, const vector<SnpInfo*> &snpVec)
    : index(index), threshold(threshold), sumPIP(sumPIP), propVar(propVar), snpVec(snpVec){
        size = snpVec.size();
        windSize = 0;
        windPropGenVar = 0.0;
        windGenVarEnrich = 0.0;
        windGenVarEnrichPP = 0.0;
        numUnconvgSNPs = -1;
    }
    
    void getNumUnconvgSNPs(const float threshold);
};

class GeneInfo {
public:
    string ensgid;
    string name;
    string type;
    string genomeBuild;
    
    int chrom;
    int start;
    int end;
        
    GeneInfo(const string &ensgid): ensgid(ensgid){
        chrom = 0;
        start = 0;
        end = 0;
    }
        
    void setFlankingWindow(const int flank);
    bool containAllSnps(vector<SnpInfo*> &snpVec);
};

struct MatrixDat
{
public:
    const vector<string> colnames;
    std::map<string, int> colname2index;
    vector<string> rownames;
    std::map<string, int> rowname2index;
    unsigned ncol;
    unsigned nrow;
    Eigen::MatrixXf values;

    MatrixDat(const vector<string> &colnames, const Eigen::MatrixXf &values)
        : colnames(colnames), ncol(int(colnames.size())), values(values)
    {
        nrow = values.rows();
        rownames.resize(nrow);
        for (unsigned j = 0; j < ncol; j++)
            colname2index.insert(pair<string, int>(colnames[j], j));
        for (unsigned j = 0; j < nrow; j++)
        {
            rownames[j] = "row" + to_string(j);
            rowname2index.insert(pair<string, int>(rownames[j], j));
        }
    }
    MatrixDat(vector<string> &rownames, const vector<string> &colnames, const Eigen::MatrixXf &values)
        : colnames(colnames), ncol(int(colnames.size())), rownames(rownames), nrow(int(rownames.size())), values(values)
    {
        for (unsigned j = 0; j < ncol; j++)
            colname2index.insert(pair<string, int>(colnames[j], j));
        for (unsigned j = 0; j < nrow; j++)
            rowname2index.insert(pair<string, int>(rownames[j], j));
    }
    Eigen::VectorXf col(string nameIdx) const { return values.col(colname2index.at(nameIdx)); }
    Eigen::VectorXf row(string nameIdx) const { return values.row(rowname2index.at(nameIdx)); }
};



class Data {
public:
    MatrixXf X;              // coefficient matrix for fixed effects
    MatrixXf W;              // coefficient matrix for random effects
    MatrixXf Z;              // coefficient matrix for SNP effects
    VectorXf D;              // 2pqn
    VectorXf y;              // phenotypes
    
    //SpMat ZPZ; // sparse Z'Z because LE is assumed for distant SNPs
    vector<VectorXf> ZPZ;
    MatrixXf ZPZmat;
    vector<SparseVector<float> > ZPZsp;
    SpMat ZPZspmat;
    SpMat ZPZinv;
    
    MatrixXf annoMat;        // annotation coefficient matrix
    MatrixXf APA;            // annotation X'X matrix

    MatrixXf XPX;            // X'X the MME lhs
    MatrixXf WPW;
    MatrixXf ZPX;            // Z'X the covariance matrix of SNPs and fixed effects
    VectorXf XPXdiag;        // X'X diagonal
    VectorXf WPWdiag;
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
    VectorXf varySnp;        // per-SNP phenotypic variance
    VectorXf scalar;         // scaling factor for GWAS b
    
    VectorXi windStart;      // leading snp position for each window
    VectorXi windSize;       // number of snps in each window

    // for Eigen dec
    VectorXi blockStarts;    // each LD block startings index in SNP included scale
    VectorXi blockSizes;     // each LD block size;
    VectorXf nGWASblock;     // mean GWAS sample size for each block in GWAS
    VectorXf numSnpsBlock;   // number of SNPs for each block
    VectorXf numEigenvalBlock;  // number of eigenvalues kept for each block
    
    VectorXf LDsamplVar;     // sum of sampling variance of LD for each SNP with all other SNPs; this is for summary-bayes methods
    VectorXf LDscore;        // sum of r^2 over SNPs in significant LD
    
    VectorXf RinverseSqrt;   // sqrt of the weights for the residuals in the individual-level model
    VectorXf Rsqrt;
    
    float ypy;               // y'y the total sum of squares adjusted for the mean
    float varGenotypic;
    float varResidual;
    float varPhenotypic;
    float varRandom;         // variance explained by random covariate effects
    
    bool reindexed;
    bool sparseLDM;
    bool shrunkLDM;
    bool readLDscore;
    bool makeWindows;
    bool weightedRes;
    
    bool lowRankModel;
    
    vector<SnpInfo*> snpInfoVec;
    vector<IndInfo*> indInfoVec;

    vector<AnnoInfo*> annoInfoVec;
    vector<string> annoNames;
    vector<string> snpAnnoPairNames;
    
    map<string, SnpInfo*> snpInfoMap;
    map<string, IndInfo*> indInfoMap;

    vector<SnpInfo*> incdSnpInfoVec;
    vector<IndInfo*> keptIndInfoVec;
    
    vector<GeneInfo*> geneInfoVec;
    
    vector<string> fixedEffectNames;
    vector<string> randomEffectNames;
    vector<string> snpEffectNames;
    
    set<int> chromosomes;
    vector<ChromInfo*> chromInfoVec;
    
    vector<bool> fullSnpFlag;
    
    vector<unsigned> numSnpMldVec;
    vector<unsigned> numSnpAnnoVec;
    VectorXf numAnnoPerSnpVec;
    
    vector<SpMat> annowiseZPZsp;
    vector<VectorXf> annowiseZPZdiag;
    
    vector<vector<unsigned> > windowSnpIdxVec;
    
    map<SnpInfo*, vector<SnpInfo*> > LDmap;
    
    //////// ld block begin ///////
     vector<LDBlockInfo *> ldBlockInfoVec;
     vector<LDBlockInfo *> keptLdBlockInfoVec;
     map<string, LDBlockInfo *> ldBlockInfoMap;
     vector<string> ldblockNames;
     vector<VectorXf> eigenValLdBlock; // store lambda  (per LD block matrix = U * diag(lambda)* V')  per gene LD
     vector<MatrixXf> eigenVecLdBlock; // store U   (per  LD block matrix = U * diag(lambda)* V')  per gene LD
     vector<VectorXf> wcorrBlocks;
     vector<MatrixXf> Qblocks;
     ///////// ld block end  ////////
    ///
    map<int, vector<int>> ldblock2gwasSnpMap;

    vector<VectorXf> gwasEffectInBlock;  // gwas marginal effect;
    vector<VectorXf> gwasPerSnpNinBlock; // gwas per-SNP sample size;
    
    vector<VectorXf> pseudoGwasEffectTrn;
    vector<VectorXf> pseudoGwasEffectVal;
    VectorXf pseudoGwasNtrnBlock;
    VectorXf pseudoGwasNValBlock;
    VectorXf b_val;

    unsigned numFixedEffects;
    unsigned numRandomEffects;
    unsigned numSnps;
    unsigned numInds;
    unsigned numIncdSnps;
    unsigned numKeptInds;
    unsigned numChroms;
    unsigned numSkeletonSnps;
    unsigned numAnnos;
    unsigned numWindows;
    unsigned numLDBlocks;
    unsigned numKeptLDBlocks;
    
    string label;
    string title;
    
    Data(){
        numFixedEffects = 0;
        numRandomEffects = 0;
        numSnps = 0;
        numInds = 0;
        numIncdSnps = 0;
        numKeptInds = 0;
        numChroms = 0;
        numSkeletonSnps = 0;
        numAnnos = 0;
        numWindows= 0;
        numLDBlocks = 0;
        numKeptLDBlocks = 0;
        
        reindexed = false;
        sparseLDM = false;
        readLDscore = false;
        makeWindows = false;
        weightedRes = false;
        lowRankModel = false;
    }
    
    void readFamFile(const string &famFile);
    void readBimFile(const string &bimFile);
    void readBedFile(const bool noscale, const string &bedFile);
    void readPhenotypeFile(const string &phenFile, const unsigned mphen);
    void readCovariateFile(const string &covarFile);
    void readRandomCovariateFile(const string &covarFile);
    void readGwasSummaryFile(const string &gwasFile, const float afDiff, const float mafmin, const float mafmax, const float pValueThreshold, const bool imputeN, const bool removeOutlierN);
    void readLDmatrixInfoFileOld(const string &ldmatrixFile);
    void readLDmatrixInfoFile(const string &ldmatrixFile);
    void readLDmatrixBinFile(const string &ldmatrixFile);
    void readLDmatrixTxtFile(const string &ldmatrixFile);
    void readGeneticMapFile(const string &freqFile);
    void readfreqFile(const string &geneticMapFile);
    void readGeneMapFile(const string &geneMapFile, const int flank, const string &genomeBuild);
    void keepMatchedInd(const string &keepIndFile, const unsigned keepIndMax);
    void includeSnp(const string &includeSnpFile);
    void excludeSnp(const string &excludeSnpFile);
    void includeChr(const unsigned chr);
    void includeBlock(const unsigned block);
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
    string partLDMatrix(const string &partParam, const string &outfilename, const string &LDmatType);
    void makeLDmatrix(const string &bedFile, const string &LDmatType, const float chisqThreshold, const float LDthreshold, const unsigned windowWidth,
                      const string &snpRange, const string &filename, const bool writeLdmTxt);
    void makeshrunkLDmatrix(const string &bedFile, const string &LDmatType, const string &snpRange, const string &filename, const bool writeLdmTxt, const float effpopNE, const float cutOff, const float genMapN);
    void resizeWindow(const vector<SnpInfo*> &incdSnpInfoVec, const VectorXi &windStartOri, const VectorXi &windSizeOri,
                      VectorXi &windStartNew, VectorXi &windSizeNew);
    void computeAlleleFreq(const MatrixXf &Z, vector<SnpInfo*> &incdSnpInfoVec, VectorXf &snp2pq);
    void reindexSnp(vector<SnpInfo*> snpInfoVec);
    void initVariances(const float heritability, const float propVarRandom);
    
    void outputSnpResults(const VectorXf &posteriorMean, const VectorXf &posteriorSqrMean, const VectorXf &pip, const bool noscale, const string &filename) const;
//    void outputFixedEffects(const MatrixXf &fixedEffects, const string &filename) const;
    void outputFixedEffects(const VectorXf &mean, const VectorXf &sd, const string &filename) const;
//    void outputRandomEffects(const MatrixXf &randomEffects, const string &filename) const;
    void outputRandomEffects(const VectorXf &mean, const VectorXf &sd, const string &filename) const;
    void outputWindowResults(const VectorXf &posteriorMean, const string &filename) const;
    void summarizeSnpResults(const SpMat &snpEffects, const string &filename) const;
    void buildSparseMME(const bool sampleOverlap, const bool noscale);
    void readMultiLDmatInfoFile(const string &mldmatFile);
    void readMultiLDmatBinFile(const string &mldmatFile);
    void outputSnpEffectSamples(const SpMat &snpEffects, const unsigned burnin, const unsigned outputFreq, const string &snpResFile, const string &filename) const;
    void resizeLDmatrix(const string &LDmatType, const float chisqThreshold, const unsigned windowWidth, const float LDthreshold, const float effpopNE, const float cutOff, const float genMapN);
    void outputLDmatrix(const string &LDmatType, const string &filename, const bool writeLdmTxt) const;
    void displayAverageWindowSize(const VectorXi &windSize);
    
    void inputMatchedSnpResults(const string &snpResFile);
    void inputSnpInfoAndResults(const string &snpResFile, const string &bayesType);
    void readLDmatrixBinFileAndShrink(const string &ldmatrixFile);
    void readMultiLDmatBinFileAndShrink(const string &mldmatFile, const float genMapN);
    void directPruneLDmatrix(const string &ldmatrixFile, const string &outLDmatType, const float chisqThreshold, const string &title, const bool writeLdmTxt);
    void jackknifeLDmatrix(const string &ldmatrixFile, const string &outLDmatType, const string &title, const bool writeLdmTxt);
    void addLDmatrixInfo(const string &ldmatrixFile);
    
    void readAnnotationFile(const string &annotationFile, const bool transpose, const bool allowMultiAnno);
    void readAnnotationFileFormat2(const string &continuousAnnoFile, const unsigned flank, const string &eQTLFile); // for continuous annotations
    void setAnnoInfoVec(void);
    void readLDscoreFile(const string &ldscFile);
    void makeAnnowiseSparseLDM(const vector<SparseVector<float> > &ZPZsp, const vector<AnnoInfo *> &annoInfoVec, const vector<SnpInfo*> &snpInfoVec);
    void imputePerSnpSampleSize(vector<SnpInfo*> &snpInfoVec, unsigned &numIncdSnps, float sd);
    void getZPZspmat(void);
    void getZPZmat(void);
    void binSnpByLDrsq(const float rsqThreshold, const string &title);
    void readWindowFile(const string &windowFile);
    void binSnpByWindowID(void);
    void filterSnpByLDrsq(const float rsqThreshold);
    void readResidualDiagFile(const string &resDiagFile);
    void makeWindowAnno(const string &annoFile, const float windowWidth);
    
    void mergeLdmInfo(const string &outLDmatType, const string &dirname, const bool print);
    void inputNewSnpResults(const string &snpResFile);
    void getOverlapWindows(const unsigned windowWidth, const unsigned stepSize);
    
    void readPlinkAFfile(const string &plinkAFfile);
    void readPlinkLDtxtfile(const string &plinkLDfile);
    void readPlinkLDbinfile(const string &plinkLDfile);
    
    void filterSnpByGelmanRubinStat(const float threshold);

    /////////// eigen decomposition for LD blocks
    void readLDBlockInfoFile(const string &ldBlockInfoFile);
    void getEigenDataFromFullLDM(const string &filename, const float eigenCutoff);

    void eigenDecomposition(const MatrixXf &X, const float &prop, VectorXf &eigenValAdjusted, MatrixXf &eigenVecAdjusted, float &sumPosEigVal);
    MatrixXf generateLDmatrixPerBlock(const string &bedFile, const vector<string> &snplists); // generate full LDM for block
    
    void makeBlockLDmatrix(const string &bedFile, const string &LDmatType, const unsigned block, const string &filename, const bool writeLdmTxt, int ldBlockRegionWind = 0);

    void readBlockLdmBinaryAndDoEigenDecomposition(const string &dirname, const unsigned block, const float eigenCutoff, const bool writeLdmTxt);
    
    void getEigenDataForLDBlock(const string &bedFile, const string &ldBlockInfoFile, int ldBlockRegionWind, const string &filename, const float eigenCutoff);
    void outputBlockLDmatrixInfo(const LDBlockInfo &block, const string &outSnpfile, const string &outldmfile) const;

    void impG(const unsigned block, double diag_mod = 0.1);

    ///////////// read LD matrix eigen-decomposition data for LD blocks
    void readEigenMatrix(const string &eigenMatrixFile, const float eigenCutoff, const bool readBinary = false, const bool writeLdmTxt = false, const string &outputDir = ".");
    void readBlockLDmatrixAndDoEigenDecomposition(const string &LDmatrixFile, const unsigned block, const float eigenCutoff, const bool writeLdmTxt);
    void readBlockLdmInfoFile(const string &dirname, const unsigned block = 0);
    void readBlockLdmSnpInfoFile(const string &dirname, const unsigned block = 0);
    void readBlockLDMbinaryFile(const string &svdLDfile, const float eigenCutoff);
    vector<LDBlockInfo *> makeKeptLDBlockInfoVec(const vector<LDBlockInfo *> &ldBlockInfoVec);
    
    void readEigenMatrixBinaryFile(const string &eigenMatrixFile, const float eigenCutoff, const bool writeLdmTxt = false, const string &outputDir = ".");
    
    void readEigenMatrixBinaryFileAndMakeWandQ(const string &dirname, const float eigenCutoff, const vector<VectorXf> &GWASeffects, const VectorXf &nGWASblock, const bool noscale, const bool makePseudoSummary);

    
    ///////////// merge eigen matrices
    void mergeMultiEigenLDMatrices(const string & infoFile, const string &filename, const string LDmatType);

    //////////// Step 2.2 Build multiple maps
    void buildMMEeigen(const string &dirname, const bool sampleOverlap, const float eigenCutoff, const bool noscale); // for eigen decomposition
    void includeMatchedBlocks(void);

    //////////// Step 2.3 build model matrix
//    void constructWandQ(const float eigenCutoff, const bool noscale);
 
    void imputeSummaryData(void);
    
    void truncateEigenMatrix(const float sumPosEigVal, const float eigenCutoff, const VectorXf &oriEigenVal, const MatrixXf &oriEigenVec, VectorXf &newEigenVal, MatrixXf &newEigenVec);
    void constructPseudoSummaryData(void);
    
    void constructWandQ(const vector<VectorXf> &GWASeffects, const float nGWAS, const bool noscale);
    
    void scaleGwasEffects(void);
    void mapSnpsToBlocks(void);
    
    void mergeBlockGwasSummary(const string &gwasSummaryFile, const string &title);

    void outputWandQ(const string &dirname);
    void readUnconvergedSnplist(const string &snplistFile);
    
    void convert(const string &eigenMatrixFile, const string &snplistFile, const string &title);
    
    void getLDfromEigenMatrix(const string &eigenMatrixFile, const float rsqThreshold, const string &title);
    
    void inputPairwiseLD(const string &ldfile, const float rsqThreshold);
    void inputLDfriends(const string &ldfriendFile);
    void getLDfriends(const string &ldfile, const float rsqThreshold, const string &title);
    
    void outputEigenMatTxt(const string &title);
    void skipSnp(const string &skipSnpFile);
    
    void readBlockLDmatrixAndMakeItSparse(const string &LDmatrixFile, const unsigned block, const float chisqThreshold, const bool writeLdmTxt);
    void readBlockLdmBinaryAndMakeItSparse(const string &dirname, const unsigned block, const float chisqThreshold, const bool writeLdmTxt);

    void readSparseBlockLDmatrixAndDoEigenDecomposition(const string &LDmatrixFile, const unsigned block, const float eigenCutoff, const bool writeLdmTxt);
    void readSparseBlockLdmBinaryAndDoEigenDecomposition(const string &dirname, const unsigned block, const float eigenCutoff, const bool writeLdmTxt);
    
    void resizeBlockLDmatrix(const string &dirname, const string &outLDmatType, const string &includeSnpFile, const string &title, const bool writeLdmTxt);

    void outputBlockLDmatrixTxt(const string &dirname, const unsigned block);
    
    void resizeBlockLDmatrixAndDoEigenDecomposition(const string &LDmatrixFile, const float eigenCutoff, const float rsqThreshold, const string &title, const bool writeLdmTxt);
    
    void readBlockLDmatrix(const string &dirname, const string &blockID, const int32_t blockSize, MatrixXf &ldm);
    void convertToBlockTriangularMatrix(const string &dirname, const bool writeLdmTxt, const string &title);
    void outputLDfriends(const MatrixXf &ldm, const LDBlockInfo *blockInfo, const string &outDirname);
    
    void calcMarginalEnrichmentPermute(const string &paramStr, const string &title);
    void calcMarignalEnrichmentJackknife(const string &paramStr, const string &title);
    void calcJointEnrichmentPermute(const string &paramStr, const string &title);
    void calcJointEnrichmentJackknife(const string &paramStr, const string &title);
    void calcJointEnrichmentJackknifeLM(const string &paramStr, const string &title);

};

#endif /* data_hpp */
