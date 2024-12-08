//
//  mcmc.hpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#ifndef mcmc_hpp
#define mcmc_hpp
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

#include <stdio.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <boost/format.hpp>
#include "model.hpp"
#include "gadgets.hpp"

using namespace std;
using namespace Eigen;


class McmcSamples {
    // rows: MCMC cycles, cols: model parameters
public:
    const string label;
    string filename;

    enum {dense, sparse, do_not_store} storageMode;
    enum {bin, txt, txt_combine_others, no_output} outputMode;
    
    unsigned chainLength;
    unsigned burnin;
    unsigned thin;
    
    unsigned nrow;
    unsigned ncol;
    unsigned nnz;  // number of non-zeros for sparse matrix
    
    MatrixXf datMat;
    SpMat datMatSp; // most of the snp effects will be zero if pi value is high
    
    MatrixXf sampleIter;
    
    VectorXf posteriorMean;
    VectorXf posteriorSqrMean;
//    VectorXf pip;  // for snp effects, will consider to remove
//    VectorXf lastSample; // save the last sample of MCMC
    
    // for multiple chains
    bool multiChain;
    unsigned numChains;
    MatrixXf perChainMean;
    MatrixXf perChainSqrMean;
    VectorXf GelmanRubinStat;
    unsigned cntPosteriorSample;

    VectorXf probGreaterThanCriticalValue;
    float criticalValue;
    
    FILE *bout;
    ofstream tout;
        
    McmcSamples(const string &label, const unsigned numChains, const unsigned chainLength, const unsigned burnin, const unsigned thin,
                const unsigned npar, const string &storage_mode, const string &output_mode, const string &title):
    label(label), numChains(numChains), chainLength(chainLength), burnin(burnin), thin(thin) {
        nrow = (chainLength/thin - burnin/thin)*numChains;
        ncol = npar;
        
        if (storage_mode == "dense") {
            storageMode = dense;
            datMat.setZero(nrow, ncol);
        } else if (storage_mode == "sparse") {
            storageMode = sparse;
            //if (myMPI::rank==0) datMatSp.reserve(VectorXi::Constant(ncol,nrow));  // for faster filling the matrix
        } else if (storage_mode == "do_not_store") {
            storageMode = do_not_store;
        } else {
            cerr << "Error: Unrecognized storage mode: " << storage_mode << ". Option is 'dense' or 'sparse'." << endl;
        }
        
        if (output_mode == "bin") {
            outputMode = bin;
            initBinFile(title);
        } else if (output_mode == "txt") {
            outputMode = txt;
            initTxtFile(title);
        } else if (output_mode == "txt_combine_others") {
            outputMode = txt_combine_others;
        } else if (output_mode == "no_output") {
            outputMode = no_output;
        } else {
            cerr << "Error: Unrecognized output mode: " << output_mode << ". Option is 'bin', 'txt', 'txt_combine_others' or 'no_output'." << endl;
        }
        
        posteriorMean.setZero(ncol);
        posteriorSqrMean.setZero(ncol);
        probGreaterThanCriticalValue.setZero(ncol);
        sampleIter.resize(npar, numChains);
        cntPosteriorSample = 0;
        
        if (numChains > 1) {
            multiChain = true;
            perChainMean.setZero(npar, numChains);
            perChainSqrMean.setZero(npar, numChains);
            GelmanRubinStat.setZero(npar);
        } else {
            multiChain = false;
        }
        
        Gadget::Tokenizer token;
        token.getTokens(label, "_");
        if (token.back() == "Enrichment") {
            criticalValue = 1;
        } else {
            criticalValue = 0;
        }
    }
    
    McmcSamples(const string &label): label(label) {}
    
    void getParSample(const unsigned iter, const Parameter* par);
    void getParSetSample(const unsigned iter, const ParamSet* parSet);
    void outputSample(const unsigned chain, ofstream &out);
    void computeGelmanRubinStat(void);
    void writeSampleBin(const unsigned iter, const VectorXf &sample, const string &title);
    void writeSampleTxt(const unsigned iter, const float sample, const string &title);
    VectorXf mean(void);
    VectorXf sd(void);
    
    void initBinFile(const string &title);
    void initTxtFile(const string &title);
    void writeDataBin(const string &title);
    void writeDataTxt(const string &title);
    void readDataBin(const string &filename);
    void readDataTxt(const string &filename);
    void readDataTxt(const string &filename, const string &label);
    void writeMatSpTxt(const string &title);
};

class MCMC {
public:
    string outfilename;
    ofstream out;
    
    enum {keep_running, restart_and_use_robust_model, stop_and_exit} action;
    
    void initTxtFile(const vector<Parameter*> &paramVec, const string &title);
    vector<McmcSamples*> initMcmcSamples(const Model &model, const unsigned numChains, const unsigned chainLength, const unsigned burnin,
                                         const unsigned thin, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior);
    void collectSamples(const Model &model, vector<McmcSamples*> &mcmcSampleVec, const unsigned iteration);
    void outputSamples(vector<McmcSamples*> &mcmcSampleVec, const unsigned numChains);
    void computeGelmanRubinStat(vector<McmcSamples*> &mcmcSampleVec);
    void printStatus(const vector<Parameter*> &paramToPrint, const unsigned thisIter, const unsigned outputFreq, const string &timeLeft);
    void printStatusR(const vector<float*> &paramToPrintR, const unsigned thisIter, const unsigned outputFreq, const string &timeLeft);
    void printSummary(const vector<Parameter*> &paramToPrint, const vector<McmcSamples*> &mcmcSampleVec, const unsigned numChains, const string &filename);
    void printSetSummary(const vector<ParamSet*> &paramSetToPrint, const vector<McmcSamples*> &mcmcSampleVec, const unsigned numChains, const string &filename, const string &enrich);
    void printSnpAnnoMembership(const vector<ParamSet*> &paramSetToPrint, const vector<McmcSamples*> &mcmcSampleVec, const string &filename);

    MCMC() {
        action = keep_running;
    }

    vector<McmcSamples*> run(Model &model, const unsigned numChains, const unsigned chainLength, const unsigned burnin, const unsigned thin, const bool print,
                             const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior);
    void convergeDiagGelmanRubin(const Model &model, vector<vector<McmcSamples*> > &mcmcSampleVecChain, const string &filename);
    void setAction(const Model &model);
};

#endif /* mcmc_hpp */
