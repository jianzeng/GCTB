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
    enum {dense, sparse} storageMode;
    
    unsigned chainLength;
    unsigned burnin;
    unsigned thin;
    
    unsigned nrow;
    unsigned ncol;
    unsigned nnz;  // number of non-zeros for sparse matrix
    
    MatrixXf datMat;
    SpMat datMatSp; // most of the snp effects will be zero if pi value is high
    
    VectorXf posteriorMean;
    VectorXf posteriorSqrMean;
    VectorXf pip;  // for snp effects, will consider to remove
    VectorXf lastSample; // save the last sample of MCMC
    
    FILE *bout;
    ofstream tout;
    
    McmcSamples(const string &label, const unsigned chainLength, const unsigned burnin, const unsigned thin,
                const unsigned npar, const string &storage_mode = "dense"):
    label(label), chainLength(chainLength), burnin(burnin), thin(thin) {
        nrow = chainLength/thin - burnin/thin;
        ncol = npar;
        if (storage_mode == "dense") {
            storageMode = dense;
            datMat.setZero(nrow, ncol);
        } else if (storage_mode == "sparse") {
            storageMode = sparse;
            //if (myMPI::rank==0) datMatSp.reserve(VectorXi::Constant(ncol,nrow));  // for faster filling the matrix
        } else {
            cerr << "Error: Unrecognized storage mode: " << storage_mode << endl;
        }
        posteriorMean.setZero(ncol);
        posteriorSqrMean.setZero(ncol);
        pip.setZero(ncol);
        lastSample.setZero(ncol);
    }
    
    McmcSamples(const string &label): label(label) {}
    
    void getSample(const unsigned iter, const VectorXf &sample, const bool writeBinPosterior, const bool writeTxtPosterior);
    void getSample(const unsigned iter, const float sample, const bool writeTxtPosterior, ofstream &out);
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
};

class MCMC {
private:
    string outfilename;
    ofstream out;
    
    void initTxtFile(const vector<Parameter*> &paramVec, const string &title);
    vector<McmcSamples*> initMcmcSamples(const Model &model, const unsigned chainLength, const unsigned burnin,
                                         const unsigned thin, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior);
    void collectSamples(const Model &model, vector<McmcSamples*> &mcmcSampleVec, const unsigned iteration, const bool writeBinPosterior, const bool writeTxtPosterior);
    void printStatus(const vector<Parameter*> &paramToPrint, const unsigned thisIter, const unsigned outputFreq, const string &timeLeft);
    void printStatusR(const vector<float*> &paramToPrintR, const unsigned thisIter, const unsigned outputFreq, const string &timeLeft);
    void printSummary(const vector<Parameter*> &paramToPrint, const vector<McmcSamples*> &mcmcSampleVec, const string &filename);
    void printSetSummary(const vector<ParamSet*> &paramSetToPrint, const vector<McmcSamples*> &mcmcSampleVec, const string &filename);
    void printSnpAnnoMembership(const vector<ParamSet*> &paramSetToPrint, const vector<McmcSamples*> &mcmcSampleVec, const string &filename);

public:
    vector<McmcSamples*> run(Model &model, const unsigned chainLength, const unsigned burnin, const unsigned thin, const bool print,
                             const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior);
    void convergeDiagGelmanRubin(const Model &model, vector<vector<McmcSamples*> > &mcmcSampleVecChain, const string &filename);
};

#endif /* mcmc_hpp */
