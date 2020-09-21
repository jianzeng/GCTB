//
//  hsq.hpp
//  gctb
//
//  Created by Jian Zeng on 20/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#ifndef hsq_hpp
#define hsq_hpp

#include <stdio.h>
#include "data.hpp"
#include "mcmc.hpp"


class Heritability {
public:
    McmcSamples varGenotypic;
    McmcSamples varResidual;
    McmcSamples hsq;
    McmcSamples numSnpVarExplained;  // No. SNPs explain more than 0.0001, 0.001, 0.01 genetic variance
    McmcSamples piVarExplained;
    McmcSamples meanEffSqVarExplained;  // the mean of effect squared for SNPs that explain more than 0.0001, 0.001, 0.01 genetic variance
    VectorXf gamma;
    
    unsigned numMcmcSamples;
    unsigned popSize;
    unsigned numIncdSnps;
    
    void getEstimate(const Data &data, const McmcSamples &snpEffects, const McmcSamples &resVar, const unsigned outputFreq);
    void writeRes(const string &filename);
    void writeMcmcSamples(const string &filename);
    
    Heritability(const unsigned nsamples): varGenotypic("GenVar"), varResidual("ResVar"), hsq("hsq"), numSnpVarExplained("numSnpVarExplained"), piVarExplained("piVarExplained"), meanEffSqVarExplained("meanEffSqVarExplained") {
        numMcmcSamples = nsamples;
        varGenotypic.storageMode = McmcSamples::dense;
        varGenotypic.datMat.resize(nsamples,1);
        varGenotypic.nrow = nsamples;
        varGenotypic.ncol = 1;
        varResidual.storageMode = McmcSamples::dense;
        varResidual.datMat.resize(nsamples,1);
        varResidual.nrow = nsamples;
        varResidual.ncol = 1;
        hsq.storageMode = McmcSamples::dense;
        hsq.datMat.resize(nsamples,1);
        hsq.nrow = nsamples;
        hsq.ncol = 1;
        gamma.resize(3);
        gamma << 1e-4, 1e-3, 1e-2;
        numSnpVarExplained.storageMode = McmcSamples::dense;
        numSnpVarExplained.datMat.resize(nsamples, gamma.size());
        numSnpVarExplained.nrow = nsamples;
        numSnpVarExplained.ncol = gamma.size();
        piVarExplained.storageMode = McmcSamples::dense;
        piVarExplained.datMat.resize(nsamples, gamma.size());
        piVarExplained.nrow = nsamples;
        piVarExplained.ncol = gamma.size();
        meanEffSqVarExplained.storageMode = McmcSamples::dense;
        meanEffSqVarExplained.datMat.resize(nsamples, gamma.size());
        meanEffSqVarExplained.nrow = nsamples;
        meanEffSqVarExplained.ncol = gamma.size();
    }
};

class Polygenicity {
public:
    McmcSamples pi;
    McmcSamples nnz;
    McmcSamples numSnpVarExplained;  // No. SNPs explain more than 0.0001, 0.001, 0.01 genetic variance
    McmcSamples piVarExplained;
    McmcSamples meanEffSqVarExplained;  // the mean of effect squared for SNPs that explain more than 0.0001, 0.001, 0.01 genetic variance
    VectorXf gamma;

    unsigned numMcmcSamples;
    unsigned numIncdSnps;
    
    void getEstimate(const Data &data, const McmcSamples &snpEffects, const McmcSamples &genVar, const unsigned outputFreq);
    void writeRes(const string &filename);
    void writeMcmcSamples(const string &filename);
    
    Polygenicity(const unsigned nsamples): pi("Pi"), nnz("NnzSnp"), numSnpVarExplained("numSnpVarExplained"), piVarExplained("piVarExplained"), meanEffSqVarExplained("meanEffSqVarExplained") {
        numMcmcSamples = nsamples;
        pi.storageMode = McmcSamples::dense;
        pi.datMat.resize(nsamples,1);
        pi.nrow = nsamples;
        pi.ncol = 1;
        nnz.storageMode = McmcSamples::dense;
        nnz.datMat.resize(nsamples,1);
        nnz.nrow = nsamples;
        nnz.ncol = 1;
        gamma.resize(4);
        gamma << 1e-5, 1e-4, 1e-3, 1e-2;
        numSnpVarExplained.storageMode = McmcSamples::dense;
        numSnpVarExplained.datMat.resize(nsamples, gamma.size());
        numSnpVarExplained.nrow = nsamples;
        numSnpVarExplained.ncol = gamma.size();
        piVarExplained.storageMode = McmcSamples::dense;
        piVarExplained.datMat.resize(nsamples, gamma.size());
        piVarExplained.nrow = nsamples;
        piVarExplained.ncol = gamma.size();
        meanEffSqVarExplained.storageMode = McmcSamples::dense;
        meanEffSqVarExplained.datMat.resize(nsamples, gamma.size());
        meanEffSqVarExplained.nrow = nsamples;
        meanEffSqVarExplained.ncol = gamma.size();
    }
};

#endif /* hsq_hpp */
