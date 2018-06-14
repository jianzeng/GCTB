//
//  hsq.cpp
//  gctb
//
//  Created by Jian Zeng on 20/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "hsq.hpp"

void Heritability::getEstimate(const Data &data, const McmcSamples &snpEffects, const float phenVar){
    // This function computes Za (Z: SNP genotypes, a:sampled SNP effects) for each MCMC cycle,
    // computes the variance of Za, and then calculate the posterior mean across cycles
    
    Gadget::Timer timer;
    timer.setTime();
    
    cout << "Estimating heritability from MCMC samples of SNP effects ..." << endl;
    
    unsigned nSamples = snpEffects.datMatSp.rows();
    
    VectorXf genVarMcmc(nSamples);
    
    hsqMcmc.resize(nSamples);
    
    for (unsigned i=0; i<nSamples; ++i) {
        VectorXf gi;
        gi.setZero(data.numKeptInds);
        for (unsigned j=0; j<data.numSnps; ++j) {
            SnpInfo *snp = data.snpInfoVec[j];
            if (snp->included)
                gi += data.Z.col(snp->index) * snpEffects.datMatSp.coeff(i, j);
        }
        genVarMcmc[i] = Gadget::calcVariance(gi);
        if (!(i%100)) cout << " iter " << i << " " << genVarMcmc[i] << endl;
    }

    hsqMcmc = genVarMcmc.array()/phenVar;
    
    varGenotypic = genVarMcmc.mean();
    varResidual  = phenVar - genVarMcmc.mean();
    hsq          = hsqMcmc.mean();
    popSize      = data.numKeptInds;
    
    timer.getTime();
    
    cout << endl;
    cout << "Genotypic variance :  " << varGenotypic  << endl;
    cout << "Residual variance  :  " << varResidual   << endl;
    cout << "Heritability       :  " << hsq           << endl;
    cout << "Population size    :  " << popSize       << endl;
    cout << "Computational time:  " << timer.format(timer.getElapse()) << endl << endl;
}

void Heritability::getEstimate(const Data &data, const McmcSamples &snpEffects, const McmcSamples &resVar){
    // essensially, it computes Za (Z: SNP genotypes, a:sampled SNP effects) for each MCMC cycle,
    // computes the variance of Za, and then calculate the posterior mean across cycles
    
    Gadget::Timer timer;
    timer.setTime();
    
    cout << "Estimating heritability from MCMC samples of SNP effects ..." << endl;

    unsigned nSamples = snpEffects.datMatSp.rows();
    
    VectorXf genVarMcmc(nSamples);
    VectorXf resVarMcmc = resVar.datMat.col(0);
    
    hsqMcmc.resize(nSamples);
    
    MatrixXf gmcmc = data.Z * snpEffects.datMatSp.transpose();
    
    for (unsigned i=0; i<nSamples; ++i) {
        genVarMcmc[i] = Gadget::calcVariance(gmcmc.col(i));
        hsqMcmc[i] = genVarMcmc[i]/(genVarMcmc[i] + resVarMcmc[i]);
    }
    
    varGenotypic = genVarMcmc.mean();
    varResidual  = resVarMcmc.mean();
    hsq          = hsqMcmc.mean();
    popSize      = data.numKeptInds;
    
    timer.getTime();

    cout << "Genotypic variance :  " << varGenotypic  << endl;
    cout << "Residual variance  :  " << varResidual   << endl;
    cout << "Heritability       :  " << hsq           << endl;
    cout << "Population size    :  " << popSize       << endl;
    cout << "Computational time:  " << timer.format(timer.getElapse()) << endl << endl;
}

void Heritability::writeRes(const string &filename){
    string file = filename + ".hsq";
    ofstream out(file.c_str());
    if (!out) {
        throw("Error: cannot open file " + file);
    }
    out << "Genotypic variance :  " << varGenotypic  << endl;
    out << "Residual variance  :  " << varResidual   << endl;
    out << "Heritability       :  " << hsq           << endl;
    out << "Population size    :  " << popSize       << endl;
    out.close();
}

void Heritability::writeMcmcSamples(const string &filename){
    string file = filename + ".hsq.mcmc";
    ofstream out(file.c_str());
    if (!out) {
        throw("Error: cannot open file " + file);
    }
    out << hsqMcmc;
    out.close();
}