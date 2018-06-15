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
    
    varGenotypic.datMat.resize(nSamples,1);
    
    for (unsigned i=0; i<nSamples; ++i) {
        VectorXf gi;
        gi.setZero(data.numKeptInds);
        for (unsigned j=0; j<data.numSnps; ++j) {
            SnpInfo *snp = data.snpInfoVec[j];
            if (snp->included)
                gi += data.Z.col(snp->index) * snpEffects.datMatSp.coeff(i, j);
        }
        varGenotypic.datMat(i,0) = Gadget::calcVariance(gi);
        if (!(i%100)) cout << " iter " << i << " " << varGenotypic.datMat(i,0) << endl;
    }

    varResidual.datMat = phenVar - varGenotypic.datMat.array();
    hsq.datMat = varGenotypic.datMat.array()/phenVar;
    popSize = data.numKeptInds;
    
    timer.getTime();
    
    cout << endl;
    cout << "Genotypic variance :  " << varGenotypic.mean() << "  " << varGenotypic.sd() << endl;
    cout << "Residual variance  :  " << varResidual.mean()  << "  " << varResidual.sd()  << endl;
    cout << "Heritability       :  " << hsq.mean()          << "  " << hsq.sd()          << endl;
    cout << "Population size    :  " << popSize << endl;
    cout << "Computational time :  " << timer.format(timer.getElapse()) << endl << endl;
}

void Heritability::writeRes(const string &filename){
    string file = filename + ".hsq";
    ofstream out(file.c_str());
    if (!out) {
        throw("Error: cannot open file " + file);
    }
    out << "Genotypic variance :  " << varGenotypic.mean() << "  " << varGenotypic.sd() << endl;
    out << "Residual variance  :  " << varResidual.mean()  << "  " << varResidual.sd()  << endl;
    out << "Heritability       :  " << hsq.mean()          << "  " << hsq.sd()          << endl;
    out << "Population size    :  " << popSize << endl;
    out.close();
}

void Heritability::writeMcmcSamples(const string &filename){
    string file = filename + ".hsq.mcmc";
    ofstream out(file.c_str());
    if (!out) {
        throw("Error: cannot open file " + file);
    }
    out << hsq.datMat;
    out.close();
}