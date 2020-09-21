//
//  hsq.cpp
//  gctb
//
//  Created by Jian Zeng on 20/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "hsq.hpp"

void Heritability::getEstimate(const Data &data, const McmcSamples &snpEffects, const McmcSamples &resVar, const unsigned outputFreq){
    // This function computes Za (Z: SNP genotypes, a:sampled SNP effects) for each MCMC cycle,
    // computes the variance of Za, and then calculate the posterior mean across cycles
    
    Gadget::Timer timer;
    timer.setTime();
    
    cout << "Estimating heritability from MCMC samples of SNP effects ..." << endl;
        
    popSize = data.numKeptInds;
    numIncdSnps = data.numIncdSnps;
    
    VectorXf perSnpGV(data.numIncdSnps);
    VectorXf gammaVarExplained(gamma.size());
    float betaj = 0.0;
    
    for (unsigned i=0; i<numMcmcSamples; ++i) {
        perSnpGV.setZero(data.numIncdSnps);
        if (data.Z.size()) {  // using individual genotype data
            VectorXf gi;
            gi.setZero(data.numKeptInds);
            for (unsigned j=0; j<data.numIncdSnps; ++j) {
                SnpInfo *snp = data.incdSnpInfoVec[j];
                betaj = snpEffects.datMatSp.coeff(i, j);
                if (betaj) {
                    gi += data.Z.col(j) * betaj;
                    perSnpGV[j] = data.snp2pq[j] * betaj * betaj;
                }
            }
            varGenotypic.datMat(i,0) = Gadget::calcVariance(gi);
        }
        else if (data.sparseLDM) {  // using sparse LD matrix
            VectorXf beta(snpEffects.datMatSp.row(i));
            varGenotypic.datMat(i,0) = beta.transpose() * data.ZPZspmat * beta;
            varGenotypic.datMat(i,0) /= float(data.numKeptInds);
            perSnpGV = data.snp2pq.cwiseProduct(beta.cwiseProduct(beta));
        }
        else {  // using full LD matrix
            //cout << snpEffects.datMatSp.row(i).size() << " " << data.ZPZmat.rows() << " " << data.ZPZmat.cols() << endl;
            VectorXf beta(snpEffects.datMatSp.row(i));
            varGenotypic.datMat(i,0) = beta.transpose() * data.ZPZmat * beta;
            varGenotypic.datMat(i,0) /= float(data.numKeptInds);
            perSnpGV = data.snp2pq.cwiseProduct(beta.cwiseProduct(beta));
        }
        varResidual.datMat(i,0) = resVar.datMat(i,0);
        hsq.datMat(i,0) = varGenotypic.datMat(i,0)/(varGenotypic.datMat(i,0) + varResidual.datMat(i,0));
        numSnpVarExplained.datMat.row(i).setZero();
        meanEffSqVarExplained.datMat.row(i).setZero();
        gammaVarExplained = gamma * varGenotypic.datMat(i,0);
        for (unsigned t=0; t<gamma.size(); ++t) {
            for (unsigned j=0; j<data.numIncdSnps; ++j) {
                if (perSnpGV[j] > gammaVarExplained[t]) {
                    numSnpVarExplained.datMat(i,t)++;
                    meanEffSqVarExplained.datMat(i,t) += snpEffects.datMatSp.coeff(i,j)*snpEffects.datMatSp.coeff(i,j);
                }
            }
        }
        for (unsigned t=0; t<gamma.size()-1; ++t) {
            numSnpVarExplained.datMat(i,t) -= numSnpVarExplained.datMat(i,t+1);
            meanEffSqVarExplained.datMat(i,t) -= meanEffSqVarExplained.datMat(i,t+1);
        }
        meanEffSqVarExplained.datMat.row(i).array() /= numSnpVarExplained.datMat.row(i).array();
        for (unsigned t=0; t<gamma.size(); ++t) {
            if (!numSnpVarExplained.datMat(i,t)) meanEffSqVarExplained.datMat(i,t) = 0;
        }
        piVarExplained.datMat.row(i) = numSnpVarExplained.datMat.row(i)/float(numIncdSnps);
        if (!(i%outputFreq)) cout << " iter " << i << " hsq " << hsq.datMat(i,0) << endl;
    }

    timer.getTime();
    
    VectorXf numSnpVarExplainedMean = numSnpVarExplained.mean();
    VectorXf numSnpVarExplainedSD = numSnpVarExplained.sd();
    VectorXf piVarExplainedMean = piVarExplained.mean();
    VectorXf piVarExplainedSD = piVarExplained.sd();
    VectorXf meanEffSqVarExplainedMean = meanEffSqVarExplained.mean();
    VectorXf meanEffSqVarExplainedSD = meanEffSqVarExplained.sd();

    cout << endl;
    cout << boost::format("%-20s %-10.6f %-10.6f\n") % "Genotypic variance" % varGenotypic.mean() % varGenotypic.sd();
    cout << boost::format("%-20s %-10.6f %-10.6f\n") % "Residual variance"  % varResidual.mean()  % varResidual.sd();
    cout << boost::format("%-20s %-10.6f %-10.6f\n") % "Heritability"       % hsq.mean()          % hsq.sd();
    for (unsigned t=0; t<gamma.size(); ++t) {
        stringstream ss;
        ss << "NumSnp" << gamma[t] << "hsq";
        cout << boost::format("%-20s %-10.6f %-10.6f\n") % ss.str() % numSnpVarExplainedMean[t] % numSnpVarExplainedSD[t];
    }
    for (unsigned t=0; t<gamma.size(); ++t) {
        stringstream ss;
        ss << "Pi" << gamma[t] << "hsq";
        cout << boost::format("%-20s %-10.6f %-10.6f\n") % ss.str() % piVarExplainedMean[t] % piVarExplainedSD[t];
    }
    for (unsigned t=0; t<gamma.size(); ++t) {
        stringstream ss;
        ss << "MeanEffSq" << gamma[t] << "hsq";
        cout << boost::format("%-20s %-10.6f %-10.6f\n") % ss.str() % meanEffSqVarExplainedMean[t] % meanEffSqVarExplainedSD[t];
    }
    cout << boost::format("%-20s %-10s\n") % "Population size" % popSize;
    cout << boost::format("%-20s %-10s\n") % "Total number of SNPs" % numIncdSnps;
    //cout << "Computational time :  " << timer.format(timer.getElapse()) << endl << endl;
}

void Heritability::writeRes(const string &filename){
    string file = filename + ".hsq";
    ofstream out(file.c_str());
    if (!out) {
        throw("Error: cannot open file " + file);
    }
    
    VectorXf numSnpVarExplainedMean = numSnpVarExplained.mean();
    VectorXf numSnpVarExplainedSD = numSnpVarExplained.sd();
    VectorXf piVarExplainedMean = piVarExplained.mean();
    VectorXf piVarExplainedSD = piVarExplained.sd();
    VectorXf meanEffSqVarExplainedMean = meanEffSqVarExplained.mean();
    VectorXf meanEffSqVarExplainedSD = meanEffSqVarExplained.sd();

    out << boost::format("%-20s %-10.6f %-10.6f\n") % "Genotypic variance" % varGenotypic.mean() % varGenotypic.sd();
    out << boost::format("%-20s %-10.6f %-10.6f\n") % "Residual variance"  % varResidual.mean()  % varResidual.sd();
    out << boost::format("%-20s %-10.6f %-10.6f\n") % "Heritability"       % hsq.mean()          % hsq.sd();
//    for (unsigned t=0; t<gamma.size(); ++t) {
//        out << boost::format("%6s %6s %-6s %-10.6f %-10.6f\n") % "NumSNP" % gamma[t] % "hsq" % numSnpVarExplainedMean[t] % numSnpVarExplainedSD[t];
//    }
    for (unsigned t=0; t<gamma.size(); ++t) {
        stringstream ss;
        ss << "NumSnp" << gamma[t] << "hsq";
        out << boost::format("%-20s %-10.6f %-10.6f\n") % ss.str() % numSnpVarExplainedMean[t] % numSnpVarExplainedSD[t];
    }
    for (unsigned t=0; t<gamma.size(); ++t) {
        stringstream ss;
        ss << "Pi" << gamma[t] << "hsq";
        out << boost::format("%-20s %-10.6f %-10.6f\n") % ss.str() % piVarExplainedMean[t] % piVarExplainedSD[t];
    }
    for (unsigned t=0; t<gamma.size(); ++t) {
        stringstream ss;
        ss << "MeanEffSq" << gamma[t] << "hsq";
        out << boost::format("%-20s %-10.6f %-10.6f\n") % ss.str() % meanEffSqVarExplainedMean[t] % meanEffSqVarExplainedSD[t];
    }
    out << boost::format("%-20s %-10s\n") % "Population size" % popSize;
    out << boost::format("%-20s %-10s\n") % "Total number of SNPs" % numIncdSnps;
    out.close();
}

void Heritability::writeMcmcSamples(const string &filename){
    string file1 = filename + ".hsq.mcmc";
    string file2 = filename + ".NumSnpVarExplained.mcmc";
    string file3 = filename + ".PiVarExplained.mcmc";
    ofstream out1(file1.c_str());
    ofstream out2(file2.c_str());
    ofstream out3(file3.c_str());
    if (!out1) throw("Error: cannot open file " + file1);
    if (!out2) throw("Error: cannot open file " + file2);
    if (!out3) throw("Error: cannot open file " + file3);
    out1 << "hsq" << endl << hsq.datMat;
    for (unsigned i=0; i<numSnpVarExplained.ncol; ++i) out2 << "NumSNP" << gamma[i] << "hsq  ";
    out2 << endl << numSnpVarExplained.datMat;
    for (unsigned i=0; i<piVarExplained.ncol; ++i) out3 << "Pi" << gamma[i] << "hsq  ";
    out3 << endl << piVarExplained.datMat;
    out1.close();
    out2.close();
    out3.close();
}

void Polygenicity::getEstimate(const Data &data, const McmcSamples &snpEffects, const McmcSamples &genVar, const unsigned outputFreq){
    cout << "Estimating polygenicity from MCMC samples of SNP effects ..." << endl;
        
    numIncdSnps = data.numIncdSnps;
    
    VectorXf perSnpGV(data.numIncdSnps);
    VectorXf gammaVarExplained(gamma.size());
    
    for (unsigned i=0; i<numMcmcSamples; ++i) {
        VectorXf beta(snpEffects.datMatSp.row(i));
        nnz.datMat(i,0) = 0;
        for (unsigned j=0; j<numIncdSnps; ++j) {
            if (beta[j]) nnz.datMat(i,0)++;
        }
        pi.datMat(i,0) = nnz.datMat(i,0)/float(numIncdSnps);
        perSnpGV = data.snp2pq.cwiseProduct(beta.cwiseProduct(beta));
        numSnpVarExplained.datMat.row(i).setZero();
        meanEffSqVarExplained.datMat.row(i).setZero();
        gammaVarExplained = gamma * genVar.datMat(i,0);
        for (unsigned t=0; t<gamma.size(); ++t) {
            for (unsigned j=0; j<data.numIncdSnps; ++j) {
                if (perSnpGV[j] > gammaVarExplained[t]) {
                    numSnpVarExplained.datMat(i,t)++;
                    meanEffSqVarExplained.datMat(i,t) += beta[j]*beta[j];
                }
            }
        }
        for (unsigned t=0; t<gamma.size()-1; ++t) {
            numSnpVarExplained.datMat(i,t) -= numSnpVarExplained.datMat(i,t+1);
            meanEffSqVarExplained.datMat(i,t) -= meanEffSqVarExplained.datMat(i,t+1);
        }
        meanEffSqVarExplained.datMat.row(i).array() /= numSnpVarExplained.datMat.row(i).array();
        for (unsigned t=0; t<gamma.size(); ++t) {
            if (!numSnpVarExplained.datMat(i,t)) meanEffSqVarExplained.datMat(i,t) = 0;
        }
        piVarExplained.datMat.row(i) = numSnpVarExplained.datMat.row(i)/nnz.datMat(i,0);
    }

    VectorXf numSnpVarExplainedMean = numSnpVarExplained.mean();
    VectorXf numSnpVarExplainedSD = numSnpVarExplained.sd();
    VectorXf piVarExplainedMean = piVarExplained.mean();
    VectorXf piVarExplainedSD = piVarExplained.sd();
    VectorXf meanEffSqVarExplainedMean = meanEffSqVarExplained.mean();
    VectorXf meanEffSqVarExplainedSD = meanEffSqVarExplained.sd();

    cout << endl;
    cout << boost::format("%-20s %-10.6f %-10.6f\n") % "Pi" % pi.mean() % pi.sd();
    cout << boost::format("%-20s %-10.6f %-10.6f\n") % "NnzSnp" % nnz.mean()  % nnz.sd();
    for (unsigned t=0; t<gamma.size(); ++t) {
        stringstream ss;
        ss << "NumSnp" << gamma[t] << "hsq";
        cout << boost::format("%-20s %-10.6f %-10.6f\n") % ss.str() % numSnpVarExplainedMean[t] % numSnpVarExplainedSD[t];
    }
    for (unsigned t=0; t<gamma.size(); ++t) {
        stringstream ss;
        ss << "Pi" << gamma[t] << "hsq";
        cout << boost::format("%-20s %-10.6f %-10.6f\n") % ss.str() % piVarExplainedMean[t] % piVarExplainedSD[t];
    }
    for (unsigned t=0; t<gamma.size(); ++t) {
        stringstream ss;
        ss << "MeanEffSq" << gamma[t] << "hsq";
        cout << boost::format("%-20s %-10.6f %-10.6f\n") % ss.str() % meanEffSqVarExplainedMean[t] % meanEffSqVarExplainedSD[t];
    }
    cout << boost::format("%-20s %-10s\n") % "Total number of SNPs" % numIncdSnps;
    //cout << "Computational time :  " << timer.format(timer.getElapse()) << endl << endl;
}

void Polygenicity::writeRes(const string &filename){
    string file = filename + ".Pi";
    ofstream out(file.c_str());
    if (!out) {
        throw("Error: cannot open file " + file);
    }
    
    VectorXf numSnpVarExplainedMean = numSnpVarExplained.mean();
    VectorXf numSnpVarExplainedSD = numSnpVarExplained.sd();
    VectorXf piVarExplainedMean = piVarExplained.mean();
    VectorXf piVarExplainedSD = piVarExplained.sd();
    VectorXf meanEffSqVarExplainedMean = meanEffSqVarExplained.mean();
    VectorXf meanEffSqVarExplainedSD = meanEffSqVarExplained.sd();

    out << boost::format("%-20s %-10.6f %-10.6f\n") % "Pi" % pi.mean() % pi.sd();
    out << boost::format("%-20s %-10.6f %-10.6f\n") % "NnzSnp" % nnz.mean()  % nnz.sd();
    for (unsigned t=0; t<gamma.size(); ++t) {
        stringstream ss;
        ss << "NumSnp" << gamma[t] << "hsq";
        out << boost::format("%-20s %-10.6f %-10.6f\n") % ss.str() % numSnpVarExplainedMean[t] % numSnpVarExplainedSD[t];
    }
    for (unsigned t=0; t<gamma.size(); ++t) {
        stringstream ss;
        ss << "Pi" << gamma[t] << "hsq";
        out << boost::format("%-20s %-10.6f %-10.6f\n") % ss.str() % piVarExplainedMean[t] % piVarExplainedSD[t];
    }
    for (unsigned t=0; t<gamma.size(); ++t) {
        stringstream ss;
        ss << "MeanEffSq" << gamma[t] << "hsq";
        out << boost::format("%-20s %-10.6f %-10.6f\n") % ss.str() % meanEffSqVarExplainedMean[t] % meanEffSqVarExplainedSD[t];
    }
    out << boost::format("%-20s %-10s\n") % "Total number of SNPs" % numIncdSnps;
    out.close();
}

void Polygenicity::writeMcmcSamples(const string &filename){
    string file1 = filename + ".Pi.mcmc";
    string file2 = filename + ".NumSnpVarExplained.mcmc";
    string file3 = filename + ".PiVarExplained.mcmc";
    ofstream out1(file1.c_str());
    ofstream out2(file2.c_str());
    ofstream out3(file3.c_str());
    if (!out1) throw("Error: cannot open file " + file1);
    if (!out2) throw("Error: cannot open file " + file2);
    if (!out3) throw("Error: cannot open file " + file3);
    out1 << "Pi" << endl << pi.datMat;
    for (unsigned i=0; i<numSnpVarExplained.ncol; ++i) out2 << "NumSNP" << gamma[i] << "hsq  ";
    out2 << endl << numSnpVarExplained.datMat;
    for (unsigned i=0; i<piVarExplained.ncol; ++i) out3 << "Pi" << gamma[i] << "hsq  ";
    out3 << endl << piVarExplained.datMat;
    out1.close();
    out2.close();
    out3.close();
}
