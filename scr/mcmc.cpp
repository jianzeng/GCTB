//
//  mcmc.cpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "mcmc.hpp"
#include "stat.hpp"
#include <Eigen/Core>
#include <cstdio>
#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#endif

static inline bool skipPosteriorTextOutput(const string &label) {
    return (label.find("Delta") != string::npos && label.find("AnnoDelta") != 0) ||
           label.find("PEP") != string::npos ||
           label.find("PIP") != string::npos;
}

void deleteMcmcSamples(vector<McmcSamples*> &mcmcSampleVec) {
    for (unsigned i = 0; i < mcmcSampleVec.size(); ++i) {
        delete mcmcSampleVec[i];
    }
    mcmcSampleVec.clear();
}

void McmcSamples::getParSample(const unsigned iter, const Parameter* par){
    if (iter % thin || iter <= burnin) return;
    if (multiChain) sampleIter.row(0) = par->perChainValue.transpose();
    else sampleIter(0,0) = par->value;
    accumulateIteration();
}

void McmcSamples::getParSetSample(const unsigned iter, const ParamSet* parSet){
    if (iter % thin || iter <= burnin) return;
    if (multiChain) sampleIter = parSet->perChainValues;
    else sampleIter.col(0) = parSet->values;
    accumulateIteration();
}

void McmcSamples::accumulateIteration(){
    ++cntPosteriorSample;
    const unsigned chains = jointStatistic ? 1 : numChains;
    const double t = cntPosteriorSample;
    const double count = t * chains;
    for (unsigned c=0; c<chains; ++c) {
        // Welford's centered moments, retaining chain identity for R-hat.
        const VectorXd delta = sampleIter.col(c).cast<double>() - perChainMean.col(c);
        perChainMean.col(c) += delta / t;
        perChainM2.col(c).array() += delta.array()
            * (sampleIter.col(c).cast<double>() - perChainMean.col(c)).array();
        thresholdCount.array() += (sampleIter.col(c).array() > criticalValue).cast<double>();
    }
    for (unsigned p=0; p<ncol; ++p) {
        const double mu = perChainMean.row(p).mean();
        const double m2 = perChainM2.row(p).sum()
            + t * (perChainMean.row(p).array() - mu).square().sum();
        posteriorMean[p] = mu;
        // Retain E[x^2] for expected squared effects / variance explained.
        posteriorSqrMean[p] = m2 / count + mu * mu;
        probGreaterThanCriticalValue[p] = jointStatistic
            ? std::numeric_limits<float>::quiet_NaN() : thresholdCount[p] / count;
    }
    if (storageMode == dense) {
        for (unsigned c=0; c<numChains; ++c)
            datMat.row((cntPosteriorSample-1)*numChains+c) = sampleIter.col(c).transpose();
    }
}

void McmcSamples::outputSample(const unsigned chain, ofstream &out){
    if (outputMode == txt) {
        tout << (burnin/thin + cntPosteriorSample)*thin << "\t" << chain+1;
        for (unsigned j=0; j<ncol; ++j) tout << "\t" << sampleIter(j,chain);
        tout << endl;
    }
    else if (outputMode == txt_combine_others) {
        out << "\t" << sampleIter.col(chain).transpose();
    }
    else if (outputMode == bin) {
        if (binarySampleRows.is_open())
            binarySampleRows << (cntPosteriorSample-1)*numChains + chain+1 << "\t"
                             << (burnin/thin + cntPosteriorSample)*thin << "\t" << chain+1 << endl;
        SparseVector<float> spvec = sampleIter.col(chain).sparseView();
        for (SparseVector<float>::InnerIterator it(spvec); it; ++it) {
            unsigned rc[2] = {(cntPosteriorSample-1)*numChains + chain, (unsigned)it.index()};
            fwrite(rc, sizeof(unsigned), 2, bout);
            float val = it.value();
            fwrite(&val, sizeof(float), 1, bout);
            ++nnz;
            //cout << it.index() << " " << it.value() << endl;
        }
        //cout << MatrixXf(spvec).transpose() << endl;
        //tout << sample.transpose() << endl;
    }
}

void McmcSamples::computeGelmanRubinStat(){
    if (!multiChain) return;
    GelmanRubinStat.setConstant(ncol, std::numeric_limits<float>::quiet_NaN());
    if (jointStatistic || cntPosteriorSample < 2) return;
    const double t = cntPosteriorSample;
    for (unsigned p=0; p<ncol; ++p) {
        const double mu = perChainMean.row(p).mean();
        const double W = perChainM2.row(p).mean() / (t-1);
        const double B = t * (perChainMean.row(p).array()-mu).square().sum() / (numChains-1);
        if (W > 0)
            GelmanRubinStat[p] = sqrt(((t-1)*W/t + B/t)/W);
        else
            GelmanRubinStat[p] = B > 0 ? std::numeric_limits<float>::infinity() : 1.0f;
    }
}

VectorXf McmcSamples::mean(){
    if (cntPosteriorSample) return perChainMean.rowwise().mean().cast<float>();
    if (storageMode == dense) {
        return VectorXf::Ones(nrow).transpose()*datMat/nrow;
    } else if (storageMode == sparse) {
        return VectorXf::Ones(nrow).transpose()*datMatSp/nrow;
    } else {
        return posteriorMean;
    }
}

VectorXf McmcSamples::sd(){
    VectorXf res(ncol);
    if (cntPosteriorSample) {
        const double t = cntPosteriorSample;
        const double count = t * perChainMean.cols();
        res.setConstant(std::numeric_limits<float>::quiet_NaN());
        if (jointStatistic || count < 2) return res;
        for (unsigned p=0; p<ncol; ++p) {
            const double mu = perChainMean.row(p).mean();
            const double m2 = perChainM2.row(p).sum()
                + t * (perChainMean.row(p).array()-mu).square().sum();
            res[p] = sqrt(std::max(0.0, m2) / (count-1));
        }
        return res;
    }
    if (storageMode == dense) {
        for (unsigned i=0; i<ncol; ++i) {
            res[i] = std::sqrt(Gadget::calcVariance(datMat.col(i)));
        }
    } else if (storageMode == sparse) {
        for (unsigned i=0; i<ncol; ++i) {
            res[i] = std::sqrt(Gadget::calcVariance(datMatSp.col(i)));
        }
    } else {
        res = (posteriorSqrMean.array() - posteriorMean.array().square()).sqrt();
    }
    return res;
}

void McmcSamples::initBinFile(const string &title){
    string dirname = title + ".mcmcsamples";
    if (!Gadget::directoryExist(dirname)) {
        throw("Error: cannot find directory " + dirname);
    }
    filename = dirname + "/" + label + ".mcmcsamples.bin";
    bout = fopen(filename.c_str(), "wb");
    if (!bout) {
        throw("Error: cannot open file " + filename);
    }
    nnz = 0;
    unsigned xyn[3] = {nrow, ncol, nnz};
    fwrite(xyn, sizeof(unsigned), 3, bout);
}

void McmcSamples::initTxtFile(const string &title){
    string dirname = title + ".mcmcsamples";
    if (!Gadget::directoryExist(dirname)) {
        throw("Error: cannot find directory " + dirname);
    }
    filename = dirname + "/" + label + ".mcmcsamples.txt";
    tout.open(filename.c_str());
    if (!tout) {
        throw("Error: cannot open file " + filename);
    }
}

void McmcSamples::writeDataBin(const string &title){
    string dirname = title + ".mcmcsamples";
    if (!Gadget::directoryExist(dirname)) {
        throw("Error: cannot find directory " + dirname);
    }
    filename = dirname + "/" + label + ".mcmcsamples.bin";
    FILE *out = fopen(filename.c_str(), "wb");
    if (!out) {
        throw("Error: cannot open file " + filename);
    }
    
    int xyn[3] = {static_cast<int>(datMatSp.rows()), static_cast<int>(datMatSp.cols()), static_cast<int>(datMatSp.nonZeros())};
    fwrite(xyn, sizeof(unsigned), 3, out);
    
    for (int i=0; i < datMatSp.outerSize(); ++i) {
        SpMat::InnerIterator it(datMatSp, i);
        for (; it; ++it) {
            unsigned rc[2] = {(unsigned)it.row(), (unsigned)it.col()};
            fwrite(rc, sizeof(unsigned), 2, out);
            float v = it.value();
            fwrite(&v, sizeof(float), 1, out);
        }
    }
    fclose(out);
}

void McmcSamples::readDataBin(const string &title){
    string dirname = title + ".mcmcsamples";
    if (!Gadget::directoryExist(dirname)) {
        throw("Error: cannot find directory " + dirname);
    }
    filename = dirname + "/" + label + ".mcmcsamples.bin";
    FILE *in = fopen(filename.c_str(), "rb");
    if (!in) {
        throw("Error: cannot open file " + filename);
    }
    
    unsigned xyn[3];
    fread(xyn, sizeof(unsigned), 3, in);
    
    nrow = xyn[0];
    ncol = xyn[1];

    datMatSp.resize(xyn[0], xyn[1]);

    if (xyn[2]) {
        typedef SpMat::StorageIndex SparseIndex;
        vector<SparseIndex> colCounts(ncol, 0);
        unsigned validNnz = 0;

        for (unsigned i = 0; i < xyn[2]; ++i) {
            unsigned rc[2];
            if (fread(rc, sizeof(unsigned), 2, in) != 2) break;
            float v;
            if (fread(&v, sizeof(float), 1, in) != 1) break;

            if(rc[0]>=xyn[0] || rc[1]>=xyn[1]) continue;

            ++colCounts[rc[1]];
            ++validNnz;
        }

        vector<SparseIndex> outer(ncol + 1, 0);
        for (unsigned col = 0; col < ncol; ++col)
            outer[col + 1] = outer[col] + colCounts[col];

        datMatSp.makeCompressed();
        datMatSp.resizeNonZeros(validNnz);
        std::copy(outer.begin(), outer.end(), datMatSp.outerIndexPtr());

        vector<SparseIndex> next = outer;
        fseek(in, sizeof(unsigned) * 3, SEEK_SET);
        for (unsigned i = 0; i < xyn[2]; ++i) {
            unsigned rc[2];
            if (fread(rc, sizeof(unsigned), 2, in) != 2) break;
            float v;
            if (fread(&v, sizeof(float), 1, in) != 1) break;

            if(rc[0]>=xyn[0] || rc[1]>=xyn[1]) continue;

            SparseIndex pos = next[rc[1]]++;
            datMatSp.innerIndexPtr()[pos] = static_cast<SparseIndex>(rc[0]);
            datMatSp.valuePtr()[pos] = v;
        }
        nnz = validNnz;
    } else {
        vector<Triplet<float> > trips;
        while (true) {
            unsigned rc[2];
            if (fread(rc, sizeof(unsigned), 2, in) != 2) break;
            float v;
            if (fread(&v, sizeof(float), 1, in) != 1) break;

            if(rc[0]>=xyn[0] || rc[1]>=xyn[1]) continue;

            trips.push_back(Triplet<float>(rc[0], rc[1], v));
        }
        nnz = trips.size();
        datMatSp.setFromTriplets(trips.begin(), trips.end());
        datMatSp.makeCompressed();
    }
    fclose(in);
    
    //cout << "nrow: " << nrow << " ncol: " << ncol << " nonzeros: " << datMatSp.nonZeros() << " " << nnz << endl;
    //cout << MatrixXf(datMatSp) << endl;
    
    storageMode = sparse;
}

void McmcSamples::closeOutputFile(void){
    if (outputMode == bin && bout) {
        unsigned xyn[3] = {nrow, ncol, nnz};
        fseek(bout, 0, SEEK_SET);
        fwrite(xyn, sizeof(unsigned), 3, bout);
        fclose(bout);
        bout = NULL;
    } else if (outputMode == txt && tout.is_open()) {
        tout.close();
    }
}

void McmcSamples::readDataTxt(const string &title){
    string dirname = title + ".mcmcsamples";
    if (!Gadget::directoryExist(dirname)) {
        throw("Error: cannot find directory " + dirname);
    }
    filename = dirname + "/" + label + ".mcmcsamples.txt";
    ifstream in(filename.c_str());
    string inputStr;
    vector<float> tmp;
    while (in >> inputStr) {
        tmp.push_back(stof(inputStr));
    }
    in.close();
    nrow = tmp.size();
    datMat.resize(nrow, 1);
    datMat.col(0) = Eigen::Map<VectorXf>(&tmp[0], nrow);
    storageMode = dense;
}

void McmcSamples::readDataTxt(const string &title, const string &label){
    string dirname = title + ".mcmcsamples";
    if (!Gadget::directoryExist(dirname)) {
        throw("Error: cannot find directory " + dirname);
    }
    filename = dirname + "/" + label + ".mcmcsamples.txt";
    ifstream in(filename.c_str());
    Gadget::Tokenizer colData;
    Gadget::Tokenizer header;
    string inputStr;
    string sep(" \t");
    vector<float> tmp;
    unsigned line = 0;
    
    std::getline(in, inputStr);
    header.getTokens(inputStr, sep);
    int idx = header.getIndex(label);
    
    if (idx==-1) throw("Error: Cannot find " + label + " in file [" + filename + "].");
    
    while (getline(in, inputStr)) {
        ++line;
        colData.getTokens(inputStr, sep);
        tmp.push_back(stof(colData[idx]));
    }
    in.close();
    nrow = tmp.size();
    datMat.resize(nrow, 1);
    datMat.col(0) = Eigen::Map<VectorXf>(&tmp[0], nrow);
    storageMode = dense;
}

void McmcSamples::writeDataTxt(const string &title){
    string dirname = title + ".mcmcsamples";
    if (!Gadget::directoryExist(dirname)) {
        throw("Error: cannot find directory " + dirname);
    }
    filename = dirname + "/" + label + ".mcmcsamples.txt";
    ofstream out(filename);
    out << datMat << endl;
    out.close();
}

void McmcSamples::writeMatSpTxt(const string &title){
    string dirname = title + ".mcmcsamples";
    if (!Gadget::directoryExist(dirname)) {
        throw("Error: cannot find directory " + dirname);
    }
    filename = dirname + "/" + label + ".mcmcsamples.txt";
    ofstream out(filename);
    out << boost::format("%8s %8s %12s\n") % "Row" % "Col" % "Value";
    
    for (int k=0; k<datMatSp.outerSize(); ++k) {
        for (SpMat::InnerIterator it(datMatSp,k); it; ++it) {
            out << boost::format("%8s %8s %12s\n")
            % it.row()
            % it.col()
            % it.value();
        }
    }
    out.close();
}

void MCMC::initTxtFile(const vector<Parameter*> &paramVec, const string &title){
    string dirname = title + ".mcmcsamples";
    if (!Gadget::directoryExist(dirname)) {
        throw("Error: cannot find directory " + dirname);
    }
    outfilename = dirname + "/CoreParameters.mcmcsamples.txt";
    out.open(outfilename.c_str());
    if (!out) {
        throw("Error: cannot open file " + outfilename);
    }
    out << "Iteration\tChain";
    for (unsigned i=0; i<paramVec.size(); ++i) {
        Parameter *par = paramVec[i];
        if (skipPosteriorTextOutput(par->label)) continue;
        out << "\t" << par->label;
    }
    out << endl;
}

vector<McmcSamples*> MCMC::initMcmcSamples(const Model &model, const unsigned numChains, const unsigned chainLength, const unsigned burnin, const unsigned thin,
                                           const string &title, const bool writeSnpMcmcBin, const bool writeParamMcmcTxt){
    vector<McmcSamples*> mcmcSampleVec;
    for (unsigned i=0; i<model.paramSetVec.size(); ++i) {
        ParamSet *parSet = model.paramSetVec[i];
        McmcSamples *mcmcSamples;
        if (parSet->label.find("SnpEffects") != string::npos) {
            if (writeSnpMcmcBin) {
                mcmcSamples = new McmcSamples(parSet->label, numChains, chainLength, burnin, thin, parSet->size, "do_not_store", "bin", title);
            } else {
                mcmcSamples = new McmcSamples(parSet->label, numChains, chainLength, burnin, thin, parSet->size, "do_not_store", "no_output", title);
            }
        } else {
            if (writeParamMcmcTxt && !skipPosteriorTextOutput(parSet->label)) {
                mcmcSamples = new McmcSamples(parSet->label, numChains, chainLength, burnin, thin, parSet->size, "do_not_store", "txt", title);
            } else {
                mcmcSamples = new McmcSamples(parSet->label, numChains, chainLength, burnin, thin, parSet->size, "do_not_store", "no_output", title);
            }
        }
        if (mcmcSamples->outputMode == McmcSamples::txt) {
            mcmcSamples->tout << "Iteration\tChain";
            for (const auto &name : parSet->header) mcmcSamples->tout << "\t" << name;
            mcmcSamples->tout << endl;
        }
        if (mcmcSamples->outputMode == McmcSamples::bin) {
            ofstream columns(mcmcSamples->filename + ".columns.txt");
            if (!columns) throw string("Error: cannot write binary sample column names.");
            columns << "Column\tName" << endl;
            for (unsigned j=0; j<parSet->header.size(); ++j)
                columns << j+1 << "\t" << parSet->header[j] << endl;
            mcmcSamples->binarySampleRows.open(mcmcSamples->filename + ".rows.txt");
            if (!mcmcSamples->binarySampleRows)
                throw string("Error: cannot write binary sample row names.");
            mcmcSamples->binarySampleRows << "Row\tIteration\tChain" << endl;
        }
        mcmcSampleVec.push_back(mcmcSamples);
        if (mcmcSamples->jointStatistic)
            cout << "Note: Marginal_Predictability_Enrichment is a joint-chain statistic; "
                 << "its mean is reported, but posterior SD, probability and R-hat are unavailable (NaN)." << endl;
    }
    for (unsigned i=0; i<model.paramVec.size(); ++i) {
        Parameter *par = model.paramVec[i];
        const string outputMode = (writeParamMcmcTxt && !skipPosteriorTextOutput(par->label)) ? "txt_combine_others" : "no_output";
        McmcSamples *mcmcSamples = new McmcSamples(par->label, numChains, chainLength, burnin, thin, 1, "do_not_store", outputMode, title);
        mcmcSampleVec.push_back(mcmcSamples);
    }
    if (writeParamMcmcTxt) initTxtFile(model.paramVec, title);
    return mcmcSampleVec;
}

void MCMC::collectSamples(const Model &model, vector<McmcSamples*> &mcmcSampleVec, const unsigned iteration){
    unsigned i=0;
    for (unsigned j=0; j<model.paramSetVec.size(); ++j) {
        McmcSamples *mcmcSamples = mcmcSampleVec[i++];
        ParamSet *parSet = model.paramSetVec[j];
        mcmcSamples->getParSetSample(iteration, parSet);
    }
    for (unsigned j=0; j<model.paramVec.size(); ++j) {
        McmcSamples *mcmcSamples = mcmcSampleVec[i++];
        Parameter *par = model.paramVec[j];
        mcmcSamples->getParSample(iteration, par);
    }
}

void MCMC::outputSamples(vector<McmcSamples*> &mcmcSampleVec, const unsigned numChains){
    for (unsigned chain=0; chain<numChains; ++chain) {
        if (out.is_open() && !mcmcSampleVec.empty()) {
            const McmcSamples *sample = mcmcSampleVec.front();
            out << (sample->burnin/sample->thin + sample->cntPosteriorSample)*sample->thin
                << "\t" << chain+1;
        }
        for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
            McmcSamples *mcmcSamples = mcmcSampleVec[i];
            if (skipPosteriorTextOutput(mcmcSamples->label)) continue;
            mcmcSamples->outputSample(chain, out);
        }
        if (out.is_open()) out << endl;
    }
}

void MCMC::computeGelmanRubinStat(vector<McmcSamples*> &mcmcSampleVec){
    for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
        McmcSamples *mcmcSamples = mcmcSampleVec[i];
        mcmcSamples->computeGelmanRubinStat();
    }
}

void MCMC::printStatus(const vector<Parameter*> &paramToPrint, const unsigned thisIter, const unsigned outputFreq, const string &timeLeft){
    if (thisIter==outputFreq) {
        cout << boost::format("%=10s ") % "Iter";
        for (unsigned i=0; i<paramToPrint.size(); ++i) {
            cout << boost::format("%=12s ") % paramToPrint[i]->label;
        }
        cout << boost::format("%=12s\n") % "TimeLeft";
    }
    cout << boost::format("%=10s ") % thisIter;
    for (unsigned i=0; i<paramToPrint.size(); ++i) {
        Parameter *par = paramToPrint[i];
        if (par->label[0] == 'N')
            cout << boost::format("%=12.0f ") % par->value;
        else
            cout << boost::format("%=12.4f ") % paramToPrint[i]->value;
    }
    cout << boost::format("%=12s\n") % timeLeft;
    
    cout.flush();
}



void MCMC::printSummary(const vector<Parameter*> &paramToPrint, const vector<McmcSamples*> &mcmcSampleVec, const unsigned numChains, const string &filename){
    if (!paramToPrint.size()) return;
    ofstream out;
    out.open(filename.c_str());
    if (!out) {
        throw("Error: cannot open file " + filename);
    }
    cout << "\nPosterior statistics from MCMC samples:\n\n";
    if (numChains > 1) {
        cout << boost::format("%10s %2s %-15s %-15s %-15s\n") %"Parameter" % "" % "Mean" % "SD " % "GelmanRubin_R";
        //out << "Posterior statistics from MCMC samples:\n\n";
        out << boost::format("%10s %2s %-15s %-15s %-15s\n") %"Parameter" % "" % "Mean" % "SD " % "GelmanRubin_R";

    } else {
        cout << boost::format("%10s %2s %-15s %-15s\n") %"Parameter" % "" % "Mean" % "SD ";
        //out << "Posterior statistics from MCMC samples:\n\n";
        out << boost::format("%10s %2s %-15s %-15s\n") %"Parameter" % "" % "Mean" % "SD ";
    }
    for (unsigned i=0; i<paramToPrint.size(); ++i) {
        Parameter *par = paramToPrint[i];
        for (unsigned j=0; j<mcmcSampleVec.size(); ++j) {
            McmcSamples *mcmcSamples = mcmcSampleVec[j];
            if (mcmcSamples->label == par->label) {
                if (mcmcSamples->numChains > 1) {
                    cout << boost::format("%10s %2s %-15.6f %-15.6f %-15.4f\n")
                    % par->label
                    % ""
                    % mcmcSamples->mean()
                    % mcmcSamples->sd()
                    % mcmcSamples->GelmanRubinStat;
                    out << boost::format("%10s %2s %-15.6f %-15.6f %-15.4f\n")
                    % par->label
                    % ""
                    % mcmcSamples->mean()
                    % mcmcSamples->sd()
                    % mcmcSamples->GelmanRubinStat;
                } else {
                    cout << boost::format("%10s %2s %-15.6f %-15.6f\n")
                    % par->label
                    % ""
                    % mcmcSamples->mean()
                    % mcmcSamples->sd();
                    out << boost::format("%10s %2s %-15.6f %-15.6f\n")
                    % par->label
                    % ""
                    % mcmcSamples->mean()
                    % mcmcSamples->sd();
                }
                break;
            }
        }
    }
    out.close();
}

void MCMC::printSetSummary(const vector<ParamSet*> &paramSetToPrint, const vector<McmcSamples*> &mcmcSampleVec, const unsigned numChains, const string &filename, const string &enrich){
    if (!paramSetToPrint.size()) return;
    ofstream out;
    out.open(filename.c_str());
    if (!out) {
        throw("Error: cannot open file " + filename);
    }
    ofstream out2;
    out2.open(enrich.c_str());
    if (!out2) {
        throw("Error: cannot open file " + enrich);
    }
    if (numChains > 1) {
        out << boost::format("%40s %40s %2s %-15s %-15s %-18s %-12s\n") % "Parameter" % "Annotation" % "" % "Mean" % "SD " % "PostProbAboveZero" % "GelmanRubin_R";
        out2 << boost::format("%40s %40s %2s %-15s %-15s %-18s %-12s\n") % "Parameter" % "Annotation" % "" % "Mean" % "SD " % "PostProbAboveOne" % "GelmanRubin_R";

    } else {
        out << boost::format("%40s %40s %2s %-15s %-15s %-18s\n") % "Parameter" % "Annotation" % "" % "Mean" % "SD " % "PostProbAboveZero";
        out2 << boost::format("%40s %40s %2s %-15s %-15s %-18s\n") % "Parameter" % "Annotation" % "" % "Mean" % "SD " % "PostProbAboveOne";
    }
    for (unsigned i=0; i<paramSetToPrint.size(); ++i) {
        ParamSet *parset = paramSetToPrint[i];
        if (parset->label == "SnpAnnoMembershipDelta") continue;
        for (unsigned j=0; j<mcmcSampleVec.size(); ++j) {
            McmcSamples *mcmcSamples = mcmcSampleVec[j];
            if (mcmcSamples->label == parset->label) {
                VectorXf mean = mcmcSamples->mean();
                VectorXf sd   = mcmcSamples->sd();
                Gadget::Tokenizer token;
                token.getTokens(parset->label, "_");
                if (token.back() == "Enrichment") {
                    for (unsigned col=0; col<parset->size; ++col) {
                        out2 << boost::format("%40s %40s %2s %-15.6f %-15.6f %-18.6f ")
                        % parset->label //token.front()
                        % parset->header[col]
                        % ""
                        % mean[col]
                        % sd[col]
                        % mcmcSamples->probGreaterThanCriticalValue[col];
                        if (mcmcSamples->numChains > 1) out2 << boost::format("%-12.4f ") % mcmcSamples->GelmanRubinStat[col];
                        out2 << endl;
                    }
                } else {
                    for (unsigned col=0; col<parset->size; ++col) {
                        out << boost::format("%40s %40s %2s %-15.6f %-15.6f %-18.6f ")
                        % parset->label
                        % parset->header[col]
                        % ""
                        % mean[col]
                        % sd[col]
                        % mcmcSamples->probGreaterThanCriticalValue[col];
                        if (mcmcSamples->numChains > 1) out << boost::format("%-12.4f ") % mcmcSamples->GelmanRubinStat[col];
                        out << endl;
                    }
                }
                break;
            }
        }
    }
    out.close();
    out2.close();
}

void MCMC::printSnpAnnoMembership(const vector<ParamSet *> &paramSetToPrint, const vector<McmcSamples *> &mcmcSampleVec, const string &filename) {
    if (!paramSetToPrint.size()) return;
    int idx = -9;
    for (unsigned i=0; i<paramSetToPrint.size(); ++i) {
        ParamSet *parset = paramSetToPrint[i];
        if (parset->label == "SnpAnnoMembershipDelta") idx = i;
    }
    if (idx == -9) return;
    ofstream out;
    out.open(filename.c_str());
    if (!out) {
        throw("Error: cannot open file " + filename);
    }
    ParamSet *parset = paramSetToPrint[idx];
    for (unsigned j=0; j<mcmcSampleVec.size(); ++j) {
        McmcSamples *mcmcSamples = mcmcSampleVec[j];
        if (mcmcSamples->label == parset->label) {
            for (unsigned col=0; col<parset->size; ++col) {
                out << boost::format("%-40s %-15.6f\n")
                % parset->header[col]
                % mcmcSamples->posteriorMean[col];
            }
            break;
        }
    }
}

vector<McmcSamples*> MCMC::run(Model &model, const unsigned numChains, const unsigned chainLength, const unsigned burnin, const unsigned thin, const bool print,
                               const unsigned outputFreq, const string &title, const bool writeSnpMcmcBin, const bool writeParamMcmcTxt,
                               const bool printFinalSummary, const unsigned rngChain){
    int savedEigenThreads = Eigen::nbThreads();
#ifdef _OPENMP
    if (omp_get_max_threads() > 1) {
        Eigen::setNbThreads(1);
    }
#endif
    Stat::logRngCheckpoint("MCMC_run_enter");
    if (print) {
        cout << "MCMC launched ..." << endl;
        cout << "  Number of chains: " << numChains << endl;
        cout << "  Chain length: " << chainLength << " iterations" << endl;
        cout << "  Burn-in: " << burnin << " iterations" << endl << endl;
    }
    
    if (writeSnpMcmcBin || writeParamMcmcTxt) {
        if (!Gadget::directoryExist(title + ".mcmcsamples")){
            Gadget::createDirectory(title + ".mcmcsamples");
            if (print) cout << "  Created directory [" << title << ".mcmcsamples] to store MCMC samples.\n\n";
        }
    }

    vector<McmcSamples*> mcmcSampleVec = initMcmcSamples(model, numChains, chainLength, burnin, thin, title, writeSnpMcmcBin, writeParamMcmcTxt);
    Stat::logRngCheckpoint("MCMC_after_init_samples");
    
    Gadget::Timer timer;
    timer.setTime();
    
    for (unsigned iteration=0; iteration<chainLength; ++iteration) {
        unsigned thisIter = iteration + 1;
        Stat::setMcmcContext(rngChain, thisIter);

        // Summaries are needed for kept samples (getParSetSample's rule) and the last iteration.
        model.summaryRequired = (thisIter > burnin && !(thisIter % thin))
                                || thisIter == chainLength;
        model.sampleUnknowns(thisIter);
        if (std::getenv("GCTB_DEBUG_RNG") &&
            (thisIter <= 5u || (thisIter % 10u) == 0u || thisIter == chainLength)) {
            char buf[96];
            std::snprintf(buf, sizeof(buf), "MCMC_after_sampleUnknowns_iter_%u", thisIter);
            Stat::logRngCheckpoint(buf);
        }
        
        collectSamples(model, mcmcSampleVec, thisIter);
        
        if ((writeSnpMcmcBin || writeParamMcmcTxt) && thisIter > burnin && !(thisIter % thin))
            outputSamples(mcmcSampleVec, numChains);
        if (numChains > 1 && thisIter == chainLength) computeGelmanRubinStat(mcmcSampleVec);
        
        setAction(model);
        if (action != keep_running) {
            cout << "\nMCMC sampling disrupted at iteration " << thisIter << endl;
            for (unsigned i = 0; i < mcmcSampleVec.size(); ++i)
                mcmcSampleVec[i]->closeOutputFile();
            Eigen::setNbThreads(savedEigenThreads);
            return mcmcSampleVec;
        }
        
        if (!(thisIter % outputFreq)) {
            timer.getTime();
            time_t timeToFinish = (chainLength-thisIter)*timer.getElapse()/thisIter; // remaining iterations multiplied by average time per iteration in seconds
            if (print) {
                printStatus(model.paramToPrint, thisIter, outputFreq, timer.format(timeToFinish));
            }
        }
    }
    
    // save the samples in the last iteration for potential continual run
    
    if (print) {
        cout << "\nMCMC cycles completed." << endl;
        if (printFinalSummary) {
            printSummary(model.paramToPrint, mcmcSampleVec, numChains, title + ".parRes");
            printSetSummary(model.paramSetToPrint, mcmcSampleVec, numChains, title + ".parSetRes", title + ".enrich");
            printSnpAnnoMembership(model.paramSetToPrint, mcmcSampleVec, title + ".snpAnnoMembership");
        }
    }

    Stat::logRngCheckpoint("MCMC_run_exit");

    for (unsigned i = 0; i < mcmcSampleVec.size(); ++i)
        mcmcSampleVec[i]->closeOutputFile();

    ///TMP
//    ofstream tmpOut;
//    tmpOut.open(title + ".varei");
//    tmpOut << static_cast<ApproxBayesS*>(&model)->vareiMean/(chainLength/100) << endl;

    Eigen::setNbThreads(savedEigenThreads);
    return mcmcSampleVec;
}

//vector<McmcSamples*> MCMC::run_multi_chains(Model model, const unsigned numChains, const unsigned chainLength, const unsigned burnin, const unsigned thin,
//                                            const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior){
//    if (myMPI::rank==0) cout << numChains << " ";
//
//    vector<vector<McmcSamples*> > mcmcSampleVecChain;
//    mcmcSampleVecChain.resize(numChains);
//
//}

void MCMC::setAction(const Model &model){
    if (model.status.empty()){
        action = keep_running;
    } else if (model.status == "Annealing") {
        action = keep_running;
    } else if (model.status == "Negative residual variance") {
        cout << "\nError: Residual variance is negative. This may indicate that effect sizes are \"blowing up\" likely due to a convergence problem. If SigmaSq variable is increasing with MCMC iterations, then this further indicates MCMC may not converge." << endl;
        action = restart_and_use_robust_model;
    } else if (model.status == "Unknown error") {
        action = stop_and_exit;
    }
}

void MCMC::convergeDiagGelmanRubin(const Model &model, vector<vector<McmcSamples *> > &mcmcSampleVecChain, const string &filename){
    if (!model.paramToPrint.size()) return;
    ofstream out;
    out.open((filename + ".parRes").c_str());
    cout << "\nPosterior statistics from multiple chains:\n\n";
    cout << boost::format("%13s %-15s %-15s %-12s\n") %"" % "Mean" % "SD " % "R_GelmanRubin ";
    //out << "Posterior statistics from multiple chains:\n\n";
    out << boost::format("%13s %-15s %-15s %-12s\n") %"" % "Mean" % "SD " % "R_GelmanRubin ";
    long numChains = mcmcSampleVecChain.size();
    VectorXf meanVec(numChains);
    VectorXf varVec(numChains);
    for (unsigned i=0; i<model.paramToPrint.size(); ++i) {
        Parameter *par = model.paramToPrint[i];
        for (unsigned j=0; j<mcmcSampleVecChain[0].size(); ++j) {
            McmcSamples *mcmcSamples = mcmcSampleVecChain[0][j];
            if (mcmcSamples->label == par->label) {
                float nsample = mcmcSamples->nrow;
                for (unsigned m=0; m<numChains; ++m) {
                    mcmcSamples = mcmcSampleVecChain[m][j];
                    meanVec[m] = mcmcSamples->mean()[0];
                    varVec[m]  = mcmcSamples->sd()[0];
                    varVec[m] *= varVec[m];
                }
                float posteriorMean = meanVec.mean();
                float B = (meanVec.array() - posteriorMean).matrix().squaredNorm()*nsample/float(numChains-1);
                float W = varVec.mean();
                float posteriorVar = (nsample-1.0)*W/nsample + B/nsample;
                float R = sqrt(posteriorVar/W);
                posteriorVar = ((nsample-1)*varVec.sum()
                    + nsample*(meanVec.array()-posteriorMean).square().sum())
                    / (nsample*numChains-1);
                
                cout << boost::format("%10s %2s %-15.6f %-15.6f %-12.3f\n")
                % par->label
                % ""
                % posteriorMean
                % sqrt(posteriorVar)
                % R;
                out << boost::format("%10s %2s %-15.6f %-15.6f %-12.3f\n")
                % par->label
                % ""
                % posteriorMean
                % sqrt(posteriorVar)
                % R;
                break;
            }
        }
    }
    out.close();
    
    if (model.paramSetToPrint.size()) {
        ofstream out2;
        out2.open((filename + ".parSetRes").c_str());
        
        for (unsigned i=0; i<model.paramSetToPrint.size(); ++i) {
            ParamSet *parset = model.paramSetToPrint[i];
            for (unsigned j=0; j<mcmcSampleVecChain[0].size(); ++j) {
                McmcSamples *mcmcSamples = mcmcSampleVecChain[0][j];
                if (mcmcSamples->label == parset->label) {
                    MatrixXf meanMat(numChains, parset->size);
                    MatrixXf varMat(numChains, parset->size);
                    float nsample = mcmcSamples->nrow;
                    for (unsigned m=0; m<numChains; ++m) {
                        mcmcSamples = mcmcSampleVecChain[m][j];
                        meanMat.row(m) = mcmcSamples->mean();
                        varMat.row(m)  = mcmcSamples->sd();
                        varMat.row(m) = varMat.row(m).array().square().matrix();
                    }
                    VectorXf posteriorMean = meanMat.colwise().mean();
                    VectorXf B = (meanMat.rowwise() - posteriorMean.transpose()).colwise().squaredNorm()*nsample/float(numChains-1);
                    VectorXf W = varMat.colwise().mean();
                    VectorXf posteriorVar = (nsample-1.0)*W/nsample + B/nsample;
                    VectorXf R = sqrt(posteriorVar.array()/W.array());
                    posteriorVar = ((nsample-1)*varMat.colwise().sum().transpose()
                        + nsample*(meanMat.rowwise()-posteriorMean.transpose()).colwise().squaredNorm().transpose())
                        / (nsample*numChains-1);
                    
                    for (unsigned col=0; col<parset->size; ++col) {
                        
                        out2 << boost::format("%25s %20s %2s %-15.6f %-15.6f ")
                        % parset->label
                        % parset->header[col]
                        % ""
                        % posteriorMean[col]
                        % sqrt(posteriorVar[col]);
                        Gadget::Tokenizer token;
                        token.getTokens(parset->label, "_");
                        float postprob = 0;
                        for (unsigned c=0; c<numChains; ++c)
                            postprob += mcmcSampleVecChain[c][j]->probGreaterThanCriticalValue[col];
                        postprob /= float(numChains);
                        out2 << boost::format("%-15.6f %-12.3f\n") % postprob % R[col];
                    }
                    break;
                }
            }
        }
        
    }

}
