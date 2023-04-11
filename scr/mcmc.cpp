//
//  mcmc.cpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "mcmc.hpp"

void McmcSamples::getSample(const unsigned iter, const VectorXf &sample, const bool writeBinPosterior, const bool writeTxtPosterior){
    if (storageMode == dense) {
        if (writeTxtPosterior) tout << sample.transpose() << endl;
    }
    if (iter % thin) return;
//    if (!sample.size()) return;
    unsigned thin_iter = iter/thin;
    unsigned thin_iter_post_burnin = thin_iter - burnin/thin;
    if (storageMode == dense) {
//        if (writeTxtPosterior) tout << sample.transpose() << endl;
        if (iter >= burnin) {
            datMat.row(thin_iter_post_burnin) = sample;
            posteriorMean.array() += (sample - posteriorMean).array()/(thin_iter_post_burnin+1);
            posteriorSqrMean.array() += (sample.array().square() - posteriorSqrMean.array())/(thin_iter_post_burnin+1);
        }
    } else if (storageMode == sparse) {
        lastSample = sample;
        //SparseVector<float>::InnerIterator it(sample.sparseView());
        SparseVector<float> spvec = sample.sparseView();
        if (writeBinPosterior) {
            for (SparseVector<float>::InnerIterator it(spvec); it; ++it) {
                unsigned rc[2] = {thin_iter, (unsigned)it.index()};
                fwrite(rc, sizeof(unsigned), 2, bout);
                float val = it.value();
                fwrite(&val, sizeof(float), 1, bout);
                //cout << it.index() << " " << it.value() << endl;
            }
            //cout << MatrixXf(spvec).transpose() << endl;
//            tout << sample.transpose() << endl;
        }
        //nnz += sample.sparseView().nonZeros();
        //cout << nnz << " " << sample.sparseView().nonZeros() <<endl;
        ArrayXf delta;
        delta.setZero(sample.size());
        for (SparseVector<float>::InnerIterator it(spvec); it; ++it) {
            delta[it.index()] = 1;
        }
        if (iter >= burnin) {
            pip.array() += (delta - pip.array())/(thin_iter_post_burnin+1);
            posteriorMean.array() += (sample - posteriorMean).array()/(thin_iter_post_burnin+1);
            posteriorSqrMean.array() += (sample.array().square() - posteriorSqrMean.array())/(thin_iter_post_burnin+1);            
        }
    }
}

void McmcSamples::getSample(const unsigned iter, const float sample, const bool writeTxtPosterior, ofstream &out){
    if (writeTxtPosterior) out << boost::format("%12s ") %sample;
    if (iter % thin) return;
    unsigned thin_iter_post_burnin = iter/thin - burnin/thin;
    if (iter >= burnin) {
        datMat(thin_iter_post_burnin,0) = sample;
        posteriorMean.array() += (sample - posteriorMean.array())/(thin_iter_post_burnin+1);
        posteriorSqrMean.array() += (sample*sample - posteriorSqrMean.array())/(thin_iter_post_burnin+1);
    }
}

VectorXf McmcSamples::mean(){
    if (storageMode == dense) {
        return VectorXf::Ones(nrow).transpose()*datMat/nrow;
    } else {
        return VectorXf::Ones(nrow).transpose()*datMatSp/nrow;
    }
}

VectorXf McmcSamples::sd(){
    VectorXf res(ncol);
    if (storageMode == dense) {
        for (unsigned i=0; i<ncol; ++i) {
            res[i] = std::sqrt(Gadget::calcVariance(datMat.col(i)));
        }
    } else {
        for (unsigned i=0; i<ncol; ++i) {
            res[i] = std::sqrt(Gadget::calcVariance(datMatSp.col(i)));
        }
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
    unsigned xyn[3] = {chainLength/thin, ncol, nnz};
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
//    vector<Triplet<float>> trips(xyn[2]);
    vector<Triplet<float> > trips;
    
    //for (int i=0; i < trips.size(); ++i){
    while (!feof(in)) {
        unsigned rc[2];
        fread(rc, sizeof(unsigned), 2, in);
        float v;
        fread(&v, sizeof(float), 1, in);
        
        if(rc[0]>xyn[0] || rc[1]>xyn[1]) continue;
        
        //trips[i] = Triplet<float>(rc[0], rc[1], v);
        trips.push_back(Triplet<float>(rc[0], rc[1], v));
    }
    fclose(in);
    
    datMatSp.setFromTriplets(trips.begin(), trips.end());
    datMatSp.makeCompressed();
    
    //cout << "nrow: " << nrow << " ncol: " << ncol << " nonzeros: " << datMatSp.nonZeros() << " " << nnz << endl;
    //cout << MatrixXf(datMatSp) << endl;
    
    storageMode = sparse;
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
    for (unsigned i=0; i<paramVec.size(); ++i) {
        Parameter *par = paramVec[i];
        out << boost::format("%12s ") %par->label;
    }
    out << endl;
}

vector<McmcSamples*> MCMC::initMcmcSamples(const Model &model, const unsigned chainLength, const unsigned burnin, const unsigned thin,
                                           const string &title, const bool writeBinPosterior, const bool writeTxtPosterior){
    vector<McmcSamples*> mcmcSampleVec;
    for (unsigned i=0; i<model.paramSetVec.size(); ++i) {
        ParamSet *parSet = model.paramSetVec[i];
        McmcSamples *mcmcSamples;
        if (parSet->label.find("SnpEffects") != string::npos) {
            mcmcSamples = new McmcSamples(parSet->label, chainLength, burnin, thin, parSet->size, "sparse");
            if (writeBinPosterior) mcmcSamples->initBinFile(title);
//            mcmcSamples->initTxtFile(title);
        } else if (parSet->label.find("Delta") != string::npos) {
            mcmcSamples = new McmcSamples(parSet->label, chainLength, burnin, thin, parSet->size, "sparse");
            if (writeBinPosterior) mcmcSamples->initBinFile(title);
        } else {
            mcmcSamples = new McmcSamples(parSet->label, chainLength, burnin, thin, parSet->size);
            if (writeTxtPosterior) mcmcSamples->initTxtFile(title);
        }
        mcmcSampleVec.push_back(mcmcSamples);
    }
    for (unsigned i=0; i<model.paramVec.size(); ++i) {
        Parameter *par = model.paramVec[i];
        McmcSamples *mcmcSamples = new McmcSamples(par->label, chainLength, burnin, thin, 1);
        //mcmcSamples->initTxtFile(title);
        mcmcSampleVec.push_back(mcmcSamples);
    }
    if (writeTxtPosterior) initTxtFile(model.paramVec, title);
    return mcmcSampleVec;
}

void MCMC::collectSamples(const Model &model, vector<McmcSamples*> &mcmcSampleVec, const unsigned iteration, const bool writeBinPosterior, const bool writeTxtPosterior){
    unsigned i = 0;
    for (unsigned j=0; j<model.paramSetVec.size(); ++j) {
        McmcSamples *mcmcSamples = mcmcSampleVec[i++];
        ParamSet *parSet = model.paramSetVec[j];
        mcmcSamples->getSample(iteration, parSet->values, writeBinPosterior, writeTxtPosterior);
    }
    for (unsigned j=0; j<model.paramVec.size(); ++j) {
        McmcSamples *mcmcSamples = mcmcSampleVec[i++];
        Parameter *par = model.paramVec[j];
        mcmcSamples->getSample(iteration, par->value, writeTxtPosterior, out);
    }
    out << endl;
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



void MCMC::printSummary(const vector<Parameter*> &paramToPrint, const vector<McmcSamples*> &mcmcSampleVec, const string &filename){
    if (!paramToPrint.size()) return;
    ofstream out;
    out.open(filename.c_str());
    if (!out) {
        throw("Error: cannot open file " + filename);
    }
    cout << "\nPosterior statistics from MCMC samples:\n\n";
    cout << boost::format("%13s %-15s %-15s\n") %"" % "Mean" % "SD ";
    //out << "Posterior statistics from MCMC samples:\n\n";
    out << boost::format("%13s %-15s %-15s\n") %"" % "Mean" % "SD ";
    for (unsigned i=0; i<paramToPrint.size(); ++i) {
        Parameter *par = paramToPrint[i];
        for (unsigned j=0; j<mcmcSampleVec.size(); ++j) {
            McmcSamples *mcmcSamples = mcmcSampleVec[j];
            if (mcmcSamples->label == par->label) {
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
                break;
            }
        }
    }
    out.close();
}

void MCMC::printSetSummary(const vector<ParamSet*> &paramSetToPrint, const vector<McmcSamples*> &mcmcSampleVec, const string &filename){
    if (!paramSetToPrint.size()) return;
    ofstream out;
    out.open(filename.c_str());
    if (!out) {
        throw("Error: cannot open file " + filename);
    }
//    cout << "\nPosterior statistics from MCMC samples:\n\n";
//    cout << boost::format("%13s %-15s %-15s\n") %"" % "Mean" % "SD ";
//    out << "Posterior statistics from MCMC samples:\n\n";
//    out << boost::format("%13s %-15s %-15s\n") %"" % "Mean" % "SD ";
    for (unsigned i=0; i<paramSetToPrint.size(); ++i) {
        ParamSet *parset = paramSetToPrint[i];
        if (parset->label == "SnpAnnoMembershipDelta") continue;
        for (unsigned j=0; j<mcmcSampleVec.size(); ++j) {
            McmcSamples *mcmcSamples = mcmcSampleVec[j];
            if (mcmcSamples->label == parset->label) {
                for (unsigned col=0; col<parset->size; ++col) {
//                    cout << boost::format("%20s %10s %2s %-15.6f %-15.6f\n")
//                    % parset->label
//                    % parset->header[col]
//                    % ""
//                    % mcmcSamples->posteriorMean[col]
//                    % sqrt(mcmcSamples->posteriorSqrMean[col]-mcmcSamples->posteriorMean[col]*mcmcSamples->posteriorMean[col]);
                    out << boost::format("%25s %20s %2s %-15.6f %-15.6f ")
                    % parset->label
                    % parset->header[col]
                    % ""
                    % mcmcSamples->posteriorMean[col]
                    % sqrt(mcmcSamples->posteriorSqrMean[col]-mcmcSamples->posteriorMean[col]*mcmcSamples->posteriorMean[col]);
                    Gadget::Tokenizer token;
                    token.getTokens(parset->label, "_");
                    float postprob = 0;
                    if (token.back() == "Enrichment") {
                        for (unsigned row=0; row<mcmcSamples->nrow; ++row) {
                            if (mcmcSamples->datMat(row, col) > 1) ++postprob;
                        }
                    } else {
                        for (unsigned row=0; row<mcmcSamples->nrow; ++row) {
                            if (mcmcSamples->datMat(row, col) > 0) ++postprob;
                        }
                    }
                    postprob /= float(mcmcSamples->nrow);
                    out << boost::format("%-15.6f\n") % postprob;
                }
                break;
            }
        }
    }
    out.close();
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

vector<McmcSamples*> MCMC::run(Model &model, const unsigned chainLength, const unsigned burnin, const unsigned thin, const bool print,
                               const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior){
    if (print) {
        cout << "MCMC launched ..." << endl;
        cout << "  Chain length: " << chainLength << " iterations" << endl;
        cout << "  Burn-in: " << burnin << " iterations" << endl << endl;
    }
    
    if (writeBinPosterior || writeTxtPosterior) {
        if (!Gadget::directoryExist(title + ".mcmcsamples")){
            Gadget::createDirectory(title + ".mcmcsamples");
            if (print) cout << "  Created directory [" << title << ".mcmcsamples] to store MCMC samples.\n\n";
        }
    }

    vector<McmcSamples*> mcmcSampleVec = initMcmcSamples(model, chainLength, burnin, thin, title, writeBinPosterior, writeTxtPosterior);
    
    Gadget::Timer timer;
    timer.setTime();
    
    for (unsigned iteration=0; iteration<chainLength; ++iteration) {
        unsigned thisIter = iteration + 1;
        
        model.sampleUnknowns();
        collectSamples(model, mcmcSampleVec, iteration, writeBinPosterior, writeTxtPosterior);
        
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
        printSummary(model.paramToPrint, mcmcSampleVec, title + ".parRes");
        printSetSummary(model.paramSetToPrint, mcmcSampleVec, title + ".parSetRes");
        printSnpAnnoMembership(model.paramSetToPrint, mcmcSampleVec, title + ".snpAnnoMembership");
    }

    ///TMP
//    ofstream tmpOut;
//    tmpOut.open(title + ".varei");
//    tmpOut << static_cast<ApproxBayesS*>(&model)->vareiMean/(chainLength/100) << endl;

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
                        varMat.row(m) *= varMat.row(m);
                    }
                    VectorXf posteriorMean = meanMat.colwise().mean();
                    VectorXf B = (meanMat.rowwise() - posteriorMean.transpose()).colwise().squaredNorm()*nsample/float(numChains-1);
                    VectorXf W = varMat.colwise().mean();
                    VectorXf posteriorVar = (nsample-1.0)*W/nsample + B/nsample;
                    VectorXf R = sqrt(posteriorVar.array()/W.array());
                    
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
                        if (token.back() == "Enrichment") {
                            for (unsigned row=0; row<mcmcSamples->nrow; ++row) {
                                if (mcmcSamples->datMat(row, col) > 1) ++postprob;
                            }
                        } else {
                            for (unsigned row=0; row<mcmcSamples->nrow; ++row) {
                                if (mcmcSamples->datMat(row, col) > 0) ++postprob;
                            }
                        }
                        postprob /= float(mcmcSamples->nrow);
                        out2 << boost::format("%-15.6f %-12.3f\n") % postprob % R[col];
                    }
                    break;
                }
            }
        }
        
    }

}
