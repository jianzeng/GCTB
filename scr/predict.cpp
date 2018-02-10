//
//  predict.cpp
//  gctb
//
//  Created by Jian Zeng on 20/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "predict.hpp"

void Predict::getAccuracy(const Data &data, const string &filename) {
    ofstream out(filename.c_str());
    if (!out) {
        throw("Error: cannot open file " + filename);
    }
    ghat.setZero(data.numKeptInds);
    SnpInfo *snp;
    for (unsigned i=0; i<data.numIncdSnps; ++i) {
        snp = data.incdSnpInfoVec[i];
        if (snp->effect)
            ghat += data.Z.col(i) * snp->effect;
    }
    cor = Gadget::calcCorrelation(data.y, ghat);
    reg = Gadget::calcRegression(data.y, ghat);
    
    cout << " Accuracy of prediction (correlation between y and ghat) : " << boost::format("%-10.4f") % cor << endl;
    cout << " Bias of prediction (one minus regression of y on  ghat) : " << boost::format("%-10.4f") % (1.0f-reg) << endl;
    
    out << " Accuracy of prediction (correlation between y and ghat) : " << boost::format("%-10.4f") % cor << endl;
    out << " Bias of prediction (one minus regression of y on  ghat) : " << boost::format("%-10.4f") % (1.0f-reg) << endl;
}

void Predict::writeRes(const Data &data, const string &filename) {
    ofstream out(filename.c_str());
    if (!out) {
        throw("Error: cannot open file " + filename);
    }
    IndInfo *ind;
    out << boost::format("%20s %20s %20s %20s\n")
    % "FID" % "PID" % "Phenotype" % "Predicted";
    for (unsigned i=0; i<data.numKeptInds; ++i) {
        ind = data.keptIndInfoVec[i];
        out << boost::format("%20s %20s %20s %20s\n")
        % ind->famID
        % ind->indID
        % ind->phenotype
        % ghat[i];
    }
    out.close();
}
