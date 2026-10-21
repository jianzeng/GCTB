#include "../scr/options.hpp"
#include <cassert>
#include <iostream>

int main(int argc, char **argv) {
    Options options;
    assert(!options.jointAnno);
    options.validateJointAnno();
    options.jointAnno=true;
    options.analysisType="SBayes";
    options.annotationFile="annotations.txt";
    options.estimatePi=true;
    for (const std::string &model : {"RC","RD","RD_indep"}) {
        options.bayesType=model;
        options.validateJointAnno();
    }
    options.bayesType="R";
    bool rejected=false;
    try { options.validateJointAnno(); } catch (const std::string &) { rejected=true; }
    assert(rejected);
    Options parsed;
    const char *args[] = {"gctb", "--sbayes", "RC", "--annot", "annotations.txt", "--joint-annot"};
    parsed.inputOptions(6,args);
    assert(parsed.jointAnno);
    assert(parsed.optionsSummary.find("--joint-annot")!=std::string::npos);
    Options defaults;
    defaults.inputOptions(5,args);
    assert(!defaults.jointAnno);
    if(argc>1) {
        Options fileOptions;
        const char *fileArgs[] = {"gctb", "--inp-file", argv[1]};
        fileOptions.inputOptions(3,fileArgs);
        assert(fileOptions.jointAnno);
    }
    options.bayesType="RC";
    options.estimatePi=false;
    rejected=false;
    try { options.validateJointAnno(); } catch (const std::string &) { rejected=true; }
    assert(rejected);
    std::cout << "Joint annotation defaults and validation passed.\n";
}
