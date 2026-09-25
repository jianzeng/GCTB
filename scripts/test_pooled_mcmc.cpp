#include "../scr/mcmc.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    Parameter p("test");
    p.perChainValue.resize(2);
    McmcSamples streaming("test",2,2,0,1,1,"do_not_store","no_output","");
    McmcSamples dense("test",2,2,0,1,1,"dense","no_output","");
    for (unsigned t=1;t<=2;++t) {
        p.perChainValue << (t==1 ? 0 : 2), (t==1 ? 4 : 6);
        p.value=p.perChainValue.mean();
        streaming.getParSample(t,&p);
        dense.getParSample(t,&p);
    }
    assert(std::abs(streaming.mean()[0]-3)<1e-6);
    assert(std::abs(streaming.sd()[0]-std::sqrt(20.0/3))<1e-6);
    assert(dense.datMat(0,0)==0 && dense.datMat(1,0)==4);
    assert(dense.datMat(2,0)==2 && dense.datMat(3,0)==6);
    assert(dense.sd()[0]==streaming.sd()[0]);
    streaming.computeGelmanRubinStat();
    assert(std::abs(streaming.GelmanRubinStat[0]-std::sqrt(4.5))<1e-6);

    std::vector<std::string> names(1,"annotation");
    ParamSet a("Test_Enrichment",names);
    a.perChainValues.resize(1,2);
    McmcSamples enrich(a.label,2,2,0,1,1,"do_not_store","no_output","");
    for(unsigned t=1;t<=2;++t) {
        a.perChainValues << 0,4;
        a.values=a.perChainValues.rowwise().mean();
        enrich.getParSetSample(t,&a);
    }
    assert(enrich.probGreaterThanCriticalValue[0]==0.5f);
    enrich.computeGelmanRubinStat();
    assert(std::isinf(enrich.GelmanRubinStat[0]));

    McmcSamples stable("large",1,3,0,1,1,"do_not_store","no_output","");
    for(unsigned t=1;t<=3;++t) {
        p.value=10000000.0f+t;
        stable.getParSample(t,&p);
    }
    assert(stable.sd()[0]==1);
    McmcSamples joint("Marginal_Predictability_Enrichment",2,2,0,1,1,"do_not_store","no_output","");
    a.perChainValues.setConstant(2);
    a.values.setConstant(2);
    joint.getParSetSample(1,&a);
    joint.getParSetSample(2,&a);
    joint.computeGelmanRubinStat();
    assert(std::isnan(joint.sd()[0]));
    assert(std::isnan(joint.probGreaterThanCriticalValue[0]));
    assert(std::isnan(joint.GelmanRubinStat[0]));
    assert(joint.mean()[0]==2);
    std::cout << "Pooled moments, dense storage, probabilities, R-hat and stability passed.\n";
}
