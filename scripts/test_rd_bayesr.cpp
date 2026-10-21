#include "../scr/model.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    Stat::seedEngine(913);
    VectorXf weights(5); weights << .6f,.1f,.1f,.1f,.1f;
    // No observations: posterior assignment probabilities equal prior weights.
    auto p = ApproxBayesRD::AnnoEffects::mixtureProbabilities(
        VectorXf::Zero(2),VectorXf::Zero(2),weights);
    assert((p-weights).norm()<1e-6);
    // Analytic integrated likelihood for a null residual score, two levels.
    VectorXf d(2); d << 100,50;
    p = ApproxBayesRD::AnnoEffects::mixtureProbabilities(d,VectorXf::Zero(2),weights);
    VectorXf expected(5); expected[0]=.6;
    const double v[5]={0,.0001,.001,.01,.1};
    for(int g=1;g<5;++g) expected[g]=.1/std::sqrt((1+100*v[g])*(1+50*v[g]));
    expected/=expected.sum(); assert((p-expected).norm()<1e-6);
    p = ApproxBayesRD::AnnoEffects::mixtureProbabilities(d,VectorXf::Constant(2,10000),weights);
    assert(p.allFinite() && p[4]>.999f);

    for(unsigned annotations : {1u,5u,25u}) for(bool independent : {false,true}) {
        MatrixXf a(12,annotations); a.col(0).setOnes();
        std::vector<std::string> names(annotations,"annotation");
        for(unsigned j=1;j<annotations;++j) for(int r=0;r<12;++r)
            a(r,j)=((r+j)%3==0); // duplicate/collinear annotations
        ApproxBayesRD::AnnoEffects effects(names,4,a);
        VectorXf pis(4); pis << .7,.15,.1,.05; effects.initIntercept(pis);
        MatrixXf probabilities(12,3);
        for(unsigned it=0;it<100;++it) {
            VectorXi membership(12);
            for(int r=0;r<12;++r) membership[r]=it%2 ? r%4 : 0;
            Stat::setMcmcContext(0,it);
            if(independent) effects.sampleFromFC_indep(membership,a,VectorXf::Ones(3),
                VectorXf::Constant(3,.1),probabilities,true,true);
            else effects.sampleFromFC_joint(membership,a,VectorXf::Ones(3),.1,
                probabilities,true,true);
            assert(probabilities.allFinite());
            assert((probabilities.array()>=0).all() && (probabilities.array()<=1).all());
            assert(effects.values.allFinite());
            assert(effects.tgsProposals==0 && effects.swapProposals==0);
            for(unsigned c=0;c<3;++c) {
                auto alpha=effects[c]->values.tail(annotations-1);
                assert(std::abs(effects.ssq[c]-alpha.squaredNorm())<1e-5);
                assert(effects.numNonZeros[c]==(alpha.array()!=0).count());
                if(independent) assert(effects.mixtureCounts.row(c).sum()==annotations-1);
            }
            if(!independent) assert(effects.mixtureCounts.row(0).sum()==annotations-1);
            effects.updateMixtureWeights(independent);
            for(unsigned c=0;c<(independent ? 3u : 1u);++c) {
                auto w=effects.mixtureWeights[c]->values;
                assert(w.size()==5 && effects.mixtureCounts.cols()==5);
                assert(w.allFinite() && (w.array()>0).all());
                assert(std::abs(w.sum()-1)<1e-6);
            }
        }
    }
    // Prior-only draw moments verify both mixture allocation and Gaussian sampling.
    std::vector<std::string> names={"intercept","annotation"};
    ApproxBayesRD::AnnoEffects effects(names,2,MatrixXf::Ones(1,2));
    effects.mixtureWeights[0]->values=weights;
    double sumsq=0; unsigned counts[5]={};
    for(int it=0;it<30000;++it) {
        VectorXf draw; float pip;
        unsigned group=effects.drawMixture(VectorXf::Zero(1),VectorXf::Zero(1),0,draw,pip);
        ++counts[group]; sumsq+=draw.squaredNorm();
        assert(std::abs(pip-.4f)<1e-6);
        if(!group) assert(draw[0]==0);
    }
    for(int g=0;g<5;++g) assert(std::abs(counts[g]/30000.0-weights[g])<.015);
    assert(std::abs(sumsq/30000.0-.01111)<.001);
    std::cout << "RD BayesR conditional, integration and prior-moment checks passed\n";
}
