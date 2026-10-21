#include "../scr/model.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

void checkIntegration(unsigned annotations, bool full) {
    const unsigned n=12, components=3;
    MatrixXf a(n,annotations);
    a.col(0).setOnes();
    std::vector<std::string> names(annotations);
    for(unsigned j=0;j<annotations;++j) {
        names[j]="A"+std::to_string(j);
        if(j) for(unsigned r=0;r<n;++r)
            a(r,j)=((r+j)%3==0) ? 1.0f : 0.0f;
    }
    // Deliberately duplicated columns and more columns than observations.
    VectorXf pis(4); pis << 0.7,0.15,0.1,0.05;
    VectorXf variance=VectorXf::Constant(components,0.4);
    ApproxBayesRC::AnnoEffects rc(names,4,a);
    assert(!rc.blockedGaussian);
    rc.blockedGaussian = full;
    rc.initIntercept_probit(pis);
    ApproxBayesRD::AnnoEffects indep(names,4,a), joint(names,4,a);
    indep.bayesR = joint.bayesR = false; // regression coverage for the legacy kernels
    indep.initIntercept(pis); joint.initIntercept(pis);
    indep.initializeSwapNeighbors(a.transpose()*a);
    joint.initializeSwapNeighbors(a.transpose()*a);
    VectorXi membership(n);
    MatrixXf probabilities(n,components);
    for(unsigned iteration=0;iteration<30;++iteration) {
        Stat::setMcmcContext(0,iteration);
        for(unsigned r=0;r<n;++r)
            membership[r]=iteration%3 ? (r+iteration)%4 : 0;
        rc.sampleFromFC_Gibbs(membership,a,variance,probabilities);
        assert(probabilities.allFinite());
        assert((probabilities.array()>=0).all() && (probabilities.array()<=1).all());
        for(unsigned c=0;c<components;++c) {
            assert(rc[c]->values.allFinite());
            assert(std::abs(rc.ssq[c]-rc[c]->values.tail(annotations-1).squaredNorm())<1e-4);
        }
        indep.sampleFromFC_indep(membership,a,variance,
            VectorXf::Constant(components,0.3),probabilities,false,full);
        joint.sampleFromFC_joint(membership,a,variance,0.3,probabilities,false,full);
        for(ApproxBayesRD::AnnoEffects *effects : {&indep,&joint}) {
            assert(effects->values.allFinite());
            for(unsigned c=0;c<components;++c) {
                const VectorXf &alpha=(*effects)[c]->values;
                assert(std::abs(effects->ssq[c]-alpha.tail(annotations-1).squaredNorm())<1e-4);
                assert(effects->numNonZeros[c]==(alpha.tail(annotations-1).array()!=0).count());
                assert((effects->values.col(c)-alpha).norm()==0);
            }
        }
        assert(probabilities.allFinite());
    }
    if(full && annotations>1) {
        assert(indep.tgsProposals>0 && joint.tgsProposals>0);
        assert(indep.tgsAccepted>0 && joint.tgsAccepted>0);
    } else {
        assert(indep.tgsProposals==0 && joint.tgsProposals==0);
    }
    if(full) {
        assert(rc.gaussianGramCache.size()==1);
        assert(rc.gaussianGramCache[0].rows()==annotations);
        if(annotations>1) {
            assert(indep.tgsGramCache.size()==1 && joint.tgsGramCache.size()==1);
            assert(indep.tgsGramCache[0].rows()==annotations-1);
            assert(joint.tgsGramCache[0].rows()==annotations-1);
        }
    } else assert(rc.gaussianGramCache.empty());
    std::cout << "Integration passed: " << annotations << " annotations, full=" << full << '\n';
}

int main() {
    Stat::seedEngine(12345);
    for(bool full : {false,true}) {
        checkIntegration(1,full); checkIntegration(4,full); checkIntegration(35,full);
    }
}
