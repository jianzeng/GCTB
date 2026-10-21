#include "../scr/model.hpp"
#include <cassert>
#include <iostream>

using namespace Eigen;

void gaussianTest(int blockSize) {
    MatrixXf a(6,3);
    a << 1,-1,-0.9, 1,-0.5,-0.6, 1,0,0.2,
         1,0.4,0.5, 1,0.8,0.7, 1,1,1.1;
    VectorXf y(6); y << -1,0,0.5,0.8,1,2;
    const double priorVariance=0.7;
    MatrixXd precision=a.cast<double>().transpose()*a.cast<double>();
    precision(1,1)+=1/priorVariance; precision(2,2)+=1/priorVariance;
    MatrixXd covariance=precision.inverse();
    VectorXd truth=covariance*a.cast<double>().transpose()*y.cast<double>();
    VectorXf alpha=VectorXf::Zero(3), residual=y;
    VectorXd sum=VectorXd::Zero(3);
    MatrixXd sumOuter=MatrixXd::Zero(3,3);
    AnnotationSampler::GramCache cache;
    const int draws=80000, burn=1000;
    for(int t=0;t<draws+burn;++t) {
        AnnotationSampler::gaussian(a,residual,alpha,priorVariance,&cache,blockSize);
        assert((residual-(y-a*alpha)).norm()<0.001);
        if(t>=burn) {
            VectorXd x=alpha.cast<double>(); sum+=x; sumOuter+=x*x.transpose();
        }
        // Production refreshes residuals with every new latent response.
        residual=y-a*alpha;
    }
    VectorXd mean=sum/draws;
    MatrixXd cov=sumOuter/draws-mean*mean.transpose();
    assert((mean-truth).cwiseAbs().maxCoeff()<0.02);
    assert((cov-covariance).cwiseAbs().maxCoeff()<0.025);
    std::cout << "Gaussian block " << blockSize << ": mean/covariance correct\n";
}

// Enumerate all four inclusion configurations and integrate Gaussian slabs.
// This checks the TGS kernel alone, without a Gibbs sweep concealing bias.
void tgsTest(int components, double pi) {
    std::vector<MatrixXf> a(components);
    std::vector<VectorXf> y(components), residual(components);
    std::vector<const MatrixXf*> designs;
    std::vector<VectorXf*> residuals;
    VectorXf variances(components);
    for(int c=0;c<components;++c) {
        a[c].resize(6-c,3);
        for(int r=0;r<a[c].rows();++r) {
            a[c](r,0)=1;
            a[c](r,1)=(r-2.5)*0.6;
            a[c](r,2)=a[c](r,1)+(r%2 ? 0.2f : -0.2f);
        }
        y[c]=0.7*a[c].col(1)-0.2*a[c].col(2);
        y[c][0]+=0.4;
        residual[c]=y[c]; variances[c]=0.7+0.3*c;
        designs.push_back(&a[c]); residuals.push_back(&residual[c]);
    }
    VectorXd logMass(4);
    std::vector<MatrixXd> conditionalMean(4), conditionalSecond(4);
    for(int mask=0;mask<4;++mask) {
        std::vector<int> active;
        for(int j=0;j<2;++j) if(mask&(1<<j)) active.push_back(j+1);
        int k=active.size();
        logMass[mask]=k*std::log(pi)+(2-k)*std::log1p(-pi);
        conditionalMean[mask]=MatrixXd::Zero(2,components);
        conditionalSecond[mask]=MatrixXd::Zero(2,components);
        if(!k) continue;
        for(int c=0;c<components;++c) {
            MatrixXd design(a[c].rows(),k);
            for(int j=0;j<k;++j) design.col(j)=a[c].col(active[j]).cast<double>();
            MatrixXd h=design.transpose()*design;
            h.diagonal().array()+=1.0/variances[c];
            VectorXd rhs=design.transpose()*y[c].cast<double>();
            MatrixXd v=h.inverse(); VectorXd m=v*rhs;
            logMass[mask]+=0.5*(rhs.dot(m)-std::log(h.determinant())-k*std::log(variances[c]));
            for(int j=0;j<k;++j) {
                conditionalMean[mask](active[j]-1,c)=m[j];
                conditionalSecond[mask](active[j]-1,c)=m[j]*m[j]+v(j,j);
            }
        }
    }
    VectorXd mass=(logMass.array()-AnnotationSampler::logSumExp(logMass)).exp();
    VectorXd expectedPip=VectorXd::Zero(2), observedPip=VectorXd::Zero(2);
    MatrixXd expectedMean=MatrixXd::Zero(2,components), expectedSecond=expectedMean;
    MatrixXd sum=expectedMean, squares=expectedMean;
    for(int mask=0;mask<4;++mask) {
        for(int j=0;j<2;++j) if(mask&(1<<j)) expectedPip[j]+=mass[mask];
        expectedMean+=mass[mask]*conditionalMean[mask];
        expectedSecond+=mass[mask]*conditionalSecond[mask];
    }
    MatrixXf alpha=MatrixXf::Zero(3,components);
    AnnotationSampler::GramCache cache;
    unsigned long long proposals=0, accepted=0;
    const int draws=150000, burn=2000;
    for(int t=0;t<draws+burn;++t) {
        AnnotationSampler::correctedTgs(designs,residuals,alpha,variances,pi,&cache,&proposals,&accepted);
        for(int c=0;c<components;++c) {
            assert((residual[c]-(y[c]-a[c]*alpha.col(c))).norm()<0.001);
            assert(alpha(0,c)==0); // intercept excluded from TGS
            residual[c]=y[c]-a[c]*alpha.col(c);
        }
        if(t>=burn) {
            for(int j=0;j<2;++j) observedPip[j]+=(alpha.row(j+1).array()!=0).any();
            sum+=alpha.bottomRows(2).cast<double>();
            squares+=alpha.bottomRows(2).cast<double>().array().square().matrix();
        }
    }
    assert((observedPip/draws-expectedPip).cwiseAbs().maxCoeff()<0.015);
    assert((sum/draws-expectedMean).cwiseAbs().maxCoeff()<0.02);
    assert((squares/draws-expectedSecond).cwiseAbs().maxCoeff()<0.025);
    assert(accepted>0 && accepted<proposals);
    std::cout << "Corrected TGS " << components << " components, pi=" << pi
              << ": exact PIPs/moments matched; acceptance=" << double(accepted)/proposals << '\n';
}

int main() {
    Stat::seedEngine(192837);
    gaussianTest(0); gaussianTest(2); gaussianTest(1);
    tgsTest(1,0.25); tgsTest(2,0.25); tgsTest(2,0.01); tgsTest(1,0.95);
    assert(std::isfinite(AnnotationSampler::logProbability(10000,false)));
    assert(std::isfinite(AnnotationSampler::logProbability(-10000,true)));
    std::cout << "Annotation sampler tests passed.\n";
}
