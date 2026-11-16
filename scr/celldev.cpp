#include "celldev.hpp"
#include "stat.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <exception>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <iomanip>
#include <numeric>
#include <chrono>
#include <unordered_map>
#include <unordered_set>

namespace CellDev {
using namespace Eigen;
namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::invalid_argument(message);
}
MatrixXd helmert(int n) {
    MatrixXd q = MatrixXd::Zero(n, std::max(0,n-1));
    for (int j=0; j<n-1; ++j) {
        const double scale=std::sqrt(double((j+1)*(j+2)));
        for (int i=0; i<=j; ++i) q(i,j)=1.0/scale;
        q(j+1,j)=-(j+1)/scale;
    }
    return q;
}
double openUniform() {
    double u;
    do { u=Stat::ranf(); } while (u<=0.0 || u>=1.0);
    return u;
}

double wallSeconds() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

MatrixXd mainWeights(const Metadata& m) {
    MatrixXd weights=MatrixXd::Zero(m.assayCells.size(),m.mainColumns());
    require(m.assayCells.size()==m.assayStages.size(),"Invalid CellDev assay mapping.");
    std::set<std::pair<int,int> > seen;
    for (int q=0;q<int(m.assayCells.size());++q) {
        const int c=m.assayCells[q], t=m.assayStages[q];
        require(c>=0 && c<m.cells && t>=0 && t<m.stages && seen.insert({c,t}).second,
                "Invalid or duplicate CellDev assay coordinates.");
        weights(q,m.parent[c])=1;
        weights.block(q,m.parents,1,m.cells-m.parents)=m.H.row(c);
        weights.block(q,m.cells,1,m.stages-1)=m.QT.row(t);
        if (m.structured()) {
            const int start=m.cells+m.stages-1;
            weights.block(q,start,1,m.stages-1)=m.groupCoding[c]*m.QT.row(t);
            for (int r=0;r<m.QJ.cols();++r)
                weights.block(q,start+(r+1)*(m.stages-1),1,m.stages-1)=
                    m.QJ(m.parent[c],r)*m.QT.row(t);
        }
    }
    return weights;
}

template<typename Matrix>
void continuationLogs(const MatrixXd& eta,Matrix& result) {
    require(eta.allFinite(),"CellDev predictor must be finite.");
    const int splits=eta.cols();
    result.resize(eta.rows(),splits+1);
    #pragma omp parallel for schedule(static) if(eta.rows()>=16384)
    for (int j=0;j<eta.rows();++j) {
        double prefix=0;
        for (int k=0;k<splits;++k) {
            const double x=eta(j,k), absolute=std::abs(x);
            // One erfc supplies both log Phi(x) and log Phi(-x).
            const double small=absolute<30?std::log(0.5*std::erfc(absolute/std::sqrt(2.0))):logPhi(-absolute);
            const double large=std::log1p(-std::exp(small));
            result(j,k)=prefix+(x>=0?small:large);
            prefix+=x>=0?large:small;
        }
        result(j,splits)=prefix;
    }
}

// Explicit row tiling provides parallel GEMM even when production MCMC has
// disabled Eigen's internal threading. Temporaries are bounded per worker.
template<typename Derived>
void rowMultiply(Eigen::Ref<MatrixXd> output, const MatrixBase<Derived>& input,
                 const MatrixXd& weights) {
    const int tileSize=4096, tiles=(input.rows()+tileSize-1)/tileSize;
    #pragma omp parallel for schedule(static) if(input.rows()>=16384)
    for (int tile=0;tile<tiles;++tile) {
        const int first=tile*tileSize, n=std::min(tileSize,int(input.rows())-first);
        output.middleRows(first,n).noalias()=input.middleRows(first,n).template cast<double>()*weights;
    }
}

VectorXd productColumnSquares(const MatrixXf& input,const MatrixXd& weights) {
    const int tileSize=4096, tiles=(input.rows()+tileSize-1)/tileSize;
    MatrixXd partial(weights.cols(),tiles);
    #pragma omp parallel for schedule(static) if(input.rows()>=16384)
    for (int tile=0;tile<tiles;++tile) {
        const int first=tile*tileSize,n=std::min(tileSize,int(input.rows())-first);
        MatrixXd block=input.middleRows(first,n).cast<double>()*weights;
        partial.col(tile)=block.colwise().squaredNorm().transpose();
    }
    return partial.rowwise().sum();
}

VectorXd vectorProduct(const MatrixXd& design,const VectorXd& value) {
    VectorXd result(design.rows());
    const int tileSize=4096,tiles=(design.rows()+tileSize-1)/tileSize;
    #pragma omp parallel for schedule(static) if(design.rows()>=16384)
    for (int tile=0;tile<tiles;++tile) {
        const int first=tile*tileSize,n=std::min(tileSize,int(design.rows())-first);
        result.segment(first,n).noalias()=design.middleRows(first,n)*value;
    }
    return result;
}

void applyChange(Eigen::Ref<VectorXd> predictor,const MatrixXd& design,const VectorXd& change) {
    const int tileSize=4096,tiles=(design.rows()+tileSize-1)/tileSize;
    #pragma omp parallel for schedule(static) if(design.rows()>=16384)
    for (int tile=0;tile<tiles;++tile) {
        const int first=tile*tileSize,n=std::min(tileSize,int(design.rows())-first);
        predictor.segment(first,n).noalias()+=design.middleRows(first,n)*change;
    }
}

// Ascending fixed-group reduction is independent of worker scheduling/count.
// Gather changing eligible subsets once per tile, not one heap vector per SNP.
void sufficientStatistics(const MatrixXd& design, const std::vector<int>& eligible,
                          const VectorXd& residual, bool needGram,
                          MatrixXd& gram, VectorXd& information) {
    const int p=design.cols(), n=eligible.size(), tileSize=4096;
    const int groups=std::max(1,std::min(32,(n+tileSize-1)/tileSize));
    const bool full=n==design.rows(); // eligible indices are sorted and unique
    std::vector<MatrixXd> gs(groups);
    std::vector<VectorXd> hs(groups);
    #pragma omp parallel for schedule(static) if(n>=16384)
    for (int group=0;group<groups;++group) {
        hs[group]=VectorXd::Zero(p);
        if (needGram) gs[group]=MatrixXd::Zero(p,p);
        const int begin=(static_cast<long long>(n)*group)/groups;
        const int end=(static_cast<long long>(n)*(group+1))/groups;
        MatrixXd block; VectorXd response;
        for (int first=begin;first<end;first+=tileSize) {
            const int count=std::min(tileSize,end-first);
            if (full && !needGram) {
                hs[group].noalias()+=design.middleRows(first,count).transpose()*residual.segment(first,count);
                continue;
            }
            block.resize(count,p); response.resize(count);
            if (full) { block=design.middleRows(first,count); response=residual.segment(first,count); }
            else {
                for (int c=0;c<p;++c) for (int i=0;i<count;++i) block(i,c)=design(eligible[first+i],c);
                for (int i=0;i<count;++i) response[i]=residual[eligible[first+i]];
            }
            if (needGram) gs[group].noalias()+=block.transpose()*block;
            hs[group].noalias()+=block.transpose()*response;
        }
    }
    gram=needGram?MatrixXd::Zero(p,p):MatrixXd(0,0);
    information=VectorXd::Zero(p);
    for (int group=0;group<groups;++group) {
        if (needGram) gram+=gs[group];
        information+=hs[group];
    }
}
}

Metadata::Metadata(const std::vector<int>& parentOfCell, int numStages)
    : cells(parentOfCell.size()), stages(numStages), parents(0), parent(parentOfCell) {
    require(cells>0 && stages>0,"CellDev requires cells and stages.");
    for (int g : parent) {
        require(g>=0 && g<cells,"CellDev parent IDs must be contiguous from zero.");
        parents=std::max(parents,g+1);
    }
    for (int c=0;c<cells;++c) for (int t=0;t<stages;++t) {
        assayCells.push_back(c); assayStages.push_back(t);
    }
    P=MatrixXd::Zero(cells,parents);
    std::vector<std::vector<int> > groups(parents);
    for (int c=0;c<cells;++c) { P(c,parent[c])=1; groups[parent[c]].push_back(c); }
    H=MatrixXd::Zero(cells,cells-parents);
    VectorXd n(parents);
    for (int g=0,offset=0;g<parents;++g) {
        require(!groups[g].empty(),"CellDev parent IDs must be contiguous from zero.");
        n[g]=groups[g].size();
        MatrixXd hg=helmert(groups[g].size());
        for (int i=0;i<hg.rows();++i)
            H.block(groups[g][i],offset,1,hg.cols())=hg.row(i);
        offset+=hg.cols();
    }
    // Deterministic Householder basis of null(n^T), including unequal families.
    MatrixXd house=MatrixXd::Identity(parents,parents);
    VectorXd v=n; v[0]+=n.norm();
    house.noalias()-=2.0*(v*v.transpose())/v.squaredNorm();
    QG=house.rightCols(parents-1);
    QT=helmert(stages);
    QU.resize(cells,cells-1);
    QU.leftCols(parents-1)=P*QG;
    QU.rightCols(cells-parents)=H;
}

void Metadata::setNeuronGroups(const std::vector<int>& groups) {
    require(int(groups.size())==cells,"Neuron groups must label every subtype.");
    require(std::all_of(groups.begin(),groups.end(),[](int x){return x==0 || x==1;}),
            "Neuron groups must be 0 or 1.");
    const int count=std::accumulate(groups.begin(),groups.end(),0);
    require(count>0 && count<cells,"Stage hierarchy requires both neuron and non-neuron subtypes.");
    neuron=groups;
    groupCoding.resize(cells);
    VectorXd n=VectorXd::Zero(parents), r=VectorXd::Zero(parents);
    for (int c=0;c<cells;++c) {
        groupCoding[c]=groups[c]-double(count)/cells;
        ++n[parent[c]]; r[parent[c]]+=groupCoding[c];
    }
    parentGroupCoding=r.cwiseQuotient(n);
    // Orthogonal to shared stage and group-stage designs on the subtype grid.
    MatrixXd constraints(2,parents); constraints.row(0)=n.transpose(); constraints.row(1)=r.transpose();
    JacobiSVD<MatrixXd> svd(constraints,ComputeFullV);
    svd.setThreshold(1e-12);
    QJ=svd.matrixV().rightCols(parents-svd.rank());
}

CoordinateLabels coordinateLabels(const Metadata& m,int baselines,int splits,int rank) {
    CoordinateLabels result;
    auto name=[](const std::vector<std::string>& names,int i,const std::string& prefix) {
        return i<int(names.size())?names[i]:prefix+std::to_string(i+1);
    };
    result.linear.push_back("Intercept");
    for (int i=0;i<baselines;++i) result.linear.push_back("Baseline:"+name(m.baselineNames,i,"baseline"));
    for (int g=0;g<m.parents;++g) {
        const std::string parent=name(m.parentNames,g,"parent");
        result.linear.push_back("Parent:"+parent);
        const int count=std::count(m.parent.begin(),m.parent.end(),g);
        for (int h=0;h<count-1;++h) result.subtypeContrasts.push_back("H:"+parent+":"+std::to_string(h+1));
    }
    for (const auto& h:result.subtypeContrasts) result.linear.push_back("SubtypeContrast:"+h);
    for (int t=0;t<m.stages-1;++t) {
        result.stageContrasts.push_back("QT:"+std::to_string(t+1));
        result.linear.push_back("StageContrast:"+result.stageContrasts.back());
    }
    if (m.structured()) {
        for (int t=0;t<m.stages-1;++t)
            result.linear.push_back("NeuronStageContrast:QT:"+std::to_string(t+1));
        for (int r=0;r<m.QJ.cols();++r) for (int t=0;t<m.stages-1;++t)
            result.linear.push_back("ParentStageContrast:QJ:"+std::to_string(r+1)+":QT:"+std::to_string(t+1));
        rank=0;
    }
    for (int g=0;g<m.parents-1;++g) result.parentContrasts.push_back("QG:"+std::to_string(g+1));
    if (m.cells==1 || m.stages==1) rank=0;
    for (int r=0;r<rank;++r) {
        const std::string factor="factor"+std::to_string(r+1)+":";
        for (const auto& q:result.parentContrasts) result.cellFactor.push_back(factor+"ParentContrast:"+q);
        for (const auto& h:result.subtypeContrasts) result.cellFactor.push_back(factor+"SubtypeContrast:"+h);
        for (const auto& t:result.stageContrasts) result.stageFactor.push_back(factor+"StageContrast:"+t);
    }
    for (int k=0;k<splits;++k) for (int r=0;r<rank;++r)
        result.strengths.push_back("split"+std::to_string(k+1)+":factor"+std::to_string(r+1));
    return result;
}

double logPhi(double x) {
    require(std::isfinite(x),"CellDev predictor must be finite.");
    if (x>=0) return std::log1p(-0.5*std::erfc(x/std::sqrt(2.0)));
    if (x>-30) return std::log(0.5*std::erfc(-x/std::sqrt(2.0)));
    // Mills-ratio asymptotic expansion. No probability clipping.
    const double inv=1.0/(x*x);
    double term=1, sum=1, previous=1;
    for (int k=1;k<100;++k) {
        term*=-(2*k-1)*inv;
        if (std::abs(term)>previous) break;
        sum+=term;
        if (std::abs(term)<std::abs(sum)*1e-16) break;
        previous=std::abs(term);
    }
    return -0.5*x*x-std::log(-x)-0.5*std::log(2*std::acos(-1.0))+std::log(sum);
}

double truncatedProbit(double mean, bool success) {
    require(std::isfinite(mean),"CellDev truncated normal mean must be finite.");
    if (!success) return -truncatedProbit(-mean,true);
    if (mean < -0.45) {
        const double bound=-mean;
        const double rate=0.5*bound+0.5*std::hypot(bound,2.0);
        for (;;) {
            const double excess=-std::log(openUniform())/rate;
            const double difference=excess-1.0/rate;
            if (std::log(openUniform())<=-0.5*difference*difference)
                return excess; // avoids cancellation of mean + bound
        }
    }
    double value;
    do { value=mean+Stat::snorm(); } while (value<=0);
    return value;
}

VectorXd gaussianDraw(const MatrixXd& precision, const VectorXd& information) {
    require(precision.rows()==precision.cols() && precision.rows()==information.size()
            && precision.allFinite() && information.allFinite(),"Invalid CellDev Gaussian system.");
    if (!information.size()) return VectorXd(0);
    LLT<MatrixXd> chol(precision);
    if (chol.info()!=Success) throw std::runtime_error("CellDev precision is not positive definite.");
    VectorXd noise(information.size());
    for (int i=0;i<noise.size();++i) noise[i]=Stat::snorm();
    return chol.solve(information)+chol.matrixU().solve(noise);
}

MatrixXd Sampler::buildMainDesign(const Metadata& m, const MatrixXf& annotation,
                                  const MatrixXf& baseline) {
    require(annotation.cols()==int(m.assayCells.size()) && annotation.rows()>0
            && baseline.rows()==annotation.rows() && annotation.allFinite() && baseline.allFinite(),
            "Invalid CellDev annotation/baseline dimensions or values.");
    const int nb=baseline.cols();
    MatrixXd x=MatrixXd::Zero(annotation.rows(),1+nb+m.mainColumns());
    x.col(0).setOnes();
    x.middleCols(1,nb)=baseline.cast<double>();
    const MatrixXd weights=mainWeights(m);
    rowMultiply(x.rightCols(weights.cols()),annotation,weights);
    return x;
}

Sampler::Sampler(const Metadata& m, const MatrixXf& annotation, const MatrixXf& baseline,
                 int k, int r, const Priors& prior, const Sampler* sharedDesigns)
    : annotationStorage(sharedDesigns?sharedDesigns->annotationStorage:
          std::make_shared<const MatrixXd>(annotation.cast<double>())),
      designStorage(sharedDesigns?sharedDesigns->designStorage:
          std::make_shared<const MatrixXd>(buildMainDesign(m,annotation,baseline))),
      metadata(m), splits(k), rank((m.structured() || m.cells==1 || m.stages==1)?0:r),
      annotations(*annotationStorage), mainDesign(*designStorage),
      priors(prior) {
    require(k>0 && r>=0,"Invalid CellDev split count or rank.");
    require(!m.structured() || !priors.rd,"Explicit stage hierarchy requires --celldev-prior Gaussian.");
    require(priors.linearMean.rows()==mainDesign.cols() && priors.linearMean.cols()==splits
            && priors.linearVariance.rows()==mainDesign.cols() && priors.linearVariance.cols()==splits
            && priors.linearMean.allFinite() && priors.linearVariance.allFinite()
            && (priors.linearVariance.array()>0).all(),"Invalid CellDev linear Gaussian priors.");
    require(priors.subtypeFactorSD>0 && std::isfinite(priors.subtypeFactorSD),
            "Invalid CellDev subtype factor scale.");
    if (rank) require(priors.interactionVariance.size()==splits
            && priors.interactionVariance.allFinite() && (priors.interactionVariance.array()>0).all(),
            "Invalid CellDev interaction prior variances.");
    if (priors.rd) require(priors.linearMean.bottomRows(priors.linearMean.rows()-1).isZero(),
                           "CellDev RD requires zero non-intercept prior means.");
    linear=priors.linearMean;
    rdWeights=MatrixXd::Constant(splits,5,0.2);
    rdInclusion=MatrixXd::Zero(linear.rows(),splits);
    rdAllocation=MatrixXi::Zero(linear.rows(),splits);
    if (priors.rd) linear.bottomRows(linear.rows()-1).setZero();
    a.resize(m.cells-1,rank); f.resize(m.stages-1,rank); lambda.resize(rank,splits);
    for (int r0=0;r0<rank;++r0) {
        for (int i=0;i<a.rows();++i) a(i,r0)=0.1*Stat::snorm();
        for (int i=0;i<f.rows();++i) f(i,r0)=0.1*Stat::snorm();
        for (int i=0;i<splits;++i) lambda(r0,i)=0.1*Stat::snorm();
    }
    eligible.resize(splits);
    z=MatrixXd::Constant(annotation.rows(),splits,std::numeric_limits<double>::quiet_NaN());
    rebuildEta();
}

MatrixXd Sampler::coefficients() const {
    const Metadata& m=metadata;
    const int offset=1+baselines();
    MatrixXd cell=m.P*linear.middleRows(offset,m.parents)
                  +m.H*linear.middleRows(offset+m.parents,m.cells-m.parents);
    MatrixXd stage=stageEffects();
    const MatrixXd difference=stageDifference(), interaction=parentInteractions();
    MatrixXd u=m.QU*a, v=m.QT*f;
    MatrixXd alpha(m.cells*m.stages,splits);
    for (int c=0;c<m.cells;++c) for (int t=0;t<m.stages;++t) {
        alpha.row(c*m.stages+t)=cell.row(c)+stage.row(t);
        if (m.structured()) alpha.row(c*m.stages+t)+=m.groupCoding[c]*difference.row(t)
            +interaction.row(m.parent[c]*m.stages+t);
        for (int r=0;r<rank;++r) alpha.row(c*m.stages+t)+=u(c,r)*v(t,r)*lambda.row(r);
    }
    return alpha;
}

MatrixXd Sampler::parentEffects() const {
    return linear.middleRows(1+baselines(),metadata.parents);
}
MatrixXd Sampler::subtypeEffects() const {
    return metadata.H*linear.middleRows(1+baselines()+metadata.parents,metadata.cells-metadata.parents);
}
MatrixXd Sampler::stageEffects() const {
    return metadata.QT*linear.middleRows(1+baselines()+metadata.cells,metadata.stages-1);
}
MatrixXd Sampler::stageDifference() const {
    if (!metadata.structured()) return MatrixXd::Zero(metadata.stages,splits);
    return metadata.QT*linear.middleRows(1+baselines()+metadata.cells+metadata.stages-1,metadata.stages-1);
}
MatrixXd Sampler::parentInteractions() const {
    const auto& m=metadata;
    MatrixXd result=MatrixXd::Zero(m.parents*m.stages,splits);
    if (!m.structured()) return result;
    const int start=1+baselines()+m.cells+2*(m.stages-1);
    for (int r=0;r<m.QJ.cols();++r) {
        const MatrixXd profile=m.QT*linear.middleRows(start+r*(m.stages-1),m.stages-1);
        for (int g=0;g<m.parents;++g) for (int t=0;t<m.stages;++t)
            result.row(g*m.stages+t)+=m.QJ(g,r)*profile.row(t);
    }
    return result;
}

std::vector<std::string> Sampler::componentNames() const {
    const auto& m=metadata;
    auto label=[&](const std::vector<std::string>& ids,int i) {
        return i<int(ids.size())?ids[i]:std::to_string(i+1);
    };
    const std::string prefix=m.modality.empty()?"":m.modality+":";
    std::vector<std::string> names;
    for (int g=0;g<m.parents;++g) names.push_back("Parent:"+prefix+label(m.parentNames,g));
    for (int c=0;c<m.cells;++c) names.push_back("Subtype:"+prefix+label(m.cellNames,c));
    for (int t=0;t<m.stages;++t) names.push_back("Stage:"+prefix+label(m.stageNames,t));
    if (!m.structured()) return names;
    for (int c=0;c<m.cells;++c) names.push_back("SubtypeTotal:"+prefix+label(m.cellNames,c));
    for (const std::string kind:{"StageNeuron:","StageNonNeuron:","StageDifference:"})
        for (int t=0;t<m.stages;++t) names.push_back(kind+prefix+label(m.stageNames,t));
    for (const std::string kind:{"ParentStageInteraction:","ParentStage:"})
        for (int g=0;g<m.parents;++g) for (int t=0;t<m.stages;++t)
            names.push_back(kind+prefix+label(m.parentNames,g)+":"+label(m.stageNames,t));
    for (int c=0;c<m.cells;++c) for (int t=0;t<m.stages;++t)
        names.push_back("SubtypeStage:"+prefix+label(m.cellNames,c)+":"+label(m.stageNames,t));
    return names;
}

MatrixXd Sampler::componentContrasts() const {
    const auto& m=metadata;
    const MatrixXd parent=parentEffects(),subtype=subtypeEffects(),stage=stageEffects();
    MatrixXd result(componentNames().size(),splits);
    int row=0;
    for (const MatrixXd* block:{&parent,&subtype,&stage}) {
        result.middleRows(row,block->rows())=*block; row+=block->rows();
    }
    if (!m.structured()) return result;
    for (int c=0;c<m.cells;++c) result.row(row++)=parent.row(m.parent[c])+subtype.row(c);
    const MatrixXd difference=stageDifference(), interaction=parentInteractions();
    const double p=double(std::accumulate(m.neuron.begin(),m.neuron.end(),0))/m.cells;
    result.middleRows(row,m.stages)=stage+(1-p)*difference; row+=m.stages;
    result.middleRows(row,m.stages)=stage-p*difference; row+=m.stages;
    result.middleRows(row,m.stages)=difference; row+=m.stages;
    result.middleRows(row,interaction.rows())=interaction; row+=interaction.rows();
    for (int g=0;g<m.parents;++g) for (int t=0;t<m.stages;++t)
        result.row(row++)=parent.row(g)+stage.row(t)+m.parentGroupCoding[g]*difference.row(t)
            +interaction.row(g*m.stages+t);
    result.bottomRows(m.cells*m.stages)=coefficients();
    return result;
}

MatrixXd Sampler::observedCoefficients() const {
    const MatrixXd full=coefficients();
    MatrixXd observed(metadata.assayCells.size(),splits);
    for (int q=0;q<observed.rows();++q)
        observed.row(q)=full.row(metadata.assayCells[q]*metadata.stages+metadata.assayStages[q]);
    return observed;
}

MatrixXd Sampler::predictor() const {
    MatrixXd result(mainDesign.rows(),splits);
    rowMultiply(result,mainDesign,linear);
    result.noalias()+=factorScores()*lambda;
    return result;
}
void Sampler::rebuildEta() { eta=predictor(); }

MatrixXd Sampler::logProbabilities() const {
    MatrixXd result;
    continuationLogs(eta,result);
    return result;
}
void Sampler::logProbabilities(MatrixXf& result) const {
    continuationLogs(eta,result);
    require(result.allFinite(),"CellDev log probabilities exceed the SNP kernel's floating-point range.");
}

MatrixXd Sampler::cellDesign(int r) const {
    require(r>=0 && r<rank,"Invalid CellDev factor index.");
    VectorXd v=metadata.QT*f.col(r);
    MatrixXd weights(annotations.cols(),a.rows());
    for (int q=0;q<weights.rows();++q)
        weights.row(q)=v[metadata.assayStages[q]]*metadata.QU.row(metadata.assayCells[q]);
    MatrixXd design(annotations.rows(),a.rows());
    rowMultiply(design,annotations,weights);
    return design;
}
MatrixXd Sampler::stageDesign(int r) const {
    require(r>=0 && r<rank,"Invalid CellDev factor index.");
    VectorXd u=metadata.QU*a.col(r);
    MatrixXd weights(annotations.cols(),f.rows());
    for (int q=0;q<weights.rows();++q)
        weights.row(q)=u[metadata.assayCells[q]]*metadata.QT.row(metadata.assayStages[q]);
    MatrixXd design(annotations.rows(),f.rows());
    rowMultiply(design,annotations,weights);
    return design;
}
const MatrixXd& Sampler::factorScores() const {
    if (scoresCache.rows()==annotations.rows() && scoresCache.cols()==rank
        && scoresA.rows()==a.rows() && scoresA.cols()==a.cols()
        && scoresF.rows()==f.rows() && scoresF.cols()==f.cols()
        && (scoresA.array()==a.array()).all() && (scoresF.array()==f.array()).all()) return scoresCache;
    MatrixXd weights(annotations.cols(),rank);
    const MatrixXd u=metadata.QU*a, v=metadata.QT*f;
    for (int q=0;q<weights.rows();++q)
        weights.row(q)=u.row(metadata.assayCells[q]).cwiseProduct(v.row(metadata.assayStages[q]));
    scoresCache.resize(annotations.rows(),rank);
    if (rank) rowMultiply(scoresCache,annotations,weights);
    scoresA=a; scoresF=f;
    return scoresCache;
}

void Sampler::invalidateZ() { augmented=false; }
void Sampler::requireAugmentation() const {
    if (!augmented) throw std::logic_error("Refresh CellDev z after SNP allocation before annotation updates.");
}
void Sampler::sampleZ(const VectorXi& membership) {
    invalidateZ();
    require(membership.size()==eta.rows(),"Invalid CellDev membership dimensions.");
    require(eta.allFinite(),"CellDev predictor must be finite.");
    require((membership.array()>=0).all() && (membership.array()<=splits).all(),
            "Invalid CellDev mixture membership.");
    z.setConstant(std::numeric_limits<double>::quiet_NaN());
    ++augmentationSweep;
    for (int k=0;k<splits;++k) {
        eligible[k].clear();
        if (k==0) { eligible[k].resize(eta.rows()); std::iota(eligible[k].begin(),eligible[k].end(),0); }
        else for (int j=0;j<eta.rows();++j) if (membership[j]>=k) eligible[k].push_back(j);
        const int tileSize=4096, n=eligible[k].size(), tiles=(n+tileSize-1)/tileSize;
        const bool parallelStreams=n>=16384;
        #pragma omp parallel for schedule(static) if(parallelStreams)
        for (int tile=0;tile<tiles;++tile) {
            // Logical tile streams are independent of worker count. Include an
            // augmentation counter so repeated sweeps at one MCMC iteration
            // still draw fresh z. Preserve the caller's scalar RNG stream.
            if (parallelStreams) {
                Stat::ParallelTaskRngScope scope(tile*splits+k,0x43445A00u^(augmentationSweep*0x9E3779B9u));
                for (int i=tile*tileSize;i<std::min(n,(tile+1)*tileSize);++i) {
                    const int j=eligible[k][i]; z(j,k)=truncatedProbit(eta(j,k),membership[j]>k);
                }
            } else for (int i=tile*tileSize;i<std::min(n,(tile+1)*tileSize);++i) {
                const int j=eligible[k][i]; z(j,k)=truncatedProbit(eta(j,k),membership[j]>k);
            }
        }
    }
    augmented=true;
}

void Sampler::updateSplit(int k, const MatrixXd& design, VectorXd& value,
                          const VectorXd& mean, const VectorXd& variance, bool linearBlock) {
    requireAugmentation();
    require(k>=0 && k<splits,"Invalid CellDev split index.");
    VectorXd inverse=variance.cwiseInverse();
    MatrixXd precision=inverse.asDiagonal();
    VectorXd information=inverse.cwiseProduct(mean);
    // The full linear block's offset is exactly the interaction contribution.
    // Avoid recomputing X*old_linear for every split merely to cancel it.
    VectorXd residual;
    if (linearBlock) residual=z.col(k)-vectorProduct(factorScores(),lambda.col(k));
    else residual=z.col(k)-eta.col(k)+vectorProduct(design,value);
    MatrixXd gram; VectorXd score;
    const bool cacheGram=linearBlock && k==0;
    const bool cached=cacheGram && firstSplitGram.rows()==design.cols();
    sufficientStatistics(design,eligible[k],residual,!cached,gram,score);
    if (cacheGram && !cached) firstSplitGram=gram;
    precision+=cached?firstSplitGram:gram;
    information+=score;
    VectorXd next=gaussianDraw(precision,information);
    applyChange(eta.col(k),design,next-value);
    value=next;
}
void Sampler::sampleLinear(int k) {
    if (priors.rd) { sampleLinearRD(k); return; }
    require(k>=0 && k<splits,"Invalid CellDev split index.");
    if (fitIntercept) {
        VectorXd value=linear.col(k);
        updateSplit(k,mainDesign,value,priors.linearMean.col(k),priors.linearVariance.col(k),true);
        linear.col(k)=value;
    } else {
        const int n=mainDesign.cols()-1;
        VectorXd value=linear.col(k).tail(n);
        updateSplit(k,mainDesign.rightCols(n),value,priors.linearMean.col(k).tail(n),
                    priors.linearVariance.col(k).tail(n),true);
        linear.col(k).tail(n)=value;
    }
}

// Relative slab multipliers normalised to preserve the calibrated variance
// under the initial uniform mixture (mean of 0,.0001,.001,.01,.1).
static const double rdMultipliers[5]={0, .0001/.02222, .001/.02222, .01/.02222, .1/.02222};
VectorXd Sampler::rdProbabilities(double diagonal, double score, double scale,
                                  const VectorXd& weights) {
    require(diagonal>=0 && std::isfinite(diagonal) && std::isfinite(score)
            && scale>0 && std::isfinite(scale) && weights.size()==5
            && weights.allFinite() && (weights.array()>0).all(),"Invalid CellDev RD conditional.");
    VectorXd logp(5);
    for (int g=0;g<5;++g) {
        const double v=scale*rdMultipliers[g];
        logp[g]=std::log(weights[g])-.5*std::log1p(diagonal*v)
                +.5*score*(v/(1+diagonal*v))*score;
    }
    VectorXd p=(logp.array()-logp.maxCoeff()).exp();
    return p/p.sum();
}
void Sampler::sampleLinearRD(int k) {
    requireAugmentation();
    require(k>=0 && k<splits,"Invalid CellDev split index.");
    VectorXf counts=VectorXf::Ones(5);
    for (int q=fitIntercept?0:1;q<linear.rows();++q) {
        double diagonal=0, score=0;
        const double old=linear(q,k);
        for (int j:eligible[k]) {
            const double x=mainDesign(j,q);
            diagonal+=x*x;
            score+=x*(z(j,k)-eta(j,k)+x*old);
        }
        double next;
        if (q==0) {
            const double v=1/(diagonal+1/priors.linearVariance(q,k));
            next=v*(score+priors.linearMean(q,k)/priors.linearVariance(q,k))+std::sqrt(v)*Stat::snorm();
        } else {
            const VectorXd p=rdProbabilities(diagonal,score,priors.linearVariance(q,k),rdWeights.row(k).transpose());
            double u=Stat::ranf(); int g=0;
            while (g<4 && u>p[g]) { u-=p[g]; ++g; }
            rdAllocation(q,k)=g; rdInclusion(q,k)=1-p[0]; ++counts[g];
            next=0;
            if (g) {
                const double v=1/(diagonal+1/(priors.linearVariance(q,k)*rdMultipliers[g]));
                next=v*score+std::sqrt(v)*Stat::snorm();
            }
        }
        eta.col(k).noalias()+=mainDesign.col(q)*(next-old);
        linear(q,k)=next;
    }
    rdWeights.row(k)=Stat::Dirichlet().sample(5,counts).cast<double>().transpose();
}

void Sampler::updateSharedFactor(int r, bool cell) {
    requireAugmentation();
    MatrixXd design=cell?cellDesign(r):stageDesign(r);
    VectorXd value=cell?VectorXd(a.col(r)):VectorXd(f.col(r));
    VectorXd inverse=VectorXd::Ones(value.size());
    if (cell) inverse.tail(metadata.cells-metadata.parents).setConstant(
        1.0/(priors.subtypeFactorSD*priors.subtypeFactorSD));
    MatrixXd precision=inverse.asDiagonal();
    VectorXd information=VectorXd::Zero(value.size());
    const VectorXd contribution=vectorProduct(design,value);
    for (int k=0;k<splits;++k) {
        VectorXd residual=z.col(k)-eta.col(k)+lambda(r,k)*contribution;
        MatrixXd gram; VectorXd score;
        sufficientStatistics(design,eligible[k],residual,true,gram,score);
        precision.noalias()+=lambda(r,k)*lambda(r,k)*gram;
        information.noalias()+=lambda(r,k)*score;
    }
    VectorXd next=gaussianDraw(precision,information);
    VectorXd delta=vectorProduct(design,next-value);
    eta.noalias()+=delta*lambda.row(r);
    if (cell) a.col(r)=next; else f.col(r)=next;
}
void Sampler::sampleA(int r) { updateSharedFactor(r,true); }
void Sampler::sampleF(int r) { updateSharedFactor(r,false); }
void Sampler::sampleLambda(int k) {
    require(k>=0 && k<splits,"Invalid CellDev split index.");
    requireAugmentation();
    if (!rank) return;
    VectorXd value=lambda.col(k);
    updateSplit(k,factorScores(),value,VectorXd::Zero(rank),
                VectorXd::Constant(rank,priors.interactionVariance[k]/rank));
    lambda.col(k)=value;
}
void Sampler::sweep(const VectorXi& membership) {
    const double start=wallSeconds();
    sampleZ(membership);
    const double latent=wallSeconds();
    for (int k=0;k<splits;++k) sampleLinear(k);
    const double main=wallSeconds();
    for (int r=0;r<rank;++r) { sampleA(r); sampleF(r); }
    const double factors=wallSeconds();
    for (int k=0;k<splits;++k) sampleLambda(k);
    const double end=wallSeconds();
    timing.z=latent-start; timing.linear=main-latent;
    timing.factors=factors-main; timing.lambda=end-factors;
}

void Sampler::sweepConditional(const Sampler& augmentation, const MatrixXd& offset, bool intercept) {
    require(augmentation.augmented,"Joint CellDev requires shared augmentation.");
    const double start=wallSeconds();
    fitIntercept=intercept;
    if (!intercept) require(linear.row(0).isZero(),"Extra modality intercept must be zero.");
    eligible=augmentation.eligible;
    if (&augmentation==this) z-=offset;
    else z=augmentation.z-offset;
    augmented=true;
    for (int k=0;k<splits;++k) sampleLinear(k);
    const double main=wallSeconds();
    for (int r=0;r<rank;++r) { sampleA(r); sampleF(r); }
    const double factors=wallSeconds();
    for (int k=0;k<splits;++k) sampleLambda(k);
    timing.z=0; timing.linear=main-start;
    timing.factors=factors-main; timing.lambda=wallSeconds()-factors;
}

namespace {
std::vector<std::string> enrichmentNames(const Sampler& s) {
    std::vector<std::string> names={"Intercept"};
    const int baselines=s.baselines();
    if (!s.metadata.enrichmentColumns.empty()) {
        require(int(s.metadata.enrichmentColumns.size())==baselines+s.annotations.cols(),
                "Invalid CellDev enrichment annotation mapping.");
        for (const auto& info:s.metadata.enrichmentColumns) names.push_back(info.name);
    } else {
        for (int i=0;i<baselines;++i) names.push_back(i<int(s.metadata.baselineNames.size())?
            s.metadata.baselineNames[i]:"Baseline"+std::to_string(i+1));
        const auto& m=s.metadata;
        for (int q=0;q<int(m.assayCells.size());++q) {
            const int c=m.assayCells[q], t=m.assayStages[q];
            names.push_back((c<int(m.cellNames.size())?m.cellNames[c]:"cell"+std::to_string(c+1))+
                "_"+(t<int(m.stageNames.size())?m.stageNames[t]:"stage"+std::to_string(t+1)));
        }
    }
    const auto& m=s.metadata;
    if (int(m.cellNames.size())==m.cells && int(m.stageNames.size())==m.stages)
        for (int q=0;q<int(m.assayCells.size());++q)
            names[1+baselines+q]=(m.modality.empty()?"":m.modality+":")+
                m.cellNames[m.assayCells[q]]+"_"+m.stageNames[m.assayStages[q]];
    require(std::set<std::string>(names.begin(),names.end()).size()==names.size(),
            "CellDev biological output names must be unique across baseline and cell-stage annotations.");
    return names;
}
std::vector<std::string> conditionalEnrichmentNames(const Sampler& s) {
    auto names=enrichmentNames(s);
    const auto components=s.componentNames();
    names.insert(names.end(),components.begin(),components.end());
    require(std::set<std::string>(names.begin(),names.end()).size()==names.size(),
            "CellDev conditional enrichment names must be unique.");
    return names;
}
double logExpectedMixtureVariance(const VectorXd& predictor,const VectorXf& gamma) {
    double prefix=0, result=-std::numeric_limits<double>::infinity();
    for (int k=0;k<predictor.size();++k) {
        // Stopping at component k (k=0 is the point mass at zero).
        if (k) {
            const double term=prefix+logPhi(-predictor[k])+std::log(double(gamma[k]));
            if (!std::isfinite(result)) result=term;
            else { const double high=std::max(result,term); result=high+std::log1p(std::exp(std::min(result,term)-high)); }
        }
        prefix+=logPhi(predictor[k]);
    }
    const double last=prefix+std::log(double(gamma[predictor.size()]));
    if (!std::isfinite(result)) return last;
    const double high=std::max(result,last);
    return high+std::log1p(std::exp(std::min(result,last)-high));
}
}

Enrichment::Enrichment(const Sampler& s)
    : names(enrichmentNames(s)), conditionalNames(conditionalEnrichmentNames(s)),
      marginal("Marginal_Heritability_Enrichment",names),
      conditional("Conditional_Heritability_Enrichment",conditionalNames) {
    const int baselines=s.baselines();
    columns=s.metadata.enrichmentColumns;
    const bool infer=columns.empty();
    if (infer) columns.resize(names.size()-1);
    means.resize(columns.size()); ssq.resize(columns.size()); sd.resize(columns.size());
    for (int i=0;i<int(columns.size());++i) {
        const VectorXd x=i<baselines?VectorXd(s.mainDesign.col(i+1)):VectorXd(s.annotations.col(i-baselines));
        means[i]=x.mean(); ssq[i]=x.squaredNorm();
        sd[i]=std::sqrt(std::max(0.0,ssq[i]/x.size()-means[i]*means[i]));
        if (infer) {
            columns[i].name=names[i+1];
            columns[i].binary=(x.array()==0 || x.array()==1).all();
            columns[i].fraction=means[i];
        }
        require(std::isfinite(columns[i].offset) && std::isfinite(columns[i].scale)
                && columns[i].scale>=0 && (!columns[i].binary ||
                (columns[i].fraction>=0 && columns[i].fraction<=1)),
                "Invalid CellDev enrichment column transformations.");
    }
    marginal.values.setOnes(); conditional.values.setOnes();
}

void Enrichment::compute(const Sampler& s,const VectorXf& effects,const VectorXf& gamma, const VectorXd* referenceOverride) {
    require(effects.size()==s.annotations.rows() && effects.allFinite() && gamma.allFinite()
            && gamma.size()==s.splits+1 && gamma[0]==0
            && (gamma.tail(s.splits).array()>0).all(),"Invalid CellDev enrichment inputs.");
    const int baselines=s.baselines();
    std::vector<int> active; std::vector<double> variance;
    double total=0;
    for (int j=0;j<effects.size();++j) if (effects[j]!=0) {
        const double v=double(effects[j])*effects[j];
        active.push_back(j); variance.push_back(v); total+=v;
    }
    marginal.values.setZero(); marginal.values[0]=1;
    // Sum-of-squared-SNP-effect convention, identical to RC's active kernel.
    #pragma omp parallel for schedule(static) if(columns.size()>=16)
    for (int i=0;i<int(columns.size());++i) {
        if (!(total>0)) continue; // RC's all-zero SNP draw convention
        double weighted=0;
        if (i<baselines) {
            const double* x=s.mainDesign.col(i+1).data();
            for (int n=0;n<int(active.size());++n) weighted+=x[active[n]]*variance[n];
        } else {
            const double* x=s.annotations.col(i-baselines).data();
            for (int n=0;n<int(active.size());++n) weighted+=x[active[n]]*variance[n];
        }
        const auto& info=columns[i];
        if (info.binary) {
            marginal.values[i+1]=info.fraction>0?
                (info.offset*total+info.scale*weighted)/(total*info.fraction):1.0;
        } else {
            const double denominator=ssq[i]-effects.size()*means[i]*means[i];
            marginal.values[i+1]=sd[i]>0 && denominator>0?
                1+(weighted-means[i]*total)/denominator*sd[i]/(total/effects.size()):1.0;
        }
    }
    MatrixXd coefficients(columns.size(),s.splits);
    coefficients.topRows(baselines)=s.linear.middleRows(1,baselines);
    coefficients.bottomRows(s.annotations.cols())=s.observedCoefficients();
    VectorXd reference=s.linear.row(0).transpose();
    for (int i=0;i<int(columns.size());++i) if (columns[i].binary && columns[i].scale>0)
        reference.noalias()-=columns[i].offset/columns[i].scale*coefficients.row(i).transpose();
    if (referenceOverride) reference=*referenceOverride;
    const double baseline=logExpectedMixtureVariance(reference,gamma);
    conditional.values.setOnes();
    for (int i=0;i<int(columns.size());++i) {
        const auto& info=columns[i];
        const double change=info.binary?(info.scale>0 && info.fraction>0 && info.fraction<1?1/info.scale:0):sd[i];
        if (!change) continue;
        VectorXd shifted=reference+change*coefficients.row(i).transpose();
        conditional.values[i+1]=std::exp(logExpectedMixtureVariance(shifted,gamma)-baseline);
    }
    // Main effects are already on the fitted predictor scale. Reconstruct
    // biological subtype/stage effects rather than reporting basis coordinates.
    const MatrixXd components=s.componentContrasts();
    for (int row=0;row<components.rows();++row)
        conditional.values[names.size()+row]=std::exp(logExpectedMixtureVariance(
            reference+components.row(row).transpose(),gamma)-baseline);
    require(marginal.values.allFinite() && conditional.values.allFinite(),
            "CellDev enrichment exceeds the supported floating-point range.");
    // Undefined annotation contrasts remain NaN in traces and become NA in
    // summaries. Main-effect rows are separate model contrasts, not SNP sets.
    for (int i=0;i<int(columns.size());++i)
        if (!(sd[i]>0) || (columns[i].binary &&
            !(columns[i].fraction>0 && columns[i].fraction<1))) {
            marginal.values[i+1]=conditional.values[i+1]=std::numeric_limits<float>::quiet_NaN();
        }

}

namespace {
std::vector<std::string> fields(const std::string& line) {
    std::istringstream input(line);
    std::vector<std::string> result;
    std::string word;
    while (input>>word) result.push_back(word);
    return result;
}
void inputError(const std::string& message) {
    throw std::runtime_error("Error: CellDev "+message);
}
double number(const std::string& text, const std::string& id,const std::string& column) {
    size_t used=0;
    double value;
    try { value=std::stod(text,&used); }
    catch (const std::exception&) { inputError("invalid numeric value at "+id+"/"+column+": "+text); return 0; }
    if (used!=text.size() || !std::isfinite(value)) inputError("nonfinite/invalid value at "+id+"/"+column+": "+text);
    return value;
}
}

std::vector<std::string> Input::selectModalities(const std::string& mapFile, const std::string& selection) {
    std::ifstream input(mapFile);
    if (!input) inputError("cannot read metadata file "+mapFile);
    std::set<std::string> available;
    std::string line;
    std::getline(input,line);
    while (std::getline(input,line)) {
        auto row=fields(line);
        if ((row.size()==6 || row.size()==7) && row[1]=="cell_stage") available.insert(row[2]);
    }
    if (available.empty()) inputError("metadata contains no cell_stage annotations.");
    if (selection.empty() || selection=="all") return {available.begin(),available.end()};
    std::vector<std::string> selected;
    std::set<std::string> seen;
    std::istringstream list(selection);
    while (std::getline(list,line,',')) {
        if (!available.count(line)) inputError("unknown modality "+line);
        if (!seen.insert(line).second) inputError("duplicate selected modality "+line);
        selected.push_back(line);
    }
    if (selected.empty()) inputError("empty modality selection.");
    return selected;
}

Input Input::read(const std::string& annotationFile, const std::string& mapFile,
                  const std::vector<std::string>& snpIds, const std::string& requestedModality, bool includeBaseline) {
    std::cout<<"CellDev: reading annotation metadata from ["<<mapFile<<"]."<<std::endl;
    struct Row { std::string role,modality,parent,cell,stage; };
    std::ifstream map(mapFile);
    if (!map) inputError("cannot read metadata file "+mapFile);
    std::string line;
    std::getline(map,line);
    const auto header=fields(line);
    const std::vector<std::string> required={"annotation_id","role","modality","parent_id","subtype_id","stage_id"};
    auto groupedHeader=required; groupedHeader.push_back("neuron_group");
    const bool grouped=header==groupedHeader;
    if (header!=required && !grouped) inputError("metadata header must be: annotation_id role modality parent_id subtype_id stage_id [neuron_group].");
    std::map<std::string,Row> rows;
    std::set<std::string> modalities;
    std::map<std::string,std::string> cellParents;
    std::map<std::string,int> cellGroups;
    std::set<std::vector<std::string> > combinations;
    while (std::getline(map,line)) {
        auto x=fields(line); if (x.empty()) continue;
        if (x.size()!=header.size()) inputError("metadata rows must match the header field count.");
        if (!rows.emplace(x[0],Row{x[1],x[2],x[3],x[4],x[5]}).second)
            inputError("duplicate annotation ID "+x[0]);
        if (x[1]=="cell_stage") {
            for (int i=2;i<6;++i) if (x[i]=="NA") inputError("missing cell-stage metadata for "+x[0]);
            modalities.insert(x[2]);
            if (!combinations.insert({x[2],x[4],x[5]}).second) inputError("duplicate modality/subtype/stage combination for "+x[0]);
            auto inserted=cellParents.emplace(x[4],x[3]);
            if (!inserted.second && inserted.first->second!=x[3]) inputError("inconsistent parent for subtype "+x[4]);
            if (grouped) {
                if (x[6]!="neuron" && x[6]!="non_neuron") inputError("neuron_group must be neuron or non_neuron for "+x[0]);
                auto group=cellGroups.emplace(x[4],x[6]=="neuron"?1:0);
                if (!group.second && group.first->second!=(x[6]=="neuron"?1:0))
                    inputError("inconsistent neuron_group for subtype "+x[4]);
            }
        } else if (x[1]=="intercept" || x[1]=="baseline") {
            for (int i=2;i<6;++i) if (x[i]!="NA") inputError("baseline/intercept hierarchy fields must be NA for "+x[0]);
            if (grouped && x[6]!="NA") inputError("baseline/intercept neuron_group must be NA for "+x[0]);
        } else inputError("unknown annotation role "+x[1]);
    }
    if (modalities.empty()) inputError("metadata contains no cell_stage annotations.");
    Input result;
    if (requestedModality.empty()) {
        if (modalities.size()!=1) inputError("single-block reader received multiple modalities; select a block explicitly (the CLI fits multiple blocks jointly).");
        result.modality=*modalities.begin();
    } else {
        if (!modalities.count(requestedModality)) inputError("unknown modality "+requestedModality);
        result.modality=requestedModality;
    }
    std::set<std::string> cells,stages,parents;
    for (const auto& item:rows) {
        const Row& r=item.second;
        if (r.role=="cell_stage" && r.modality==result.modality) {
            cells.insert(r.cell); stages.insert(r.stage); parents.insert(r.parent);
        }
    }
    result.cells.assign(cells.begin(),cells.end());
    result.stages.assign(stages.begin(),stages.end());
    result.parents.assign(parents.begin(),parents.end());
    std::map<std::string,int> ci,ti,pi;
    for (int i=0;i<int(result.cells.size());++i) ci[result.cells[i]]=i;
    for (int i=0;i<int(result.stages.size());++i) ti[result.stages[i]]=i;
    for (int i=0;i<int(result.parents.size());++i) pi[result.parents[i]]=i;
    for (const auto& cell:result.cells) result.parentOfCell.push_back(pi.at(cellParents.at(cell)));
    if (grouped) {
        for (const auto& cell:result.cells) result.neuron.push_back(cellGroups.at(cell));
        result.metadata(); // validate both groups within each selected modality
        std::cout<<"CellDev: explicit Gaussian neuron/non-neuron stage hierarchy; low-rank interaction replaced by constrained parent-stage effects."<<std::endl;
    }
    std::vector<std::string> grid(cells.size()*stages.size());
    for (const auto& item:rows) {
        const Row& r=item.second;
        if (r.role=="cell_stage" && r.modality==result.modality)
            grid[ci.at(r.cell)*stages.size()+ti.at(r.stage)]=item.first;
    }
    int missing=0;
    for (int c=0;c<int(cells.size());++c) for (int t=0;t<int(stages.size());++t) {
        const auto& id=grid[c*stages.size()+t];
        if (id.empty()) { ++missing; continue; }
        result.assayCells.push_back(c); result.assayStages.push_back(t);
        result.cellColumnIds.push_back(id);
    }
    if (missing) std::cout<<"CellDev: "<<result.modality<<" has "<<result.cellColumnIds.size()
        <<" observed assays and "<<missing<<" missing subtype-stage assays; missing assays "
        <<"are excluded from the predictor (no zero filling)."<<std::endl;

    std::ifstream input(annotationFile);
    if (!input) inputError("cannot read annotation matrix "+annotationFile);
    std::cout<<"CellDev: reading annotation matrix for "<<snpIds.size()<<" retained SNPs, modality "<<result.modality<<"."<<std::endl;
    std::getline(input,line);
    auto columns=fields(line);
    if (columns.size()<2 || columns[0]!="SNP") inputError("annotation matrix must start with SNP and have a header.");
    std::map<std::string,int> col;
    int intercept=-1;
    for (int i=1;i<int(columns.size());++i) {
        if (!col.emplace(columns[i],i).second) inputError("duplicate matrix column "+columns[i]);
        auto it=rows.find(columns[i]);
        if (it==rows.end()) inputError("no metadata for annotation "+columns[i]);
        if (includeBaseline && it->second.role=="baseline") result.baselineIds.push_back(columns[i]);
        if (it->second.role=="intercept") {
            if (intercept!=-1) inputError("requires exactly one intercept column.");
            intercept=i;
        }
    }
    if (intercept<0) inputError("requires one intercept column.");
    for (const auto& item:rows) if (!col.count(item.first)) inputError("metadata annotation absent from matrix: "+item.first);
    if (snpIds.empty()) inputError("no retained SNPs.");
    std::unordered_map<std::string,int> snpIndex;
    snpIndex.reserve(snpIds.size());
    for (int i=0;i<int(snpIds.size());++i)
        if (!snpIndex.emplace(snpIds[i],i).second) inputError("duplicate retained SNP "+snpIds[i]);
    result.annotations.resize(snpIds.size(),result.cellColumnIds.size());
    result.baseline.resize(snpIds.size(),result.baselineIds.size());
    std::vector<int> ac,bc;
    for (const auto& name:result.cellColumnIds) ac.push_back(col.at(name));
    for (const auto& name:result.baselineIds) bc.push_back(col.at(name));
    std::unordered_set<std::string> seen;
    seen.reserve(snpIds.size());
    std::vector<bool> matched(snpIds.size(),false);
    size_t ignoredRows=0;
    size_t matchedRows=0;
    while (std::getline(input,line)) {
        auto x=fields(line); if (x.empty()) continue;
        if (x.size()!=columns.size()) inputError("inconsistent matrix row width at SNP "+x[0]);
        if (!seen.insert(x[0]).second) inputError("duplicate SNP row "+x[0]);
        // Numeric annotations are needed only for retained SNPs. In particular,
        // an unmatched metadata row must not be interpreted as missing SNP data.
        auto found=snpIndex.find(x[0]);
        if (found==snpIndex.end()) { ++ignoredRows; continue; }
        std::vector<double> values(x.size());
        for (int k=1;k<int(x.size());++k) values[k]=number(x[k],x[0],columns[k]);
        if (values[intercept]!=1) inputError("intercept must equal one at SNP "+x[0]);
        const int j=found->second; matched[j]=true;
        for (int k=0;k<int(ac.size());++k) result.annotations(j,k)=values[ac[k]];
        for (int k=0;k<int(bc.size());++k) result.baseline(j,k)=values[bc[k]];
        if (!(++matchedRows%1000000))
            std::cout<<"CellDev: read annotations for "<<matchedRows<<" retained SNPs."<<std::endl;
    }
    if (ignoredRows)
        std::cout<<"CellDev: ignored "<<ignoredRows<<" annotation rows whose IDs are not retained SNPs."<<std::endl;
    for (int j=0;j<int(matched.size());++j) if (!matched[j]) inputError("missing retained SNP row "+snpIds[j]);
    if (!result.annotations.allFinite() || !result.baseline.allFinite())
        inputError("annotation values exceed the supported floating-point range.");
    std::map<std::string,EnrichmentColumn> enrichmentInfo;
    auto standardise=[&](MatrixXf& matrix,const std::vector<std::string>& ids,bool baseline) {
        for (int k=0;k<matrix.cols();++k) {
            VectorXd column=matrix.col(k).cast<double>();
            const double mean=column.mean(); column.array()-=mean;
            const double sd=std::sqrt(column.squaredNorm()/column.size());
            const bool binary=(matrix.col(k).array()==0 || matrix.col(k).array()==1).all();
            EnrichmentColumn info;
            info.name=ids[k]; info.binary=binary; info.fraction=binary?mean:0;
            info.offset=binary && sd>0?0:mean;
            info.scale=binary && sd>0?1:sd;
            enrichmentInfo.emplace(info.name,info);
            if (!(sd>0)) {
                matrix.col(k).setZero();
                result.transforms.push_back({ids[k],"constant_zero_design",mean,0});
            } else if (binary)
                result.transforms.push_back({ids[k],"binary_unscaled",0,1});
            else {
                matrix.col(k)=(column/sd).cast<float>();
                result.transforms.push_back({ids[k],"standardised",mean,sd});
            }
        }
    };
    std::cout<<"CellDev: preparing "<<result.annotations.cols()<<" cell-stage and "<<result.baseline.cols()<<" baseline columns."<<std::endl;
    standardise(result.annotations,result.cellColumnIds,false);
    standardise(result.baseline,result.baselineIds,true);
    for (const auto& id:result.baselineIds) result.enrichmentColumns.push_back(enrichmentInfo.at(id));
    for (const auto& id:result.cellColumnIds) result.enrichmentColumns.push_back(enrichmentInfo.at(id));
    return result;
}

Priors Input::calibratedPriors(const VectorXf& pis,int rank) const {
    std::cout<<"CellDev: calibrating fixed annotation prior scales."<<std::endl;
    Metadata m=metadata();
    const VectorXd norms=productColumnSquares(annotations,mainWeights(m));
    Priors p;
    p.linearMean=MatrixXd::Zero(1+baseline.cols()+norms.size(),pis.size()-1);
    p.linearVariance=MatrixXd::Ones(p.linearMean.rows(),pis.size()-1);
    if (pis.size()<2 || (pis.array()<=0).any() || std::abs(pis.sum()-1)>1e-5)
        inputError("requires positive initial mixture probabilities summing to one.");
    for (int k=0;k<pis.size()-1;++k) {
        double prob=pis.tail(pis.size()-k-1).cast<double>().sum()/pis.tail(pis.size()-k).cast<double>().sum();
        p.linearMean(0,k)=boost::math::quantile(boost::math::normal_distribution<double>(),prob);
    }
    int offset=1+baseline.cols(),designOffset=0;
    for (int width:{m.parents,m.cells-m.parents,m.stages-1,
                   m.structured()?m.stages-1:0,m.structured()?int(m.QJ.cols())*(m.stages-1):0}) {
        if (width) {
            double h=norms.segment(designOffset,width).sum()/annotations.rows();
            p.linearVariance.middleRows(offset,width).setConstant(h>0?0.0625/h:1.0);
        }
        offset+=width;
        designOffset+=width;
    }
    double h=0;
    if (!m.structured() && rank>0 && m.cells>1 && m.stages>1) {
        MatrixXd weights(annotations.cols(),(m.cells-1)*(m.stages-1));
        for (int q=0;q<weights.rows();++q)
            for (int u=0;u<m.cells-1;++u) for (int v=0;v<m.stages-1;++v)
                weights(q,u*(m.stages-1)+v)=m.QU(m.assayCells[q],u)*m.QT(m.assayStages[q],v);
        h=productColumnSquares(annotations,weights).sum()/annotations.rows();
    }
    p.interactionVariance=VectorXd::Constant(pis.size()-1,h>0?0.0625/h:1.0);
    return p;
}

void Input::writeManifest(const std::string& prefix,const Priors& prior,int rank,const std::string& parameterPrefix) const {
    Metadata m=metadata();
    if (m.structured()) rank=0;
    const auto labels=coordinateLabels(m,baseline.cols(),prior.linearMean.cols(),rank);
    std::ofstream out(prefix+".celldev.input.txt");
    if (!out) inputError("cannot write input manifest for "+prefix);
    out<<"annotation_id\trole\tmodality\tparent_id\tsubtype_id\tstage_id\tmean\tsd\tstatus\tis_binary\traw_fraction\toutput_name\n"<<std::setprecision(17);
    std::map<std::string,Transform> transform;
    for (const auto& t:transforms) transform.emplace(t.id,t);
    std::map<std::string,EnrichmentColumn> info;
    for (const auto& column:enrichmentColumns) info.emplace(column.name,column);
    auto writeBinaryInfo=[&](const std::string& id,const std::string& outputName) {
        const auto& column=info.at(id);
        out<<'\t'<<column.binary<<'\t';
        if (column.binary) out<<column.fraction; else out<<"NA";
        out<<'\t'<<outputName<<'\n';
    };
    for (int i=0;i<int(cellColumnIds.size());++i) {
        const auto& t=transform.at(cellColumnIds[i]); const int c=m.assayCells[i], stage=m.assayStages[i];
        out<<t.id<<"\tcell_stage\t"<<modality<<'\t'<<parents[parentOfCell[c]]<<'\t'<<cells[c]<<'\t'<<stages[stage]<<'\t'<<t.mean<<'\t'<<t.sd<<'\t'<<t.status;
        writeBinaryInfo(t.id,(joint?modality+":":"")+cells[c]+"_"+stages[stage]);
    }
    for (const auto& id:baselineIds) {
        const auto& t=transform.at(id);
        out<<id<<"\tbaseline\tNA\tNA\tNA\tNA\t"<<t.mean<<'\t'<<t.sd<<'\t'<<t.status;
        writeBinaryInfo(id,id);
    }
    std::ofstream assays(prefix+".celldev.assays.txt");
    if (!assays) inputError("cannot write assay manifest for "+prefix);
    assays<<"modality\tparent_id\tsubtype_id\tstage_id\tannotation_id\tassay_status\toutput_name\n";
    std::map<std::pair<int,int>,std::string> assayed;
    for (int q=0;q<int(cellColumnIds.size());++q)
        assayed[{m.assayCells[q],m.assayStages[q]}]=cellColumnIds[q];
    for (int c=0;c<m.cells;++c) for (int t=0;t<m.stages;++t) {
        auto found=assayed.find({c,t});
        assays<<modality<<'\t'<<parents[parentOfCell[c]]<<'\t'<<cells[c]<<'\t'<<stages[t]<<'\t'
              <<(found==assayed.end()?"NA":found->second)<<'\t'
              <<(found==assayed.end()?"missing_predicted":"observed")<<'\t'
              <<(joint?modality+":":"")<<cells[c]<<'_'<<stages[t]<<'\n';
    }
    std::ofstream scales(prefix+".celldev.priors.txt");
    if (!scales) inputError("cannot write prior manifest for "+prefix);
    if (prior.rd) scales<<"# RD: zero-centred mixture; relative variances 0,.0001,.001,.01,.1 divided by .02222; Dirichlet(1,1,1,1,1); intercept Gaussian\n";
    scales<<"split\tcoordinate\tmean\tvariance\tname\n"<<std::setprecision(17);
    for (int k=0;k<prior.linearMean.cols();++k) {
        for (int i=0;i<prior.linearMean.rows();++i)
            scales<<k+1<<'\t'<<i+1<<'\t'<<prior.linearMean(i,k)<<'\t'<<(i==0 && !parameterPrefix.empty()?0:prior.linearVariance(i,k))<<'\t'<<labels.linear[i]<<'\n';
        if (!m.structured()) scales<<k+1<<"\tinteraction_total\t0\t"<<prior.interactionVariance[k]<<"\tInteractionTotalVarianceBeforeRankDivision\n";
    }
    std::ofstream contrasts(prefix+".celldev.contrasts.txt");
    if (!contrasts) inputError("cannot write contrast manifest for "+prefix);
    auto writeMatrix=[&](const std::string& id,const MatrixXd& matrix,const std::vector<std::string>& rows,
                         const std::vector<std::string>& columns,const std::string& description) {
        contrasts<<"# "<<id<<": "<<description<<"\n# dimensions: "<<matrix.rows()<<" x "<<matrix.cols()<<"\n# rows:";
        for (const auto& row:rows) contrasts<<'\t'<<row;
        contrasts<<"\n# columns:";
        for (const auto& col:columns) contrasts<<'\t'<<col;
        contrasts<<'\n'<<id<<'\n'<<std::setprecision(17)<<matrix<<'\n';
    };
    writeMatrix("P",m.P,cells,parents,"subtype-to-parent membership; 1 means subtype belongs to parent");
    writeMatrix("H",m.H,cells,labels.subtypeContrasts,"within-parent orthonormal Helmert contrasts; P^T H = 0");
    writeMatrix("QG",m.QG,parents,labels.parentContrasts,"parent-factor contrasts; n^T QG = 0, n = subtype counts per parent");
    writeMatrix("QT",m.QT,stages,labels.stageContrasts,"stage Helmert contrasts; 1^T QT = 0");
    auto qu=labels.parentContrasts; qu.insert(qu.end(),labels.subtypeContrasts.begin(),labels.subtypeContrasts.end());
    writeMatrix("QU",m.QU,cells,qu,"centred cell-factor basis [P QG, H]; not generally orthonormal");
    if (m.structured()) {
        std::vector<std::string> interactionLabels;
        for (int r=0;r<m.QJ.cols();++r) interactionLabels.push_back("QJ:"+std::to_string(r+1));
        writeMatrix("QJ",m.QJ,parents,interactionLabels,"residual parent-stage cell basis; n^T QJ = r^T QJ = 0");
        std::ofstream groups(prefix+".celldev.groups.txt");
        if (!groups) inputError("cannot write group manifest for "+prefix);
        groups<<"subtype_id\tneuron_group\tq\n"<<std::setprecision(17);
        for (int c=0;c<m.cells;++c) groups<<cells[c]<<'\t'<<(m.neuron[c]?"neuron":"non_neuron")<<'\t'<<m.groupCoding[c]<<'\n';
    }
    std::ofstream coordinates(prefix+".celldev.coordinates.txt");
    if (!coordinates) inputError("cannot write coordinate dictionary for "+prefix);
    coordinates<<"parameter\tindex\tlegacy_name\tname\tmeaning\n";
    auto writeCoordinates=[&](const std::string& parameter,const std::vector<std::string>& names,const std::string& meaning) {
        for (int i=0;i<int(names.size());++i)
            coordinates<<parameterPrefix+parameter<<'\t'<<i+1<<"\tcoordinate"<<i+1<<'\t'<<names[i]<<'\t'<<meaning<<'\n';
    };
    for (int k=0;k<prior.linearMean.cols();++k)
        writeCoordinates("CellDevLinear_p"+std::to_string(k+2),labels.linear,
            "split "+std::to_string(k+1)+": Pr(d >= "+std::to_string(k+1)+" | d >= "+std::to_string(k)+"); intercept/baseline/parent/H/QT coordinates");
    writeCoordinates("CellDevA",labels.cellFactor,"a coordinates; u = QU a; factor-major order");
    writeCoordinates("CellDevF",labels.stageFactor,"f coordinates; v = QT f; factor-major order");
    writeCoordinates("CellDevLambda",labels.strengths,"split-specific interaction strength; split-major order");
}
}

namespace {
const Data& checkedEigenData(const Data& d, const VectorXf& pis,
                            const VectorXf& gamma, const MatrixXf& a,
                            const MatrixXf& b) {
    auto require=[](bool ok, const char* message) {
        if (!ok) throw std::invalid_argument(message);
    };
    require(d.lowRankModel && !d.sparseLDM && !d.coreContextLdmEigen
            && !d.integrateOutContext
            && d.coreContextLdBlocks.empty(),
            "CellDev requires ordinary eigen LD; quantized, sparse and tiled LD are unsupported.");
    // Ordinary eigen loading also allocates inactive quantization placeholders.
    for (const QuantizedEigenBlock& block:d.quantizedEigenBlocks)
        require(!block.active(),"CellDev does not support quantized eigen LD.");
    require(!d.numFixedEffectSnps && !d.removeSkepticalSnps,
            "Core CellDev does not support fixed-effect SNPs or skeptical-SNP removal.");
    require(pis.size()>1 && gamma.size()==pis.size() && pis.allFinite()
            && gamma.allFinite() && (pis.array()>0).all()
            && std::abs(pis.sum()-1)<1e-5 && gamma[0]==0
            && (gamma.tail(gamma.size()-1).array()>0).all(),
            "Invalid CellDev SNP mixture probabilities or variances.");
    for (int k=2;k<gamma.size();++k)
        require(gamma[k]>gamma[k-1],"CellDev mixture variances must be increasing.");
    require(a.rows()==d.numIncdSnps && b.rows()==d.numIncdSnps,
            "CellDev annotations must follow retained SNP order.");
    const unsigned blocks=d.keptLdBlockInfoVec.size();
    require(blocks && d.Qblocks.size()==blocks && d.wcorrBlocks.size()==blocks
            && d.nGWASblock.size()==blocks && d.numEigenvalBlock.size()==blocks,
            "Incomplete CellDev eigen LD data.");
    unsigned next=0;
    for (unsigned i=0;i<blocks;++i) {
        const LDBlockInfo* block=d.keptLdBlockInfoVec[i];
        require(block && !block->hasCoreContextTileMetadata() && block->eigenColRemap.size()==0,
                "CellDev requires full ordinary eigen blocks without tile metadata or SNP remapping.");
        require(block->startSnpIdx==int(next) && block->endSnpIdx>=block->startSnpIdx,
                "CellDev eigen blocks must cover retained SNPs once in order.");
        const int n=block->endSnpIdx-block->startSnpIdx+1;
        const bool valid=d.Qblocks[i].cols()==n && d.Qblocks[i].rows()>0
                && d.Qblocks[i].rows()==d.wcorrBlocks[i].size()
                && d.Qblocks[i].rows()==d.numEigenvalBlock[i]
                && block->qDiag.size()==n && block->snpInfoVec.size()>=unsigned(n)
                && d.Qblocks[i].allFinite() && d.wcorrBlocks[i].allFinite()
                && block->qDiag.allFinite() && (block->qDiag.array()>0).all()
                && d.nGWASblock[i]>0;
        if (!valid) {
            std::ostringstream detail;
            detail<<"Error: invalid CellDev eigen block "<<block->ID
                  <<": SNPs="<<n<<", Q="<<d.Qblocks[i].rows()<<"x"<<d.Qblocks[i].cols()
                  <<", wcorr="<<d.wcorrBlocks[i].size()<<", rank="<<d.numEigenvalBlock[i]
                  <<", qDiag="<<block->qDiag.size()<<", SNP map="<<block->snpInfoVec.size()
                  <<", N="<<d.nGWASblock[i]<<".";
            throw std::invalid_argument(detail.str());
        }
        next+=n;
    }
    require(next==d.numIncdSnps && d.incdSnpInfoVec.size()==next,
            "CellDev eigen block coverage does not match retained SNPs.");
    for (const LDBlockInfo* block:d.keptLdBlockInfoVec)
        for (int i=0;i<=block->endSnpIdx-block->startSnpIdx;++i)
            require(block->snpInfoVec[i]==d.incdSnpInfoVec[block->startSnpIdx+i],
                    "CellDev eigen SNP maps must match retained SNP order.");
    for (const SnpInfo* snp:d.incdSnpInfoVec)
        require(snp && !snp->skip && !snp->fitAsFixedEffect,
                "Core CellDev requires retained random-effect SNPs without skipped entries.");
    return d;
}
}

ApproxBayesRCellDev::ApproxBayesRCellDev(
        const Data& d, const VectorXf& pis, const VectorXf& piPar, const VectorXf& gamma,
        bool noscale, bool hsqPercModel, bool robustMode, const CellDev::Metadata& metadata,
        const MatrixXf& a, const MatrixXf& b, int rank, const CellDev::Priors& priors,
        const ApproxBayesRCellDev* sharedDesigns)
    : ApproxBayesR(checkedEigenData(d,pis,gamma,a,b),true,d.varGenotypic,d.varResidual,
                   pis,piPar,gamma,true,noscale,hsqPercModel,robustMode,"Gibbs",false),
      snpEffects(d.snpEffectNames,pis), annotation(metadata,a,b,pis.size()-1,rank,priors,sharedDesigns?&sharedDesigns->annotation:nullptr),
      enrichment(annotation) {
    const char* profile=std::getenv("GCTB_PROFILE_CELLDEV");
    profileTiming=profile && std::string(profile)!="0";
    bayesType="RCellDev";
    paramSetVec={&snpEffects,&snpPip,&snpHsqPep,&vargBlk,&vareBlk};
    paramSetVec.insert(paramSetVec.end(),deltaPi.begin(),deltaPi.end());
    paramSetToPrint.clear();
    reports.reset(new CellDev::Reports(annotation,enrichment));
    reports->registerSets(paramSetVec,paramSetToPrint);
    paramSetVec.push_back(&enrichment.marginal);
    paramSetVec.push_back(&enrichment.conditional);
    paramSetToPrint.push_back(&enrichment.marginal);
    paramSetToPrint.push_back(&enrichment.conditional);
    // Bad-SNP removal is outside the fixed model target in this core version.
    paramToPrint.erase(std::remove(paramToPrint.begin(),paramToPrint.end(),&nBadSnps),paramToPrint.end());
    annotation.logProbabilities(snpLogPi);
    collectAnnotationParameters();
}

void ApproxBayesRCellDev::collectAnnotationParameters() {
    reports->collect(annotation);
    for (auto& block:otherModalities) block.reports->collect(*block.sampler);
}

CellDev::Reports::Reports(const Sampler& annotation, const Enrichment& enrichment, const std::string& prefix) {
    // Register free-coordinate draws using the existing MCMC output machinery.
    annotationHeaders.reserve(annotation.splits+3);
    annotationParameters.reserve(annotation.splits+3);
    const auto labels=CellDev::coordinateLabels(annotation.metadata,annotation.baselines(),annotation.splits,annotation.rank);
    for (int k=0;k<annotation.splits;++k) {
        annotationHeaders.push_back(labels.linear);
        annotationParameters.emplace_back("CellDevLinear_p"+std::to_string(k+2),annotationHeaders.back());
    }
    for (const auto& item:std::vector<std::pair<std::string,std::vector<std::string> > >{
            {"CellDevA",labels.cellFactor},{"CellDevF",labels.stageFactor},
            {"CellDevLambda",labels.strengths}}) {
        annotationHeaders.push_back(item.second);
        annotationParameters.emplace_back(item.first,annotationHeaders.back());
    }
    // Publish reconstructed biological quantities in .parSetRes. Free contrast
    // coordinates remain registered for draw-level traces and the dictionary.
    if (annotation.priors.rd) {
        rdHeaders.reserve(2*annotation.splits); rdParameters.reserve(2*annotation.splits);
        for (int k=0;k<annotation.splits;++k) {
            rdHeaders.emplace_back(labels.linear.begin()+1,labels.linear.end());
            rdParameters.emplace_back(prefix+"CellDevRDInclusion_p"+std::to_string(k+2),rdHeaders.back());
            rdHeaders.push_back({"zero","small","medium","large","largest"});
            rdParameters.emplace_back(prefix+"CellDevRDWeights_p"+std::to_string(k+2),rdHeaders.back());
        }
    }
    const int sets=6*annotation.splits+2*annotation.rank;
    biologicalHeaders.reserve(sets); biologicalParameters.reserve(sets);
    auto entities=[&](const std::vector<std::string>& ids,int count,const std::string& prefix) {
        if (int(ids.size())==count) {
            auto result=ids;
            if (!annotation.metadata.modality.empty())
                for (auto& name:result) name=annotation.metadata.modality+":"+name;
            return result;
        }
        std::vector<std::string> names;
        for (int i=0;i<count;++i) names.push_back(prefix+std::to_string(i+1));
        return names;
    };
    const auto cells=entities(annotation.metadata.cellNames,annotation.metadata.cells,"cell");
    const auto stages=entities(annotation.metadata.stageNames,annotation.metadata.stages,"stage");
    const auto parents=entities(annotation.metadata.parentNames,annotation.metadata.parents,"parent");
    std::vector<std::string> factors;
    for (int r=0;r<annotation.rank;++r) factors.push_back("factor"+std::to_string(r+1));
    auto addReport=[&](const std::string& label,const std::vector<std::string>& names) {
        biologicalHeaders.push_back(names);
        biologicalParameters.emplace_back(label,biologicalHeaders.back());
    };
    const std::vector<std::string> cellStages(enrichment.names.begin()+1+annotation.baselines(),enrichment.names.end());
    for (int k=0;k<annotation.splits;++k) {
        const auto suffix="_p"+std::to_string(k+2);
        addReport("CellDevEffects"+suffix,enrichment.names);
        addReport("CellDevParent"+suffix,parents);
        addReport("CellDevSubtype"+suffix,cells);
        addReport("CellDevStage"+suffix,stages);
        addReport("CellDevInteraction"+suffix,cellStages);
        addReport("CellDevStrength"+suffix,factors);
    }
    for (int r=0;r<annotation.rank;++r) {
        addReport("CellDevCellLoadings_factor"+std::to_string(r+1),cells);
        addReport("CellDevStageProfile_factor"+std::to_string(r+1),stages);
    }
    const auto& m=annotation.metadata;
    std::set<int> assayed;
    for (int q=0;q<int(m.assayCells.size());++q) assayed.insert(m.assayCells[q]*m.stages+m.assayStages[q]);
    std::vector<std::string> predicted;
    for (int c=0;c<m.cells;++c) for (int t=0;t<m.stages;++t)
        if (!assayed.count(c*m.stages+t)) {
            missingGridRows.push_back(c*m.stages+t);
            const std::string cell=c<int(m.cellNames.size())?m.cellNames[c]:"cell"+std::to_string(c+1);
            const std::string stage=t<int(m.stageNames.size())?m.stageNames[t]:"stage"+std::to_string(t+1);
            predicted.push_back((m.modality.empty()?"":m.modality+":")+cell+"_"+stage);
        }
    predictionHeaders.reserve(annotation.splits); predictionParameters.reserve(annotation.splits);
    if (!predicted.empty()) for (int k=0;k<annotation.splits;++k) {
        predictionHeaders.push_back(predicted);
        predictionParameters.emplace_back(prefix+"CellDevPredictedEffects_p"+std::to_string(k+2),predictionHeaders.back());
    }
    for (auto& p:annotationParameters) p.label=prefix+p.label;
    for (auto& p:biologicalParameters) p.label=prefix+p.label;
    if (annotation.metadata.structured()) {
        componentHeaders.reserve(annotation.splits); componentParameters.reserve(annotation.splits);
        for (int k=0;k<annotation.splits;++k) {
            componentHeaders.push_back(annotation.componentNames());
            componentParameters.emplace_back(prefix+"CellDevComponents_p"+std::to_string(k+2),componentHeaders.back());
        }
    }
}

void CellDev::Reports::registerSets(std::vector<ParamSet*>& all, std::vector<ParamSet*>& printed) {
    for (auto& p:annotationParameters) if (p.size) all.push_back(&p);
    for (auto& p:biologicalParameters) if (p.size) { all.push_back(&p); printed.push_back(&p); }
    for (auto& p:predictionParameters) { all.push_back(&p); printed.push_back(&p); }
    for (auto& p:rdParameters) { all.push_back(&p); printed.push_back(&p); }
    for (auto& p:componentParameters) { all.push_back(&p); printed.push_back(&p); }
}

void CellDev::Reports::collect(const Sampler& annotation) {
    if (annotation.priors.rd) for (int k=0;k<annotation.splits;++k) {
        rdParameters[2*k].values=annotation.rdInclusion.col(k).tail(annotation.linear.rows()-1).cast<float>();
        rdParameters[2*k+1].values=annotation.rdWeights.row(k).transpose().cast<float>();
    }
    for (int k=0;k<annotation.splits;++k)
        annotationParameters[k].values=annotation.linear.col(k).cast<float>();
    const MatrixXd* matrices[]={&annotation.a,&annotation.f,&annotation.lambda};
    for (int i=0;i<3;++i) if (matrices[i]->size())
        annotationParameters[annotation.splits+i].values=
            Eigen::Map<const VectorXd>(matrices[i]->data(),matrices[i]->size()).cast<float>();
    const auto& m=annotation.metadata;
    const int nb=annotation.baselines(), parentStart=1+nb;
    const MatrixXd parent=annotation.parentEffects();
    const MatrixXd subtype=annotation.subtypeEffects();
    const MatrixXd stage=annotation.stageEffects();
    const MatrixXd coefficients=annotation.coefficients(),u=m.QU*annotation.a,v=m.QT*annotation.f;
    const MatrixXd parentInteraction=annotation.parentInteractions();
    if (!componentParameters.empty()) {
        const MatrixXd components=annotation.componentContrasts();
        for (int k=0;k<annotation.splits;++k) componentParameters[k].values=components.col(k).cast<float>();
    }
    for (int k=0;k<annotation.splits;++k) {
        auto& effect=biologicalParameters[6*k].values;
        effect.head(1+nb)=annotation.linear.col(k).head(1+nb).cast<float>();
        for (int q=0;q<int(m.assayCells.size());++q)
            effect[1+nb+q]=coefficients(m.assayCells[q]*m.stages+m.assayStages[q],k);
        biologicalParameters[6*k+1].values=parent.col(k).cast<float>();
        biologicalParameters[6*k+2].values=subtype.col(k).cast<float>();
        biologicalParameters[6*k+3].values=stage.col(k).cast<float>();
        auto& interaction=biologicalParameters[6*k+4].values;
        interaction.setZero();
        for (int q=0;q<int(m.assayCells.size());++q) {
            const int c=m.assayCells[q], t=m.assayStages[q];
            double value=0;
            for (int r=0;r<annotation.rank;++r) value+=u(c,r)*v(t,r)*annotation.lambda(r,k);
            if (m.structured()) value=parentInteraction(m.parent[c]*m.stages+t,k);
            interaction[q]=value;
        }
        biologicalParameters[6*k+5].values=annotation.lambda.col(k).cast<float>();
    }
    for (int k=0;k<int(predictionParameters.size());++k)
        for (int q=0;q<int(missingGridRows.size());++q)
            predictionParameters[k].values[q]=coefficients(missingGridRows[q],k);
    for (int r=0;r<annotation.rank;++r) {
        biologicalParameters[6*annotation.splits+2*r].values=u.col(r).cast<float>();
        biologicalParameters[6*annotation.splits+2*r+1].values=v.col(r).cast<float>();
    }
}

void ApproxBayesRCellDev::addModality(const CellDev::Input& input, const CellDev::Priors& priors, int rank,
                                    const ApproxBayesRCellDev* sharedDesigns) {
    if (input.baseline.cols()!=0)
        throw std::invalid_argument("Additional CellDev modalities must omit the shared baseline.");
    ModalityBlock block;
    block.sampler.reset(new CellDev::Sampler(input.metadata(),input.annotations,input.baseline,
        annotation.splits,rank,priors,sharedDesigns?
            sharedDesigns->otherModalities.at(otherModalities.size()).sampler.get():nullptr));
    block.sampler->linear.row(0).setZero();
    block.sampler->rebuildEta();
    block.enrichment.reset(new CellDev::Enrichment(*block.sampler));
    block.enrichment->marginal.label=input.modality+"_"+block.enrichment->marginal.label;
    block.enrichment->conditional.label=input.modality+"_"+block.enrichment->conditional.label;
    block.reports.reset(new CellDev::Reports(*block.sampler,*block.enrichment,input.modality+"_"));
    block.reports->registerSets(paramSetVec,paramSetToPrint);
    paramSetVec.push_back(&block.enrichment->marginal); paramSetVec.push_back(&block.enrichment->conditional);
    paramSetToPrint.push_back(&block.enrichment->marginal); paramSetToPrint.push_back(&block.enrichment->conditional);
    otherModalities.push_back(std::move(block));
    collectAnnotationParameters();
}

MatrixXd ApproxBayesRCellDev::jointPredictor() const {
    MatrixXd result=annotation.eta;
    for (const auto& block:otherModalities) result+=block.sampler->eta;
    return result;
}

void ApproxBayesRCellDev::sampleAnnotations() {
    if (otherModalities.empty()) { annotation.sweep(snpEffects.membership); return; }
    const double start=CellDev::wallSeconds();
    MatrixXd total=jointPredictor();
    annotation.eta.swap(total);
    annotation.sampleZ(snpEffects.membership);
    annotation.eta.swap(total);
    // Keep one immutable shared augmentation throughout the block Gibbs sweep.
    const MatrixXd sharedZ=annotation.z;
    const double latent=CellDev::wallSeconds();
    total-=annotation.eta;
    annotation.sweepConditional(annotation,total,true);
    total+=annotation.eta;
    annotation.z=sharedZ;
    for (auto& block:otherModalities) {
        total-=block.sampler->eta;
        block.sampler->sweepConditional(annotation,total,false);
        total+=block.sampler->eta;
        annotation.timing.linear+=block.sampler->timing.linear;
        annotation.timing.factors+=block.sampler->timing.factors;
        annotation.timing.lambda+=block.sampler->timing.lambda;
    }
    annotation.timing.z=latent-start;
}

void ApproxBayesRCellDev::sampleStartVal() {
    throw std::logic_error("Construct an independent CellDev model for each fresh chain.");
}

void ApproxBayesRCellDev::sampleUnknowns(unsigned iter) {
    const double start=CellDev::wallSeconds();
    annotation.invalidateZ();
    if (otherModalities.empty()) annotation.logProbabilities(snpLogPi);
    else {
        MatrixXd total=jointPredictor();
        annotation.eta.swap(total);
        annotation.logProbabilities(snpLogPi);
        annotation.eta.swap(total);
    }
    const double logEnd=CellDev::wallSeconds();
    snpEffects.sampleFromFC_eigen(wcorrBlocks,data.Qblocks,whatBlocks,data.keptLdBlockInfoVec,
        data.nGWASblock,vareBlk.values,snpLogPi,gamma.values,vargRandom.value,deltaPi,
        hsqPercModel,sigmaSq.value,data.incdSnpInfoVec,data.numFixedEffectSnps);
    const double snpEnd=CellDev::wallSeconds();
    snpEffects.computePosteriorMean(iter);
    snpPip.getValues(snpEffects.pip);
    nnzSnp.getValue(snpEffects.numNonZeros);
    numSnps.getValues(snpEffects.numSnpMix);
    if (robustMode)
        sigmaSq.value=vargRandom.value/(data.numIncdSnps*gamma.values.dot(Pis.values));
    else sigmaSq.sampleFromFC(snpEffects.wtdSumSq,snpEffects.numNonZeros);
    vargBlk.compute(whatBlocks);
    vareBlk.sampleFromFC(wcorrBlocks,vargBlk.values,snpEffects.ssqBlocks,
                        data.nGWASblock,data.numEigenvalBlock);
    vargRandom.value=varg.value=vargBlk.total;
    vare.value=vareBlk.mean;
    hsqRandom.compute(vargRandom.value,data.varPhenotypic);
    hsq.value=hsqRandom.value;
    const double varianceEnd=CellDev::wallSeconds();
    sampleAnnotations();
    Vgs.compute(snpEffects.values,snpEffects.snpset);
    if (summaryRequired) snpHsqPep.compute(snpEffects.values,vargRandom.value);
    if (!(iter%100)) {
        annotation.rebuildEta();
        for (auto& block:otherModalities) block.sampler->rebuildEta();
        rounding.computeWcorr_eigen(data.wcorrBlocks,data.Qblocks,data.keptLdBlockInfoVec,
                                    snpEffects.values,wcorrBlocks);
    }
    const double enrichmentStart=CellDev::wallSeconds();
    if (summaryRequired) {
        collectAnnotationParameters();
        enrichment.compute(annotation,snpEffects.values,gamma.values);
        if (!otherModalities.empty()) {
            // File-facing binary inputs are raw; continuous inputs have mean zero.
            const VectorXd reference=annotation.linear.row(0).transpose();
            for (auto& block:otherModalities)
                block.enrichment->compute(*block.sampler,snpEffects.values,gamma.values,&reference);
        }
    }
    const double enrichmentEnd=CellDev::wallSeconds();
    if (profileTiming || iter==1 || !(iter%10)) {
        std::ostringstream line;
        line<<std::fixed<<std::setprecision(3)<<"CellDev ";
        if (chainIndex>=0) line<<"chain "<<chainIndex+1<<" ";
        line<<"iteration "<<iter
            <<": total="<<CellDev::wallSeconds()-start<<"s, logPi="<<logEnd-start
            <<"s, SNP="<<snpEnd-logEnd<<"s, variance/summaries="<<varianceEnd-snpEnd
            <<"s, z="<<annotation.timing.z<<"s, linear="<<annotation.timing.linear
            <<"s, factors="<<annotation.timing.factors<<"s, lambda="<<annotation.timing.lambda
            <<"s, enrichment="<<enrichmentEnd-enrichmentStart<<"s";
        if (iter==1) line<<" (includes first-split Gram cache construction)";
        std::cout<<line.str()<<std::endl;
    }
}

MultiChainCellDev::MultiChainCellDev(std::vector<std::unique_ptr<ApproxBayesRCellDev> > models)
    : chains(std::move(models)) {
    if (chains.size()<2) throw std::invalid_argument("Multi-chain CellDev requires at least two chains.");
    bayesType="RCellDev";
    const unsigned count=chains.size();
    for (unsigned c=0;c<count;++c) chains[c]->chainIndex=c;
    outerThreads=std::min(count,unsigned(omp_get_max_threads()));
    innerThreads=std::max(1u,unsigned(omp_get_max_threads())/outerThreads);
    omp_set_max_active_levels(2);
    const Model& first=*chains.front();
    std::map<const Parameter*,Parameter*> scalarMap;
    std::vector<Parameter*> originals=first.paramVec;
    for (auto* p:first.paramToPrint)
        if (std::find(originals.begin(),originals.end(),p)==originals.end()) originals.push_back(p);
    for (auto* original:originals) {
        Scalar item;
        item.pooled.reset(new Parameter(original->label));
        item.pooled->numChains=count; item.pooled->perChainValue.setZero(count);
        for (const auto& chain:chains) {
            Parameter* match=nullptr;
            for (auto* p:chain->paramVec) if (p->label==original->label) match=p;
            if (!match) for (auto* p:chain->paramToPrint) if (p->label==original->label) match=p;
            if (!match) throw std::logic_error("Inconsistent CellDev chain scalar registration.");
            item.source.push_back(match);
        }
        scalarMap[original]=item.pooled.get();
        paramVec.push_back(item.pooled.get());
        scalars.push_back(std::move(item));
    }
    for (auto* p:first.paramToPrint) paramToPrint.push_back(scalarMap.at(p));
    for (auto* original:first.paramSetVec) {
        Set item;
        item.pooled.reset(new ParamSet(original->label,original->header));
        item.pooled->numChains=count;
        item.pooled->perChainValues.setZero(original->size,count);
        for (const auto& chain:chains) {
            ParamSet* match=nullptr;
            for (auto* p:chain->paramSetVec) if (p->label==original->label) match=p;
            if (!match || match->header!=original->header)
                throw std::logic_error("Inconsistent CellDev chain parameter-set registration.");
            item.source.push_back(match);
        }
        paramSetVec.push_back(item.pooled.get());
        const bool factorSpecific=original->label.find("CellDevStrength")!=std::string::npos ||
            original->label.find("CellDevCellLoadings")!=std::string::npos ||
            original->label.find("CellDevStageProfile")!=std::string::npos;
        if (!factorSpecific && std::find(first.paramSetToPrint.begin(),first.paramSetToPrint.end(),original)
            !=first.paramSetToPrint.end()) paramSetToPrint.push_back(item.pooled.get());
        sets.push_back(std::move(item));
    }
    collect();
    std::cout<<"Multi-chain CellDev: "<<count<<" independent chains; "<<outerThreads
             <<" chains in parallel, "<<innerThreads<<" threads per chain.\n"
             <<"Factor-specific loadings and strengths are retained in chain-labelled traces; "
             <<"pooled summaries report reconstructed effects and total interactions.\n";
}

void MultiChainCellDev::collect() {
    for (auto& item:scalars) {
        for (unsigned c=0;c<chains.size();++c) item.pooled->perChainValue[c]=item.source[c]->value;
        item.pooled->value=item.pooled->perChainValue.mean();
    }
    for (auto& item:sets) {
        for (unsigned c=0;c<chains.size();++c) item.pooled->perChainValues.col(c)=item.source[c]->values;
        item.pooled->values=item.pooled->perChainValues.rowwise().mean();
    }
}

void MultiChainCellDev::sampleUnknowns(unsigned iter) {
    Stat::prepareParallelChainRng(outerThreads);
    std::exception_ptr failure;
    #pragma omp parallel for num_threads(outerThreads) schedule(static)
    for (unsigned c=0;c<chains.size();++c) {
        omp_set_num_threads(innerThreads);
        Stat::McmcChainRngScope rng(c);
        try {
            chains[c]->summaryRequired = summaryRequired;
            chains[c]->sampleUnknowns(iter);
        }
        catch (...) {
            #pragma omp critical(gctb_celldev_chain_exception)
            { if (!failure) failure=std::current_exception(); }
        }
    }
    Stat::finishParallelChainRng();
    if (failure) std::rethrow_exception(failure);
    collect();
}

void MultiChainCellDev::sampleStartVal() {
    throw std::logic_error("Construct independent CellDev models for a fresh multi-chain run.");
}
