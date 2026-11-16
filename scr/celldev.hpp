#ifndef GCTB_CELLDEV_HPP
#define GCTB_CELLDEV_HPP

#include <Eigen/Dense>
#include "model.hpp"
#include <vector>
#include <memory>

// Per-modality CellDev kernels for single or joint fits. Inputs are already processed and in
// retained-SNP order. Cell-stage columns follow the observed assay mapping;
// baseline excludes intercept. Default metadata retains complete c*T+t order.
// File parsing, standardisation, output and checkpoint policy belong to callers.
namespace CellDev {

struct EnrichmentColumn {
    std::string name;
    bool binary = false;
    // Raw binary indicator = offset + scale * fitted annotation column.
    double offset = 0, scale = 1, fraction = 0;
};

struct Metadata {
    int cells, stages, parents;
    std::string modality; // nonempty for modality-qualified joint output labels
    std::vector<int> parent;
    // Annotation-column mapping; only assayed subtype/stage pairs appear here.
    std::vector<int> assayCells, assayStages;
    std::vector<std::string> cellNames, parentNames, stageNames, baselineNames;
    std::vector<EnrichmentColumn> enrichmentColumns; // baseline, then cell-stage
    Eigen::MatrixXd P, H, QG, QT, QU;
    // Optional explicit Gaussian stage hierarchy; labels are per subtype.
    std::vector<int> neuron; // 1 neuron, 0 non-neuron; empty selects legacy model
    Eigen::VectorXd groupCoding, parentGroupCoding;
    Eigen::MatrixXd QJ; // orthonormal residual parent-interaction cell basis
    Metadata(const std::vector<int>& parentOfCell, int numStages);
    void setNeuronGroups(const std::vector<int>& groups);
    bool structured() const { return !neuron.empty(); }
    int extraColumns() const { return structured()?(stages-1)*(1+QJ.cols()):0; }
    int mainColumns() const { return cells+stages-1+extraColumns(); }
};

struct CoordinateLabels {
    std::vector<std::string> linear, cellFactor, stageFactor, strengths;
    std::vector<std::string> subtypeContrasts, parentContrasts, stageContrasts;
};
CoordinateLabels coordinateLabels(const Metadata& metadata, int baselines, int splits, int rank);

struct Priors {
    bool rd = false;
    // Rows: intercept, baseline, parent, subtype contrasts, stage contrasts.
    // Columns: continuation splits. All variances must be finite and positive.
    Eigen::MatrixXd linearMean, linearVariance;
    // Total strength variance before division by rank, one per split.
    Eigen::VectorXd interactionVariance;
    double subtypeFactorSD = 1.0;
};

// File-facing preparation for one modality block. No stage-order file.
struct Input {
    struct Transform { std::string id, status; double mean, sd; };
    std::string modality;
    static std::vector<std::string> selectModalities(const std::string& mapFile, const std::string& selection);
    bool joint = false;
    std::vector<std::string> cells, stages, parents, cellColumnIds, baselineIds;
    std::vector<int> parentOfCell, assayCells, assayStages;
    std::vector<int> neuron;
    std::vector<Transform> transforms;
    std::vector<EnrichmentColumn> enrichmentColumns;
    Eigen::MatrixXf annotations, baseline;
    Metadata metadata() const {
        Metadata m(parentOfCell,stages.size());
        if (!neuron.empty()) m.setNeuronGroups(neuron);
        if (joint) m.modality=modality;
        m.cellNames=cells; m.parentNames=parents; m.stageNames=stages; m.baselineNames=baselineIds;
        m.enrichmentColumns=enrichmentColumns;
        if (!assayCells.empty()) { m.assayCells=assayCells; m.assayStages=assayStages; }
        return m;
    }
    static Input read(const std::string& annotationFile, const std::string& mapFile,
                      const std::vector<std::string>& snpIds, const std::string& modality = "", bool includeBaseline = true);
    Priors calibratedPriors(const Eigen::VectorXf& pis, int rank) const;
    void writeManifest(const std::string& prefix, const Priors& priors, int rank = 0,
                       const std::string& parameterPrefix = "") const;
};

double logPhi(double x);
double truncatedProbit(double mean, bool success);
Eigen::VectorXd gaussianDraw(const Eigen::MatrixXd& precision,
                             const Eigen::VectorXd& information);

class Sampler {
    std::shared_ptr<const Eigen::MatrixXd> annotationStorage, designStorage;
public:
    const Metadata metadata;
    const int splits, rank;
    const Eigen::MatrixXd &annotations, &mainDesign;
    const Priors priors;
    // Free coordinates and caches are exposed for posterior collection and tests.
    Eigen::MatrixXd linear, a, f, lambda, eta, z;
    std::vector<std::vector<int> > eligible;
    Eigen::MatrixXd rdWeights, rdInclusion; // split x mixture; coordinate x split
    Eigen::MatrixXi rdAllocation;
    static Eigen::VectorXd rdProbabilities(double diagonal, double score, double scale,
                                            const Eigen::VectorXd& weights);

    Sampler(const Metadata& metadata, const Eigen::MatrixXf& annotations,
            const Eigen::MatrixXf& baseline, int splits, int rank,
            const Priors& priors, const Sampler* sharedDesigns = nullptr);
    static Eigen::MatrixXd buildMainDesign(const Metadata& metadata,
                                          const Eigen::MatrixXf& annotations,
                                          const Eigen::MatrixXf& baseline);
    Eigen::MatrixXd coefficients() const; // full (c*T+t) grid, including predictions
    int baselines() const { return mainDesign.cols()-1-metadata.mainColumns(); }
    Eigen::MatrixXd parentEffects() const;
    Eigen::MatrixXd subtypeEffects() const;
    Eigen::MatrixXd stageEffects() const;
    Eigen::MatrixXd stageDifference() const;
    Eigen::MatrixXd parentInteractions() const; // full parent x stage grid
    Eigen::MatrixXd componentContrasts() const; // rows match componentNames()
    std::vector<std::string> componentNames() const;
    Eigen::MatrixXd observedCoefficients() const; // annotation-column order
    Eigen::MatrixXd predictor() const;
    void rebuildEta();
    Eigen::MatrixXd logProbabilities() const;
    void logProbabilities(Eigen::MatrixXf& result) const;
    Eigen::MatrixXd cellDesign(int factor) const; // without split strength
    Eigen::MatrixXd stageDesign(int factor) const;
    const Eigen::MatrixXd& factorScores() const;
    struct Timing { double z=0, linear=0, factors=0, lambda=0; } timing;
    void sampleZ(const Eigen::VectorXi& membership);
    void sampleLinear(int split);
    void sampleLinearRD(int split);
    void sampleA(int factor);
    void sampleF(int factor);
    void sampleLambda(int split);
    // Called only AFTER a collapsed SNP allocation sweep. Refreshes all z.
    void sweep(const Eigen::VectorXi& membership);
    void invalidateZ();
    // Conditional block update with one shared probit augmentation.
    void sweepConditional(const Sampler& augmentation, const Eigen::MatrixXd& offset,
                          bool fitIntercept);


private:
    bool augmented = false;
    bool fitIntercept = true;
    unsigned augmentationSweep = 0;
    Eigen::MatrixXd firstSplitGram;
    mutable Eigen::MatrixXd scoresCache, scoresA, scoresF;
    void requireAugmentation() const;
    void updateSharedFactor(int factor, bool cell);
    void updateSplit(int split, const Eigen::MatrixXd& design,
                     Eigen::VectorXd& value, const Eigen::VectorXd& mean,
                     const Eigen::VectorXd& variance, bool linearBlock = false);
};

class Enrichment {
public:
    const std::vector<std::string> names, conditionalNames;
    ParamSet marginal, conditional;
    explicit Enrichment(const Sampler& sampler);
    void compute(const Sampler& sampler, const Eigen::VectorXf& effects,
                 const Eigen::VectorXf& gamma, const Eigen::VectorXd* reference = nullptr);
private:
    std::vector<EnrichmentColumn> columns;
    Eigen::VectorXd means, ssq, sd;
};
// Reusable posterior registration for a modality block.
class Reports {
public:
    std::vector<std::vector<std::string> > annotationHeaders, biologicalHeaders;
    std::vector<ParamSet> annotationParameters, biologicalParameters;
    std::vector<std::vector<std::string> > predictionHeaders;
    std::vector<ParamSet> predictionParameters;
    std::vector<std::vector<std::string> > rdHeaders;
    std::vector<ParamSet> rdParameters;
    std::vector<std::vector<std::string> > componentHeaders;
    std::vector<ParamSet> componentParameters;
    std::vector<int> missingGridRows;
    Reports(const Sampler& annotation, const Enrichment& enrichment, const std::string& prefix="");
    void collect(const Sampler& annotation);
    void registerSets(std::vector<ParamSet*>& all, std::vector<ParamSet*>& printed);
};
}

// Programmatic core model: ordinary, complete, non-overlapping eigen LD only.
class ApproxBayesRCellDev : public ApproxBayesR {
public:
    ApproxBayesRC::SnpEffects snpEffects;
    CellDev::Sampler annotation;
    CellDev::Enrichment enrichment;
    MatrixXf snpLogPi;
    int chainIndex = -1; // set by the multi-chain wrapper for progress messages

    ApproxBayesRCellDev(const Data& data, const VectorXf& pis,
                       const VectorXf& piPar, const VectorXf& gamma,
                       bool noscale, bool hsqPercModel, bool robustMode,
                       const CellDev::Metadata& metadata,
                       const MatrixXf& cellStageAnnotations,
                       const MatrixXf& baselineAnnotations,
                       int rank, const CellDev::Priors& priors,
                       const ApproxBayesRCellDev* sharedDesigns = nullptr);
    void addModality(const CellDev::Input& input, const CellDev::Priors& priors, int rank,
                     const ApproxBayesRCellDev* sharedDesigns = nullptr);
    void sampleUnknowns(unsigned iter) override;
    // Core starts are set in the constructor; inherited BayesC start logic
    // operates on a different SNP member and must not be used here.
    void sampleStartVal() override;

private:
    struct ModalityBlock {
        std::unique_ptr<CellDev::Sampler> sampler;
        std::unique_ptr<CellDev::Enrichment> enrichment;
        std::unique_ptr<CellDev::Reports> reports;
    };
    std::unique_ptr<CellDev::Reports> reports;
    std::vector<ModalityBlock> otherModalities;
    Eigen::MatrixXd jointPredictor() const;
    void sampleAnnotations();
    bool profileTiming = false;
    void collectAnnotationParameters();
    ApproxBayesRCellDev(const ApproxBayesRCellDev&) = delete;
    ApproxBayesRCellDev& operator=(const ApproxBayesRCellDev&) = delete;
};
// Generic aggregation of independently sampled CellDev chains.
class MultiChainCellDev : public Model {
public:
    std::vector<std::unique_ptr<ApproxBayesRCellDev> > chains;
    explicit MultiChainCellDev(std::vector<std::unique_ptr<ApproxBayesRCellDev> > chains);
    void sampleUnknowns(unsigned iter) override;
    void sampleStartVal() override;
private:
    struct Scalar {
        std::unique_ptr<Parameter> pooled;
        std::vector<Parameter*> source;
    };
    struct Set {
        std::unique_ptr<ParamSet> pooled;
        std::vector<ParamSet*> source;
    };
    std::vector<Scalar> scalars;
    std::vector<Set> sets;
    unsigned outerThreads, innerThreads;
    void collect();
};
#endif
