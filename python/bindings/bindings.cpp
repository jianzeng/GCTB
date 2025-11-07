/*
 * GCTB Python Bindings
 * 
 * This file contains pybind11 bindings for the GCTB C++ core library
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>

// Include GCTB headers
#include "data.hpp"
#include "gctb.hpp"
#include "model.hpp"
#include "mcmc.hpp"
#include "options.hpp"
#include "stat.hpp"
#include "gadgets.hpp"
#include "hsq.hpp"
#include "predict.hpp"
#include "stratify.hpp"
#include "xci.hpp"
#include "vgmaf.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = "GCTB: Genome-wide Complex Trait Bayesian Analysis - C++ Core";

    // ============================================================================
    // SnpInfo Class
    // ============================================================================
    py::class_<SnpInfo>(m, "SnpInfo", "SNP information structure")
        .def(py::init<const int, const std::string&, const std::string&, 
                      const std::string&, const int, const float, const int>(),
             py::arg("idx"), py::arg("id"), py::arg("allele1"), 
             py::arg("allele2"), py::arg("chr"), py::arg("gpos"), py::arg("ppos"))
        .def_readonly("ID", &SnpInfo::ID, "SNP identifier")
        .def_readonly("a1", &SnpInfo::a1, "Reference allele")
        .def_readonly("a2", &SnpInfo::a2, "Coded allele")
        .def_readonly("chrom", &SnpInfo::chrom, "Chromosome")
        .def_readonly("physPos", &SnpInfo::physPos, "Physical position")
        .def_readwrite("genPos", &SnpInfo::genPos, "Genetic position")
        .def_readwrite("index", &SnpInfo::index, "Index")
        .def_readwrite("af", &SnpInfo::af, "Allele frequency")
        .def_readwrite("twopq", &SnpInfo::twopq, "2pq")
        .def_readwrite("included", &SnpInfo::included, "Inclusion flag")
        .def_readwrite("effect", &SnpInfo::effect, "Estimated effect")
        .def_readwrite("pip", &SnpInfo::pip, "Posterior inclusion probability")
        .def_readwrite("varExplained", &SnpInfo::varExplained, "Variance explained")
        .def_readwrite("gwas_b", &SnpInfo::gwas_b, "GWAS beta")
        .def_readwrite("gwas_se", &SnpInfo::gwas_se, "GWAS standard error")
        .def_readwrite("gwas_n", &SnpInfo::gwas_n, "GWAS sample size")
        .def_readwrite("gwas_pvalue", &SnpInfo::gwas_pvalue, "GWAS p-value");

    // ============================================================================
    // IndInfo Class
    // ============================================================================
    py::class_<IndInfo>(m, "IndInfo", "Individual information structure")
        .def(py::init<const int, const std::string&, const std::string&, 
                      const std::string&, const std::string&, const int>(),
             py::arg("idx"), py::arg("fid"), py::arg("pid"), 
             py::arg("dad"), py::arg("mom"), py::arg("sex"))
        .def_readonly("famID", &IndInfo::famID, "Family ID")
        .def_readonly("indID", &IndInfo::indID, "Individual ID")
        .def_readonly("catID", &IndInfo::catID, "Concatenated ID")
        .def_readonly("fatherID", &IndInfo::fatherID, "Father ID")
        .def_readonly("motherID", &IndInfo::motherID, "Mother ID")
        .def_readonly("sex", &IndInfo::sex, "Sex (1=male, 2=female)")
        .def_readwrite("index", &IndInfo::index, "Index")
        .def_readwrite("kept", &IndInfo::kept, "Keep flag")
        .def_readwrite("phenotype", &IndInfo::phenotype, "Phenotype value");

    // ============================================================================
    // Data Class
    // ============================================================================
    py::class_<Data>(m, "Data", "Main data container for GCTB")
        .def(py::init<>())
        // Data reading methods with exception handling
        .def("read_fam_file", [](Data& data, const std::string& fam_file) {
            try {
                data.readFamFile(fam_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("Error reading FAM file: ") + e.what());
            }
        }, py::arg("fam_file"), "Read PLINK .fam file")
        
        .def("read_bim_file", [](Data& data, const std::string& bim_file) {
            try {
                data.readBimFile(bim_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("Error reading BIM file: ") + e.what());
            }
        }, py::arg("bim_file"), "Read PLINK .bim file")
        
        .def("read_bed_file", [](Data& data, bool noscale, const std::string& bed_file) {
            try {
                data.readBedFile(noscale, bed_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("Error reading BED file: ") + e.what());
            }
        }, py::arg("noscale"), py::arg("bed_file"), "Read PLINK .bed file")
        .def("read_phenotype_file", [](Data& data, const std::string& pheno_file, unsigned mphen) {
            try {
                data.readPhenotypeFile(pheno_file, mphen);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("Error reading phenotype file: ") + e.what());
            }
        }, py::arg("pheno_file"), py::arg("mphen"), "Read phenotype file")
        
        .def("read_covariate_file", [](Data& data, const std::string& covar_file) {
            try {
                data.readCovariateFile(covar_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("Error reading covariate file: ") + e.what());
            }
        }, py::arg("covar_file"), "Read covariate file")
        
        .def("read_gwas_summary_file", [](Data& data, const std::string& gwas_file, 
                                          float af_diff, float maf_min, float maf_max,
                                          float pvalue_threshold, bool impute_n, bool remove_outlier_n) {
            try {
                data.readGwasSummaryFile(gwas_file, af_diff, maf_min, maf_max, 
                                        pvalue_threshold, impute_n, remove_outlier_n);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("Error reading GWAS summary file: ") + e.what());
            }
        }, py::arg("gwas_file"), py::arg("af_diff"), py::arg("maf_min"),
           py::arg("maf_max"), py::arg("pvalue_threshold"), 
           py::arg("impute_n"), py::arg("remove_outlier_n"),
           "Read GWAS summary statistics file")
        .def("read_ld_matrix_info_file", &Data::readLDmatrixInfoFile,
             py::arg("ldm_file"),
             "Read LD matrix info file")
        .def("read_ld_matrix_bin_file", &Data::readLDmatrixBinFile,
             py::arg("ldm_file"),
             "Read LD matrix binary file")
        // SNP/Individual selection methods
        .def("include_matched_snp", &Data::includeMatchedSnp,
             "Build the list of included SNPs (must call after reading BIM)")
        .def("build_kept_individuals", [](Data& data) {
            // Build kept individuals list (keeps all by default)
            data.keptIndInfoVec = data.makeKeptIndInfoVec(data.indInfoVec);
            data.numKeptInds = (unsigned)data.keptIndInfoVec.size();
        }, "Build the list of kept individuals (call after reading FAM)")
        
        // Data properties
        .def_readonly("num_snps", &Data::numSnps, "Number of SNPs")
        .def_readonly("num_inds", &Data::numInds, "Number of individuals")
        .def_readonly("num_incd_snps", &Data::numIncdSnps, "Number of included SNPs")
        .def_readonly("num_kept_inds", &Data::numKeptInds, "Number of kept individuals")
        .def_readonly("var_phenotypic", &Data::varPhenotypic, "Phenotypic variance")
        .def_readonly("var_genotypic", &Data::varGenotypic, "Genotypic variance")
        .def_readonly("var_residual", &Data::varResidual, "Residual variance")
        .def_readwrite("title", &Data::title, "Title/label for output")
        // Data vectors (read-only access)
        .def("get_snp_info_vec", [](const Data& data) {
            return data.snpInfoVec;
        }, py::return_value_policy::reference_internal, "Get SNP info vector")
        .def("get_ind_info_vec", [](const Data& data) {
            return data.indInfoVec;
        }, py::return_value_policy::reference_internal, "Get individual info vector")
        .def("get_incd_snp_info_vec", [](const Data& data) {
            return data.incdSnpInfoVec;
        }, py::return_value_policy::reference_internal, "Get included SNP info vector")
        
        // Data initialization methods
        .def("init_variances", &Data::initVariances,
             py::arg("heritability"), py::arg("prop_var_random"),
             "Initialize variance components (must call before building model)")
        
        // Matrix building methods (for summary stats)
        .def("build_sparse_mme", [](Data& data, bool sample_overlap, bool noscale) {
            try {
                data.buildSparseMME(sample_overlap, noscale);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("sample_overlap") = false, py::arg("noscale") = false,
           "Build sparse MME for summary statistics");

    // ============================================================================
    // Model Class (Opaque - don't expose internals)
    // ============================================================================
    py::class_<Model>(m, "Model", "Bayesian model (opaque - use via build_model)")
        .def_readonly("num_snps", &Model::numSnps, "Number of SNPs in model");

    // ============================================================================
    // GCTB Class & Factory Functions
    // ============================================================================
    py::class_<GCTB>(m, "GCTB", "Main GCTB controller")
        .def(py::init<Options&>());
    
    // Factory function for summary-stats models (ApproxBayes* - no genotypes needed)
    m.def("build_model_summary", [](Data& data, const std::string& bayes_type,
                                    float heritability, float pi) -> Model* {
        try {
            // Initialize variances
            data.initVariances(heritability, 0.05f);
            
            // Create dummy options
            Options dummy_opts;
            GCTB gctb(dummy_opts);
            
            // Default parameters
            VectorXf pis(4);
            pis << 0.95f, 0.02f, 0.02f, 0.01f;
            VectorXf gamma(4);  
            gamma << 0.0f, 0.01f, 0.1f, 1.0f;
            VectorXf piPar = VectorXf::Ones(4);
            std::vector<float> S = {0.0f};
            
            // Pass "dummy" for gwasFile to force summary-stats path
            return gctb.buildModel(data, "", "dummy", bayes_type, 0,
                                  heritability, 0.05f, pi, 1.0f, 1.0f, 
                                  true, false, pis, piPar, gamma, true,
                                  0.0f, 10.0f, "Gibbs", 2, 1.0f, S, 0.0f, false,
                                  0.0f, 0.0f, false, true, false, false, false);
        } catch (const std::string& e) {
            throw std::runtime_error(e);
        } catch (const char* e) {
            throw std::runtime_error(e);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Error building model: ") + e.what());
        }
    }, py::arg("data"), py::arg("bayes_type"),
       py::arg("heritability") = 0.5f, py::arg("pi") = 0.01f,
       py::return_value_policy::take_ownership,
       "Build a summary-stats Bayesian model (ApproxBayes* - for summary data)");

    // ============================================================================
    // McmcSamples Class
    // ============================================================================
    py::class_<McmcSamples>(m, "McmcSamples", "MCMC sample storage")
        .def_readonly("label", &McmcSamples::label, "Parameter label")
        .def_readonly("chain_length", &McmcSamples::chainLength, "Chain length")
        .def_readonly("burnin", &McmcSamples::burnin, "Burn-in iterations")
        .def_readonly("thin", &McmcSamples::thin, "Thinning interval")
        .def_readonly("nrow", &McmcSamples::nrow, "Number of rows (samples)")
        .def_readonly("ncol", &McmcSamples::ncol, "Number of columns (parameters)")
        .def_readonly("posterior_mean", &McmcSamples::posteriorMean, "Posterior mean")
        .def_readonly("posterior_sqr_mean", &McmcSamples::posteriorSqrMean, "Posterior squared mean")
        .def_readonly("pip", &McmcSamples::pip, "Posterior inclusion probability")
        .def_readonly("last_sample", &McmcSamples::lastSample, "Last MCMC sample")
        .def("to_dict", [](const McmcSamples& samples) {
            py::dict d;
            d["label"] = samples.label;
            d["posterior_mean"] = samples.posteriorMean;
            d["pip"] = samples.pip;
            d["last_sample"] = samples.lastSample;
            return d;
        }, "Convert to Python dictionary");

    // ============================================================================
    // MCMC Class
    // ============================================================================
    py::class_<MCMC>(m, "MCMC", "MCMC sampler")
        .def(py::init<>())
        .def("run", [](MCMC& mcmc, Model& model, unsigned chain_length,
                      unsigned burnin, unsigned thin, unsigned output_freq,
                      const std::string& title) -> std::vector<McmcSamples*> {
            try {
                py::gil_scoped_release release;  // Release GIL during computation
                return mcmc.run(model, chain_length, burnin, thin, true,
                               output_freq, title, false, false);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("Error in MCMC: ") + e.what());
            }
        }, py::arg("model"), py::arg("chain_length") = 3000,
           py::arg("burnin") = 1000, py::arg("thin") = 10,
           py::arg("output_freq") = 100, py::arg("title") = "gctb",
           py::return_value_policy::take_ownership,
           "Run MCMC sampler");
    
    // Convenience function combining GCTB.runMcmc
    m.def("run_mcmc", [](Model& model, unsigned chain_length, unsigned burnin,
                        unsigned thin, unsigned output_freq, const std::string& title) {
        try {
            MCMC mcmc;
            py::gil_scoped_release release;  // Release GIL  
            return mcmc.run(model, chain_length, burnin, thin, true,
                           output_freq, title, false, false);
        } catch (const std::string& e) {
            throw std::runtime_error(e);
        } catch (const char* e) {
            throw std::runtime_error(e);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Error running MCMC: ") + e.what());
        }
    }, py::arg("model"), py::arg("chain_length") = 3000,
       py::arg("burnin") = 1000, py::arg("thin") = 10,
       py::arg("output_freq") = 100, py::arg("title") = "gctb",
       py::return_value_policy::take_ownership,
       "Run MCMC sampler (convenience function)");

    // ============================================================================
    // Gadget::Timer Class
    // ============================================================================
    py::class_<Gadget::Timer>(m, "Timer", "Timer utility")
        .def(py::init<>())
        .def("set_time", &Gadget::Timer::setTime, "Set start time")
        .def("get_time", &Gadget::Timer::getTime, "Get current time")
        .def("get_date", &Gadget::Timer::getDate, "Get current date")
        .def("get_elapse", &Gadget::Timer::getElapse, "Get elapsed time")
        .def("format", &Gadget::Timer::format, 
             py::arg("seconds"), "Format seconds to readable string");

    // ============================================================================
    // Version Info
    // ============================================================================
    m.attr("__version__") = "3.0.0";
    m.attr("__author__") = "Jian Zeng, Luke Lloyd-Jones, Zhili Zheng, Shouye Liu";
}

