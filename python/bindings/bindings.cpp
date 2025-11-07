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
        // Data reading methods
        .def("read_fam_file", &Data::readFamFile, 
             py::arg("fam_file"),
             "Read PLINK .fam file")
        .def("read_bim_file", &Data::readBimFile, 
             py::arg("bim_file"),
             "Read PLINK .bim file")
        .def("read_bed_file", &Data::readBedFile, 
             py::arg("noscale"), py::arg("bed_file"),
             "Read PLINK .bed file")
        .def("read_phenotype_file", &Data::readPhenotypeFile, 
             py::arg("pheno_file"), py::arg("mphen"),
             "Read phenotype file")
        .def("read_covariate_file", &Data::readCovariateFile, 
             py::arg("covar_file"),
             "Read covariate file")
        .def("read_gwas_summary_file", &Data::readGwasSummaryFile,
             py::arg("gwas_file"), py::arg("af_diff"), py::arg("maf_min"),
             py::arg("maf_max"), py::arg("pvalue_threshold"), 
             py::arg("impute_n"), py::arg("remove_outlier_n"),
             "Read GWAS summary statistics file")
        .def("read_ld_matrix_info_file", &Data::readLDmatrixInfoFile,
             py::arg("ldm_file"),
             "Read LD matrix info file")
        .def("read_ld_matrix_bin_file", &Data::readLDmatrixBinFile,
             py::arg("ldm_file"),
             "Read LD matrix binary file")
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
        }, py::return_value_policy::reference_internal, "Get included SNP info vector");

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

