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
        .def("read_random_covariate_file", [](Data& data, const std::string& covar_file) {
            try {
                data.readRandomCovariateFile(covar_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("Error reading random covariate file: ") + e.what());
            }
        }, py::arg("covar_file"),
           "Read random covariate file")
        .def("read_residual_diag_file", [](Data& data, const std::string& resdiag_file) {
            try {
                data.readResidualDiagFile(resdiag_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("Error reading residual diag file: ") + e.what());
            }
        }, py::arg("resdiag_file"),
           "Read residual diagonal file")
        
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
        .def("read_genetic_map_file", [](Data& data, const std::string& map_file) {
            try {
                data.readGeneticMapFile(map_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("map_file"),
           "Read genetic map file")
        .def("read_ld_block_info_file", [](Data& data, const std::string& ldb_file) {
            try {
                data.readLDBlockInfoFile(ldb_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("ldb_file"),
           "Read LD block info file")
        .def("include_snp", [](Data& data, const std::string& include_snp_file) {
            try {
                data.includeSnp(include_snp_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("include_snp_file"),
           "Include SNPs based on file")
        .def("exclude_snp", [](Data& data, const std::string& exclude_snp_file) {
            try {
                data.excludeSnp(exclude_snp_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("exclude_snp_file"),
           "Exclude SNPs based on file")
        .def("include_chr", [](Data& data, unsigned chr) {
            try {
                data.includeChr(chr);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("chromosome"),
           "Restrict analysis to a chromosome")
        .def("include_block", [](Data& data, unsigned block) {
            try {
                data.includeBlock(block);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("block"),
           "Restrict analysis to an LD block")
        .def("exclude_region", [](Data& data, const std::string& region_file) {
            try {
                data.excludeRegion(region_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("region_file"),
           "Exclude SNPs in specified regions")
        .def("exclude_mhc", [](Data& data) {
            try {
                data.excludeMHC();
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, "Exclude SNPs in the MHC region")
        .def("exclude_ambiguous_snp", [](Data& data) {
            try {
                data.excludeAmbiguousSNP();
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, "Exclude ambiguous SNPs (A/T, C/G)")
        .def("include_skeleton_snp", [](Data& data, const std::string& skeleton_file) {
            try {
                data.includeSkeletonSnp(skeleton_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("skeleton_file"),
           "Include only skeleton SNPs")
        .def("read_annotation_file", [](Data& data, const std::string& annotation_file,
                                        bool transpose, bool allow_multi_anno) {
            try {
                data.readAnnotationFile(annotation_file, transpose, allow_multi_anno);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("annotation_file"), py::arg("transpose") = false,
           py::arg("allow_multi_anno") = true,
           "Read categorical annotation file")
        .def("read_annotation_file_format2", [](Data& data, const std::string& continuous_file,
                                                unsigned flank, const std::string& eqtl_file) {
            try {
                data.readAnnotationFileFormat2(continuous_file, flank, eqtl_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("continuous_annotation_file"), py::arg("flank"),
           py::arg("eqtl_file") = "",
           "Read continuous annotation file (format 2)")
        .def("set_annotation_info", [](Data& data) {
            try {
                data.setAnnoInfoVec();
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, "Finalize annotation info structures")
        .def("read_ldscore_file", [](Data& data, const std::string& ldsc_file) {
            try {
                data.readLDscoreFile(ldsc_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("ldsc_file"),
           "Read LD score file")
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
        
        .def("keep_matched_ind", [](Data& data, const std::string& keep_ind_file, 
                                   unsigned keep_ind_max) {
            try {
                data.keepMatchedInd(keep_ind_file, keep_ind_max);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("keep_ind_file") = "", py::arg("keep_ind_max") = UINT_MAX,
           "Process individuals - initializes matrices and filters based on phenotypes/covariates")
        
        // Data properties
        .def_readonly("num_snps", &Data::numSnps, "Number of SNPs")
        .def_readonly("num_inds", &Data::numInds, "Number of individuals")
        .def_readonly("num_incd_snps", &Data::numIncdSnps, "Number of included SNPs")
        .def_readonly("num_kept_inds", &Data::numKeptInds, "Number of kept individuals")
        .def_readonly("num_annos", &Data::numAnnos, "Number of annotations")
        .def_readonly("num_ld_blocks", &Data::numLDBlocks, "Number of LD blocks")
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
           "Build sparse MME for summary statistics")
        .def("read_multi_ld_matrix_info_file", [](Data& data, const std::string& mldm_file) {
            try {
                data.readMultiLDmatInfoFile(mldm_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("mldm_file"),
           "Read multi-chromosome LD matrix info file")
        .def("read_multi_ld_matrix_bin_file", [](Data& data, const std::string& mldm_file) {
            try {
                data.readMultiLDmatBinFile(mldm_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("mldm_file"),
           "Read multi-chromosome LD matrix binary file")
        .def("read_ld_matrix_txt_file", [](Data& data, const std::string& ldm_file) {
            try {
                data.readLDmatrixTxtFile(ldm_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("ldm_file"),
           "Read LD matrix text file")
        .def("read_plink_af_file", [](Data& data, const std::string& af_file) {
            try {
                data.readPlinkAFfile(af_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("af_file"),
           "Read PLINK allele frequency file")
        .def("read_plink_ld_txt_file", [](Data& data, const std::string& ld_txt_file) {
            try {
                data.readPlinkLDtxtfile(ld_txt_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("ld_txt_file"),
           "Read PLINK LD text file")
        .def("read_plink_ld_bin_file", [](Data& data, const std::string& ld_bin_file) {
            try {
                data.readPlinkLDbinfile(ld_bin_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("ld_bin_file"),
           "Read PLINK LD binary file")
        .def("part_ld_matrix", [](Data& data, const std::string& part_param,
                                  const std::string& out_filename,
                                  const std::string& ldmat_type) {
            try {
                return data.partLDMatrix(part_param, out_filename, ldmat_type);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("part_param"), py::arg("out_filename"), py::arg("ldmat_type"),
           "Partition LD matrix according to parameter file")
        .def("make_ld_matrix", [](Data& data, const std::string& bed_file,
                                  const std::string& ldmat_type,
                                  float chisq_threshold, float ld_threshold,
                                  unsigned window_width, const std::string& snp_range,
                                  const std::string& filename, bool write_ldm_txt) {
            try {
                data.makeLDmatrix(bed_file, ldmat_type, chisq_threshold, ld_threshold,
                                  window_width, snp_range, filename, write_ldm_txt);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("bed_file"), py::arg("ldmat_type"), py::arg("chisq_threshold"),
           py::arg("ld_threshold"), py::arg("window_width"), py::arg("snp_range"),
           py::arg("filename"), py::arg("write_ldm_txt") = false,
           "Construct LD matrix from genotype data")
        .def("make_shrunk_ld_matrix", [](Data& data, const std::string& bed_file,
                                         const std::string& ldmat_type,
                                         const std::string& snp_range,
                                         const std::string& filename, bool write_ldm_txt,
                                         float effpop_ne, float cutoff, float gen_map_n) {
            try {
                data.makeshrunkLDmatrix(bed_file, ldmat_type, snp_range, filename,
                                        write_ldm_txt, effpop_ne, cutoff, gen_map_n);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("bed_file"), py::arg("ldmat_type"), py::arg("snp_range"),
           py::arg("filename"), py::arg("write_ldm_txt") = false,
           py::arg("effpop_ne") = 11400.0f, py::arg("cutoff") = 1.0f,
           py::arg("gen_map_n") = 0.0f,
           "Construct shrunk LD matrix")
        .def("make_block_ld_matrix", [](Data& data, const std::string& bed_file,
                                        const std::string& ldmat_type, unsigned block,
                                        const std::string& filename, bool write_ldm_txt,
                                        int ld_block_region_wind) {
            try {
                data.makeBlockLDmatrix(bed_file, ldmat_type, block, filename,
                                       write_ldm_txt, ld_block_region_wind);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("bed_file"), py::arg("ldmat_type"), py::arg("block"),
           py::arg("filename"), py::arg("write_ldm_txt") = false,
           py::arg("ld_block_region_wind") = 0,
           "Construct block-specific LD matrix")
        .def("resize_ld_matrix", [](Data& data, const std::string& ldmat_type,
                                    float chisq_threshold, unsigned window_width,
                                    float ld_threshold, float effpop_ne, float cutoff,
                                    float gen_map_n) {
            try {
                data.resizeLDmatrix(ldmat_type, chisq_threshold, window_width,
                                    ld_threshold, effpop_ne, cutoff, gen_map_n);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("ldmat_type"), py::arg("chisq_threshold"), py::arg("window_width"),
           py::arg("ld_threshold"), py::arg("effpop_ne"), py::arg("cutoff"),
           py::arg("gen_map_n"),
           "Resize existing LD matrix according to thresholds")
        .def("output_ld_matrix", [](Data& data, const std::string& ldmat_type,
                                    const std::string& filename, bool write_ldm_txt) {
            try {
                data.outputLDmatrix(ldmat_type, filename, write_ldm_txt);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("ldmat_type"), py::arg("filename"),
           py::arg("write_ldm_txt") = false,
           "Write LD matrix to disk")
        .def("direct_prune_ld_matrix", [](Data& data, const std::string& ldm_file,
                                          const std::string& out_ldmat_type,
                                          float chisq_threshold, const std::string& title,
                                          bool write_ldm_txt) {
            try {
                data.directPruneLDmatrix(ldm_file, out_ldmat_type, chisq_threshold,
                                         title, write_ldm_txt);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("ldm_file"), py::arg("out_ldmat_type"),
           py::arg("chisq_threshold"), py::arg("title"),
           py::arg("write_ldm_txt") = false,
           "Directly prune LD matrix")
        .def("jackknife_ld_matrix", [](Data& data, const std::string& ldm_file,
                                       const std::string& out_ldmat_type,
                                       const std::string& title, bool write_ldm_txt) {
            try {
                data.jackknifeLDmatrix(ldm_file, out_ldmat_type, title, write_ldm_txt);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("ldm_file"), py::arg("out_ldmat_type"),
           py::arg("title"), py::arg("write_ldm_txt") = false,
           "Perform jackknife on LD matrix")
        .def("merge_ldm_info", [](Data& data, const std::string& out_ldmat_type,
                                  const std::string& dirname) {
            try {
                data.mergeLdmInfo(out_ldmat_type, dirname);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("out_ldmat_type"), py::arg("dirname"),
           "Merge LD matrix info files")
        .def("filter_snp_by_ld_rsq", [](Data& data, float rsq_threshold) {
            try {
                data.filterSnpByLDrsq(rsq_threshold);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("rsq_threshold"),
           "Filter SNPs by LD r-squared threshold")
        .def("bin_snp_by_ld_rsq", [](Data& data, float rsq_threshold,
                                     const std::string& title) {
            try {
                data.binSnpByLDrsq(rsq_threshold, title);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("rsq_threshold"), py::arg("title"),
           "Bin SNPs by LD r-squared threshold")
        .def("read_window_file", [](Data& data, const std::string& window_file) {
            try {
                data.readWindowFile(window_file);
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, py::arg("window_file"),
           "Read window definition file")
        .def("bin_snp_by_window_id", [](Data& data) {
            try {
                data.binSnpByWindowID();
            } catch (const std::string& e) {
                throw std::runtime_error(e);
            } catch (const char* e) {
                throw std::runtime_error(e);
            }
        }, "Assign SNPs to windows by ID");

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
    
    // Model factory for individual-level data (Bayes*)
    m.def("build_model", [](Data& data, const std::string& bayes_type,
                           float heritability, float pi) -> Model* {
        try {
            data.initVariances(heritability, 0.05f);
            
            if (bayes_type == "C") {
                return new BayesC(data, data.varGenotypic, data.varResidual, data.varRandom,
                                 pi, 1.0f, 1.0f, true, false, "Gibbs");
            } 
            else if (bayes_type == "B") {
                return new BayesB(data, data.varGenotypic, data.varResidual, data.varRandom,
                                 pi, 1.0f, 1.0f, true, false);
            }
            else if (bayes_type == "R") {
                VectorXf pis(4);
                pis << 0.95f, 0.02f, 0.02f, 0.01f;
                VectorXf gamma(4);
                gamma << 0.0f, 0.01f, 0.1f, 1.0f;
                VectorXf piPar = VectorXf::Ones(4);
                return new BayesR(data, data.varGenotypic, data.varResidual, data.varRandom,
                                 pis, piPar, gamma, true, false, true, "Gibbs");
            }
            else if (bayes_type == "S") {
                std::vector<float> S = {0.0f};
                return new BayesS(data, data.varGenotypic, data.varResidual, data.varRandom,
                                 pi, 1.0f, 1.0f, true, 1.0f, S, "Gibbs");
            }
            else {
                throw std::runtime_error("Unknown bayes_type: " + bayes_type + 
                                       ". Supported: C, B, R, S");
            }
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
       "Build Bayesian model for individual-level data (requires genotypes)");
    
    // Model factory for summary statistics (ApproxBayes*)
    m.def("build_model_summary", [](Data& data, const std::string& sbayes_type,
                                   float heritability, float pi) -> Model* {
        try {
            data.initVariances(heritability, 0.05f);
            
            // Default parameters
            VectorXf pis(4);
            pis << 0.95f, 0.02f, 0.02f, 0.01f;
            VectorXf gamma(4);
            gamma << 0.0f, 0.01f, 0.1f, 1.0f;
            VectorXf piPar = VectorXf::Ones(4);
            std::vector<float> S = {0.0f};
            
            // Create ApproxBayes* models for summary statistics
            if (sbayes_type == "C") {
                return new ApproxBayesC(data, data.lowRankModel, data.varGenotypic, 
                                       data.varResidual, data.varRandom, pi, 1.0f, 1.0f, 
                                       true, false, 0.0f, 0.0f, false, 0.0f, 0.0f, false, false);
            }
            else if (sbayes_type == "R") {
                return new ApproxBayesR(data, data.lowRankModel, data.varGenotypic,
                                       data.varResidual, pis, piPar, gamma, true, true,
                                       false, true, 0.0f, false, 0.0f, false, false, "Gibbs", false);
            }
            else if (sbayes_type == "S") {
                return new ApproxBayesS(data, data.lowRankModel, data.varGenotypic,
                                       data.varResidual, pi, 1.0f, 1.0f, true, 0.0f, 0.0f,
                                       false, 0.0f, 0.0f, 1.0f, S, "Gibbs", false, false);
            }
            else {
                throw std::runtime_error("Unknown sbayes_type: " + sbayes_type +
                                       ". Supported: C, R, S");
            }
        } catch (const std::string& e) {
            throw std::runtime_error(e);
        } catch (const char* e) {
            throw std::runtime_error(e);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Error building summary-stats model: ") + e.what());
        }
    }, py::arg("data"), py::arg("sbayes_type"),
       py::arg("heritability") = 0.5f, py::arg("pi") = 0.01f,
       py::return_value_policy::take_ownership,
       "Build ApproxBayes model for summary statistics (requires LD matrix)");

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
    m.attr("__version__") = "1.1.0";
    m.attr("__author__") = "Jian Zeng, Luke Lloyd-Jones, Zhili Zheng, Shouye Liu";
}

