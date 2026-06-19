#ifndef EIGEN_QUANTIZE_QUANTIZER_HPP
#define EIGEN_QUANTIZE_QUANTIZER_HPP

#include <cstdint>
#include <string>

namespace eigen_quantize {

using namespace std;

struct QuantizationSummary {
    uint64_t total_original_bytes = 0;
    uint64_t total_quantized_bytes = 0;
    uint64_t num_files = 0;
};

struct QuantizationOptions {
    int bits = 16;
    /** If true and bits==8, write .eigen.q8e.bin with zlib-compressed int8 matrix. */
    bool entropy_coding = false;
    /**
     * If true, quantize Q = diag(sqrt(lambda)) * U' (k x m) with one scale per SNP column
     * (max abs in that column). Writes .eigen.q8qc.bin/.eigen.q16qc.bin
     * (and .eigen.q8eqc.bin with --entropy).
     */
    bool q_per_snp_column = false;
};

QuantizationSummary quantize_directory(const string& input_dir,
                                       const string& output_dir,
                                       const QuantizationOptions& options);

}  // namespace eigen_quantize

#endif
