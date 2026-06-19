#include "quantizer.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

#include <zlib.h>

namespace eigen_quantize {

const char* kInputSuffix = ".eigen.bin";

static string join_path(const string& dir, const string& file) {
    if (dir.empty()) return file;
    if (dir[dir.size() - 1] == '/') return dir + file;
    return dir + "/" + file;
}

static bool is_directory(const string& path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

static bool is_regular_file(const string& path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0 && S_ISREG(st.st_mode);
}

static uint64_t file_size(const string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) {
        throw runtime_error("Cannot stat file " + path);
    }
    return static_cast<uint64_t>(st.st_size);
}

static void create_directories(const string& path) {
    if (path.empty() || is_directory(path)) return;

    string current;
    size_t pos = 0;
    if (path[0] == '/') {
        current = "/";
        pos = 1;
    }
    while (pos <= path.size()) {
        size_t next = path.find('/', pos);
        string part = path.substr(pos, next == string::npos ? string::npos : next - pos);
        if (!part.empty()) {
            if (!current.empty() && current[current.size() - 1] != '/') current += "/";
            current += part;
            if (!is_directory(current)) {
                if (mkdir(current.c_str(), 0755) != 0 && errno != EEXIST) {
                    throw runtime_error("Cannot create folder " + current + ": " + strerror(errno));
                }
            }
        }
        if (next == string::npos) break;
        pos = next + 1;
    }
}

static void copy_file_if_exists(const string& input_dir, const string& output_dir, const string& filename) {
    const string input_path = join_path(input_dir, filename);
    if (!is_regular_file(input_path)) return;
    const string output_path = join_path(output_dir, filename);

    ifstream in(input_path.c_str(), ios::binary);
    if (!in) throw runtime_error("Cannot open input file " + input_path);
    ofstream out(output_path.c_str(), ios::binary);
    if (!out) throw runtime_error("Cannot open output file " + output_path);
    out << in.rdbuf();
    if (!out) throw runtime_error("Write error copying " + output_path);
}

const char* output_suffix(const QuantizationOptions& options) {
    if (options.q_per_snp_column) {
        throw runtime_error("Q-column quantization is not enabled in this GCTB integration.");
    }
    if (options.bits == 8 && options.entropy_coding) return ".eigen.q8e.bin";
    if (options.bits == 8) return ".eigen.q8.bin";
    return ".eigen.q16.bin";
}

static float quantization_bound(int bits) {
    if (bits == 8) return 127.0f;
    return 32767.0f;
}

static int16_t quantize_value(float value, float scale, int bits) {
    if (scale <= 0.0f) return 0;

    const float bound = quantization_bound(bits);
    const float scaled = value / scale * bound;
    const long rounded = lround(scaled);
    const long limit = static_cast<long>(bound);
    const long clamped = std::max<long>(-limit, std::min<long>(limit, rounded));
    return static_cast<int16_t>(clamped);
}

static void write_zlib_matrix(FILE* out, const vector<int8_t>& raw, const string& output_path) {
    const uLong src_len = static_cast<uLong>(raw.size());
    uLong dst_cap = compressBound(src_len);
    vector<Bytef> compressed(dst_cap);
    uLong dst_len = dst_cap;
    const int zrc = compress(compressed.data(), &dst_len,
                             reinterpret_cast<const Bytef*>(raw.data()), src_len);
    if (zrc != Z_OK) {
        throw runtime_error("zlib compress failed for " + output_path);
    }

    const uint64_t uncompressed_size = static_cast<uint64_t>(raw.size());
    const uint64_t compressed_size = static_cast<uint64_t>(dst_len);
    if (fwrite(&uncompressed_size, sizeof(uint64_t), 1, out) != 1) {
        throw runtime_error("Write error (uncompressed size) in " + output_path);
    }
    if (fwrite(&compressed_size, sizeof(uint64_t), 1, out) != 1) {
        throw runtime_error("Write error (compressed size) in " + output_path);
    }
    if (fwrite(compressed.data(), 1, dst_len, out) != dst_len) {
        throw runtime_error("Write error (compressed matrix) in " + output_path);
    }
}

static void quantize_file(const string& input_dir,
                          const string& output_dir,
                          const string& name,
                          const QuantizationOptions& options,
                          QuantizationSummary& summary) {
    const string input_path = join_path(input_dir, name + kInputSuffix);
    const string output_path = join_path(output_dir, name + output_suffix(options));

    FILE* fp = fopen(input_path.c_str(), "rb");
    if (!fp) {
        throw runtime_error("Cannot open input file " + input_path);
    }

    FILE* out = fopen(output_path.c_str(), "wb");
    if (!out) {
        fclose(fp);
        throw runtime_error("Cannot open output file " + output_path);
    }

    try {
        int32_t num_snps = 0;
        int32_t num_eigenvalues = 0;
        float sum_pos_eigval = 0.0f;
        float eigen_cutoff = 0.0f;

        if (fread(&num_snps, sizeof(int32_t), 1, fp) != 1) {
            throw runtime_error("Read error (m) in " + input_path);
        }
        if (fread(&num_eigenvalues, sizeof(int32_t), 1, fp) != 1) {
            throw runtime_error("Read error (k) in " + input_path);
        }
        if (fread(&sum_pos_eigval, sizeof(float), 1, fp) != 1) {
            throw runtime_error("Read error (sumPosEigVal) in " + input_path);
        }
        if (fread(&eigen_cutoff, sizeof(float), 1, fp) != 1) {
            throw runtime_error("Read error (eigenCutoff) in " + input_path);
        }
        if (num_snps < 0 || num_eigenvalues < 0) {
            throw runtime_error("Negative matrix dimensions in " + input_path);
        }

        vector<float> eigenvalues(num_eigenvalues);
        if (fread(eigenvalues.data(), sizeof(float), num_eigenvalues, fp) !=
            static_cast<size_t>(num_eigenvalues)) {
            throw runtime_error("Read error (eigenvalues) in " + input_path);
        }

        const long matrix_offset = ftell(fp);
        if (matrix_offset < 0) {
            throw runtime_error("Failed to record matrix offset in " + input_path);
        }

        vector<float> column(num_snps);
        vector<float> scales(num_eigenvalues, 0.0f);

        for (int32_t col = 0; col < num_eigenvalues; ++col) {
            if (fread(column.data(), sizeof(float), num_snps, fp) !=
                static_cast<size_t>(num_snps)) {
                throw runtime_error("Read error (eigenvector column) in " + input_path);
            }

            float scale = 0.0f;
            for (int32_t row = 0; row < num_snps; ++row) {
                scale = std::max(scale, fabs(column[row]));
            }
            scales[col] = scale;
        }

        if (fseek(fp, matrix_offset, SEEK_SET) != 0) {
            throw runtime_error("Failed to rewind eigenvectors in " + input_path);
        }

        if (fwrite(&num_snps, sizeof(int32_t), 1, out) != 1) {
            throw runtime_error("Write error (m) in " + output_path);
        }
        if (fwrite(&num_eigenvalues, sizeof(int32_t), 1, out) != 1) {
            throw runtime_error("Write error (k) in " + output_path);
        }
        if (fwrite(&sum_pos_eigval, sizeof(float), 1, out) != 1) {
            throw runtime_error("Write error (sumPosEigVal) in " + output_path);
        }
        if (fwrite(&eigen_cutoff, sizeof(float), 1, out) != 1) {
            throw runtime_error("Write error (eigenCutoff) in " + output_path);
        }
        if (fwrite(eigenvalues.data(), sizeof(float), num_eigenvalues, out) !=
            static_cast<size_t>(num_eigenvalues)) {
            throw runtime_error("Write error (eigenvalues) in " + output_path);
        }
        if (fwrite(scales.data(), sizeof(float), num_eigenvalues, out) !=
            static_cast<size_t>(num_eigenvalues)) {
            throw runtime_error("Write error (column scales) in " + output_path);
        }

        const bool q8_entropy = (options.bits == 8 && options.entropy_coding);
        vector<int8_t> matrix_q8;
        if (q8_entropy) {
            const uint64_t total = static_cast<uint64_t>(num_snps) *
                                   static_cast<uint64_t>(num_eigenvalues);
            matrix_q8.resize(total);
        }

        for (int32_t col = 0; col < num_eigenvalues; ++col) {
            if (fread(column.data(), sizeof(float), num_snps, fp) !=
                static_cast<size_t>(num_snps)) {
                throw runtime_error("Read error (eigenvector column) in " + input_path);
            }

            const float scale = scales[col];
            if (options.bits == 8) {
                if (q8_entropy) {
                    int8_t* col_dst = matrix_q8.data() +
                                      static_cast<size_t>(col) * static_cast<size_t>(num_snps);
                    for (int32_t row = 0; row < num_snps; ++row) {
                        col_dst[row] =
                            static_cast<int8_t>(quantize_value(column[row], scale, options.bits));
                    }
                } else {
                    vector<int8_t> quantized_column(num_snps);
                    for (int32_t row = 0; row < num_snps; ++row) {
                        quantized_column[row] =
                            static_cast<int8_t>(quantize_value(column[row], scale, options.bits));
                    }
                    if (fwrite(quantized_column.data(), sizeof(int8_t), num_snps, out) !=
                        static_cast<size_t>(num_snps)) {
                        throw runtime_error("Write error (quantized eigenvector column) in " +
                                            output_path);
                    }
                }
            } else {
                vector<int16_t> quantized_column(num_snps);
                for (int32_t row = 0; row < num_snps; ++row) {
                    quantized_column[row] = quantize_value(column[row], scale, options.bits);
                }
                if (fwrite(quantized_column.data(), sizeof(int16_t), num_snps, out) !=
                    static_cast<size_t>(num_snps)) {
                    throw runtime_error("Write error (quantized eigenvector column) in " +
                                        output_path);
                }
            }
        }

        if (q8_entropy) {
            write_zlib_matrix(out, matrix_q8, output_path);
        }

        fclose(fp);
        fclose(out);

        const uint64_t original_bytes = file_size(input_path);
        const uint64_t quantized_bytes = file_size(output_path);

        cout << "Quantized " << name << ".eigen.bin"
             << " -> " << name << output_suffix(options)
             << " (m=" << num_snps
             << ", k=" << num_eigenvalues
             << ", q=" << options.bits
             << (options.entropy_coding ? ", zlib" : "")
             << ", bytes " << original_bytes
             << " -> " << quantized_bytes << ")\n" << flush;

        summary.total_original_bytes += original_bytes;
        summary.total_quantized_bytes += quantized_bytes;
        summary.num_files += 1;
    } catch (...) {
        fclose(fp);
        fclose(out);
        throw;
    }
}

static vector<string> list_input_names(const string& input_dir) {
    vector<string> names;
    DIR* dir = opendir(input_dir.c_str());
    if (!dir) {
        throw runtime_error("Cannot open input folder " + input_dir);
    }

    const size_t suffix_len = string(kInputSuffix).size();
    struct dirent* entry = NULL;
    while ((entry = readdir(dir)) != NULL) {
        const string filename(entry->d_name);
        if (filename.size() < suffix_len) continue;
        if (filename.substr(filename.size() - suffix_len) != kInputSuffix) continue;

        const string full_path = join_path(input_dir, filename);
        if (!is_regular_file(full_path)) continue;
        names.push_back(filename.substr(0, filename.size() - suffix_len));
    }
    closedir(dir);

    sort(names.begin(), names.end());
    return names;
}

QuantizationSummary quantize_directory(const string& input_dir,
                                       const string& output_dir,
                                       const QuantizationOptions& options) {
    if (!is_directory(input_dir)) {
        throw runtime_error("Input folder does not exist: " + input_dir);
    }
    if (options.entropy_coding && options.bits != 8) {
        throw runtime_error("--entropy is only supported with 8-bit quantization.");
    }
    if (options.bits != 8 && options.bits != 16) {
        throw runtime_error("Only q8 and q16 are supported.");
    }
    if (options.q_per_snp_column) {
        throw runtime_error("Q-column quantization is not enabled in this GCTB integration.");
    }

    create_directories(output_dir);
    copy_file_if_exists(input_dir, output_dir, "ldm.info");
    copy_file_if_exists(input_dir, output_dir, "snp.info");

    const vector<string> names = list_input_names(input_dir);
    if (names.empty()) {
        throw runtime_error("No .eigen.bin files found in " + input_dir);
    }

    QuantizationSummary summary;
    for (size_t i = 0; i < names.size(); ++i) {
        quantize_file(input_dir, output_dir, names[i], options, summary);
    }

    return summary;
}

}  // namespace eigen_quantize
