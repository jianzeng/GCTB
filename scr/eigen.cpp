//
//  eigen.cpp
//  gctb
//
//  Created by Shouye Liu and Jian Zeng on 05/01/2023.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "data.hpp"
#include "stat.hpp"
#include <Eigen/Core>
#include <sys/stat.h>
#include <dirent.h>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <stdint.h>
#include <unordered_map>
#include <unordered_set>
#include <cstring>
#include <iomanip>
#include <limits.h>
#include <zlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

bool fileExistsAndNonEmptyLocal(const string &filename) {
    ifstream in(filename.c_str(), ios::ate | ios::binary);
    return in.is_open() && static_cast<long long>(in.tellg()) > 0;
}

bool pairwiseLdHasRecordsLocal(const string &filename);

void updateFnv1a(uint64_t &hash, const string &value) {
    for (unsigned i = 0; i < value.size(); ++i) {
        hash ^= (unsigned char)value[i];
        hash *= 1099511628211ULL;
    }
}

void updateFnv1aInt(uint64_t &hash, const long long value) {
    std::ostringstream oss;
    oss << value;
    updateFnv1a(hash, oss.str());
}

bool blockInfoSignatureLocal(const string &dirname, const vector<string> &blockFileIds, uint64_t &signature) {
    uint64_t hash = 1469598103934665603ULL;
    for (vector<string>::const_iterator it = blockFileIds.begin(); it != blockFileIds.end(); ++it) {
        updateFnv1a(hash, *it);
        const string files[2] = {
            dirname + "/block" + *it + ".snp.info",
            dirname + "/block" + *it + ".ldm.info"
        };
        for (unsigned f = 0; f < 2; ++f) {
            struct stat sb;
            if (stat(files[f].c_str(), &sb) != 0 || !S_ISREG(sb.st_mode)) return false;
            updateFnv1aInt(hash, (long long)sb.st_size);
            updateFnv1aInt(hash, (long long)sb.st_mtime);
        }
    }
    signature = hash;
    return true;
}

bool combinedMergeOutputsReadyLocal(const string &dirname, const bool requirePairwiseLD) {
    const string snpInfoFile = dirname + "/snp.info";
    const string ldmInfoFile = dirname + "/ldm.info";
    if (!fileExistsAndNonEmptyLocal(snpInfoFile) || !fileExistsAndNonEmptyLocal(ldmInfoFile)) return false;

    ifstream snpIn(snpInfoFile.c_str());
    ifstream ldmIn(ldmInfoFile.c_str());
    string header, firstRecord;
    if (!getline(snpIn, header) || !getline(snpIn, firstRecord) || firstRecord.empty()) return false;
    if (!getline(ldmIn, header) || !getline(ldmIn, firstRecord) || firstRecord.empty()) return false;

    if (!requirePairwiseLD) return true;
    return pairwiseLdHasRecordsLocal(dirname + "/rsq0.5.pwld.bin")
        || pairwiseLdHasRecordsLocal(dirname + "/rsq0.5.pwld")
        || pairwiseLdHasRecordsLocal(dirname + "/rsq0.5.pwld.gz");
}

bool readMergeMarkerLocal(const string &filename, unsigned &blockCount, uint64_t &signature, bool &pairwiseLDReady) {
    ifstream in(filename.c_str());
    if (!in) return false;
    string key;
    blockCount = 0;
    signature = 0;
    pairwiseLDReady = false;
    while (in >> key) {
        if (key == "BlockInfoCount") in >> blockCount;
        else if (key == "BlockInfoSignature") in >> signature;
        else if (key == "PairwiseLDReady") {
            int ready = 0;
            in >> ready;
            pairwiseLDReady = ready != 0;
        } else {
            string rest;
            getline(in, rest);
        }
    }
    return blockCount > 0 && signature != 0;
}

void writeMergeMarkerLocal(const string &filename, const unsigned blockCount, const uint64_t signature, const bool pairwiseLDReady) {
    ofstream out(filename.c_str());
    if (!out) throw("Error: can not write merged LD metadata marker [" + filename + "].");
    out << "Version 1\n";
    out << "BlockInfoCount " << blockCount << "\n";
    out << "BlockInfoSignature " << signature << "\n";
    out << "PairwiseLDReady " << (pairwiseLDReady ? 1 : 0) << "\n";
}

const char PWLD_BIN_MAGIC[8] = {'G','C','T','B','P','W','L','D'};
const uint32_t PWLD_BIN_VERSION = 2;
const char CORE_CONTEXT_BIN_MAGIC[8] = {'G','C','T','B','R','C','F','1'};
const uint32_t CORE_CONTEXT_BIN_VERSION = 1;
const char CORE_CONTEXT_META_MAGIC[8] = {'G','C','T','B','C','T','X','1'};
const char SOURCE_SPARSE_LD_MAGIC_V1[8] = {'G','C','T','B','S','L','D','1'};
const uint32_t SOURCE_SPARSE_LD_VERSION_V1 = 1;
const char SOURCE_SPARSE_LD_MAGIC_V2[8] = {'G','C','T','B','S','L','D','2'};
const uint32_t SOURCE_SPARSE_LD_VERSION_V2 = 2;
const char SOURCE_SPARSE_LD_MAGIC_V3[8] = {'G','C','T','B','S','L','D','3'};
const uint32_t SOURCE_SPARSE_LD_VERSION_V3 = 3;
const float SOURCE_SPARSE_LD_Q16_SCALE = 32767.0f;
const uint32_t SOURCE_SPARSE_LD_MISSING_INDEX = 0xffffffffu;
const char TILED_EIGEN_BUNDLE_MAGIC[8] = {'G','C','T','B','T','E','G','1'};
const uint32_t TILED_EIGEN_BUNDLE_VERSION = 1;
const char MARGINAL_CONTEXT_BUNDLE_MAGIC[8] = {'G','C','T','B','M','C','X','1'};
const uint32_t MARGINAL_CONTEXT_BUNDLE_VERSION = 1;

void cacheRealizedQdiag(LDBlockInfo *block, const MatrixXf &Q) {
    block->qDiag = Q.colwise().squaredNorm().transpose();
    for (int j = 0; j < block->qDiag.size(); ++j) {
        if (!(block->qDiag[j] >= 0.0f) || !isfinite(block->qDiag[j]))
            throw("Error: invalid realized Q'Q diagonal in LD block " + block->ID + ".");
    }
}

void cacheRealizedQdiag(LDBlockInfo *block, const QuantizedEigenBlock &Q) {
    block->qDiag.resize(Q.m);
    for (int j = 0; j < Q.m; ++j) {
        block->qDiag[j] = Q.dotQColumns(j, j);
        if (!(block->qDiag[j] >= 0.0f) || !isfinite(block->qDiag[j]))
            throw("Error: invalid realized quantized Q'Q diagonal in LD block " + block->ID + ".");
    }
}

struct PwldBinRecord {
    uint32_t snp1;
    uint32_t snp2;
    float ldcor;
};

struct CoreContextBinRecord {
    uint32_t coreLocal;
    uint32_t flankGlobal;
    float ldcor;
};

struct SparseSourceLdRecord {
    uint32_t snp1;
    uint32_t snp2;
    float ldcor;
};

void writeBinaryString(FILE *out, const string &value);

void assignTileByLocalIndex(const vector<vector<unsigned> > &tileLocalIdx,
                            const int numSnps, vector<int> &tileByLocalIndex) {
    tileByLocalIndex.assign(numSnps, -1);
    for (unsigned t = 0; t < tileLocalIdx.size(); ++t) {
        for (unsigned j = 0; j < tileLocalIdx[t].size(); ++j)
            tileByLocalIndex[(int)tileLocalIdx[t][j]] = (int)t;
    }
}

vector<vector<unsigned> > makeBalancedSnpCountTiles(const int numSnps,
                                                    const unsigned maxSnps,
                                                    const unsigned minSnps) {
    if (numSnps <= 0) return vector<vector<unsigned> >();
    if (!maxSnps) throw("Error: --tile-snps must be >= 1.");
    if (minSnps > maxSnps) throw("Error: --tile-min-snps cannot be greater than --tile-snps.");

    unsigned numTiles = ((unsigned)numSnps + maxSnps - 1) / maxSnps;
    if (!numTiles) numTiles = 1;
    while (numTiles > 1) {
        const unsigned nextNumTiles = numTiles - 1;
        const unsigned nextMaxSize = ((unsigned)numSnps + nextNumTiles - 1) / nextNumTiles;
        if (nextMaxSize > maxSnps) break;
        const unsigned nextMinSize = (unsigned)numSnps / nextNumTiles;
        if (minSnps && nextMinSize < minSnps) break;
        numTiles = nextNumTiles;
    }

    vector<vector<unsigned> > tiles;
    tiles.reserve(numTiles);
    const unsigned baseSize = (unsigned)numSnps / numTiles;
    const unsigned extra = (unsigned)numSnps % numTiles;
    unsigned idx = 0;
    for (unsigned t = 0; t < numTiles; ++t) {
        const unsigned tileSize = baseSize + (t < extra ? 1 : 0);
        vector<unsigned> tile;
        tile.reserve(tileSize);
        for (unsigned j = 0; j < tileSize; ++j) tile.push_back(idx++);
        tiles.push_back(tile);
    }
    return tiles;
}

vector<uint32_t> makeSparseSourceLocalToGlobal(const vector<SnpInfo*> &sourceSnps,
                                               const unordered_map<int, int> &newIndexByOldIndex) {
    vector<uint32_t> localToGlobal(sourceSnps.size(), SOURCE_SPARSE_LD_MISSING_INDEX);
    for (unsigned j = 0; j < sourceSnps.size(); ++j) {
        SnpInfo *snp = sourceSnps[j];
        if (!snp) continue;
        unordered_map<int, int>::const_iterator idxIt = newIndexByOldIndex.find(snp->index);
        if (idxIt == newIndexByOldIndex.end() || idxIt->second < 0) continue;
        localToGlobal[j] = (uint32_t)idxIt->second;
    }
    return localToGlobal;
}

long writeSparseSourceLdHeader(FILE *out, const string &sourceBlockID,
                               const vector<uint32_t> &localToGlobal,
                               const bool useUint16Index,
                               const bool q16Correlation) {
    if (q16Correlation) {
        fwrite(SOURCE_SPARSE_LD_MAGIC_V2, sizeof(char), 8, out);
        fwrite(&SOURCE_SPARSE_LD_VERSION_V2, sizeof(uint32_t), 1, out);
    } else {
        fwrite(SOURCE_SPARSE_LD_MAGIC_V3, sizeof(char), 8, out);
        fwrite(&SOURCE_SPARSE_LD_VERSION_V3, sizeof(uint32_t), 1, out);
    }
    writeBinaryString(out, sourceBlockID);
    uint32_t numSnps = (uint32_t)localToGlobal.size();
    uint64_t numRecords = 0;
    uint8_t indexWidth = useUint16Index ? 2 : 4;
    uint8_t corrWidth = q16Correlation ? 2 : 4;
    float corrScale = q16Correlation ? SOURCE_SPARSE_LD_Q16_SCALE : 0.0f;
    fwrite(&numSnps, sizeof(uint32_t), 1, out);
    const long countPos = ftell(out);
    fwrite(&numRecords, sizeof(uint64_t), 1, out);
    fwrite(&indexWidth, sizeof(uint8_t), 1, out);
    if (!q16Correlation)
        fwrite(&corrWidth, sizeof(uint8_t), 1, out);
    fwrite(&corrScale, sizeof(float), 1, out);
    if (numSnps && fwrite(localToGlobal.data(), sizeof(uint32_t), numSnps, out) != numSnps) {
        fclose(out);
        throw("Error: failed writing sparse source LD local-to-global index map.");
    }
    return countPos;
}

int16_t quantizeSparseLdCorrQ16(float r) {
    if (r > 1.0f) r = 1.0f;
    else if (r < -1.0f) r = -1.0f;
    int q = (int)(r >= 0.0f ? r * SOURCE_SPARSE_LD_Q16_SCALE + 0.5f
                             : r * SOURCE_SPARSE_LD_Q16_SCALE - 0.5f);
    if (q > 32767) q = 32767;
    else if (q < -32767) q = -32767;
    return (int16_t)q;
}

void writeSparseSourceLdV2Record(FILE *out, const uint32_t snp1, const uint32_t snp2,
                                 const float r, const bool useUint16Index,
                                 const bool q16Correlation, uint64_t &numRecords) {
    float clipped = r;
    if (clipped > 1.0f) clipped = 1.0f;
    else if (clipped < -1.0f) clipped = -1.0f;
    if (useUint16Index) {
        uint16_t i16 = (uint16_t)snp1;
        uint16_t j16 = (uint16_t)snp2;
        fwrite(&i16, sizeof(uint16_t), 1, out);
        fwrite(&j16, sizeof(uint16_t), 1, out);
    } else {
        fwrite(&snp1, sizeof(uint32_t), 1, out);
        fwrite(&snp2, sizeof(uint32_t), 1, out);
    }
    if (q16Correlation) {
        int16_t q = quantizeSparseLdCorrQ16(clipped);
        fwrite(&q, sizeof(int16_t), 1, out);
    } else {
        fwrite(&clipped, sizeof(float), 1, out);
    }
    ++numRecords;
}

void finishSparseSourceLdV2(FILE *out, const long countPos, const uint64_t numRecords) {
    fseek(out, countPos, SEEK_SET);
    fwrite(&numRecords, sizeof(uint64_t), 1, out);
    fclose(out);
}

void writeBinaryString(FILE *out, const string &value) {
    uint32_t len = (uint32_t)value.size();
    fwrite(&len, sizeof(uint32_t), 1, out);
    if (len && fwrite(value.data(), sizeof(char), len, out) != len) {
        throw("Error: failed writing binary string.");
    }
}

string readBinaryString(FILE *in, const string &filename) {
    uint32_t len = 0;
    if (fread(&len, sizeof(uint32_t), 1, in) != 1) {
        throw("Error: incomplete binary string length in [" + filename + "].");
    }
    string value(len, '\0');
    if (len && fread(&value[0], sizeof(char), len, in) != len) {
        throw("Error: incomplete binary string in [" + filename + "].");
    }
    return value;
}

void prepareCleanOutputDirectory(const string &dirname, const string &label) {
    struct stat sb;
    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        if (mkdir(dirname.c_str(), 0755) != 0) {
            throw("Error: can not create folder [" + dirname + "] to store " + label + ".");
        }
        cout << "Created folder [" << dirname << "] to store " << label << "." << endl;
        return;
    }

    DIR *existingDir = opendir(dirname.c_str());
    if (!existingDir) throw("Error: can not open folder [" + dirname + "].");
    dirent *entry;
    while ((entry = readdir(existingDir)) != NULL) {
        string fname = entry->d_name;
        if (fname == "." || fname == "..") continue;
        string path = dirname + "/" + fname;
        struct stat fstat;
        if (stat(path.c_str(), &fstat) != 0) continue;
        if (S_ISDIR(fstat.st_mode)) {
            closedir(existingDir);
            throw("Error: output folder [" + dirname + "] contains subfolder [" + fname + "]. "
                  "Please remove it or use a different --out folder.");
        }
        if (remove(path.c_str()) != 0) {
            closedir(existingDir);
            throw("Error: failed to remove existing file [" + path + "] before rewriting output.");
        }
    }
    closedir(existingDir);
    cout << "Cleared existing files in folder [" << dirname << "] before writing " << label << "." << endl;
}

void prepareOutputDirectoryPreserveExisting(const string &dirname, const string &label) {
    struct stat sb;
    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        if (mkdir(dirname.c_str(), 0755) != 0) {
            throw("Error: can not create folder [" + dirname + "] to store " + label + ".");
        }
        cout << "Created folder [" << dirname << "] to store " << label << "." << endl;
        return;
    }
    cout << "Using existing folder [" << dirname << "] and preserving files from other blocks." << endl;
}

class SourceBlockOutputLock {
    string lockDir;
    bool locked;
public:
    SourceBlockOutputLock(const string &dirname, const string &sourceBlockID, const bool enabled)
    : lockDir(dirname + "/block" + sourceBlockID + ".lock"), locked(false) {
        if (!enabled) return;
        if (mkdir(lockDir.c_str(), 0755) != 0) {
            throw("Error: block-specific output lock [" + lockDir + "] already exists. "
                  "Another job may be writing this block, or a previous job failed. "
                  "Remove this lock directory only after confirming no job is still running for this block.");
        }
        locked = true;
    }

    ~SourceBlockOutputLock() {
        if (locked) remove(lockDir.c_str());
    }
};

void removeFileIfExistsLocal(const string &filename) {
    struct stat sb;
    if (stat(filename.c_str(), &sb) == 0 && S_ISREG(sb.st_mode)) remove(filename.c_str());
}

bool hasSuffixLocal(const string &value, const string &suffix) {
    return value.size() >= suffix.size()
        && value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool findBlockOutputLocksLocal(const string &dirname, vector<string> &locks) {
    locks.clear();
    DIR *dirp = opendir(dirname.c_str());
    if (!dirp) throw("Error: can not open folder [" + dirname + "] to check block output locks.");
    dirent *dp;
    while ((dp = readdir(dirp)) != NULL) {
        string fname = dp->d_name;
        if (fname.find("block") != 0 || !hasSuffixLocal(fname, ".lock")) continue;
        string path = dirname + "/" + fname;
        struct stat sb;
        if (stat(path.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) locks.push_back(fname);
    }
    closedir(dirp);
    std::sort(locks.begin(), locks.end());
    return !locks.empty();
}

bool isBlockSpecificMergeInputLocal(const string &fname) {
    if (fname.find("block") != 0) return false;
    return hasSuffixLocal(fname, ".snp.info")
        || hasSuffixLocal(fname, ".ldm.info")
        || hasSuffixLocal(fname, ".rsq0.5.pwld")
        || hasSuffixLocal(fname, ".rsq0.5.pwld.gz")
        || hasSuffixLocal(fname, ".rsq0.5.pwld.bin")
        || hasSuffixLocal(fname, ".pwld.source");
}

unsigned removeBlockSpecificMergeInputsLocal(const string &dirname) {
    DIR *dirp = opendir(dirname.c_str());
    if (!dirp) throw("Error: can not open folder [" + dirname + "] to clean block-specific merge inputs.");
    vector<string> files;
    dirent *dp;
    while ((dp = readdir(dirp)) != NULL) {
        string fname = dp->d_name;
        if (isBlockSpecificMergeInputLocal(fname)) files.push_back(fname);
    }
    closedir(dirp);

    unsigned removed = 0;
    for (unsigned i = 0; i < files.size(); ++i) {
        string path = dirname + "/" + files[i];
        struct stat sb;
        if (stat(path.c_str(), &sb) == 0 && S_ISREG(sb.st_mode) && remove(path.c_str()) == 0)
            ++removed;
    }
    return removed;
}

void removeSourceBlockTiledOutputs(const string &dirname, const string &sourceBlockID) {
    DIR *dirp = opendir(dirname.c_str());
    if (!dirp) throw("Error: can not open folder [" + dirname + "] before rewriting source block " + sourceBlockID + ".");
    const string tilePrefix = "block" + sourceBlockID + ".tile";
    dirent *dp;
    while ((dp = readdir(dirp)) != NULL) {
        string fname = dp->d_name;
        if (fname.find(tilePrefix) != 0) continue;
        if (fname.size() >= 10 && fname.substr(fname.size() - 10) == ".eigen.bin")
            remove((dirname + "/" + fname).c_str());
        else if (fname.size() >= 14 && fname.substr(fname.size() - 14) == ".eigen.q16.bin")
            remove((dirname + "/" + fname).c_str());
        else if (fname.size() >= 10 && fname.substr(fname.size() - 10) == ".eigen.txt")
            remove((dirname + "/" + fname).c_str());
    }
    closedir(dirp);

    const string prefix = dirname + "/block" + sourceBlockID;
    removeFileIfExistsLocal(prefix + ".ldm.info");
    removeFileIfExistsLocal(prefix + ".snp.info");
    removeFileIfExistsLocal(prefix + ".sparse.bin");
    removeFileIfExistsLocal(prefix + ".sparse.q16.bin");
    removeFileIfExistsLocal(prefix + ".tiles.eigen.bin");
    removeFileIfExistsLocal(prefix + ".tiles.eigen.q16.bin");
    removeFileIfExistsLocal(prefix + ".tiles.mctx.bin");
    removeFileIfExistsLocal(prefix + ".rsq0.5.pwld");
    removeFileIfExistsLocal(prefix + ".rsq0.5.pwld.gz");
    removeFileIfExistsLocal(prefix + ".rsq0.5.pwld.bin");
    removeFileIfExistsLocal(prefix + ".pwld.source");
    removeFileIfExistsLocal(prefix + ".done");
}

void writeEigenBlockPayload(FILE *outbin, int32_t numSnpInBlock, const VectorXf &eigenVal,
                            const MatrixXf &eigenVec, float sumPosEigVal, float eigenCutoff) {
    int32_t numEigenValue = eigenVal.size();
    fwrite(&numSnpInBlock, sizeof(int32_t), 1, outbin);
    fwrite(&numEigenValue, sizeof(int32_t), 1, outbin);
    fwrite(&sumPosEigVal, sizeof(float), 1, outbin);
    fwrite(&eigenCutoff, sizeof(float), 1, outbin);
    fwrite(eigenVal.data(), sizeof(float), numEigenValue, outbin);
    uint64_t nElements = (uint64_t)numSnpInBlock * (uint64_t)numEigenValue;
    fwrite(eigenVec.data(), sizeof(float), nElements, outbin);
}

void writeEigenBlockBinary(const string &filename, int32_t numSnpInBlock, const VectorXf &eigenVal,
                           const MatrixXf &eigenVec, float sumPosEigVal, float eigenCutoff) {
    FILE *outbin = fopen(filename.c_str(), "wb");
    if (!outbin) throw("Error: can not open the file [" + filename + "] to write.");
    writeEigenBlockPayload(outbin, numSnpInBlock, eigenVal, eigenVec, sumPosEigVal, eigenCutoff);
    fclose(outbin);
}

void writeEigenBlockPayloadQ16(FILE *outbin, int32_t numSnpInBlock, const VectorXf &eigenVal,
                               const MatrixXf &eigenVec, float sumPosEigVal, float eigenCutoff) {
    int32_t numEigenValue = eigenVal.size();
    fwrite(&numSnpInBlock, sizeof(int32_t), 1, outbin);
    fwrite(&numEigenValue, sizeof(int32_t), 1, outbin);
    fwrite(&sumPosEigVal, sizeof(float), 1, outbin);
    fwrite(&eigenCutoff, sizeof(float), 1, outbin);
    fwrite(eigenVal.data(), sizeof(float), numEigenValue, outbin);

    VectorXf scales(numEigenValue);
    for (int col = 0; col < numEigenValue; ++col) {
        float maxAbs = eigenVec.col(col).cwiseAbs().maxCoeff();
        if (maxAbs <= 0.0f) maxAbs = 1.0f;
        scales[col] = maxAbs;
    }
    fwrite(scales.data(), sizeof(float), numEigenValue, outbin);

    vector<int16_t> colRaw(numSnpInBlock);
    for (int col = 0; col < numEigenValue; ++col) {
        const float invScale = SOURCE_SPARSE_LD_Q16_SCALE / scales[col];
        for (int row = 0; row < numSnpInBlock; ++row) {
            int q = (int)lrintf(eigenVec(row, col) * invScale);
            if (q > 32767) q = 32767;
            else if (q < -32767) q = -32767;
            colRaw[row] = (int16_t)q;
        }
        fwrite(colRaw.data(), sizeof(int16_t), numSnpInBlock, outbin);
    }
}

void writeEigenBlockBinaryQ16(const string &filename, int32_t numSnpInBlock, const VectorXf &eigenVal,
                              const MatrixXf &eigenVec, float sumPosEigVal, float eigenCutoff) {
    FILE *outbin = fopen(filename.c_str(), "wb");
    if (!outbin) throw("Error: can not open the file [" + filename + "] to write.");
    writeEigenBlockPayloadQ16(outbin, numSnpInBlock, eigenVal, eigenVec, sumPosEigVal, eigenCutoff);
    fclose(outbin);
}

string tiledEigenBundleFilename(const string &dirname, const string &sourceBlockID, const int outputQuantBits) {
    if (outputQuantBits == 0) return dirname + "/block" + sourceBlockID + ".tiles.eigen.bin";
    if (outputQuantBits == 16) return dirname + "/block" + sourceBlockID + ".tiles.eigen.q16.bin";
    throw("Error: tiled LD eigen bundles currently support only float or --bits 16.");
}

string marginalizedContextBundleFilename(const string &dirname, const string &sourceBlockID) {
    return dirname + "/block" + sourceBlockID + ".tiles.mctx.bin";
}

string sparseSourceLdFilename(const string &dirname, const string &sourceBlockID, const bool q16Correlation) {
    return dirname + "/block" + sourceBlockID + (q16Correlation ? ".sparse.q16.bin" : ".sparse.bin");
}

string existingSparseSourceLdFilename(const string &dirname, const string &sourceBlockID) {
    string q16File = sparseSourceLdFilename(dirname, sourceBlockID, true);
    if (fileExistsAndNonEmptyLocal(q16File)) return q16File;
    string floatFile = sparseSourceLdFilename(dirname, sourceBlockID, false);
    if (fileExistsAndNonEmptyLocal(floatFile)) return floatFile;
    return floatFile;
}

FILE *openTiledEigenBundleForWrite(const string &dirname, const string &sourceBlockID,
                                   const unsigned numRecords, const int outputQuantBits) {
    const string filename = tiledEigenBundleFilename(dirname, sourceBlockID, outputQuantBits);
    FILE *outbin = fopen(filename.c_str(), "wb");
    if (!outbin) throw("Error: can not open the file [" + filename + "] to write.");
    fwrite(TILED_EIGEN_BUNDLE_MAGIC, sizeof(char), 8, outbin);
    fwrite(&TILED_EIGEN_BUNDLE_VERSION, sizeof(uint32_t), 1, outbin);
    uint32_t bits = (uint32_t)outputQuantBits;
    uint32_t nrec = (uint32_t)numRecords;
    fwrite(&bits, sizeof(uint32_t), 1, outbin);
    fwrite(&nrec, sizeof(uint32_t), 1, outbin);
    return outbin;
}

void writeTiledEigenBundleRecord(FILE *outbin, const string &blockID, const int32_t numSnpInBlock,
                                 const VectorXf &eigenVal, const MatrixXf &eigenVec,
                                 const float sumPosEigVal, const float eigenCutoff,
                                 const int outputQuantBits) {
    writeBinaryString(outbin, blockID);
    uint64_t payloadBytes = 0;
    const long payloadSizePos = ftell(outbin);
    fwrite(&payloadBytes, sizeof(uint64_t), 1, outbin);
    const long payloadStart = ftell(outbin);
    if (outputQuantBits == 0)
        writeEigenBlockPayload(outbin, numSnpInBlock, eigenVal, eigenVec, sumPosEigVal, eigenCutoff);
    else if (outputQuantBits == 16)
        writeEigenBlockPayloadQ16(outbin, numSnpInBlock, eigenVal, eigenVec, sumPosEigVal, eigenCutoff);
    else
        throw("Error: tiled LD eigen bundles currently support only float or --bits 16.");
    const long payloadEnd = ftell(outbin);
    payloadBytes = (uint64_t)(payloadEnd - payloadStart);
    fseek(outbin, payloadSizePos, SEEK_SET);
    fwrite(&payloadBytes, sizeof(uint64_t), 1, outbin);
    fseek(outbin, payloadEnd, SEEK_SET);
}

FILE *openMarginalizedContextBundleForWrite(const string &dirname, const string &sourceBlockID,
                                             const unsigned numRecords) {
    const string filename = marginalizedContextBundleFilename(dirname, sourceBlockID);
    FILE *outbin = fopen(filename.c_str(), "wb");
    if (!outbin) throw("Error: can not open the file [" + filename + "] to write.");
    fwrite(MARGINAL_CONTEXT_BUNDLE_MAGIC, sizeof(char), 8, outbin);
    fwrite(&MARGINAL_CONTEXT_BUNDLE_VERSION, sizeof(uint32_t), 1, outbin);
    const uint32_t nrec = (uint32_t)numRecords;
    fwrite(&nrec, sizeof(uint32_t), 1, outbin);
    return outbin;
}

void writeMarginalizedContextBundleRecord(FILE *outbin, const string &blockID,
                                           const VectorXf &lambda, const VectorXf &contextEigenvalues,
                                           const MatrixXf &rotatedWOperator, const MatrixXf &rotatedQ,
                                           const uint64_t numLinks, const uint32_t numFlankSnps) {
    const int32_t rank = (int32_t)lambda.size();
    const int32_t numSnps = (int32_t)rotatedQ.cols();
    if (rank <= 0 || rotatedQ.rows() != rank || contextEigenvalues.size() != rank
        || rotatedWOperator.rows() != rank || rotatedWOperator.cols() != rank) {
        throw("Error: inconsistent marginalized-context dimensions for tile " + blockID + ".");
    }

    writeBinaryString(outbin, blockID);
    uint64_t payloadBytes = 0;
    const long payloadSizePos = ftell(outbin);
    fwrite(&payloadBytes, sizeof(uint64_t), 1, outbin);
    const long payloadStart = ftell(outbin);

    fwrite(&numSnps, sizeof(int32_t), 1, outbin);
    fwrite(&rank, sizeof(int32_t), 1, outbin);
    fwrite(&numLinks, sizeof(uint64_t), 1, outbin);
    fwrite(&numFlankSnps, sizeof(uint32_t), 1, outbin);
    fwrite(lambda.data(), sizeof(float), rank, outbin);
    fwrite(contextEigenvalues.data(), sizeof(float), rank, outbin);
    fwrite(rotatedWOperator.data(), sizeof(float), (uint64_t)rank * rank, outbin);
    fwrite(rotatedQ.data(), sizeof(float), (uint64_t)rank * numSnps, outbin);

    const long payloadEnd = ftell(outbin);
    payloadBytes = (uint64_t)(payloadEnd - payloadStart);
    fseek(outbin, payloadSizePos, SEEK_SET);
    fwrite(&payloadBytes, sizeof(uint64_t), 1, outbin);
    fseek(outbin, payloadEnd, SEEK_SET);
}

void writeTiledEigenBlockBinary(const string &dirname, const string &blockID, const int32_t numSnpInBlock,
                                const VectorXf &eigenVal, const MatrixXf &eigenVec,
                                const float sumPosEigVal, const float eigenCutoff,
                                const int outputQuantBits) {
    if (outputQuantBits == 0) {
        writeEigenBlockBinary(dirname + "/block" + blockID + ".eigen.bin",
                              numSnpInBlock, eigenVal, eigenVec, sumPosEigVal, eigenCutoff);
    } else if (outputQuantBits == 16) {
        writeEigenBlockBinaryQ16(dirname + "/block" + blockID + ".eigen.q16.bin",
                                 numSnpInBlock, eigenVal, eigenVec, sumPosEigVal, eigenCutoff);
    } else {
        throw("Error: tiled LD eigen output currently supports only float or --bits 16.");
    }
}

void writeCoreContextLdBlockBinary(const string &filename, const CoreContextLdBlock &block) {
    FILE *outbin = fopen(filename.c_str(), "wb");
    if (!outbin) throw("Error: can not open the file [" + filename + "] to write.");
    fwrite(CORE_CONTEXT_BIN_MAGIC, sizeof(char), 8, outbin);
    fwrite(&CORE_CONTEXT_BIN_VERSION, sizeof(uint32_t), 1, outbin);
    int32_t numCore = block.numCoreSnps;
    uint64_t numRecords = block.coreLocalIndex.size();
    fwrite(&numCore, sizeof(int32_t), 1, outbin);
    fwrite(&numRecords, sizeof(uint64_t), 1, outbin);
    for (uint64_t r = 0; r < numRecords; ++r) {
        CoreContextBinRecord rec;
        rec.coreLocal = (uint32_t)block.coreLocalIndex[r];
        rec.flankGlobal = (uint32_t)block.flankGlobalIndex[r];
        rec.ldcor = block.ldCorr[r];
        fwrite(&rec, sizeof(CoreContextBinRecord), 1, outbin);
    }
    fclose(outbin);
}

void writePwldGzFromDense(const string &filename, const vector<SnpInfo*> &sourceSnps,
                          const MatrixXf &ldm, const float r2Threshold) {
    gzFile out = gzopen(filename.c_str(), "wb");
    if (!out) throw("Error: can not open the file [" + filename + "] to write.");
    gzprintf(out, "%12s %12s %12s\n", "SNP1", "SNP2", "LDcorrelation");
    const int m = ldm.cols();
    if (ldm.rows() != m || m != (int)sourceSnps.size()) {
        gzclose(out);
        throw("Error: inconsistent dimensions while writing pairwise LD file [" + filename + "].");
    }
    const float threshold = std::max(0.0f, r2Threshold);
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < j; ++i) {
            float r = ldm(i, j);
            if (r > 1.0f) r = 1.0f;
            else if (r < -1.0f) r = -1.0f;
            if (r * r < threshold) continue;
            gzprintf(out, "%12s %12s %12.6f\n",
                     sourceSnps[i]->ID.c_str(), sourceSnps[j]->ID.c_str(), r);
        }
    }
    gzclose(out);
}

bool readCoreContextLdBlockBinary(const string &filename, CoreContextLdBlock &block) {
    block.clear();
    FILE *in = fopen(filename.c_str(), "rb");
    if (!in) return false;
    char magic[8];
    uint32_t version = 0;
    int32_t numCore = 0;
    uint64_t numRecords = 0;
    if (fread(magic, sizeof(char), 8, in) != 8 ||
        memcmp(magic, CORE_CONTEXT_BIN_MAGIC, 8) != 0 ||
        fread(&version, sizeof(uint32_t), 1, in) != 1 ||
        fread(&numCore, sizeof(int32_t), 1, in) != 1 ||
        fread(&numRecords, sizeof(uint64_t), 1, in) != 1) {
        fclose(in);
        throw("Error: invalid legacy core-context sparse LD file [" + filename + "].");
    }
    if (version != CORE_CONTEXT_BIN_VERSION) {
        fclose(in);
        throw("Error: unsupported legacy core-context sparse LD file version in [" + filename + "].");
    }
    block.numCoreSnps = numCore;
    block.coreLocalIndex.reserve(numRecords);
    block.flankGlobalIndex.reserve(numRecords);
    block.ldCorr.reserve(numRecords);
    for (uint64_t r = 0; r < numRecords; ++r) {
        CoreContextBinRecord rec;
        if (fread(&rec, sizeof(CoreContextBinRecord), 1, in) != 1) {
            fclose(in);
            throw("Error: incomplete core-context LD record in [" + filename + "].");
        }
        block.coreLocalIndex.push_back((unsigned)rec.coreLocal);
        block.flankGlobalIndex.push_back((unsigned)rec.flankGlobal);
        block.ldCorr.push_back(rec.ldcor);
    }
    fclose(in);
    return true;
}

void writeCoreContextTileMetaBinary(const string &filename, const string &sourceBlockID) {
    FILE *out = fopen(filename.c_str(), "wb");
    if (!out) throw("Error: can not open the file [" + filename + "] to write.");
    fwrite(CORE_CONTEXT_META_MAGIC, sizeof(char), 8, out);
    fwrite(&CORE_CONTEXT_BIN_VERSION, sizeof(uint32_t), 1, out);
    writeBinaryString(out, sourceBlockID);
    fclose(out);
}

bool readCoreContextTileMetaBinary(const string &filename, string &sourceBlockID) {
    FILE *in = fopen(filename.c_str(), "rb");
    if (!in) return false;
    char magic[8];
    uint32_t version = 0;
    if (fread(magic, sizeof(char), 8, in) != 8 ||
        memcmp(magic, CORE_CONTEXT_META_MAGIC, 8) != 0 ||
        fread(&version, sizeof(uint32_t), 1, in) != 1) {
        fclose(in);
        throw("Error: invalid core-context tile metadata [" + filename + "].");
    }
    if (version != CORE_CONTEXT_BIN_VERSION) {
        fclose(in);
        throw("Error: unsupported core-context tile metadata version in [" + filename + "].");
    }
    sourceBlockID = readBinaryString(in, filename);
    fclose(in);
    return true;
}

uint64_t writeSparseSourceLdFromQ(const string &filename, const string &sourceBlockID,
                                  const vector<SnpInfo*> &sourceSnps,
                                  const unordered_map<int, int> &newIndexByOldIndex,
                                  const vector<int> &tileByLocalIndex,
                                  const MatrixXf &Qsource, const VectorXf &qDiag,
                                  const float r2Threshold, const bool q16Correlation) {
    FILE *out = fopen(filename.c_str(), "wb");
    if (!out) throw("Error: can not open the file [" + filename + "] to write.");

    const int m = Qsource.cols();
    if (m != (int)sourceSnps.size() || qDiag.size() != m || (int)tileByLocalIndex.size() != m) {
        fclose(out);
        throw("Error: inconsistent dimensions while writing sparse source LD for block " + sourceBlockID + ".");
    }
    vector<uint32_t> localToGlobal = makeSparseSourceLocalToGlobal(sourceSnps, newIndexByOldIndex);
    const bool useUint16Index = m <= 65536;
    uint64_t numRecords = 0;
    const long countPos = writeSparseSourceLdHeader(out, sourceBlockID, localToGlobal, useUint16Index, q16Correlation);

    MatrixXf Qnorm(Qsource.rows(), Qsource.cols());
    for (int j = 0; j < m; ++j) {
        const float d = qDiag[j];
        if (d > 0.0f) Qnorm.col(j) = Qsource.col(j) / sqrtf(d);
        else Qnorm.col(j).setZero();
    }

    const int chunk = 256;
    const float threshold = std::max(0.0f, r2Threshold);
    for (int c0 = 0; c0 < m; c0 += chunk) {
        const int width = std::min(chunk, m - c0);
        MatrixXf R = Qnorm.transpose() * Qnorm.middleCols(c0, width);
        for (int jj = 0; jj < width; ++jj) {
            const int j = c0 + jj;
            if (localToGlobal[j] == SOURCE_SPARSE_LD_MISSING_INDEX) continue;
            for (int i = 0; i < j; ++i) {
                if (tileByLocalIndex[i] == tileByLocalIndex[j]) continue;
                if (localToGlobal[i] == SOURCE_SPARSE_LD_MISSING_INDEX) continue;
                float r = R(i, jj);
                if (r > 1.0f) r = 1.0f;
                else if (r < -1.0f) r = -1.0f;
                if (r * r < threshold) continue;
                writeSparseSourceLdV2Record(out, (uint32_t)i, (uint32_t)j, r, useUint16Index, q16Correlation, numRecords);
            }
        }
    }

    finishSparseSourceLdV2(out, countPos, numRecords);
    return numRecords;
}

uint64_t writeSparseSourceLdFromDense(const string &filename, const string &sourceBlockID,
                                      const vector<SnpInfo*> &sourceSnps,
                                      const unordered_map<int, int> &newIndexByOldIndex,
                                      const vector<int> &tileByLocalIndex,
                                      const MatrixXf &ldm,
                                      const float r2Threshold, const bool q16Correlation) {
    FILE *out = fopen(filename.c_str(), "wb");
    if (!out) throw("Error: can not open the file [" + filename + "] to write.");

    const int m = ldm.cols();
    if (ldm.rows() != m || m != (int)sourceSnps.size() || (int)tileByLocalIndex.size() != m) {
        fclose(out);
        throw("Error: inconsistent dimensions while writing sparse source LD for block " + sourceBlockID + ".");
    }
    vector<uint32_t> localToGlobal = makeSparseSourceLocalToGlobal(sourceSnps, newIndexByOldIndex);
    const bool useUint16Index = m <= 65536;
    uint64_t numRecords = 0;
    const long countPos = writeSparseSourceLdHeader(out, sourceBlockID, localToGlobal, useUint16Index, q16Correlation);

    const float threshold = std::max(0.0f, r2Threshold);
    for (int j = 0; j < m; ++j) {
        if (localToGlobal[j] == SOURCE_SPARSE_LD_MISSING_INDEX) continue;
        for (int i = 0; i < j; ++i) {
            if (tileByLocalIndex[i] == tileByLocalIndex[j]) continue;
            if (localToGlobal[i] == SOURCE_SPARSE_LD_MISSING_INDEX) continue;
            float r = ldm(i, j);
            if (r > 1.0f) r = 1.0f;
            else if (r < -1.0f) r = -1.0f;
            if (r * r < threshold) continue;
            writeSparseSourceLdV2Record(out, (uint32_t)i, (uint32_t)j, r, useUint16Index, q16Correlation, numRecords);
        }
    }

    finishSparseSourceLdV2(out, countPos, numRecords);
    return numRecords;
}

struct BlockLdBinaryInfo {
    string filename;
    int quantBits;
    bool lowerTriangular;
    bool fullSquare;
    BlockLdBinaryInfo() : quantBits(0), lowerTriangular(false), fullSquare(false) {}
};

bool inspectBlockLdBinary(const string &dirname, const string &blockID, const int32_t blockSize,
                          BlockLdBinaryInfo &info) {
    const string floatFile = dirname + "/block" + blockID + ".ldm.bin";
    const string q16File = dirname + "/block" + blockID + ".ldm.q16.bin";
    const string q8File = dirname + "/block" + blockID + ".ldm.q8.bin";

    struct stat statBuf;
    if (stat(floatFile.c_str(), &statBuf) == 0) {
        info.filename = floatFile;
        info.quantBits = 0;
    } else if (stat(q16File.c_str(), &statBuf) == 0) {
        info.filename = q16File;
        info.quantBits = 16;
    } else if (stat(q8File.c_str(), &statBuf) == 0) {
        info.filename = q8File;
        info.quantBits = 8;
    } else {
        return false;
    }

    const size_t nFull = (size_t)blockSize * (size_t)blockSize;
    const size_t nLowerTri = (size_t)blockSize * (size_t)(blockSize + 1) / 2;
    const size_t bytesPerValue = info.quantBits == 0 ? sizeof(float) : (info.quantBits == 8 ? sizeof(int8_t) : sizeof(int16_t));
    const size_t fileBytes = (size_t)statBuf.st_size;
    info.lowerTriangular = fileBytes == nLowerTri * bytesPerValue;
    info.fullSquare = fileBytes == nFull * bytesPerValue;
    return true;
}

bool buildCoreTileLdAndSparseFromLowerTriLdStream(const string &ldmDirname,
                                                  const string &outDirname,
                                                  const string &sourceBlockID,
                                                  const int32_t mFull,
                                                  const vector<unsigned> &sourceLocal,
                                                  const vector<SnpInfo*> &sourceSnps,
                                                  const unordered_map<int, int> &newIndexByOldIndex,
                                                  const vector<vector<unsigned> > &tileLocalIdx,
                                                  const vector<int> &tileByLocalIndex,
                                                  const float sparseR2Threshold,
                                                  const int outputQuantBits,
                                                  vector<MatrixXf> &tileRcc,
                                                  uint64_t &sparseRecordsOut) {
    sparseRecordsOut = 0;
    BlockLdBinaryInfo info;
    if (!inspectBlockLdBinary(ldmDirname, sourceBlockID, mFull, info)) {
        throw("Error: cannot find block LD matrix file for source block " + sourceBlockID + " in [" + ldmDirname + "].");
    }
    if (!info.lowerTriangular) return false;

    FILE *in = fopen(info.filename.c_str(), "rb");
    if (!in) throw("Error: cannot open the file [" + info.filename + "] to read.");

    const int m = (int)sourceSnps.size();
    vector<int> originalToSource(mFull, -1);
    for (int s = 0; s < m; ++s) originalToSource[(int)sourceLocal[s]] = s;

    vector<int> sourceToTileLocal(m, -1);
    for (unsigned t = 0; t < tileLocalIdx.size(); ++t) {
        for (unsigned a = 0; a < tileLocalIdx[t].size(); ++a)
            sourceToTileLocal[(int)tileLocalIdx[t][a]] = (int)a;
    }
    if ((int)tileByLocalIndex.size() != m) {
        fclose(in);
        throw("Error: inconsistent tile assignment while streaming source block " + sourceBlockID + ".");
    }
    for (int s = 0; s < m; ++s) {
        if (tileByLocalIndex[s] < 0 || sourceToTileLocal[s] < 0) {
            fclose(in);
            throw("Error: SNPs are not fully assigned to core tiles in source block " + sourceBlockID + ".");
        }
    }

    tileRcc.resize(tileLocalIdx.size());
    for (unsigned t = 0; t < tileLocalIdx.size(); ++t)
        tileRcc[t].setZero((int)tileLocalIdx[t].size(), (int)tileLocalIdx[t].size());

    const bool q16Sparse = outputQuantBits == 16;
    vector<uint32_t> localToGlobal = makeSparseSourceLocalToGlobal(sourceSnps, newIndexByOldIndex);
    const bool useUint16Index = m <= 65536;
    uint64_t sparseRecords = 0;
    const string sparseFile = sparseSourceLdFilename(outDirname, sourceBlockID, q16Sparse);
    FILE *outSparse = fopen(sparseFile.c_str(), "wb");
    if (!outSparse) {
        fclose(in);
        throw("Error: can not open the file [" + sparseFile + "] to write.");
    }
    const long sparseCountPos = writeSparseSourceLdHeader(outSparse, sourceBlockID, localToGlobal, useUint16Index, q16Sparse);

    const string pwldFile = outDirname + "/block" + sourceBlockID + ".rsq0.5.pwld.gz";
    gzFile outPwld = gzopen(pwldFile.c_str(), "wb");
    if (!outPwld) {
        fclose(in);
        fclose(outSparse);
        throw("Error: can not open the file [" + pwldFile + "] to write.");
    }
    gzprintf(outPwld, "%12s %12s %12s\n", "SNP1", "SNP2", "LDcorrelation");

    const float sparseThreshold = std::max(0.0f, sparseR2Threshold);
    const float pwldThreshold = 0.5f;
    const float invQ = info.quantBits == 8 ? (1.0f / 127.0f)
                     : info.quantBits == 16 ? (1.0f / 32767.0f)
                     : 1.0f;

    vector<float> floatRow;
    vector<int16_t> q16Row;
    vector<int8_t> q8Row;
    for (int row = 0; row < mFull; ++row) {
        const int rowLen = row + 1;
        if (info.quantBits == 0) {
            floatRow.resize(rowLen);
            if (fread(floatRow.data(), sizeof(float), rowLen, in) != (size_t)rowLen) {
                gzclose(outPwld);
                fclose(outSparse);
                fclose(in);
                throw("Error: incomplete lower-triangular LD row in [" + info.filename + "].");
            }
        } else if (info.quantBits == 16) {
            q16Row.resize(rowLen);
            if (fread(q16Row.data(), sizeof(int16_t), rowLen, in) != (size_t)rowLen) {
                gzclose(outPwld);
                fclose(outSparse);
                fclose(in);
                throw("Error: incomplete q16 lower-triangular LD row in [" + info.filename + "].");
            }
        } else {
            q8Row.resize(rowLen);
            if (fread(q8Row.data(), sizeof(int8_t), rowLen, in) != (size_t)rowLen) {
                gzclose(outPwld);
                fclose(outSparse);
                fclose(in);
                throw("Error: incomplete q8 lower-triangular LD row in [" + info.filename + "].");
            }
        }

        const int si = originalToSource[row];
        if (si < 0) continue;
        const int tileI = tileByLocalIndex[si];
        const int tileLocalI = sourceToTileLocal[si];
        for (int col = 0; col <= row; ++col) {
            const int sj = originalToSource[col];
            if (sj < 0) continue;
            float r = 0.0f;
            if (info.quantBits == 0) r = floatRow[col];
            else if (info.quantBits == 16) r = (float)q16Row[col] * invQ;
            else r = (float)q8Row[col] * invQ;
            if (r > 1.0f) r = 1.0f;
            else if (r < -1.0f) r = -1.0f;

            const int tileJ = tileByLocalIndex[sj];
            const int tileLocalJ = sourceToTileLocal[sj];
            const float r2 = r * r;
            if (row != col && r2 >= pwldThreshold)
                gzprintf(outPwld, "%12s %12s %12.6f\n", sourceSnps[sj]->ID.c_str(), sourceSnps[si]->ID.c_str(), r);
            if (tileI == tileJ) {
                MatrixXf &Rcc = tileRcc[tileI];
                Rcc(tileLocalI, tileLocalJ) = r;
                Rcc(tileLocalJ, tileLocalI) = r;
                continue;
            }
            if (r2 >= sparseThreshold)
                writeSparseSourceLdV2Record(outSparse, (uint32_t)sj, (uint32_t)si, r, useUint16Index, q16Sparse, sparseRecords);
        }
    }

    gzclose(outPwld);
    sparseRecordsOut = sparseRecords;
    finishSparseSourceLdV2(outSparse, sparseCountPos, sparseRecords);
    fclose(in);

    for (unsigned t = 0; t < tileRcc.size(); ++t)
        Data::normaliseDenseCorrelation(tileRcc[t]);

    return true;
}

void outputCoreContextSourceInfo(const vector<LDBlockInfo*> &tileBlocks, const string &dirname,
                                 const string &sourceBlockID, const int numKeptInds) {
    string outSnpfile = dirname + "/block" + sourceBlockID + ".snp.info";
    string outldmfile = dirname + "/block" + sourceBlockID + ".ldm.info";

    ofstream out1(outSnpfile.c_str());
    out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %12s %10s %15s\n")
    % "Chrom"
    % "ID"
    % "Index"
    % "GenPos"
    % "PhysPos"
    % "A1"
    % "A2"
    % "A1Freq"
    % "N"
    % "Block";
    for (unsigned b = 0; b < tileBlocks.size(); ++b) {
        const LDBlockInfo *block = tileBlocks[b];
        for (unsigned i = 0; i < block->snpInfoVec.size(); ++i) {
            SnpInfo *snp = block->snpInfoVec[i];
            out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %12f %10s %15s\n")
            % snp->chrom
            % snp->ID
            % snp->index
            % snp->genPos
            % snp->physPos
            % snp->a1
            % snp->a2
            % snp->af
            % numKeptInds
            % block->ID;
        }
    }
    out1.close();

    ofstream out2(outldmfile.c_str());
    out2 << boost::format("%15s %6s %15s %15s %15s %15s %12s %15s %10s\n")
    % "Block" % "Chrom" % "StartSnpIdx" % "StartSnpID" % "EndSnpIdx" % "EndSnpID" % "NumSnps"
    % "SourceBlockID" % "TileID";
    for (unsigned b = 0; b < tileBlocks.size(); ++b) {
        const LDBlockInfo *block = tileBlocks[b];
        out2 << boost::format("%15s %6s %15s %15s %15s %15s %12s %15s %10s\n")
        % block->ID
        % block->chrom
        % block->startSnpIdx
        % block->snpInfoVec.front()->ID
        % block->endSnpIdx
        % block->snpInfoVec.back()->ID
        % block->numSnpInBlock
        % block->sourceBlockID
        % block->sourceTileID;
    }
    out2.close();
}

bool readSparseSourceLdBinary(const string &filename, vector<SparseSourceLdRecord> &records) {
    records.clear();
    FILE *in = fopen(filename.c_str(), "rb");
    if (!in) return false;
    char magic[8];
    uint32_t version = 0;
    if (fread(magic, sizeof(char), 8, in) != 8 ||
        fread(&version, sizeof(uint32_t), 1, in) != 1) {
        fclose(in);
        throw("Error: invalid sparse source LD file [" + filename + "].");
    }

    if (memcmp(magic, SOURCE_SPARSE_LD_MAGIC_V1, 8) == 0 && version == SOURCE_SPARSE_LD_VERSION_V1) {
        (void)readBinaryString(in, filename);
        uint64_t numRecords = 0;
        if (fread(&numRecords, sizeof(uint64_t), 1, in) != 1) {
            fclose(in);
            throw("Error: incomplete sparse source LD header in [" + filename + "].");
        }
        records.resize(numRecords);
        if (numRecords && fread(records.data(), sizeof(SparseSourceLdRecord), numRecords, in) != numRecords) {
            fclose(in);
            throw("Error: incomplete sparse source LD records in [" + filename + "].");
        }
        fclose(in);
        return true;
    }

    const bool sparseV2Q16 = memcmp(magic, SOURCE_SPARSE_LD_MAGIC_V2, 8) == 0 && version == SOURCE_SPARSE_LD_VERSION_V2;
    const bool sparseV3 = memcmp(magic, SOURCE_SPARSE_LD_MAGIC_V3, 8) == 0 && version == SOURCE_SPARSE_LD_VERSION_V3;
    if (!sparseV2Q16 && !sparseV3) {
        fclose(in);
        throw("Error: unsupported sparse source LD file version in [" + filename + "].");
    }

    (void)readBinaryString(in, filename);
    uint32_t numSnps = 0;
    uint64_t numRecords = 0;
    uint8_t indexWidth = 0;
    uint8_t corrWidth = 2;
    float corrScale = 0.0f;
    if (fread(&numSnps, sizeof(uint32_t), 1, in) != 1 ||
        fread(&numRecords, sizeof(uint64_t), 1, in) != 1 ||
        fread(&indexWidth, sizeof(uint8_t), 1, in) != 1) {
        fclose(in);
        throw("Error: incomplete sparse source LD header in [" + filename + "].");
    }
    if (sparseV3 && fread(&corrWidth, sizeof(uint8_t), 1, in) != 1) {
        fclose(in);
        throw("Error: incomplete sparse source LD v3 correlation-width header in [" + filename + "].");
    }
    if (fread(&corrScale, sizeof(float), 1, in) != 1) {
        fclose(in);
        throw("Error: incomplete sparse source LD correlation-scale header in [" + filename + "].");
    }
    if ((indexWidth != 2 && indexWidth != 4) || (corrWidth != 2 && corrWidth != 4)
        || (corrWidth == 2 && corrScale <= 0.0f)) {
        fclose(in);
        throw("Error: invalid sparse source LD header in [" + filename + "].");
    }
    vector<uint32_t> localToGlobal(numSnps);
    if (numSnps && fread(localToGlobal.data(), sizeof(uint32_t), numSnps, in) != numSnps) {
        fclose(in);
        throw("Error: incomplete sparse source LD v2 local-to-global index map in [" + filename + "].");
    }

    records.reserve(numRecords);
    for (uint64_t r = 0; r < numRecords; ++r) {
        uint32_t snp1 = 0;
        uint32_t snp2 = 0;
        if (indexWidth == 2) {
            uint16_t i16 = 0;
            uint16_t j16 = 0;
            if (fread(&i16, sizeof(uint16_t), 1, in) != 1 ||
                fread(&j16, sizeof(uint16_t), 1, in) != 1) {
                fclose(in);
                throw("Error: incomplete sparse source LD v2 record in [" + filename + "].");
            }
            snp1 = (uint32_t)i16;
            snp2 = (uint32_t)j16;
        } else {
            if (fread(&snp1, sizeof(uint32_t), 1, in) != 1 ||
                fread(&snp2, sizeof(uint32_t), 1, in) != 1) {
                fclose(in);
                throw("Error: incomplete sparse source LD v2 record in [" + filename + "].");
            }
        }
        float ldcor = 0.0f;
        if (corrWidth == 2) {
            int16_t q = 0;
            if (fread(&q, sizeof(int16_t), 1, in) != 1) {
                fclose(in);
                throw("Error: incomplete sparse source LD q16 correlation in [" + filename + "].");
            }
            ldcor = ((float)q) / corrScale;
        } else {
            if (fread(&ldcor, sizeof(float), 1, in) != 1) {
                fclose(in);
                throw("Error: incomplete sparse source LD float correlation in [" + filename + "].");
            }
        }
        if (snp1 >= numSnps || snp2 >= numSnps ||
            localToGlobal[snp1] == SOURCE_SPARSE_LD_MISSING_INDEX ||
            localToGlobal[snp2] == SOURCE_SPARSE_LD_MISSING_INDEX) {
            fclose(in);
            throw("Error: SNP index out of range in sparse source LD v2 file [" + filename + "].");
        }
        SparseSourceLdRecord rec;
        rec.snp1 = localToGlobal[snp1];
        rec.snp2 = localToGlobal[snp2];
        rec.ldcor = ldcor;
        records.push_back(rec);
    }
    fclose(in);
    return true;
}

void writePrecomputedMarginalizedContextBundle(
        const string &dirname, const string &sourceBlockID,
        const vector<LDBlockInfo*> &tileBlocks, const vector<VectorXf> &tileEigenvalues,
        const vector<MatrixXf> &tileQ, const string &sparseFilename,
        uint64_t &totalDirectedLinks, uint64_t &totalTileFlankSnps) {
    const unsigned numTiles = tileBlocks.size();
    if (!numTiles || tileEigenvalues.size() != numTiles || tileQ.size() != numTiles)
        throw("Error: incomplete tile data while preparing marginalized context for source block "
              + sourceBlockID + ".");

    unsigned sourceSize = 0;
    for (unsigned t = 0; t < numTiles; ++t) sourceSize += tileBlocks[t]->snpInfoVec.size();

    unordered_map<unsigned, pair<unsigned, unsigned> > tileLocalByGlobal;
    unordered_map<unsigned, unsigned> sourceColumnByGlobal;
    tileLocalByGlobal.reserve(sourceSize * 2 + 1);
    sourceColumnByGlobal.reserve(sourceSize * 2 + 1);
    unsigned sourceColumn = 0;
    for (unsigned t = 0; t < numTiles; ++t) {
        if (tileQ[t].cols() != (int)tileBlocks[t]->snpInfoVec.size()
            || tileQ[t].rows() != tileEigenvalues[t].size()) {
            throw("Error: inconsistent tile eigen data for " + tileBlocks[t]->ID + ".");
        }
        for (unsigned j = 0; j < tileBlocks[t]->snpInfoVec.size(); ++j) {
            SnpInfo *snp = tileBlocks[t]->snpInfoVec[j];
            if (!snp || snp->index < 0)
                throw("Error: invalid SNP index while preparing marginalized context for "
                      + tileBlocks[t]->ID + ".");
            tileLocalByGlobal[(unsigned)snp->index] = make_pair(t, j);
            sourceColumnByGlobal[(unsigned)snp->index] = sourceColumn++;
        }
    }

    vector<MatrixXf> transforms(numTiles);
    vector<MatrixXf> crossProjections(numTiles);
    vector<vector<unsigned char> > flankSeen(numTiles, vector<unsigned char>(sourceSize, 0));
    vector<uint64_t> linkCounts(numTiles, 0);
    for (unsigned t = 0; t < numTiles; ++t) {
        transforms[t] = tileQ[t];
        for (int r = 0; r < transforms[t].rows(); ++r) {
            if (!(tileEigenvalues[t][r] > 0.0f))
                throw("Error: non-positive tile eigenvalue while preparing marginalized context for "
                      + tileBlocks[t]->ID + ".");
            transforms[t].row(r) /= tileEigenvalues[t][r];
        }
        crossProjections[t].setZero(tileQ[t].rows(), sourceSize);
    }

    vector<SparseSourceLdRecord> records;
    if (!readSparseSourceLdBinary(sparseFilename, records))
        throw("Error: can not read sparse cross-tile LD file [" + sparseFilename
              + "] while preparing marginalized context.");
    for (unsigned r = 0; r < records.size(); ++r) {
        unordered_map<unsigned, pair<unsigned, unsigned> >::const_iterator a = tileLocalByGlobal.find(records[r].snp1);
        unordered_map<unsigned, pair<unsigned, unsigned> >::const_iterator b = tileLocalByGlobal.find(records[r].snp2);
        if (a == tileLocalByGlobal.end() || b == tileLocalByGlobal.end())
            throw("Error: sparse cross-tile LD references an unknown SNP in source block " + sourceBlockID + ".");
        if (a->second.first == b->second.first) continue;

        const unsigned tileA = a->second.first;
        const unsigned tileB = b->second.first;
        const unsigned colA = sourceColumnByGlobal[records[r].snp1];
        const unsigned colB = sourceColumnByGlobal[records[r].snp2];
        crossProjections[tileA].col(colB) += transforms[tileA].col(a->second.second) * records[r].ldcor;
        crossProjections[tileB].col(colA) += transforms[tileB].col(b->second.second) * records[r].ldcor;
        flankSeen[tileA][colB] = 1;
        flankSeen[tileB][colA] = 1;
        ++linkCounts[tileA];
        ++linkCounts[tileB];
    }
    vector<SparseSourceLdRecord>().swap(records);
    vector<MatrixXf>().swap(transforms);

    FILE *out = openMarginalizedContextBundleForWrite(dirname, sourceBlockID, numTiles);
    totalDirectedLinks = 0;
    totalTileFlankSnps = 0;
    for (unsigned t = 0; t < numTiles; ++t) {
        MatrixXf G(crossProjections[t].rows(), crossProjections[t].rows());
        G.noalias() = crossProjections[t] * crossProjections[t].transpose();
        crossProjections[t].resize(0, 0);

        SelfAdjointEigenSolver<MatrixXf> eig(G);
        if (eig.info() != Success) {
            fclose(out);
            throw("Error: context eigendecomposition failed for tile " + tileBlocks[t]->ID + ".");
        }
        VectorXf d = eig.eigenvalues();
        for (int j = 0; j < d.size(); ++j)
            if (d[j] < 0.0f) d[j] = 0.0f;
        const MatrixXf &V = eig.eigenvectors();

        MatrixXf invLambdaV = V;
        for (int r = 0; r < invLambdaV.rows(); ++r)
            invLambdaV.row(r) /= tileEigenvalues[t][r];
        MatrixXf rotatedWOperator = V.transpose() * invLambdaV;
        MatrixXf rotatedQ = V.transpose() * tileQ[t];

        uint32_t numFlankSnps = 0;
        for (unsigned j = 0; j < flankSeen[t].size(); ++j)
            if (flankSeen[t][j]) ++numFlankSnps;
        writeMarginalizedContextBundleRecord(out, tileBlocks[t]->ID, tileEigenvalues[t], d,
                                             rotatedWOperator, rotatedQ, linkCounts[t], numFlankSnps);
        totalDirectedLinks += linkCounts[t];
        totalTileFlankSnps += numFlankSnps;
    }
    fclose(out);
}

void appendSparseSourceRecordToCoreContexts(
    const uint32_t snp1,
    const uint32_t snp2,
    const float ldcor,
    const unordered_map<unsigned, pair<unsigned, unsigned> > &coreByGlobal,
    vector<CoreContextLdBlock> &coreContextLdBlocks) {

    unordered_map<unsigned, pair<unsigned, unsigned> >::const_iterator c1 = coreByGlobal.find(snp1);
    unordered_map<unsigned, pair<unsigned, unsigned> >::const_iterator c2 = coreByGlobal.find(snp2);
    if (c1 == coreByGlobal.end() && c2 == coreByGlobal.end()) return;
    if (c1 != coreByGlobal.end() && c2 != coreByGlobal.end()
        && c1->second.first == c2->second.first) {
        return;  // within-core LD is already represented by that tile's Q matrix
    }
    if (c1 != coreByGlobal.end()) {
        CoreContextLdBlock &context = coreContextLdBlocks[c1->second.first];
        context.coreLocalIndex.push_back(c1->second.second);
        context.flankGlobalIndex.push_back(snp2);
        context.ldCorr.push_back(ldcor);
    }
    if (c2 != coreByGlobal.end()) {
        CoreContextLdBlock &context = coreContextLdBlocks[c2->second.first];
        context.coreLocalIndex.push_back(c2->second.second);
        context.flankGlobalIndex.push_back(snp1);
        context.ldCorr.push_back(ldcor);
    }
}

bool appendSparseSourceLdBinaryToCoreContexts(
    const string &filename,
    const unordered_map<unsigned, pair<unsigned, unsigned> > &coreByGlobal,
    vector<CoreContextLdBlock> &coreContextLdBlocks,
    const vector<uint32_t> *sourceCurrentByLocal = NULL) {

    FILE *in = fopen(filename.c_str(), "rb");
    if (!in) return false;
    char magic[8];
    uint32_t version = 0;
    if (fread(magic, sizeof(char), 8, in) != 8 ||
        fread(&version, sizeof(uint32_t), 1, in) != 1) {
        fclose(in);
        throw("Error: invalid sparse source LD file [" + filename + "].");
    }

    if (memcmp(magic, SOURCE_SPARSE_LD_MAGIC_V1, 8) == 0 && version == SOURCE_SPARSE_LD_VERSION_V1) {
        (void)readBinaryString(in, filename);
        uint64_t numRecords = 0;
        if (fread(&numRecords, sizeof(uint64_t), 1, in) != 1) {
            fclose(in);
            throw("Error: incomplete sparse source LD header in [" + filename + "].");
        }
        for (uint64_t r = 0; r < numRecords; ++r) {
            SparseSourceLdRecord rec;
            if (fread(&rec, sizeof(SparseSourceLdRecord), 1, in) != 1) {
                fclose(in);
                throw("Error: incomplete sparse source LD records in [" + filename + "].");
            }
            appendSparseSourceRecordToCoreContexts(rec.snp1, rec.snp2, rec.ldcor, coreByGlobal, coreContextLdBlocks);
        }
        fclose(in);
        return true;
    }

    const bool sparseV2Q16 = memcmp(magic, SOURCE_SPARSE_LD_MAGIC_V2, 8) == 0 && version == SOURCE_SPARSE_LD_VERSION_V2;
    const bool sparseV3 = memcmp(magic, SOURCE_SPARSE_LD_MAGIC_V3, 8) == 0 && version == SOURCE_SPARSE_LD_VERSION_V3;
    if (!sparseV2Q16 && !sparseV3) {
        fclose(in);
        throw("Error: unsupported sparse source LD file version in [" + filename + "].");
    }

    (void)readBinaryString(in, filename);
    uint32_t numSnps = 0;
    uint64_t numRecords = 0;
    uint8_t indexWidth = 0;
    uint8_t corrWidth = 2;
    float corrScale = 0.0f;
    if (fread(&numSnps, sizeof(uint32_t), 1, in) != 1 ||
        fread(&numRecords, sizeof(uint64_t), 1, in) != 1 ||
        fread(&indexWidth, sizeof(uint8_t), 1, in) != 1) {
        fclose(in);
        throw("Error: incomplete sparse source LD header in [" + filename + "].");
    }
    if (sparseV3 && fread(&corrWidth, sizeof(uint8_t), 1, in) != 1) {
        fclose(in);
        throw("Error: incomplete sparse source LD v3 correlation-width header in [" + filename + "].");
    }
    if (fread(&corrScale, sizeof(float), 1, in) != 1) {
        fclose(in);
        throw("Error: incomplete sparse source LD correlation-scale header in [" + filename + "].");
    }
    if ((indexWidth != 2 && indexWidth != 4) || (corrWidth != 2 && corrWidth != 4)
        || (corrWidth == 2 && corrScale <= 0.0f)) {
        fclose(in);
        throw("Error: invalid sparse source LD header in [" + filename + "].");
    }
    vector<uint32_t> localToGlobal(numSnps);
    if (numSnps && fread(localToGlobal.data(), sizeof(uint32_t), numSnps, in) != numSnps) {
        fclose(in);
        throw("Error: incomplete sparse source LD v2 local-to-global index map in [" + filename + "].");
    }
    unordered_map<uint32_t, uint32_t> sidecarIndexToCurrent;
    if (sourceCurrentByLocal && sourceCurrentByLocal->size() == numSnps) {
        sidecarIndexToCurrent.reserve(numSnps * 2 + 1);
        for (uint32_t j = 0; j < numSnps; ++j) {
            if (localToGlobal[j] == SOURCE_SPARSE_LD_MISSING_INDEX ||
                (*sourceCurrentByLocal)[j] == SOURCE_SPARSE_LD_MISSING_INDEX) continue;
            sidecarIndexToCurrent[localToGlobal[j]] = (*sourceCurrentByLocal)[j];
        }
    }

    for (uint64_t r = 0; r < numRecords; ++r) {
        uint32_t snp1 = 0;
        uint32_t snp2 = 0;
        if (indexWidth == 2) {
            uint16_t i16 = 0;
            uint16_t j16 = 0;
            if (fread(&i16, sizeof(uint16_t), 1, in) != 1 ||
                fread(&j16, sizeof(uint16_t), 1, in) != 1) {
                fclose(in);
                throw("Error: incomplete sparse source LD v2 record in [" + filename + "].");
            }
            snp1 = (uint32_t)i16;
            snp2 = (uint32_t)j16;
        } else {
            if (fread(&snp1, sizeof(uint32_t), 1, in) != 1 ||
                fread(&snp2, sizeof(uint32_t), 1, in) != 1) {
                fclose(in);
                throw("Error: incomplete sparse source LD v2 record in [" + filename + "].");
            }
        }
        float ldcor = 0.0f;
        if (corrWidth == 2) {
            int16_t q = 0;
            if (fread(&q, sizeof(int16_t), 1, in) != 1) {
                fclose(in);
                throw("Error: incomplete sparse source LD q16 correlation in [" + filename + "].");
            }
            ldcor = ((float)q) / corrScale;
        } else {
            if (fread(&ldcor, sizeof(float), 1, in) != 1) {
                fclose(in);
                throw("Error: incomplete sparse source LD float correlation in [" + filename + "].");
            }
        }
        if (snp1 >= numSnps || snp2 >= numSnps ||
            localToGlobal[snp1] == SOURCE_SPARSE_LD_MISSING_INDEX ||
            localToGlobal[snp2] == SOURCE_SPARSE_LD_MISSING_INDEX) {
            fclose(in);
            throw("Error: SNP index out of range in sparse source LD v2 file [" + filename + "].");
        }
        uint32_t current1 = localToGlobal[snp1];
        uint32_t current2 = localToGlobal[snp2];
        if (!sidecarIndexToCurrent.empty()) {
            unordered_map<uint32_t, uint32_t>::const_iterator map1 = sidecarIndexToCurrent.find(current1);
            unordered_map<uint32_t, uint32_t>::const_iterator map2 = sidecarIndexToCurrent.find(current2);
            if (map1 == sidecarIndexToCurrent.end() || map2 == sidecarIndexToCurrent.end()) continue;
            current1 = map1->second;
            current2 = map2->second;
        }
        appendSparseSourceRecordToCoreContexts(current1, current2, ldcor, coreByGlobal, coreContextLdBlocks);
    }
    fclose(in);
    return true;
}

void writePwldBinHeader(FILE *out, uint64_t numSnps, uint64_t numRecords, const vector<string> &snpIDs) {
    if (snpIDs.size() != numSnps) {
        throw("Error: internal mismatch while writing binary pairwise LD SNP IDs.");
    }
    fwrite(PWLD_BIN_MAGIC, sizeof(char), 8, out);
    fwrite(&PWLD_BIN_VERSION, sizeof(uint32_t), 1, out);
    fwrite(&numSnps, sizeof(uint64_t), 1, out);
    fwrite(&numRecords, sizeof(uint64_t), 1, out);
    for (uint64_t i = 0; i < numSnps; ++i) {
        const string &id = snpIDs[i];
        uint32_t len = (uint32_t)id.size();
        fwrite(&len, sizeof(uint32_t), 1, out);
        if (len && fwrite(id.data(), sizeof(char), len, out) != len) {
            throw("Error: failed writing binary pairwise LD SNP ID.");
        }
    }
}

bool readPwldBinHeaderLocal(FILE *in, const string &filename, uint64_t &numSnps, uint64_t &numRecords, vector<string> &snpIDs) {
    char magic[8];
    uint32_t version = 0;
    if (fread(magic, sizeof(char), 8, in) != 8) return false;
    if (memcmp(magic, PWLD_BIN_MAGIC, 8) != 0) {
        throw("Error: invalid binary pairwise LD file header in [" + filename + "].");
    }
    if (fread(&version, sizeof(uint32_t), 1, in) != 1 ||
        fread(&numSnps, sizeof(uint64_t), 1, in) != 1 ||
        fread(&numRecords, sizeof(uint64_t), 1, in) != 1) {
        throw("Error: incomplete binary pairwise LD file header in [" + filename + "].");
    }
    if (version != PWLD_BIN_VERSION) {
        throw("Error: unsupported binary pairwise LD file version in [" + filename + "]. Please regenerate .pwld.bin so SNP IDs are embedded.");
    }
    snpIDs.clear();
    snpIDs.reserve(numSnps);
    for (uint64_t i = 0; i < numSnps; ++i) {
        uint32_t len = 0;
        if (fread(&len, sizeof(uint32_t), 1, in) != 1) {
            throw("Error: incomplete SNP ID table in binary pairwise LD file [" + filename + "].");
        }
        string id(len, '\0');
        if (len && fread(&id[0], sizeof(char), len, in) != len) {
            throw("Error: incomplete SNP ID in binary pairwise LD file [" + filename + "].");
        }
        snpIDs.push_back(id);
    }
    return true;
}

bool readGzLineLocal(gzFile in, string &line) {
    line.clear();
    char buffer[65536];
    bool got = false;
    while (gzgets(in, buffer, sizeof(buffer))) {
        got = true;
        line += buffer;
        if (!line.empty() && line[line.size() - 1] == '\n') break;
    }
    if (!got) return false;
    while (!line.empty() && (line[line.size() - 1] == '\n' || line[line.size() - 1] == '\r')) {
        line.erase(line.size() - 1);
    }
    return true;
}

bool pairwiseLdHasRecordsLocal(const string &filename) {
    if (!fileExistsAndNonEmptyLocal(filename)) return false;

    if (hasSuffixLocal(filename, ".bin")) {
        FILE *in = fopen(filename.c_str(), "rb");
        if (!in) return false;
        char magic[8];
        uint32_t version = 0;
        uint64_t numSnps = 0;
        uint64_t numRecords = 0;
        const bool valid = fread(magic, sizeof(char), 8, in) == 8
            && fread(&version, sizeof(uint32_t), 1, in) == 1
            && fread(&numSnps, sizeof(uint64_t), 1, in) == 1
            && fread(&numRecords, sizeof(uint64_t), 1, in) == 1
            && memcmp(magic, PWLD_BIN_MAGIC, 8) == 0
            && version == PWLD_BIN_VERSION;
        fclose(in);
        return valid && numRecords > 0;
    }

    string line;
    if (hasSuffixLocal(filename, ".gz")) {
        gzFile in = gzopen(filename.c_str(), "rb");
        if (!in) return false;
        readGzLineLocal(in, line);  // header
        while (readGzLineLocal(in, line)) {
            if (!line.empty()) {
                gzclose(in);
                return true;
            }
        }
        gzclose(in);
        return false;
    }

    ifstream in(filename.c_str());
    getline(in, line);  // header
    while (getline(in, line)) {
        if (!line.empty()) return true;
    }
    return false;
}

string findPairwiseLdWithRecordsForBaseLocal(const string &base) {
    const string candidates[3] = {base + ".bin", base, base + ".gz"};
    for (unsigned i = 0; i < 3; ++i) {
        if (pairwiseLdHasRecordsLocal(candidates[i])) return candidates[i];
    }
    return "";
}

string absolutePathLocal(const string &filename) {
    char resolved[PATH_MAX];
    if (!realpath(filename.c_str(), resolved))
        throw("Error: can not resolve pairwise LD file path [" + filename + "].");
    return string(resolved);
}

void copyFileBinaryLocal(const string &source, const string &destination) {
    ifstream in(source.c_str(), ios::binary);
    if (!in) throw("Error: can not open pairwise LD source [" + source + "] to copy.");
    const string temporary = destination + ".tmp";
    ofstream out(temporary.c_str(), ios::binary | ios::trunc);
    if (!out) throw("Error: can not open pairwise LD destination [" + temporary + "] to write.");
    out << in.rdbuf();
    out.close();
    if (in.bad() || !out)
        throw("Error: failed while copying pairwise LD file [" + source + "] to [" + destination + "].");
    if (rename(temporary.c_str(), destination.c_str()) != 0) {
        remove(temporary.c_str());
        throw("Error: can not install copied pairwise LD file [" + destination + "].");
    }
}

string copyPairwiseLdToBaseLocal(const string &source, const string &destinationBase) {
    string destination = destinationBase;
    if (hasSuffixLocal(source, ".bin")) destination += ".bin";
    else if (hasSuffixLocal(source, ".gz")) destination += ".gz";
    copyFileBinaryLocal(source, destination);
    const string variants[3] = {destinationBase, destinationBase + ".gz", destinationBase + ".bin"};
    for (unsigned i = 0; i < 3; ++i) {
        if (variants[i] != destination) removeFileIfExistsLocal(variants[i]);
    }
    return destination;
}

void writePwldSourceMarkerLocal(const string &dirname, const string &blockID, const string &source) {
    const string marker = dirname + "/block" + blockID + ".pwld.source";
    ofstream out(marker.c_str());
    if (!out) throw("Error: can not write pairwise LD source marker [" + marker + "].");
    out << absolutePathLocal(source) << '\n';
}

string installRecordedPairwiseLdSourceLocal(const string &dirname) {
    DIR *dirp = opendir(dirname.c_str());
    if (!dirp) throw("Error: can not open folder [" + dirname + "] to find pairwise LD source markers.");
    vector<string> markers;
    dirent *dp;
    while ((dp = readdir(dirp)) != NULL) {
        const string fname = dp->d_name;
        if (fname.find("block") == 0 && hasSuffixLocal(fname, ".pwld.source"))
            markers.push_back(dirname + "/" + fname);
    }
    closedir(dirp);
    if (markers.empty()) return "";

    string source;
    for (unsigned i = 0; i < markers.size(); ++i) {
        ifstream in(markers[i].c_str());
        string candidate;
        getline(in, candidate);
        if (!candidate.empty() && candidate[candidate.size() - 1] == '\r') candidate.erase(candidate.size() - 1);
        if (candidate.empty() || !pairwiseLdHasRecordsLocal(candidate))
            throw("Error: pairwise LD source recorded in [" + markers[i] + "] is missing or has no records.");
        candidate = absolutePathLocal(candidate);
        if (source.empty()) source = candidate;
        else if (source != candidate)
            throw("Error: tiled block jobs recorded different combined pairwise LD sources; cannot choose one safely.");
    }

    const string existing = findPairwiseLdWithRecordsForBaseLocal(dirname + "/rsq0.5.pwld");
    if (!existing.empty()) return existing;
    const string installed = copyPairwiseLdToBaseLocal(source, dirname + "/rsq0.5.pwld");
    cout << "Copied source pairwise LD file [" << source << "] to [" << installed << "]." << endl;
    return installed;
}

struct EigenNbThreadsGuard {
    int saved;
    EigenNbThreadsGuard() : saved(Eigen::nbThreads()) {
#ifdef _OPENMP
        if (omp_get_max_threads() > 1) {
            omp_set_max_active_levels(1);
            Eigen::setNbThreads(1);
        }
#endif
    }
    ~EigenNbThreadsGuard() { Eigen::setNbThreads(saved); }
};

/** Anchor SNPs for --impute-summary: must have finite beta and positive finite SE (z = b/SE). */
inline bool gwasUsableAnchorSnp(const SnpInfo *snp) {
    const bool finite_b = std::isfinite(snp->gwas_b);
    const bool finite_se = std::isfinite(snp->gwas_se);
    const bool pos_se = snp->gwas_se > 0.0;
    const bool ok = snp->included && finite_b && finite_se && pos_se;
    // cout << snp->ID << " included=" << snp->included
    //      << " gwas_b=" << snp->gwas_b << " finite_b=" << finite_b
    //      << " gwas_se=" << snp->gwas_se << " finite_se=" << finite_se
    //      << " se>0=" << pos_se << " usable=" << ok << endl;
    return ok;
}

[[noreturn]] inline void throwEigenReadErr(const string &msg) {
    cerr << msg << endl;
    throw std::runtime_error(msg);
}

inline string eigenBlockFilename(const string &dirname, const string &blockID, const int quantizedBits, const bool q8Entropy) {
    string suffix;
    if (quantizedBits == 0) suffix = ".eigen.bin";
    else if (quantizedBits == 8 && q8Entropy) suffix = ".eigen.q8e.bin";
    else if (quantizedBits == 8) suffix = ".eigen.q8.bin";
    else if (quantizedBits == 16) suffix = ".eigen.q16.bin";
    else throw std::runtime_error("Error: unsupported quantized eigen bit width " + to_string(quantizedBits) + ". Use 8 or 16.");
    return dirname + "/block" + blockID + suffix;
}

bool sourceBlockIDFromTileID(const string &blockID, string &sourceBlockID) {
    const string marker = ".tile";
    size_t pos = blockID.rfind(marker);
    if (pos == string::npos || pos == 0) return false;
    sourceBlockID = blockID.substr(0, pos);
    return true;
}

bool openMarginalizedContextPayloadForRead(const string &dirname, const string &blockID,
                                            FILE *&fp, string &infile) {
    fp = NULL;
    string sourceBlockID;
    if (!sourceBlockIDFromTileID(blockID, sourceBlockID)) return false;
    infile = marginalizedContextBundleFilename(dirname, sourceBlockID);
    fp = fopen(infile.c_str(), "rb");
    if (!fp) return false;

    char magic[8];
    uint32_t version = 0, nrec = 0;
    if (fread(magic, sizeof(char), 8, fp) != 8
        || memcmp(magic, MARGINAL_CONTEXT_BUNDLE_MAGIC, 8) != 0
        || fread(&version, sizeof(uint32_t), 1, fp) != 1
        || fread(&nrec, sizeof(uint32_t), 1, fp) != 1) {
        fclose(fp);
        fp = NULL;
        throw("Error: invalid marginalized-context bundle [" + infile + "].");
    }
    if (version != MARGINAL_CONTEXT_BUNDLE_VERSION) {
        fclose(fp);
        fp = NULL;
        throw("Error: unsupported marginalized-context bundle format [" + infile + "].");
    }

    for (uint32_t r = 0; r < nrec; ++r) {
        const string recordID = readBinaryString(fp, infile);
        uint64_t payloadBytes = 0;
        if (fread(&payloadBytes, sizeof(uint64_t), 1, fp) != 1) {
            fclose(fp);
            fp = NULL;
            throw("Error: incomplete marginalized-context bundle index [" + infile + "].");
        }
        if (recordID == blockID) {
            infile += "::" + blockID;
            return true;
        }
        if (fseek(fp, (long)payloadBytes, SEEK_CUR) != 0) {
            fclose(fp);
            fp = NULL;
            throw("Error: invalid marginalized-context bundle payload in [" + infile + "].");
        }
    }

    fclose(fp);
    fp = NULL;
    return false;
}

void readMarginalizedContextPayloadOrThrow(FILE *fp, const string &infile, const string &blockID,
                                            const int expectedNumSnps, VectorXf &lambda,
                                            VectorXf &contextEigenvalues, MatrixXf &rotatedWOperator,
                                            MatrixXf &rotatedQ, uint64_t &numLinks,
                                            unsigned &numFlankSnps) {
    int32_t numSnps = 0, rank = 0;
    uint32_t flankCount = 0;
    if (fread(&numSnps, sizeof(int32_t), 1, fp) != 1
        || fread(&rank, sizeof(int32_t), 1, fp) != 1
        || fread(&numLinks, sizeof(uint64_t), 1, fp) != 1
        || fread(&flankCount, sizeof(uint32_t), 1, fp) != 1) {
        throwEigenReadErr("Error: incomplete marginalized-context payload in " + infile + ".");
    }
    if (numSnps != expectedNumSnps || rank <= 0) {
        throwEigenReadErr("Error: inconsistent marginalized-context dimensions for tile " + blockID
                          + " in " + infile + ".");
    }

    lambda.resize(rank);
    contextEigenvalues.resize(rank);
    rotatedWOperator.resize(rank, rank);
    rotatedQ.resize(rank, numSnps);
    if (fread(lambda.data(), sizeof(float), rank, fp) != (size_t)rank
        || fread(contextEigenvalues.data(), sizeof(float), rank, fp) != (size_t)rank
        || fread(rotatedWOperator.data(), sizeof(float), (uint64_t)rank * rank, fp)
             != (size_t)((uint64_t)rank * rank)
        || fread(rotatedQ.data(), sizeof(float), (uint64_t)rank * numSnps, fp)
             != (size_t)((uint64_t)rank * numSnps)) {
        throwEigenReadErr("Error: incomplete marginalized-context matrices for tile " + blockID
                          + " in " + infile + ".");
    }
    numFlankSnps = (unsigned)flankCount;
}

bool openEigenBlockPayloadForRead(const string &dirname, const string &blockID,
                                  const int quantizedBits, const bool q8Entropy,
                                  FILE *&fp, string &infile) {
    fp = NULL;
    infile = eigenBlockFilename(dirname, blockID, quantizedBits, q8Entropy);
    fp = fopen(infile.c_str(), "rb");
    if (fp) return true;

    if (q8Entropy || (quantizedBits != 0 && quantizedBits != 16)) return false;
    string sourceBlockID;
    if (!sourceBlockIDFromTileID(blockID, sourceBlockID)) return false;

    infile = tiledEigenBundleFilename(dirname, sourceBlockID, quantizedBits);
    fp = fopen(infile.c_str(), "rb");
    if (!fp) return false;

    char magic[8];
    uint32_t version = 0, bits = 0, nrec = 0;
    if (fread(magic, sizeof(char), 8, fp) != 8 ||
        memcmp(magic, TILED_EIGEN_BUNDLE_MAGIC, 8) != 0 ||
        fread(&version, sizeof(uint32_t), 1, fp) != 1 ||
        fread(&bits, sizeof(uint32_t), 1, fp) != 1 ||
        fread(&nrec, sizeof(uint32_t), 1, fp) != 1) {
        fclose(fp);
        fp = NULL;
        throw("Error: invalid tiled eigen bundle [" + infile + "].");
    }
    if (version != TILED_EIGEN_BUNDLE_VERSION || (int)bits != quantizedBits) {
        fclose(fp);
        fp = NULL;
        throw("Error: unsupported tiled eigen bundle format [" + infile + "].");
    }

    for (uint32_t r = 0; r < nrec; ++r) {
        string recordID = readBinaryString(fp, infile);
        uint64_t payloadBytes = 0;
        if (fread(&payloadBytes, sizeof(uint64_t), 1, fp) != 1) {
            fclose(fp);
            fp = NULL;
            throw("Error: incomplete tiled eigen bundle index [" + infile + "].");
        }
        const long payloadStart = ftell(fp);
        if (recordID == blockID) {
            infile += "::" + blockID;
            return true;
        }
        if (fseek(fp, (long)payloadBytes, SEEK_CUR) != 0) {
            fclose(fp);
            fp = NULL;
            throw("Error: invalid tiled eigen bundle payload size in [" + infile + "].");
        }
        if (ftell(fp) < payloadStart) {
            fclose(fp);
            fp = NULL;
            throw("Error: invalid tiled eigen bundle offset in [" + infile + "].");
        }
    }

    fclose(fp);
    fp = NULL;
    return false;
}

void readEigenHeaderOrThrow(FILE *fp, const string &infile, const string &blockID, const int expectedM,
                            int32_t &cur_m, int32_t &cur_k, float &sumPosEigVal, float &oldEigenCutoff, VectorXf &lambda) {
    if (fread(&cur_m, sizeof(int32_t), 1, fp) != 1) throwEigenReadErr("Read " + infile + " error (m)");
    if (cur_m != expectedM) throwEigenReadErr("In LD block " + blockID + ", inconsistent marker number to marker information in " + infile);
    if (fread(&cur_k, sizeof(int32_t), 1, fp) != 1) throwEigenReadErr("In LD block " + blockID + ", error about number of eigenvalues in  " + infile);
    if (fread(&sumPosEigVal, sizeof(float), 1, fp) != 1) throwEigenReadErr("In LD block " + blockID + ", error about the sum of positive eigenvalues in " + infile);
    if (fread(&oldEigenCutoff, sizeof(float), 1, fp) != 1) throwEigenReadErr("In LD block " + blockID + ", error about eigen cutoff used in " + infile);
    lambda.resize(cur_k);
    if (fread(lambda.data(), sizeof(float), cur_k, fp) != (size_t)cur_k) throwEigenReadErr("In LD block " + blockID + ",size error about eigenvalues in " + infile);
}

void readEigenUOrThrow(FILE *fp, const string &infile, const string &blockID, const int32_t cur_m, const int32_t cur_k,
                       const int quantizedBits, const bool q8Entropy, MatrixXf &U) {
    U.resize(cur_m, cur_k);
    const uint64_t nElements = (uint64_t)cur_m * (uint64_t)cur_k;
    if (quantizedBits == 0) {
        size_t nReadU = fread(U.data(), sizeof(float), nElements, fp);
        if (nReadU != nElements) {
            cerr << "Expected " << nElements << " eigen matrix floats; got " << nReadU
                 << " (rows " << cur_m << ", cols " << cur_k << ")." << endl;
            throwEigenReadErr("In LD block " + blockID + ",size error about eigenvectors in " + infile);
        }
        return;
    }

    VectorXf scales(cur_k);
    if (fread(scales.data(), sizeof(float), cur_k, fp) != (size_t)cur_k) {
        throwEigenReadErr("In LD block " + blockID + ",size error about quantized eigenvector scales in " + infile);
    }
    float invBound = 0.0f;
    if (quantizedBits == 8) invBound = 1.0f / 127.0f;
    else if (quantizedBits == 16) invBound = 1.0f / 32767.0f;
    else throwEigenReadErr("Error: unsupported quantized eigen bit width " + to_string(quantizedBits) + ".");

    if (quantizedBits == 8) {
        vector<int8_t> raw(nElements);
        if (q8Entropy) {
            uint64_t uncompressedSize = 0, compressedSize = 0;
            if (fread(&uncompressedSize, sizeof(uint64_t), 1, fp) != 1)
                throwEigenReadErr("In LD block " + blockID + ", read error (zlib uncompressed size) in " + infile);
            if (fread(&compressedSize, sizeof(uint64_t), 1, fp) != 1)
                throwEigenReadErr("In LD block " + blockID + ", read error (zlib compressed size) in " + infile);
            if (uncompressedSize != nElements)
                throwEigenReadErr("In LD block " + blockID + ", zlib uncompressed size mismatch in " + infile);
            vector<Bytef> compressed(compressedSize);
            if (fread(compressed.data(), 1, compressedSize, fp) != compressedSize)
                throwEigenReadErr("In LD block " + blockID + ", size error about zlib payload in " + infile);
            uLongf destLen = (uLongf)uncompressedSize;
            int zrc = uncompress(reinterpret_cast<Bytef*>(raw.data()), &destLen, compressed.data(), (uLong)compressedSize);
            if (zrc != Z_OK || destLen != (uLongf)uncompressedSize)
                throwEigenReadErr("In LD block " + blockID + ", zlib uncompress failed for " + infile);
        } else {
            if (fread(raw.data(), sizeof(int8_t), nElements, fp) != nElements)
                throwEigenReadErr("In LD block " + blockID + ",size error about q8 eigenvectors in " + infile);
        }
        for (int col = 0; col < cur_k; ++col) {
            const float scale = scales[col] * invBound;
            for (int row = 0; row < cur_m; ++row)
                U(row, col) = scale * (float)raw[(uint64_t)col * (uint64_t)cur_m + (uint64_t)row];
        }
    } else {
        vector<int16_t> colRaw(cur_m);
        for (int col = 0; col < cur_k; ++col) {
            if (fread(colRaw.data(), sizeof(int16_t), cur_m, fp) != (size_t)cur_m)
                throwEigenReadErr("In LD block " + blockID + ",size error about q16 eigenvectors in " + infile);
            const float scale = scales[col] * invBound;
            for (int row = 0; row < cur_m; ++row)
                U(row, col) = scale * (float)colRaw[row];
        }
    }
}

void readQuantizedEigenUOrThrow(FILE *fp, const string &infile, const string &blockID, const int32_t cur_m, const int32_t cur_k,
                                const int quantizedBits, const bool q8Entropy, const VectorXf &lambda, QuantizedEigenBlock &qb) {
    if (quantizedBits == 0) throwEigenReadErr("Internal error: readQuantizedEigenUOrThrow called for unquantized eigen file.");
    qb.bits = quantizedBits;
    qb.m = cur_m;
    qb.k = cur_k;
    qb.sqrtLambda = lambda.array().sqrt();

    VectorXf scales(cur_k);
    if (fread(scales.data(), sizeof(float), cur_k, fp) != (size_t)cur_k) {
        throwEigenReadErr("In LD block " + blockID + ",size error about quantized eigenvector scales in " + infile);
    }

    float invBound = 0.0f;
    if (quantizedBits == 8) invBound = 1.0f / 127.0f;
    else if (quantizedBits == 16) invBound = 1.0f / 32767.0f;
    else throwEigenReadErr("Error: unsupported quantized eigen bit width " + to_string(quantizedBits) + ".");
    qb.dequantScale = scales * invBound;
    qb.qScale = qb.sqrtLambda.array() * qb.dequantScale.array();

    const uint64_t nElements = (uint64_t)cur_m * (uint64_t)cur_k;
    qb.q8.clear();
    qb.q16.clear();

    if (quantizedBits == 8) {
        vector<int8_t> raw(nElements);
        if (q8Entropy) {
            uint64_t uncompressedSize = 0, compressedSize = 0;
            if (fread(&uncompressedSize, sizeof(uint64_t), 1, fp) != 1)
                throwEigenReadErr("In LD block " + blockID + ", read error (zlib uncompressed size) in " + infile);
            if (fread(&compressedSize, sizeof(uint64_t), 1, fp) != 1)
                throwEigenReadErr("In LD block " + blockID + ", read error (zlib compressed size) in " + infile);
            if (uncompressedSize != nElements)
                throwEigenReadErr("In LD block " + blockID + ", zlib uncompressed size mismatch in " + infile);
            vector<Bytef> compressed(compressedSize);
            if (fread(compressed.data(), 1, compressedSize, fp) != compressedSize)
                throwEigenReadErr("In LD block " + blockID + ", size error about zlib payload in " + infile);
            uLongf destLen = (uLongf)uncompressedSize;
            int zrc = uncompress(reinterpret_cast<Bytef*>(raw.data()), &destLen, compressed.data(), (uLong)compressedSize);
            if (zrc != Z_OK || destLen != (uLongf)uncompressedSize)
                throwEigenReadErr("In LD block " + blockID + ", zlib uncompress failed for " + infile);
        } else {
            if (fread(raw.data(), sizeof(int8_t), nElements, fp) != nElements)
                throwEigenReadErr("In LD block " + blockID + ",size error about q8 eigenvectors in " + infile);
        }
        qb.q8.resize(nElements);
        for (int row = 0; row < cur_m; ++row)
            for (int col = 0; col < cur_k; ++col)
                qb.q8[(size_t)row * (size_t)cur_k + (size_t)col] = raw[(size_t)col * (size_t)cur_m + (size_t)row];
    } else {
        qb.q16.resize(nElements);
        for (int col = 0; col < cur_k; ++col) {
            vector<int16_t> colRaw(cur_m);
            if (fread(colRaw.data(), sizeof(int16_t), cur_m, fp) != (size_t)cur_m)
                throwEigenReadErr("In LD block " + blockID + ",size error about q16 eigenvectors in " + infile);
            for (int row = 0; row < cur_m; ++row)
                qb.q16[(size_t)row * (size_t)cur_k + (size_t)col] = colRaw[row];
        }
    }
}

int eigenTailStartIndex(const float sumPosEigVal, const float eigenCutoff, const VectorXf &lambda) {
    int revIdx = (int)lambda.size();
    VectorXf cumsumNonNeg(revIdx);
    cumsumNonNeg.setZero();
    revIdx = revIdx - 1;
    cumsumNonNeg(revIdx) = lambda(revIdx);
    revIdx = revIdx - 1;

    while (revIdx >= 0 && lambda(revIdx) > 1e-10f) {
        cumsumNonNeg(revIdx) = lambda(revIdx) + cumsumNonNeg(revIdx + 1);
        revIdx = revIdx - 1;
        if (revIdx < 0) break;
    }
    cumsumNonNeg = cumsumNonNeg / sumPosEigVal;

    bool haveValue = false;
    for (revIdx = revIdx + 1; revIdx < lambda.size(); revIdx++) {
        if (eigenCutoff >= cumsumNonNeg(revIdx)) {
            revIdx = revIdx - 1;
            haveValue = true;
            break;
        }
    }
    if (!haveValue) revIdx = lambda.size() - 1;
    if (revIdx < 0) revIdx = 0;
    return revIdx;
}

void keepQuantizedEigenTail(QuantizedEigenBlock &qb, const int startCol) {
    if (startCol <= 0) return;
    const int oldK = qb.k;
    const int newK = oldK - startCol;
    if (newK <= 0) return;

    qb.sqrtLambda = qb.sqrtLambda.tail(newK).eval();
    qb.dequantScale = qb.dequantScale.tail(newK).eval();
    qb.qScale = qb.qScale.tail(newK).eval();

    if (qb.bits == 8) {
        vector<int8_t> kept((size_t)qb.m * (size_t)newK);
        for (int row = 0; row < qb.m; ++row) {
            const int8_t *src = &qb.q8[(size_t)row * (size_t)oldK + (size_t)startCol];
            int8_t *dst = &kept[(size_t)row * (size_t)newK];
            std::copy(src, src + newK, dst);
        }
        qb.q8.swap(kept);
    } else if (qb.bits == 16) {
        vector<int16_t> kept((size_t)qb.m * (size_t)newK);
        for (int row = 0; row < qb.m; ++row) {
            const int16_t *src = &qb.q16[(size_t)row * (size_t)oldK + (size_t)startCol];
            int16_t *dst = &kept[(size_t)row * (size_t)newK];
            std::copy(src, src + newK, dst);
        }
        qb.q16.swap(kept);
    }
    qb.k = newK;
}

void quantizeEigenMatrixToBlock(const MatrixXf &U, const VectorXf &lambda, const int quantizedBits, QuantizedEigenBlock &qb) {
    const int m = U.rows();
    const int k = U.cols();
    qb = QuantizedEigenBlock();
    qb.bits = quantizedBits;
    qb.m = m;
    qb.k = k;
    qb.sqrtLambda = lambda.array().sqrt();
    qb.dequantScale.resize(k);

    int bound = 0;
    if (quantizedBits == 8) bound = 127;
    else if (quantizedBits == 16) bound = 32767;
    else throwEigenReadErr("Error: unsupported quantized eigen bit width " + to_string(quantizedBits) + ".");

    VectorXf maxAbs(k);
    for (int col = 0; col < k; ++col) {
        maxAbs[col] = U.col(col).cwiseAbs().maxCoeff();
        if (maxAbs[col] <= 0.0f) maxAbs[col] = 1.0f;
        qb.dequantScale[col] = maxAbs[col] / (float)bound;
    }
    qb.qScale = qb.sqrtLambda.array() * qb.dequantScale.array();

    if (quantizedBits == 8) {
        qb.q8.resize((size_t)m * (size_t)k);
        for (int col = 0; col < k; ++col) {
            const float invScale = (float)bound / maxAbs[col];
            for (int row = 0; row < m; ++row) {
                int q = (int)lrintf(U(row, col) * invScale);
                if (q > bound) q = bound;
                if (q < -bound) q = -bound;
                qb.q8[(size_t)row * (size_t)k + (size_t)col] = (int8_t)q;
            }
        }
    } else {
        qb.q16.resize((size_t)m * (size_t)k);
        for (int col = 0; col < k; ++col) {
            const float invScale = (float)bound / maxAbs[col];
            for (int row = 0; row < m; ++row) {
                int q = (int)lrintf(U(row, col) * invScale);
                if (q > bound) q = bound;
                if (q < -bound) q = -bound;
                qb.q16[(size_t)row * (size_t)k + (size_t)col] = (int16_t)q;
            }
        }
    }
}

VectorXf makeWcorrFromQuantizedEigen(const QuantizedEigenBlock &qb, const VectorXf &gwasEffects) {
    if (gwasEffects.size() != qb.m) {
        throw("Error: quantized eigen/GWAS dimension mismatch: eigen block has "
              + to_string(qb.m) + " SNPs but the GWAS-effect vector has "
              + to_string((long long)gwasEffects.size()) + ".");
    }
    VectorXf wcorr(qb.k);
    wcorr.setZero();
    for (int eigIdx = 0; eigIdx < qb.k; ++eigIdx) {
        float uDotB = 0.0f;
        for (int row = 0; row < qb.m; ++row) {
            const float u = qb.qComponent(row, eigIdx) / qb.sqrtLambda[eigIdx];
            uDotB += u * gwasEffects[row];
        }
        wcorr[eigIdx] = uDotB / qb.sqrtLambda[eigIdx];
    }
    return wcorr;
}
}

///////////////////////////////////////////////////////////////////////////////////////
////////    Step 1. perform eigen-decomposition for ld blocks                   ///////
///////////////////////////////////////////////////////////////////////////////////////


void Data::readLDBlockInfoFile(const string &ldBlockInfoFile){
    // Read bim file: recombination rate is defined between SNP i and SNP i-1
    ifstream in(ldBlockInfoFile.c_str());
    if (!in) throw ("Error: can not open the file [" + ldBlockInfoFile + "] to read.");
    cout << "Reading ld block info from file [" + ldBlockInfoFile + "]." << endl;
    ldBlockInfoVec.clear();
    ldBlockInfoMap.clear();
    string header;
    string id;
    int  chr,start, stop;
    int idx = 0;
    getline(in, header);
    while (in >> id >> chr >> start >> stop) {
        LDBlockInfo *ld = new LDBlockInfo(idx++, id, chr);
        ld->startPos    = start;
        ld->endPos     = stop;
        ldBlockInfoVec.push_back(ld);
        //chromosomes.insert(ld->chr);
        if (ldBlockInfoMap.insert(pair<string, LDBlockInfo*>(id, ld)).second == false) {
            throw ("Error: Duplicate LD block ID found: \"" + id + "\".");
        }
    }
    in.close();
    numLDBlocks = (unsigned) ldBlockInfoVec.size();
    cout << numLDBlocks << " LD Blocks to be included from [" + ldBlockInfoFile + "]." << endl;
}

void Data::eigenDecomposition( const MatrixXf &X, const float &prop, VectorXf &eigenValAdjusted, MatrixXf &eigenVecAdjusted, float &sumPosEigVal){
    // VectorXf cumsumNonNeg; // cumulative sums of non-negative values
    sumPosEigVal = 0.0;

    SelfAdjointEigenSolver<MatrixXf> eigensolver(X);
    VectorXf eigenVal = eigensolver.eigenvalues();
    MatrixXf eigenVec = eigensolver.eigenvectors();
    int n = (int)eigenVal.size();
    if (n == 0) {
        eigenValAdjusted.resize(0);
        eigenVecAdjusted.resize(0, 0);
        return;
    }
    // 1x1 matrix (e.g. sub-LD with one retained SNP): one eigenvalue explains all variance; the general loop uses revIdx-1 and breaks for n==1.
    if (n == 1) {
        if (eigenVal[0] < 0) cout << "Error, all eigenvector are negative" << endl;
        sumPosEigVal = (eigenVal[0] > 1e-10f) ? eigenVal[0] : 0.f;
        eigenValAdjusted.resize(1);
        eigenValAdjusted[0] = eigenVal[0];
        eigenVecAdjusted = eigenVec;
        return;
    }

    int revIdx = n;
    VectorXf cumsumNonNeg(revIdx);
    cumsumNonNeg.setZero();
    revIdx = revIdx -1;
    if(eigenVal[revIdx] < 0) cout << "Error, all eigenvector are negative" << endl;
    cumsumNonNeg[revIdx] = eigenVal[revIdx];
    sumPosEigVal = eigenVal[revIdx];
    revIdx = revIdx -1;
    
    int numPosEigVal = 0;
    while (revIdx >= 0 && eigenVal[revIdx] > 1e-10 ){
        sumPosEigVal = sumPosEigVal + eigenVal[revIdx];
        cumsumNonNeg[revIdx] = eigenVal[revIdx] + cumsumNonNeg[revIdx + 1];
        ++numPosEigVal;
        revIdx = revIdx - 1;
        if(revIdx < 0) break;
    }
    // cout << "revIdx: " << revIdx << endl;
    // cout << "size: " << eigenVal.size()  << " eigenVal: " << eigenVal << endl;
    // cout << "cumsumNoNeg: " << cumsumNonNeg << endl;
    cumsumNonNeg = cumsumNonNeg/sumPosEigVal;
    bool haveValue = false;
    // cout << "cumsumNonNeg: " << cumsumNonNeg << endl;
    // cout << "revIdx : " << revIdx << endl;
    // cout << "revIdx: "  << revIdx  << endl;
    for (revIdx = revIdx + 1; revIdx < eigenVal.size(); revIdx ++ ){
        // cout << "revIdx: " << revIdx << " cumsumNonNeg: " << cumsumNonNeg(revIdx) << endl;
        if(prop >= cumsumNonNeg[revIdx] ){
            revIdx = revIdx -1;
            haveValue = true;
            break;
        }
    }
    // cout << "revIdx : " << revIdx << endl;
    //if(!haveValue) revIdx = eigenVal.size() - 1;
    if(!haveValue) revIdx = numPosEigVal - 1;
    // cout << "cumsumNonNeg: " << cumsumNonNeg.size() << endl;
    // cout << "cumsumNoNeg: " << cumsumNonNeg << endl;
    // cout << "revIdx: "  << revIdx  << endl;

    eigenVecAdjusted = eigenVec.rightCols(eigenVal.size() - revIdx);
    eigenValAdjusted = eigenVal.tail(eigenVal.size() - revIdx);
    // cout << "eigenValAdjusted size: " << eigenValAdjusted.size() << " eigenValue eventually: " << eigenValAdjusted << endl;
    //eigenvalueNum = eigenVal.size() - revIdx;
    // cout << endl;
}

MatrixXf Data::generateLDmatrixPerBlock(const string &bedFile, const vector<string> &snplists){
    int numSnpInRange = snplists.size();
    IndInfo *indi = NULL;
    SnpInfo *snpj = NULL;
    SnpInfo *snpk = NULL;

    snpj = snpInfoMap.at(snplists[0]); // start
    snpk = snpInfoMap.at(snplists[numSnpInRange - 1]); // end;
    unsigned start = snpj->index;
    unsigned end = snpk->index;
    
    if (numIncdSnps == 0) throw ("Error: No SNP is retained for analysis.");
    if (numKeptInds == 0) throw ("Error: No individual is retained for analysis.");
    // if (start >= numIncdSnps) throw ("Error: Specified a SNP range of " + snpRange + " but " + to_string(static_cast<long long>(numIncdSnps)) + " SNPs are included.");
    
    // Gadget::Timer timer;
    // timer.setTime();
    
    //////////////////////////////////////////////////////
    // Step 1. read in the genotypes of SNPs in the given range
    //////////////////////////////////////////////////////
    
    const int bedToGeno[4] = {2, -9, 1, 0};
    unsigned size = (numInds+3)>>2;
    
    MatrixXf ZP(numSnpInRange, numKeptInds);  // SNP x Ind
    VectorXf Dtmp;
    Dtmp.setZero(numSnpInRange);

    if (numKeptInds < 2) throw("Error: Cannot calculate LD matrix with number of individuals < 2.");
    
    FILE *in1 = fopen(bedFile.c_str(), "rb");
    if (!in1) throw ("Error: can not open the file [" + bedFile + "] to read.");
    // cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
    char header[3];
    fread(header, sizeof(header), 1, in1);
    if (!in1 || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
        cerr << "Error: Incorrect first three bytes of bed file: " << bedFile << endl;
        exit(1);
    }

    int genoValue;
    unsigned i, j, k;
    unsigned incj, inck; // index of included SNP
    unsigned long long skipj = 0;
    unsigned nmiss;
    float mean;
    
    set<int> chromInRange;
    
    for (j = 0, incj = 0; j < numSnps; j++) {
        snpj = snpInfoVec[j];
        
        if (snpj->index < start || !snpj->included) {
            skipj += size;
            continue;
        }
        
        if (skipj) fseek(in1, skipj, SEEK_CUR);
        skipj = 0;
        
        char *bedLineIn = new char[size];
        fread(bedLineIn, sizeof(char), size, in1);
        
        chromInRange.insert(snpj->chrom);
        
        mean = 0.0;
        nmiss = 0;
        
        for (i = 0; i < numInds; i++) {
            indi = indInfoVec[i];
            if (!indi->kept) continue;
            genoValue = bedToGeno[(bedLineIn[i>>2]>>((i&3)<<1))&3];
            ZP(incj, indi->index) = genoValue;
            if (genoValue == -9) ++nmiss;
            else mean += genoValue;
        }
        delete[] bedLineIn;
        
        // fill missing values with the mean
        snpj->sampleSize = numKeptInds-nmiss;
        mean /= float(snpj->sampleSize);
        if (nmiss) {
            for (i=0; i<numKeptInds; ++i) {
                if (ZP(incj, i) == -9) ZP(incj, i) = mean;
            }
        }
        
        // compute allele frequency
        snpj->af = 0.5f*mean;
        snp2pq[incj] = snpj->twopq = 2.0f*snpj->af*(1.0f-snpj->af);
        
        if (snp2pq[incj]==0) throw ("Error: " + snpj->ID + " is a fixed SNP (MAF=0)!");
        
        Dtmp[incj] = Gadget::calcVariance(ZP.row(incj))*numKeptInds;
        
        ZP.row(incj) = (ZP.row(incj).array() - ZP.row(incj).mean())/sqrt(Dtmp[incj]);
        
        if (++incj == numSnpInRange) break;
    }
    
    fclose(in1);
    
    ZPZdiag = ZP.rowwise().squaredNorm();
    
    //////////////////////////////////////////////////////
    // Step 2. read in the bed file again to compute Z'Z
    //////////////////////////////////////////////////////
    
    MatrixXf denseZPZ;
    denseZPZ.setZero(numSnpInRange, numSnpInRange);
    VectorXf Zk(numKeptInds);
    Dtmp.setZero(numIncdSnps);
    
    FILE *in2 = fopen(bedFile.c_str(), "rb");
    fseek(in2, 3, SEEK_SET);
    unsigned long long skipk = 0;
    
    set<int>::iterator setend = chromInRange.end();

    if (numSkeletonSnps) {
        for (k = 0, inck = 0; k < numSnps; k++) {
            snpk = snpInfoVec[k];
            
            // if (!snpk->included) {
            if (snpk->index < start || !snpk->included) {
                skipk += size;
                continue;
            }

            if (chromInRange.find(snpk->chrom) == setend && !snpk->skeleton) {
                skipk += size;
                ++inck;       // ensure the index is correct
                continue;
            }
            
            if (skipk) fseek(in2, skipk, SEEK_CUR);
            skipk = 0;
            
            char *bedLineIn = new char[size];
            fread(bedLineIn, sizeof(char), size, in2);
            
            mean = 0.0;
            nmiss = 0;
            
            for (i = 0; i < numInds; i++) {
                indi = indInfoVec[i];
                if (!indi->kept) continue;
                genoValue = bedToGeno[(bedLineIn[i>>2]>>((i&3)<<1))&3];
                Zk[indi->index] = genoValue;
                if (genoValue == -9) ++nmiss;   // missing genotype
                else mean += genoValue;
            }
            delete[] bedLineIn;
            
            // fill missing values with the mean
            snpk->sampleSize = numKeptInds-nmiss;
            mean /= float(snpk->sampleSize);
            if (nmiss) {
                for (i=0; i<numKeptInds; ++i) {
                    if (Zk[i] == -9) Zk[i] = mean;
                }
            }
            
            // compute allele frequency
            snpk->af = 0.5f*mean;
            snp2pq[inck] = snpk->twopq = 2.0f*snpk->af*(1.0f-snpk->af);
            
            if (snp2pq[inck]==0) throw ("Error: " + snpk->ID + " is a fixed SNP (MAF=0)!");
            Dtmp[inck] = Gadget::calcVariance(Zk.row(inck))*numKeptInds;
            Zk = (Zk.array() - Zk.mean())/sqrt(Dtmp[inck]);
            denseZPZ.col(inck) = ZP * Zk;

            //++inck;
            if (++inck == numSnpInRange) break;
        }
    }
    else {
        for (k = 0, inck = 0; k < numSnps; k++) {
            snpk = snpInfoVec[k];
            
            // if (!snpk->included) {
            if (snpk->index < start || !snpk->included) {
                skipk += size;
                continue;
            }

            if (skipk) fseek(in2, skipk, SEEK_CUR);
            skipk = 0;
            
            char *bedLineIn = new char[size];
            fread(bedLineIn, sizeof(char), size, in2);
            
            mean = 0.0;
            nmiss = 0;
            
            for (i = 0; i < numInds; i++) {
                indi = indInfoVec[i];
                if (!indi->kept) continue;
                genoValue = bedToGeno[(bedLineIn[i>>2]>>((i&3)<<1))&3];
                Zk[indi->index] = genoValue;
                if (genoValue == -9) ++nmiss;   // missing genotype
                else mean += genoValue;
            }
            delete[] bedLineIn;
            
            // fill missing values with the mean
            snpk->sampleSize = numKeptInds-nmiss;
            mean /= float(snpk->sampleSize);
            if (nmiss) {
                for (i=0; i<numKeptInds; ++i) {
                    if (Zk[i] == -9) Zk[i] = mean;
                }
            }
            
            // compute allele frequency
            snpk->af = 0.5f*mean;
            snp2pq[inck] = snpk->twopq = 2.0f*snpk->af*(1.0f-snpk->af);
            
            if (snp2pq[inck]==0) throw ("Error: " + snpk->ID + " is a fixed SNP (MAF=0)!");

            Dtmp[inck] = Gadget::calcVariance(Zk)*numKeptInds;

            Zk = (Zk.array() - Zk.mean())/sqrt(Dtmp[inck]);
            
            denseZPZ.col(inck) = ZP * Zk;

            if (++inck == numSnpInRange) break;
        }
    }

    fclose(in2);
    // timer.getTime();
    return denseZPZ;
}

void Data::getEigenDataFromFullLDM(const string &filename, const float eigenCutoff = 0.995){
    
    string outfilename = filename + ".eigen.bin";
    FILE *out3 = fopen(outfilename.c_str(), "wb");

    for (unsigned blk=0; blk < numKeptLDBlocks; blk++){
        
        MatrixXf eigenVec;
        VectorXf eigenVal;
        float sumPosEigVal;
        // cout << "rval: " << rval << endl;
        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        eigenDecomposition(ZPZmat, eigenCutoff, eigenVal, eigenVec, sumPosEigVal);
        // cout << "rval: " << rval.row(0) << endl;
        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
        // save svd matrix
        
        int32_t numEigenValue = eigenVal.size();
        int32_t numSnpInBlock = keptLdBlockInfoVec[blk]->numSnpInBlock;  // TMP

        cout << " Generate Eigen decomposition result for LD block " << blk << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << endl;
        
        // save summary
        // 1, the number of SNPs in the block
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, out3);
        // 2, the number of eigenvalues at with the given cutoff
        fwrite(&numEigenValue, sizeof(int32_t), 1, out3);
        // 3. sum of all the positive eigenvalues
        fwrite(&sumPosEigVal, sizeof(float), 1, out3);
        // 4. eigenvalue cutoff based on the proportion of variance explained in LD
        fwrite(&eigenCutoff, sizeof(float), 1, out3);
        // 5. the selected eigenvalues
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, out3);
        // 6. the selected eigenvector;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, out3);
        
    }
    
    fclose(out3);

    //outputEigenDataForLDM(filename);
    
}

void Data::mapSnpsToBlocks(){
    //unsigned snpIdx = 0;
    int mapped = 0;
    ChromInfo* chr;
    LDBlockInfo *block;
    SnpInfo *snp;
    vector<ChromInfo*>::iterator it = chromInfoVec.begin(), end = chromInfoVec.end();
    unsigned chrCur = (*it)->id;
//    int chrStartSnpIdx = (*it)->startSnpIdx;
//    int chrEndSnpIdx = (*it)->endSnpIdx;
    
    for (unsigned i=0; i<numLDBlocks; ++i) {
        block = ldBlockInfoVec[i];
        block->snpNameVec.clear();
        block->snpInfoVec.clear();
        
//        cout << i << " " << numLDBlocks << " " << block->chrom << " " << chrCur << endl;
        
        if (block->chrom > chrCur) {  // move to a new chromosome
            if (it != end) ++it;
            if (it == end) {
                block->kept = false;
                continue;
            }
//            chrStartSnpIdx = (*it)->startSnpIdx;
//            chrEndSnpIdx = (*it)->endSnpIdx;
            chrCur = (*it)->id;
        }
        if (block->chrom != chrCur) {
            block->kept = false;
            continue;
        }
        
//        cout << i << " " << numLDBlocks << " " << block->chrom << " " << block->startPos << " " << block->endPos << endl;

        for (unsigned snpIdx = 0; snpIdx < numSnps; ++snpIdx) {
            //cout << "snpIdx " << snpIdx << " chrEndSnpIdx " << chrEndSnpIdx << endl;
            snp = snpInfoVec[snpIdx];
            if (!snp->included) continue;
            if (snp->chrom != chrCur) continue;
            if (snp->physPos < block->startPos) continue;
            else if (snp->physPos > block->endPos) break;
            block->snpNameVec.push_back(snp->ID);
            block->snpInfoVec.push_back(snp);
            snp->block = block->ID;
            snp->blockIdx = block->index;
//            cout << snpIdx << " " << snp->ID << " " << snp->chrom << " " << snp->physPos << " " << block->startPos << " " << block->endPos << " " << snp->included << endl;
        }
        block->numSnpInBlock = block->snpInfoVec.size();
        if (block->numSnpInBlock) {
            block->startSnpIdx = block->snpInfoVec[0]->index;
            block->endSnpIdx = block->snpInfoVec[block->numSnpInBlock-1]->index;
            block->kept = true;
            ++mapped;
        } else {
            block->kept = false;
        }
        //cout << i << " " << numLDBlocks << " " << block->numSnpInBlock << endl;

    }
        
    if (mapped < 1) throw(0, "No SNP can be mapped to the provided LD block list. Please check the input data regarding chromosome and bp.");
    else cout << mapped << " LD block(s) are retained." << endl;
}

void Data::makeBlockLDmatrix(const string &bedFile, const string &LDmatType, const unsigned block, const string &dirname, const bool writeLdmTxt, int ldBlockRegionWind){
    cout << "Making block LD matricies ..." << endl;
    
    struct stat sb;
    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        string create_cmd = "mkdir " + dirname;
        system(create_cmd.c_str());
        cout << "Created folder [" << dirname << "] to store LD matrices." << endl;
    }
    
    mapSnpsToBlocks();
    
    keptLdBlockInfoVec = makeKeptLDBlockInfoVec(ldBlockInfoVec);
    numKeptLDBlocks = (unsigned) keptLdBlockInfoVec.size();
    
#pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < numKeptLDBlocks; i++) {
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];

        string outBinfile = dirname + "/block" + ldblock->ID + ".ldm.bin";
        FILE *outbin = fopen(outBinfile.c_str(), "wb");
        ofstream outtxt;
        string outTxtfile;
        if (writeLdmTxt) {
            outTxtfile = dirname + "/block" + ldblock->ID + ".ldm.txt";
            outtxt.open(outTxtfile.c_str());
        }
        string outSnpfile = dirname + "/block" + ldblock->ID + ".snp.info";
        string outldmfile = dirname + "/block" + ldblock->ID + ".ldm.info";
        string pwldfilename = dirname + "/block" + ldblock->ID + ".rsq0.5.pwld";  // pairwise LD file

        if(!i) cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
            
        MatrixXf rval = generateLDmatrixPerBlock(bedFile, ldblock->snpNameVec);
        normaliseDenseCorrelation(rval);
        
//        unsigned numSnpInBlock = ldblock->numSnpInBlock;
//        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numSnpInBlock;
//        fwrite(rval.data(), sizeof(float), nElements, outbin);
        
        // save the lower trangular matrix only
        unsigned numSnpInBlock = ldblock->numSnpInBlock;
        for (unsigned row = 0; row < numSnpInBlock; ++row) {
            for (unsigned col = 0; col <= row; ++col) {
                float value = rval(row, col);
                fwrite(&value, sizeof(float), 1, outbin);
            }
        }
        
        if (writeLdmTxt) {
//            for (unsigned ii=0; ii<numSnpInBlock; ++ii){
//                for (unsigned jj=0; jj<numSnpInBlock; ++jj) {
//                    outtxt << ldblock->ID << "\t" << ldblock->snpNameVec[ii] << "\t" << ldblock->snpNameVec[jj] << "\t" << rval(ii,jj) << endl;
//                }
//            }
            for (unsigned row = 0; row < numSnpInBlock; ++row) {
                for (unsigned col = 0; col <= row; ++col) {
                    outtxt << rval(row, col) << "\t";
                }
                outtxt << endl;
            }
        }
        
        fclose(outbin);
        if (writeLdmTxt) outtxt.close();
        
        outputBlockLDmatrixInfo(*ldblock, outSnpfile, outldmfile);
        
        outputLDfriends(rval, ldblock, dirname, writeLdmTxt);
        
        
        if(!(i%1)) cout << " computed block " << ldblock->ID << "\r" << flush;

        if (block) {
            cout << "Written the LD matrix into file [" << outBinfile << "]." << endl;
            if (writeLdmTxt) cout << "Written the LD matrix into file [" << outTxtfile << "]." << endl;
            cout << "Written the LD matrix SNP info into file [" << outSnpfile << "]." << endl;
            cout << "Written the LD matrix ldm info into file [" << outldmfile << "]." << endl;
            cout << "Written the high pairwise LD correlations (rsq>0.5) into [" << pwldfilename << ".gz]." << endl;
        }
    }
    
    if (!block) {
        cout << "Written the LD matrix into folder [" << dirname << "/block*.ldm.bin]." << endl;
        if (writeLdmTxt) cout << "Written the LD matrix into text file [" << dirname << "/block*.ldm.txt]." << endl;
        
        //if (chromInfoVec.size() >= 22) {  // genome-wide build of LD matrices
        
            mergeLdmInfo(LDmatType, dirname, true);
        
        //} else {
//            cout << "Written the LD matrix into folder [" << dirname << "/block*.snp.info]." << endl;
//            cout << "Written the LD matrix into folder [" << dirname << "/block*.ldm.info]." << endl;
//            cout << "Written the high pairwise LD correlations (rsq>0.5) into [" << dirname << "/block*.rsq0.5.pwld]." << endl;
        //}
    }
}

void Data::makeCoreContextTiledBlockEigenFromEigen(const string &eigenDirname, const string &dirname, const float eigenCutoff, const unsigned tileSnps, const unsigned tileMinSnps, const float flankR2Threshold, const bool writeLdmTxt, const bool mergeCombinedInfo, const int outputQuantBits){
    if (!tileSnps) throw("Error: --tile-snps must be >= 1.");
    if (outputQuantBits != 0 && outputQuantBits != 16)
        throw("Error: --make-tiled-ldm-eigen currently supports only float or --bits 16 output.");
    if (!Gadget::directoryExist(eigenDirname)) throw("Error: cannot find the input eigen folder [" + eigenDirname + "].");
    if (!numLDBlocks) {
        readBlockLdmInfoFile(eigenDirname);
        readBlockLdmSnpInfoFile(eigenDirname);
    }
    if (coreContextLdmEigen) {
        throw("Error: --make-tiled-ldm-eigen from --ldm-eigen currently requires ordinary independent LD-block eigen input.");
    }

    const float r2Threshold = (flankR2Threshold >= 0.0f && flankR2Threshold < 1.0f) ? flankR2Threshold : 0.001f;
    cout << "Making core-context tiled LD eigen matrices from eigen data in [" << eigenDirname << "]..." << endl;
    cout << "  Core tile size: up to " << tileSnps
         << " SNPs; context: all sparse LD-connected SNPs in the original LD block; r^2 threshold: " << r2Threshold << endl;
    cout << "  Core tiles are balanced by SNP count; target minimum is " << tileMinSnps
         << " SNPs when compatible with the maximum." << endl;
    if (outputQuantBits == 16)
        cout << "  Writing q16 tiled eigen bundles and q16 sparse cross-tile LD files." << endl;
    if (r2Threshold > 0.001f) {
        cout << "Warning: sparse flank r^2 threshold > 0.001 can be too sparse for high-N GWAS and may inflate SNP effects. "
             << "Consider regenerating with --rsq 0.001 or lower." << endl;
    }

    const bool perSourceBlockOutput = !mergeCombinedInfo;
    if (perSourceBlockOutput)
        prepareOutputDirectoryPreserveExisting(dirname, "core-context tiled LD eigen matrices");
    else
        prepareCleanOutputDirectory(dirname, "core-context tiled LD eigen matrices");

    vector<float> snpNumInBlock;
    vector<LDBlockInfo*> tiledBlocks;
    unordered_map<int, int> newIndexByOldIndex;
    uint64_t totalSparseCrossTilePairs = 0;
    uint64_t totalContextDirectedLinks = 0;
    uint64_t totalContextTileFlankSnps = 0;
    int nextSnpIndex = 0;
    int blockId = 1;
    const string combinedSourcePwld = findPairwiseLdWithRecordsForBaseLocal(eigenDirname + "/rsq0.5.pwld");
    bool foundPairwiseLdSource = !combinedSourcePwld.empty();

    for (unsigned srcIdx = 0; srcIdx < ldBlockInfoVec.size(); ++srcIdx) {
        LDBlockInfo *srcBlock = ldBlockInfoVec[srcIdx];
        const int m = srcBlock->numSnpInBlock;
        if (m <= 0 || srcBlock->snpInfoVec.empty()) continue;
        SourceBlockOutputLock outputLock(dirname, srcBlock->ID, perSourceBlockOutput);
        if (perSourceBlockOutput) removeSourceBlockTiledOutputs(dirname, srcBlock->ID);

        int32_t cur_m = 0, cur_k = 0;
        float sumPosEigVal = 0.0f, oldEigenCutoff = 0.0f;
        VectorXf lambda;
        MatrixXf U;
        if (!readBlockEigenBinContents(eigenDirname, srcBlock->ID, m, cur_m, cur_k, sumPosEigVal, oldEigenCutoff, lambda, U)) {
            throw("Error: can not read float eigen data for input LD block " + srcBlock->ID + " from [" + eigenDirname + "].");
        }

        VectorXf lambdaUse;
        MatrixXf Uuse;
        if (eigenCutoff + 1e-6f < oldEigenCutoff) {
            truncateEigenMatrix(sumPosEigVal, eigenCutoff, lambda, U, lambdaUse, Uuse);
        } else {
            lambdaUse = lambda;
            Uuse = U;
        }
        if (!lambdaUse.size()) continue;

        VectorXf sqrtLambda = lambdaUse.array().sqrt();
        MatrixXf Qsource = sqrtLambda.asDiagonal() * Uuse.transpose();  // rank x original-block SNPs
        VectorXf qDiag(m);
        for (int j = 0; j < m; ++j) qDiag[j] = Qsource.col(j).squaredNorm();

        vector<vector<unsigned> > tileLocalIdx = makeBalancedSnpCountTiles(m, tileSnps, tileMinSnps);
        if (tileLocalIdx.empty()) continue;
        vector<int> tileByLocalIndex(m, -1);
        assignTileByLocalIndex(tileLocalIdx, m, tileByLocalIndex);

        for (unsigned t = 0; t < tileLocalIdx.size(); ++t) {
            for (unsigned a = 0; a < tileLocalIdx[t].size(); ++a) {
                SnpInfo *snp = srcBlock->snpInfoVec[tileLocalIdx[t][a]];
                newIndexByOldIndex[snp->index] = nextSnpIndex++;
            }
        }

        const string sparseFilename = sparseSourceLdFilename(dirname, srcBlock->ID, outputQuantBits == 16);
        uint64_t sparseCrossTilePairs = writeSparseSourceLdFromQ(sparseFilename,
                                                                 srcBlock->ID, srcBlock->snpInfoVec, newIndexByOldIndex,
                                                                 tileByLocalIndex, Qsource, qDiag, r2Threshold, outputQuantBits == 16);
        totalSparseCrossTilePairs += sparseCrossTilePairs;

        FILE *eigenBundle = openTiledEigenBundleForWrite(dirname, srcBlock->ID, (unsigned)tileLocalIdx.size(), outputQuantBits);
        vector<LDBlockInfo*> sourceTileBlocks;
        vector<VectorXf> sourceTileEigenvalues;
        vector<MatrixXf> sourceTileQ;
        for (unsigned t = 0; t < tileLocalIdx.size(); ++t) {
            const vector<unsigned> &coreLocal = tileLocalIdx[t];

            const int coreSize = (int)coreLocal.size();
            MatrixXf Qcore(Qsource.rows(), coreSize);
            for (int a = 0; a < coreSize; ++a) Qcore.col(a) = Qsource.col((int)coreLocal[a]);
            MatrixXf Rcc = Qcore.transpose() * Qcore;
            normaliseDenseCorrelation(Rcc);

            MatrixXf eigenVec;
            VectorXf eigenVal;
            float coreSumPosEigVal = 0.0f;
            eigenDecomposition(Rcc, eigenCutoff, eigenVal, eigenVec, coreSumPosEigVal);

            vector<SnpInfo*> coreSnps;
            coreSnps.reserve(coreSize);
            for (int a = 0; a < coreSize; ++a) {
                SnpInfo *old = srcBlock->snpInfoVec[coreLocal[a]];
                int newIdx = newIndexByOldIndex[old->index];
                SnpInfo *copy = new SnpInfo(newIdx, old->ID, old->a1, old->a2, old->chrom, old->genPos, old->physPos);
                copy->af = old->af;
                copy->twopq = old->twopq;
                copy->ld_n = old->ld_n;
                coreSnps.push_back(copy);
            }

            string id = srcBlock->ID + ".tile" + to_string(t + 1);
            ++blockId;
            LDBlockInfo *block = new LDBlockInfo((int)tiledBlocks.size(), id, srcBlock->chrom);
            block->startPos = coreSnps.front()->physPos;
            block->endPos = coreSnps.back()->physPos;
            block->startSnpIdx = coreSnps.front()->index;
            block->endSnpIdx = coreSnps.back()->index;
            block->numSnpInBlock = coreSize;
            block->snpInfoVec = coreSnps;
            block->sumPosEigVal = coreSumPosEigVal;
            block->sourceBlockID = srcBlock->ID;
            block->sourceTileID = (int)t + 1;
            for (unsigned a = 0; a < coreSnps.size(); ++a) block->snpNameVec.push_back(coreSnps[a]->ID);
            tiledBlocks.push_back(block);
            sourceTileBlocks.push_back(block);
            sourceTileEigenvalues.push_back(eigenVal);
            sourceTileQ.push_back(eigenVal.array().sqrt().matrix().asDiagonal() * eigenVec.transpose());

            writeTiledEigenBundleRecord(eigenBundle, block->ID, coreSize, eigenVal, eigenVec,
                                        coreSumPosEigVal, eigenCutoff, outputQuantBits);

            if (writeLdmTxt) {
                string outTxtfile = dirname + "/block" + block->ID + ".eigen.txt";
                ofstream outtxt(outTxtfile.c_str());
                outtxt << "Block " << block->ID << endl;
                outtxt << "NumSnp " << coreSize << endl;
                outtxt << "NumEigen " << eigenVal.size() << endl;
                outtxt << "EigenCutoff " << eigenCutoff << endl;
                outtxt << "Eigenvalues\n" << eigenVal.transpose() << endl;
                outtxt << "Eigenvectors\n" << eigenVec << endl;
            }

            snpNumInBlock.push_back((float)coreSize);
            const float retainedPct = coreSumPosEigVal > 0.0f ? 100.0f * eigenVal.sum() / coreSumPosEigVal : 0.0f;
            cout << "  Tile " << block->ID << ": " << coreSize << " SNPs, retained "
                 << eigenVal.size() << " eigenvalues";
            if (coreSumPosEigVal > 0.0f)
                cout << " (" << boost::format("%.2f") % retainedPct << "% LD variance)";
            cout << "." << endl;
        }
        fclose(eigenBundle);
        uint64_t sourceContextLinks = 0;
        uint64_t sourceContextFlankSnps = 0;
        writePrecomputedMarginalizedContextBundle(dirname, srcBlock->ID, sourceTileBlocks,
                                                   sourceTileEigenvalues, sourceTileQ, sparseFilename,
                                                   sourceContextLinks, sourceContextFlankSnps);
        totalContextDirectedLinks += sourceContextLinks;
        totalContextTileFlankSnps += sourceContextFlankSnps;
        outputCoreContextSourceInfo(sourceTileBlocks, dirname, srcBlock->ID, numKeptInds);
        if (!combinedSourcePwld.empty()) {
            writePwldSourceMarkerLocal(dirname, srcBlock->ID, combinedSourcePwld);
        } else {
            const string sourceBlockBase = eigenDirname + "/block" + srcBlock->ID + ".rsq0.5.pwld";
            const string sourceBlockPwld = findPairwiseLdWithRecordsForBaseLocal(sourceBlockBase);
            if (!sourceBlockPwld.empty()) {
                copyPairwiseLdToBaseLocal(sourceBlockPwld,
                    dirname + "/block" + srcBlock->ID + ".rsq0.5.pwld");
                foundPairwiseLdSource = true;
            }
        }
        cout << "  Source block " << srcBlock->ID << ": " << sparseCrossTilePairs
             << " sparse cross-tile SNP pairs with r^2 >= " << r2Threshold << "." << endl;
    }

    if (tiledBlocks.empty()) throw("Error: no core-context tiled LD eigen blocks could be generated.");
    if (!foundPairwiseLdSource) {
        cout << "Warning: no source rsq0.5.pwld file with LD records was found in [" << eigenDirname
             << "]. The tiled reference cannot provide pairwise LD for local credible sets; use --pwld-file with the original pairwise LD file." << endl;
    }

    if (perSourceBlockOutput) {
        cout << "\nWritten block-specific tiled LD eigen bundle for source block "
             << ldBlockInfoVec.front()->ID << " into folder [" << dirname << "]." << endl;
        cout << "Combined ldm.info/snp.info/pairwise LD files are not merged during --block jobs." << endl;
    } else {
        mergeLdmInfo("block", dirname, false);
    }

    cout << "\nWritten core-context tiled LD eigen bundles into folder [" << dirname
         << (outputQuantBits == 16 ? "/block*.tiles.eigen.q16.bin]." : "/block*.tiles.eigen.bin].") << endl;
    cout << "Written " << (outputQuantBits == 16 ? "q16 " : "float ")
         << "sparse cross-tile LD files into folder ["
         << dirname << (outputQuantBits == 16 ? "/block*.sparse.q16.bin]." : "/block*.sparse.bin].") << endl;
    cout << "Written source-block tile metadata into folder [" << dirname << "/block*.ldm.info]." << endl;
    cout << "Written precomputed marginalized-context bundles into folder [" << dirname
         << "/block*.tiles.mctx.bin] using " << totalContextDirectedLinks
         << " directed links to " << totalContextTileFlankSnps << " tile-specific flank SNPs." << endl;
    if (writeLdmTxt) cout << "Written tiled LD eigen text files into folder [" << dirname << "/block*.eigen.txt]." << endl;
    cout << "Generated " << snpNumInBlock.size() << " core tiles; per-tile retained eigen ranks are reported above." << endl;
    cout << "Total sparse cross-tile SNP pairs with r^2 >= " << r2Threshold << ": "
         << totalSparseCrossTilePairs << "." << endl;
}

void Data::makeCoreContextTiledBlockEigenFromBlockLD(const string &ldmDirname, const string &dirname, const float eigenCutoff, const unsigned tileSnps, const unsigned tileMinSnps, const float flankR2Threshold, const bool writeLdmTxt, const bool mergeCombinedInfo, const int outputQuantBits){
    if (!tileSnps) throw("Error: --tile-snps must be >= 1.");
    if (outputQuantBits != 0 && outputQuantBits != 16)
        throw("Error: --make-tiled-ldm-eigen currently supports only float or --bits 16 output.");
    if (!Gadget::directoryExist(ldmDirname)) throw("Error: cannot find the input block LD folder [" + ldmDirname + "].");
    if (!numLDBlocks) {
        readBlockLdmInfoFile(ldmDirname);
        readBlockLdmSnpInfoFile(ldmDirname);
    }
    if (coreContextLdmEigen) {
        throw("Error: --make-tiled-ldm-eigen from --ldm requires ordinary independent block LD input.");
    }

    const float r2Threshold = (flankR2Threshold >= 0.0f && flankR2Threshold < 1.0f) ? flankR2Threshold : 0.001f;
    cout << "Making core-context tiled LD eigen matrices from block LD matrices in [" << ldmDirname << "]..." << endl;
    cout << "  Core tile size: up to " << tileSnps
         << " SNPs; context: sparse LD-connected SNPs within each original LD block; r^2 threshold: " << r2Threshold << endl;
    cout << "  Core tiles are balanced by SNP count; target minimum is " << tileMinSnps
         << " SNPs when compatible with the maximum." << endl;
    if (outputQuantBits == 16)
        cout << "  Writing q16 tiled eigen bundles and q16 sparse cross-tile LD files." << endl;

    const bool perSourceBlockOutput = !mergeCombinedInfo;
    if (perSourceBlockOutput)
        prepareOutputDirectoryPreserveExisting(dirname, "core-context tiled LD eigen matrices");
    else
        prepareCleanOutputDirectory(dirname, "core-context tiled LD eigen matrices");

    vector<float> snpNumInBlock;
    vector<LDBlockInfo*> tiledBlocks;
    unordered_map<int, int> newIndexByOldIndex;
    uint64_t totalSparseCrossTilePairs = 0;
    uint64_t totalContextDirectedLinks = 0;
    uint64_t totalContextTileFlankSnps = 0;
    int nextSnpIndex = 0;

    for (unsigned srcIdx = 0; srcIdx < ldBlockInfoVec.size(); ++srcIdx) {
        LDBlockInfo *srcBlock = ldBlockInfoVec[srcIdx];
        if (!srcBlock || !srcBlock->kept) continue;
        const int mFull = srcBlock->numSnpInBlock;
        if (mFull <= 0 || srcBlock->snpInfoVec.empty()) continue;
        SourceBlockOutputLock outputLock(dirname, srcBlock->ID, perSourceBlockOutput);
        if (perSourceBlockOutput) removeSourceBlockTiledOutputs(dirname, srcBlock->ID);

        vector<unsigned> sourceLocal;
        vector<SnpInfo*> sourceSnps;
        for (int j = 0; j < mFull; ++j) {
            SnpInfo *snp = srcBlock->snpInfoVec[j];
            if (!snp || !snp->included) continue;
            sourceLocal.push_back((unsigned)j);
            sourceSnps.push_back(snp);
        }
        const int m = (int)sourceSnps.size();
        if (!m) continue;

        vector<vector<unsigned> > tileLocalIdx = makeBalancedSnpCountTiles(m, tileSnps, tileMinSnps);
        vector<int> tileByLocalIndex(m, -1);
        assignTileByLocalIndex(tileLocalIdx, m, tileByLocalIndex);
        if (tileLocalIdx.empty()) continue;

        for (unsigned t = 0; t < tileLocalIdx.size(); ++t) {
            for (unsigned a = 0; a < tileLocalIdx[t].size(); ++a) {
                SnpInfo *snp = sourceSnps[tileLocalIdx[t][a]];
                newIndexByOldIndex[snp->index] = nextSnpIndex++;
            }
        }

        vector<MatrixXf> tileRcc;
        uint64_t sparseCrossTilePairs = 0;
        const bool streamedLowerTri = buildCoreTileLdAndSparseFromLowerTriLdStream(
            ldmDirname, dirname, srcBlock->ID, mFull, sourceLocal, sourceSnps,
            newIndexByOldIndex, tileLocalIdx, tileByLocalIndex, r2Threshold, outputQuantBits, tileRcc,
            sparseCrossTilePairs);

        const string sparseFilename = sparseSourceLdFilename(dirname, srcBlock->ID, outputQuantBits == 16);
        MatrixXf Rsource;
        if (!streamedLowerTri) {
            cout << "\nWarning: source block " << srcBlock->ID
                 << " is not stored as a lower-triangular LD matrix; falling back to dense read." << endl;
            MatrixXf Rfull;
            readBlockLDmatrix(ldmDirname, srcBlock->ID, mFull, Rfull);
            normaliseDenseCorrelation(Rfull);
            Rsource.setZero(m, m);
            for (int j = 0; j < m; ++j) {
                for (int i = 0; i < m; ++i)
                    Rsource(i, j) = Rfull((int)sourceLocal[i], (int)sourceLocal[j]);
            }
            normaliseDenseCorrelation(Rsource);

            sparseCrossTilePairs = writeSparseSourceLdFromDense(sparseFilename,
                                                                 srcBlock->ID, sourceSnps, newIndexByOldIndex,
                                                                 tileByLocalIndex, Rsource, r2Threshold, outputQuantBits == 16);
            writePwldGzFromDense(dirname + "/block" + srcBlock->ID + ".rsq0.5.pwld.gz", sourceSnps, Rsource, 0.5f);
        }
        totalSparseCrossTilePairs += sparseCrossTilePairs;

        FILE *eigenBundle = openTiledEigenBundleForWrite(dirname, srcBlock->ID, (unsigned)tileLocalIdx.size(), outputQuantBits);
        vector<LDBlockInfo*> sourceTileBlocks;
        vector<VectorXf> sourceTileEigenvalues;
        vector<MatrixXf> sourceTileQ;
        for (unsigned t = 0; t < tileLocalIdx.size(); ++t) {
            const vector<unsigned> &coreLocal = tileLocalIdx[t];
            const int coreSize = (int)coreLocal.size();
            MatrixXf Rcc;
            if (streamedLowerTri) {
                Rcc.swap(tileRcc[t]);
            } else {
                Rcc.setZero(coreSize, coreSize);
                for (int a = 0; a < coreSize; ++a) {
                    for (int b = 0; b < coreSize; ++b)
                        Rcc(a, b) = Rsource((int)coreLocal[a], (int)coreLocal[b]);
                }
                normaliseDenseCorrelation(Rcc);
            }

            MatrixXf eigenVec;
            VectorXf eigenVal;
            float coreSumPosEigVal = 0.0f;
            eigenDecomposition(Rcc, eigenCutoff, eigenVal, eigenVec, coreSumPosEigVal);

            vector<SnpInfo*> coreSnps;
            coreSnps.reserve(coreSize);
            for (int a = 0; a < coreSize; ++a) {
                SnpInfo *old = sourceSnps[coreLocal[a]];
                int newIdx = newIndexByOldIndex[old->index];
                SnpInfo *copy = new SnpInfo(newIdx, old->ID, old->a1, old->a2, old->chrom, old->genPos, old->physPos);
                copy->af = old->af;
                copy->twopq = old->twopq;
                copy->ld_n = old->ld_n;
                coreSnps.push_back(copy);
            }

            string id = srcBlock->ID + ".tile" + to_string(t + 1);
            LDBlockInfo *block = new LDBlockInfo((int)tiledBlocks.size(), id, srcBlock->chrom);
            block->startPos = coreSnps.front()->physPos;
            block->endPos = coreSnps.back()->physPos;
            block->startSnpIdx = coreSnps.front()->index;
            block->endSnpIdx = coreSnps.back()->index;
            block->numSnpInBlock = coreSize;
            block->snpInfoVec = coreSnps;
            block->sumPosEigVal = coreSumPosEigVal;
            block->sourceBlockID = srcBlock->ID;
            block->sourceTileID = (int)t + 1;
            for (unsigned a = 0; a < coreSnps.size(); ++a) block->snpNameVec.push_back(coreSnps[a]->ID);
            tiledBlocks.push_back(block);
            sourceTileBlocks.push_back(block);
            sourceTileEigenvalues.push_back(eigenVal);
            sourceTileQ.push_back(eigenVal.array().sqrt().matrix().asDiagonal() * eigenVec.transpose());

            writeTiledEigenBundleRecord(eigenBundle, block->ID, coreSize, eigenVal, eigenVec,
                                        coreSumPosEigVal, eigenCutoff, outputQuantBits);

            if (writeLdmTxt) {
                string outTxtfile = dirname + "/block" + block->ID + ".eigen.txt";
                ofstream outtxt(outTxtfile.c_str());
                outtxt << "Block " << block->ID << endl;
                outtxt << "NumSnp " << coreSize << endl;
                outtxt << "NumEigen " << eigenVal.size() << endl;
                outtxt << "EigenCutoff " << eigenCutoff << endl;
                outtxt << "Eigenvalues\n" << eigenVal.transpose() << endl;
                outtxt << "Eigenvectors\n" << eigenVec << endl;
            }

            snpNumInBlock.push_back((float)coreSize);
            const float retainedPct = coreSumPosEigVal > 0.0f ? 100.0f * eigenVal.sum() / coreSumPosEigVal : 0.0f;
            cout << "  Tile " << block->ID << ": " << coreSize << " SNPs, retained "
                 << eigenVal.size() << " eigenvalues";
            if (coreSumPosEigVal > 0.0f)
                cout << " (" << boost::format("%.2f") % retainedPct << "% LD variance)";
            cout << "." << endl;
        }
        fclose(eigenBundle);
        uint64_t sourceContextLinks = 0;
        uint64_t sourceContextFlankSnps = 0;
        writePrecomputedMarginalizedContextBundle(dirname, srcBlock->ID, sourceTileBlocks,
                                                   sourceTileEigenvalues, sourceTileQ, sparseFilename,
                                                   sourceContextLinks, sourceContextFlankSnps);
        totalContextDirectedLinks += sourceContextLinks;
        totalContextTileFlankSnps += sourceContextFlankSnps;

        outputCoreContextSourceInfo(sourceTileBlocks, dirname, srcBlock->ID, numKeptInds);
        cout << "  Source block " << srcBlock->ID << ": " << sparseCrossTilePairs
             << " sparse cross-tile SNP pairs with r^2 >= " << r2Threshold << "." << endl;
    }

    if (tiledBlocks.empty()) throw("Error: no core-context tiled LD eigen blocks could be generated from block LD matrices.");

    if (perSourceBlockOutput) {
        cout << "\nWritten block-specific tiled LD eigen bundle for source block "
             << ldBlockInfoVec.front()->ID << " into folder [" << dirname << "]." << endl;
        cout << "Combined ldm.info/snp.info/pairwise LD files are not merged during --block jobs." << endl;
    } else {
        mergeLdmInfo("block", dirname, false);
    }

    cout << "\nWritten core-context tiled LD eigen bundles into folder [" << dirname
         << (outputQuantBits == 16 ? "/block*.tiles.eigen.q16.bin]." : "/block*.tiles.eigen.bin].") << endl;
    cout << "Written " << (outputQuantBits == 16 ? "q16 " : "float ")
         << "exact sparse cross-tile LD files into folder ["
         << dirname << (outputQuantBits == 16 ? "/block*.sparse.q16.bin]." : "/block*.sparse.bin].") << endl;
    cout << "Written source-block tile metadata into folder [" << dirname << "/block*.ldm.info]." << endl;
    cout << "Written precomputed marginalized-context bundles into folder [" << dirname
         << "/block*.tiles.mctx.bin] using " << totalContextDirectedLinks
         << " directed links to " << totalContextTileFlankSnps << " tile-specific flank SNPs." << endl;
    if (writeLdmTxt) cout << "Written tiled LD eigen text files into folder [" << dirname << "/block*.eigen.txt]." << endl;
    cout << "Generated " << snpNumInBlock.size() << " core tiles; per-tile retained eigen ranks are reported above." << endl;
    cout << "Total sparse cross-tile SNP pairs with r^2 >= " << r2Threshold << ": "
         << totalSparseCrossTilePairs << "." << endl;
}

void Data::makeBlockLDmatrixFromEigen(const string &eigenDirname, const float eigenCutoff, const string &LDmatType, const unsigned block, const string &dirname, const bool writeLdmTxt, int ldBlockRegionWind){
    (void)ldBlockRegionWind;
    cout << "Making block LD matrices from eigen decomposition in [" << eigenDirname << "] ..." << endl;
    
    struct stat sb;
    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        string create_cmd = "mkdir " + dirname;
        system(create_cmd.c_str());
        cout << "Created folder [" << dirname << "] to store LD matrices." << endl;
    }
    
    readBlockLdmInfoFile(eigenDirname, block);
    readBlockLdmSnpInfoFile(eigenDirname, block);
    readEigenMatrixBinaryFile(eigenDirname, eigenCutoff, false, ".");
    
#pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < numLDBlocks; i++) {
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (!ldblock->kept) continue;
        
        string outBinfile = dirname + "/block" + ldblock->ID + ".ldm.bin";
        FILE *outbin = fopen(outBinfile.c_str(), "wb");
        ofstream outtxt;
        string outTxtfile;
        if (writeLdmTxt) {
            outTxtfile = dirname + "/block" + ldblock->ID + ".ldm.txt";
            outtxt.open(outTxtfile.c_str());
        }
        string outSnpfile = dirname + "/block" + ldblock->ID + ".snp.info";
        string outldmfile = dirname + "/block" + ldblock->ID + ".ldm.info";
        string pwldfilename = dirname + "/block" + ldblock->ID + ".rsq0.5.pwld";  // pairwise LD file
        
        MatrixXf rval = eigenVecLdBlock[i] * eigenValLdBlock[i].asDiagonal() * eigenVecLdBlock[i].transpose();
        normaliseDenseCorrelation(rval);
        
        unsigned numSnpInBlock = ldblock->numSnpInBlock;
        for (unsigned row = 0; row < numSnpInBlock; ++row) {
            for (unsigned col = 0; col <= row; ++col) {
                float value = rval(row, col);
                fwrite(&value, sizeof(float), 1, outbin);
            }
        }
        
        if (writeLdmTxt) {
            for (unsigned row = 0; row < numSnpInBlock; ++row) {
                for (unsigned col = 0; col <= row; ++col) {
                    outtxt << rval(row, col) << "\t";
                }
                outtxt << endl;
            }
        }
        
        fclose(outbin);
        if (writeLdmTxt) outtxt.close();
        
        outputBlockLDmatrixInfo(*ldblock, outSnpfile, outldmfile);
        
        outputLDfriends(rval, ldblock, dirname, writeLdmTxt);
        
        
        if(!(i%1)) cout << " computed block " << ldblock->ID << "\r" << flush;

        if (block) {
            cout << "Written the LD matrix into file [" << outBinfile << "]." << endl;
            if (writeLdmTxt) cout << "Written the LD matrix into file [" << outTxtfile << "]." << endl;
            cout << "Written the LD matrix SNP info into file [" << outSnpfile << "]." << endl;
            cout << "Written the LD matrix ldm info into file [" << outldmfile << "]." << endl;
            cout << "Written the high pairwise LD correlations (rsq>0.5) into [" << pwldfilename << ".gz]." << endl;
        }
    }
    
    if (!block) {
        cout << "Written the LD matrix into folder [" << dirname << "/block*.ldm.bin]." << endl;
        if (writeLdmTxt) cout << "Written the LD matrix into text file [" << dirname << "/block*.ldm.txt]." << endl;
        
        mergeLdmInfo(LDmatType, dirname, true);
    }
}


void Data::outputBlockLDmatrixInfo(const LDBlockInfo &block, const string &outSnpfile, const string &outldmfile) const {
    // write snp info
    ofstream out1(outSnpfile.c_str());
    out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %12s %10s %10s\n")
    % "Chrom"
    % "ID"
    % "Index"
    % "GenPos"
    % "PhysPos"
    % "A1"
    % "A2"
    % "A1Freq"
    % "N"
    % "Block";
    SnpInfo *snp;
    for (unsigned i=0; i < block.numSnpInBlock; ++i) {
        snp = block.snpInfoVec[i];
        out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %12f %10s %10s\n")
        % snp->chrom
        % snp->ID
        % snp->index
        % snp->genPos
        % snp->physPos
        % snp->a1
        % snp->a2
        % snp->af
        % numKeptInds
        % block.ID;
    }
    out1.close();
    
    // svd matrix for ld blocks here.
    ofstream out2(outldmfile.c_str());
    out2 << boost::format("%10s %6s %15s %15s %15s %15s %12s\n")
    % "Block"
    % "Chrom"
    % "StartSnpIdx"
    % "StartSnpID"
    % "EndSnpIdx"
    % "EndSnpID"
    % "NumSnps";
    out2 << boost::format("%10s %6s %15s %15s %15s %15s %12s\n")
    % block.ID
    % block.chrom
    % block.startSnpIdx
    % block.snpInfoVec[0]->ID
    % block.endSnpIdx
    % block.snpInfoVec[block.numSnpInBlock-1]->ID
    % block.numSnpInBlock;
    out2.close();
}

//void Data::outputBlockLDmatrixInfo(const unsigned block, const string &dirname) const {
//    struct stat sb;
//    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
//        throw("Error: Folder " + dirname + " does not exist!");
//    }
//
//    string outSnpfile;
//    string outldmfile;
//
//    if (block) {
//        LDBlockInfo *blockinfo = keptLdBlockInfoVec[0];
//        outSnpfile = dirname + "/block" + blockinfo->ID + ".snp.info";
//        outldmfile = dirname + "/block" + blockinfo->ID + ".ldm.info";
//    } else {
//        outSnpfile = dirname + "/snp.info";
//        outldmfile = dirname + "/ldm.info";
//    }
//
//    // write snp info
//    ofstream out1(outSnpfile.c_str());
//    out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %12s %10s %10s\n")
//    % "Chrom"
//    % "ID"
//    % "Index"
//    % "GenPos"
//    % "PhysPos"
//    % "A1"
//    % "A2"
//    % "A1Freq"
//    % "N"
//    % "Block";
//    SnpInfo *snp;
//    for (unsigned i=0; i < numIncdSnps ; ++i) {
//        snp = incdSnpInfoVec[i];
//        out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %12f %10s %10s\n")
//        % snp->chrom
//        % snp->ID
//        % snp->index
//        % snp->genPos
//        % snp->physPos
//        % snp->a1
//        % snp->a2
//        % snp->af
//        % numKeptInds
//        % snp->block;
//    }
//    out1.close();
//
//    // svd matrix for ld blocks here.
//    ofstream out2(outldmfile.c_str());
//    out2 << boost::format("%10s %6s %15s %15s %15s %15s %12s\n")
//    % "Block"
//    % "Chrom"
//    % "StartSnpIdx"
//    % "StartSnpID"
//    % "EndSnpIdx"
//    % "EndSnpID"
//    % "NumSnps";
//    LDBlockInfo * ldblock;
//    for (unsigned i=0; i < numKeptLDBlocks ; ++i) {
//        ldblock = keptLdBlockInfoVec[i];
//        //cout << "ldblock id: " << ldblock->ID << endl;
//        out2 << boost::format("%10s %6s %15s %15s %15s %15s %12s\n")
//        % ldblock->ID
//        % ldblock->chrom
//        % ldblock->startSnpIdx
//        % incdSnpInfoVec[ldblock->startSnpIdx]->ID
//        % ldblock->endSnpIdx
//        % incdSnpInfoVec[ldblock->endSnpIdx]->ID
//        % ldblock->snpNameVec.size();
//    }
//    out2.close();
//
//    cout << "Written the LD matrix SNP info into file [" << outSnpfile << "]." << endl;
//    cout << "Written the LD matrix ldm info into file [" << outldmfile << "]." << endl;
//}


//void Data::makeBlockLDmatrix(const string &bedFile, const string &LDmatType, const unsigned block, const string &dirname, const bool writeLdmTxt, int ldBlockRegionWind){
//    int i,j;
//    vector<locus_bp> snpVec;
//    SnpInfo *snp;
//
//    map<int, string>  chrEndSnp;
//    for (i = 1; i < numIncdSnps; i++) {
//        snp = incdSnpInfoVec[i];
//        if(incdSnpInfoVec[i]->chrom != incdSnpInfoVec[i-1]->chrom){
//            chrEndSnp.insert(pair<int, string>(incdSnpInfoVec[i - 1]->chrom,incdSnpInfoVec[i - 1]->ID ));
//        }
//    }
//    chrEndSnp.insert(pair<int, string>(incdSnpInfoVec[numIncdSnps - 1]->chrom,incdSnpInfoVec[numIncdSnps - 1]->ID ));
//    //Step 1.2  Read block file
//    //readLDBlockInfoFile(ldBlockInfoFile);
//
//
//    /////////////////////////////////////////
//    // Step 2. Map snps to blocks
//    /////////////////////////////////////////
//    vector<string> block2snp_1(numLDBlocks), block2snp_2(numLDBlocks);
//    map<string,int> keptLdBlock2AllLdBlcokMap;
//    vector<locus_bp>::iterator iter;
//    map<int, string>::iterator chrIter;
//    LDBlockInfo *ldblock;
//    for (i = 0; i < numIncdSnps ; i++) {
//        snp = incdSnpInfoVec[i];
//        snpVec.push_back(locus_bp(snp->ID, snp->chrom, snp->physPos ));
//    }
//#pragma omp parallel for private(iter, chrIter)
//    for (i = 0; i < numLDBlocks; i++) {
//        // find lowest snp_name in the block
//        ldblock = ldBlockInfoVec[i];
//
//        iter = find_if(snpVec.begin(), snpVec.end(), locus_bp( ldblock->ID ,ldblock->chrom, ldblock->startPos - ldBlockRegionWind));
//        if (iter != snpVec.end()) block2snp_1[i] = iter->locusName;
//        else block2snp_1[i] = "NA";
//    }
//#pragma omp parallel for private(iter, chrIter)
//    for (i = 0; i < numLDBlocks; i++) {
//        ldblock = ldBlockInfoVec[i];
//        if (block2snp_1[i] == "NA") {
//            block2snp_2[i] = "NA";
//            continue;
//        }
//        iter = find_if(snpVec.begin(), snpVec.end(), locus_bp(ldblock->ID, ldblock->chrom, ldblock->endPos + ldBlockRegionWind));
//        if (iter != snpVec.end()){
//            if (iter->bp ==  ldblock->endPos + ldBlockRegionWind){
//                block2snp_2[i] = iter->locusName;
//            }else {
//                if(iter!=snpVec.begin()){
//                    iter--;
//                    block2snp_2[i] = iter->locusName;
//                }
//                else block2snp_2[i] = "NA";
//            }
//        }
//        else {
//            chrIter = chrEndSnp.find(ldblock->chrom);
//            if (chrIter == chrEndSnp.end()) block2snp_2[i] = "NA";
//            else block2snp_2[i] = chrIter->second;
//        }
//    }
//    int mapped = 0;
//    for (i = 0; i < numLDBlocks; i++) {
//        ldblock = ldBlockInfoVec[i];
//        if (block2snp_1[i] != "NA" && block2snp_2[i] != "NA")
//        {
//            mapped++;
//            // ldblock->kept = true;
//            keptLdBlock2AllLdBlcokMap.insert(pair<string, int>(ldblock->ID, i));
//        } else {
//            ldblock->kept = false;
//        }
//    }
//
//    if (mapped < 1) throw(0, "No SNP can be mapped to the provided LD block list. Please check the input data regarding chromosome and bp.");
//    else cout << mapped << " LD block(s) are retained." << endl;
//
//    keptLdBlockInfoVec = makeKeptLDBlockInfoVec(ldBlockInfoVec);
//    numKeptLDBlocks = (unsigned) keptLdBlockInfoVec.size();
//
//    map<string, int>::iterator iter1, iter2;
//    map<string, int> snp_name_map;
//    VectorXf snpNumInldblock(numKeptLDBlocks);
//    VectorXf eigenNumInldblock(numKeptLDBlocks);
//    vector<VectorXf> cumsumNonNeg(numKeptLDBlocks);
//
//    for (i = 0; i < numIncdSnps; i++) {
//        snp = incdSnpInfoVec[i];
//        snp_name_map.insert(pair<string,int>(snp->ID, i));
//
//    }
//    // eigenvector and eigenvalue
//    //eigenValLdBlock.resize(numKeptLDBlocks);
//    //eigenVecLdBlock.resize(numKeptLDBlocks);
//
//    //    string outfilename = filename + "." + LDmatType + ".eigen" + ".bin";
//    //    FILE *out = fopen(outfilename.c_str(), "wb");
//
//    struct stat sb;
//    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
//        // Folder doesn't exist, create it
//        string create_cmd = "mkdir " + dirname;
//        system(create_cmd.c_str());
//        cout << "Created folder [" << dirname << "] to store results." << endl;
//
//    }
//
//    string outBinfile;
//    string outTxtfile;
//    bool readBedBool = true;
//    for (i = 0; i < numKeptLDBlocks; i++) {
//        ldblock = keptLdBlockInfoVec[i];
//
//        outBinfile = dirname + "/block" + ldblock->ID + ".ldm.bin";
//        FILE *outbin = fopen(outBinfile.c_str(), "wb");
//        ofstream outtxt;
//        outTxtfile;
//        if (writeLdmTxt) {
//            outTxtfile = dirname + "/block" + ldblock->ID + ".ldm.txt";
//            outtxt.open(outTxtfile.c_str());
//        }
//
//        // cout << "ldblock id: " << ldblock->ID << endl;
//        iter1 = snp_name_map.find(block2snp_1[keptLdBlock2AllLdBlcokMap.at(ldblock->ID ) ]);
//        iter2 = snp_name_map.find(block2snp_2[keptLdBlock2AllLdBlcokMap.at(ldblock->ID ) ]);
//
//        bool skip = false;
//        if (iter1 == snp_name_map.end() || iter2 == snp_name_map.end() || iter1->second >= iter2->second) ldblock->kept = false;
//        snpNumInldblock[i] = iter2->second - iter1->second + 1;
//        // cout << "ldblock->kept: " << ldblock->kept << endl;
//        if(!ldblock->kept) continue;
//        vector<int> snp_indx;
//        SnpInfo *snp;
//        for (j = iter1->second; j <= iter2->second; j++) {
//            snp_indx.push_back(j);
//            snp = incdSnpInfoVec[j];
//            ldblock->snpNameVec.push_back(snp->ID);
//            snp->block = ldblock->ID;
//            // cout << incdSnpInfoVec[j]->ID << " ";
//        }
//        if(readBedBool) {
//            cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
//            readBedBool = false;
//        }
//        //        MatrixXf eigenVec;
//        //        VectorXf eigenVal,cumsumNonNegPerLD;
//        MatrixXf rval = generateLDmatrixPerBlock(bedFile, ldblock->snpNameVec);
//        //        // cout << "rval: " << rval << endl;
//        //        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
//        //        eigenDecomposition(rval, eigenCutoff,eigenVal, eigenVec,cumsumNonNegPerLD);
//        //        // cout << "rval: " << rval.row(0) << endl;
//        //        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
//        //        // save svd matrix
//        //        int32_t numEigenValue = eigenVal.size();
//        int32_t numSnpInBlock = ldblock->snpNameVec.size();
//
//        ldblock->startSnpIdx = snp_indx[0];
//        ldblock->endSnpIdx = snp_indx[snp_indx.size()-1];
//
//        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numSnpInBlock;
//        fwrite(rval.data(), sizeof(float), nElements, outbin);
//        //        cout << " Generate and save Eigen decomposition result for LD block " << i << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << "\r" << flush;
//
//        if (writeLdmTxt) {
//            for (unsigned ii=0; ii<numSnpInBlock; ++ii){
//                for (unsigned jj=0; jj<numSnpInBlock; ++jj) {
//                    outtxt << ldblock->ID << "\t" << ldblock->snpNameVec[ii] << "\t" << ldblock->snpNameVec[jj] << "\t" << rval(ii,jj) << endl;
//                }
//            }
//        }
//
//        fclose(outbin);
//        if (writeLdmTxt) outtxt.close();
//    }
//
//    if (block) {
//        cout << "Written the LD matrix into file [" << outBinfile << "]." << endl;
//        if (writeLdmTxt) cout << "Written the LD matrix into file [" << outTxtfile << "]." << endl;
//    }
//    else {
//        cout << "Written the LD matrix into folder [" << dirname << "/block*.ldm.bin]." << endl;
//        if (writeLdmTxt) cout << "Written the LD matrix into text file [" << dirname << "/block*.ldm.txt]." << endl;
//    }
//
//    outputBlockLDmatrixInfo(block, dirname);
//
//}


void Data::impG(const unsigned block, double diag_mod){
    VectorXi numImpSnp;
    VectorXi numValidTyp;  // included in GWAS with non-zero SE (anchors for LD imputation)
    numImpSnp.setZero(numLDBlocks);
    numValidTyp.setZero(numLDBlocks);
    for (unsigned i = 0; i < numLDBlocks; i++ ){
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (!ldblock->kept) continue;
        for (unsigned j=0; j<ldblock->numSnpInBlock; ++j) {
            SnpInfo *snp = ldblock->snpInfoVec[j];
            if (gwasUsableAnchorSnp(snp)) {
                ++numValidTyp[i];
            } else {
                // Not in GWAS, or present but SE missing/zero/non-finite — impute from LD with anchors
                ++numImpSnp[i];
            }
        }
    }
    unsigned totalNumImpSnp = numImpSnp.sum();
    
    cout << boost::format("%12s %12s %12s %12s %12s\n") % "Block" % "TotalSNPs" % "Anchors" % "ToImpute" % "PctImpute";
    for (unsigned i = 0; i < numLDBlocks; i++ ){
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (!ldblock->kept) continue;
        cout << boost::format("%12s %12s %12s %12s %12.3f\n") % (i+1) % ldblock->numSnpInBlock % numValidTyp[i] % numImpSnp[i] % (float(numImpSnp[i])/float(ldblock->numSnpInBlock));
        if (numValidTyp[i] == 0) {
            cout << "  Warning: Block " << (i+1) << " has no anchor SNP (finite beta and SE>0); block skipped." << endl;
            ldblock->kept = false;
            if (block) {
                // A header-only file is an explicit completion marker for an
                // intentionally empty block.  Without it, the merge step
                // cannot distinguish this case from a failed block job.
                string outfile = title + ".block" + ldblock->ID + ".imputed.ma";
                ofstream out(outfile.c_str());
                if (!out)
                    throw("Error: can not open the file [" + outfile + "] to write.");
                out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n")
                    % "SNP" % "A1" % "A2" % "freq" % "b" % "se" % "p" % "N";
                out.close();
                cout << "Created empty completed-block summary file [" << outfile << "]." << endl;
            }
        }
    }
        
    cout << "Imputing summary statistics for " << to_string(totalNumImpSnp) << " SNPs (not in GWAS and/or SE=0)." << endl;

    Gadget::Timer timer;
    timer.setTime();
    
#pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < numLDBlocks; i++ ){
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (!ldblock->kept) continue;
        
        if (numImpSnp[i] && numValidTyp[i]) {
            
            Stat::Normal normal;
            
            /// Step 1. construct LD
            MatrixXf LDPerBlock = eigenVecLdBlock[i] * eigenValLdBlock[i].asDiagonal() * eigenVecLdBlock[i].transpose();
            
            LDPerBlock.diagonal().array() += (float)diag_mod;
            /// Step 2. Construct the LD correlation matrix among anchor SNPs (LDtt) and between imputation targets and anchors (LDit).
            // Anchors: gwasUsableAnchorSnp (finite b, finite SE>0).  Else imputed.
            VectorXi typedSnpIdx(numValidTyp[i]);
            VectorXi untypedSnpIdx(numImpSnp[i]);
            VectorXf zTypSnp(numValidTyp[i]);
            VectorXf nTypSnp(numValidTyp[i]);
            VectorXf varyTypSnp(numValidTyp[i]);
            for(unsigned j=0, idxTyp=0, idxImp=0; j < ldblock->numSnpInBlock; j++){
                SnpInfo *snp = ldblock->snpInfoVec[j];
                if (gwasUsableAnchorSnp(snp)) {
                    typedSnpIdx[idxTyp] = j;
                    zTypSnp[idxTyp] = float(snp->gwas_b / snp->gwas_se);
                    nTypSnp[idxTyp] = snp->gwas_n;
                    float hetj = 2.0 * snp->gwas_af * (1.0 - snp->gwas_af);
                    varyTypSnp[idxTyp] = hetj * (snp->gwas_n * snp->gwas_se * snp->gwas_se + snp->gwas_b * snp->gwas_b);
                    ++idxTyp;
                } else {
                    untypedSnpIdx[idxImp] = j;
                    ++idxImp;
                }
            }
            // Step 2.2 construct LDtt and LDit and Ztt.
            MatrixXf LDtt = LDPerBlock(typedSnpIdx,typedSnpIdx);
            MatrixXf LDit = LDPerBlock(untypedSnpIdx,typedSnpIdx);
            // Step 2.3 //  The Z score for the missing SNPs;
            VectorXf LDi_Z = LDtt.ldlt().solve(zTypSnp);
            VectorXf zImpSnp = LDit * LDi_Z;
            // Step 3. re-calcualte beta and se
            // if snp is missing use median to replace N
            std::sort(nTypSnp.data(), nTypSnp.data() + nTypSnp.size());
            float nMedian = nTypSnp[nTypSnp.size()/2];  // median
            // calcuate median of phenotypic variance
            std::sort(varyTypSnp.data(), varyTypSnp.data() + varyTypSnp.size());
            float varyMedian = varyTypSnp[varyTypSnp.size()/2];  // median
            // begin impute
            for(unsigned j = 0; j < numImpSnp[i]; j++){
                SnpInfo *snp = ldblock->snpInfoVec[untypedSnpIdx[j]];
                if (snp->keepSumstatIntact) continue;
                float base = float(sqrt(2.0 * snp->af *(1.0 - snp->af) * (nMedian + zImpSnp[j] * zImpSnp[j])));
                if (!(base > 1e-12f) || !std::isfinite(base) || !std::isfinite(zImpSnp[j])) {
                    cout << "Warning: cannot impute SNP " << snp->ID << " (MAF at 0/1 in reference or non-finite Z); excluding from analysis." << endl;
                    snp->included = false;
                    snp->gwas_b = -999;
                    snp->gwas_se = -999;
                    continue;
                }
                snp->gwas_b = zImpSnp[j] * sqrt(varyMedian)/base;
                snp->gwas_se = sqrt(varyMedian) / base;
                snp->gwas_n = nMedian;
                snp->gwas_af = snp->af;
                {
                    double denom = std::max(snp->gwas_se, 1e-20);
                    double zpv = std::abs(snp->gwas_b / denom);
                    if (!std::isfinite(zpv)) zpv = 0.0;
                    snp->gwas_pvalue = float(2*(1.0-normal.cdf_01(zpv)));
                }
                snp->included = true;
                //cout << "b " << snp->gwas_b << " se " << snp->gwas_se << " z " << snp->gwas_b/snp->gwas_se << " p " << snp->gwas_pvalue << endl;
            }
        }
        
        if(!(i%10)) cout << " imputed block " << i+1 << "\r" << flush;
//        cout << " imputed block " << i << endl;

        if (block) {
            string outfile = title + ".block" + ldblock->ID + ".imputed.ma";
            ofstream out(outfile.c_str());
            out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n") % "SNP" % "A1" % "A2" % "freq" % "b" % "se" % "p" % "N";
            for (unsigned i=0; i<ldblock->numSnpInBlock; ++i) {
                SnpInfo *snp = ldblock->snpInfoVec[i];
                out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n")
                % snp->ID
                % snp->a1
                % snp->a2
                % snp->gwas_af
                % snp->gwas_b
                % snp->gwas_se
                % snp->gwas_pvalue
                % snp->gwas_n;
            }
            out.close();

            timer.getTime();
            cout << "Imputation of summary statistics is completed (time used: " << timer.format(timer.getElapse()) << ")." << endl;
            cout << "Summary statistics of all SNPs are saved into file [" + outfile + "]." << endl;

        }

    }

    if (block) return;
    
    string outfile = title + ".imputed.ma";
    ofstream out(outfile.c_str());
    out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n") % "SNP" % "A1" % "A2" % "freq" % "b" % "se" % "p" % "N";
    for (unsigned bi=0; bi<numLDBlocks; ++bi) {
        LDBlockInfo *ldblock = ldBlockInfoVec[bi];
        if (!ldblock->kept) continue;
        for (unsigned j=0; j<ldblock->numSnpInBlock; ++j) {
            SnpInfo *snp = ldblock->snpInfoVec[j];
            out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n")
            % snp->ID
            % snp->a1
            % snp->a2
            % snp->gwas_af
            % snp->gwas_b
            % snp->gwas_se
            % snp->gwas_pvalue
            % snp->gwas_n;
        }
    }
    out.close();

    timer.getTime();
    cout << "Imputation of summary statistics is completed (time used: " << timer.format(timer.getElapse()) << ")." << endl;
    cout << "Summary statistics of all SNPs are saved into file [" + outfile + "]." << endl;
}



void Data::getEigenDataForLDBlock(const string &bedFile, const string &ldBlockInfoFile, int ldBlockRegionWind, const string &filename, const float eigenCutoff){
    int i,j;
    vector<locus_bp> snpVec;
    SnpInfo *snp;

    map<int, string>  chrEndSnp;
    for (i = 1; i < numIncdSnps; i++) {
        snp = incdSnpInfoVec[i];
        if(incdSnpInfoVec[i]->chrom != incdSnpInfoVec[i-1]->chrom){
            chrEndSnp.insert(pair<int, string>(incdSnpInfoVec[i - 1]->chrom,incdSnpInfoVec[i - 1]->ID ));
        }
    }
    chrEndSnp.insert(pair<int, string>(incdSnpInfoVec[numIncdSnps - 1]->chrom,incdSnpInfoVec[numIncdSnps - 1]->ID ));
    //Step 1.2  Read block file
    readLDBlockInfoFile(ldBlockInfoFile);
    /////////////////////////////////////////
    // Step 2. Map snps to blocks
    /////////////////////////////////////////
    vector<string> block2snp_1(numLDBlocks), block2snp_2(numLDBlocks);
    map<string,int> keptLdBlock2AllLdBlcokMap;
    vector<locus_bp>::iterator iter;
    map<int, string>::iterator chrIter;
    for (i = 0; i < numIncdSnps ; i++) {
        snp = incdSnpInfoVec[i];
        snpVec.push_back(locus_bp(snp->ID, snp->chrom, snp->physPos ));
    }
#pragma omp parallel for private(iter, chrIter)
    for (i = 0; i < numLDBlocks; i++) {
        // find lowest snp_name in the block
        LDBlockInfo *ldblock = ldBlockInfoVec[i];

        iter = find_if(snpVec.begin(), snpVec.end(), locus_bp( ldblock->ID ,ldblock->chrom, ldblock->startPos - ldBlockRegionWind));
        if (iter != snpVec.end()) block2snp_1[i] = iter->locusName;
        else block2snp_1[i] = "NA";
    }
#pragma omp parallel for private(iter, chrIter)
    for (i = 0; i < numLDBlocks; i++) {
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (block2snp_1[i] == "NA") {
            block2snp_2[i] = "NA";
            continue;
        }
        iter = find_if(snpVec.begin(), snpVec.end(), locus_bp(ldblock->ID, ldblock->chrom, ldblock->endPos + ldBlockRegionWind));
        if (iter != snpVec.end()){
            if (iter->bp ==  ldblock->endPos + ldBlockRegionWind){
                block2snp_2[i] = iter->locusName;
            }else {
                if(iter!=snpVec.begin()){
                    iter--;
                    block2snp_2[i] = iter->locusName;
                }
                else block2snp_2[i] = "NA";
            }
        }
        else {
            chrIter = chrEndSnp.find(ldblock->chrom);
            if (chrIter == chrEndSnp.end()) block2snp_2[i] = "NA";
            else block2snp_2[i] = chrIter->second;
        }
    }
    int mapped = 0;
    for (i = 0; i < numLDBlocks; i++) {
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (block2snp_1[i] != "NA" && block2snp_2[i] != "NA")
        {
            mapped++;
            // ldblock->kept = true;
            keptLdBlock2AllLdBlcokMap.insert(pair<string, int>(ldblock->ID, i));
        } else {
            ldblock->kept = false;
        }
    }
    if (mapped < 1) throw(0, "No SNP can be mapped to the provided LD block list. Please check the input data regarding chromosome and bp.");
    else cout << mapped << " LD blocks have at least one SNP." << endl;

    keptLdBlockInfoVec = makeKeptLDBlockInfoVec(ldBlockInfoVec);
    numKeptLDBlocks = (unsigned) keptLdBlockInfoVec.size();

    map<string, int>::iterator iter1, iter2;
    map<string, int> snp_name_map;
    VectorXf snpNumInldblock(numKeptLDBlocks);
    VectorXf eigenNumInldblock(numKeptLDBlocks);
    vector<VectorXf> cumsumNonNeg(numKeptLDBlocks);
    
    for (i = 0; i < numIncdSnps; i++) {
        SnpInfo *snp = incdSnpInfoVec[i];
        snp_name_map.insert(pair<string,int>(snp->ID, i));
        
    }
      // eigenvector and eigenvalue
    //eigenValLdBlock.resize(numKeptLDBlocks);
    //eigenVecLdBlock.resize(numKeptLDBlocks);
    string LDmatType = "block";
    string outfilename = filename + "." + LDmatType + ".eigen" + ".bin";
    FILE *out3 = fopen(outfilename.c_str(), "wb");

    bool readBedBool = true;
    for (i = 0; i < numKeptLDBlocks; i++) {
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
        // cout << "ldblock id: " << ldblock->ID << endl;
        iter1 = snp_name_map.find(block2snp_1[keptLdBlock2AllLdBlcokMap.at(ldblock->ID ) ]);
        iter2 = snp_name_map.find(block2snp_2[keptLdBlock2AllLdBlcokMap.at(ldblock->ID ) ]);
        bool skip = false;
        if (iter1 == snp_name_map.end() || iter2 == snp_name_map.end() || iter1->second >= iter2->second) ldblock->kept = false;
        snpNumInldblock[i] = iter2->second - iter1->second + 1;
        // cout << "ldblock->kept: " << ldblock->kept << endl;
        if(!ldblock->kept) continue;
        vector<int> snp_indx;
        for (j = iter1->second; j <= iter2->second; j++) {
            snp_indx.push_back(j);
            ldblock->snpNameVec.push_back(incdSnpInfoVec[j]->ID);
           // cout << incdSnpInfoVec[j]->ID << " ";
        }
        if(readBedBool) {
            cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
            readBedBool = false;
        }
        MatrixXf eigenVec;
        VectorXf eigenVal;
        float sumPosEigVal = 0;
        MatrixXf rval = generateLDmatrixPerBlock(bedFile, ldblock->snpNameVec);
        // cout << "rval: " << rval << endl;
        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        eigenDecomposition(rval, eigenCutoff,eigenVal, eigenVec, sumPosEigVal);
        // cout << "rval: " << rval.row(0) << endl;
        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
        // save svd matrix
        int32_t numEigenValue = eigenVal.size();
        int32_t numSnpInBlock = ldblock->snpNameVec.size();
        eigenNumInldblock[i] = numEigenValue;
        // save summary
        // 1, nrow of eigenVecGene[i]
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, out3);
        //        cout << "rval: " << rval << endl;
        // cout << "eigenVec: "  << eigenVec << endl;
        // cout << "";
        // 2, ncol of eigenVecGene[i]
        fwrite(&numEigenValue, sizeof(int32_t), 1, out3);
        // 3. sum of all positive eigenvalues
        fwrite(&sumPosEigVal, sizeof(float), 1, out3);
        // 4. eigenCutoff
        fwrite(&eigenCutoff, sizeof(float), 1, out3);
        // 5. eigenvalues
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, out3);
        // 6. eigenvectors;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, out3);
        cout << " Generate and save Eigen decomposition result for LD block " << i << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << "\r" << flush;
    }
    fclose(out3);
    // cout << "size of eigenValuegene: " << eigenValGene.size() << endl;
    //outputEigenDataForLDM(filename);
    
    cout << "To explain " << eigenCutoff*100 << "% variance in LD, on average " << int(eigenNumInldblock.mean()) << " eigenvalues are selected across LD blocks (mean number of SNPs is " << int(snpNumInldblock.mean()) << ")." << endl;

}


void Data::readBlockLdmInfoFile(const string &dirname, const unsigned block){
    // Read bim file: recombination rate is defined between SNP i and SNP i-1
    string infoFile;
    if (block) {
        infoFile = dirname + "/block" + to_string(block) + ".ldm.info";
        ifstream in(infoFile.c_str());
        if (!in) infoFile = dirname + "/ldm.info";
    } else {
        infoFile = dirname + "/ldm.info";
    }
    ifstream in(infoFile.c_str());
    if (!in) throw ("Error: can not open the file [" + infoFile + "] to read.");
    cout << "Reading LDM info from file [" + infoFile + "]." << endl;
    ldBlockInfoVec.clear();
    ldBlockInfoMap.clear();
        
    string header;
    string id;
    int chr, blockStart, blockEnd, snpNum;
    int idx = 0;
    int snpCount =  1;
    string snpName;
    string startSnpID, endSnpID;
    LDBlockInfo *ldblock;
    getline(in, header);
    const bool sourceTileHeader = header.find("SourceBlockID") != string::npos && header.find("TileID") != string::npos;
    const bool phasedTileHeader = header.find("TilePhase") != string::npos || header.find("TileGroupSize") != string::npos;
    const bool overlapTileHeader = header.find("CoreStartSnpIdx") != string::npos || header.find("CoreEndSnpIdx") != string::npos;
    if (phasedTileHeader || overlapTileHeader) {
        throw("Error: overlapping and phased tiled eigen LD formats are no longer supported. Regenerate the LD reference as ordinary or core-context eigen data.");
    }
    const bool blockSpecificFile = block && infoFile == dirname + "/block" + to_string(block) + ".ldm.info";
    coreContextLdmEigen = false;
    string line;
    while (getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        if (!(iss >> id >> chr >> blockStart >> startSnpID >> blockEnd >> endSnpID >> snpNum)) continue;
        string idstr = id;
        ldblock = new LDBlockInfo(idx++, idstr, chr);
        ldblock->startSnpIdx = blockStart;
        ldblock->endSnpIdx   = blockEnd;
        ldblock->numSnpInBlock = snpNum;
        vector<string> extra;
        string token;
        while (iss >> token) extra.push_back(token);
        if (sourceTileHeader && extra.size() >= 2) {
            ldblock->sourceBlockID = extra[0];
            ldblock->sourceTileID = stoi(extra[1]);
            coreContextLdmEigen = true;
        }

        if (block && !blockSpecificFile) {
            const string blockStr = to_string(block);
            if (idstr != blockStr && ldblock->sourceBlockID != blockStr) {
                delete ldblock;
                --idx;
                continue;
            }
        }
                
        ldBlockInfoVec.push_back(ldblock);
        ldBlockInfoMap[idstr] = ldblock;

    }
    in.close();
    numLDBlocks = (unsigned) ldBlockInfoVec.size();
    coreContextLdBlocks.clear();
    cout << numLDBlocks << " LD Blocks to be included from [" + infoFile + "]." << endl;
}

void Data::readBlockLdmSnpInfoFile(const string &dirname, const unsigned block){
    string snpInfoFile;
    bool blockSpecificFile = false;
    if (block) {
        snpInfoFile = dirname + "/block" + to_string(block) + ".snp.info";
        ifstream in(snpInfoFile.c_str());
        if (!in) snpInfoFile = dirname + "/snp.info";
        else blockSpecificFile = true;
    } else {
        snpInfoFile = dirname + "/snp.info";
    }
    ifstream in(snpInfoFile.c_str());
    if (!in) throw ("Error: can not open the file [" + snpInfoFile + "] to read.");
    cout << "Reading LDM SNP info from file [" + snpInfoFile + "]." << endl;
    snpInfoVec.clear();
    snpInfoMap.clear();
    for (unsigned i=0; i<ldBlockInfoVec.size(); ++i) {
        ldBlockInfoVec[i]->snpNameVec.clear();
        ldBlockInfoVec[i]->snpInfoVec.clear();
    }
    
    map<string,int> ld2snpMap;
    string header;
    string id, allele1, allele2;
    int chr, physPos,ld_n;
    float genPos;
    float allele1Freq;
    int index;
    string blockID;
    LDBlockInfo *ldblock;
    getline(in, header);
    while (in >> chr >> id >> index >> genPos >> physPos >> allele1 >> allele2 >> allele1Freq >> ld_n >> blockID) {
        if (block && !blockSpecificFile) {
            const string blockStr = to_string(block);
            if (blockID != blockStr && blockID.find(blockStr + ".tile") != 0) continue;
        }
        
        string blockIDstr = blockID;

        map<string, LDBlockInfo*>::iterator blockIt = ldBlockInfoMap.find(blockIDstr);
        if (blockIt == ldBlockInfoMap.end()) {
            throw("Error: SNP " + id + " refers to LD block [" + blockIDstr + "] not found in ldm.info.");
        }
        ldblock = blockIt->second;

        SnpInfo *snp = NULL;
        map<string, SnpInfo*>::iterator snpIt = snpInfoMap.find(id);
        if (snpIt != snpInfoMap.end()) {
            throw ("Error: Duplicate SNP ID found: \"" + id + "\".");
        } else {
            snp = new SnpInfo(index, id, allele1, allele2, chr, genPos, physPos);
            snp->af = allele1Freq;
            snp->ld_n = ld_n;
            snpInfoVec.push_back(snp);
            chromosomes.insert(snp->chrom);
            snpInfoMap.insert(pair<string, SnpInfo*>(id, snp));
        }
        
//        cout << idx << " " << id << " " << blockID << endl;
        
        snp->af = allele1Freq;
        snp->ld_n = ld_n;
        snp->block = blockIDstr;
        
        if (ld2snpMap.insert(pair<string, int>(blockIDstr + "_" + id, index)).second == false) {
            throw ("Error: Duplicate LDBlock-SNP pair found: \"" + blockIDstr + "_" + id + "\".");
        } else{
            ldblock->snpNameVec.push_back(id);
            ldblock->snpInfoVec.push_back(snp);
        }

    }
    in.close();
    numSnps = (unsigned) snpInfoVec.size();

    if (coreContextLdmEigen && !block && numSnps) {
        vector<char> seenIndex(numSnps, 0);
        unsigned badIndexCount = 0;
        for (unsigned i = 0; i < snpInfoVec.size(); ++i) {
            SnpInfo *snp = snpInfoVec[i];
            if (!snp || snp->index < 0 || (unsigned)snp->index >= numSnps || seenIndex[(unsigned)snp->index]) {
                ++badIndexCount;
                continue;
            }
            seenIndex[(unsigned)snp->index] = 1;
        }
        if (badIndexCount) {
            throw("Error: core-context tiled LD SNP indices are not unique and sequential in [" + snpInfoFile + "]. "
                  "This can happen if block-specific tiled outputs were merged by an older beta version. "
                  "Please regenerate the tiled LD folder or rebuild merged snp.info/ldm.info with the current version.");
        }
    }
    
    for (unsigned i=0; i<numLDBlocks; ++i) {
        LDBlockInfo *b = ldBlockInfoVec[i];
        if (block) {
            // Filtered snp.info: snpInfoVec is only this block's SNPs; startSnpIdx/endSnpIdx in ldm.info are global indices.
            if (b->snpInfoVec.empty()) {
                throw ("Error: no SNPs read for LD block " + b->ID + "; check snp.info block IDs match ldm.info.");
            }
            b->startPos = b->snpInfoVec.front()->physPos;
            b->endPos = b->snpInfoVec.back()->physPos;
        } else {
            if ((unsigned) b->startSnpIdx >= snpInfoVec.size() || (unsigned) b->endSnpIdx >= snpInfoVec.size()) {
                throw ("Error: LD block " + b->ID + " startSnpIdx/endSnpIdx out of range for snp.info (size " + to_string(snpInfoVec.size()) + ").");
            }
            b->startPos = snpInfoVec[b->startSnpIdx]->physPos;
            b->endPos = snpInfoVec[b->endSnpIdx]->physPos;
        }
    }
    
    numKeptInds = ld_n;
    
    cout << numSnps << " SNPs to be included from [" + snpInfoFile + "]." << endl;
}

void Data::readBlockLdmBinaryAndMakeItSparse(const string &dirname, const unsigned block, const float chisqThreshold, const bool writeLdmTxt){
    struct stat sb;
    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        throw("Error: cannot find the folder [" + dirname + "]");
    }
    
    vector<int> numSnpInRegion;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        LDBlockInfo *block = ldBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }
        
    keptLdBlockInfoVec = ldBlockInfoVec;
    
    cout << "Making sparse matrices using a chisq threshold of " << chisqThreshold << endl;

#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numLDBlocks; i++){
        //if (block && i != block - 1) continue;
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[i];

        string infile = dirname + "/block" + blockInfo->ID + ".ldm.bin";
        FILE *fp = fopen(infile.c_str(), "rb");
        if(!fp){throw ("Error: can not open the file [" + infile + "] to read.");}

        string outBinfile = dirname + "/block" + blockInfo->ID + ".sparse.bin";
//        FILE *outbin = fopen(outBinfile.c_str(), "wb");

        string outTxtfile;
        ofstream outtxt;
        if (writeLdmTxt) {
            outTxtfile = dirname + "/block" + blockInfo->ID + ".sparse.txt";
//            outtxt.open(outTxtfile.c_str());
        }
        
        int32_t blockSize = blockInfo->numSnpInBlock;
        
//        MatrixXf ldm(blockSize, blockSize);
//        uint64_t nElements = (uint64_t)blockSize * (uint64_t)blockSize;
//                
//        if(fread(ldm.data(), sizeof(float), nElements, fp) != nElements){
//            cout << "fread(U.data(), sizeof(float), nElements, fp): " << fread(ldm.data(), sizeof(float), nElements, fp) << endl;
//            cout << "nEle: " << nElements << " ldm.size: " << ldm.size() <<  " ldm.col: " << ldm.cols() << " row: " << ldm.rows() << endl;
//            throw("In LD block " + blockInfo->ID + ",size error in " + outBinfile);
//            // cout << "Read " << svdLDfile << " error (U)" << endl;
//            // throw("read file error");
//        }
        
        MatrixXf ldm(blockSize, blockSize);
        readBlockLDmatrix(dirname, blockInfo->ID, blockSize, ldm);

        
//        cout << " set zero " << blockInfo->ID << endl;
        SparseMatrix<float> ldmSp = ldm.sparseView();
        ldmSp.makeCompressed();
        
//        cout << " computed block " << blockInfo->ID << endl;

        // output sparse LDM
        Gadget::writeSparseMatrixBinary(ldmSp, outBinfile);
        
//        cout << " output block " << blockInfo->ID << endl;

        if (writeLdmTxt) {
            Gadget::writeSparseMatrixToText(ldmSp, outTxtfile);
        }
        
        
        fclose(fp);
//        fclose(outbin);
//        if (writeLdmTxt) outtxt.close();

        if(!(i%1)) cout << " computed block " << blockInfo->ID << "\r" << flush;
        
        if (block) {
            cout << "Written the sparse matrix data for block LD matrix into file [" << outBinfile << "]." << endl;
            if (writeLdmTxt) cout << "Written the sparse matrix data for block LD matrix into file [" << outTxtfile << "]." << endl;
        }
    }

    if (!block) {
        cout << "Written the sparse matrix data for block LD matrix into file [" << dirname << "/block*.sparse.bin]." << endl;
        if (writeLdmTxt) cout << "Written the sparse matrix data for block LD matrix into file [" << dirname << "/block*.sparse.txt]." << endl;
    }

}

void Data::readBlockLdmBinaryAndDoEigenDecomposition(const string &dirname, const unsigned block, const float eigenCutoff, const bool writeLdmTxt){
    struct stat sb;
    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        throw("Error: cannot find the folder [" + dirname + "]");
    }
    
    vector<int> numSnpInRegion;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        LDBlockInfo *block = ldBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }
        
    keptLdBlockInfoVec = ldBlockInfoVec;
    
#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numLDBlocks; i++){
        //if (block && i != block - 1) continue;
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[i];

        string outBinfile = dirname + "/block" + blockInfo->ID + ".eigen.bin";
        FILE *outbin = fopen(outBinfile.c_str(), "wb");

        string outTxtfile;
        ofstream outtxt;
        if (writeLdmTxt) {
            outTxtfile = dirname + "/block" + blockInfo->ID + ".eigen.txt";
            outtxt.open(outTxtfile.c_str());
        }
        
        
//        MatrixXf ldm(blockSize, blockSize);
//        uint64_t nElements = (uint64_t)blockSize * (uint64_t)blockSize;
//                
//        if(fread(ldm.data(), sizeof(float), nElements, fp) != nElements){
//            cout << "fread(ldm.data(), sizeof(float), nElements, fp): " << fread(ldm.data(), sizeof(float), nElements, fp) << endl;
//            cout << "nEle: " << nElements << " ldm.size: " << ldm.size() <<  " ldm.col: " << ldm.cols() << " row: " << ldm.rows() << endl;
//            throw("In LD block " + blockInfo->ID + ",size error in " + outBinfile);
//            // cout << "Read " << svdLDfile << " error (U)" << endl;
//            // throw("read file error");
//        }
        
        MatrixXf ldm;
        readBlockLDmatrix(dirname, blockInfo->ID, blockInfo->numSnpInBlock, ldm);
                
        
        MatrixXf eigenVec;
        VectorXf eigenVal;
        float sumPosEigVal;  // sum of all positive eigenvalues
        // cout << "rval: " << rval << endl;
        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        eigenDecomposition(ldm, eigenCutoff, eigenVal, eigenVec, sumPosEigVal);
        // cout << "rval: " << rval.row(0) << endl;
        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
        // save svd matrix

        blockInfo->sumPosEigVal = sumPosEigVal;
        
        int32_t numEigenValue = eigenVal.size();
        int32_t numSnpInBlock = blockInfo->numSnpInBlock;

        //cout << " Generate Eigen decomposition result for LD block " << i << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << endl;
        
        // save summary
        // 1, the number of SNPs in the block
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, outbin);
        // 2, the number of eigenvalues at with the given cutoff
        fwrite(&numEigenValue, sizeof(int32_t), 1, outbin);
        // 3. sum of all the positive eigenvalues
        fwrite(&sumPosEigVal, sizeof(float), 1, outbin);
        // 4. eigenvalue cutoff based on the proportion of variance explained in LD
        fwrite(&eigenCutoff, sizeof(float), 1, outbin);
        // 5. the selected eigenvalues
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, outbin);
        // 6. the selected eigenvector;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, outbin);
        
        if (writeLdmTxt) {
            outtxt << "Block " << blockInfo->ID << endl;
            outtxt << "numSnps " << blockInfo->numSnpInBlock << endl;
            outtxt << "numEigenvalues " << numEigenValue << endl;
            outtxt << "SumPositiveEigenvalues " << sumPosEigVal << endl;
            outtxt << "EigenCutoff " << eigenCutoff << endl;
            outtxt << "Eigenvalues\n" << eigenVal.transpose() << endl;
            outtxt << "Eigenvectors\n" << eigenVec << endl;
            outtxt << endl;
        }
        
        fclose(outbin);
        if (writeLdmTxt) outtxt.close();

        if(!(i%1)) cout << " computed block " << blockInfo->ID << "\r" << flush;
        
        if (block) {
            cout << "Written the eigen data for block LD matrix into file [" << outBinfile << "]." << endl;
            if (writeLdmTxt) cout << "Written the eigen data for block LD matrix into file [" << outTxtfile << "]." << endl;
        }
    }

    if (!block) {
        cout << "Written the eigen data for block LD matrix into file [" << dirname << "/block*.eigen.bin]." << endl;
        if (writeLdmTxt) cout << "Written the eigen data for block LD matrix into file [" << dirname << "/block*.eigen.txt]." << endl;
    }

}

void Data::readSparseBlockLdmBinaryAndDoEigenDecomposition(const string &dirname, const unsigned block, const float eigenCutoff, const bool writeLdmTxt){
    struct stat sb;
    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        throw("Error: cannot find the folder [" + dirname + "]");
    }
    
    vector<int> numSnpInRegion;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        LDBlockInfo *block = ldBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }
        
    keptLdBlockInfoVec = ldBlockInfoVec;
    
#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numLDBlocks; i++){
        //if (block && i != block - 1) continue;
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[i];

        string infile = dirname + "/block" + blockInfo->ID + ".sparse.bin";
        FILE *fp = fopen(infile.c_str(), "rb");
        if(!fp){throw ("Error: can not open the file [" + infile + "] to read.");}

        string outBinfile = dirname + "/block" + blockInfo->ID + ".eigen.bin";
        FILE *outbin = fopen(outBinfile.c_str(), "wb");

        string outTxtfile;
        ofstream outtxt;
        if (writeLdmTxt) {
            outTxtfile = dirname + "/block" + blockInfo->ID + ".eigen.txt";
            outtxt.open(outTxtfile.c_str());
        }
        
        int32_t blockSize = blockInfo->numSnpInBlock;
        
//        MatrixXf ldm(blockSize, blockSize);
        uint64_t nElements = (uint64_t)blockSize * (uint64_t)blockSize;
                
        MatrixXf ldm = MatrixXf(Gadget::readSparseMatrixBinary(infile));
        
//        if(fread(ldm.data(), sizeof(float), nElements, fp) != nElements){
//            cout << "fread(U.data(), sizeof(float), nElements, fp): " << fread(ldm.data(), sizeof(float), nElements, fp) << endl;
//            cout << "nEle: " << nElements << " ldm.size: " << ldm.size() <<  " ldm.col: " << ldm.cols() << " row: " << ldm.rows() << endl;
//            throw("In LD block " + blockInfo->ID + ",size error in " + outBinfile);
//            // cout << "Read " << svdLDfile << " error (U)" << endl;
//            // throw("read file error");
//        }
        
        MatrixXf eigenVec;
        VectorXf eigenVal;
        float sumPosEigVal;  // sum of all positive eigenvalues
        // cout << "rval: " << rval << endl;
        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        eigenDecomposition(ldm, eigenCutoff, eigenVal, eigenVec, sumPosEigVal);
        // cout << "rval: " << rval.row(0) << endl;
        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
        // save svd matrix

        blockInfo->sumPosEigVal = sumPosEigVal;
        
        int32_t numEigenValue = eigenVal.size();
        int32_t numSnpInBlock = blockSize;

        //cout << " Generate Eigen decomposition result for LD block " << i << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << endl;
        
        // save summary
        // 1, the number of SNPs in the block
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, outbin);
        // 2, the number of eigenvalues at with the given cutoff
        fwrite(&numEigenValue, sizeof(int32_t), 1, outbin);
        // 3. sum of all the positive eigenvalues
        fwrite(&sumPosEigVal, sizeof(float), 1, outbin);
        // 4. eigenvalue cutoff based on the proportion of variance explained in LD
        fwrite(&eigenCutoff, sizeof(float), 1, outbin);
        // 5. the selected eigenvalues
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, outbin);
        // 6. the selected eigenvector;
        nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, outbin);
        
        if (writeLdmTxt) {
            outtxt << "Block " << blockInfo->ID << endl;
            outtxt << "numSnps " << numSnpInBlock << endl;
            outtxt << "numEigenvalues " << numEigenValue << endl;
            outtxt << "SumPositiveEigenvalues " << sumPosEigVal << endl;
            outtxt << "EigenCutoff " << eigenCutoff << endl;
            outtxt << "Eigenvalues\n" << eigenVal.transpose() << endl;
            outtxt << "Eigenvectors\n" << eigenVec << endl;
            outtxt << endl;
        }
        
        fclose(fp);
        fclose(outbin);
        if (writeLdmTxt) outtxt.close();

        if(!(i%1)) cout << " computed block " << blockInfo->ID << "\r" << flush;
        
        if (block) {
            cout << "Written the eigen data for block LD matrix into file [" << outBinfile << "]." << endl;
            if (writeLdmTxt) cout << "Written the eigen data for block LD matrix into file [" << outTxtfile << "]." << endl;
        }
    }

    if (!block) {
        cout << "Written the eigen data for block LD matrix into file [" << dirname << "/block*.eigen.bin]." << endl;
        if (writeLdmTxt) cout << "Written the eigen data for block LD matrix into file [" << dirname << "/block*.eigen.txt]." << endl;
    }

}

void Data::readEigenMatrixBinaryFile(const string &dirname, const float eigenCutoff, const bool writeLdmTxt, const string &outputDir, const int quantizedBits, const bool q8Entropy){
    if (!Gadget::directoryExist(dirname)) {
        throwEigenReadErr("Error: cannot find the folder [" + dirname + "]");
    }
    
    vector<int>numSnpInRegion;
    LDBlockInfo * block;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        block = ldBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }
    eigenValLdBlock.resize(numLDBlocks);
    eigenVecLdBlock.resize(numLDBlocks);
    cout << "Reading per-block eigen binaries under [" << dirname << "] (" << numLDBlocks << " blocks, serial) ..." << endl;
    // Serial loop: exceptions must not be thrown from OpenMP parallel regions (undefined behavior with libgomp;
    // often terminates with std::string and no message on stderr).
    for(int i = 0; i < (int)numLDBlocks; i++){
        LDBlockInfo * block;
        block = ldBlockInfoVec[i];
        
        if (!block->kept) continue;
        
        int32_t cur_m = 0;
        int32_t cur_k = 0;
        float sumPosEigVal = 0;
        float oldEigenCutoff =0;
        
        string infile;
        FILE *fp = NULL;
        if (!openEigenBlockPayloadForRead(dirname, block->ID, quantizedBits, q8Entropy, fp, infile))
            throwEigenReadErr("Error: can not open eigen data for block [" + block->ID + "] in [" + dirname + "] to read.");

        VectorXf lambda;
        MatrixXf U;
        readEigenHeaderOrThrow(fp, infile, block->ID, numSnpInRegion[i], cur_m, cur_k, sumPosEigVal, oldEigenCutoff, lambda);
        readEigenUOrThrow(fp, infile, block->ID, cur_m, cur_k, quantizedBits, q8Entropy, U);
        
        fclose(fp);
        
        bool haveValue = false;
        int revIdx = 0;
        if(oldEigenCutoff < eigenCutoff & i == 0){
            cout << "Warning: current proportion of variance in LD block is set as " + to_string(eigenCutoff)+ ". But the proportion of variance is set as "<< to_string(oldEigenCutoff) + " in "  + infile + ".\n";
            // throw("");
        }
        // cout << "lambda: " << lambda << endl;
        // cout << "U: " << U << endl;
        // eigenVecLdBlock[i] = U;
        // eigenValLdBlock[i] = lambda;
        
        if (eigenCutoff < oldEigenCutoff) {
            truncateEigenMatrix(sumPosEigVal, eigenCutoff, lambda, U, eigenValLdBlock[i], eigenVecLdBlock[i]);
        } else {
            eigenValLdBlock[i] = lambda;
            eigenVecLdBlock[i] = U;
        }
        block->sumPosEigVal = sumPosEigVal;
        block->eigenvalues = lambda;
        
        string outTxtfile;
        ofstream outtxt;
        if (writeLdmTxt) {
            Gadget::createDirectory(outputDir);
            outTxtfile = outputDir + "/block" + block->ID + ".eigen.txt";
            outtxt.open(outTxtfile.c_str());
            outtxt << "Block " << block->ID << endl;
            outtxt << "numSnps " << cur_m << endl;
            outtxt << "numEigenvalues " << cur_k << endl;
            outtxt << "SumPositiveEigenvalues " << sumPosEigVal << endl;
            outtxt << "EigenCutoff " << eigenCutoff << endl;
            outtxt << "Eigenvalues\n" << lambda.transpose() << endl;
            outtxt << "Eigenvectors\n" << U << endl;
            outtxt.close();
            cout << "Output eigen matrix for block " << block->ID << " in file [" << outTxtfile << "]." << endl;
        }
    }
    
    lowRankModel = true;
}

bool Data::blockLdmBinFileExists(const string &dirname, const string &blockID) const {
    string path = dirname + "/block" + blockID + ".ldm.bin";
    struct stat sb;
    return stat(path.c_str(), &sb) == 0 && S_ISREG(sb.st_mode);
}

bool Data::blockEigenBinFileExists(const string &dirname, const string &blockID) const {
    string path = dirname + "/block" + blockID + ".eigen.bin";
    struct stat sb;
    if (stat(path.c_str(), &sb) == 0 && S_ISREG(sb.st_mode)) return true;
    string sourceBlockID;
    if (!sourceBlockIDFromTileID(blockID, sourceBlockID)) return false;
    path = tiledEigenBundleFilename(dirname, sourceBlockID, 0);
    return stat(path.c_str(), &sb) == 0 && S_ISREG(sb.st_mode);
}

bool Data::readBlockEigenBinContents(
    const string &dirname,
    const string &blockID,
    int32_t expect_m,
    int32_t &cur_m,
    int32_t &cur_k,
    float &sumPosEigVal,
    float &oldEigenCutoff,
    VectorXf &lambda,
    MatrixXf &U
) {
    string infile;
    FILE *fp = NULL;
    if (!openEigenBlockPayloadForRead(dirname, blockID, 0, false, fp, infile))
        return false;
    if (fread(&cur_m, sizeof(int32_t), 1, fp) != 1) {
        fclose(fp);
        return false;
    }
    if (cur_m != expect_m) {
        fclose(fp);
        return false;
    }
    if (fread(&cur_k, sizeof(int32_t), 1, fp) != 1) {
        fclose(fp);
        return false;
    }
    if (fread(&sumPosEigVal, sizeof(float), 1, fp) != 1) {
        fclose(fp);
        return false;
    }
    if (fread(&oldEigenCutoff, sizeof(float), 1, fp) != 1) {
        fclose(fp);
        return false;
    }
    lambda.resize(cur_k);
    if (fread(lambda.data(), sizeof(float), cur_k, fp) != (size_t)cur_k) {
        fclose(fp);
        return false;
    }
    U.resize(cur_m, cur_k);
    uint64_t nElements = (uint64_t)cur_m * (uint64_t)cur_k;
    if (fread(U.data(), sizeof(float), nElements, fp) != nElements) {
        fclose(fp);
        return false;
    }
    fclose(fp);
    return true;
}

bool Data::buildSubmatrixLdFromEigenFactors(
    const float eigenCutoff,
    float sumPosEigVal,
    float oldEigenCutoff,
    const VectorXf &lambda,
    const MatrixXf &U,
    const vector<unsigned> &kept,
    MatrixXf &RssOut
) {
    unsigned m2 = (unsigned)kept.size();
    if (m2 == 0) return false;

    VectorXf lambdaUse;
    MatrixXf Uuse;
    if (eigenCutoff < oldEigenCutoff) {
        truncateEigenMatrix(sumPosEigVal, eigenCutoff, lambda, U, lambdaUse, Uuse);
    } else {
        lambdaUse = lambda;
        Uuse = U;
    }
    int k2 = (int)lambdaUse.size();
    if (k2 <= 0) return false;

    MatrixXf U_K((int)m2, k2);
    for (unsigned a = 0; a < m2; ++a)
        U_K.row((int)a) = Uuse.row((int)kept[a]);

    RssOut = U_K * lambdaUse.asDiagonal() * U_K.transpose();
    normaliseDenseCorrelation(RssOut);
    return true;
}

bool Data::trySubLdEigenFromFullLdm(
    LDBlockInfo *block,
    const string &dirname,
    const float eigenCutoff,
    const VectorXf &gwasBlock,
    VectorXf &eigenValOut,
    MatrixXf &eigenVecOut,
    VectorXf &wcorrOut,
    MatrixXf &Qout,
    float &sumPosEigValOut,
    VectorXi &remapOut,
    vector<unsigned> &keptLocalOut,
    bool *outUsedEigenRecon,
    const int quantizedBits,
    const bool q8Entropy
) {
    if (outUsedEigenRecon) *outUsedEigenRecon = false;

    int m = block->numSnpInBlock;
    if (m <= 1) return false;
    if (!recomputeEigen) return false;
    unsigned n_skip = 0;
    for (int j = 0; j < m; ++j)
        if (block->snpInfoVec[j]->skip) ++n_skip;
    if (n_skip == 0 || n_skip >= (unsigned)m) return false;

    vector<unsigned> kept;
    kept.reserve(m - n_skip);
    for (int j = 0; j < m; ++j)
        if (!block->snpInfoVec[j]->skip) kept.push_back((unsigned)j);

    unsigned m2 = (unsigned)kept.size();
    if (m2 == 0) return false;

    MatrixXf Rss(m2, m2);
    bool haveRss = false;

    if (blockLdmBinFileExists(dirname, block->ID)) {
        try {
            MatrixXf L;
            readBlockLDmatrix(dirname, block->ID, m, L);
            for (unsigned a = 0; a < m2; ++a)
                for (unsigned b = 0; b < m2; ++b)
                    Rss(a, b) = L((int)kept[a], (int)kept[b]);
            haveRss = true;
        } catch (...) {
            haveRss = false;
        }
    }

    if (!haveRss) {
        int32_t cur_m = 0, cur_k = 0;
        float sumPosEigVal = 0, oldEigenCutoff = 0;
        VectorXf lambda;
        MatrixXf U;
        if (readBlockEigenBinContents(dirname, block->ID, m, cur_m, cur_k, sumPosEigVal, oldEigenCutoff, lambda, U)) {
            if (!buildSubmatrixLdFromEigenFactors(eigenCutoff, sumPosEigVal, oldEigenCutoff, lambda, U, kept, Rss))
                return false;
            if (outUsedEigenRecon) *outUsedEigenRecon = true;
        } else if (quantizedBits > 0) {
            string infile;
            FILE *fp = NULL;
            if (!openEigenBlockPayloadForRead(dirname, block->ID, quantizedBits, q8Entropy, fp, infile))
                return false;
            QuantizedEigenBlock qb;
            try {
                readEigenHeaderOrThrow(fp, infile, block->ID, m, cur_m, cur_k, sumPosEigVal, oldEigenCutoff, lambda);
                readQuantizedEigenUOrThrow(fp, infile, block->ID, cur_m, cur_k, quantizedBits, q8Entropy, lambda, qb);
            } catch (...) {
                fclose(fp);
                return false;
            }
            fclose(fp);
            if (eigenCutoff + 1e-6f < oldEigenCutoff) {
                const int startCol = eigenTailStartIndex(sumPosEigVal, eigenCutoff, lambda);
                keepQuantizedEigenTail(qb, startCol);
            }
            Rss.resize(m2, m2);
            Rss.setZero();
            for (unsigned a = 0; a < m2; ++a) {
                for (unsigned b = 0; b <= a; ++b) {
                    float val = 0.0f;
                    for (int eigIdx = 0; eigIdx < qb.k; ++eigIdx)
                        val += qb.qComponent((int)kept[a], eigIdx) * qb.qComponent((int)kept[b], eigIdx);
                    Rss((int)a, (int)b) = val;
                    Rss((int)b, (int)a) = val;
                }
            }
            normaliseDenseCorrelation(Rss);
            if (outUsedEigenRecon) *outUsedEigenRecon = true;
        } else {
            return false;
        }
    }

    eigenDecomposition(Rss, eigenCutoff, eigenValOut, eigenVecOut, sumPosEigValOut);

    remapOut.resize(m);
    int lc = 0;
    for (int j = 0; j < m; ++j) {
        if (block->snpInfoVec[j]->skip) remapOut[j] = -1;
        else remapOut[j] = lc++;
    }
    keptLocalOut.swap(kept);

    VectorXf sqrtLambda = eigenValOut.array().sqrt();
    VectorXf b_kept(m2);
    for (unsigned a = 0; a < m2; ++a)
        b_kept[a] = gwasBlock[(int)keptLocalOut[a]];
    wcorrOut = (1.0f / sqrtLambda.array()).matrix().asDiagonal() * (eigenVecOut.transpose() * b_kept);
    Qout = sqrtLambda.asDiagonal() * eigenVecOut.transpose();
    return true;
}

bool Data::trySubLdEigenFromFullLdmBivariate(
    LDBlockInfo *block,
    const string &dirname,
    const float eigenCutoff,
    const VectorXf &gwas1,
    const VectorXf &gwas2,
    VectorXf &wcorrOut,
    MatrixXf &Qout,
    VectorXf &eigenValOut,
    MatrixXf &eigenVecOut,
    float &sumPosEigValOut,
    VectorXi &remapOut,
    vector<unsigned> &keptLocalOut,
    bool *outUsedEigenRecon
) {
    if (outUsedEigenRecon) *outUsedEigenRecon = false;

    int m = block->numSnpInBlock;
    if (m <= 1) return false;
    if (!recomputeEigen) return false;
    unsigned n_skip = 0;
    for (int j = 0; j < m; ++j)
        if (block->snpInfoVec[j]->skip) ++n_skip;
    if (n_skip == 0 || n_skip >= (unsigned)m) return false;

    vector<unsigned> kept;
    kept.reserve(m - n_skip);
    for (int j = 0; j < m; ++j)
        if (!block->snpInfoVec[j]->skip) kept.push_back((unsigned)j);

    unsigned m2 = (unsigned)kept.size();
    if (m2 == 0) return false;

    MatrixXf Rss(m2, m2);
    bool haveRss = false;

    if (blockLdmBinFileExists(dirname, block->ID)) {
        try {
            MatrixXf L;
            readBlockLDmatrix(dirname, block->ID, m, L);
            for (unsigned a = 0; a < m2; ++a)
                for (unsigned b = 0; b < m2; ++b)
                    Rss(a, b) = L((int)kept[a], (int)kept[b]);
            haveRss = true;
        } catch (...) {
            haveRss = false;
        }
    }

    if (!haveRss) {
        int32_t cur_m = 0, cur_k = 0;
        float sumPosEigVal = 0, oldEigenCutoff = 0;
        VectorXf lambda;
        MatrixXf U;
        if (!readBlockEigenBinContents(dirname, block->ID, m, cur_m, cur_k, sumPosEigVal, oldEigenCutoff, lambda, U))
            return false;
        if (!buildSubmatrixLdFromEigenFactors(eigenCutoff, sumPosEigVal, oldEigenCutoff, lambda, U, kept, Rss))
            return false;
        if (outUsedEigenRecon) *outUsedEigenRecon = true;
    }

    eigenDecomposition(Rss, eigenCutoff, eigenValOut, eigenVecOut, sumPosEigValOut);

    remapOut.resize(m);
    int lc = 0;
    for (int j = 0; j < m; ++j) {
        if (block->snpInfoVec[j]->skip) remapOut[j] = -1;
        else remapOut[j] = lc++;
    }
    keptLocalOut.swap(kept);

    VectorXf sqrtLambda = eigenValOut.array().sqrt();
    VectorXf b1_kept(m2), b2_kept(m2);
    for (unsigned a = 0; a < m2; ++a) {
        b1_kept[a] = gwas1[(int)keptLocalOut[a]];
        b2_kept[a] = gwas2[(int)keptLocalOut[a]];
    }
    VectorXf w1 = (1.0f / sqrtLambda.array()).matrix().asDiagonal() * (eigenVecOut.transpose() * b1_kept);
    VectorXf w2 = (1.0f / sqrtLambda.array()).matrix().asDiagonal() * (eigenVecOut.transpose() * b2_kept);
    wcorrOut.resize(w1.size() + w2.size());
    wcorrOut.head(w1.size()) = w1;
    wcorrOut.tail(w2.size()) = w2;
    Qout = sqrtLambda.asDiagonal() * eigenVecOut.transpose();
    return true;
}

void Data::readEigenMatrixBinaryFileAndMakeWandQ(const string &dirname, const float eigenCutoff, const vector<VectorXf> &GWASeffects, const VectorXf &nGWASblock, const bool noscale, const bool makePseudoSummary, const int quantizedBits, const bool q8Entropy){
    const bool firstEigenRead = !wcorrBlocks.size();
    Gadget::Timer eigenReadTimer;
    if (firstEigenRead) {  // only print for the first time reading the data
        cout << "Reading eigenvectors from binary file and making W and Q matrices..." << endl;
    }
    if (!Gadget::directoryExist(dirname)) {
        throw("Error: cannot find the folder [" + dirname + "]");
    }

    precomputedMarginalizedContext = integrateOutContext && coreContextLdmEigen;
    if (precomputedMarginalizedContext) {
        set<string> checkedSources;
        for (unsigned i = 0; i < keptLdBlockInfoVec.size(); ++i) {
            const string &sourceBlockID = keptLdBlockInfoVec[i]->sourceBlockID;
            if (sourceBlockID.empty()
                || (checkedSources.insert(sourceBlockID).second
                    && !fileExistsAndNonEmptyLocal(marginalizedContextBundleFilename(dirname, sourceBlockID)))) {
                precomputedMarginalizedContext = false;
                break;
            }
        }
        if (precomputedMarginalizedContext) {
            cout << "Reading precomputed marginalized-context tile matrices; sparse cross-tile LD links will not be loaded." << endl;
        } else {
            cout << "Precomputed marginalized-context bundles were not found; using the legacy sparse-link preparation path." << endl;
        }
    }
    marginalizedContextPrepared = false;
    marginalizedContextMaxRatio = 0.0f;
    marginalizedContextTotalLinks = 0;
    marginalizedContextTotalFlankSnps = 0;
    
    vector<int>numSnpInRegion(numKeptLDBlocks);

    if (GWASeffects.size() != numKeptLDBlocks || nGWASblock.size() != numKeptLDBlocks) {
        throw("Error: incomplete block-level GWAS data while constructing low-rank W/Q matrices (expected "
              + to_string(numKeptLDBlocks) + " blocks, found "
              + to_string(GWASeffects.size()) + " effect vectors and "
              + to_string((long long)nGWASblock.size()) + " sample-size values).");
    }
    
    for(int i = 0; i < numKeptLDBlocks;i++){
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
        if (GWASeffects[i].size() != block->numSnpInBlock) {
            throw("Error: GWAS-effect vector for LD block " + block->ID + " has "
                  + to_string((long long)GWASeffects[i].size()) + " entries; expected "
                  + to_string(block->numSnpInBlock) + ".");
        }
    }

    unsigned subLdPlanned = 0;
    vector<string> subLdMissingLdBlockIds;
    for (int i = 0; i < numKeptLDBlocks; ++i) {
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        int m = block->numSnpInBlock;
        if (m <= 1) continue;
        if (!recomputeEigen) continue;
        unsigned n_skip = 0;
        for (int j = 0; j < m; ++j)
            if (block->snpInfoVec[j]->skip) ++n_skip;
        if (n_skip == 0 || n_skip >= (unsigned)m) continue;
        bool haveQuantEigen = false;
        if (quantizedBits > 0) {
            string qpath;
            FILE *qfp = NULL;
            haveQuantEigen = openEigenBlockPayloadForRead(dirname, block->ID, quantizedBits, q8Entropy, qfp, qpath);
            if (qfp) fclose(qfp);
        }
        if (blockLdmBinFileExists(dirname, block->ID) || blockEigenBinFileExists(dirname, block->ID) || haveQuantEigen)
            ++subLdPlanned;
        else
            subLdMissingLdBlockIds.push_back(block->ID);
    }
    if (!subLdMissingLdBlockIds.empty() && !wcorrBlocks.size()) {
        cout << "Warning: --recompute-eigen with --skip: partial skip in " << subLdMissingLdBlockIds.size()
             << " LD block(s); need block<ID>.ldm.bin, block<ID>.eigen.bin, or matching block<ID>.eigen.q*.bin in [" << dirname << "]. Missing for block ID(s):";
        for (size_t t = 0; t < subLdMissingLdBlockIds.size(); ++t)
            cout << " " << subLdMissingLdBlockIds[t];
        cout << ". Falling back to full-block block*.eigen.bin for those blocks." << endl;
    }
    if (subLdPlanned && !wcorrBlocks.size()) {
        cout << "For " << subLdPlanned << " LD block(s), the LD submatrix of retained SNPs will be reconstructed (--recompute-eigen) and the GWAS effects of skipped SNPs will not be used." << endl;
    }

    eigenValLdBlock.resize(numLDBlocks);
    eigenVecLdBlock.resize(numLDBlocks);
    wcorrBlocks.resize(numKeptLDBlocks);
    numSnpsBlock.resize(numKeptLDBlocks);
    numEigenvalBlock.resize(numKeptLDBlocks);
    Qblocks.resize(numKeptLDBlocks);
    quantizedEigenBlocks.assign(numKeptLDBlocks, QuantizedEigenBlock());
    marginalizedContextEigenvalues.assign(numKeptLDBlocks, VectorXf());
    marginalizedContextOriginalNormWeights.assign(numKeptLDBlocks, VectorXf());
    marginalizedContextLinkCounts.assign(numKeptLDBlocks, 0);
    marginalizedContextFlankCounts.assign(numKeptLDBlocks, 0);

    EigenNbThreadsGuard eigenNbThreadsGuard;

    VectorXf n_trn(numKeptLDBlocks);
    VectorXf n_val(numKeptLDBlocks);
    if (makePseudoSummary) {
        pseudoGwasNtrnBlock.resize(numKeptLDBlocks);
        pseudoGwasNValBlock.resize(numKeptLDBlocks);
        for(int i = 0; i < numKeptLDBlocks;i++){
            n_trn[i] = 0.9f * nGWASblock[i];
            n_val[i] = nGWASblock[i] - n_trn[i];
            pseudoGwasNtrnBlock[i] = n_trn[i];
            pseudoGwasNValBlock[i] = n_val[i];
        }

        pseudoGwasEffectTrn.resize(numKeptLDBlocks);
        pseudoGwasEffectVal.resize(numKeptLDBlocks);
        b_val.setZero(numIncdSnps);
    }
    
#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        block->eigenColRemap.resize(0);
        block->subLdKeptLocalIdx.clear();
        block->qDiag.resize(0);

        if (precomputedMarginalizedContext) {
            string infile;
            FILE *fp = NULL;
            if (!openMarginalizedContextPayloadForRead(dirname, block->ID, fp, infile))
                throw("Error: can not open precomputed marginalized-context data for tile " + block->ID + ".");

            VectorXf lambda;
            MatrixXf rotatedWOperator;
            MatrixXf rotatedQ;
            uint64_t numLinks = 0;
            unsigned numFlankSnps = 0;
            readMarginalizedContextPayloadOrThrow(fp, infile, block->ID, block->numSnpInBlock,
                                                   lambda, marginalizedContextEigenvalues[i],
                                                   rotatedWOperator, rotatedQ, numLinks, numFlankSnps);
            fclose(fp);

            eigenValLdBlock[i] = lambda;
            block->eigenvalues = lambda;
            block->sumPosEigVal = lambda.sum();
            wcorrBlocks[i].noalias() = rotatedWOperator * (rotatedQ * GWASeffects[i]);
            Qblocks[i].swap(rotatedQ);
            eigenVecLdBlock[i].resize(0, 0);
            quantizedEigenBlocks[i] = QuantizedEigenBlock();
            numSnpsBlock[i] = Qblocks[i].cols();
            numEigenvalBlock[i] = Qblocks[i].rows();
            marginalizedContextLinkCounts[i] = numLinks;
            marginalizedContextFlankCounts[i] = numFlankSnps;

            if (makePseudoSummary) {
                Stat::ParallelTaskRngScope pseudoRng((unsigned)i, Stat::pseudoSummaryRngSalt);
                VectorXf rnd(Qblocks[i].rows());
                for (int j = 0; j < rnd.size(); ++j) rnd[j] = Stat::snorm();
                pseudoGwasEffectTrn[i] = gwasEffectInBlock[i]
                    + sqrt(1.0f/n_trn[i] - 1.0f/nGWASblock[i]) * Qblocks[i].transpose() * rnd;
                pseudoGwasEffectVal[i] = nGWASblock[i]/n_val[i] * gwasEffectInBlock[i]
                    - n_trn[i]/n_val[i] * pseudoGwasEffectTrn[i];
                b_val.segment(block->startSnpIdx, block->numSnpInBlock) = pseudoGwasEffectVal[i];
            }
            continue;
        }

        VectorXf eigenValSub;
        MatrixXf eigenVecSub;
        VectorXf wcorrSub;
        MatrixXf Qsub;
        float sumPosSub = 0;
        VectorXi remapSub;
        vector<unsigned> keptSub;

        bool didSubLd = trySubLdEigenFromFullLdm(block, dirname, eigenCutoff, GWASeffects[i], eigenValSub, eigenVecSub, wcorrSub, Qsub, sumPosSub, remapSub, keptSub, nullptr, quantizedBits, q8Entropy);

        if (didSubLd) {
#pragma omp critical
            {
                cout << "LD block " << block->ID << ": " << block->numSnpInBlock - (int)keptSub.size()
                     << " SNPs skipped, " << keptSub.size() << " retained." << endl;
            }
            eigenValLdBlock[i] = eigenValSub;
            wcorrBlocks[i] = wcorrSub;
            block->sumPosEigVal = sumPosSub;
            block->eigenvalues = eigenValSub;
            block->eigenColRemap = remapSub;
            block->subLdKeptLocalIdx = keptSub;
            if (quantizedBits > 0) {
                quantizeEigenMatrixToBlock(eigenVecSub, eigenValSub, quantizedBits, quantizedEigenBlocks[i]);
                eigenVecLdBlock[i].resize(0, 0);
                Qblocks[i].resize(0, 0);
                numSnpsBlock[i] = quantizedEigenBlocks[i].m;
                numEigenvalBlock[i] = quantizedEigenBlocks[i].k;
            } else {
                eigenVecLdBlock[i] = eigenVecSub;
                Qblocks[i] = Qsub;
                quantizedEigenBlocks[i] = QuantizedEigenBlock();
                numSnpsBlock[i] = Qblocks[i].cols();
                numEigenvalBlock[i] = Qblocks[i].rows();
            }

            if (makePseudoSummary) {
                Stat::ParallelTaskRngScope pseudoRng((unsigned)i, Stat::pseudoSummaryRngSalt);
                VectorXf rnd(eigenValLdBlock[i].size());
                for (int j = 0; j < rnd.size(); ++j) rnd[j] = Stat::snorm();
                VectorXf noise;
                if (quantizedBits > 0)
                    noise = quantizedEigenBlocks[i].transposeMultiply(rnd);
                else
                    noise = eigenVecSub
                        * (eigenValLdBlock[i].array().sqrt().matrix().asDiagonal() * rnd);
                pseudoGwasEffectTrn[i] = gwasEffectInBlock[i]
                    + sqrt(1.0f/n_trn[i] - 1.0f/nGWASblock[i]) * noise;
                pseudoGwasEffectVal[i] = nGWASblock[i]/n_val[i] * gwasEffectInBlock[i]
                    - n_trn[i]/n_val[i] * pseudoGwasEffectTrn[i];
                b_val.segment(block->startSnpIdx, block->numSnpInBlock) = pseudoGwasEffectVal[i];
            }

            eigenVecLdBlock[i].resize(0,0);
            continue;
        }

        int32_t cur_m = 0;
        int32_t cur_k = 0;
        float sumPosEigVal = 0;
        float oldEigenCutoff =0;
        
        string infile;
        FILE *fp = NULL;
        if (!openEigenBlockPayloadForRead(dirname, block->ID, quantizedBits, q8Entropy, fp, infile)) {
#pragma omp critical
            {
                cout << "Error: can not open eigen data for block [" + block->ID + "] in [" + dirname + "] to read." << endl;
            }
            throw("Error: can not open eigen data for block [" + block->ID + "] in [" + dirname + "] to read.");
        }

        VectorXf lambda;
        MatrixXf U;
        readEigenHeaderOrThrow(fp, infile, block->ID, numSnpInRegion[i], cur_m, cur_k, sumPosEigVal, oldEigenCutoff, lambda);
        bool haveValue = false;
        int revIdx = 0;
        if(oldEigenCutoff + 1e-6f < eigenCutoff & i == 0){
            cout << "Warning: current proportion of variance in LD block is set as " + to_string(eigenCutoff)+ ". But the proportion of variance is set as "<< to_string(oldEigenCutoff) + " in "  + infile + ".\n";
            // throw("");
        }

        const bool keepQuantizedInMemory = quantizedBits > 0;
        if (keepQuantizedInMemory) {
            if (i == 0 && firstEigenRead)
                cout << "Keeping quantized eigenvectors in memory; Q columns will be de-quantized on demand during sampling." << endl;
            readQuantizedEigenUOrThrow(fp, infile, block->ID, cur_m, cur_k, quantizedBits, q8Entropy, lambda, quantizedEigenBlocks[i]);
            fclose(fp);   // this is important!

            if (eigenCutoff + 1e-6f < oldEigenCutoff) {
                const int startCol = eigenTailStartIndex(sumPosEigVal, eigenCutoff, lambda);
                keepQuantizedEigenTail(quantizedEigenBlocks[i], startCol);
                lambda = lambda.tail(quantizedEigenBlocks[i].k).eval();
            }

            eigenValLdBlock[i] = lambda;
            eigenVecLdBlock[i].resize(0, 0);
            block->sumPosEigVal = sumPosEigVal;
            block->eigenvalues = lambda;
            wcorrBlocks[i] = makeWcorrFromQuantizedEigen(quantizedEigenBlocks[i], GWASeffects[i]);
            Qblocks[i].resize(0, 0);
            numSnpsBlock[i] = cur_m;
            numEigenvalBlock[i] = quantizedEigenBlocks[i].k;

            if (makePseudoSummary) {
                Stat::ParallelTaskRngScope pseudoRng((unsigned)i, Stat::pseudoSummaryRngSalt);
                VectorXf rnd(eigenValLdBlock[i].size());
                for (int j = 0; j < rnd.size(); ++j) rnd[j] = Stat::snorm();
                VectorXf noise = quantizedEigenBlocks[i].transposeMultiply(rnd);
                pseudoGwasEffectTrn[i] = gwasEffectInBlock[i]
                    + sqrt(1.0f/n_trn[i] - 1.0f/nGWASblock[i]) * noise;
                pseudoGwasEffectVal[i] = nGWASblock[i]/n_val[i] * gwasEffectInBlock[i]
                    - n_trn[i]/n_val[i] * pseudoGwasEffectTrn[i];
                b_val.segment(block->startSnpIdx, block->numSnpInBlock) = pseudoGwasEffectVal[i];
            }

            continue;
        }

        readEigenUOrThrow(fp, infile, block->ID, cur_m, cur_k, quantizedBits, q8Entropy, U);
        
        fclose(fp);   // this is important!

        // cout << "lambda: " << lambda << endl;
        // cout << "U: " << U << endl;
        // eigenVecLdBlock[i] = U;
        // eigenValLdBlock[i] = lambda;
        
        if (eigenCutoff + 1e-6f < oldEigenCutoff) {
            truncateEigenMatrix(sumPosEigVal, eigenCutoff, lambda, U, eigenValLdBlock[i], eigenVecLdBlock[i]);
        } else {
            eigenValLdBlock[i] = lambda;
            eigenVecLdBlock[i] = U;
        }
        quantizedEigenBlocks[i] = QuantizedEigenBlock();
        block->sumPosEigVal = sumPosEigVal;
        block->eigenvalues = lambda;
                        
        // make w and Q
        VectorXf sqrtLambda = eigenValLdBlock[i].array().sqrt();
        wcorrBlocks[i] = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (eigenVecLdBlock[i].transpose() * GWASeffects[i] );
        //MatrixXf tmpQblocks = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
        //MatrixDat matrixDat = MatrixDat(block->snpNameVec, tmpQblocks);
        Qblocks[i] = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
        
        // if (noscale) {
        //     VectorXf Dsqrt(block->numSnpInBlock);
        //     for (unsigned j=0; j<block->numSnpInBlock; ++j) {
        //         SnpInfo *snp = block->snpInfoVec[j];
        //         Dsqrt[j] = sqrt(snp->twopq)*snp->gwas_scalar;
        //         //Dsqrt[j] = snp->gwas_scalar;
        //     }
        //     Qblocks[i] = Qblocks[i] * Dsqrt.asDiagonal();
        // }

        numSnpsBlock[i] = Qblocks[i].cols();
        numEigenvalBlock[i] = Qblocks[i].rows();

        if (makePseudoSummary) {
            Stat::ParallelTaskRngScope pseudoRng((unsigned)i, Stat::pseudoSummaryRngSalt);
            VectorXf rnd(eigenValLdBlock[i].size());
            for (int j = 0; j < rnd.size(); ++j) rnd[j] = Stat::snorm();
            pseudoGwasEffectTrn[i] = gwasEffectInBlock[i]
                + sqrt(1.0f/n_trn[i] - 1.0f/nGWASblock[i]) * eigenVecLdBlock[i]
                * (eigenValLdBlock[i].array().sqrt().matrix().asDiagonal() * rnd);
            pseudoGwasEffectVal[i] = nGWASblock[i]/n_val[i] * gwasEffectInBlock[i]
                - n_trn[i]/n_val[i] * pseudoGwasEffectTrn[i];
            b_val.segment(block->startSnpIdx, block->numSnpInBlock) = pseudoGwasEffectVal[i];
        }
                
        eigenVecLdBlock[i].resize(0,0);
    }

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < numKeptLDBlocks; ++i) {
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        if (Qblocks[i].size())
            cacheRealizedQdiag(block, Qblocks[i]);
        else if (quantizedEigenBlocks[i].active())
            cacheRealizedQdiag(block, quantizedEigenBlocks[i]);
        else
            throw("Error: missing Q matrix while caching the realized Q'Q diagonal for LD block " + block->ID + ".");
    }

    if (firstEigenRead) {
        eigenReadTimer.getTime();
        cout << "Loaded eigenvectors and constructed W/Q matrices (time used: "
             << eigenReadTimer.format(eigenReadTimer.getElapse()) << ")." << endl;
    }

    readCoreContextLdBlocks(dirname);
}

void Data::readCoreContextLdBlocks(const string &dirname) {
    coreContextLdBlocks.clear();
    if (precomputedMarginalizedContext) {
        coreContextLdmEigen = true;
        marginalizedContextTotalLinks = 0;
        marginalizedContextTotalFlankSnps = 0;
        for (unsigned i = 0; i < marginalizedContextLinkCounts.size(); ++i) {
            marginalizedContextTotalLinks += marginalizedContextLinkCounts[i];
            marginalizedContextTotalFlankSnps += marginalizedContextFlankCounts[i];
        }
        cout << "Using precomputed marginalized-context tile matrices from [" << dirname << "]." << endl;
        cout << "Skipped " << marginalizedContextTotalLinks
             << " in-memory sparse cross-tile LD links to "
             << marginalizedContextTotalFlankSnps << " tile-specific flank SNPs." << endl;
        return;
    }
    coreContextLdmEigen = false;
    if (!keptLdBlockInfoVec.size()) return;

    unsigned foundMeta = 0;
    uint64_t totalCoreSnps = 0;
    uint64_t totalFlankLinks = 0;
    coreContextLdBlocks.resize(keptLdBlockInfoVec.size());
    bool likelyCoreContext = false;
    for (unsigned i = 0; i < keptLdBlockInfoVec.size(); ++i) {
        if (!keptLdBlockInfoVec[i]->sourceBlockID.empty()) {
            likelyCoreContext = true;
            break;
        }
    }
    Gadget::Timer sparseSidecarTimer;
    if (likelyCoreContext)
        cout << "Loading sparse cross-tile LD files..." << endl;

    if (likelyCoreContext) {
        map<string, vector<unsigned> > blocksBySource;
        for (unsigned i = 0; i < keptLdBlockInfoVec.size(); ++i) {
            LDBlockInfo *block = keptLdBlockInfoVec[i];
            if (block->sourceBlockID.empty()) {
                throw("Error: core-context tiled ldm.info is missing SourceBlockID for block " + block->ID + ".");
            }
            blocksBySource[block->sourceBlockID].push_back(i);
        }

        for (map<string, vector<unsigned> >::iterator src = blocksBySource.begin(); src != blocksBySource.end(); ++src) {
            unordered_map<unsigned, pair<unsigned, unsigned> > coreByGlobal;
            unsigned totalCoreInSource = 0;
            vector<unsigned> sourceBlockIdx = src->second;
            std::sort(sourceBlockIdx.begin(), sourceBlockIdx.end(), [&](const unsigned a, const unsigned b) {
                LDBlockInfo *ba = keptLdBlockInfoVec[a];
                LDBlockInfo *bb = keptLdBlockInfoVec[b];
                if (ba->sourceTileID != bb->sourceTileID) return ba->sourceTileID < bb->sourceTileID;
                return ba->startSnpIdx < bb->startSnpIdx;
            });
            for (unsigned b = 0; b < src->second.size(); ++b) {
                const unsigned blockIdx = src->second[b];
                LDBlockInfo *block = keptLdBlockInfoVec[blockIdx];
                coreContextLdBlocks[blockIdx].numCoreSnps = block->numSnpInBlock;
                totalCoreInSource += block->snpInfoVec.size();
            }
            coreByGlobal.reserve(totalCoreInSource * 2 + 1);
            vector<uint32_t> sourceCurrentByLocal;
            sourceCurrentByLocal.reserve(totalCoreInSource);
            for (unsigned b = 0; b < sourceBlockIdx.size(); ++b) {
                const unsigned blockIdx = sourceBlockIdx[b];
                LDBlockInfo *block = keptLdBlockInfoVec[blockIdx];
                for (unsigned j = 0; j < block->snpInfoVec.size(); ++j) {
                    SnpInfo *snp = block->snpInfoVec[j];
                    if (!snp || snp->index < 0) {
                        sourceCurrentByLocal.push_back(SOURCE_SPARSE_LD_MISSING_INDEX);
                        continue;
                    }
                    coreByGlobal[(unsigned)snp->index] = make_pair(blockIdx, j);
                    sourceCurrentByLocal.push_back((uint32_t)snp->index);
                }
            }

            string sourceFile = existingSparseSourceLdFilename(dirname, src->first);
            if (!appendSparseSourceLdBinaryToCoreContexts(sourceFile, coreByGlobal, coreContextLdBlocks, &sourceCurrentByLocal)) {
                sourceFile = dirname + "/sourceblock" + src->first + ".sld.bin";
                if (!appendSparseSourceLdBinaryToCoreContexts(sourceFile, coreByGlobal, coreContextLdBlocks, &sourceCurrentByLocal)) {
                    throw("Error: can not open sparse cross-tile LD file for source block [" + src->first + "] in [" + dirname + "].");
                }
            }
        }

        totalCoreSnps = 0;
        totalFlankLinks = 0;
        for (unsigned i = 0; i < coreContextLdBlocks.size(); ++i) {
            totalCoreSnps += (uint64_t)coreContextLdBlocks[i].numCoreSnps;
            totalFlankLinks += (uint64_t)coreContextLdBlocks[i].ldCorr.size();
        }
        coreContextLdmEigen = true;
        sparseSidecarTimer.getTime();
        cout << "Loaded sparse cross-tile LD files (time used: "
             << sparseSidecarTimer.format(sparseSidecarTimer.getElapse()) << ")." << endl;
        cout << "Using core-context tiled eigen blocks; flank effects will be accounted for by sparse cross-tile LD offsets during SBayesR sampling." << endl;
        if (totalCoreSnps) {
            cout << "Sparse cross-tile LD files contain on average "
                 << double(totalFlankLinks) / double(totalCoreSnps)
                 << " LD link(s) per core SNP." << endl;
        }
        return;
    }

    map<string, vector<SparseSourceLdRecord> > sourceRecordCache;
    for (unsigned i = 0; i < keptLdBlockInfoVec.size(); ++i) {
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        string sourceBlockID = block->sourceBlockID;
        if (sourceBlockID.empty()) {
            string metaFile = dirname + "/block" + block->ID + ".ctx.bin";
            if (!readCoreContextTileMetaBinary(metaFile, sourceBlockID)) continue;
        }

        map<string, vector<SparseSourceLdRecord> >::iterator srcIt = sourceRecordCache.find(sourceBlockID);
        if (srcIt == sourceRecordCache.end()) {
            vector<SparseSourceLdRecord> records;
            string sourceFile = existingSparseSourceLdFilename(dirname, sourceBlockID);
            if (!readSparseSourceLdBinary(sourceFile, records)) {
                sourceFile = dirname + "/sourceblock" + sourceBlockID + ".sld.bin";
                if (!readSparseSourceLdBinary(sourceFile, records)) {
                    throw("Error: can not open sparse cross-tile LD file for source block [" + sourceBlockID + "] in [" + dirname + "] for core-context tiled eigen block " + block->ID + ".");
                }
            }
            srcIt = sourceRecordCache.insert(pair<string, vector<SparseSourceLdRecord> >(sourceBlockID, records)).first;
        }

        unordered_map<unsigned, unsigned> coreLocalByGlobal;
        coreLocalByGlobal.reserve(block->snpInfoVec.size() * 2 + 1);
        for (unsigned j = 0; j < block->snpInfoVec.size(); ++j) {
            SnpInfo *snp = block->snpInfoVec[j];
            if (!snp || snp->index < 0) continue;
            coreLocalByGlobal[(unsigned)snp->index] = j;
        }

        CoreContextLdBlock context;
        context.numCoreSnps = block->numSnpInBlock;
        const vector<SparseSourceLdRecord> &records = srcIt->second;
        for (unsigned r = 0; r < records.size(); ++r) {
            const unsigned snp1 = records[r].snp1;
            const unsigned snp2 = records[r].snp2;
            unordered_map<unsigned, unsigned>::iterator c1 = coreLocalByGlobal.find(snp1);
            unordered_map<unsigned, unsigned>::iterator c2 = coreLocalByGlobal.find(snp2);
            const bool snp1Core = c1 != coreLocalByGlobal.end();
            const bool snp2Core = c2 != coreLocalByGlobal.end();
            if (snp1Core == snp2Core) continue;
            if (snp1Core) {
                context.coreLocalIndex.push_back(c1->second);
                context.flankGlobalIndex.push_back(snp2);
            } else {
                context.coreLocalIndex.push_back(c2->second);
                context.flankGlobalIndex.push_back(snp1);
            }
            context.ldCorr.push_back(records[r].ldcor);
        }
        coreContextLdBlocks[i] = context;
        totalCoreSnps += (uint64_t)context.numCoreSnps;
        totalFlankLinks += (uint64_t)context.ldCorr.size();
        ++foundMeta;
    }

    if (foundMeta) {
        if (foundMeta != keptLdBlockInfoVec.size()) {
            throw("Error: found core-context tile metadata for only " + to_string(foundMeta) + " of "
                  + to_string(keptLdBlockInfoVec.size()) + " retained LD block(s). Please regenerate the tiled eigen folder.");
        }
        coreContextLdmEigen = true;
        if (likelyCoreContext) {
            sparseSidecarTimer.getTime();
            cout << "Loaded sparse cross-tile LD files (time used: "
                 << sparseSidecarTimer.format(sparseSidecarTimer.getElapse()) << ")." << endl;
        }
        cout << "Using core-context tiled eigen blocks; flank effects will be accounted for by sparse cross-tile LD offsets during SBayesR sampling." << endl;
        if (totalCoreSnps) {
            cout << "Sparse cross-tile LD files contain on average "
                 << double(totalFlankLinks) / double(totalCoreSnps)
                 << " LD link(s) per core SNP." << endl;
        }
        return;
    }

    unsigned foundOldRcf = 0;
    for (unsigned i = 0; i < keptLdBlockInfoVec.size(); ++i) {
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        string filename = dirname + "/block" + block->ID + ".rcf.bin";
        CoreContextLdBlock sidecar;
        if (!readCoreContextLdBlockBinary(filename, sidecar)) continue;
        if (sidecar.numCoreSnps != block->numSnpInBlock) {
            throw("Error: legacy core-context sparse LD file [" + filename + "] has " + to_string(sidecar.numCoreSnps)
                  + " core SNPs but ldm.info reports " + to_string(block->numSnpInBlock) + ".");
        }
        coreContextLdBlocks[i] = sidecar;
        ++foundOldRcf;
    }

    if (!foundOldRcf) {
        coreContextLdBlocks.clear();
        return;
    }
    if (foundOldRcf != keptLdBlockInfoVec.size()) {
        throw("Error: found legacy core-context sparse LD files for only " + to_string(foundOldRcf) + " of "
              + to_string(keptLdBlockInfoVec.size()) + " retained LD block(s). Please regenerate the tiled eigen folder.");
    }
    coreContextLdmEigen = true;
    cout << "Using legacy core-context tiled eigen blocks; flank effects will be accounted for by sparse cross-tile LD offsets during SBayesR sampling." << endl;
}

void Data::preparePrecomputedMarginalizedContext(void) {
    if (!precomputedMarginalizedContext || marginalizedContextPrepared) return;
    const unsigned nBlocks = keptLdBlockInfoVec.size();
    if (Qblocks.size() != nBlocks || wcorrBlocks.size() != nBlocks
        || marginalizedContextEigenvalues.size() != nBlocks
        || nGWASblock.size() != (int)nBlocks) {
        throw("Error: incomplete precomputed marginalized-context data.");
    }
    if (!(contextMarginalVar > 0.0f) || !isfinite(contextMarginalVar)
        || !(varPhenotypic > 0.0f) || !isfinite(varPhenotypic)) {
        throw("Error: invalid context or residual variance while preparing marginalized-context tiles.");
    }

    marginalizedContextOriginalNormWeights.resize(nBlocks);
    marginalizedContextMaxRatio = 0.0f;
    for (unsigned blk = 0; blk < nBlocks; ++blk) {
        MatrixXf &Q = Qblocks[blk];
        VectorXf &w = wcorrBlocks[blk];
        const VectorXf &d = marginalizedContextEigenvalues[blk];
        if (Q.rows() != d.size() || w.size() != d.size() || !(nGWASblock[blk] > 0.0f)) {
            throw("Error: inconsistent precomputed marginalized-context dimensions for tile "
                  + keptLdBlockInfoVec[blk]->ID + ".");
        }

        const float rho = nGWASblock[blk] * contextMarginalVar / varPhenotypic;
        VectorXf &weights = marginalizedContextOriginalNormWeights[blk];
        weights = 1.0f + rho * d.array();
        for (int row = 0; row < Q.rows(); ++row) {
            if (!(weights[row] > 0.0f) || !isfinite(weights[row]))
                throw("Error: invalid marginalized-context weight for tile "
                      + keptLdBlockInfoVec[blk]->ID + ".");
            const float shrink = 1.0f / sqrtf(weights[row]);
            Q.row(row) *= shrink;
            w[row] *= shrink;
            marginalizedContextMaxRatio = std::max(marginalizedContextMaxRatio, rho * d[row]);
        }
        cacheRealizedQdiag(keptLdBlockInfoVec[blk], Q);
    }
    marginalizedContextPrepared = true;
}

void Data::releasePseudoSummaryData(void){
    for (unsigned i = 0; i < pseudoGwasEffectTrn.size(); ++i) {
        VectorXf empty;
        pseudoGwasEffectTrn[i].swap(empty);
    }
    for (unsigned i = 0; i < pseudoGwasEffectVal.size(); ++i) {
        VectorXf empty;
        pseudoGwasEffectVal[i].swap(empty);
    }
    vector<VectorXf>().swap(pseudoGwasEffectTrn);
    vector<VectorXf>().swap(pseudoGwasEffectVal);
    VectorXf emptyNtrn;
    pseudoGwasNtrnBlock.swap(emptyNtrn);
    VectorXf emptyNval;
    pseudoGwasNValBlock.swap(emptyNval);
    VectorXf emptyBval;
    b_val.swap(emptyBval);
    VectorXf emptyBvalMask;
    b_val_mask.swap(emptyBvalMask);
}

void Data::truncateEigenMatrix(const float sumPosEigVal, const float eigenCutoff, const VectorXf &oriEigenVal, const MatrixXf &oriEigenVec, VectorXf &newEigenVal, MatrixXf &newEigenVec){
    int revIdx = (int)oriEigenVal.size();
    VectorXf cumsumNonNeg(revIdx);
    cumsumNonNeg.setZero();
    revIdx = revIdx -1;
    if(oriEigenVal(revIdx) < 0) cout << "Error, all eigenvector are negative" << endl;
    cumsumNonNeg(revIdx) = oriEigenVal(revIdx);
    revIdx = revIdx -1;
    
    while (revIdx >= 0 && oriEigenVal(revIdx) > 1e-10 ){
        cumsumNonNeg(revIdx) = oriEigenVal(revIdx) + cumsumNonNeg(revIdx + 1);
        revIdx =revIdx - 1;
        if(revIdx < 0) break;
    }
    // cout << "revIdx: " << revIdx << endl;
    // cout << "size: " << eigenVal.size()  << " eigenVal: " << eigenVal << endl;
    // cout << "cumsumNoNeg: " << cumsumNonNeg << endl;
    cumsumNonNeg = cumsumNonNeg/sumPosEigVal;  // calcualte the cumulative proportion of variance explained
    bool haveValue = false;
    // cout << "cumsumNonNeg: " << cumsumNonNeg << endl;
    // cout << "revIdx : " << revIdx << endl;
    // cout << "revIdx: "  << revIdx  << endl;
    for (revIdx = revIdx + 1; revIdx < oriEigenVal.size(); revIdx ++ ){
        // cout << "revIdx: " << revIdx << " cumsumNonNeg: " << cumsumNonNeg(revIdx) << endl;
        if(eigenCutoff >= cumsumNonNeg(revIdx) ){
            revIdx = revIdx -1;
            haveValue = true;
            break;
        }
    }
    // cout << "revIdx : " << revIdx << endl;
    if(!haveValue) revIdx = oriEigenVal.size() - 1;
    // cout << "cumsumNonNeg: " << cumsumNonNeg.size() << endl;
    // cout << "cumsumNoNeg: " << cumsumNonNeg << endl;
    // cout << "revIdx: "  << revIdx  << endl;

    newEigenVec = oriEigenVec.rightCols(oriEigenVal.size() - revIdx);
    newEigenVal = oriEigenVal.tail(oriEigenVal.size() - revIdx);
    // cout << "eigenValAdjusted size: " << eigenValAdjusted.size() << " eigenValue eventually: " << eigenValAdjusted << endl;
    //eigenvalueNum = eigenVal.size() - revIdx;
    // cout << endl;

}

void Data::readBlockLDmatrixAndMakeItSparse(const string &dirname, const unsigned block, const float chisqThreshold, const bool writeLdmTxt){
    readBlockLdmInfoFile(dirname, block);
    readBlockLdmSnpInfoFile(dirname, block);
    readBlockLdmBinaryAndMakeItSparse(dirname, block, chisqThreshold, writeLdmTxt);
}

void Data::readBlockLDmatrixAndDoEigenDecomposition(const string &dirname, const unsigned block, const float eigenCutoff, const bool writeLdmTxt){
    readBlockLdmInfoFile(dirname, block);
    readBlockLdmSnpInfoFile(dirname, block);
    readBlockLdmBinaryAndDoEigenDecomposition(dirname, block, eigenCutoff, writeLdmTxt);
}

void Data::readSparseBlockLDmatrixAndDoEigenDecomposition(const string &dirname, const unsigned block, const float eigenCutoff, const bool writeLdmTxt){
    readBlockLdmInfoFile(dirname, block);
    readBlockLdmSnpInfoFile(dirname, block);
    readSparseBlockLdmBinaryAndDoEigenDecomposition(dirname, block, eigenCutoff, writeLdmTxt);
}

void Data::readEigenMatrix(const string &dirname, const float eigenCutoff, const bool readBinary, const bool writeLdmTxt, const string &outputDir, const int quantizedBits, const bool q8Entropy){
    cout << "Reading LD matrix eigen-decomposition data..." << endl;
    //Gadget::Timer timer;
    //timer.setTime();
    
    readBlockLdmInfoFile(dirname);
    readBlockLdmSnpInfoFile(dirname);
    if (readBinary) readEigenMatrixBinaryFile(dirname, eigenCutoff, writeLdmTxt, outputDir, quantizedBits, q8Entropy);
    
    //timer.getTime();
    //cout << "Read LD data completed (time used: " << timer.format(timer.getElapse()) << ")." << endl;
}

vector<LDBlockInfo*> Data::makeKeptLDBlockInfoVec(const vector<LDBlockInfo*> &ldBlockInfoVec){
    vector<LDBlockInfo*> keptLDBlock;
    ldblockNames.clear();
    LDBlockInfo * ldblock = NULL;
    for (unsigned i=0, j=0; i< numLDBlocks; ++i) {
        ldblock = ldBlockInfoVec[i];
        if(ldblock->kept) {
            ldblock->index = j++;  // reindex inds
            keptLDBlock.push_back(ldblock);
            ldblockNames.push_back(ldblock->ID);
        }
    }
    return keptLDBlock;
}

void Data::buildMMEeigen(const string &dirname, const bool sampleOverlap, const float eigenCutoff, const bool noscale, const int quantizedBits, const bool q8Entropy, const bool makeHeldoutGwasValidation){
    includeMatchedBlocks();
    pruneLdBlocksWithNoGwasOnAllReferenceSnps(false);
    if (!numKeptLDBlocks) {
        string msg = "Error: no LD blocks remain after matching GWAS to LD. Check --gwas / --impute-summary.";
        cerr << msg << endl;
        throw std::runtime_error(msg);
    }

//    for (unsigned i=0; i<numIncdSnps; ++i) {
//        SnpInfo *snp = incdSnpInfoVec[i];
//        if (snp->gwas_b == -999) {
//            throw("Error: SNP " + snp->ID + " in the LD reference has no summary data. Run --impute-summary first.");
//        }
//    }
    unsigned nmiss = 0;
    for (unsigned i=0; i<numKeptLDBlocks; ++i) {
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        for (unsigned j=0; j<block->numSnpInBlock; ++j) {
            SnpInfo *snp = block->snpInfoVec[j];
            if (snp->gwas_b == -999) {
//                cout << "Error: SNP " + snp->ID + " in the LD reference has no summary data." << endl;
                ++nmiss;
            }
        }
    }
    if (nmiss) {
        string msg = "Error: " + to_string(nmiss) + " SNP(s) in the LD reference still have no summary data within a retained LD block (incomplete coverage). Run --impute-summary first, or check GWAS SE (finite, >0).";
        cerr << msg << endl;
        throw std::runtime_error(msg);
    }
    
    scaleGwasEffects();
    
    ldBlockInfoMap.clear();
//    snpInfoMap.clear();  // free memory from snpInfoMap which is no longer needed
    
    readEigenMatrixBinaryFileAndMakeWandQ(dirname, eigenCutoff, gwasEffectInBlock, nGWASblock, noscale, true, quantizedBits, q8Entropy);
    if (coreContextLdmEigen && quantizedBits && quantizedBits != 16) {
        throw("Error: core-context tiled eigen LD currently supports only float or q16 eigen data.");
    }
    //if (numIncdSnps!=0) constructWandQ(gwasEffectInBlock, numKeptInds);
    //if (numIncdSnps!=0) constructWandQ(eigenCutoff, noscale);

    cout << "\nData summary:" << endl;
    cout << boost::format("%40s %8s %8s\n") %"" %"mean" %"sd";
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP Phenotypic variance" %Gadget::calcMean(varySnp) %sqrt(Gadget::calcVariance(varySnp));
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP heterozygosity" %Gadget::calcMean(snp2pq) %sqrt(Gadget::calcVariance(snp2pq));
    cout << boost::format("%40s %8.0f %8.0f\n") %"GWAS SNP sample size" %Gadget::calcMean(n) %sqrt(Gadget::calcVariance(n));
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP effect (in genotype SD unit)" %Gadget::calcMean(b) %sqrt(Gadget::calcVariance(b));
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP SE" %Gadget::calcMean(se) %sqrt(Gadget::calcVariance(se));
    cout << boost::format("%40s %8.0f %8.0f\n") %"LD block size" %Gadget::calcMean(numSnpsBlock) %sqrt(Gadget::calcVariance(numSnpsBlock));
    cout << boost::format("%40s %8.0f %8.0f\n") %"LD block rank" %Gadget::calcMean(numEigenvalBlock) %sqrt(Gadget::calcVariance(numEigenvalBlock));
    cout << endl;

    if (numAnnos) setAnnoInfoVec();
    lowRankModel = true;

}

void Data::buildMMEeigenBivariate(const string &dirname, const bool sampleOverlap, const float eigenCutoff, const bool noscale){
    includeMatchedBlocks();
    pruneLdBlocksWithNoGwasOnAllReferenceSnps(true);
    if (!numKeptLDBlocks) {
        string msg = "Error: no LD blocks remain after matching GWAS to LD. Check --gwas / --impute-summary.";
        cerr << msg << endl;
        throw std::runtime_error(msg);
    }

    // Check for missing data for both traits (within retained blocks only)
    unsigned nmiss1 = 0, nmiss2 = 0;
    for (unsigned i=0; i<numKeptLDBlocks; ++i) {
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        for (unsigned j=0; j<block->numSnpInBlock; ++j) {
            SnpInfo *snp = block->snpInfoVec[j];
            if (snp->gwas_b == -999) ++nmiss1;
            if (snp->gwas_b2 == -999) ++nmiss2;
        }
    }
    if (nmiss1) {
        string msg = "Error: " + to_string(nmiss1) + " SNP(s) in the LD reference still have no summary data for trait 1 within a retained LD block (incomplete coverage). Run --impute-summary first, or check GWAS SE (finite, >0).";
        cerr << msg << endl;
        throw std::runtime_error(msg);
    }
    if (nmiss2) {
        string msg = "Error: " + to_string(nmiss2) + " SNP(s) in the LD reference still have no summary data for trait 2 within a retained LD block (incomplete coverage). Run --impute-summary first, or check GWAS SE (finite, >0).";
        cerr << msg << endl;
        throw std::runtime_error(msg);
    }
    
    scaleBivariateGwasEffects();
    
    ldBlockInfoMap.clear();
    
    // Create gwasEffectInBlock for both traits
    vector<VectorXf> gwasEffectInBlock1, gwasEffectInBlock2;
    vector<Vector2f> nGWASblockBivariate;
    
    if (numKeptLDBlocks) {
        gwasEffectInBlock1.resize(numKeptLDBlocks);
        gwasEffectInBlock2.resize(numKeptLDBlocks);
        nGWASblockBivariate.resize(numKeptLDBlocks);
        
        for (unsigned i = 0; i < numKeptLDBlocks; i++){
            LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
            gwasEffectInBlock1[i] = b(ldblock->block2GwasSnpVec);
            gwasEffectInBlock2[i] = b2(ldblock->block2GwasSnpVec);
            
            VectorXf nBlock1 = n(ldblock->block2GwasSnpVec);
            VectorXf nBlock2 = n2(ldblock->block2GwasSnpVec);
            VectorXf nBlock1Srt = nBlock1;
            VectorXf nBlock2Srt = nBlock2;
            std::sort(nBlock1Srt.data(), nBlock1Srt.data() + nBlock1Srt.size());
            std::sort(nBlock2Srt.data(), nBlock2Srt.data() + nBlock2Srt.size());
            nGWASblockBivariate[i] << nBlock1Srt[nBlock1Srt.size()/2], nBlock2Srt[nBlock2Srt.size()/2]; // median for each trait
        }
    }
    
    readEigenMatrixBinaryFileAndMakeWandQBivariate(dirname, eigenCutoff, gwasEffectInBlock1, gwasEffectInBlock2, nGWASblockBivariate, noscale);
    
    // Store nGWASblockBivariate for later use (member variable)
    this->nGWASblockBivariate = nGWASblockBivariate;
    // Also set nGWASblock to trait 1 values for compatibility
    nGWASblock.resize(nGWASblockBivariate.size());
    for (unsigned i = 0; i < nGWASblockBivariate.size(); ++i) {
        nGWASblock[i] = nGWASblockBivariate[i][0];
    }

    cout << "\nBivariate data summary:" << endl;
    cout << boost::format("%40s %8s %8s\n") %"" %"mean" %"sd";
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP Phenotypic variance (Trait 1)" %Gadget::calcMean(varySnp)  %sqrt(Gadget::calcVariance(varySnp));
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP Phenotypic variance (Trait 2)" %Gadget::calcMean(varySnp2) %sqrt(Gadget::calcVariance(varySnp2));
    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP heterozygosity"               %Gadget::calcMean(snp2pq)   %sqrt(Gadget::calcVariance(snp2pq));
    cout << boost::format("%40s %8.0f %8.0f\n") %"GWAS SNP sample size (Trait 1)"        %Gadget::calcMean(n)        %sqrt(Gadget::calcVariance(n));
    cout << boost::format("%40s %8.0f %8.0f\n") %"GWAS SNP sample size (Trait 2)"        %Gadget::calcMean(n2)       %sqrt(Gadget::calcVariance(n2));
    cout << boost::format("%40s %8.0f %8.0f\n") %"LD block size"                         %Gadget::calcMean(numSnpsBlock)    %sqrt(Gadget::calcVariance(numSnpsBlock));
    cout << boost::format("%40s %8.0f %8.0f\n") %"LD block rank"                         %Gadget::calcMean(numEigenvalBlock) %sqrt(Gadget::calcVariance(numEigenvalBlock));
    cout << endl;

    if (numAnnos) setAnnoInfoVec();
    lowRankModel = true;
}

void Data::includeMatchedBlocks(){
    // this step is to construct gwasSnp2geneVec
//    cout << "Matching blocks..." << endl;
    SnpInfo * snp;
    LDBlockInfo * ldblock;
    
    for (unsigned i=0; i<numLDBlocks; ++i){
        ldblock = ldBlockInfoVec[i];
        ldblock->block2GwasSnpVec.clear();
        if (coreContextLdmEigen) {
            for (unsigned j=0; j<ldblock->snpInfoVec.size(); ++j) {
                snp = ldblock->snpInfoVec[j];
                if (snp->included && snp->index >= 0) ldblock->block2GwasSnpVec.push_back(snp->index);
            }
            if (ldblock->block2GwasSnpVec.size() == 0) {
                cout << "Warning: Skipping LD block " << ldblock->ID
                     << ": no GWAS summary SNPs matched any SNP in this core-context block." << endl;
                ldblock->kept = false;
            } else {
                if (ldblock->block2GwasSnpVec.size() != ldblock->snpInfoVec.size()) {
                    throw("Error: core-context eigen block " + ldblock->ID + " has SNPs removed by analysis filters. V1 requires the analysis SNP set to match the generated eigen reference.");
                }
                ldblock->startSnpIdx = ldblock->block2GwasSnpVec.front();
                ldblock->endSnpIdx = ldblock->block2GwasSnpVec.back();
                ldblock->kept = true;
            }
        }
    }
    if (coreContextLdmEigen) {
        refreshKeptLdBlockMaps();
        cout << numKeptLDBlocks << " core-context tiled LD blocks are included." << endl;
        return;
    }
    for (unsigned j=0; j<numIncdSnps; ++j){
        snp = incdSnpInfoVec[j];
        ldblock = ldBlockInfoMap[snp->block];
//        cout << j << " " << snp->ID << " " << snp->block << endl;
        ldblock->block2GwasSnpVec.push_back(j);
        ldblock->snpInfoVec.push_back(snp);
    }
    for (unsigned i=0; i<numLDBlocks; ++i){
        ldblock = ldBlockInfoVec[i];
        if(ldblock->block2GwasSnpVec.size() == 0){
            cout << "Warning: Skipping LD block " << ldblock->ID
                 << ": no GWAS summary SNPs matched any SNP in this LD reference block." << endl;
            ldblock->kept = false;
        } else {
            ldblock->startSnpIdx = ldblock->block2GwasSnpVec[0];
            ldblock->endSnpIdx = ldblock->block2GwasSnpVec[ldblock->numSnpInBlock-1];
            ldblock->kept = true;
        }
    }

    refreshKeptLdBlockMaps();
    cout << numKeptLDBlocks << " LD blocks are included." << endl;
}

void Data::refreshKeptLdBlockMaps(void){
    keptLdBlockInfoVec = makeKeptLDBlockInfoVec(ldBlockInfoVec);
    numKeptLDBlocks = (unsigned) keptLdBlockInfoVec.size();

    ldblock2gwasSnpMap.clear();
    for(unsigned i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
        ldblock2gwasSnpMap.insert(pair<int, vector<int> > (i,ldblock->block2GwasSnpVec));
    }

    for (unsigned j=0; j<numIncdSnps; ++j){
        SnpInfo *snp = incdSnpInfoVec[j];
        LDBlockInfo *ldblock = ldBlockInfoMap[snp->block];
        snp->blockIdx = ldblock->index;
    }
}

unsigned Data::pruneLdBlocksWithNoGwasOnAllReferenceSnps(bool bivariate){
    unsigned pruned = 0;
    for (unsigned i=0; i<numLDBlocks; ++i){
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (!ldblock->kept) continue;
        const int m = ldblock->numSnpInBlock;
        if (m <= 0) continue;
        unsigned missing1 = 0, missing2 = 0;
        for (int j=0; j<m; ++j){
            SnpInfo *s = ldblock->snpInfoVec[j];
            if (s->gwas_b == -999) ++missing1;
            if (bivariate && s->gwas_b2 == -999) ++missing2;
        }
        const bool noTrait1 = (missing1 == (unsigned)m);
        const bool noTrait2 = bivariate && (missing2 == (unsigned)m);
        if (noTrait1 || noTrait2){
            cout << "Warning: Skipping LD block " << ldblock->ID << " (" << m << " SNPs in LD reference): ";
            if (bivariate) {
                if (noTrait1 && noTrait2)
                    cout << "all SNPs lack GWAS summary data for both traits (every SNP missing)";
                else if (noTrait1)
                    cout << "all SNPs lack GWAS summary data for trait 1 (every SNP has gwas_b == -999)";
                else
                    cout << "all SNPs lack GWAS summary data for trait 2 (every SNP has gwas_b2 == -999)";
            } else {
                cout << "all SNPs lack GWAS summary data (every SNP has gwas_b == -999)";
            }
            cout << "; excluding this block from the low-rank model." << endl;
            ldblock->kept = false;
            ++pruned;
        }
    }
    if (pruned) refreshKeptLdBlockMaps();
    return pruned;
}

//void Data::constructWandQ(const float eigenCutoff, const bool noscale){
//    VectorXf nMinusOne;
//    snp2pq.resize(numIncdSnps);
//    D.resize(numIncdSnps);
//   // ZPZdiag.resize(numIncdSnps);
//    ZPy.resize(numIncdSnps);
//    b.resize(numIncdSnps);
//    n.resize(numIncdSnps);
//    nMinusOne.resize(numIncdSnps);
//    se.resize(numIncdSnps);
//    tss.resize(numIncdSnps);
//    SnpInfo *snp;
//    for (unsigned i=0; i<numIncdSnps; ++i) {
//        snp = incdSnpInfoVec[i];
//        snp->af = snp->gwas_af;
//        snp2pq[i] = snp->twopq = 2.0f*snp->gwas_af*(1.0f-snp->gwas_af);
//        if(snp2pq[i]==0) cout << "Error: SNP " << snp->ID << " af " << snp->af << " has 2pq = 0." << endl;
//        D[i] = snp2pq[i]*snp->gwas_n;
//        b[i] = snp->gwas_b * sqrt(snp2pq[i]); // scale the marginal effect so that it's in per genotype SD unit
//        n[i] = snp->gwas_n;
//        nMinusOne[i] = snp->gwas_n - 1;
//        se[i]= snp->gwas_se * sqrt(snp2pq[i]);
//        tss[i] = D[i]*(n[i]*se[i]*se[i] + b[i]*b[i]);
//        ZPy[i] = n[i]*b[i];
//        //  D[i] = 1.0/(se[i]*se[i]+b[i]*b[i]/snp->gwas_n);  // NEW!
//        //  snp2pq[i] = snp->twopq = D[i]/snp->gwas_n;       // NEW!
//    }
//    cout << endl;
//    LDBlockInfo * ldblock;
//    wcorrBlocks.resize(numKeptLDBlocks);
//    // Qblocks.resize(numKeptLDBlocks);
//    numSnpsBlock.resize(numKeptLDBlocks);
//    numEigenvalBlock.resize(numKeptLDBlocks);
//    Qblocks.clear();
//    VectorXf sqrtLambda;
//    // save gwas marginal effect into block
//    gwasEffectInBlock.resize(numKeptLDBlocks);
//    for (unsigned i = 0; i < numKeptLDBlocks; i++){
//        ldblock = keptLdBlockInfoVec[i];
//        gwasEffectInBlock[i] = b(ldblock->block2GwasSnpVec);
//        // calculate wbcorr and Qblocks
//        sqrtLambda = eigenValLdBlock[i].array().sqrt();
//        //cout << eigenVecLdBlock[i].transpose().rows() << " " << eigenVecLdBlock[i].transpose().cols() << " " << gwasEffectInBlock[i].size() << endl;
//        wcorrBlocks[i] = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (eigenVecLdBlock[i].transpose() * gwasEffectInBlock[i] );
//        // cout << "eigenVecLdBlock[i]: " << eigenVecLdBlock[i] << endl;
//        // cout << "wcorrBlocks[i]: " << wcorrBlocks[i] << endl;
//        // cout << "sqrtLambda: " << sqrtLambda << endl;
//        // cout << gwasEffectInBlock[i] << endl;
//        MatrixXf tmpQblocks = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
//        MatrixDat matrixDat = MatrixDat(ldblock->snpNameVec,tmpQblocks );
//        // cout << "Qblock: " << endl;
//        // cout << matrixDat.values << endl;
//        Qblocks.push_back(matrixDat);
//        numSnpsBlock[i] = Qblocks[i].ncol;
//        numEigenvalBlock[i] = Qblocks[i].nrow;
//    }
//    
//    //b.array() -= b.mean();  // DO NOT CENTER b
//    // estimate phenotypic variance based on the input allele frequencies in GWAS
//
//    //  Vp_buf = h_buf * N_buf * se_buf * se_buf + h_buf * b_buf * b_buf * N_buf / (N_buf - 1.0);
//    //VectorXf ypySrt = D.array()*n.array()*se.array().square() + D.array() * b.array().square() * n.array() / nMinusOne.array();
//    VectorXf ypySrt = D.array()*(n.array()*se.array().square()+b.array().square());
//    VectorXf varpSrt = ypySrt.array()/n.array();
//    std::sort(ypySrt.data(), ypySrt.data() + ypySrt.size());
//    std::sort(varpSrt.data(), varpSrt.data() + varpSrt.size());
//    ypy = ypySrt[ypySrt.size()/2];  // median
//    varPhenotypic = varpSrt[varpSrt.size()/2];
//    //cout << "varPhenotypic: " << varPhenotypic << endl;
//    VectorXf nSrt = n;
//    std::sort(nSrt.data(), nSrt.data() + nSrt.size());
//    numKeptInds = nSrt[nSrt.size()/2]; // median
//    
//    nGWASblock.resize(numKeptLDBlocks);
//    for (unsigned i=0; i<numKeptLDBlocks; ++i) {
//        nGWASblock[i] = numKeptInds;
//    }
//
//    for (unsigned i=0; i<numIncdSnps; ++i) {
//        snp = incdSnpInfoVec[i];
//        D[i] = varPhenotypic/(se[i]*se[i]+b[i]*b[i]/snp->gwas_n);  // NEW!
//        snp2pq[i] = snp->twopq = D[i]/snp->gwas_n;       // NEW!
//        tss[i] = D[i]*(n[i]*se[i]*se[i] + b[i]*b[i]);
//        // Need to adjust R and C models X'X matrix depending scale of genotypes or not
//        if (noscale == true) {
//            D[i] = snp2pq[i]*snp->gwas_n;
//        } else {
//            D[i] = snp->gwas_n;
//        }
//    }
//    // cout << endl << "snp2pq: " << endl << snp2pq << endl;
//    //ypy = numKeptInds;
//    // NEW END
//    // data summary
//    cout << "\nData summary:" << endl;
//    cout << boost::format("%40s %8s %8s\n") %"" %"mean" %"sd";
//    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP Phenotypic variance" %Gadget::calcMean(varpSrt) %sqrt(Gadget::calcVariance(varpSrt));
//    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP heterozygosity" %Gadget::calcMean(snp2pq) %sqrt(Gadget::calcVariance(snp2pq));
//    cout << boost::format("%40s %8.0f %8.0f\n") %"GWAS SNP sample size" %Gadget::calcMean(n) %sqrt(Gadget::calcVariance(n));
//    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP effect (in genotype SD unit)" %Gadget::calcMean(b) %sqrt(Gadget::calcVariance(b));
//    cout << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP SE" %Gadget::calcMean(se) %sqrt(Gadget::calcVariance(se));
//    cout << boost::format("%40s %8.3f %8.3f\n") %"LD block size" %Gadget::calcMean(numSnpsBlock) %sqrt(Gadget::calcVariance(numSnpsBlock));
//    cout << boost::format("%40s %8.3f %8.3f\n") %"LD block rank" %Gadget::calcMean(numEigenvalBlock) %sqrt(Gadget::calcVariance(numEigenvalBlock));
//}


void Data::mergeLdmInfo(const string &outLDmatType, const string &dirname, const bool print, const bool mergePairwiseLD, const bool cleanupBlockMergeFiles) {
    if (outLDmatType != "block") throw("Error: --merge-ldm-info only works for block LD matrices at the moment!");
    
    string dir_path = dirname; // Replace with your folder path
    string search_str = "snp.info";
    string outSnpInfoFile = dirname + "/snp.info";
    string outldmInfoFile = dirname + "/ldm.info";
    string outpwldFile    = dirname + "/rsq0.5.pwld";
    string outpwldGzFile  = outpwldFile + ".gz";
    string mergeMarkerFile = dirname + "/ldm.merge.done";
    DIR* dirp = opendir(dir_path.c_str());
    
    if (dirp == NULL) {
        throw("Error opening directory [" + dirname + "]");
    }
    
    // find out all ldm in the folder
    vector<string> file_list;
    
    dirent* dp;
    while ((dp = readdir(dirp)) != NULL) {
        string file_name = dp->d_name;
        if (file_name.find(search_str) != string::npos) {
            file_list.push_back(file_name);
        }
    }
    
    closedir(dirp);
    
    search_str = "block";

    vector<string> filtered_list;
    for (vector<string>::iterator it = file_list.begin(); it != file_list.end(); ++it) {
        if (it->find("block") != string::npos) {
            filtered_list.push_back(*it);
        }
    }
    file_list.swap(filtered_list); // replace original with filtered

    vector<string> activeLocks;
    if (findBlockOutputLocksLocal(dirname, activeLocks)) {
        string msg = "block-specific output lock(s) still exist in [" + dirname + "]; block jobs may still be running";
        for (unsigned i = 0; i < activeLocks.size() && i < 20; ++i)
            msg += " " + activeLocks[i];
        if (activeLocks.size() > 20) msg += " ...";
        msg += ". Please wait for all jobs to finish before merging.";
        if (print) throw("Error: " + msg);
        else {
            cout << "Warning: " << msg << endl;
            return;
        }
    }

    if (mergePairwiseLD) installRecordedPairwiseLdSourceLocal(dirname);

    if (file_list.size() == 0) {
        if (combinedMergeOutputsReadyLocal(dirname, mergePairwiseLD)) {
            cout << "Using existing merged LD metadata in [" << dirname << "]." << endl;
            return;
        }
        string msg = "there are no block-specific SNP info files (block*.snp.info) in folder [" + dirname + "] to merge.";
        if (print) throw ("Error: " + msg);
        else return;
    }
    if (file_list.size() == 1) {
        if (print) throw ("Error: there is only one block-specific info file in folder [" + dirname + "], so no need to merge.");
    }

    
    vector<string> blockFileIds;
    set<unsigned> numericBlockIdxSet;

    for (vector<string>::iterator it = file_list.begin(); it != file_list.end(); ++it) {
        size_t block_pos = it->find(search_str);
        const string suffix = ".snp.info";
        size_t suffix_pos = it->rfind(suffix);
        if (block_pos != string::npos && suffix_pos != string::npos && suffix_pos > block_pos + search_str.size()) {
            string block_id_str = it->substr(block_pos + search_str.size(), suffix_pos - block_pos - search_str.size());
            blockFileIds.push_back(block_id_str);
            bool numeric = !block_id_str.empty();
            for (unsigned c = 0; c < block_id_str.size(); ++c) {
                if (!isdigit(block_id_str[c])) {
                    numeric = false;
                    break;
                }
            }
            if (numeric) numericBlockIdxSet.insert((unsigned)atoi(block_id_str.c_str()));
        }
    }
    std::sort(blockFileIds.begin(), blockFileIds.end(), [](const string &a, const string &b) {
        char *enda = NULL;
        char *endb = NULL;
        long na = strtol(a.c_str(), &enda, 10);
        long nb = strtol(b.c_str(), &endb, 10);
        bool anum = enda && *enda == '\0';
        bool bnum = endb && *endb == '\0';
        if (anum && bnum) return na < nb;
        if (anum != bnum) return anum;
        return a < b;
    });

    bool mergeSourceTile = false;
    {
        string firstLdmInfoFile = dirname + "/block" + blockFileIds.front() + ".ldm.info";
        ifstream probe(firstLdmInfoFile.c_str());
        string headerProbe;
        getline(probe, headerProbe);
        mergeSourceTile = headerProbe.find("SourceBlockID") != string::npos && headerProbe.find("TileID") != string::npos;
        if (headerProbe.find("TilePhase") != string::npos || headerProbe.find("TileGroupSize") != string::npos
            || headerProbe.find("CoreStartSnpIdx") != string::npos || headerProbe.find("CoreEndSnpIdx") != string::npos) {
            throw("Error: overlapping and phased tiled eigen LD formats are no longer supported. Regenerate the LD reference as ordinary or core-context eigen data.");
        }
    }

    unsigned nldm = blockFileIds.size();
    if (!nldm) {
        string msg = "no valid block IDs were found in block-specific SNP info files in folder [" + dirname + "].";
        if (print) throw("Error: " + msg);
        else {
            cout << "Warning: " << msg << endl;
            return;
        }
    }

    uint64_t blockInfoSignature = 0;
    bool haveBlockInfoSignature = blockInfoSignatureLocal(dirname, blockFileIds, blockInfoSignature);
    if (haveBlockInfoSignature && combinedMergeOutputsReadyLocal(dirname, mergePairwiseLD)) {
        unsigned markerBlockCount = 0;
        uint64_t markerSignature = 0;
        bool markerPairwiseReady = false;
        if (readMergeMarkerLocal(mergeMarkerFile, markerBlockCount, markerSignature, markerPairwiseReady)
            && markerBlockCount == nldm
            && markerSignature == blockInfoSignature
            && (!mergePairwiseLD || markerPairwiseReady)) {
            cout << "Using cached merged LD metadata in [" << dirname << "]." << endl;
            return;
        }
    }

    if (!mergeSourceTile && numericBlockIdxSet.size() == blockFileIds.size()) {
        unsigned firstBlock = *numericBlockIdxSet.begin();
        unsigned lastBlock = *numericBlockIdxSet.rbegin();
        vector<unsigned> missingBlocks;
        for (unsigned blockID = firstBlock; blockID <= lastBlock; ++blockID) {
            if (!numericBlockIdxSet.count(blockID)) missingBlocks.push_back(blockID);
        }
        if (!missingBlocks.empty()) {
            string msg = "block-specific info files are incomplete in [" + dirname + "]; missing block ID(s)";
            for (unsigned i = 0; i < missingBlocks.size() && i < 20; ++i)
                msg += " " + to_string(missingBlocks[i]);
            if (missingBlocks.size() > 20)
                msg += " ...";
            msg += ". Please wait until all per-block jobs finish before merging.";
            if (print) throw("Error: " + msg);
            else {
                cout << "Warning: " << msg << endl;
                return;
            }
        }
    }

    if (nldm > 100)
        cout << "Checking completeness of " << nldm << " block-specific info files..." << endl;
    unsigned checkIdx = 0;
    for (vector<string>::iterator chk = blockFileIds.begin(); chk != blockFileIds.end(); ++chk, ++checkIdx) {
        string snpInfoFile = dirname + "/block" + *chk + ".snp.info";
        string ldmInfoFile = dirname + "/block" + *chk + ".ldm.info";
        ifstream snpIn(snpInfoFile.c_str());
        ifstream ldmIn(ldmInfoFile.c_str());
        if (!snpIn || !ldmIn) {
            string msg = "missing paired block-specific info file(s) for block " + *chk
                + " in [" + dirname + "]. Expected both block*.snp.info and block*.ldm.info.";
            if (print) throw("Error: " + msg);
            else {
                cout << "Warning: " << msg << endl;
                return;
            }
        }

        string line;
        unsigned snpLines = 0;
        getline(snpIn, line);
        while (getline(snpIn, line)) {
            if (!line.empty()) ++snpLines;
        }
        if (!snpLines) {
            string msg = "block" + *chk + ".snp.info has no SNP records; block job may still be writing.";
            if (print) throw("Error: " + msg);
            else {
                cout << "Warning: " << msg << endl;
                return;
            }
        }

        unsigned ldmLines = 0;
        int snpNum = -1;
        unsigned snpNumSum = 0;
        string header;
        getline(ldmIn, header);
        bool sourceTileHeader = header.find("SourceBlockID") != string::npos && header.find("TileID") != string::npos;
        while (getline(ldmIn, line)) {
            if (line.empty()) continue;
            ++ldmLines;
            std::istringstream iss(line);
            string id, startSnpID, endSnpID;
            int chr, blockStart, blockEnd;
            if (!(iss >> id >> chr >> blockStart >> startSnpID >> blockEnd >> endSnpID >> snpNum)) {
                string msg = "cannot parse ldm info for block " + *chk + " in [" + ldmInfoFile + "].";
                if (print) throw("Error: " + msg);
                else {
                    cout << "Warning: " << msg << endl;
                    return;
                }
            }
            if (sourceTileHeader) {
                string eigenFile = dirname + "/block" + id + ".eigen.bin";
                string eigenQ16File = dirname + "/block" + id + ".eigen.q16.bin";
                string eigenBundleFile = dirname + "/block" + *chk + ".tiles.eigen.bin";
                string eigenQ16BundleFile = dirname + "/block" + *chk + ".tiles.eigen.q16.bin";
                if (!fileExistsAndNonEmptyLocal(eigenFile) && !fileExistsAndNonEmptyLocal(eigenQ16File)
                    && !fileExistsAndNonEmptyLocal(eigenBundleFile) && !fileExistsAndNonEmptyLocal(eigenQ16BundleFile)) {
                    string msg = "incomplete tiled eigen output for source block " + *chk
                        + ": ldm.info lists tile [" + id + "] but neither per-tile eigen files nor source-block eigen bundle exists. "
                        + "Block job may still be running or was interrupted.";
                    if (print) throw("Error: " + msg);
                    else {
                        cout << "Warning: " << msg << endl;
                        return;
                    }
                }
            }
            if (snpNum > 0) snpNumSum += (unsigned)snpNum;
        }
        if (sourceTileHeader) {
            string sparseFile = sparseSourceLdFilename(dirname, *chk, false);
            string sparseQ16File = sparseSourceLdFilename(dirname, *chk, true);
            if (!fileExistsAndNonEmptyLocal(sparseFile) && !fileExistsAndNonEmptyLocal(sparseQ16File)) {
                string msg = "incomplete tiled eigen output for source block " + *chk
                    + ": sparse cross-tile LD file is missing or empty.";
                if (print) throw("Error: " + msg);
                else {
                    cout << "Warning: " << msg << endl;
                    return;
                }
            }
        }
        bool inconsistent = false;
        if (sourceTileHeader) inconsistent = ldmLines == 0 || snpNumSum != snpLines;
        else inconsistent = ldmLines != 1 || snpNum < 0 || (unsigned)snpNum != snpLines;
        if (inconsistent) {
            string msg = "inconsistent block-specific info for block " + *chk
                + ": snp.info has " + to_string(snpLines)
                + " SNP records but ldm.info reports " + to_string(sourceTileHeader ? (long long)snpNumSum : (long long)snpNum)
                + ". Block job may still be writing.";
            if (print) throw("Error: " + msg);
            else {
                cout << "Warning: " << msg << endl;
                return;
            }
        }
    }
    
    cout << "Merging " << nldm << " LD matrices info files" << endl;
    
    vector<string>::iterator it = blockFileIds.begin();

    {
        ifstream existingLdm(outldmInfoFile.c_str());
        if (existingLdm) {
            string line;
            unsigned existingBlockCount = 0;
            set<string> existingBlockIdxSet;
            set<string> existingSourceIdxSet;
            getline(existingLdm, line);
            while (getline(existingLdm, line)) {
                if (line.empty()) continue;
                ++existingBlockCount;
                std::istringstream iss(line);
                string id;
                if (iss >> id) {
                    existingBlockIdxSet.insert(id);
                    string sourceID = id;
                    size_t tilePos = sourceID.find(".tile");
                    if (tilePos != string::npos) sourceID = sourceID.substr(0, tilePos);
                    existingSourceIdxSet.insert(sourceID);
                }
            }
            bool newCoversExisting = true;
            set<string> visibleIDs(blockFileIds.begin(), blockFileIds.end());
            for (set<string>::iterator oldIt = existingBlockIdxSet.begin(); oldIt != existingBlockIdxSet.end(); ++oldIt) {
                string sourceID = *oldIt;
                size_t tilePos = sourceID.find(".tile");
                if (tilePos != string::npos) sourceID = sourceID.substr(0, tilePos);
                if (!visibleIDs.count(sourceID)) {
                    newCoversExisting = false;
                    break;
                }
            }
            unsigned comparableExistingCount = mergeSourceTile ? (unsigned)existingSourceIdxSet.size() : existingBlockCount;
            if (comparableExistingCount > nldm || !newCoversExisting) {
                string msg = "existing combined ldm.info has " + to_string(existingBlockCount)
                    + " blocks, but the currently visible block-specific info files cover " + to_string(nldm)
                    + " block(s) and do not contain all existing merged block IDs; keeping the existing combined info to avoid dropping blocks.";
                if (print) throw("Error: " + msg);
                else {
                    cout << "Warning: " << msg << endl;
                    return;
                }
            }
        }
    }

    // If a combined rsq0.5.pwld already exists (e.g. copied from the source reference),
    // merge only snp.info / ldm.info from per-block files and keep the existing pwld source
    // (do not require block*.rsq0.5.pwld).
    string existingCombinedPwldFile = mergePairwiseLD
        ? findPairwiseLdWithRecordsForBaseLocal(outpwldFile) : "";
    const bool useExistingCombinedPwld = !existingCombinedPwldFile.empty();
    if (mergePairwiseLD && useExistingCombinedPwld) {
        cout << "Found existing combined pairwise LD file [" + existingCombinedPwldFile + "]; merging snp/ldm info only and keeping this pwld file." << endl;
    }
    const bool writeCombinedPwld = mergePairwiseLD && !useExistingCombinedPwld;
    if (!mergePairwiseLD) {
        cout << "Skipping combined pairwise LD merge for this analysis." << endl;
    } else if (writeCombinedPwld && nldm > 100) {
        cout << "Building combined pairwise LD file; this can take a while for large references." << endl;
    }

    ofstream out1(outSnpInfoFile.c_str());
    out1 << "Chrom ID Index GenPos PhysPos A1 A2 A1Freq N Block\n";
    out1 << fixed << setprecision(6);
    
    ofstream out2(outldmInfoFile.c_str());
    if (mergeSourceTile) {
        out2 << "Block Chrom StartSnpIdx StartSnpID EndSnpIdx EndSnpID NumSnps SourceBlockID TileID\n";
    } else {
        out2 << "Block Chrom StartSnpIdx StartSnpID EndSnpIdx EndSnpID NumSnps\n";
    }

    gzFile out3gz = NULL;
    if (writeCombinedPwld) {
        out3gz = gzopen(outpwldGzFile.c_str(), "wb");
        if (!out3gz) throw("Error: cannot open file [" + outpwldGzFile + "] to write.");
        gzprintf(out3gz, "%12s %12s %12s\n", "SNP1", "SNP2", "LDcorrelation");
        remove(outpwldFile.c_str());
        remove((outpwldFile + ".bin").c_str());
    }

    unsigned snpIdx = 0;
    unsigned ldmIdx = 0;
    
    for (unsigned i=0; i<nldm; ++i) {
        
        string snpInfoFile = dirname + "/block" + *it + ".snp.info";
        string ldmInfoFile = dirname + "/block" + *it + ".ldm.info";
        string pwldFile = dirname + "/block" + *it + ".rsq0.5.pwld";

        // read snp info file
        ifstream in1(snpInfoFile.c_str());
        if (!in1) throw ("Error: can not open the file [" + snpInfoFile + "] to read.");
        //cout << "Reading SNP info from file [" + snpInfoFile + "]." << endl;
        string header;
        string id, allele1, allele2;
        int chr, physPos,ld_n;
        float genPos;
        float allele1Freq;
        int idx = 0;
        int index;
        string blockID;
        unordered_map<string, int> snpID2index;
        getline(in1, header);
        while (in1 >> chr >> id >> index >> genPos >> physPos >> allele1 >> allele2 >> allele1Freq >> ld_n >> blockID) {
            unsigned outIndex = snpIdx;
            out1 << chr << ' ' << id << ' ' << outIndex << ' ' << genPos << ' ' << physPos << ' '
                 << allele1 << ' ' << allele2 << ' ' << allele1Freq << ' ' << ld_n << ' ' << blockID << '\n';
            snpID2index[id] = outIndex;
            ++snpIdx;
        }
        in1.close();
        
        
        // read ldm info file
        ifstream in2(ldmInfoFile.c_str());
        if (!in2) throw ("Error: can not open the file [" + ldmInfoFile + "] to read.");
        //cout << "Reading LDM info from file [" + ldmInfoFile + "]." << endl;
        
        int blockStart, blockEnd, snpNum;
        string startSnpID, endSnpID;
        getline(in2, header);
        string ldmLine;
        while (getline(in2, ldmLine)) {
            if (ldmLine.empty()) continue;
            std::istringstream iss(ldmLine);
            if (!(iss >> id >> chr >> blockStart >> startSnpID >> blockEnd >> endSnpID >> snpNum)) continue;
            if (mergeSourceTile) {
                string sourceBlockID;
                int tileID;
                if (!(iss >> sourceBlockID >> tileID)) {
                    throw("Error: inconsistent core-context tiled ldm.info line in [" + ldmInfoFile + "].");
                }
                out2 << id << ' ' << chr << ' ' << snpID2index[startSnpID] << ' ' << startSnpID << ' '
                     << snpID2index[endSnpID] << ' ' << endSnpID << ' ' << snpNum << ' '
                     << sourceBlockID << ' ' << tileID << '\n';
            } else {
                out2 << id << ' ' << chr << ' ' << snpID2index[startSnpID] << ' ' << startSnpID << ' '
                     << snpID2index[endSnpID] << ' ' << endSnpID << ' ' << snpNum << '\n';
            }
            
            ++ldmIdx;
        }
        in2.close();
        
        
        // read per-block pwld only when building merged rsq0.5.pwld (not when keeping existing combined file)
        if (writeCombinedPwld) {
            string pwldGzFile = pwldFile + ".gz";
            string pwldBinFile = pwldFile + ".bin";
            if (pairwiseLdHasRecordsLocal(pwldGzFile)) {
                gzFile in3gz = gzopen(pwldGzFile.c_str(), "rb");
                if (!in3gz) throw("Error: can not open the file [" + pwldGzFile + "] to read.");
                string line;
                readGzLineLocal(in3gz, line);
                while (readGzLineLocal(in3gz, line)) {
                    if (line.empty()) continue;
                    gzputs(out3gz, line.c_str());
                    gzputc(out3gz, '\n');
                }
                gzclose(in3gz);
            } else if (pairwiseLdHasRecordsLocal(pwldBinFile)) {
                FILE *in3bin = fopen(pwldBinFile.c_str(), "rb");
                if (!in3bin) throw("Error: can not open the file [" + pwldBinFile + "] to read.");
                uint64_t fileNumSnps = 0, fileNumRecords = 0;
                vector<string> fileSnpIDs;
                readPwldBinHeaderLocal(in3bin, pwldBinFile, fileNumSnps, fileNumRecords, fileSnpIDs);
                for (uint64_t r = 0; r < fileNumRecords; ++r) {
                    PwldBinRecord inrec;
                    if (fread(&inrec, sizeof(PwldBinRecord), 1, in3bin) != 1) {
                        fclose(in3bin);
                        throw("Error: incomplete binary pairwise LD record in [" + pwldBinFile + "].");
                    }
                    if (inrec.snp1 >= fileSnpIDs.size() || inrec.snp2 >= fileSnpIDs.size()) {
                        fclose(in3bin);
                        throw("Error: SNP index out of range in binary pairwise LD file [" + pwldBinFile + "].");
                    }
                    const string &snp1ID = fileSnpIDs[inrec.snp1];
                    const string &snp2ID = fileSnpIDs[inrec.snp2];
                    unordered_map<string, int>::iterator idx1 = snpID2index.find(snp1ID);
                    unordered_map<string, int>::iterator idx2 = snpID2index.find(snp2ID);
                    if (idx1 == snpID2index.end() || idx2 == snpID2index.end()) continue;
                    gzprintf(out3gz, "%12s %12s %12.6f\n", snp1ID.c_str(), snp2ID.c_str(), inrec.ldcor);
                }
                fclose(in3bin);
            } else {
                ifstream in3(pwldFile.c_str());
                if (!in3) {
                    cout << "Warning: can not open the file [" + pwldFile + "] or [" + pwldBinFile + "] to read; skipping pairwise LD for this block." << endl;
                } else {
                    getline(in3, header);
                    string line;
                    while (getline(in3, line)) {
                        if (line.empty()) continue;
                        gzputs(out3gz, line.c_str());
                        gzputc(out3gz, '\n');
                    }
                    in3.close();
                }
            }
        }

        ++it;
    }
    
    out1.close();
    out2.close();
    if (writeCombinedPwld) {
        gzclose(out3gz);
    }

    cout << "Written " << snpIdx << " SNPs info into file [" + outSnpInfoFile + "]." << endl;
    cout << "Written " << ldmIdx << " LDMs info into file [" + outldmInfoFile + "]." << endl;
    if (!mergePairwiseLD) {
        cout << "Pairwise LD merge skipped." << endl;
    } else if (useExistingCombinedPwld) {
        cout << "Left pairwise LD file unchanged: [" + existingCombinedPwldFile + "]." << endl;
    } else {
        cout << "Written pairwise LD correlations (rsq>0.5) into file [" + outpwldGzFile + "]." << endl;
    }
    const bool pairwiseLDReady = mergePairwiseLD
        && !findPairwiseLdWithRecordsForBaseLocal(outpwldFile).empty();
    if (mergePairwiseLD && !pairwiseLDReady) {
        cout << "Warning: the merged pairwise LD output contains no LD records. "
             << "Keeping block-specific merge inputs so pairwise LD can be repaired or merged again." << endl;
    }
    if (haveBlockInfoSignature)
        writeMergeMarkerLocal(mergeMarkerFile, nldm, blockInfoSignature, pairwiseLDReady);
    if (cleanupBlockMergeFiles && pairwiseLDReady) {
        unsigned removed = removeBlockSpecificMergeInputsLocal(dirname);
        if (removed)
            cout << "Removed " << removed << " block-specific info/pairwise LD file(s) after successful merge." << endl;
    } else if (!cleanupBlockMergeFiles && pairwiseLDReady) {
        cout << "Kept block-specific info/pairwise LD files as requested." << endl;
    }
}

void Data::mergeBlockGwasSummary(const string &gwasSummaryFile, const string &title) {
    string dir_path = "."; // Replace with your folder path
    string search_str1 = gwasSummaryFile + ".block";
    string search_str2 = ".ma";
    DIR* dirp = opendir(dir_path.c_str());
    
    if (dirp == NULL) {
        throw("Error opening directory [" + dir_path + "]");
    }
    
    // find out all .ma files in the folder
    vector<string> file_list;
    
    dirent* dp;
    while ((dp = readdir(dirp)) != NULL) {
        string file_name = dp->d_name;
        if (file_name.find(search_str1) != string::npos && file_name.find(search_str2) != string::npos) {
            file_list.push_back(file_name);
        }
    }
    
    closedir(dirp);
    
    map<unsigned, string> blockIdxMap;
    blockIdxMap.clear();
    
    for (vector<string>::iterator it = file_list.begin(); it != file_list.end(); ++it) {
        size_t block_pos = it->find(search_str1);
        if (block_pos == string::npos) continue;
        const size_t number_pos = block_pos + search_str1.size();
        const size_t dot_pos = it->find('.', number_pos);
        if (dot_pos == string::npos || dot_pos == number_pos) continue;
        const string block_num_str = it->substr(number_pos, dot_pos - number_pos);
        bool numeric = true;
        for (unsigned c = 0; c < block_num_str.size(); ++c) {
            if (!isdigit(block_num_str[c])) {
                numeric = false;
                break;
            }
        }
        if (!numeric) continue;
        unsigned block_num = (unsigned)strtoul(block_num_str.c_str(), NULL, 10);
        if (blockIdxMap.count(block_num))
            throw("Error: multiple block summary files found for block " + block_num_str + ".");
        blockIdxMap[block_num] = *it;
    }
    
    if (blockIdxMap.size() == 0) {
        throw ("Error: there is no info file to merge in folder [" + dir_path + "].");
    }

    vector<unsigned> missingBlocks;
    uint64_t missingBlockCount = 0;
    unsigned previousBlock = blockIdxMap.begin()->first;
    map<unsigned, string>::const_iterator blockIt = blockIdxMap.begin();
    ++blockIt;
    for (; blockIt != blockIdxMap.end(); ++blockIt) {
        const uint64_t gap = (uint64_t)blockIt->first - (uint64_t)previousBlock - 1u;
        missingBlockCount += gap;
        for (uint64_t offset = 1; offset <= gap && missingBlocks.size() < 20; ++offset)
            missingBlocks.push_back((unsigned)((uint64_t)previousBlock + offset));
        previousBlock = blockIt->first;
    }
    if (missingBlockCount) {
        string msg = "block summary files are incomplete; missing block ID(s)";
        for (unsigned i = 0; i < missingBlocks.size() && i < 20; ++i)
            msg += " " + to_string(missingBlocks[i]);
        if (missingBlockCount > missingBlocks.size()) msg += " ...";
        msg += ". A zero-match block should have a header-only completed-block file; otherwise rerun or wait for the missing block job.";
        throw("Error: " + msg);
    }
    
    unsigned nBlk = blockIdxMap.size();
    
    cout << "Merging GWAS summary statistics files across " + to_string(nBlk) + " blocks..." << endl;
    
    map<unsigned, string>::iterator it = blockIdxMap.begin();
    
    string outMaFile = dir_path + "/" + title + ".ma";
    
    ofstream out(outMaFile.c_str());
    out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n") % "SNP" % "A1" % "A2" % "freq" % "b" % "se" % "p" % "N";
    
    unsigned snpIdx = 0;
    
    for (unsigned i=0; i<nBlk; ++i) {
        
        string mafile = it->second;
        
        // read snp info file
        ifstream in(mafile.c_str());
        if (!in) throw ("Error: can not open the file [" + mafile + "] to read.");
        cout << "Reading summary statistics from file [" + mafile + "]." << endl;

        string header;
        getline(in, header);

        string id, allele1, allele2, freq, b, se, pval, n;
        while (in >> id >> allele1 >> allele2 >> freq >> b >> se >> pval >> n) {
            out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n")
            % id
            % allele1
            % allele2
            % freq
            % b
            % se
            % pval
            % n;
            ++snpIdx;
        }
        in.close();
        
        ++it;
    }
    
    out.close();
    
    cout << "Written " << snpIdx << " SNPs info into file [" + outMaFile + "]." << endl;
}

void Data::constructPseudoSummaryData(){
    cout << "Constructing pseudo summary statistics for training and validation data sets, with 90% sample size for training and 10% for validation." << endl;
    
    pseudoGwasEffectTrn.resize(numKeptLDBlocks);
    pseudoGwasEffectVal.resize(numKeptLDBlocks);
    
    VectorXf n_trn(numKeptLDBlocks);
    VectorXf n_val(numKeptLDBlocks);
    pseudoGwasNtrnBlock.resize(numKeptLDBlocks);
    b_val.resize(numIncdSnps);
    VectorXf emptyMask;
    b_val_mask.swap(emptyMask);

    for (unsigned i=0; i<numKeptLDBlocks; ++i) {
        LDBlockInfo* block = keptLdBlockInfoVec[i];
        
        n_trn[i] = 0.9*nGWASblock[i];
        n_val[i] = nGWASblock[i] - n_trn[i];
        pseudoGwasNtrnBlock[i] = n_trn[i];

        long size = eigenValLdBlock[i].size();
        VectorXf rnd(size);
        for (unsigned j=0; j<size; ++j) {
            rnd[j] = Stat::snorm();
        }

        VectorXf noise;
        if (i < Qblocks.size() && Qblocks[i].size()) {
            noise = Qblocks[i].transpose() * rnd;
        } else if (i < quantizedEigenBlocks.size() && quantizedEigenBlocks[i].active()) {
            MatrixXf Qtmp;
            quantizedEigenBlocks[i].materializeQ(Qtmp);
            noise = Qtmp.transpose() * rnd;
        } else if (i < eigenVecLdBlock.size() && eigenVecLdBlock[i].size()) {
            noise = eigenVecLdBlock[i] * (eigenValLdBlock[i].array().sqrt().matrix().asDiagonal() * rnd);
        } else {
            throw("Error: cannot construct pseudo summary data because neither Q nor eigenvectors are available for LD block " + block->ID + ".");
        }

        pseudoGwasEffectTrn[i] = gwasEffectInBlock[i] + sqrt(1.0/n_trn[i] - 1.0/nGWASblock[i]) * noise;

        pseudoGwasEffectVal[i] = nGWASblock[i]/n_val[i] * gwasEffectInBlock[i] - n_trn[i]/n_val[i] * pseudoGwasEffectTrn[i];
        b_val.segment(block->startSnpIdx, block->numSnpInBlock) = pseudoGwasEffectVal[i];
    }
}

void Data::readEigenMatrixBinaryFileAndMakeWandQBivariate(const string &dirname, const float eigenCutoff, const vector<VectorXf> &GWASeffects1, const vector<VectorXf> &GWASeffects2, const vector<Vector2f> &nGWASblock, const bool noscale){
    if (!wcorrBlocks.size()) {  // only print for the first time reading the data
        cout << "Reading eigenvectors from binary file and making W and Q matrices for bivariate analysis..." << endl;
    }
    if (!Gadget::directoryExist(dirname)) {
        throw("Error: cannot find the folder [" + dirname + "]");
    }
    
    vector<int>numSnpInRegion(numKeptLDBlocks);
    
    for(int i = 0; i < numKeptLDBlocks;i++){
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }

    unsigned subLdPlannedBiv = 0;
    vector<string> subLdMissingLdBlockIdsBiv;
    for (int i = 0; i < numKeptLDBlocks; ++i) {
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        int m = block->numSnpInBlock;
        if (m <= 1) continue;
        if (!recomputeEigen) continue;
        unsigned n_skip = 0;
        for (int j = 0; j < m; ++j)
            if (block->snpInfoVec[j]->skip) ++n_skip;
        if (n_skip == 0 || n_skip >= (unsigned)m) continue;
        if (blockLdmBinFileExists(dirname, block->ID) || blockEigenBinFileExists(dirname, block->ID))
            ++subLdPlannedBiv;
        else
            subLdMissingLdBlockIdsBiv.push_back(block->ID);
    }
    if (!subLdMissingLdBlockIdsBiv.empty() && !wcorrBlocks.size()) {
        cout << "Warning (bivariate): --recompute-eigen with --skip: partial skip in " << subLdMissingLdBlockIdsBiv.size()
             << " LD block(s); need block<ID>.ldm.bin or block<ID>.eigen.bin in [" << dirname << "]. Missing for block ID(s):";
        for (size_t t = 0; t < subLdMissingLdBlockIdsBiv.size(); ++t)
            cout << " " << subLdMissingLdBlockIdsBiv[t];
        cout << ". Falling back to full-block block*.eigen.bin for those blocks." << endl;
    }
    if (subLdPlannedBiv && !wcorrBlocks.size()) {
        cout << "For " << subLdPlannedBiv << " LD block(s) (bivariate), the LD submatrix of retained SNPs will be reconstructed (--recompute-eigen) and the GWAS effects of skipped SNPs will not be used." << endl;
    }

    eigenValLdBlock.resize(numLDBlocks);
    eigenVecLdBlock.resize(numLDBlocks);
    wcorrBlocks.resize(numKeptLDBlocks);
    numSnpsBlock.resize(numKeptLDBlocks);
    numEigenvalBlock.resize(numKeptLDBlocks);
    Qblocks.resize(numKeptLDBlocks);

#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        block->eigenColRemap.resize(0);
        block->subLdKeptLocalIdx.clear();
        block->qDiag.resize(0);

        VectorXf eigenValSub;
        MatrixXf eigenVecSub;
        VectorXf wcorrSub;
        MatrixXf Qsub;
        float sumPosSub = 0;
        VectorXi remapSub;
        vector<unsigned> keptSub;

        if (trySubLdEigenFromFullLdmBivariate(block, dirname, eigenCutoff, GWASeffects1[i], GWASeffects2[i], wcorrSub, Qsub, eigenValSub, eigenVecSub, sumPosSub, remapSub, keptSub, nullptr)) {
#pragma omp critical
            {
                cout << "LD block " << block->ID << " (bivariate): " << block->numSnpInBlock - (int)keptSub.size()
                     << " SNPs skipped, " << keptSub.size() << " retained." << endl;
            }
            eigenValLdBlock[i] = eigenValSub;
            eigenVecLdBlock[i] = eigenVecSub;
            wcorrBlocks[i] = wcorrSub;
            Qblocks[i] = Qsub;
            block->sumPosEigVal = sumPosSub;
            block->eigenvalues = eigenValSub;
            block->eigenColRemap = remapSub;
            block->subLdKeptLocalIdx = keptSub;
            numSnpsBlock[i] = Qblocks[i].cols();
            numEigenvalBlock[i] = Qblocks[i].rows();
            eigenVecLdBlock[i].resize(0,0);
            continue;
        }

        int32_t cur_m = 0;
        int32_t cur_k = 0;
        float sumPosEigVal = 0;
        float oldEigenCutoff =0;
        
        string infile = dirname + "/block" + block->ID + ".eigen.bin";
        FILE *fp = fopen(infile.c_str(), "rb");
        if(!fp){cout << "Error: can not open the file [" + infile + "] to read." << endl;}
        if(!fp){throw ("Error: can not open the file [" + infile + "] to read.");}

        // 1. marker number
        if(fread(&cur_m, sizeof(int32_t), 1, fp) != 1){
            throw("Read " + infile + " error (m)");
        }

        if(cur_m != numSnpInRegion[i]){
            throw("In LD block " + block->ID + ", inconsistent marker number to marker information in " + infile);
        }

        // 2. ncol of eigenVec (number of eigenvalues)
        if(fread(&cur_k, sizeof(int32_t), 1, fp) != 1){
            throw("In LD block " + block->ID + ", error about number of eigenvalues in  " + infile);
        }

        // 3. sum of all positive eigenvalues
        if(fread(&sumPosEigVal, sizeof(float), 1, fp) != 1){
            throw("In LD block " + block->ID + ", error about the sum of positive eigenvalues in " + infile);
        }

        // 4. eigenCutoff
        if(fread(&oldEigenCutoff, sizeof(float), 1, fp) != 1){
            throw("In LD block " + block->ID + ", error about eigen cutoff used in " + infile);
        }

        // 5. eigenvalues
        VectorXf lambda(cur_k);
        if(fread(lambda.data(), sizeof(float), cur_k, fp) != cur_k){
            throw("In LD block " + block->ID + ",size error about eigenvalues in " + infile);
        }

        // 6. eigenvector
        MatrixXf U(cur_m, cur_k);
        uint64_t nElements = (uint64_t)cur_m * (uint64_t)cur_k;
        if(fread(U.data(), sizeof(float), nElements, fp) != nElements){
            throw("In LD block " + block->ID + ",size error about eigenvectors in " + infile);
        }
        
        fclose(fp);
        
        if(oldEigenCutoff < eigenCutoff & i == 0){
            cout << "Warning: current proportion of variance in LD block is set as " + to_string(eigenCutoff)+ ". But the proportion of variance is set as "<< to_string(oldEigenCutoff) + " in "  + infile + ".\n";
        }
        
        if (eigenCutoff < oldEigenCutoff) {
            truncateEigenMatrix(sumPosEigVal, eigenCutoff, lambda, U, eigenValLdBlock[i], eigenVecLdBlock[i]);
        } else {
            eigenValLdBlock[i] = lambda;
            eigenVecLdBlock[i] = U;
        }
        block->sumPosEigVal = sumPosEigVal;
        block->eigenvalues = lambda;
                        
        // make w and Q for bivariate: concatenate w1 and w2
        VectorXf sqrtLambda = eigenValLdBlock[i].array().sqrt();
        VectorXf w1 = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (eigenVecLdBlock[i].transpose() * GWASeffects1[i]);
        VectorXf w2 = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (eigenVecLdBlock[i].transpose() * GWASeffects2[i]);
        
        // Concatenate w1 and w2: [w1; w2]
        wcorrBlocks[i].resize(w1.size() + w2.size());
        wcorrBlocks[i].head(w1.size()) = w1;
        wcorrBlocks[i].tail(w2.size()) = w2;
        
        Qblocks[i] = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
        
        numSnpsBlock[i] = Qblocks[i].cols();
        numEigenvalBlock[i] = Qblocks[i].rows();
    }

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < numKeptLDBlocks; ++i)
        cacheRealizedQdiag(keptLdBlockInfoVec[i], Qblocks[i]);
}

void Data::constructWandQ(const vector<VectorXf> &GWASeffects, const float nGWAS, const bool noscale) {
    wcorrBlocks.resize(numKeptLDBlocks);
    numSnpsBlock.resize(numKeptLDBlocks);
    numEigenvalBlock.resize(numKeptLDBlocks);
    Qblocks.resize(numKeptLDBlocks);

    for (unsigned i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
        // calculate wbcorr and Qblocks
        VectorXf sqrtLambda = eigenValLdBlock[i].array().sqrt();
        //cout << eigenVecLdBlock[i].transpose().rows() << " " << eigenVecLdBlock[i].transpose().cols() << " " << gwasEffectInBlock[i].size() << endl;
        wcorrBlocks[i] = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (eigenVecLdBlock[i].transpose() * GWASeffects[i] );
        // cout << "eigenVecLdBlock[i]: " << eigenVecLdBlock[i] << endl;
        // cout << "wcorrBlocks[i]: " << wcorrBlocks[i] << endl;
        // cout << "sqrtLambda: " << sqrtLambda << endl;
        // cout << gwasEffectInBlock[i] << endl;
//        MatrixXf tmpQblocks = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
//        MatrixDat matrixDat = MatrixDat(ldblock->snpNameVec, tmpQblocks);
//        // cout << "Qblock: " << endl;
//        // cout << matrixDat.values << endl;
//        Qblocks.push_back(matrixDat);
        Qblocks[i] = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
        
        if (noscale) {
            VectorXf Dsqrt(ldblock->numSnpInBlock);
            for (unsigned j=0; j<ldblock->numSnpInBlock; ++j) {
                SnpInfo *snp = ldblock->snpInfoVec[j];
                //Dsqrt[j] = sqrt(snp->twopq);
                Dsqrt[j] = snp->gwas_scalar;
            }
            Qblocks[i] = Qblocks[i] * Dsqrt.asDiagonal();
        }
        
        numSnpsBlock[i] = Qblocks[i].cols();
        numEigenvalBlock[i] = Qblocks[i].rows();
        cacheRealizedQdiag(ldblock, Qblocks[i]);
        
        eigenVecLdBlock[i].resize(0,0);
    }
    
    nGWASblock.resize(numKeptLDBlocks);
    for (unsigned i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
        nGWASblock[i] = nGWAS;
    }
}

void Data::outputWandQ(const string &dirname) {
    // output w and Q in text format
    if (!Gadget::directoryExist(dirname)){
        Gadget::createDirectory(dirname);
        cout << "  Created directory [" << dirname << "] to store w and Q matrices.\n\n";
    }
#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[i];
        
        string outTxtfile1 = dirname + "/block" + blockInfo->ID + ".w.txt";
        string outTxtfile2 = dirname + "/block" + blockInfo->ID + ".Q.txt";
        ofstream outtxt1(outTxtfile1.c_str());
        ofstream outtxt2(outTxtfile2.c_str());

        outtxt1 << wcorrBlocks[i] << endl;
        outtxt2 << Qblocks[i] << endl;
        
        outtxt1.close();
        outtxt2.close();
    }
}

void Data::scaleGwasEffects(){
    snp2pq.resize(numIncdSnps);
    b.resize(numIncdSnps);
    n.resize(numIncdSnps);
    se.resize(numIncdSnps);
    scalar.resize(numIncdSnps);
    tss.resize(numIncdSnps); // only used in SBayesC
    ZPy.resize(numIncdSnps);
    const char *gwasScalarLabel = "estimated by (1/sqrt(n*SE^2+b^2) from summary stats";
    if (gwasScalarMode == GwasScalarMode::Gwas2pq)
        gwasScalarLabel = "sqrt(2pq) based on allele frequency in GWAS)";
    else if (gwasScalarMode == GwasScalarMode::Ref2pq)
        gwasScalarLabel = "(sqrt(2pq) based on allele frequency in LD reference";
    cout << "GWAS effect scaling: " << gwasScalarLabel << "." << endl;

    SnpInfo *snp;
    unsigned ref2pqFallback = 0;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snp = incdSnpInfoVec[i];
        const float refAf = snp->af;
        snp->af = snp->gwas_af;
        snp2pq[i] = snp->twopq = 2.0f*snp->gwas_af*(1.0f-snp->gwas_af);
        if(snp2pq[i]==0) cout << "Error: SNP " << snp->ID << " af " << snp->af << " has 2pq = 0." << endl;
        b[i] = snp->gwas_b;
        n[i] = snp->gwas_n;
        se[i]= snp->gwas_se;
        double sc = 0.0;
        if (gwasScalarMode == GwasScalarMode::LeastSquares) {
            const double denom = double(n[i])*double(se[i])*double(se[i]) + double(b[i])*double(b[i]);
            if (!(denom > 0.0))
                throw("Error: SNP " + snp->ID + " has non-positive n*SE^2+b^2 for --gwas-scalar ls.");
            sc = 1.0/std::sqrt(denom);
        } else if (gwasScalarMode == GwasScalarMode::Ref2pq) {
            if (refAf > 0.0f && refAf < 1.0f)
                sc = std::sqrt(double(2.0f*refAf*(1.0f-refAf)));
            else {
                sc = std::sqrt(double(snp2pq[i]));
                ++ref2pqFallback;
            }
        } else {
            sc = std::sqrt(double(snp2pq[i]));
        }
        snp->gwas_scalar = sc;
        scalar[i] = snp->gwas_scalar;
//        b[i] = snp->gwas_b * sqrt(snp2pq[i]); // scale the marginal effect so that it's in per genotype SD unit
//        n[i] = snp->gwas_n;
//        se[i]= snp->gwas_se * sqrt(snp2pq[i]);
//        tss[i] = n[i]*(n[i]*se[i]*se[i] + b[i]*b[i]);
//        ZPy[i] = n[i]*b[i];
    }
    if (ref2pqFallback)
        cout << "Warning: --gwas-scalar ref2pq used GWAS sqrt(2pq) for " << ref2pqFallback << " SNPs (LD reference AF missing or not in (0,1))." << endl;

    // estimate phenotypic variance
    varySnp = snp2pq.array()*(n.array()*se.array().square()+b.array().square());
    VectorXf varpSrt = varySnp;
    std::sort(varpSrt.data(), varpSrt.data() + varpSrt.size());
    float obsVarPhenotypic = varpSrt[varpSrt.size()/2];
    varPhenotypic = 1.0;
    //cout << "varPhenotypic: " << varPhenotypic << endl;
    VectorXf nSrt = n;
    std::sort(nSrt.data(), nSrt.data() + nSrt.size());
    numKeptInds = nSrt[nSrt.size()/2]; // median
    
    // estimate per-SNP 2pq using the estimated phenotypic variance
    // for (unsigned i=0; i<numIncdSnps; ++i) {
    //     snp = incdSnpInfoVec[i];
    //     snp2pq[i] = snp->twopq = obsVarPhenotypic/(snp->gwas_n*se[i]*se[i]+b[i]*b[i]);       // NEW!
    // }
    
    // scale GWAS effects
    b.array() *= scalar.array();
    se.array() *= scalar.array();
    tss.array() = n.array()*(n.array()*se.array().square() + b.array().square());
    ZPy.array() = n.array()*b.array();
    
    // map to blocks
    if (numKeptLDBlocks) {
        gwasEffectInBlock.resize(numKeptLDBlocks);
        gwasPerSnpNinBlock.resize(numKeptLDBlocks);
        nGWASblock.resize(numKeptLDBlocks);
        for (unsigned i = 0; i < numKeptLDBlocks; i++){
            LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
            gwasEffectInBlock[i] = b(ldblock->block2GwasSnpVec);
            gwasPerSnpNinBlock[i] = n(ldblock->block2GwasSnpVec);
//            nGWASblock[i] = gwasPerSnpNinBlock[i].mean();
            VectorXf nBlockSrt = gwasPerSnpNinBlock[i];
            std::sort(nBlockSrt.data(), nBlockSrt.data() + nBlockSrt.size());
            nGWASblock[i] = nBlockSrt[nBlockSrt.size()/2]; // median
//            nGWASblock[i] = numKeptInds;
        }
    }
    
    // calculate ypy (total sum of squares) for summary statistics
    // ypy = mean phenotypic variance which is one after scaling * number of individuals
    ypy = numKeptInds;


    // output GWAS data
}

void Data::scaleBivariateGwasEffects(){
    // Scale GWAS effects for both traits (--gwas-scalar matches univariate Data::scaleGwasEffects)
    snp2pq.resize(numIncdSnps);
    b.resize(numIncdSnps);
    n.resize(numIncdSnps);
    se.resize(numIncdSnps);
    scalar.resize(numIncdSnps);
    b2.resize(numIncdSnps);
    se2.resize(numIncdSnps);
    scalar2.resize(numIncdSnps);
    n2.resize(numIncdSnps);

    const char *gwasScalarLabel = "ls (1/sqrt(n*SE^2+b^2) from summary stats)";
    if (gwasScalarMode == GwasScalarMode::Gwas2pq)
        gwasScalarLabel = "gwas2pq (sqrt(2pq) from GWAS AF)";
    else if (gwasScalarMode == GwasScalarMode::Ref2pq)
        gwasScalarLabel = "ref2pq (sqrt(2pq) from LD reference AF before GWAS AF overwrite)";
    cout << "GWAS effect scaling (bivariate): " << gwasScalarLabel << "." << endl;
    
    SnpInfo *snp;
    unsigned ref2pqFallback = 0;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snp = incdSnpInfoVec[i];
        const float refAf = snp->af;
        snp->af = snp->gwas_af;  // trait 1 GWAS AF on SnpInfo::af (same convention as univariate)
        snp2pq[i] = snp->twopq = 2.0f*snp->gwas_af*(1.0f-snp->gwas_af);
        if(snp2pq[i]==0) cout << "Error: SNP " << snp->ID << " af " << snp->af << " has 2pq = 0." << endl;
        
        b[i] = snp->gwas_b;
        n[i] = snp->gwas_n;
        se[i] = snp->gwas_se;
        b2[i] = snp->gwas_b2;
        n2[i] = snp->gwas_n2;
        se2[i] = snp->gwas_se2;

        double sc1 = 0.0, sc2 = 0.0;
        if (gwasScalarMode == GwasScalarMode::LeastSquares) {
            const double denom1 = double(n[i])*double(se[i])*double(se[i]) + double(b[i])*double(b[i]);
            if (!(denom1 > 0.0))
                throw("Error: SNP " + snp->ID + " trait 1: non-positive n*SE^2+b^2 for --gwas-scalar ls.");
            sc1 = 1.0/std::sqrt(denom1);
            const double denom2 = double(n2[i])*double(se2[i])*double(se2[i]) + double(b2[i])*double(b2[i]);
            if (!(denom2 > 0.0))
                throw("Error: SNP " + snp->ID + " trait 2: non-positive n*SE^2+b^2 for --gwas-scalar ls.");
            sc2 = 1.0/std::sqrt(denom2);
        } else if (gwasScalarMode == GwasScalarMode::Ref2pq) {
            if (refAf > 0.0f && refAf < 1.0f) {
                const double t = double(2.0f*refAf*(1.0f-refAf));
                sc1 = sc2 = std::sqrt(t);
            } else {
                sc1 = std::sqrt(double(snp2pq[i]));
                const float tpq2 = 2.0f*snp->gwas_af2*(1.0f-snp->gwas_af2);
                if (tpq2 <= 0.f)
                    throw("Error: SNP " + snp->ID + " trait 2 GWAS AF gives 2pq = 0; cannot complete ref2pq fallback.");
                sc2 = std::sqrt(double(tpq2));
                ++ref2pqFallback;
            }
        } else { // Gwas2pq
            sc1 = std::sqrt(double(snp2pq[i]));
            const float tpq2 = 2.0f*snp->gwas_af2*(1.0f-snp->gwas_af2);
            if (tpq2 <= 0.f)
                throw("Error: SNP " + snp->ID + " trait 2 GWAS AF gives 2pq = 0 for --gwas-scalar gwas2pq.");
            sc2 = std::sqrt(double(tpq2));
        }
        snp->gwas_scalar = sc1;
        scalar[i] = snp->gwas_scalar;
        snp->gwas_scalar2 = scalar2[i] = sc2;
    }
    if (ref2pqFallback)
        cout << "Warning: bivariate --gwas-scalar ref2pq used per-trait GWAS sqrt(2pq) for " << ref2pqFallback << " SNPs (LD reference AF missing or not in (0,1))." << endl;

    // estimate phenotypic variance for both traits (before scaling)
    varySnp  = snp2pq.array()*(n.array() *se.array().square() +b.array().square());
    varySnp2 = snp2pq.array()*(n2.array()*se2.array().square()+b2.array().square());
    VectorXf varpSrt = varySnp;
    std::sort(varpSrt.data(), varpSrt.data() + varpSrt.size());
    float obsVarPhenotypic = varpSrt[varpSrt.size()/2];
    varPhenotypic = 1.0;
    VectorXf nSrt = n;
    std::sort(nSrt.data(), nSrt.data() + nSrt.size());
    numKeptInds = nSrt[nSrt.size()/2]; // median sample size trait 1
    VectorXf n2Srt = n2;
    std::sort(n2Srt.data(), n2Srt.data() + n2Srt.size());
    numKeptInds2 = n2Srt[n2Srt.size()/2]; // median sample size trait 2
    
    // Implied 2pq from trait 1 summary stats (legacy bivariate LS path only)
    if (gwasScalarMode == GwasScalarMode::LeastSquares) {
        for (unsigned i=0; i<numIncdSnps; ++i) {
            snp = incdSnpInfoVec[i];
            snp2pq[i] = snp->twopq = obsVarPhenotypic/(snp->gwas_n*se[i]*se[i]+b[i]*b[i]);
        }
    }
    
    // scale GWAS effects for both traits
    b.array() *= scalar.array();
    se.array() *= scalar.array();
    b2.array() *= scalar2.array();
    se2.array() *= scalar2.array();
    
    // calculate ypy (total sum of squares) for summary statistics
    ypy = numKeptInds;
    ypy2 = numKeptInds2;
//    string outfile = "ma.txt";
//    ofstream out(outfile.c_str());
//    for (unsigned i=0; i<numIncdSnps; ++i) {
//        snp = incdSnpInfoVec[i];
//        out << snp->ID << " " << snp->a1 << " " << snp->a2 << " " << snp->af << " " << snp->gwas_b << " " << snp->gwas_se << " " << snp->gwas_pvalue << " " << snp->gwas_n << " " << snp2pq[i]*n[i] << " " << varySnp[i] << " " << b[i] << endl;
//    }
//    out.close();

}


void Data::resizeBlockLDmatrix(const string &inDirname, const string &outLDmatType, const string &includeSnpFile, const string &outDirname, const bool writeLdmTxt){
    readBlockLdmInfoFile(inDirname);
    readBlockLdmSnpInfoFile(inDirname);

    includeSnp(includeSnpFile);
    //includeMatchedSnp();
    //mapSnpsToBlocks();
    //includeMatchedBlocks();
    
    
//    ldBlockInfoVec = keptLdBlockInfoVec;
//    numLDBlocks = numKeptLDBlocks;
    
    struct stat sb;
    if (stat(inDirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        throw("Error: cannot find the folder [" + inDirname + "]");
    }
    
    if (stat(outDirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        string create_cmd = "mkdir " + outDirname;
        system(create_cmd.c_str());
        cout << "Created folder [" << outDirname << "] to store resized LD matrices." << endl;
    }

    vector<int> numSnpInRegion;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        LDBlockInfo *block = ldBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }
        
    keptLdBlockInfoVec = ldBlockInfoVec;

    cout << "Resizing block LD matrices... " << endl;

#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numLDBlocks; i++){
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[i];
                
        int32_t blockSize = blockInfo->numSnpInBlock;
        
        unsigned newBlockSize = 0;
        for (unsigned j=0; j<blockSize; ++j) {
            SnpInfo *snpj = blockInfo->snpInfoVec[j];
            if (snpj->included) {
                ++newBlockSize;
            }
        }
        
        if (!newBlockSize) continue;
        
        
        if (outLDmatType == "block") {
            
            MatrixXf ldm(blockSize, blockSize);
            readBlockLDmatrix(inDirname, blockInfo->ID, blockSize, ldm);

            //        cout << " read block " << blockInfo->ID << endl;
                        
            vector<string> newSnpNameVec;
            vector<SnpInfo*> newSnpInfoVec;
            
            MatrixXf ldm_resized(newBlockSize, newBlockSize);
            
            unsigned jj=0;
            for (unsigned j=0; j<blockSize; ++j) {
                SnpInfo *snpj = blockInfo->snpInfoVec[j];
                if (!snpj->included) continue;
                newSnpNameVec.push_back(snpj->ID);
                newSnpInfoVec.push_back(snpj);
                unsigned kk=0;
                for (unsigned k=0; k<blockSize; ++k) {
                    SnpInfo *snpk = blockInfo->snpInfoVec[k];
                    if (!snpk->included) continue;
                    ldm_resized(jj,kk) = ldm(j,k);
                    ++kk;
                }
                ++jj;
            }
            
            blockInfo->numSnpInBlock = newBlockSize;
            blockInfo->snpNameVec = newSnpNameVec;
            blockInfo->snpInfoVec = newSnpInfoVec;
            
            string outBinfile = outDirname + "/block" + blockInfo->ID + ".ldm.bin";
            FILE *outbin = fopen(outBinfile.c_str(), "wb");
            uint64_t nElements = (uint64_t) newBlockSize * (uint64_t) newBlockSize;
            fwrite(ldm_resized.data(), sizeof(float), nElements, outbin);
            fclose(outbin);
            
            if (writeLdmTxt) {
                string outTxtfile = outDirname + "/block" + blockInfo->ID + ".ldm.txt";
                ofstream outtxt;
                outtxt.open(outTxtfile.c_str());
                for (unsigned ii=0; ii<newBlockSize; ++ii){
                    for (unsigned jj=0; jj<newBlockSize; ++jj) {
                        outtxt << blockInfo->ID << "\t" << blockInfo->snpNameVec[ii] << "\t" << blockInfo->snpNameVec[jj] << "\t" << ldm_resized(ii,jj) << endl;
                    }
                }
                outtxt.close();
            }
            
        }
        else if (outLDmatType == "blockSparse") {
            string infile = inDirname + "/block" + blockInfo->ID + ".sparse.bin";
            FILE *fp = fopen(infile.c_str(), "rb");
            if(!fp){throw ("Error: can not open the file [" + infile + "] to read.");}
            
            SparseMatrix<float> ldmSp = Gadget::readSparseMatrixBinary(infile);
            
            //        cout << " read block " << blockInfo->ID << endl;
            
            vector<string> newSnpNameVec;
            vector<SnpInfo*> newSnpInfoVec;
            
            vector<Triplet<float> > tripletList;
            tripletList.reserve(newBlockSize * newBlockSize);  // preallocate for speed
            
            unsigned jj=0;
            for (unsigned j=0; j<blockSize; ++j) {
                SnpInfo *snpj = blockInfo->snpInfoVec[j];
                if (!snpj->included) continue;
                newSnpNameVec.push_back(snpj->ID);
                newSnpInfoVec.push_back(snpj);
                unsigned kk=0;
                for (unsigned k=0; k<blockSize; ++k) {
                    SnpInfo *snpk = blockInfo->snpInfoVec[k];
                    if (!snpk->included) continue;
                    float val = ldmSp.coeff(j, k);  // safely access sparse matrix element
                    if (val != 0.0f) {
                        tripletList.emplace_back(jj, kk, val);
                    }
                    ++kk;
                }
                ++jj;
            }
            
            SparseMatrix<float> ldmSp_resized(newBlockSize, newBlockSize);
            ldmSp_resized.setFromTriplets(tripletList.begin(), tripletList.end());
            
            blockInfo->numSnpInBlock = newBlockSize;
            blockInfo->snpNameVec = newSnpNameVec;
            blockInfo->snpInfoVec = newSnpInfoVec;
            
            string outBinfile = outDirname + "/block" + blockInfo->ID + ".sparse.bin";
            Gadget::writeSparseMatrixBinary(ldmSp_resized, outBinfile);
            
            if (writeLdmTxt) {
                string outTxtfile = outDirname + "/block" + blockInfo->ID + ".sparse.txt";
                Gadget::writeSparseMatrixToText(ldmSp_resized, outTxtfile);
            }
        }
            
        string outSnpfile = outDirname + "/block" + blockInfo->ID + ".snp.info";
        string outldmfile = outDirname + "/block" + blockInfo->ID + ".ldm.info";
        
        outputBlockLDmatrixInfo(*blockInfo, outSnpfile, outldmfile);
        
        if(!(i%1)) cout << " computed block " << blockInfo->ID << "\r" << flush;
    }
    
    if (outLDmatType == "block") {
        cout << "Written the LD matrix data into file [" << outDirname << "/block*.ldm.bin]." << endl;
        if (writeLdmTxt) cout << "Written the LD matrix data into file [" << outDirname << "/block*.ldm.txt]." << endl;
    } else if (outLDmatType == "blockSparse") {
        cout << "Written the LD matrix data into file [" << outDirname << "/block*.sparse.bin]." << endl;
        if (writeLdmTxt) cout << "Written the LD matrix data into file [" << outDirname << "/block*.sparse.txt]." << endl;
    }
}

void Data::outputBlockLDmatrixTxt(const string &dirname, const unsigned int block){
    readBlockLdmInfoFile(dirname, block);
    readBlockLdmSnpInfoFile(dirname, block);

    struct stat sb;
    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        throw("Error: cannot find the folder [" + dirname + "]");
    }

    vector<int> numSnpInRegion;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        LDBlockInfo *block = ldBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }
            
#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numLDBlocks; i++){
//        if (block && i != block - 1) continue;
        
        LDBlockInfo *blockInfo = ldBlockInfoVec[i];
        
        string infile = dirname + "/block" + blockInfo->ID + ".ldm.bin";
        FILE *fp = fopen(infile.c_str(), "rb");
        if(!fp){throw ("Error: can not open the file [" + infile + "] to read.");}
        
        string outTxtfile = dirname + "/block" + blockInfo->ID + ".ldm.txt";
        ofstream outtxt;
        outtxt.open(outTxtfile.c_str());
        
        int32_t blockSize = blockInfo->numSnpInBlock;
        
        MatrixXf ldm(blockSize, blockSize);
        readBlockLDmatrix(dirname, blockInfo->ID, blockSize, ldm);
        
//        uint64_t nElements = (uint64_t)blockSize * (uint64_t)blockSize;
//        
//        if(fread(ldm.data(), sizeof(float), nElements, fp) != nElements){
//            cout << "fread(U.data(), sizeof(float), nElements, fp): " << fread(ldm.data(), sizeof(float), nElements, fp) << endl;
//            cout << "nEle: " << nElements << " ldm.size: " << ldm.size() <<  " ldm.col: " << ldm.cols() << " row: " << ldm.rows() << endl;
//            throw("In LD block " + blockInfo->ID + ",size error in " + infile);
//            // cout << "Read " << svdLDfile << " error (U)" << endl;
//            // throw("read file error");
//        }
        
//        for (unsigned ii=0; ii<blockSize; ++ii){
//            for (unsigned jj=0; jj<blockSize; ++jj) {
//                outtxt << blockInfo->ID << "\t" << blockInfo->snpNameVec[ii] << "\t" << blockInfo->snpNameVec[jj] << "\t" << ldm(ii,jj) << endl;
//            }
//        }

        //        outtxt << ldm << endl;
        
        for (unsigned row = 0; row < blockSize; ++row) {
            for (unsigned col = 0; col <= row; ++col) {
                outtxt << ldm(row, col) << "\t";
            }
            outtxt << endl;
        }

        outtxt.close();

    }
    
    if (block) {
        cout << "Written the LD matrix data into file [" << dirname << "/block" << block << ".ldm.txt]." << endl;
    } else {
        cout << "Written the LD matrix data into file [" << dirname << "/block*.ldm.txt]." << endl;
    }

}

void Data::resizeBlockLDmatrixAndDoEigenDecomposition(const string &inDirname, const float eigenCutoff, const float rsqThreshold, const string &outDirname, const bool writeLdmTxt){
    //readBlockLdmInfoFile(inDirname + "/ldm.info");
    //readBlockLdmSnpInfoFile(inDirname + "/snp.info");

    struct stat sb;
    if (stat(inDirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        throw("Error: cannot find the folder [" + inDirname + "]");
    }
    
    if (stat(outDirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        string create_cmd = "mkdir " + outDirname;
        system(create_cmd.c_str());
        cout << "Created folder [" << outDirname << "] to store resized LD matrices." << endl;
    }

    vector<int> numSnpInRegion;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        LDBlockInfo *block = ldBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }


    keptLdBlockInfoVec.resize(0);

    cout << "Resizing block LD matrices... " << endl;

#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numLDBlocks; i++){
        
        LDBlockInfo *blockInfo = ldBlockInfoVec[i];
                
        int32_t blockSize = blockInfo->numSnpInBlock;
        
        unsigned newBlockSize = 0;
        for (unsigned j=0; j<blockSize; ++j) {
            SnpInfo *snpj = blockInfo->snpInfoVec[j];
            if (snpj->included) {
                ++newBlockSize;
            }
        }
        
        if (!newBlockSize) {
            blockInfo->kept = false;
            continue;
        } else {
            blockInfo->kept = true;
            keptLdBlockInfoVec.push_back(blockInfo);
        }
        
        MatrixXf ldm(blockSize, blockSize);
        readBlockLDmatrix(inDirname, blockInfo->ID, blockSize, ldm);
        
        //        cout << " read block " << blockInfo->ID << endl;
        
        vector<string> newSnpNameVec;
        vector<SnpInfo*> newSnpInfoVec;
        
        MatrixXf ldm_resized(newBlockSize, newBlockSize);
        
        unsigned jj=0;
        for (unsigned j=0; j<blockSize; ++j) {
            SnpInfo *snpj = blockInfo->snpInfoVec[j];
            if (!snpj->included) continue;
            newSnpNameVec.push_back(snpj->ID);
            newSnpInfoVec.push_back(snpj);
            unsigned kk=0;
            for (unsigned k=0; k<blockSize; ++k) {
                SnpInfo *snpk = blockInfo->snpInfoVec[k];
                if (!snpk->included) continue;
                ldm_resized(jj,kk) = ldm(j,k);
                ++kk;
            }
            ++jj;
        }
        
        blockInfo->numSnpInBlock = newBlockSize;
        blockInfo->snpNameVec = newSnpNameVec;
        blockInfo->snpInfoVec = newSnpInfoVec;
        
        MatrixXf eigenVec;
        VectorXf eigenVal;
        float sumPosEigVal;  // sum of all positive eigenvalues
        // cout << "rval: " << rval << endl;
        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        eigenDecomposition(ldm_resized, eigenCutoff, eigenVal, eigenVec, sumPosEigVal);

        blockInfo->sumPosEigVal = sumPosEigVal;
        
        int32_t numEigenValue = eigenVal.size();
        int32_t numSnpInBlock = newBlockSize;


        
        string outBinfile = outDirname + "/block" + blockInfo->ID + ".eigen.bin";
        FILE *outbin = fopen(outBinfile.c_str(), "wb");
                
        // save summary
        // 1, the number of SNPs in the block
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, outbin);
        // 2, the number of eigenvalues at with the given cutoff
        fwrite(&numEigenValue, sizeof(int32_t), 1, outbin);
        // 3. sum of all the positive eigenvalues
        fwrite(&sumPosEigVal, sizeof(float), 1, outbin);
        // 4. eigenvalue cutoff based on the proportion of variance explained in LD
        fwrite(&eigenCutoff, sizeof(float), 1, outbin);
        // 5. the selected eigenvalues
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, outbin);
        // 6. the selected eigenvector;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, outbin);
        
        fclose(outbin);


        if (writeLdmTxt) {
            string outTxtfile = outDirname + "/block" + blockInfo->ID + ".eigen.txt";
            ofstream outtxt;
            outtxt.open(outTxtfile.c_str());
            outtxt << "Block " << blockInfo->ID << endl;
            outtxt << "numSnps " << numSnpInBlock << endl;
            outtxt << "numEigenvalues " << numEigenValue << endl;
            outtxt << "SumPositiveEigenvalues " << sumPosEigVal << endl;
            outtxt << "EigenCutoff " << eigenCutoff << endl;
            outtxt << "Eigenvalues\n" << eigenVal.transpose() << endl;
            outtxt << "Eigenvectors\n" << eigenVec << endl;
            outtxt << endl;
            outtxt.close();
        }
        
            
        string outSnpfile = outDirname + "/block" + blockInfo->ID + ".snp.info";
        string outldmfile = outDirname + "/block" + blockInfo->ID + ".ldm.info";
        
        outputBlockLDmatrixInfo(*blockInfo, outSnpfile, outldmfile);
        
        outputLDfriends(ldm_resized, blockInfo, outDirname, writeLdmTxt);
        
        if(!(i%1)) cout << " computed block " << blockInfo->ID << "\r" << flush;
    }
    
    
    numKeptLDBlocks = keptLdBlockInfoVec.size();
    
    if (numKeptLDBlocks == 1) {
        cout << "Written the LD matrix data into file [" << outDirname << "/block" << keptLdBlockInfoVec[0]->ID << ".eigen.bin]." << endl;
        if (writeLdmTxt) cout << "Written the LD matrix data into file [" << outDirname << "/block" << keptLdBlockInfoVec[0]->ID << ".eigen.txt]." << endl;
        cout << "Written the high pairwise LD correlations (rsq>0.5) into [" << outDirname << "/block" << keptLdBlockInfoVec[0]->ID << ".rsq0.5.pwld.gz]." << endl;
    } else {
        cout << "Written the LD matrix data into file [" << outDirname << "/block*.eigen.bin]." << endl;
        if (writeLdmTxt) cout << "Written the LD matrix data into file [" << outDirname << "/block*.eigen.txt]." << endl;
        mergeLdmInfo("block", outDirname, true);
    }
}


void Data::readBlockLDmatrix(const string &dirname, const string &blockID, const int32_t blockSize, MatrixXf &ldm) {
    string infile = dirname + "/block" + blockID + ".ldm.bin";
    int quantBits = 0;

    // Get file size
    struct stat stat_buf;
    if (stat(infile.c_str(), &stat_buf) != 0) {
        string q16file = dirname + "/block" + blockID + ".ldm.q16.bin";
        string q8file = dirname + "/block" + blockID + ".ldm.q8.bin";
        if (stat(q16file.c_str(), &stat_buf) == 0) {
            infile = q16file;
            quantBits = 16;
        } else if (stat(q8file.c_str(), &stat_buf) == 0) {
            infile = q8file;
            quantBits = 8;
        } else {
            throw("Error: cannot stat the file [" + infile + "] or quantized LD files [" + q16file + ", " + q8file + "]");
        }
    }
    size_t file_bytes = stat_buf.st_size;

    // Expected sizes
    size_t n_full = (size_t)blockSize * blockSize;
    size_t n_lower_tri = (size_t)blockSize * (blockSize + 1) / 2;

    FILE *fp = fopen(infile.c_str(), "rb");
    if (!fp) {
        throw("Error: cannot open the file [" + infile + "] to read.");
    }

    ldm.setZero(blockSize, blockSize);

    if (quantBits) {
        const size_t expected_bytes = n_lower_tri * (quantBits == 8 ? sizeof(int8_t) : sizeof(int16_t));
        if (file_bytes != expected_bytes) {
            fclose(fp);
            throw("Unrecognized quantized LD matrix format in " + infile +
                  ": got " + to_string((long long)file_bytes) + " bytes, but expected " +
                  to_string((long long)expected_bytes) + ".");
        }
        const float invBound = quantBits == 8 ? (1.0f / 127.0f) : (1.0f / 32767.0f);
        if (quantBits == 8) {
            vector<int8_t> buffer(n_lower_tri);
            if (fread(buffer.data(), sizeof(int8_t), n_lower_tri, fp) != n_lower_tri) {
                fclose(fp);
                throw("Error reading quantized LD matrix from " + infile);
            }
            size_t idx = 0;
            for (int i = 0; i < blockSize; ++i) {
                for (int j = 0; j <= i; ++j) {
                    float value = float(buffer[idx]) * invBound;
                    ldm(i, j) = value;
                    if (i != j) ldm(j, i) = value;
                    ++idx;
                }
            }
        } else {
            vector<int16_t> buffer(n_lower_tri);
            if (fread(buffer.data(), sizeof(int16_t), n_lower_tri, fp) != n_lower_tri) {
                fclose(fp);
                throw("Error reading quantized LD matrix from " + infile);
            }
            size_t idx = 0;
            for (int i = 0; i < blockSize; ++i) {
                for (int j = 0; j <= i; ++j) {
                    float value = float(buffer[idx]) * invBound;
                    ldm(i, j) = value;
                    if (i != j) ldm(j, i) = value;
                    ++idx;
                }
            }
        }
    } else {
        if (file_bytes % sizeof(float) != 0) {
            fclose(fp);
            throw("Unrecognized matrix format in " + infile + ": file size is not a multiple of float size.");
        }
        size_t n_floats = file_bytes / sizeof(float);

        if (n_floats == n_full) {
        // Full square matrix
            if (fread(ldm.data(), sizeof(float), n_full, fp) != n_full) {
                throw("Error reading full matrix from " + infile);
            }
        } else if (n_floats == n_lower_tri) {
            // Lower triangular matrix
            vector<float> buffer(n_lower_tri);
            if (fread(buffer.data(), sizeof(float), n_lower_tri, fp) != n_lower_tri) {
                throw("Error reading lower-triangular matrix from " + infile);
            }

            size_t idx = 0;
            for (int i = 0; i < blockSize; ++i) {
                for (int j = 0; j <= i; ++j) {
                    ldm(i, j) = buffer[idx];
                    if (i != j) ldm(j, i) = buffer[idx];
                    ++idx;
                }
            }
        } else {
            fclose(fp);
            throw("Unrecognized matrix format in " + infile +
                  ": got " + to_string(n_floats) + " elements, but expected " +
                  to_string(n_full) + " (full) or " + to_string(n_lower_tri) + " (lower-triangular).");
        }
    }

    fclose(fp);
    normaliseDenseCorrelation(ldm);
}

void Data::convertToBlockTriangularMatrix(const string &inDirname, const bool writeLdmTxt, const string &outDirname){
    readBlockLdmInfoFile(inDirname);
    readBlockLdmSnpInfoFile(inDirname);

    struct stat sb;
    if (stat(outDirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        string create_cmd = "mkdir " + outDirname;
        system(create_cmd.c_str());
        cout << "Created folder [" << outDirname << "] to store LD matrices." << endl;
    }

#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numLDBlocks; i++){
        LDBlockInfo *blockInfo = ldBlockInfoVec[i];
                
        int32_t blockSize = blockInfo->numSnpInBlock;
        
        MatrixXf ldm(blockSize, blockSize);
        readBlockLDmatrix(inDirname, blockInfo->ID, blockSize, ldm);

        string outBinfile = outDirname + "/block" + blockInfo->ID + ".ldm.bin";
        FILE *outbin = fopen(outBinfile.c_str(), "wb");
        ofstream outtxt;
        string outTxtfile;
        if (writeLdmTxt) {
            outTxtfile = outDirname + "/block" + blockInfo->ID + ".ldm.txt";
            outtxt.open(outTxtfile.c_str());
        }

        // save the lower trangular matrix only
        unsigned numSnpInBlock = blockInfo->numSnpInBlock;
        for (unsigned row = 0; row < numSnpInBlock; ++row) {
            for (unsigned col = 0; col <= row; ++col) {
                float value = ldm(row, col);
                fwrite(&value, sizeof(float), 1, outbin);
            }
        }
        
        if (writeLdmTxt) {
            for (unsigned row = 0; row < numSnpInBlock; ++row) {
                for (unsigned col = 0; col <= row; ++col) {
                    outtxt << ldm(row, col) << "\t";
                }
                outtxt << endl;
            }
        }
        
        fclose(outbin);
        if (writeLdmTxt) outtxt.close();
        
        string outSnpfile = outDirname + "/block" + blockInfo->ID + ".snp.info";
        string outldmfile = outDirname + "/block" + blockInfo->ID + ".ldm.info";

        outputBlockLDmatrixInfo(*blockInfo, outSnpfile, outldmfile);
        
        outputLDfriends(ldm, blockInfo, outDirname, writeLdmTxt);
       
        if(!(i%1)) cout << " computed block " << blockInfo->ID << "\r" << flush;

    }
    
    cout << "Written the LD matrix into folder [" << outDirname << "/block*.ldm.bin]." << endl;
    if (writeLdmTxt) cout << "Written the LD matrix into text file [" << outDirname << "/block*.ldm.txt]." << endl;
    
    mergeLdmInfo("block", outDirname, true);
}


void Data::outputLDfriends(const MatrixXf &ldm, const LDBlockInfo *blockInfo, const string &outDirname, const bool writeText){
    // find LD friends for each SNP with LD rsq > threshold
    float rsqThreshold = 0.5;
    int blockSize = ldm.rows();
    vector<vector<int> > ldSnpIdxBlk(blockSize);
    vector<vector<float> > ldcorBlk(blockSize);
    for(unsigned j=0; j < blockSize; j++){
        SnpInfo *snp = blockInfo->snpInfoVec[j];
        for(unsigned k=0; k < j; k++){
            float rsq = ldm(j,k)*ldm(j,k);
            if(rsq > rsqThreshold){
                ldSnpIdxBlk[j].push_back(k);
                ldcorBlk[j].push_back(ldm(j,k));
            }
        }
    }
    
    string pwldfilename = outDirname + "/block" + blockInfo->ID + ".rsq0.5.pwld";  // pairwise LD file
    string gzPwldFilename = pwldfilename + ".gz";
    gzFile gzout = gzopen(gzPwldFilename.c_str(), "wb");
    if (!gzout) throw("Error: cannot open file [" + gzPwldFilename + "] to write.");
    gzprintf(gzout, "%12s %12s %12s\n", "SNP1", "SNP2", "LDcorrelation");

    remove(pwldfilename.c_str());
    remove((pwldfilename + ".bin").c_str());
    
    for (unsigned j=0; j<blockInfo->numSnpInBlock; ++j) {
        unsigned numLDfrd = ldSnpIdxBlk[j].size();
        if (!numLDfrd) continue;
        SnpInfo *snpj = blockInfo->snpInfoVec[j];
        for (unsigned k=0; k<numLDfrd; ++k) {
            SnpInfo *snpk = blockInfo->snpInfoVec[ldSnpIdxBlk[j][k]];
            
            gzprintf(gzout, "%12s %12s %12.6f\n", snpj->ID.c_str(), snpk->ID.c_str(), ldcorBlk[j][k]);
        }
    }
    gzclose(gzout);

}
