//
//  gadgets.cpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "gadgets.hpp"

void Gadget::Tokenizer::getTokens(const string &str, const string &sep){
    clear();
    string::size_type begidx,endidx;
    begidx = str.find_first_not_of(sep);
    while (begidx != string::npos) {
        endidx = str.find_first_of(sep,begidx);
        if (endidx == string::npos) endidx = str.length();
        push_back(str.substr(begidx,endidx - begidx));
        begidx = str.find_first_not_of(sep,endidx);
    }
}

int Gadget::Tokenizer::getIndex(const string &str, const bool err){
    for (unsigned i=0; i<size(); i++){
        if((*this)[i]==str){
            return i;
        }
    }
    if (err) {
        throw ("Error: can not find '" + str + "'.");
    } else {
        return -1;
    }
}

void Gadget::Timer::setTime(){
    prev = curr = time(0);
}

time_t Gadget::Timer::getTime(){
    return curr = time(0);
}

time_t Gadget::Timer::getElapse(){
    return curr - prev;
}

string Gadget::Timer::format(const time_t time){
    return to_string((long long)(time/3600)) + ":" + to_string((long long)((time % 3600)/60)) + ":" + to_string((long long)(time % 60));
}

string Gadget::Timer::getDate(){
    return ctime(&curr);
}

void Gadget::Timer::printElapse(){
    getTime();
    cout << "Time elapse: " << format(getElapse()) << endl;
}

string Gadget::getFileName(const string &file){
    size_t start = file.rfind('/');
    size_t end   = file.rfind('.');
    start = start==string::npos ? 0 : start+1;
    return file.substr(start, end-start);
}

string Gadget::getFileSuffix(const string &file){
    size_t start = file.rfind('.');
    return file.substr(start);
}

void Gadget::fileExist(const string &filename){
    ifstream file(filename.c_str());
    if(!file) throw("Error: can not open the file ["+filename+"] to read.");
}

bool Gadget::directoryExist(const string& dirname)
{
    struct stat info;

    if (stat(dirname.c_str(), &info) != 0)
        return false;
    else if (info.st_mode & S_IFDIR)
        return true;
    else
        return false;
}

bool Gadget::createDirectory(const string& dirname)
{
    int status = mkdir(dirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    if (status == -1)
        return false;
    else
        return true;
}

float Gadget::calcMean(const VectorXf &vec){
    VectorXd vec_double = vec.cast<double>();
    return vec_double.mean();
}

float Gadget::calcVariance(const VectorXf &vec){
    VectorXd vec_double = vec.cast<double>();
    return (vec_double.array() - vec_double.mean()).square().sum()/vec_double.size();
}

float Gadget::calcCovariance(const VectorXf &vec1, const VectorXf &vec2){
    if (vec1.size() != vec2.size()) {
        throw("Error: Gadget::calcCovariance: the two vectors have different sizes.");
    }
    VectorXd vec1_double = vec1.cast<double>();
    VectorXd vec2_double = vec2.cast<double>();
    return (vec1_double.array()-vec1_double.mean()).cwiseProduct(vec2_double.array()-vec2_double.mean()).sum()/vec1_double.size();
}

float Gadget::calcCorrelation(const VectorXf &vec1, const VectorXf &vec2){
    float cov = calcCovariance(vec1, vec2);
    float var1 = calcVariance(vec1);
    float var2 = calcVariance(vec2);
    return cov/sqrt(var1*var2);
}

float Gadget::calcRegression(const VectorXf &y, const VectorXf &x){
    float cov = calcCovariance(y, x);
    float varx = calcVariance(x);
    return cov/varx;
}

float Gadget::findMedian(const VectorXf &vec){
    VectorXf tmp = vec;
    std::sort(tmp.data(), tmp.data() + tmp.size());
    // return tmp[tmp.size()/2];
    return tmp.size() % 2 == 0 ? tmp.segment( (tmp.size()-2)/2, 2 ).mean() : tmp( tmp.size()/2 );
}

vector<int> Gadget::shuffle_index(const int start, const int end){
    vector<int> vec;
    for (unsigned i = start; i <= end; i++) {
        vec.push_back(i);
    }
    
    Gadget::shuffle_vector(vec);

    return vec;
}

void Gadget::shuffle_vector(vector<int> &vec){
    // Get a random integer using Boost random generator (which has been seeded)
    constexpr int max_integer = std::numeric_limits<int>::max();
    int random_integer = static_cast<int>(Stat::ranf() * max_integer);
    
    // Create a thread-local random number generator
    thread_local std::mt19937 rng(random_integer);

    // Shuffle using the thread-local RNG
    std::shuffle(vec.begin(), vec.end(), rng);
}

void Gadget::removeSecondElement(VectorXf &vec){
    // Erase the second element (index 1)
    vec.segment(1, vec.size() - 1) = vec.segment(2, vec.size() - 2);

    // Resize the vector to one less element
    vec.conservativeResize(vec.size() - 1);
}

void Gadget::writeSparseMatrixBinary(const SparseMatrix<float>& mat, const std::string& filename) {
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Failed to open file for writing.");
    }

    // Ensure matrix is compressed
    SparseMatrix<float> m = mat;
    m.makeCompressed();

    int rows = m.rows();
    int cols = m.cols();
    int nnz  = m.nonZeros();  // number of non-zero entries

    // Write metadata
    out.write(reinterpret_cast<const char*>(&rows), sizeof(int));
    out.write(reinterpret_cast<const char*>(&cols), sizeof(int));
    out.write(reinterpret_cast<const char*>(&nnz),  sizeof(int));

    // Write outer index (size: cols+1 for col-major, rows+1 for row-major)
    int outerSize = m.outerSize();
    out.write(reinterpret_cast<const char*>(m.outerIndexPtr()), sizeof(int) * (outerSize + 1));

    // Write inner indices (row or column indices of non-zeros)
    out.write(reinterpret_cast<const char*>(m.innerIndexPtr()), sizeof(int) * nnz);

    // Write values
    out.write(reinterpret_cast<const char*>(m.valuePtr()), sizeof(float) * nnz);

    out.close();
}

SparseMatrix<float> Gadget::readSparseMatrixBinary(const std::string& filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) throw std::runtime_error("Failed to open file.");

    int rows, cols, nnz;
    in.read(reinterpret_cast<char*>(&rows), sizeof(int));
    in.read(reinterpret_cast<char*>(&cols), sizeof(int));
    in.read(reinterpret_cast<char*>(&nnz),  sizeof(int));

    std::vector<int> outer(rows + 1);
    std::vector<int> inner(nnz);
    std::vector<float> values(nnz);

    in.read(reinterpret_cast<char*>(outer.data()), sizeof(int) * (rows + 1));
    in.read(reinterpret_cast<char*>(inner.data()), sizeof(int) * nnz);
    in.read(reinterpret_cast<char*>(values.data()), sizeof(float) * nnz);

    SparseMatrix<float> mat(rows, cols);
    mat.reserve(nnz);

    std::copy(outer.begin(), outer.end(), const_cast<int*>(mat.outerIndexPtr()));
    std::copy(inner.begin(), inner.end(), const_cast<int*>(mat.innerIndexPtr()));
    std::copy(values.begin(), values.end(), const_cast<float*>(mat.valuePtr()));

    mat.finalize();  // update internal structure

    return mat;
}

void Gadget::writeSparseMatrixToText(const SparseMatrix<float>& mat, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        throw std::runtime_error("Failed to open file for writing.");
    }

    out << "# Rows: " << mat.rows() << ", Cols: " << mat.cols() << ", Nonzeros: " << mat.nonZeros() << "\n";
    out << "# Format: row_index col_index value\n";

    for (int k = 0; k < mat.outerSize(); ++k) {
        for (SparseMatrix<float>::InnerIterator it(mat, k); it; ++it) {
            out << it.row() << " " << it.col() << " " << it.value() << "\n";
        }
    }

    out.close();
}
