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

int Gadget::Tokenizer::getIndex(const string &str){
    for (unsigned i=0; i<size(); i++){
        if((*this)[i]==str){
            return i;
        }
    }
    return -1;
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
    return tmp[tmp.size()/2];
}

