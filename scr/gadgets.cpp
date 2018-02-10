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

float Gadget::calcVariance(const VectorXf &vec){
    return (vec.array() - vec.mean()).square().sum()/vec.size();
}

float Gadget::calcCovariance(const VectorXf &vec1, const VectorXf &vec2){
    if (vec1.size() != vec2.size()) {
        throw("Error: Gadget::calcCovariance: the two vectors have different sizes.");
    }
    return (vec1.array()-vec1.mean()).cwiseProduct(vec2.array()-vec2.mean()).sum()/vec1.size();
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
