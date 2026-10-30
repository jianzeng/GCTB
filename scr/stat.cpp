//
//  stat.cpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "stat.hpp"
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

unsigned Stat::masterSeed = 11415u;
unsigned Stat::mcmcIterationForRng = 0u;
unsigned Stat::mcmcChainForRng = 0u;

thread_local unsigned long long Stat::rngDrawCount = 0;

#ifdef GCTB_COUNT_RNG_DRAWS
Stat::CountedMt19937::CountedMt19937() : impl_() {}
Stat::CountedMt19937::CountedMt19937(result_type s) : impl_(s) {}
void Stat::CountedMt19937::seed(result_type s) { impl_.seed(s); }
Stat::CountedMt19937::result_type Stat::CountedMt19937::operator()() {
    ++rngDrawCount;
    return impl_();
}
#endif

thread_local Stat::random_engine Stat::engine;
thread_local Stat::uniform01_generator Stat::ranf(Stat::engine, Stat::uniform_01());
thread_local Stat::normal_generator Stat::snorm(Stat::engine, Stat::normal_distribution(0,1));

namespace {
    std::vector<unsigned> parallelChainByOuterThread;
    bool parallelChainRngActive = false;
    thread_local bool chainOverrideActive = false;
    thread_local unsigned chainOverride = 0u;

    uint64_t mixRngKey(uint64_t value) {
        value += UINT64_C(0x9E3779B97F4A7C15);
        value = (value ^ (value >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
        value = (value ^ (value >> 27)) * UINT64_C(0x94D049BB133111EB);
        return value ^ (value >> 31);
    }

    Stat::random_engine::result_type contextSeed(unsigned chain, unsigned iteration,
                                                  unsigned taskKey, unsigned salt) {
        uint64_t key = mixRngKey(Stat::masterSeed);
        key = mixRngKey(key ^ (uint64_t(chain) << 32) ^ iteration);
        key = mixRngKey(key ^ (uint64_t(taskKey) << 32) ^ salt);
        return static_cast<Stat::random_engine::result_type>(key ^ (key >> 32));
    }
}

void Stat::seedEngine(const int seed){
    rngDrawCount = 0;
    if (seed) {
        masterSeed = (unsigned)seed;
        srand(seed);
        engine.seed(seed);
    } else {
        masterSeed = (unsigned)time(NULL);
        srand((int)masterSeed);
        engine.seed(masterSeed);
    }
    snorm.distribution() = normal_distribution(0,1);
    mcmcChainForRng = 0u;
    mcmcIterationForRng = 0u;
}

void Stat::setMcmcIteration(unsigned iter) {
    mcmcIterationForRng = iter;
}

void Stat::setMcmcContext(unsigned chain, unsigned iter) {
    mcmcChainForRng = chain;
    mcmcIterationForRng = iter;
}

unsigned Stat::currentMcmcChain() {
    if (chainOverrideActive) return chainOverride;
#ifdef _OPENMP
    if (parallelChainRngActive && omp_get_level() > 1) {
        int outerThread = omp_get_ancestor_thread_num(1);
        if (outerThread >= 0 && (unsigned)outerThread < parallelChainByOuterThread.size())
            return parallelChainByOuterThread[outerThread];
    }
#endif
    return mcmcChainForRng;
}

void Stat::prepareParallelChainRng(unsigned numOuterThreads) {
    parallelChainByOuterThread.assign(numOuterThreads ? numOuterThreads : 1u, mcmcChainForRng);
    parallelChainRngActive = true;
}

void Stat::finishParallelChainRng() {
    parallelChainRngActive = false;
}

Stat::random_engine::result_type Stat::parallelTaskSeed(unsigned taskKey, unsigned salt) {
    return contextSeed(currentMcmcChain(), mcmcIterationForRng, taskKey, salt);
}

void Stat::seedEngineForParallelTask(unsigned taskKey, unsigned salt) {
    engine.seed(parallelTaskSeed(taskKey, salt));
    snorm.distribution() = normal_distribution(0,1);
}

void Stat::seedEngineForParallelTaskIfNeeded(unsigned taskKey, unsigned salt) {
#ifdef _OPENMP
    if (omp_get_max_threads() <= 1) return;
#endif
    seedEngineForParallelTask(taskKey, salt);
}

Stat::McmcChainRngScope::McmcChainRngScope(unsigned chain) :
    savedEngine(engine),
    savedNormalDistribution(snorm.distribution()),
    savedChainOverrideActive(chainOverrideActive),
    savedChainOverride(chainOverride) {
    chainOverrideActive = true;
    chainOverride = chain;
#ifdef _OPENMP
    if (parallelChainRngActive && omp_get_level() == 1) {
        int outerThread = omp_get_thread_num();
        if (outerThread >= 0 && (unsigned)outerThread < parallelChainByOuterThread.size())
            parallelChainByOuterThread[outerThread] = chain;
    }
#endif
    engine.seed(contextSeed(chain, mcmcIterationForRng, 0xFFFFFFFFu, 0x43484149u));
    snorm.distribution() = normal_distribution(0,1);
}

Stat::McmcChainRngScope::~McmcChainRngScope() {
    engine = savedEngine;
    snorm.distribution() = savedNormalDistribution;
    chainOverrideActive = savedChainOverrideActive;
    chainOverride = savedChainOverride;
}

Stat::ParallelTaskRngScope::ParallelTaskRngScope(unsigned taskKey, unsigned salt) :
    savedEngine(engine),
    savedNormalDistribution(snorm.distribution()) {
    seedEngineForParallelTask(taskKey, salt);
}

Stat::ParallelTaskRngScope::~ParallelTaskRngScope() {
    engine = savedEngine;
    snorm.distribution() = savedNormalDistribution;
}

void Stat::logRngCheckpoint(const char *label) {
    if (!std::getenv("GCTB_DEBUG_RNG")) return;
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    std::cerr << "[GCTB_DEBUG_RNG] " << (label ? label : "(null)")
              << " omp_tid=" << tid << " rng_draws=";
#ifdef GCTB_COUNT_RNG_DRAWS
    std::cerr << rngDrawCount;
#else
    std::cerr << "disabled (rebuild with -DGCTB_COUNT_RNG_DRAWS)";
#endif
    std::cerr << std::endl;
}

double Stat::Normal::sample(const double mean, const double variance){
    return mean + snorm()*sqrtf(variance);
}

double Stat::Normal::cdf_01(const double value){
    // Phi(x) = erfc(-x/sqrt(2))/2. NaN is passed to Boost, which raises its domain error.
    if (std::isnan(value)) return cdf(d, value);
    return std::erfc(-value/1.41421356237309504880)/2;
}

namespace {
    // Inverse standard normal CDF for 0 < p < 1: Wichura (1988), Algorithm AS 241 (PPND16),
    // Applied Statistics 37:477-484; accurate to about 16 significant digits.
    double normalQuantileAS241(const double p) {
        const double q=p-0.5;
        if (fabs(q)<=0.425) {
            const double r=0.180625-q*q;
            return q*(((((((r*2509.0809287301226727+33430.575583588128105)*r+67265.770927008700853)*r+45921.953931549871457)*r+13731.693765509461125)*r+1971.5909503065514427)*r+133.14166789178437745)*r+3.387132872796366608)
                   /(((((((r*5226.495278852545925+28729.085735721942674)*r+39307.89580009271061)*r+21213.794301586595867)*r+5394.1960214247511077)*r+687.1870074920579083)*r+42.313330701600911252)*r+1.0);
        }
        double r=sqrt(-log(q<0 ? p : 1.0-p)), quantileResult;
        if (r<=5.0) {
            r-=1.6;
            quantileResult=(((((((r*7.7454501427834140764e-4+0.0227238449892691845833)*r+0.24178072517745061177)*r+1.27045825245236838258)*r+3.64784832476320460504)*r+5.7694972214606914055)*r+4.6303378461565452959)*r+1.42343711074968357734)
                 /(((((((r*1.05075007164441684324e-9+5.475938084995344946e-4)*r+0.0151986665636164571966)*r+0.14810397642748007459)*r+0.68976733498510000455)*r+1.6763848301838038494)*r+2.05319162663775882187)*r+1.0);
        } else {
            r-=5.0;
            quantileResult=(((((((r*2.01033439929228813265e-7+2.71155556874348757815e-5)*r+0.0012426609473880784386)*r+0.026532189526576123093)*r+0.29656057182850489123)*r+1.7848265399172913358)*r+5.4637849111641143699)*r+6.6579046435011037772)
                 /(((((((r*2.04426310338993978564e-15+1.4215117583164458887e-7)*r+1.8463183175100546818e-5)*r+7.868691311456132591e-4)*r+0.0148753612908506148525)*r+0.13692988092273580531)*r+0.59983220655588793769)*r+1.0);
        }
        return q<0 ? -quantileResult : quantileResult;
    }
}

double Stat::Normal::quantile_01(const double value){
    double clamped_value = std::max(0.0, std::min(1.0, value)); // Clamp to [0, 1]
    if (clamped_value == 0.0) return -std::numeric_limits<float>::infinity();
    if (clamped_value == 1.0) return std::numeric_limits<float>::infinity();
    return normalQuantileAS241(clamped_value);
}

double Stat::Normal::pdf_01(const double value) {
    return pdf(d, value);
//    double SQRT_2PI = std::sqrt(2 * M_PI); // Precompute sqrt(2π)
//    return std::exp(-0.5 * value * value) / SQRT_2PI;
}

float Stat::InvChiSq::sample(const float df, const float scale){
    //inverse_chi_squared_distribution invchisq(df, scale);
    //return boost::math::quantile(invchisq, ranf());   // don't know why this is not correct
    
    gamma_generator sgamma(engine, gamma_distribution(0.5f*df, 1));
    return scale/(2.0f*sgamma());
}

float Stat::Gamma::sample(const float shape, const float scale){
    gamma_generator sgamma(engine, gamma_distribution(shape, scale));
    return sgamma();
}

float Stat::Beta::sample(const float a, const float b){
    beta_distribution beta(a,b);
    return boost::math::quantile(beta,ranf());
}

unsigned Stat::Bernoulli::sample(const float p){
    if (std::isnan(p)) return ranf() < 0.5 ? 1:0;
    else return ranf() < p ? 1:0;
}

unsigned Stat::Bernoulli::sample(const VectorXf &p){
    float cum = 0;
    float rnd = ranf();
    long size = p.size();
    unsigned ret = 0;
    for (unsigned i=0; i<size; ++i) {
        if (!std::isnan(p[i])) cum += p[i];
        if (rnd < cum) {
            ret = i;
            break;
        }
    }
    return ret;
}

unsigned Stat::Bernoulli::sample(const VectorXf &p, const float rnd){
    // sampler with a given random number
    float cum = 0;
    //float rnd = ranf();
    long size = p.size();
    unsigned ret = 0;
    for (unsigned i=0; i<size; ++i) {
        if (!std::isnan(p[i])) cum += p[i];
        if (rnd < cum) {
            ret = i;
            break;
        }
    }
    return ret;
}

float Stat::NormalZeroMixture::sample(const float mean, const float variance, const float p){
    return bernoulli.sample(p) ? normal.sample(mean, variance) : 0;
}

// Sample Dirichlet

VectorXf Stat::Dirichlet::sample(const int n, const VectorXf &irx){
    VectorXf ps(n);
    double sx = 0.0;
    for (int i = 0; i < n; i++)
    {
        ps[i] = gamma.sample(irx(i), 1.0);
        sx += ps[i];
    }
    ps = ps / sx;
    return ps;
}

float Stat::TruncatedNormal::sample_tail_01_rejection(const float a){  // a < x < inf
    float b = 0.5 * a*a;
    float w, v;
    do {
        float u = ranf();
        w = b - log(u);
        v = ranf();
        // Marsaglia (1964): accept x = sqrt(2w) with probability a/x, i.e. v*v <= b/w.
        // The previous test, v > w/b, never rejected because w >= b.
    } while (v * v > b / w);
    return sqrt(2.0 * w);
}

float Stat::TruncatedNormal::sample_lower_truncated(const float mean, const float sd, const float a){  // a < x < inf
//    int seed = rand();
//    float x = truncated_normal_a_sample (mean, sd, a, seed);
//    return x;
    
    if (a - mean > 5*sd) {
        return mean + sd * sample_tail_01_rejection((a - mean)/sd);
    } else {
        float alpha = (a - mean)/sd;
        // double: in float, a CDF near 1 can round below Phi(alpha), which allowed draws below a.
        double alpha_cdf = cdf_01(alpha);
        double u;
        // Exclude u = 0: if alpha_cdf underflows to 0, x = 0 would give quantile_01(0) = -inf.
        do {
            u = ranf();
        } while (!u);
        double x = alpha_cdf + (1.0 - alpha_cdf) * u;  // x ~ Uniform(alpha_cdf, 1)
        if (x <= 0 || x >= 1) std::cout << "alpha " << alpha << " alpha_cdf " << alpha_cdf << " u " << u << " x " << x << std::endl;
        return mean + sd * quantile_01(x);
    }
}

float Stat::TruncatedNormal::sample_upper_truncated(const float mean, const float sd, const float b){  // -inf < x < b
//    int seed = rand();
//    float x = truncated_normal_b_sample (mean, sd, b, seed);
//    return x;

    if (mean - b > 5*sd) {
        return mean - sd * sample_tail_01_rejection((mean - b)/sd);
    } else{
        float beta = (b - mean)/sd;
        // Keep the CDF in double so rounding near 1 does not allow draws above b.
        double beta_cdf = cdf_01(beta);
        double u;
        do {
            u = ranf();
        } while (!u);
        double x = beta_cdf * u;  // x ~ Uniform(0, beta_cdf);
        if (x <= 0 || x >= 1) std::cout << "beta " << beta << " beta_cdf " << beta_cdf << " u " << u << " x " << x << std::endl;
        return mean + sd * quantile_01(x);
    }
}
