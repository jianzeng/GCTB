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
    return cdf(d, value);
}

double Stat::Normal::quantile_01(const double value){
    double clamped_value = std::max(0.0, std::min(1.0, value)); // Clamp to [0, 1]
    if (clamped_value == 0.0) return -std::numeric_limits<float>::infinity();
    if (clamped_value == 1.0) return std::numeric_limits<float>::infinity();
    return quantile(d, clamped_value);
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
    } while (v > w/b);
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
        float alpha_cdf = cdf_01(alpha);
        double u = ranf();
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
        float beta_cdf = cdf_01(beta);
        double u;
        do {
            u = ranf();
        } while (!u);
        double x = beta_cdf * u;  // x ~ Uniform(0, beta_cdf);
        if (x <= 0 || x >= 1) std::cout << "beta " << beta << " beta_cdf " << beta_cdf << " u " << u << " x " << x << std::endl;
        return mean + sd * quantile_01(x);
    }
}
