//
//  stat.hpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#ifndef stat_hpp
#define stat_hpp

#include <stdio.h>
//#include <random>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <Eigen/Eigen>

using namespace Eigen;
//using namespace std;


namespace Stat {

    /// Build with GCTB_COUNT_RNG_DRAWS to count engine draws for reproducibility debugging.
#ifdef GCTB_COUNT_RNG_DRAWS
    class CountedMt19937 {
        boost::mt19937 impl_;
    public:
        typedef boost::mt19937::result_type result_type;
        CountedMt19937();
        explicit CountedMt19937(result_type s);
        static result_type min BOOST_PREVENT_MACRO_SUBSTITUTION () { return (boost::mt19937::min)(); }
        static result_type max BOOST_PREVENT_MACRO_SUBSTITUTION () { return (boost::mt19937::max)(); }
        void seed(result_type s);
        result_type operator()();
    };

    typedef CountedMt19937 random_engine;
#else
    typedef boost::mt19937 random_engine;
#endif
    typedef boost::uniform_01<> uniform_01;
    typedef boost::normal_distribution<> normal_distribution;
    typedef boost::gamma_distribution<> gamma_distribution;
    typedef boost::math::inverse_chi_squared_distribution<> inverse_chi_squared_distribution;
    typedef boost::math::beta_distribution<> beta_distribution;
    
    typedef boost::variate_generator<random_engine&, uniform_01> uniform01_generator;
    typedef boost::variate_generator<random_engine&, normal_distribution> normal_generator;
    typedef boost::variate_generator<random_engine&, gamma_distribution> gamma_generator;
    
    extern thread_local random_engine engine;
    extern thread_local uniform01_generator ranf;
    extern thread_local normal_generator snorm;  // standard normal
    
    void seedEngine(const int seed);
    /// Set the current MCMC context used to derive deterministic parallel streams.
    void setMcmcIteration(unsigned iter);
    void setMcmcContext(unsigned chain, unsigned iter);
    unsigned currentMcmcChain();
    /// Register chain ids for nested OpenMP workers before/after the outer chain loop.
    void prepareParallelChainRng(unsigned numOuterThreads);
    void finishParallelChainRng();
    /// Mix master seed with chain, MCMC iteration, task id, and purpose salt.
    random_engine::result_type parallelTaskSeed(unsigned taskKey, unsigned salt = 0u);
    void seedEngineForParallelTask(unsigned taskKey, unsigned salt = 0u);
    /// No-op when OpenMP uses a single thread (preserves serial RNG stream).
    void seedEngineForParallelTaskIfNeeded(unsigned taskKey, unsigned salt = 0u);
    /// Purpose salts mixed into deterministic streams.
    static const unsigned pseudoSummaryRngSalt = 0x50534441u;
    static const unsigned blockShuffleRngSalt = 0x9E3779B9u;
    static const unsigned tgsBlockRngSalt = 0x54475342u;
    static const unsigned ldSwapRngSalt = 0x4C44534Du;
    static const unsigned annotationSwapRngSalt = 0x414E5357u;

    extern unsigned masterSeed;
    extern unsigned mcmcIterationForRng;
    extern unsigned mcmcChainForRng;

    /// Give one model/chain an independent scalar stream for this iteration.
    class McmcChainRngScope {
        random_engine savedEngine;
        normal_distribution savedNormalDistribution;
        bool savedChainOverrideActive;
        unsigned savedChainOverride;
    public:
        explicit McmcChainRngScope(unsigned chain);
        ~McmcChainRngScope();
    private:
        McmcChainRngScope(const McmcChainRngScope&);
        McmcChainRngScope& operator=(const McmcChainRngScope&);
    };

    /// Isolate random draws made by one parallel task from its worker's scalar stream.
    class ParallelTaskRngScope {
        random_engine savedEngine;
        normal_distribution savedNormalDistribution;
    public:
        ParallelTaskRngScope(unsigned taskKey, unsigned salt = 0u);
        ~ParallelTaskRngScope();
    private:
        ParallelTaskRngScope(const ParallelTaskRngScope&);
        ParallelTaskRngScope& operator=(const ParallelTaskRngScope&);
    };

    extern thread_local unsigned long long rngDrawCount;
    void logRngCheckpoint(const char *label);
    
    class Normal {
    public:
        boost::math::normal_distribution <> d;
        double SQRT_2;
        double Inv_SQRT_2;
        
        Normal(){
            d = boost::math::normal_distribution <> (0 ,1);
            SQRT_2 = sqrt(2);
            Inv_SQRT_2 = 1.0/SQRT_2;
        }
        
        double sample(const double mean, const double variance);
        double cdf_01(const double value);
        double quantile_01(const double value);
        double pdf_01(const double value);
    };
    
    class Flat : public Normal {
    public:
        // flat prior is a normal with infinit variance
    };
    
    class InvChiSq {
    public:
        float sample(const float df, const float scale);
    };
    
    class Gamma {
    public:
        float sample(const float shape, const float scale);
    };
    
    class Beta {
    public:
        float sample(const float a, const float b);
    };
    
    class Bernoulli {
    public:
        unsigned sample(const float p);
        unsigned sample(const VectorXf &p); // multivariate sampling, return the index of component.
        unsigned sample(const VectorXf &p, const float rnd); // sampling with a given random number.
    };

    class Dirichlet {
    public:
        Gamma gamma;
        VectorXf sample(const int n, const VectorXf &irx);
    };
    
    class NormalZeroMixture {
    public:
        Normal normal;
        Bernoulli bernoulli;
        
        float sample(const float mean, const float variance, const float p);
    };
    
    class MixtureNormals {
    public:
        
    };

    class TruncatedNormal : public Normal {
    public:
        
        float sample_tail_01_rejection(const float a);
        float sample_lower_truncated(const float mean, const float sd, const float a);  // a < x < inf
        float sample_upper_truncated(const float mean, const float sd, const float b);  // -inf < x < b
    };
}

#endif /* stat_hpp */
