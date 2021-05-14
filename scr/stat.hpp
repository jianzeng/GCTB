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
#include <Eigen/Eigen>

using namespace Eigen;
//using namespace std;


namespace Stat {
    
    typedef boost::mt19937 random_engine;
    typedef boost::uniform_01<> uniform_01;
    typedef boost::normal_distribution<> normal_distribution;
    typedef boost::gamma_distribution<> gamma_distribution;
    typedef boost::math::inverse_chi_squared_distribution<> inverse_chi_squared_distribution;
    typedef boost::math::beta_distribution<> beta_distribution;
    
    typedef boost::variate_generator<random_engine&, uniform_01> uniform01_generator;
    typedef boost::variate_generator<random_engine&, normal_distribution> normal_generator;
    typedef boost::variate_generator<random_engine&, gamma_distribution> gamma_generator;
    
    static random_engine engine;
    static uniform01_generator ranf(engine, uniform_01());
    static normal_generator snorm(engine, normal_distribution(0,1));  // standard normal
    
    void seedEngine(const int seed);
    
    class Normal {
    public:
        boost::math::normal_distribution <> d;
        float SQRT_2;
        float Inv_SQRT_2;
        
        Normal(){
            d = boost::math::normal_distribution <> (0 ,1);
            SQRT_2 = sqrt(2);
            Inv_SQRT_2 = 1.0/SQRT_2;
        }
        
        float sample(const float mean, const float variance);
        float cdf_01(const double value);
        float quantile_01(const double value);
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
