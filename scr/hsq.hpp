//
//  hsq.hpp
//  gctb
//
//  Created by Jian Zeng on 20/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#ifndef hsq_hpp
#define hsq_hpp

#include <stdio.h>
#include "data.hpp"
#include "mcmc.hpp"


class Heritability {
public:
    McmcSamples varGenotypic;
    McmcSamples varResidual;
    McmcSamples hsq;
    
    unsigned popSize;
    
    void getEstimate(const Data &data, const McmcSamples &snpEffects, const float phenVar);
    void writeRes(const string &filename);
    void writeMcmcSamples(const string &filename);
    
    Heritability(): varGenotypic("GenVar"), varResidual("ResVar"), hsq("hsq"){}
};

#endif /* hsq_hpp */
