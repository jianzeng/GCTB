//
//  predict.hpp
//  gctb
//
//  Created by Jian Zeng on 20/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#ifndef predict_hpp
#define predict_hpp

#include <stdio.h>
#include "data.hpp"

class Predict {
public:
    VectorXf ghat;
    float cor, reg;
    
    void getAccuracy(const Data &data, const string &filename);
    void writeRes(const Data &data, const string &filename);
};

#endif /* predict_hpp */
