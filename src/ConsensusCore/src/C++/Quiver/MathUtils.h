//
//  MathUtils.h
//  ConsensusCoreXCode
//
//  Created by Nigel Delaney on 2/23/15.
//  Copyright (c) 2015 Pacific Biosciences. All rights reserved.
//

#ifndef ConsensusCoreXCode_MathUtils_h
#define ConsensusCoreXCode_MathUtils_h

#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <math.h>
#include <cmath>

#define NEG_INF - INFINITY
using namespace std;

const double log_one_third = -1.098612289;

// computes log(a + b) given log(a) and log(b)
inline double logadd(double lna, double lnb)
{
    auto max_val = std::max(lna, lnb);
    if (max_val == NEG_INF) {
        return  max_val;
    }
    lna -= max_val;
    lnb -= max_val;
    auto sum = std::exp(lna) + std::exp(lnb);
    sum = std::log(sum) + max_val;
    return sum;
}
// Computes the log of the summed exp values for these.
inline double logsumlog(double v1, double v2, double v3, double v4)
{
    auto max_val = std::max(std::max(std::max(v1, v2), v3),v4);
    if (max_val == NEG_INF)
    {
        return max_val;
    }
    v1 -= max_val;
    v2 -= max_val;
    v3 -= max_val;
    v4 -= max_val;
    auto sum = exp(v1) + exp (v2) + exp(v3) + exp (v4);
    return log (sum) + max_val;
}
#endif