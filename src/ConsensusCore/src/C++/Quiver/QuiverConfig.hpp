// Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: David Alexander

#pragma once

#include <list>
#include <string>
#include <utility>
#include <vector>
#include <math.h>
#include "Utils.hpp"
#include "Quiver/MathUtils.h"
#include "ContextParameters.hpp"


/* Hard coded mismatch probability for now.
   This is derived as the mean in PlotBinnedTraining.R */
#define MISMATCH_PROBABILITY 0.00742185448528985
#define PMF_BINS 21
namespace ConsensusCore
{
    // private anonymous parameters
    namespace
    {

       const double INSERT_IQV_PMF[]= {0, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
       const double MATCH_IQV_PMF[]  = {0, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
        
       //const double INSERT_IQV_PMF[]= {0, 0.397798348635352, 0.115911699112264, 0.132138696992897, 0.0864300540208446, 0.0527942777746241, 0.0354284604620736, 0.0251104624063461, 0.0185084357783954, 0.017607162054922, 0.0115809987762472, 0.0111757260461004, 0.00924314791350931, 0.0097807015736225, 0.00776735793606597, 0.00500761578734957, 0.00500922413276454, 0.00411890672085751, 0.00452905901990536, 0.00380091335895079, 0.0462587514969085};
       //const double MATCH_IQV_PMF[]  = {0, 0.0402222420548246, 0.0376480138102037, 0.0974564068792295, 0.109946602873532, 0.0980000091454914, 0.0828500771501694, 0.0703516958856884, 0.0597266541262199, 0.0506139988753366, 0.0429203549490854, 0.0356372895401746, 0.0306063836981209, 0.0262898394063799, 0.0226236494380033, 0.0198178894955275, 0.0172403905468526, 0.0150612896294664, 0.0132756440949005, 0.0117154978478068, 0.117996070552987};
        
        
        const double MATCH_IQV_PMF[] =  {0.00982273914584329, 0.00148653036235552, 0.0145314464846872, 0.031052607672955, 0.0597920593735699, 0.0968491657730054, 0.13267875890945, 0.147126197659831, 0.137630883321425, 0.109607179458889, 0.0746914877493785, 0.0466715133009795, 0.0268424677580686, 0.0139536924849399, 0.00632313480332473, 0.00189745041791118, 0.000447075922997731, 7.95894084457873e-05, 2.51102250576693e-05, 0, 0.0884909097668853};
    }

    /// \brief The banding optimizations to be used by a recursor
    struct BandingOptions
    {
        double ScoreDiff;

        BandingOptions(double scoreDiff)
            : ScoreDiff(scoreDiff)
        {
            if(scoreDiff < 0) {
                throw InvalidInputError("ScoreDiff must be positive!");
            }
        }

    };


    /// \brief A parameter vector for analysis using the QV model
    struct ModelParams
    {
        double MatchIqvPmf[PMF_BINS];
        double InsertIqvPmf[PMF_BINS];
        double PrMiscall;
        double PrNotMiscall;
        double PrThirdOfMiscall;
        double MatchScalingFactor;
        //
        // Constructor for single merge rate and merge rate slope
        //
        ModelParams(const double matchIqvPmf[]  = MATCH_IQV_PMF,
                    const double insertIqvPmf[] = INSERT_IQV_PMF,
                    double mismatch             = MISMATCH_PROBABILITY)
            : PrMiscall(mismatch)
            , PrNotMiscall(1.0 - mismatch)
            , PrThirdOfMiscall(mismatch / 3.0)
        {
            memcpy(MatchIqvPmf,  matchIqvPmf,  sizeof(MatchIqvPmf));
            memcpy(InsertIqvPmf, insertIqvPmf, sizeof(InsertIqvPmf));
        }
        
        // Copy constructor
        ModelParams(const ModelParams& src) = default;
        
        // Move constructor
        ModelParams(ModelParams&& src) = default;
        
        // Copy Assignment operator
        ModelParams& operator=(const ModelParams& rhs) = default;
    };


    class QuiverConfig
    {
        public:
            ModelParams QvParams;
            ContextParameters CtxParams;
            BandingOptions Banding;
            double FastScoreThreshold;
            double AddThreshold;
            /* This scaling factor is needed to avoid the banding of the 
             recursion matrices only going along the top row.  The problem is that 
             with two emissions the probability for a match, emission and emission, 
             is nearly equivalent to a deletion.  This can lead to an "all deletion"
             path being selected, which goes to the top right of the matrix, but has 
             no hope of having substantial probability after that.  To avoid this, 
             we "reward" every emitted base being removed by multiplying it by the
             average value for an emission.  We then cancel this scaling factor out in the likelihood
             */
            double MatchScalingFactor;
            QuiverConfig(const ContextParameters& ctxParams,
                         const BandingOptions& bandingOptions,
                         double fastScoreThreshold = -12.5,
                         double addThreshold = 1.0f);

            QuiverConfig(const QuiverConfig& qvConfig);
        
            // Assuming compiler generated destructor is sufficient.
    };

}
