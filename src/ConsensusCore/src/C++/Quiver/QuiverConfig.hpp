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
#define MISMATCH_PROBABILITY 0.00426009351993522
#define PMF_BINS 21
namespace ConsensusCore
{
    // private anonymous parameters
    namespace
    {

        const double INSERT_IQV_PMF[]= {0.0072684562651915, 0.0240552650751861, 0.0237846848538813, 0.029340665076333, 0.034201839702731, 0.0386722852646867, 0.0437507891849498, 0.0480601329196798, 0.048676216487064, 0.0434572016967914, 0.0324746636513489, 0.0225080398050568, 0.0164014664673396, 0.00943317342172977, 0.00252650564181938, 0.000413787960651576, 5.89243179825966e-05, 2.59571927009966e-14, 4.03060292021505e-14, 8.16338579512612e-13, 0.574915902206694};
        
        const double MATCH_IQV_PMF[]  = {0.00129442722911149, 0.014661557601951, 0.0295917604059299, 0.0569817671501941, 0.0945487132632137, 0.130061943434755, 0.147276818474665, 0.140239360284842, 0.112053751572591, 0.080113558106073, 0.0506341590620082, 0.0290391202860332, 0.0163367616127287, 0.0077230580035009, 0.00268746849468852, 0.000719037588632366, 0.000118088918985325, 5.37310083566748e-06, 1.07462016724177e-05, 2.14924032762575e-05, 0.0858810368043126};
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
