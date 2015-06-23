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

        const double INSERT_IQV_PMF[] =  {0.501821163640924, 0.00283884596666172, 0.00933423183225264, 0.00988386173870922, 0.0136500321346603, 0.0164142513151896, 0.0201900339403832, 0.0225358175440831, 0.0251487240942705, 0.0235731709030832, 0.0205337679747805, 0.0147962549014448, 0.0111277088854302, 0.00746319061677649, 0.00395329827432197, 0.000872383753013508, 0.000311370423255691, 1.31024952373585e-05, 9.64599207305972e-08, 0, 0.295538693105601};
        
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
