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
#define MISMATCH_PROBABILITY 0.00416541391862461


namespace ConsensusCore
{
    // private anonymous parameters
    namespace
    {
        const double INSERT_IQV_PMF[]= {0, 0.398450512344855, 0.1159414219823, 0.131649010580163, 0.0858029945309225, 0.0537532970859549, 0.034406135234615, 0.0249021814886181, 0.0185778274137733, 0.0172867594420362, 0.0113596854045226, 0.0113338026349943, 0.00925507600451598, 0.00977233056331712, 0.00766698518308035, 0.00515244493156259, 0.00480540943157229, 0.00412033607840733, 0.00455217476472953, 0.00375502498302623, 0.0474565899170333};
        
        const double MATCH_IQV_PMF[]  = {0, 0.040102840489951, 0.0376322290309637, 0.0974954582368313, 0.110007947413278, 0.0979197378060836, 0.0829516058979775, 0.0703783192772605, 0.0597271797647386, 0.050648847109636, 0.0429458513349509, 0.0356268974922942, 0.0306088572131205, 0.0262933614919783, 0.022635320175671, 0.0198070995987446, 0.0172610952234939, 0.0150629851495138, 0.0132749873057055, 0.0117210208439848, 0.117898359143822};
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
        double MatchIqvPmf[21];
        double InsertIqvPmf[21];
        double PrMiscall;
        double PrNotMiscall;
        double PrThirdOfMiscall;
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

            QuiverConfig(const ContextParameters& ctxParams,
                         const BandingOptions& bandingOptions,
                         double fastScoreThreshold = -12.5,
                         double addThreshold = 1.0f);

            QuiverConfig(const QuiverConfig& qvConfig);
        
            // Assuming compiler generated destructor is sufficient.
    };

}
