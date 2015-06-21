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

#include <stdexcept>
#include <string>
#include <vector>

#include "Quiver/QuiverConfig.hpp"

namespace ConsensusCore {
    QuiverConfig::QuiverConfig(const ContextParameters& ctxParams,
                               const BandingOptions& bandingOptions,
                               double fastScoreThreshold,
                               double addThreshold)
        : QvParams()
        , CtxParams(ctxParams)
        , Banding(bandingOptions)
        , FastScoreThreshold(fastScoreThreshold)
        , AddThreshold(addThreshold)
    {
        
        /* Calculate a match scaling factor, going for somthing near the inverse
           of the expected score for an emission. (But note this is not the expected
           value, just related to it 
         */
        double prob =0;
        double prior = 1 / (double) ctxParams.contexts.size();
        for (auto &ctx : ctxParams.contexts) {
            auto probs = ctxParams.GetParametersForContext(ctx[0], ctx[1]);
            //TODO: Revisit this, note that I should be integrating over the bin probabilities, but because these
            // may be set at runtime, this will require some refactoring.  In the future this could be changed.
            // just seeing how it works for now...
            prob += prior * (1.0 / PMF_BINS) * (probs.Match + probs.Branch + probs.Stick);
        }
        MatchScalingFactor = 1 / prob;

    }

    QuiverConfig::QuiverConfig(const QuiverConfig& qvConfig)
        : QvParams(qvConfig.QvParams),
          CtxParams(qvConfig.CtxParams),
          Banding(qvConfig.Banding),
          FastScoreThreshold(qvConfig.FastScoreThreshold),
          AddThreshold(qvConfig.AddThreshold),
          MatchScalingFactor(qvConfig.MatchScalingFactor)
    {}


    
}
