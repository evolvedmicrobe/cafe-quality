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
#define MISMATH_PROBABILITY 0.002671256

namespace ConsensusCore
{

    /// \brief The banding optimizations to be used by a recursor
    struct BandingOptions
    {
        double ScoreDiff;

        BandingOptions(int diagonalCross, double scoreDiff)
            : ScoreDiff(scoreDiff)
        {}

        BandingOptions(int diagonalCross, double scoreDiff,
                       double dynamicAdjustFactor, double dynamicAdjustOffset)
            : ScoreDiff(scoreDiff)
        {}
    };


    /// \brief A parameter vector for analysis using the QV model
    struct ModelParams
    {
        double log_miscallprobability;
        double log_one_minus_miscallprobability;
        double log_miscall_probability_times_one_third;
        //
        // Constructor for single merge rate and merge rate slope
        //
        ModelParams()
        : log_miscallprobability(log(MISMATH_PROBABILITY))
        {
            log_one_minus_miscallprobability = log(1 - MISMATH_PROBABILITY);
            log_miscall_probability_times_one_third = log_one_third + log_miscallprobability;
        }
    };


    class QuiverConfig
    {
        public:
            ModelParams QvParams;
            ContextParameters Ctx_params;
            BandingOptions Banding;
            double FastScoreThreshold;
            double AddThreshold;

            QuiverConfig(const ContextParameters& dinucleotide_params,
                         const BandingOptions& bandingOptions,
                         double fastScoreThreshold = -12.5,
                         double addThreshold = 1.0f);

            QuiverConfig(const QuiverConfig& qvConfig);
        
            // Assuming compiler generated destructor is sufficient.
    };



    class QuiverConfigTable
    {
    private:
        typedef std::pair<const std::string, const QuiverConfig> QuiverConfigTableEntry;
        std::list<QuiverConfigTableEntry> table;

    public:
        typedef std::list<QuiverConfigTableEntry>::const_iterator const_iterator;

        QuiverConfigTable();

        bool Insert(const std::string& name, const QuiverConfig& config);
        int Size() const;

        const QuiverConfig& At(const std::string& name) const throw(InvalidInputError);

        std::vector<std::string> Keys() const;

#ifndef SWIG
        const_iterator begin() const;
        const_iterator end() const;
#endif
    };
}
