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


#include <xmmintrin.h>
#include <pmmintrin.h>

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <climits>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>


#include "Quiver/QuiverConfig.hpp"
#include "Quiver/MathUtils.h"
#include "Features.hpp"
#include "Types.hpp"
#include "Utils.hpp"
#include "Read.hpp"

#ifndef SWIG
using std::min;
using std::max;
#endif  // SWIG



namespace ConsensusCore
{
    //
    // Utility functions
    //
    static inline char encodeTplBase(char base)
    {
        switch (base) {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            case 'M': return 4;  // For testing
            case 'N': return 5;  // For testing
            default:  ShouldNotReachHere();
        }
    }

    //
    // Evaluator classes
    //

    /// \brief An Evaluator that can compute move scores using a QvSequenceFeatures
    class QvEvaluator
    {
    public:
        typedef ModelParams      ParamsType;

    public:
        QvEvaluator(const Read& read,
                    const std::string& tpl,
                    const ModelParams& params)
            : read_(read),
              params_(params),
              tpl_(tpl),
              pinStart_(true),
              pinEnd_(true)
        {}

        ~QvEvaluator()
        {}

        std::string ReadName() const
        {
            return read_.Name;
        }

        std::string Template() const
        {
            return tpl_;
        }

        void Template(std::string tpl)
        {
            tpl_ = tpl;
        }


        int ReadLength() const
        {
            return read_.Length();
        }

        int TemplateLength() const
        {
            return tpl_.length();
        }

        bool PinEnd() const
        {
            return pinEnd_;
        }

        bool PinStart() const
        {
            return pinStart_;
        }

        bool IsMatch(int i, int j) const
        {
            assert(0 <= i && i < ReadLength());
            assert (0 <= j && j < TemplateLength());
            return (read_.Sequence[i] == tpl_[j]);
        }
        /* This calculates the incorporation (Match) score by comparing read base i
           and template base j, this corresponds to positions [(i+1),(j+1)] in
           the alpha matrix.
        */
        double Match(int i, int j) const
        {
            assert(0 <= j && j < TemplateLength() &&
                   0 <= i && i < ReadLength() );
            auto emission_prob = (IsMatch(i, j)) ?
                params_.log_one_minus_miscallprobability :
                params_.log_miscallprobability + log_one_third;
            return read_.trans_probs[i-1].Match + emission_prob;
        }
        /* This calculates the deletion score for a move out of template position j
           corresponding to alpha column j+1 -> j+2
         */
        double Deletion(int j) const
        {
            // Can't delete the last template base
            assert(0 <= j && j < TemplateLength());
            return read_.trans_probs[j].Deletion;
        }

        /* This calculates the insertion score by comparing read base i
         and template base j, this corresponds to positions [(i+1),(j+1)] in
         the alpha matrix.
         
         Note that for this comparison, the transition parameters are determined by
         template base J, but whether it is a branch or stick is determined by template
         base J+1
         Typically, this will be called as [(i-1), j]
         */
        double Insertion(int i, int j) const
        {
            assert(1 <= j && j <= (TemplateLength() -1) &&
                   0 <= i && i < (ReadLength() -1) );
            
            return (j < TemplateLength() && IsMatch(i, j)) ?
            read_.trans_probs[j-1].Branch : read_.trans_probs[j-1].Stick + log_one_third;
        }
    



    protected:
        Read read_;
        ModelParams params_;
        std::string tpl_;
        bool pinStart_;
        bool pinEnd_;
    };
}
