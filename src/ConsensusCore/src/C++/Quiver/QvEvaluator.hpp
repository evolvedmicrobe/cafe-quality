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



#pragma once



#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "Quiver/QuiverConfig.hpp"
#include "Quiver/MathUtils.h"
#include "Features.hpp"
#include "Types.hpp"
#include "Utils.hpp"
#include "Read.hpp"
#include "TransitionParameters.hpp"
#include "TemplateParameterPair.hpp"

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
        typedef ModelParams ParamsType;

    public:
        QvEvaluator(const Read& read,
                    const TemplateParameterPair& tpl,
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

        TemplateParameterPair Template() const
        {
            return tpl_;
        }

        void Template(TemplateParameterPair tpl)
        {
            tpl_ = tpl;
        }


        int ReadLength() const
        {
            return read_.Length();
        }

        int TemplateLength() const
        {
            return (int)tpl_.tpl.length();
        }

        bool PinEnd() const
        {
            return pinEnd_;
        }

        bool PinStart() const
        {
            return pinStart_;
        }
        
        /**
         Returns true if the read bp at position i matches the template basepair
         at position j.
         @param i read location
         @param j template location
         @returns tru if matching
         */
        bool IsMatch(int i, int j) const
        {
            assert(0 <= i && i < ReadLength());
            assert (0 <= j && j < TemplateLength());
            return (read_.Sequence[i] == tpl_.tpl[j]);
        }
        
        /**
         This calculates the incorporation (Match) score by comparing read base i
         and template base j, this corresponds to positions [(i+1),(j+1)] in
         the alpha matrix.  The transition probability to a match is determined by the
         j -1 base.
         @param i the read position (0 - indexed)
         @param j the template position (j-indexed), the likelihood of transitioning to a match is stored in the preceding template position (j-1)
         @returns the probability of transitioning to a match state and emitting that base.
         */
        /*
        */
        double Match(int i, int j) const
        {
            assert(0 <= j && j < TemplateLength() &&
                   0 <= i && i < ReadLength() );
            auto emission_prob = (IsMatch(i, j)) ?
                params_.log_one_minus_miscallprobability :
                params_.log_miscall_probability_times_one_third;
            return tpl_.trans_probs[j-1].Match + emission_prob;
        }
        /**
         This is a special edge case match score, for the first position, which
         does not have a dinucleotide context to
         @param <#parameter#>
         @returns <#retval#>
         @exception <#throws#>
         */
        double Match_Just_Emission(int i, int j) const
        {
            assert( (j==0 && i==0) || (j == (TemplateLength()-1) && i == (ReadLength()-1)) ); // The arguments are just for readability.
            auto emission_prob = (IsMatch(i, j)) ?
            params_.log_one_minus_miscallprobability :
            params_.log_miscall_probability_times_one_third;
            return emission_prob;
        }
        
        
        /**
            This calculates the deletion score for a move out of template position j
            corresponding to alpha column j+1 -> j+2.
            The emission probability for a deletion is 1, and the transition probability 
            is determined by the parameters at the jth base.
         @param j template position to get probability from
         @returns Log scale transition probability.
         */
        double Deletion(int j) const
        {
            // Can't delete the last template base
            assert(0 <= j && j < TemplateLength());
            return tpl_.trans_probs[j].Deletion;
        }

        /** This calculates the insertion score by comparing read base i
         and template base j, this corresponds to positions [(i+1),(j+1)] in
         the alpha matrix.
         
         Note that for this comparison, the transition parameters are determined by
         template base J, but whether it is a branch or stick is determined by template
         base J+1
         @param i The index of the base that will be inserted
         @param j The index of the template base that determines the insertion probability
         @returns The probability of transitioning to insertion and emitting read.
         */
        double Insertion(int i, int j) const
        {
            assert(0 <= j && j < (TemplateLength()-1) &&
                   0 <= i && i < (ReadLength() -1) );
            
            return IsMatch(i, j + 1) ? tpl_.trans_probs[j].Branch : tpl_.trans_probs[j].Stick + log_one_third;
        }
    



    protected:
        Read read_;
        ModelParams params_;
        TemplateParameterPair tpl_;
        bool pinStart_;
        bool pinEnd_;
    };
}
