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

//TODO: There is not a single comment in this file.

#include <iostream>
#include <string>

#include "PairwiseAlignment.hpp"
#include "Quiver/QuiverConfig.hpp"
#include "Quiver/ReadScorer.hpp"
#include "Matrix/SparseMatrix.hpp"
#include "Quiver/SimpleRecursor.hpp"

using std::string;
using std::cout;
using std::endl;

namespace ConsensusCore
{
    template<typename R>
    ReadScorer<R>::ReadScorer(const QuiverConfig& config)
        : _quiverConfig(config)
    {}

    template<typename R>
    double ReadScorer<R>::Score(const WrappedTemplateParameterPair& tpl, const Read& read) const
        throw(AlphaBetaMismatchException)
    {
        int I, J;
        R r(_quiverConfig.QvParams, read, tpl, _quiverConfig.Banding);
        
        I = read.Length();
        J = (int)tpl.Length();
        SparseMatrix alpha(I+1, J+1), beta(I+1, J+1);
        r.FillAlphaBeta(alpha, beta);

        return std::log(beta(0, 0)) + beta.GetLogProdScales();
    }
    
    template<typename R>
    const SparseMatrix*
    ReadScorer<R>::Alpha(const WrappedTemplateParameterPair& tpl, const Read& read) const
        throw(AlphaBetaMismatchException)
    {
        int I, J;
        R r(_quiverConfig.QvParams, read, tpl, _quiverConfig.Banding);

        I = read.Length();
        J = tpl.Length();
        SparseMatrix *alpha = new SparseMatrix(I+1, J+1);
        SparseMatrix *beta  = new SparseMatrix(I+1, J+1);
        r.FillAlphaBeta(*alpha, *beta);
        return alpha;
    }

    template<typename R>
    const SparseMatrix*
    ReadScorer<R>::Beta(const WrappedTemplateParameterPair& tpl, const Read& read) const
        throw(AlphaBetaMismatchException)
    {
        int I, J;
        R r(_quiverConfig.QvParams, read, tpl, _quiverConfig.Banding);
        
        I = read.Length();
        J = tpl.Length();
        SparseMatrix *alpha = new SparseMatrix(I+1, J+1);
        SparseMatrix *beta  = new SparseMatrix(I+1, J+1);
        r.FillAlphaBeta(*alpha, *beta);
        return beta;
    }

    typedef ReadScorer<SparseSimpleQvRecursor> \
      SparseSimpleQvReadScorer;
    typedef ReadScorer<SparseSimpleQvSumProductRecursor> \
      SparseSimpleQvSumProductReadScorer;
    
    template class ReadScorer<SparseSimpleQvSumProductRecursor>;
    template class ReadScorer<SparseSimpleQvRecursor>;
}
