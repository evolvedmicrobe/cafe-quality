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


#include "Quiver/MutationScorer.hpp"

#include <string>

#include "Matrix/DenseMatrix.hpp"
#include "Matrix/SparseMatrix.hpp"
#include "Quiver/QvEvaluator.hpp"
#include "Quiver/SimpleRecursor.hpp"
#include "Mutation.hpp"
#include "TemplateParameterPair.hpp"

#define EXTEND_BUFFER_COLUMNS 8

namespace ConsensusCore
{
    template<typename R>
    MutationScorer<R>::MutationScorer(const EvaluatorType& evaluator, const R& recursor)
        throw(AlphaBetaMismatchException)
        : evaluator_(new EvaluatorType(evaluator)),
          recursor_(new R(recursor))
    {
        // Allocate alpha and beta
        alpha_ = new MatrixType(evaluator.ReadLength() + 1,
                                evaluator.TemplateLength() + 1);
        beta_ = new MatrixType(evaluator.ReadLength() + 1,
                               evaluator.TemplateLength() + 1);
        // Buffer where we extend into
        extendBuffer_ = new MatrixType(evaluator.ReadLength() + 1, EXTEND_BUFFER_COLUMNS);
        // Initial alpha and beta
        numFlipFlops_ = recursor.FillAlphaBeta(*evaluator_, *alpha_, *beta_);
    }

    template<typename R>
    MutationScorer<R>::MutationScorer(const MutationScorer<R>& other)
    {
        evaluator_ = new EvaluatorType(*other.evaluator_);
        recursor_ = new R(*other.recursor_);

        // Copy alpha and beta
        alpha_ = new MatrixType(*other.alpha_);
        beta_ = new MatrixType(*other.beta_);
        // Buffer where we extend into
        extendBuffer_ = new MatrixType(*other.extendBuffer_);
        numFlipFlops_ = other.numFlipFlops_;
    }

    template<typename R>
    double
    MutationScorer<R>::Score() const
    {
        return (*beta_)(0, 0);
    }

    template<typename R> TemplateParameterPair
    MutationScorer<R>::Template() const
    {
        return evaluator_->Template();
    }

    template<typename R>
    void MutationScorer<R>::Template(TemplateParameterPair tpl)
        throw(AlphaBetaMismatchException)
    {
        delete alpha_;
        delete beta_;
        evaluator_->Template(tpl);
        alpha_ = new MatrixType(evaluator_->ReadLength() + 1,
                                evaluator_->TemplateLength() + 1);
        beta_  = new MatrixType(evaluator_->ReadLength() + 1,
                                evaluator_->TemplateLength() + 1);
        recursor_->FillAlphaBeta(*evaluator_, *alpha_, *beta_);
    }

    template<typename R>
    const typename R::MatrixType* MutationScorer<R>::Alpha() const
    {
        return alpha_;
    }

    template<typename R>
    const typename R::MatrixType* MutationScorer<R>::Beta() const
    {
        return beta_;
    }

    template<typename R>
    const typename R::EvaluatorType* MutationScorer<R>::Evaluator() const
    {
        return evaluator_;
    }

//    template<typename R>
//    const PairwiseAlignment* MutationScorer<R>::Alignment() const
//    {
//        return recursor_->Alignment(*evaluator_, *alpha_);
//    }

    template<typename R>
    double
    MutationScorer<R>::ScoreMutation(const Mutation& m, const ContextParameters& ctx_params) const
    {
        int betaLinkCol = 1 + m.End();
        int absoluteLinkColumn = 1 + m.End() + m.LengthDiff();
        TemplateParameterPair old_Tpl = evaluator_->Template();
        TemplateParameterPair new_tpl = ApplyMutation(m, old_Tpl, ctx_params);
     
        
        //TODO: Add logic here to update the parameters for all reads....
        
        double score;

        bool atBegin = (m.Start() < 3); 
        bool atEnd   = (m.End() > (int)old_Tpl.tpl.length() - 2);

        if (!atBegin && !atEnd)
        {
            // Install mutated template
            evaluator_->Template(new_tpl);

            int extendStartCol, extendLength;

            if (m.Type() == DELETION)
            {
                // Future thought: If we revise the semantic of Extra,
                // we can remove the extend and just link alpha and
                // beta directly.
                extendStartCol = m.Start() - 1;
                extendLength = 2;
            }
            else
            {
                extendStartCol = m.Start();
                extendLength   = 1 + m.NewBases().length();
                assert(extendLength <= EXTEND_BUFFER_COLUMNS);
            }

            recursor_->ExtendAlpha(*evaluator_, *alpha_,
                                   extendStartCol, *extendBuffer_, extendLength);
            score = recursor_->LinkAlphaBeta(*evaluator_,
                                             *extendBuffer_, extendLength,
                                             *beta_, betaLinkCol,
                                             absoluteLinkColumn);
        }
        else if (!atBegin && atEnd)
        {
            //
            // Extend alpha to end
            //
            evaluator_->Template(new_tpl);

            int extendStartCol = m.Start() - 1;
            int extendLength = new_tpl.tpl.length() - extendStartCol + 1;

            recursor_->ExtendAlpha(*evaluator_, *alpha_,
                                   extendStartCol, *extendBuffer_, extendLength);
            score = (*extendBuffer_)(evaluator_->ReadLength(), extendLength - 1);


        }
        else if (atBegin && !atEnd)
        {
            //
            // Extend beta back
            //
            evaluator_->Template(new_tpl);

            int extendLastCol = m.End();
            int extendLength = m.End() + m.LengthDiff() + 1;

            recursor_->ExtendBeta(*evaluator_, *beta_,
                                  extendLastCol, *extendBuffer_, extendLength,
                                  m.LengthDiff());
            score = (*extendBuffer_)(0, 0);
        }
        else
        {
            assert(atBegin && atEnd);
            //
            // Just do the whole fill
            //
            MatrixType alphaP(evaluator_->ReadLength() + 1,
                              new_tpl.tpl.length() + 1);
            evaluator_->Template(new_tpl);
            recursor_->FillAlpha(*evaluator_, MatrixType::Null(), alphaP);
            score = alphaP(evaluator_->ReadLength(), new_tpl.tpl.length());
        }

        // Restore the original template.
        evaluator_->Template(old_Tpl);

        // if (fabs(score - Score()) > 50) { Breakpoint(); }

        return score;
    }


    template<typename R>
    MutationScorer<R>::~MutationScorer()
    {
        delete extendBuffer_;
        delete beta_;
        delete alpha_;
        delete recursor_;
        delete evaluator_;
    }

    template class MutationScorer<SimpleQvRecursor>;
    template class MutationScorer<SimpleQvSumProductRecursor>;
    template class MutationScorer<SparseSimpleQvRecursor>;
    template class MutationScorer<SparseSimpleQvSumProductRecursor>;

}
