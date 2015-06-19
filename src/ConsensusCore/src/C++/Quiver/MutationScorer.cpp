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

#include <fstream>
#include <iostream>
#include <string>

#include "Matrix/DenseMatrix.hpp"
#include "Matrix/SparseMatrix.hpp"
#include "Quiver/SimpleRecursor.hpp"
#include "Mutation.hpp"
#include "TemplateParameterPair.hpp"

#define EXTEND_BUFFER_COLUMNS 8

namespace ConsensusCore
{
    template<typename R>
    MutationScorer<R>::MutationScorer(const R& recursor)
        throw(AlphaBetaMismatchException)
        : recursor_(new R(recursor))
    {
        try {
            cachedMatchScalingCorrectionFactor = -std::log(recursor_->MatchScalingFactor) * recursor_->read_.Length();
            int I = (int)recursor_->read_.Sequence.length() + 1;
            int J = recursor_->tpl_.Length() + 1;
            // Allocate alpha and beta
            alpha_ = new MatrixType(I, J);
            beta_ = new MatrixType(I, J);
            // Buffer where we extend into
            extendBuffer_ = new MatrixType(I, EXTEND_BUFFER_COLUMNS);
            // Initial alpha and beta
            numFlipFlops_ = recursor.FillAlphaBeta(*alpha_, *beta_);
        }
        catch(AlphaBetaMismatchException e) {
            delete alpha_;
            delete beta_;
            delete extendBuffer_;
            delete recursor_;
            throw;
        }
    }

    template<typename R>
    MutationScorer<R>::MutationScorer(const MutationScorer<R>& other)
    {
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
        
        return std::log((*beta_)(0, 0)) + beta_->GetLogProdScales() + MatchScalingFactorCorrection();
    }

    template<typename R>
    double
    MutationScorer<R>::MatchScalingFactorCorrection() const
    {
        return cachedMatchScalingCorrectionFactor;
    }
   
    template<typename R>
    void MutationScorer<R>::DumpAlphaMatrix() const
    {
        DumpMatrix(*alpha_, "/Users/nigel/git/cafe-quality/src/Tests/bin/Release/Alpha.csv");
    }
    
    template<typename R>
    void MutationScorer<R>::DumpBetaMatrix() const
    {
        DumpMatrix(*beta_, "Beta.csv");
    }

    template<typename R> WrappedTemplateParameterPair
    MutationScorer<R>::Template() const
    {
        return recursor_->tpl_;
    }

    template<typename R>
    void MutationScorer<R>::Template(WrappedTemplateParameterPair tpl)
        throw(AlphaBetaMismatchException)
    {
        delete alpha_;
        delete beta_;
        recursor_->tpl_ = tpl;
        int I = (int)recursor_->read_.Sequence.length() + 1;
        int J = recursor_->tpl_.Length() + 1;
        alpha_ = new MatrixType(I, J);
        beta_  = new MatrixType(I,J);
        recursor_->FillAlphaBeta(*alpha_, *beta_);
    }
    
    template<typename R>
    void MutationScorer<R>::DumpMatrix(const MatrixType& mat, const std::string& fname) const
    {
        if (mat.Rows() == 0 || mat.Columns() == 0) return;
        std::ofstream myfile;
        myfile.open(fname);
        for(int i=0; i< mat.Rows(); i++)
        {
            myfile << mat(i,0);
            for(int j=1; j<mat.Columns(); j++)
            {
                myfile << "," << mat(i,j);
            }
            myfile << std::endl;
        }
        myfile << mat.GetLoggedScale(0);
        for (int j=1; j < mat.Columns(); j++)
        {
            myfile << "," << mat.GetLoggedScale(j);
        }
        myfile << std::endl;
        myfile.close();
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
    double
    MutationScorer<R>::ScoreMutation(const Mutation& m) const
    {
        /*  This is a very weak guard on someone trying to score a
            mutation without first applying a virtual mutation to an underlying
            template.  In the future, I should just ensure this method can only be
            called from a MultiReadMutationScorer
         */
        if (!(recursor_->tpl_.VirtualMutationActive())) {
            throw BadExecutionOrderException();
        }

        if(std::abs(m.LengthDiff()) > 1) {
            throw new InvalidInputError("Only mutations of size 1 allowed");
        }

        int betaLinkCol = 1 + m.End();
        int absoluteLinkColumn = 1 + m.End() + m.LengthDiff();
        
        double score;

        bool atBegin = (m.Start() < 3);
        bool atEnd   = (m.End() > (int)recursor_->tpl_.Length() - 2);

        if (!atBegin && !atEnd)
        {
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
                extendLength   = 1 + (int)m.NewBases().length();
                assert(extendLength <= EXTEND_BUFFER_COLUMNS);
            }
            
            recursor_->ExtendAlpha(*alpha_,
                                   extendStartCol, *extendBuffer_, extendLength);
            score = recursor_->LinkAlphaBeta(*extendBuffer_, extendLength,
                                             *beta_, betaLinkCol,
                                             absoluteLinkColumn);
            score += alpha_->GetLogProdScales(0, extendStartCol) + MatchScalingFactorCorrection();
        }
        else if (!atBegin && atEnd)
        {
            //
            // Extend alpha to end
            //
            int extendStartCol = m.Start() - 1;
            int extendLength = (int)recursor_->tpl_.VirtualLength() - extendStartCol + 1;

            recursor_->ExtendAlpha(*alpha_,
                                   extendStartCol, *extendBuffer_, extendLength);
            score = (std::log((*extendBuffer_)((int)recursor_->read_.Sequence.length(), extendLength - 1))
                    + alpha_->GetLogProdScales(0, extendStartCol)
                    + extendBuffer_->GetLogProdScales(0, extendLength)
                    + MatchScalingFactorCorrection());
        }
        else if (atBegin && !atEnd)
        {
            // If the mutation occurs at positions 0 - 2
            int extendLastCol = m.End();
            // We duplicate this math inside the function
            int extendLength = m.End() + m.LengthDiff() + 1;
            
            recursor_->ExtendBeta(*beta_, extendLastCol,
                                  *extendBuffer_, m.LengthDiff());
            score = ( std::log((*extendBuffer_)(0, 0))
                    + beta_->GetLogProdScales(extendLastCol + 1, beta_->Columns())
                    + extendBuffer_->GetLogProdScales(0, extendLength)
                     + MatchScalingFactorCorrection());
        }
        else
        {
            assert(atBegin && atEnd);
            /* This should basically never happen...
               and is a total disaster if it does.  The basic idea is that
               FillAlpha and FillBeta use the real "template" while we test 
               mutations using "virtual" template positions and the Extend/Link 
               methods.  Trying to call FillAlpha to calculate the likelihood of a virtual
               mutation is therefore going to fail, as it calculates using the
               "real" template.
             */
            throw TooSmallTemplateException();
            
            //
            // Just do the whole fill
            //
//            MatrixType alphaP(recursor_->read_.Length() + 1,
//                              new_tpl->tpl.length() + 1);
//            evaluator_->Template(*new_tpl);
//            recursor_->FillAlpha(*evaluator_, MatrixType::Null(), alphaP);
//            score = ( std::log(alphaP(evaluator_->ReadLength(), new_tpl->tpl.length()))
//                    + alphaP.GetLogProdScales() );
        }

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
    }

    template class MutationScorer<SimpleQvRecursor>;
    template class MutationScorer<SimpleQvSumProductRecursor>;
    template class MutationScorer<SparseSimpleQvRecursor>;
    template class MutationScorer<SparseSimpleQvSumProductRecursor>;

}
