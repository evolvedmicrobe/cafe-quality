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


#include "Quiver/SimpleRecursor.hpp"

#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <climits>
#include <utility>

#include "Matrix/DenseMatrix.hpp"
#include "Matrix/SparseMatrix.hpp"
#include "Quiver/detail/Combiner.hpp"
#include "Quiver/detail/RecursorBase.hpp"
#include "Quiver/QvEvaluator.hpp"
#include "Interval.hpp"
#include "Utils.hpp"


using std::min;
using std::max;


namespace ConsensusCore {
    
    /**
     Fill in the alpha matrix.  This matrix has the read run along the rows, and the 
     template run along the columns.  The first row and column do not correspond to
     a template position.  Therefore the match represented at position (i,j) corresponds
     to a match between template positions (i+1, j+1).  
     
     The alpha matrix is the "Forward" matrix used in the forward/backward algorithm.
     The i,j position of the matrix represents the probability of all paths up 
     to the point where the ith read position and jth template have been "emitted."
     The matrix is calculated recursively by examining all possible transitions 
     into (i,j), and calculating the probability we were in the previous state, 
     times the probability of a transition into (i,j) times the probability of 
     emitting the observation that corresponds to (i,j). All probabilities are 
     calculated and stored as LOG values.
     
     Note that in doing this calculation, in order to work with di-nucleotide contexts, we
     require that the first and last transition be a match.  In other words the start and end of
     the read and template are "pinned" to each other.
     
     @param e An Evaluator type instance, such as QvEvaluator.  This type must be able to 
     calculate Match, Insertion and Deletion values at given parameters.  Used to populate the matrix.
     @param guide An object that helps inform how to select the size of "bands" for the 
     banded algorithm used.  This is typically the beta matrix if we are "repopulating" the matrix.
     @param alpha The matrix to be filled.
     */
    template<typename M, typename E, typename C>
    void
    SimpleRecursor<M, E, C>::FillAlpha(const E& e, const M& guide, M& alpha) const
    {
        // We are pinning, so should never go all the way to the end of the read/template
        int I = e.ReadLength();
        int J = e.TemplateLength();
        
        assert(alpha.Rows() == I + 1 && alpha.Columns() == J + 1);
        assert(guide.IsNull() ||
               (guide.Rows() == alpha.Rows() && guide.Columns() == alpha.Columns()));

        
        // Initial condition, we always start with a match
        alpha.StartEditingColumn(0, 0, 1);
        alpha.Set(0, 0, 1.0);
        alpha.FinishEditingColumn(0, 0, 1);
        // End initial conditions
        
        //
        int hintBeginRow = 1, hintEndRow = I - 1;
        
        for (int j = 1; j < J; ++j)
        {            
            this->RangeGuide(j, guide, alpha, &hintBeginRow, &hintEndRow);

            int requiredEndRow = min(I , hintEndRow);
            int i;
            double thresholdScore = 0.0;
            double maxScore = 0.0;
            double score = 0.0;
            alpha.StartEditingColumn(j, hintBeginRow, hintEndRow);

            int beginRow = hintBeginRow, endRow;
            for (i = beginRow;
                 i < I && (score >= thresholdScore || i < requiredEndRow);
                 ++i)
            {
                double thisMoveScore;
                score = 0.0;

                 if(j==0 || i ==0) {
                     assert (j !=0 && i!=0);
                    /* The first move must be a match so all others are 0.0
                       along first column/row */
                    //score = 0.0;
                }
                
                else {
                    
                    // Match:
                    if (i > 0 && j > 0) // Always true, left for clarity.
                    {
                        /* Important!  Note that because we require the initial state to be a match,
                           when i = 1 and j = 1 the match transition probability must be 1, since no other options
                           are allowed.  Similarly, the probability for the match probability to the end base should be 1.
                         
                           Note that for the first "match" between a read and template, we have no choice but to
                           hard code it to 1, as there is no defined transition probability for a dinucleotide context.
                            
                          ***********  EDGE_CONDITION ************
                         */
                        if (i == 1 && j == 1) { //TODO: Remove this branch bottleneck...
                            thisMoveScore = alpha(i-1, j-1) * e.Match_Just_Emission(0,0);
                            score = C::Combine(score, thisMoveScore);
                        }
                        else if (i != 1 && j != 1) {
                            thisMoveScore = alpha(i - 1, j - 1) * e.Match(i - 1, j - 1);
                            score = C::Combine(score, thisMoveScore);
                        }                        
                    }
                    // Stick or Branch:
                    if (i > 1) // Due to pinning, can't "insert" first or last read bpbase
                    {
                        thisMoveScore = alpha(i - 1, j) * e.Insertion(i - 1, j - 1);
                        score = C::Combine(score, thisMoveScore);
                    }

                    // Deletion:
                    if (j > 1) // Due to pinning, can't "delete" first or last template bp
                    {
                        thisMoveScore = alpha(i, j - 1) * e.Deletion(j - 2);
                        score = C::Combine(score, thisMoveScore);
                    }
                }
                //  Save score
                alpha.Set(i, j, score);
                if (score > maxScore)
                {
                    maxScore = score;
                    thresholdScore = maxScore / std::exp(this->bandingOptions_.ScoreDiff);
                }

            }
            
            endRow = i;
            alpha.FinishEditingColumn(j, beginRow, endRow);
            
            // Now, revise the hints to tell the caller where the mass of the
            // distribution really lived in this column.
            hintEndRow = endRow;
            for (i = beginRow; i < endRow && alpha(i, j) < thresholdScore; ++i);
            hintBeginRow = i;
        }
        /* Now fill out the probability in the last pinned position.
         * We require that we end in a match. 
         * search for the term EDGE_CONDITION to find a comment with more information */
        auto likelihood = alpha(I - 1, J - 1) * e.Match_Just_Emission(I - 1, J - 1);
        alpha.StartEditingColumn(J, I, I + 1);
        alpha.Set(I, J, likelihood);
        alpha.FinishEditingColumn(J, I, I + 1);
        //DumpAlphaMatrix(alpha);
        
    }


    /**
     Fill the Beta matrix, the backwards half of the forward-backward algorithm.  
     This represents the probability that starting from the (i,j) state, the combined
     probability of transitioning out and following all paths through to the end.  
     
     In combination with the Alpha matrix, this allows us to calculate all paths that 
     pass through the (i,j) element, as exp(Alpha(i,j) + Beta(i,j))
     
     All probabilities stored in the matrix are stored as LOG probabilities.
     
     @param e The evaluator, such as QvEvaluator
     @param M the guide matrix for banding (this needs more documentation)
     @param beta The Beta matrix, stored as either a DenseMatrix or a SparseMatrix.
     */

    template<typename M, typename E, typename C>
    void
    SimpleRecursor<M, E, C>::FillBeta(const E& e, const M& guide, M& beta) const
    {
        int I = e.ReadLength();
        int J = e.TemplateLength();

        assert(beta.Rows() == I + 1 && beta.Columns() == J + 1);
        assert(guide.IsNull() ||
               (guide.Rows() == beta.Rows() && guide.Columns() == beta.Columns()));

        //Setup initial condition, at the end we are one
        beta.StartEditingColumn(J, I, I+1);
        beta.Set(I, J, 1.0);
        beta.FinishEditingColumn(J, I, I + 1);
        
        // Totally arbitray decision here...
        int hintBeginRow = std::max(I - 45, 0), hintEndRow = I;

        for (int j = (J-1); j > 0; --j)
        {
            this->RangeGuide(j, guide, beta, &hintBeginRow, &hintEndRow);

            int requiredBeginRow = max(0, hintBeginRow);

            beta.StartEditingColumn(j, hintBeginRow, hintEndRow);

            int i;
            double score = 0.0;
            double thresholdScore = 0.0;
            double maxScore = 0.0;
            
            int beginRow, endRow = hintEndRow;
            for (i = endRow - 1;
                 i > 0 && (score >= thresholdScore || i >= requiredBeginRow);
                 --i)
            {
                double thisMoveScore;
                score = 0.0;

                // Match:
                // TODO: Right now it assumes the probability of a match transition is 1.
                // Which might introduce problems in homopolymers on boundaries.
                
                if (i == (I-1) && (j == (J-1)))
                {
                    thisMoveScore = beta(i+1, j+1) * e.Match_Just_Emission(i,j);
                    score = C::Combine(score, thisMoveScore); // TODO: Redundant on first pass?
                }
                else if (i < (I - 1))
                {
                    thisMoveScore = beta(i + 1, j + 1) * e.Match(i, j);
                    score = C::Combine(score, thisMoveScore);
                }

                // Stick or Branch:
                if (i < (I - 1) && i > 0) // Can only transition to an insertion for the 2nd to last read base
                {
                    thisMoveScore = beta(i + 1, j) * e.Insertion(i, j - 1);
                    score = C::Combine(score, thisMoveScore);
                }

                // Deletion:
                if (j < (J - 1) && j > 0)
                {
                    thisMoveScore = beta(i, j + 1) * e.Deletion(j - 1);
                    score = C::Combine(score, thisMoveScore);
                }


                //  Save score
                beta.Set(i, j, score);

                if (score > maxScore)
                {
                    maxScore = score;
                    thresholdScore = maxScore / std::exp(this->bandingOptions_.ScoreDiff);
                }
            }

            beginRow = i + 1;
            beta.FinishEditingColumn(j, beginRow, endRow);
            //DumpBetaMatrix(beta);
            // Now, revise the hints to tell the caller where the mass of the
            // distribution really lived in this column.
            hintBeginRow = beginRow;
            for (i = endRow;
                 i > beginRow && beta(i - 1, j) < thresholdScore;
                 --i);
            hintEndRow = i;
        }
        
        beta.StartEditingColumn(0, 0, 1);
        /* Now to fill the top row which must be a match
         * search for the term EDGE_CONDITION to find a comment with more information */
        beta.Set(0, 0, e.Match_Just_Emission(0,0) * beta(1, 1));
        beta.FinishEditingColumn(0, 0, 1);
    }

    /// Calculate the recursion score by "stitching" together partial
    /// alpha and beta matrices.  alphaColumn, betaColumn, and
    /// absoluteColumn all refer to the same logical position in the
    /// template, but may have different values if, for instance,
    /// alpha here is a sub-range of the columns of the full alpha
    /// matrix.  Columns betaColumn and betaColumn + 1 of beta will be
    /// read; columns alphaColumn - 1 and alphaColumn - 2 of alpha
    /// will be read.
    template<typename M, typename E, typename C>
    double
    SimpleRecursor<M, E, C>::LinkAlphaBeta(const E& e,
                                           const M& alpha, int alphaColumn,
                                           const M& beta, int betaColumn,
                                           int absoluteColumn) const
    {
        const int I = e.ReadLength();

        assert(alphaColumn > 1 && absoluteColumn > 1);
        assert(absoluteColumn < e.TemplateLength());

        int usedBegin, usedEnd;
        boost::tie(usedBegin, usedEnd) = \
            RangeUnion(alpha.UsedRowRange(alphaColumn - 2),
                       alpha.UsedRowRange(alphaColumn - 1),
                       beta.UsedRowRange(betaColumn),
                       beta.UsedRowRange(betaColumn + 1));

        double v = 0.0, thisMoveScore;

        for (int i = usedBegin; i < usedEnd; i++)
        {
            if (i < I)
            {
                // Incorporate
                thisMoveScore = alpha(i, alphaColumn - 1) *
                                e.Match(i, absoluteColumn - 1) *
                                beta(i + 1, betaColumn);
                v = C::Combine(v, thisMoveScore);

                }

            // Delete:
            thisMoveScore = alpha(i, alphaColumn - 1) *
                            e.Deletion(absoluteColumn - 2) *
                            beta(i, betaColumn);
            v = C::Combine(v, thisMoveScore);
        }

        return v;
    }


    /**
     This method extends that Alpha matrix into a temporary matrix given by 
     ext.  It extends the region [beginColumn, beginColumn + numExtColumns)
     
     All new data is placed in the extension matrix.  The guesses for start/end
     rows in the banding are determined by evaluating neighbors of each position.
     @param <#parameter#>
     @returns <#retval#>
     @exception <#throws#>
     */
    template<typename M, typename E, typename C>
    void
    SimpleRecursor<M, E, C>::ExtendAlpha(const E& e,
                                         const M& alpha, int beginColumn,
                                         M& ext, int numExtColumns) const
    {
        assert(numExtColumns >= 2); // We have to fill at least one
        assert(alpha.Rows() == e.ReadLength() + 1 &&
               ext.Rows() == e.ReadLength() + 1); // The read never mutates

        // The new template may not be the same length as the old template.
        // Just make sure that we have anough room to fill out the extend buffer
        assert(beginColumn + 1 < e.TemplateLength() + 1);
        assert(ext.Columns() >= numExtColumns);
        assert(beginColumn >= 2);

        for (int extCol = 0; extCol < numExtColumns; extCol++)
        {
            int j = beginColumn + extCol;
            int beginRow, endRow;

            //
            // If this extend is contained within the column bounds of
            // the original alpha, we use the row range that was
            // previously determined.  Otherwise start at alpha's last
            // UsedRow beginRow and go to the end.
            //
            // BULLSHIT! If there was a deletion or insertion, the row range for the previous
            // column, not the column of interest will be used.
            // TODO: ERROR! Fix this. Temporary hack is to merge the columns in front and behind.
            // Still totally broken.
            if (j < e.TemplateLength())
            {
                boost::tie(beginRow, endRow) = alpha.UsedRowRange(j);
                int pBegin, pEnd, nBegin, nEnd;
                if ( (j-1) >= 0) {
                    boost::tie(pBegin, pEnd) = alpha.UsedRowRange(j-1);
                    beginRow = std::min(beginRow, pBegin);
                    endRow = std::max(endRow, pEnd);
                }
                if ( (j+1) < e.TemplateLength() )  {
                    boost::tie(nBegin, nEnd) = alpha.UsedRowRange(j+1);
                    beginRow = std::min(beginRow, nBegin);
                    endRow = std::max(endRow, nEnd);
                }
            }
            else
            {
                beginRow = alpha.UsedRowRange(alpha.Columns() - 1).Begin;
                endRow = alpha.Rows();
            }
            
       

            ext.StartEditingColumn(extCol, beginRow, endRow);

            int i;
            double score;
            // Due to pinning at the end, moves are only possible if less than these positions.
            int maxLeftMovePossible = e.TemplateLength();
            int maxDownMovePossible = e.ReadLength();

            for (i = beginRow; i < endRow; i++)
            {
                double thisMoveScore = 0.0;
                score = 0.0;
                

                // Match:
                if (i > 0 && j > 0)
                {
                    double prev = extCol == 0 ?
                            alpha(i - 1, j - 1) :
                            ext(i - 1, extCol - 1);
                    if (i == 1 && j == 1) { //TODO: Remove this branch bottleneck...
                        thisMoveScore = e.Match_Just_Emission(0,0); // prev should be zero, so no need for explicit prev + e.Match_Just_Emission
                    }
                    else if (i == maxDownMovePossible && j == maxLeftMovePossible) {
                        thisMoveScore = prev * e.Match_Just_Emission(i -1 ,j - 1);
                    }
                    else if (i < maxDownMovePossible && j < maxLeftMovePossible) {
                        thisMoveScore = prev * e.Match(i - 1, j - 1);
                    }
                    score = C::Combine(score, thisMoveScore);
                }

                // Stick or Branch:
                if (i > 1 && i < maxDownMovePossible && j != maxLeftMovePossible)
                {
                    thisMoveScore = ext(i - 1, extCol) * e.Insertion(i - 1, j - 1);
                    score = C::Combine(score, thisMoveScore);
                }

                // Delete:
                if (j > 1 && j < maxLeftMovePossible && i != maxDownMovePossible)
                {
                    double prev = extCol == 0 ?
                            alpha(i, j - 1) :
                            ext(i, extCol - 1);
                    thisMoveScore = prev * e.Deletion(j - 2);
                    score = C::Combine(score, thisMoveScore);
                }
              
                ext.Set(i, extCol, score);
            }
            assert (i == endRow);
            ext.FinishEditingColumn(extCol, beginRow, endRow);
        }
    }


    // Semantic: After ExtendBeta(B, j), we have
    //    ext(:, numExtColumns-1) = B'(:,j)
    //    ext(:, numExtColumns-2) = B'(:,j-1) ...
    //
    // Note: lastColumn is the numerically largest column number that
    // will be filled, but it is filled first since beta fill is done
    // backwards.
    //
    // Accesses B(:, ..(j+2))
    template<typename M, typename E, typename C>
    void
    SimpleRecursor<M, E, C>::ExtendBeta(const E& e,
                                        const M& beta, int lastColumn,
                                        M& ext, int numExtColumns,
                                        int lengthDiff) const
    {
        int I = e.ReadLength();
        int J = e.TemplateLength();
        
        int firstColumn = 0 - lengthDiff;
        int lastExtColumn = numExtColumns - 1;

        assert(beta.Rows() == I + 1 && ext.Rows() == I + 1);

        // The new template may not be the same length as the old template.
        // Just make sure that we have enough room to fill out the extend buffer
        assert(lastColumn + 2 <= J);
        assert(lastColumn >= 0);
        assert(ext.Columns() >= numExtColumns);

        for (int j = lastColumn; j > lastColumn - numExtColumns; j--)
        {
            int jp = j + lengthDiff; // Convert from old template to new template coordinates.
            int extCol = lastExtColumn - (lastColumn - j);
            int beginRow, endRow;

            if (j < 0)
            {
                beginRow = 0;
                endRow = beta.UsedRowRange(0).End;
            }
            else
            {
                boost::tie(beginRow, endRow) = beta.UsedRowRange(j);
                int pBegin, pEnd, nBegin, nEnd;
                if ( (j-1) >= 0) {
                    boost::tie(pBegin, pEnd) = beta.UsedRowRange(j-1);
                    beginRow = std::min(beginRow, pBegin);
                    endRow = std::max(endRow, pEnd);
                }
                if ( (j+1) < e.TemplateLength() ) {
                    boost::tie(nBegin, nEnd) = beta.UsedRowRange(j+1);
                    beginRow = std::min(beginRow, nBegin);
                    endRow = std::max(endRow, nEnd);
                }
            }

            ext.StartEditingColumn(extCol, beginRow, endRow);

            int i;
            double score;

            for (i = endRow - 1;
                 i >= beginRow;
                 i--)
            {
                double thisMoveScore = 0.0;
                score = 0.0;

                // Incorporation:
                if (i < I && j < J)
                {
                    double prev = (extCol == lastExtColumn) ?
                        beta(i + 1, j + 1) :
                        ext(i + 1, extCol + 1);
                    if ((i == (I-1) && jp == (J-1)) || (i==0 && j == firstColumn)) {
                        thisMoveScore = prev * e.Match_Just_Emission(i, jp);
                    }
                    else if ( j > firstColumn && i > 0) {
                        thisMoveScore = prev * e.Match(i, jp);
                    }
                    score = C::Combine(score, thisMoveScore);
                }

                // Stick or branch
                if (i < (I-1) && i > 0 && j > firstColumn)
                {
                    thisMoveScore = ext(i + 1, extCol) * e.Insertion(i, jp - 1);
                    score = C::Combine(score, thisMoveScore);
                }

                // Deletion
                if (j < (J-1) && j > firstColumn && i > 0)
                {
                    double prev = (extCol == lastExtColumn) ?
                        beta(i, j + 1) :
                        ext(i, extCol + 1);
                    thisMoveScore = prev * e.Deletion(jp - 1);
                    score = C::Combine(score, thisMoveScore);
                }

                ext.Set(i, extCol, score);
            }
            ext.FinishEditingColumn(extCol, beginRow, endRow);
        }
    }


    template<typename M, typename E, typename C>
    SimpleRecursor<M, E, C>::SimpleRecursor(const BandingOptions& banding)
        : detail::RecursorBase<M, E, C>(banding)
    {}


    template class SimpleRecursor<DenseMatrix, QvEvaluator, detail::ViterbiCombiner>;
    template class SimpleRecursor<DenseMatrix, QvEvaluator, detail::SumProductCombiner>;
    template class SimpleRecursor<SparseMatrix, QvEvaluator, detail::ViterbiCombiner>;
    template class SimpleRecursor<SparseMatrix, QvEvaluator, detail::SumProductCombiner>;

}
