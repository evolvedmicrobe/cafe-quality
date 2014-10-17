// Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
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

// Author: David Alexander, Lance Hepler

#include <algorithm>
#include <boost/functional/hash.hpp>
#include <boost/tuple/tuple.hpp>
#include <cmath>
#include <utility>

#include "Quiver/MultiReadMutationScorer.hpp"
#include "Quiver/QuiverConsensus.hpp"
#include "Quiver/MutationEnumerator.hpp"
#include "Mutation.hpp"

#include "Utils.hpp"

#include "Logging/Logging.hpp"


namespace ConsensusCore
{
    using std::vector;

    namespace { // PRIVATE
    using std::max_element;

    struct RefineDinucleotideRepeatOptions : RefineOptions
    {
        RefineDinucleotideRepeatOptions(int minDinucleotideRepeatElements)
            : MinDinucleotideRepeatElements(minDinucleotideRepeatElements)
        {
            MaximumIterations = 1;
        }

        int MinDinucleotideRepeatElements;
    };

    vector<ScoredMutation>
    DeleteRange(vector<ScoredMutation> input, int rStart, int rEnd)
    {
        vector<ScoredMutation> output;
        foreach (ScoredMutation s, input)
        {
            int pos = s.Start();
            if (!(rStart <= pos && pos <= rEnd))
            {
                output.push_back(s);
            }
        }
        return output;
    }

    bool ScoreComparer(const ScoredMutation& i, const ScoredMutation& j)
    {
        return i.Score() < j.Score();
    }

    //    Given a list of (mutation, score) tuples, this utility method
    //    greedily chooses the highest scoring well-separated elements.  We
    //    use this to avoid applying adjacent high scoring mutations, which
    //    are the rule, not the exception.  We only apply the best scoring one
    //    in each neighborhood, and then revisit the neighborhoods after
    //    applying the mutations.
    //
    //    This is highly unoptimized.  It is not in the critical path.
    vector<ScoredMutation>
    BestSubset(vector<ScoredMutation> input, int mutationSeparation)
    {
        if (mutationSeparation == 0)
            return input;

        vector<ScoredMutation> output;

        while (!input.empty())
        {
            ScoredMutation& best = *max_element(input.begin(), input.end(), ScoreComparer);
            output.push_back(best);
            int nStart = best.Start() - mutationSeparation;
            int nEnd = best.Start() + mutationSeparation;
            input = DeleteRange(input, nStart, nEnd);
        }

        return output;
    }


    // Sadly and annoyingly there is no covariance on std::vector in C++, so we have
    // to explicitly project back down to the superclass type to use the APIs as written.
    vector<Mutation>
    ProjectDown(const vector<ScoredMutation>& smuts)
    {
        return vector<Mutation>(smuts.begin(), smuts.end());
    }


    int ProbabilityToQV(double probability, int cap = 93)
    {
        using std::min;

        if (probability <= 0.0)
            return cap;

        return min(cap, static_cast<int>(round(-10.0 * log10(probability))));
    }

    template <typename E, typename O>
    E MutationEnumerator(const std::string& tpl, const O& opts)
    {
        return E(tpl);
    }

    //
    // this MUST go last to properly specialize the MutationEnumerator
    //
    template <>
    DinucleotideRepeatMutationEnumerator
    MutationEnumerator<>(const std::string& tpl, const RefineDinucleotideRepeatOptions& opts)
    {
        return DinucleotideRepeatMutationEnumerator(tpl, opts.MinDinucleotideRepeatElements);
    }

    template <typename E, typename O>
    bool AbstractRefineConsensus(AbstractMultiReadMutationScorer& mms, const O& opts)
    {
        bool isConverged = false;
        float score = mms.BaselineScore();
        boost::hash<std::string> hash;
        std::set<size_t> tplHistory;

        vector<ScoredMutation> favorableMutsAndScores;

        for (int iter = 0; iter < opts.MaximumIterations; iter++)
        {
            LDEBUG << "Round " << iter;
            LDEBUG << "State of MMS: " << std::endl << mms.ToString();

            if (tplHistory.find(hash(mms.Template())) != tplHistory.end())
            {
                LDEBUG << "Cycle detected!";
            }

            if (mms.BaselineScore() < score)
            {
                LDEBUG << "Score decrease"; // Usually recoverable, so we allow iteration to continue
            }
            score = mms.BaselineScore();

            //
            // Try all mutations in iteration 0.  In subsequent iterations, try mutations
            // nearby those used in previous iteration.
            //
            E mutationEnumerator = MutationEnumerator<E, O>(mms.Template(), opts);
            vector<Mutation> mutationsToTry;
            if (iter == 0) {
                mutationsToTry = mutationEnumerator.Mutations();
            }
            else {
                mutationsToTry = UniqueNearbyMutations(mutationEnumerator,
                                                       ProjectDown(favorableMutsAndScores),
                                                       opts.MutationNeighborhood);
            }

            //
            // Screen for favorable mutations.  If none, we are done (converged).
            //
            favorableMutsAndScores.clear();
            foreach (const Mutation& m, mutationsToTry)
            {
                if (mms.FastIsFavorable(m)) {
                    float mutScore = mms.Score(m);
                    favorableMutsAndScores.push_back(m.WithScore(mutScore));
                }
            }
            if (favorableMutsAndScores.empty())
            {
                isConverged = true;
                break;
            }

            //
            // Go with the "best" subset of well-separated high scoring mutations
            //
            vector<ScoredMutation> bestSubset = BestSubset(favorableMutsAndScores, opts.MutationSeparation);

            //
            // Attempt to avoid cycling.  We could do a better job here.
            //
            if (bestSubset.size() > 1)
            {
                std::string nextTpl = ApplyMutations(ProjectDown(bestSubset), mms.Template());
                if (tplHistory.find(hash(nextTpl)) != tplHistory.end())
                {
                    LDEBUG << "Attempting to avoid cycle";
                    bestSubset = std::vector<ScoredMutation>(bestSubset.begin(), bestSubset.begin() + 1);
                }
            }

            LDEBUG << "Applying mutations:";
            foreach (const ScoredMutation& smut, bestSubset)
            {
                LDEBUG << "\t" << smut;
            }

            tplHistory.insert(hash(mms.Template()));
            mms.ApplyMutations(ProjectDown(bestSubset));
        }

        return isConverged;
    }
    } // PRIVATE


    bool RefineConsensus(AbstractMultiReadMutationScorer& mms, const RefineOptions& opts)
    {
        return AbstractRefineConsensus<UniqueSingleBaseMutationEnumerator>(mms, opts);
    }


    void RefineDinucleotideRepeats(AbstractMultiReadMutationScorer& mms, int minDinucleotideRepeatElements)
    {
        RefineDinucleotideRepeatOptions opts(minDinucleotideRepeatElements);
        AbstractRefineConsensus<DinucleotideRepeatMutationEnumerator>(mms, opts);
    }


    std::vector<int> ConsensusQVs(AbstractMultiReadMutationScorer& mms)
    {
        std::vector<int> QVs;
        UniqueSingleBaseMutationEnumerator mutationEnumerator(mms.Template());
        for (size_t pos = 0; pos < mms.Template().length(); pos++)
        {
            double scoreSum = 0.0;
            foreach (const Mutation& m, mutationEnumerator.Mutations(pos, pos + 1))
            {
                scoreSum += exp(mms.FastScore(m));
            }
            QVs.push_back(ProbabilityToQV(1.0 - 1.0 / (1.0 + scoreSum)));
        }
        return QVs;
    }


#if 0
    Matrix<float> MutationScoresMatrix(mms)
    {
        NotYetImplemented();
    }


    Matrix<float> MutationScoresMatrix(mms, mutationsToScore)
    {
        NotYetImplemented();
    }
#endif
}
