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

// Author: Patrick Marks, David Alexander

#include "Mutation.hpp"
#include <cassert>


using std::max;
using namespace std;

namespace ConsensusCore
{
    std::string
    Mutation::ToString() const
    {
        using boost::str;
        using boost::format;

        switch (Type())
        {
            case INSERTION:
                return str(format("Insertion (%s) @%d") % newBases_ % start_);
            case DELETION:
                return str(format("Deletion @%d:%d") % start_ % end_);
            case SUBSTITUTION:
                return str(format("Substitution (%s) @%d:%d") % newBases_ % start_ % end_);
            default: ShouldNotReachHere();
        }
    }

    std::ostream& operator<<(std::ostream& out, const Mutation& m)
    {
        out << m.ToString();
        return out;
    }


    ScoredMutation Mutation::WithScore(double score) const
    {
        return ScoredMutation(*this, score);
    }

    // NOTE: Start is not just equal to Mut.Start() here because the location of a mutation can change
    // as earlier mutations are applied.
    static void
    _ApplyMutationInPlace(const Mutation& mut, int start, TemplateParameterPair& tpl, const ContextParameters& ctx_params)
    {
        if (mut.IsSubstitution())
        {
            tpl.tpl.replace(start, mut.End() - mut.Start(), mut.NewBases());
            if ((start + 1) < tpl.tpl.length()) {
                tpl.trans_probs[start] = ctx_params.GetParametersForContext(tpl.tpl.at(start), tpl.tpl.at(start+1));
            }
            if (start > 0) {
                tpl.trans_probs[start-1] = ctx_params.GetParametersForContext(tpl.tpl.at(start -1), tpl.tpl.at(start));
            }
        }
        else if (mut.IsDeletion())
        {
            tpl.tpl.erase(start, mut.End() - mut.Start());
            auto maxEnd = tpl.tpl.length() - 1; // Only the second to last base can have an update
            // Only update if not first base though
            if(start > 0 && start < maxEnd) {
                tpl.trans_probs[start-1] = ctx_params.GetParametersForContext(tpl.tpl.at(start-1), tpl.tpl.at(start));
            }
            if (start < maxEnd ) {
                tpl.trans_probs.erase(tpl.trans_probs.begin() + start, tpl.trans_probs.begin() + start + ( mut.End()- mut.Start()));                
            }
        }
        else if (mut.IsInsertion())
        {
            tpl.tpl.insert(start, mut.NewBases());
            if (start > 0) {
                tpl.trans_probs[start-1] = ctx_params.GetParametersForContext(tpl.tpl.at(start-1), tpl.tpl.at(start));
            }
            if(start < (tpl.tpl.length()-1))
            {
                auto new_params = ctx_params.GetParametersForContext(tpl.tpl.at(start), tpl.tpl.at(start+1));
                tpl.trans_probs.insert(tpl.trans_probs.begin() + start, new_params);
            }
        }
    }

    TemplateParameterPair
    ApplyMutation(const Mutation& mut, const TemplateParameterPair& tpl, const ContextParameters& ctx_params)
    {
        auto tplCopy = string(tpl.tpl);
        auto new_probs = vector<TransitionParameters>(tpl.trans_probs);
        TemplateParameterPair new_tpl;
        new_tpl.tpl = tplCopy;
        new_tpl.trans_probs = new_probs;
        _ApplyMutationInPlace(mut, mut.Start(), new_tpl, ctx_params);
        return new_tpl;
    }

    TemplateParameterPair
    ApplyMutations(const std::vector<Mutation>& muts, const TemplateParameterPair& tpl,  const ContextParameters& ctx_params)
    {
        // Make a new copy
        auto tplCopy = string(tpl.tpl);
        auto new_probs = vector<TransitionParameters>(tpl.trans_probs);
        TemplateParameterPair new_tpl (tplCopy, new_probs);
        
        // Apply mutations
        std::vector<Mutation> sortedMuts(muts);
        std::sort(sortedMuts.begin(), sortedMuts.end());
        int runningLengthDiff = 0;
        foreach (const Mutation& mut, sortedMuts)
        {
            _ApplyMutationInPlace(mut, mut.Start() + runningLengthDiff, new_tpl, ctx_params);
            runningLengthDiff += mut.LengthDiff();
        }
        return new_tpl;
    }

    std::string MutationsToTranscript(const std::vector<Mutation>& mutations,
                                      const std::string& tpl)
    {
        std::vector<Mutation> sortedMuts(mutations);
        std::sort(sortedMuts.begin(), sortedMuts.end());

        // Build an alignnment transcript corresponding to these mutations.
        int tpos = 0;
        std::string transcript = "";
        foreach (const Mutation& m, sortedMuts)
        {
            for (; tpos < m.Start(); ++tpos)
            {
                transcript.push_back('M');
            }

            if (m.IsInsertion())
            {
                transcript += std::string(m.LengthDiff(), 'I');
            }
            else if (m.IsDeletion())
            {
                transcript += std::string(-m.LengthDiff(), 'D');
                tpos += -m.LengthDiff();
            }
            else if (m.IsSubstitution())
            {
                int len = m.End() - m.Start();
                transcript += std::string(len, 'R');
                tpos += len;
            }
            else
            {
                ShouldNotReachHere();
            }
        }
        for (; tpos < (int)tpl.length(); ++tpos)
        {
            transcript.push_back('M');
        }
        return transcript;
    }

    // MutatedTemplatePositions:
    //  * Returns a vector of length (tpl.length()+1), which, roughly speaking,
    //    indicates the positions in the mutated template tpl' of the characters
    //    in tpl.
    //  * More precisely, for any slice [s, e) of tpl, letting:
    //      - t[s, e) denote the induced substring of the template;
    //      - m[s, e) denote the subvector of mutations with Position
    //        in [s, e);
    //      - t' denote the mutated template; and
    //      - t[s, e)' denote the result of applying mutation m[s, e) to t[s, e),
    //    the resultant vector mtp satisfies t'[mtp[s], mtp[e]) == t[s,e)'.
    //  * Example:
    //               01234567                           0123456
    //              "GATTACA" -> (Del T@2, Ins C@5) -> "GATACCA";
    //    here mtp = 01223567, which makes sense, because for instance:
    //      - t[0,3)=="GAT" has become t'[0,2)=="GA";
    //      - t[0,2)=="GA"  remains "GA"==t'[0,2);
    //      - t[4,7)=="ACA" has become t[3,7)=="ACCA",
    //      - t[5,7)=="CA"  remains "CA"==t'[5,7).
    //
    std::vector<int> TargetToQueryPositions(const std::vector<Mutation>& mutations,
                                            const std::string& tpl)
    {
        return TargetToQueryPositions(MutationsToTranscript(mutations, tpl));
    }


    ScoredMutation::ScoredMutation(const Mutation& m, double score)
        : Mutation(m),
          score_(score)
    {}

    ScoredMutation::ScoredMutation()
        : Mutation(),
          score_(0)
    {}

    double ScoredMutation::Score() const
    {
        return score_;
    }

    std::string ScoredMutation::ToString() const
    {
        std::stringstream ss;
        ss << *this;
        return ss.str();
    }

    std::ostream& operator<<(std::ostream& out, const ScoredMutation& m)
    {
        out << m.Mutation::ToString() << " " << boost::format("%0.2f") % m.Score();
        return out;
    }
}
