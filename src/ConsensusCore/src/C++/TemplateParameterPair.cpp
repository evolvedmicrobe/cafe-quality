//
//  TemplateParameterPair.cpp
//  ConsensusCoreXCode
//
//  Created by Nigel Delaney on 3/10/15.
//  Copyright (c) 2015 Pacific Biosciences. All rights reserved.
//

#include <stdio.h>
#include "TemplateParameterPair.hpp"

namespace ConsensusCore {
    
    
    TemplateParameterPair::TemplateParameterPair(const std::string& tpl_,
                          const std::vector<TransitionParameters>& trans_probs_)
    : tpl(tpl_), trans_probs(trans_probs_)
    {
        assert(tpl.size() == (trans_probs.size() + 1));
    }
    
    TemplateParameterPair::TemplateParameterPair()
    : tpl()
    , trans_probs()
    {}
    
    TemplateParameterPair::TemplateParameterPair(const TemplateParameterPair& other) :
    tpl(other.tpl), trans_probs(other.trans_probs)
    {
        assert(tpl.size() == (trans_probs.size() + 1));
    }
    
    TemplateParameterPair::TemplateParameterPair(const std::string& tpl_,
                          const ContextParameters& ctx) : tpl(tpl_), trans_probs (tpl_.size()-1)
    {
        // Initialize the dinucleotide context values.
        for(int i = 0; i < (tpl.size() -1); i ++)
        {
            auto ps = ctx.GetParametersForContext(tpl.at(i), tpl.at(i+1));
            trans_probs[i] = ps;
        }
        assert(tpl.size() == (trans_probs.size() + 1));
    }
    
    
    
    
    TemplateParameterPair
    TemplateParameterPair::GetSubSection(int start, int len) const {
        auto starti = trans_probs.begin() + start;
        auto endi = trans_probs.begin() + start + len;
        // Remove one as the transition parameters need to be one smaller.
        return TemplateParameterPair(tpl.substr(start, len), std::vector<TransitionParameters>(starti, endi - 1 ));
    }
    
    

    
    
    
}