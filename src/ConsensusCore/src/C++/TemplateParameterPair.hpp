//
//  TemplateParameterPair.h
//  ConsensusCoreXCode
//
//  Created by Nigel Delaney on 3/3/15.
//  Copyright (c) 2015 Pacific Biosciences. All rights reserved.
//

#ifndef __ConsensusCoreXCode__TemplateParameterPair__
#define __ConsensusCoreXCode__TemplateParameterPair__

#include <stdio.h>
#include <string>
#include <vector>
#include <memory>

#include "TransitionParameters.hpp"
#include "ContextParameters.hpp"

namespace ConsensusCore {
    
    struct TemplateParameterPair {
    public:
        std::string tpl;
        std::vector<TransitionParameters> trans_probs;
        
        TemplateParameterPair(const std::string& tpl_,
                              const std::vector<TransitionParameters>& trans_probs_)
        : tpl(tpl_)
        , trans_probs(trans_probs_)
        {}
        
        TemplateParameterPair()
        : tpl()
        , trans_probs()
        {}
        
        TemplateParameterPair(const std::string& tpl_,
                              const ContextParameters& ctx) : tpl(tpl_), trans_probs (tpl.size()-1)
        {
            // Initialize the dinucleotide context values.
            for(int i = 0; i < (tpl.size() -1); i ++)
            {
                auto ps = ctx.GetParametersForContext(tpl.at(i), tpl.at(i+1));
                trans_probs[i] = ps;
            }
        }
        
        /* Get a subsection of this template parameter pair
         TODO: This is a brutal copy operation */
        TemplateParameterPair GetSubSection(int start, int end) const {
            auto starti = trans_probs.begin() + start;
            auto endi = trans_probs.begin() + end;
            return TemplateParameterPair(tpl.substr(start, end), std::vector<TransitionParameters>(starti, endi));
        }
        
    };
}

#endif /* defined(__ConsensusCoreXCode__TemplateParameterPair__) */
